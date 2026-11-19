# =============================================================================
# STEP 6 — type I error at the BE boundary, PBPK vs PopPK.
#
# In : data/pbpk_endpoints.csv, results/<model>/step3_poppk_endpoints.csv
# Out: results/<model>/step6_boundaries.csv    T50 where the true GMR hits 0.80 / 1.25
#      results/<model>/step6_type1_error.csv   pass rate at that T50 (target 5%)
#      results/<model>/step6_tost_control.csv  pass rate with the reference x 0.80
#      results/<model>/fig_step6_type1_bars.png, fig_step6_oc.png
#
# 1. Boundary: walk out from the reference T50 on a 0.005 h grid until the
#    noise-free GMR leaves 0.80–1.25 (Cmax and AUC separately; Joint = the
#    one reached first).
# 2. Type I error = share of 20,000 simulated BE trials that pass AT the
#    boundary. A correct test gives about 5%.
# 3. Cross-check: the pass rate each method gives at the boundary the OTHER
#    method found. PBPK is the truth, so "PBPK at the PopPK boundary" is the
#    type I error a PopPK-based safe space would really have.
# 4. Control: test = reference x 0.80 exactly, which checks the TOST itself.
# =============================================================================

using CSV, DataFrames, CairoMakie
include("settings.jl")
include("be_functions.jl")

pbpk = CSV.read(joinpath(DATA, "pbpk_endpoints.csv"), DataFrame)
const GRID = collect(range(minimum(pbpk.T50), maximum(pbpk.T50); step = T50_STEP_TYPE1))
COLOR = Dict("PBPK" => :firebrick, "PopPK" => :royalblue)

"Wilson 95% interval of a simulated proportion (Monte Carlo uncertainty)."
function wilson(p, n; z = 1.96)
    d = 1 + z^2 / n; c = (p + z^2 / (2n)) / d
    h = z * sqrt(p * (1 - p) / n + z^2 / (4n^2)) / d
    (max(0.0, c - h), min(1.0, c + h))
end

for model_name in MODELS
    dir   = joinpath(RESULTS, model_name)
    poppk = CSV.read(joinpath(dir, "step3_poppk_endpoints.csv"), DataFrame)
    sources = ["PBPK" => pbpk, "PopPK" => poppk]

    # ---- 1. boundaries ------------------------------------------------------
    bnd = DataFrame()
    for (source, pk) in sources
        gmr = Dict(ep => [gmr_at(pk, ep, t) for t in GRID] for ep in ENDPOINTS)
        for side in (:low, :high)
            hits = Dict(ep => first_crossing(GRID, gmr[ep], side) for ep in ENDPOINTS)
            for ep in ENDPOINTS
                push!(bnd, (endpoint = ep, source = source, side = String(side),
                            T50 = hits[ep].t50, limit = hits[ep].limit, binding = ep); promote = true)
            end
            reached = [ep for ep in ENDPOINTS if !ismissing(hits[ep].t50)]
            if isempty(reached)
                push!(bnd, (endpoint = "Joint", source = source, side = String(side),
                            T50 = missing, limit = "not reached", binding = "none"); promote = true)
            else
                first_ep = reached[argmin([abs(hits[ep].t50 - REF_T50) for ep in reached])]
                push!(bnd, (endpoint = "Joint", source = source, side = String(side),
                            T50 = hits[first_ep].t50, limit = hits[first_ep].limit, binding = first_ep); promote = true)
            end
        end
    end
    CSV.write(joinpath(dir, "step6_boundaries.csv"), bnd)

    # ---- 2–3. type I error at own and at the other method's boundary --------
    t1 = DataFrame()
    for ep in [ENDPOINTS; "Joint"], side in ("low", "high"), n in N_SUBJ, (source, pk) in sources
        own   = only(bnd[(bnd.endpoint .== ep) .& (bnd.side .== side) .& (bnd.source .== source), :])
        other = only(bnd[(bnd.endpoint .== ep) .& (bnd.side .== side) .& (bnd.source .!= source), :])
        ismissing(own.T50) && continue
        t1_own = pass_rate_at(pk, ep, own.T50, n, N_TRIAL_TYPE1)
        lo, hi = wilson(t1_own, N_TRIAL_TYPE1)
        t1_other  = ismissing(other.T50) ? missing : pass_rate_at(pk, ep, other.T50, n, N_TRIAL_TYPE1)
        gmr_other = ismissing(other.T50) ? missing : gmr_at(pk, own.binding, other.T50)
        push!(t1, (endpoint = ep, side = side, n_subj = n, judged_by = source,
                   own_boundary = own.T50, limit = own.limit, binding = own.binding,
                   type1_own = t1_own, ci_low = lo, ci_high = hi,
                   other_boundary = other.T50, type1_at_other_boundary = t1_other,
                   true_GMR_at_other_boundary = gmr_other); promote = true)
    end
    CSV.write(joinpath(dir, "step6_type1_error.csv"), t1)

    # ---- 4. TOST control: test = reference x 0.80 ---------------------------
    ctl = DataFrame()
    for ep in [ENDPOINTS; "Joint"], (source, pk) in sources, n in N_SUBJ
        rC, rA = values_at(pk, "Cmax", REF_T50), values_at(pk, "AUC_tEnd", REF_T50)
        p = ep == "Joint" ?
            pass_rate_joint(rC .* BE_LOWER, rC, rA .* BE_LOWER, rA; n_subj = n, n_trial = N_TRIAL_TYPE1,
                            cvC = WSV_CV["Cmax"], cvA = WSV_CV["AUC_tEnd"]) :
            pass_rate(values_at(pk, ep, REF_T50) .* BE_LOWER, values_at(pk, ep, REF_T50);
                      n_subj = n, n_trial = N_TRIAL_TYPE1, cv = WSV_CV[ep])
        push!(ctl, (endpoint = ep, source = source, n_subj = n, pass_rate = p))
    end
    CSV.write(joinpath(dir, "step6_tost_control.csv"), ctl)

    # ---- figures ------------------------------------------------------------
    d = t1[(t1.judged_by .== "PBPK") .& (t1.endpoint .== "Joint") .& .!ismissing.(t1.type1_at_other_boundary), :]
    fig = Figure(size = (900, 450))
    ax = Axis(fig[1, 1]; title = "Type I error, joint ABE, judged by PBPK ($model_name PopPK)",
              ylabel = "P(conclude BE)", xticks = (1:nrow(d), ["$(r.side)\nn=$(r.n_subj)" for r in eachrow(d)]))
    barplot!(ax, (1:nrow(d)) .- 0.2, d.type1_own; width = 0.4, color = COLOR["PBPK"], label = "at PBPK boundary")
    barplot!(ax, (1:nrow(d)) .+ 0.2, Float64.(d.type1_at_other_boundary); width = 0.4,
             color = COLOR["PopPK"], label = "at PopPK boundary")
    hlines!(ax, [0.05]; color = :black, linestyle = :dash)
    axislegend(ax; position = :lt)
    save(joinpath(dir, "fig_step6_type1_bars.png"), fig)

    # operating characteristic: joint pass rate vs the true Cmax GMR
    fig = Figure(size = (1100, 450))
    for (j, n) in enumerate(N_SUBJ)
        ax = Axis(fig[1, j]; title = "Joint ABE, n = $n ($model_name PopPK)",
                  xlabel = "true Cmax GMR", ylabel = "P(conclude BE)")
        vlines!(ax, [BE_LOWER, BE_UPPER]; color = :gray60, linestyle = :dash)
        hlines!(ax, [0.05]; color = :gray60, linestyle = :dot)
        for (source, pk) in sources
            pts = [t for t in GRID[1:8:end] if 0.70 <= gmr_at(pk, "Cmax", t) <= 1.35]
            g = [gmr_at(pk, "Cmax", t) for t in pts]
            p = [pass_rate_at(pk, "Joint", t, n, N_TRIAL_OC) for t in pts]
            o = sortperm(g)
            lines!(ax, g[o], p[o]; color = COLOR[source], label = source)
        end
        j == 1 && axislegend(ax; position = :rt)
    end
    save(joinpath(dir, "fig_step6_oc.png"), fig)

    println("\n$model_name — boundaries:"); show(bnd; allrows = true); println()
    println("\n$model_name — type I error:")
    show(t1[:, [:endpoint, :side, :n_subj, :judged_by, :own_boundary, :type1_own,
                :other_boundary, :type1_at_other_boundary]]; allrows = true); println()
end
