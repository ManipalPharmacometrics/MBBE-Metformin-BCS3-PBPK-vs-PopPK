# =============================================================================
# 05_gamma_comparison.jl  — compare the two γ (Weibull shape) choices in the
# JOINT PopPK model, against the PBPK truth:
#     (open) γ estimated as ONE shared value across arms
#     (γ=1)  γ fixed at 1.0 (= the PBPK dissolution β)
# Both use shared, formulation-independent disposition (CL, Vc) per physiology;
# only the release scale α varies with T50.
#
# Produces:
#   output/gamma_comparison_boundaries.csv   safe-space bounds: PBPK vs both modes
#   output/gamma_comparison_ivivc.csv        fitted α vs T50 for both modes
#   output/fig_gamma_comparison_<ep>.png     power curves PBPK vs open vs γ=1
#
# License-free (open packages). Run:
#   julia --project=nopumas_env scripts/05_gamma_comparison.jl
# =============================================================================

using CSV, DataFrames, Statistics, CairoMakie, Optim

const HERE = @__DIR__
const ROOT = normpath(joinpath(HERE, ".."))
const OUT  = joinpath(ROOT, "output"); mkpath(OUT)

include(joinpath(ROOT, "src", "weibull_absorption.jl"))   # simulate_weibull, fit_physiology_joint
include(joinpath(ROOT, "src", "nca_lite.jl"))             # nca_metrics
include(joinpath(ROOT, "src", "be_utils.jl")); using .BEUtils

const DENSE = collect(0.0:0.05:24.0)
const LAG   = 0.25

data = CSV.read(joinpath(OUT, "sparse_be_dataset.csv"), DataFrame)
dose = maximum(skipmissing(data.amt))
ids  = sort(unique(data.IndividualId))
t50s = sort(unique(data.T50))

"Joint per-physiology fit for a γ mode -> (endpoints DataFrame, α-by-T50 DataFrame)."
function fit_mode(gamma_fixed)
    endp = DataFrame(); alphas = Dict{Float64,Vector{Float64}}(t => Float64[] for t in t50s)
    for id in ids
        arms = [(t50 = t, tobs = sort(filter(r -> r.IndividualId==id && r.T50==t && r.evid==0, data), :time).time,
                 cobs = collect(skipmissing(sort(filter(r -> r.IndividualId==id && r.T50==t && r.evid==0, data), :time).dv)))
                for t in t50s]
        f = fit_physiology_joint(arms, dose; lag = LAG, gamma_fixed = gamma_fixed)
        for t in t50s
            push!(alphas[t], f.alpha[t])
            m = nca_metrics(DENSE, simulate_weibull(DENSE, f.alpha[t], f.gamma, f.CL, f.Vc, LAG, dose))
            push!(endp, (T50 = t, IndividualId = id, Cmax = m.Cmax,
                         AUC_inf = m.AUC_inf, AUC_tEnd = m.AUC_tEnd); cols = :union)
        end
    end
    gm(x) = exp(mean(log.(x)))
    ivivc = DataFrame(T50 = t50s, gm_alpha = [gm(alphas[t]) for t in t50s])
    endp, ivivc
end

# ---- safe-space helpers (same math as 04) -----------------------------------
metric_matrix(pk, metric, anchors, ids) = begin
    M = Matrix{Union{Missing,Float64}}(missing, length(ids), length(anchors))
    for (j, t) in enumerate(anchors)
        d = sort(filter(:T50 => ==(t), pk), :IndividualId)
        for (i, id) in enumerate(ids)
            r = findfirst(==(id), d.IndividualId); r !== nothing && (M[i, j] = d[r, metric])
        end
    end; M
end
interp_log(vals, anchors, t) = begin
    ok = .!ismissing.(vals); sum(ok) < 2 && return missing
    x = anchors[ok]; v = Float64.(vals[ok])
    t <= x[1] && return v[1]; t >= x[end] && return v[end]
    k = searchsortedlast(x, t); (t == x[k]) && return v[k]
    w = (log(t)-log(x[k]))/(log(x[k+1])-log(x[k])); exp((1-w)*log(v[k]) + w*log(v[k+1]))
end
safe_bounds(grid, power, ptgt) = begin
    above = power .>= ptgt; ri = argmin(abs.(grid .- REF_T50)); lo = ri; hi = ri
    while lo > 1 && above[lo-1]; lo -= 1; end
    while hi < length(grid) && above[hi+1]; hi += 1; end
    (grid[lo], grid[hi])
end
function powcurve(pk, metric, cv, grid, anchors, ids, nss)
    M = metric_matrix(pk, metric, anchors, ids); rj = findfirst(==(REF_T50), anchors); rv = M[:, rj]
    [BEUtils.vbe_replicate([interp_log(M[i,:], anchors, t) for i in eachindex(ids)], rv;
                           n_subj=nss, n_trial=N_TRIAL_DEFAULT, wsv_cv=cv).pass_rate for t in grid]
end
function powcurve_joint(pk, cvC, cvA, grid, anchors, ids, nss)
    MC = metric_matrix(pk, "Cmax", anchors, ids); MA = metric_matrix(pk, "AUC_tEnd", anchors, ids)
    rj = findfirst(==(REF_T50), anchors); rC = MC[:,rj]; rA = MA[:,rj]
    [BEUtils.vbe_replicate_joint([interp_log(MC[i,:],anchors,t) for i in eachindex(ids)], rC,
                                 [interp_log(MA[i,:],anchors,t) for i in eachindex(ids)], rA;
                                 n_subj=nss, n_trial=N_TRIAL_DEFAULT, cv_C=cvC, cv_A=cvA).pass_rate for t in grid]
end

# ---- run both modes ---------------------------------------------------------
println("Fitting joint model — γ ESTIMATED (shared) ...")
endp_open, ivivc_open = fit_mode(nothing)
println("Fitting joint model — γ FIXED = 1.0 ...")
endp_g1,   ivivc_g1   = fit_mode(1.0)

pbpk = CSV.read(joinpath(ROOT, "data", "pbpk_pk_parameters_individual.csv"), DataFrame)
SOURCES = ["PBPK"=>pbpk, "PopPK_gammaOpen"=>endp_open, "PopPK_gamma1"=>endp_g1]
COL = Dict("PBPK"=>:black, "PopPK_gammaOpen"=>:royalblue, "PopPK_gamma1"=>:firebrick)

anchors = sort(unique(pbpk.T50)); grid = collect(range(anchors[1], anchors[end]; step=0.05))
endpoints = [("Cmax",:single,WSV_CV["Cmax"]), ("AUC_tEnd",:single,WSV_CV["AUC_tEnd"]),
             ("AUC_inf",:single,WSV_CV["AUC_inf"]),
             ("Joint (Cmax & AUC_0-t)",:joint,(WSV_CV["Cmax"],WSV_CV["AUC_tEnd"]))]

bnd = DataFrame()
for (label, kind, cv) in endpoints
    token = replace(label, r"[^A-Za-z0-9]"=>"")
    fig = Figure(size=(1200,520), fontsize=16)
    Label(fig[0,1:2], "γ comparison — $label — PBPK vs PopPK(γ open) vs PopPK(γ=1)"; fontsize=18, font=:bold)
    for (pcol, nss) in enumerate([48, 86])
        ax = Axis(fig[1,pcol]; xlabel="Test T50 (h)", ylabel=pcol==1 ? "power" : "", title="n=$nss", xticks=anchors)
        xlims!(ax, anchors[1], anchors[end]); ylims!(ax, -0.02, 1.02)
        hlines!(ax, [0.05,0.80]; linestyle=:dot, color=:gray50); vlines!(ax, [REF_T50]; linestyle=:dash, color=:gray40)
        for (src, pk) in SOURCES
            idl = sort(unique(pk.IndividualId))
            pw = kind===:joint ? powcurve_joint(pk, cv[1], cv[2], grid, anchors, idl, nss) :
                                  powcurve(pk, label=="AUC_tEnd" ? "AUC_tEnd" : label=="AUC_inf" ? "AUC_inf" : "Cmax", cv, grid, anchors, idl, nss)
            lines!(ax, grid, pw; color=COL[src], linewidth=2.5, label=src)
            lo, hi = safe_bounds(grid, pw, 0.80)
            push!(bnd, (endpoint=label, source=src, n_subj=nss, ss_lo=lo, ss_hi=hi); cols=:union)
        end
        pcol==2 && axislegend(ax; position=:rb, framevisible=false)
    end
    save(joinpath(OUT, "fig_gamma_comparison_$(token).png"), fig, px_per_unit=2)
end

# ---- IVIVC comparison table -------------------------------------------------
ivivc = DataFrame(T50=t50s, alpha_gammaOpen=ivivc_open.gm_alpha, alpha_gamma1=ivivc_g1.gm_alpha)
CSV.write(joinpath(OUT, "gamma_comparison_ivivc.csv"), ivivc)
CSV.write(joinpath(OUT, "gamma_comparison_boundaries.csv"), bnd)

println("\n=== IVIVC: fitted geometric-mean α vs T50 ===")
show(ivivc; allrows=true, allcols=true); println()
println("\n=== Safe-space boundaries (power ≥ 0.80): PBPK vs γ-open vs γ=1 ===")
show(sort(bnd, [:endpoint,:n_subj,:source]); allrows=true, allcols=true); println()
println("\nWrote gamma_comparison_boundaries.csv, gamma_comparison_ivivc.csv, fig_gamma_comparison_*.png -> $OUT")
