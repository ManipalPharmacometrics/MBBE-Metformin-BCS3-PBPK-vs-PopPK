# =============================================================================
# 04_be_safespace_overlay.jl  — PBPK vs PopPK decision-equivalence overlay.
#
# THE HEADLINE DELIVERABLE. Loads the per-subject BE endpoints from BOTH
# methods and runs the *same* TOST / VBE / safe-space machinery (BEUtils, the
# exact code used on the PBPK side) on each, then overlays them:
#   1. GMR(T50) with 90% CI, PBPK vs PopPK, per endpoint           -> fig_overlay_gmr.png
#   2. Safe-space power curves, PBPK vs PopPK, per endpoint, n=48/86 -> fig_overlay_safespace_<ep>.png
#   3. Boundary comparison table                                    -> safe_space_overlay_boundaries.csv
#
# "Decision equivalence" = do the two methods put the safe-space boundaries in
# the same place and reach the same BE conclusions?
#
# Inputs : data/pbpk_pk_parameters_individual.csv     (PBPK, from the PK-Sim pipeline)
#          output/poppk_pk_parameters_individual.csv  (PopPK, from 03)
# =============================================================================

using CSV, DataFrames, Statistics, CairoMakie

const HERE = @__DIR__
const ROOT = normpath(joinpath(HERE, ".."))
const OUT  = joinpath(ROOT, "output");  mkpath(OUT)

include(joinpath(ROOT, "src", "be_utils.jl"));  using .BEUtils

pbpk  = CSV.read(joinpath(ROOT, "data", "pbpk_pk_parameters_individual.csv"), DataFrame)
poppk = CSV.read(joinpath(OUT, "poppk_pk_parameters_individual.csv"), DataFrame)
SOURCES = ["PBPK" => pbpk, "PopPK" => poppk]
COL = Dict("PBPK" => :firebrick, "PopPK" => :royalblue)

# ---- 1. GMR(T50) + 90% CI table for both sources ----------------------------
gmr_rows = DataFrame()
for (src, pk) in SOURCES
    ref = sort(filter(:T50 => ==(REF_T50), pk), :IndividualId)
    for t50 in sort(unique(pk.T50))
        tst = sort(filter(:T50 => ==(t50), pk), :IndividualId)
        for (m, _) in BE_METRICS
            r = BEUtils.be_paired(tst[!, m], ref[!, m])
            push!(gmr_rows, (source = src, T50 = t50, metric = m, n = r.n,
                             GMR = r.GMR, CI_lower = r.CI_lower, CI_upper = r.CI_upper,
                             BE = r.BE); cols = :union)
        end
    end
end
CSV.write(joinpath(OUT, "be_overlay_gmr.csv"), gmr_rows)

# ---- GMR overlay figure (one column per endpoint) ---------------------------
let
    eps = first.(BE_METRICS)
    fig = Figure(size = (1500, 480), fontsize = 16)
    Label(fig[0, 1:length(eps)],
          "GMR vs T50 — PBPK vs PopPK  (ref T50 = $(REF_T50) h, 90% CI)";
          fontsize = 19, font = :bold)
    for (col, m) in enumerate(eps)
        ax = Axis(fig[1, col]; title = m, xlabel = "Test T50 (h)",
                  ylabel = col == 1 ? "GMR (test/ref)" : "",
                  xticks = sort(unique(pbpk.T50)))
        hspan!(ax, BE_LOWER, BE_UPPER; color = (:seagreen, 0.10))
        hlines!(ax, [BE_LOWER, BE_UPPER]; color = :gray60, linestyle = :dash)
        hlines!(ax, [1.0]; color = :gray60, linestyle = :dot)
        for (src, _) in SOURCES
            sub = sort(filter(r -> r.source == src && r.metric == m, gmr_rows), :T50)
            off = src == "PBPK" ? -0.04 : 0.04          # dodge so CIs don't overlap
            xs  = sub.T50 .+ off
            errorbars!(ax, xs, sub.GMR, sub.GMR .- sub.CI_lower, sub.CI_upper .- sub.GMR;
                       color = COL[src], whiskerwidth = 8, linewidth = 1.5)
            scatterlines!(ax, xs, sub.GMR; color = COL[src], markersize = 9,
                          label = src, linewidth = 1.5)
        end
        col == length(eps) && axislegend(ax; position = :rt, framevisible = false)
    end
    save(joinpath(OUT, "fig_overlay_gmr.png"), fig, px_per_unit = 2)
end

# ---- 2. Safe-space power curves, PBPK vs PopPK ------------------------------
const GRID_DX = 0.05
function metric_matrix(pk, metric, anchors, ids)
    M = Matrix{Union{Missing,Float64}}(missing, length(ids), length(anchors))
    for (j, t) in enumerate(anchors)
        d = sort(filter(:T50 => ==(t), pk), :IndividualId)
        for (i, id) in enumerate(ids)
            r = findfirst(==(id), d.IndividualId)
            r !== nothing && (M[i, j] = d[r, metric])
        end
    end
    M
end
function interp_log(vals, anchors, t)
    ok = .!ismissing.(vals); sum(ok) < 2 && return missing
    x = anchors[ok]; v = Float64.(vals[ok])
    t <= x[1]   && return v[1]
    t >= x[end] && return v[end]
    k = searchsortedlast(x, t); (t == x[k]) && return v[k]
    w = (log(t) - log(x[k])) / (log(x[k+1]) - log(x[k]))
    exp((1 - w) * log(v[k]) + w * log(v[k+1]))
end
safe_bounds(grid, power, ptgt) = begin
    above = power .>= ptgt
    ri = argmin(abs.(grid .- REF_T50)); lo = ri; hi = ri
    while lo > 1 && above[lo-1]; lo -= 1; end
    while hi < length(grid) && above[hi+1]; hi += 1; end
    (grid[lo], grid[hi])
end

function power_curve(pk, metric, cv, grid, anchors, ids, nss; n_trial = N_TRIAL_DEFAULT)
    M  = metric_matrix(pk, metric, anchors, ids)
    rj = findfirst(==(REF_T50), anchors); rv = M[:, rj]
    [BEUtils.vbe_replicate([interp_log(M[i, :], anchors, t) for i in eachindex(ids)], rv;
                           n_subj = nss, n_trial = n_trial, wsv_cv = cv).pass_rate for t in grid]
end
function power_curve_joint(pk, cvC, cvA, grid, anchors, ids, nss; n_trial = N_TRIAL_DEFAULT)
    MC = metric_matrix(pk, "Cmax", anchors, ids); MA = metric_matrix(pk, "AUC_tEnd", anchors, ids)
    rj = findfirst(==(REF_T50), anchors); rC = MC[:, rj]; rA = MA[:, rj]
    [BEUtils.vbe_replicate_joint([interp_log(MC[i, :], anchors, t) for i in eachindex(ids)], rC,
                                 [interp_log(MA[i, :], anchors, t) for i in eachindex(ids)], rA;
                                 n_subj = nss, n_trial = n_trial, cv_C = cvC, cv_A = cvA).pass_rate
     for t in grid]
end

const PTGT = 0.80
bnd_rows = DataFrame()
# endpoint => (kind, cv-or-cvs)
endpoints = [("Cmax", :single, WSV_CV["Cmax"]),
             ("AUC_tEnd", :single, WSV_CV["AUC_tEnd"]),
             ("AUC_inf", :single, WSV_CV["AUC_inf"]),
             ("Joint (Cmax & AUC_0-t)", :joint, (WSV_CV["Cmax"], WSV_CV["AUC_tEnd"]))]

for (label, kind, cv) in endpoints
    anchors = sort(unique(pbpk.T50))                    # same 8 anchors both sources
    grid    = collect(range(anchors[1], anchors[end]; step = GRID_DX))
    token   = replace(label, r"[^A-Za-z0-9]" => "")
    fig = Figure(size = (1200, 520), fontsize = 16)
    Label(fig[0, 1:2], "Safe space — $label — PBPK vs PopPK  (ref $(REF_T50) h, $(N_TRIAL_DEFAULT) trials)";
          fontsize = 18, font = :bold)
    for (pcol, nss) in enumerate(N_SUBJ_LIST)
        ax = Axis(fig[1, pcol]; xlabel = "Test T50 (h)",
                  ylabel = pcol == 1 ? "P(conclude BE) — power" : "",
                  title = "n = $nss", xticks = anchors)
        xlims!(ax, anchors[1], anchors[end]); ylims!(ax, -0.02, 1.02)
        hlines!(ax, [0.05, PTGT]; linestyle = :dot, color = :gray50)
        vlines!(ax, [REF_T50]; linestyle = :dash, color = :gray40)
        for (src, pk) in SOURCES
            ids = sort(unique(pk.IndividualId))
            pw  = kind === :joint ?
                  power_curve_joint(pk, cv[1], cv[2], grid, anchors, ids, nss) :
                  power_curve(pk, label == "AUC_tEnd" ? "AUC_tEnd" :
                                  label == "AUC_inf"  ? "AUC_inf"  : "Cmax",
                              cv, grid, anchors, ids, nss)
            lines!(ax, grid, pw; color = COL[src], linewidth = 2.5, label = src)
            lo, hi = safe_bounds(grid, pw, PTGT)
            vlines!(ax, [lo, hi]; color = COL[src], linestyle = :dot, linewidth = 1.5)
            push!(bnd_rows, (endpoint = label, source = src, n_subj = nss,
                             ss_lo = lo, ss_hi = hi); cols = :union)
        end
        pcol == 2 && axislegend(ax; position = :rb, framevisible = false)
    end
    save(joinpath(OUT, "fig_overlay_safespace_$(token).png"), fig, px_per_unit = 2)
end

CSV.write(joinpath(OUT, "safe_space_overlay_boundaries.csv"), bnd_rows)

println("\nSafe-space boundaries (power >= $(PTGT)) — PBPK vs PopPK:")
show(sort(bnd_rows, [:endpoint, :n_subj, :source]); allrows = true, allcols = true); println()
println("\nWrote: be_overlay_gmr.csv, safe_space_overlay_boundaries.csv,")
println("       fig_overlay_gmr.png, fig_overlay_safespace_*.png -> $OUT")
