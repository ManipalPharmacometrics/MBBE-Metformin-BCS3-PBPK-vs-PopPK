# =============================================================================
# 02_be_safe_space.jl   (Julia primary pipeline -- Steps 2-3)
#
# BE assessment + VBE replicate-trial pass-rate / error rate.
# Reference formulation: T50 = 3.0 h. Every other T50 is a TEST against it.
#
# Outputs under output/:
#   be_point_estimates.csv      paired BE (GMR, 90% CI, pass) per endpoint  (CI-based)
#   be_table_publication.csv    publication wide BE table (one row per formulation)
#   vbe_trial_results.csv        VBE pass-rate (power / Type I error), n=48 & 86, 500 trials
#   safe_space_summary.csv       per-formulation overall BE verdict
#   fig_safe_space_gmr.png       GMR +/- 90% CI vs T50, with 0.80-1.25 band
#
# The Simcyp-style safe-space POWER curves (per endpoint) are produced by
# 03_safe_space_plot.jl.
#
# VBE settings (as specified): 500 trials; study sizes n = 48 and 86.
# Within-subject CV is endpoint-specific and sourced (BEUtils.WSV_CV):
#   metformin ER, FDA NDA 212595 -> Cmax 12.2%, AUC0-t 11.2%, AUC_inf 11.4%.
#
# Run from programs/:  julia --project=. 02_be_safe_space.jl
# =============================================================================

using CSV, DataFrames, Statistics, Printf, CairoMakie
include(joinpath(@__DIR__, "julia", "be_utils.jl"))
using .BEUtils

const OUT = joinpath(@__DIR__, "output")
mkpath(OUT)

pk = CSV.read(joinpath(OUT, "pk_parameters_individual.csv"), DataFrame)
nrow(pk) == 0 && error("Run 01_extract_pk_parameters.jl first.")

ref_df   = filter(:T50 => ==(REF_T50), pk)
nrow(ref_df) == 0 && error("Reference T50 = $REF_T50 not in data.")
test_t50 = sort(setdiff(unique(pk.T50), [REF_T50]))

const N_TRIAL = N_TRIAL_DEFAULT   # 500

# --- 1. Point (single-trial) paired BE, every test formulation ---------------
# Deterministic CI-based BE on the as-simulated n=24 population.
point = DataFrame()
for tt in test_t50
    res = be_all_metrics(filter(:T50 => ==(tt), pk), ref_df)
    res[!, :T50]     .= tt
    res[!, :ref_T50] .= REF_T50
    append!(point, res)
end
select!(point, [:T50, :ref_T50, :metric, :n, :GMR, :CI_lower, :CI_upper, :BE])
CSV.write(joinpath(OUT, "be_point_estimates.csv"), point)

# --- 1b. Publication wide BE table (one row per formulation) -----------------
fmt_gmr(x) = @sprintf("%.3f", x)
fmt_ci(lo, hi) = @sprintf("%.2f–%.2f", lo, hi)
getrow(tt, m) = only(filter(r -> r.T50 == tt && r.metric == m, point))
pub = DataFrame()
for tt in test_t50
    c = getrow(tt, "Cmax"); a = getrow(tt, "AUC_tEnd"); ai = getrow(tt, "AUC_inf")
    push!(pub, (
        var"Test T50 (h)"        = tt,
        var"Cmax GMR"            = fmt_gmr(c.GMR),
        var"Cmax 90% CI (%)"     = fmt_ci(100c.CI_lower, 100c.CI_upper),
        var"Cmax"                = c.BE ? "Pass" : "Fail",
        var"AUC0-t GMR"          = fmt_gmr(a.GMR),
        var"AUC0-t 90% CI (%)"   = fmt_ci(100a.CI_lower, 100a.CI_upper),
        var"AUC0-t"              = a.BE ? "Pass" : "Fail",
        var"AUCinf GMR"          = fmt_gmr(ai.GMR),
        var"AUCinf 90% CI (%)"   = fmt_ci(100ai.CI_lower, 100ai.CI_upper),
        var"AUCinf"              = ai.BE ? "Pass" : "Fail",
        var"Overall BE"          = (c.BE && a.BE) ? "Pass" : "Fail",
    ))
end
CSV.write(joinpath(OUT, "be_table_publication.csv"), pub)

# --- 2. VBE replicate trials: pass-rate (power / Type I error) ---------------
# Endpoint-specific sourced CV; study sizes n = 48 and 86; 500 trials each.
vbe = DataFrame()
for nss in N_SUBJ_LIST, tt in test_t50
    tdf = sort(filter(:T50 => ==(tt), pk), :IndividualId)
    rdf = sort(ref_df, :IndividualId)
    for m in first.(BE_METRICS)
        cv = WSV_CV[m]
        r = vbe_replicate(tdf[!, m], rdf[!, m]; n_subj = nss, n_trial = N_TRIAL, wsv_cv = cv)
        tv = tdf[!, m]; rv = rdf[!, m]
        ok = @. !ismissing(tv) & !ismissing(rv) & isfinite(tv) & isfinite(rv) & (tv > 0) & (rv > 0)
        true_gmr = exp(mean(log.(tv[ok])) - mean(log.(rv[ok])))
        truly_be = BE_LOWER <= true_gmr <= BE_UPPER
        push!(vbe, (T50 = tt, metric = m, n_subj = nss, wsv_cv = cv,
                    true_GMR = true_gmr, truly_BE = truly_be,
                    pass_rate = r.pass_rate,
                    risk_type = truly_be ? "power" : "type_I_error"))
    end
end
CSV.write(joinpath(OUT, "vbe_trial_results.csv"), vbe)

# --- 3. Safe-space summary (overall BE = Cmax AND AUC_0-t both pass) ----------
summ = DataFrame()
for tt in test_t50
    c = getrow(tt, "Cmax"); a = getrow(tt, "AUC_tEnd")
    push!(summ, (T50 = tt,
                 Cmax_GMR = c.GMR,
                 Cmax_CI  = @sprintf("[%.3f, %.3f]", c.CI_lower, c.CI_upper),
                 AUCt_GMR = a.GMR,
                 AUCt_CI  = @sprintf("[%.3f, %.3f]", a.CI_lower, a.CI_upper),
                 overall_BE = c.BE && a.BE))
end
CSV.write(joinpath(OUT, "safe_space_summary.csv"), summ)

# --- 4. GMR +/- 90% CI forest plot vs T50, with acceptance band --------------
metric_cols = ["Cmax" => :firebrick, "AUC_inf" => :steelblue, "AUC_tEnd" => :seagreen]
all_t50 = sort(unique(pk.T50))
fig = Figure(size = (1000, 700), fontsize = 18)
ax  = Axis(fig[1, 1], xlabel = "Test formulation T50 (h)",
           ylabel = "GMR vs reference (T50 = 3.0 h)",
           title = "Metformin XR BE point estimates (paired PBPK crossover, n=24)",
           xlabelsize = 19, ylabelsize = 19, titlesize = 20,
           xticklabelsize = 16, yticklabelsize = 16,
           xticks = all_t50)   # a tick at every simulated T50
band!(ax, [minimum(all_t50)-0.3, maximum(all_t50)+0.3], [BE_LOWER, BE_LOWER],
      [BE_UPPER, BE_UPPER], color = (:seagreen, 0.12))
hlines!(ax, [BE_LOWER, BE_UPPER]; linestyle = :dash, color = :gray40)
hlines!(ax, [1.0]; linestyle = :solid, color = :gray60)
vlines!(ax, [REF_T50]; linestyle = :dot, color = :gray50)
for (m, cl) in metric_cols
    d = sort(filter(:metric => ==(m), point), :T50)
    rangebars!(ax, d.T50, d.CI_lower, d.CI_upper, color = cl, whiskerwidth = 8)
    scatterlines!(ax, d.T50, d.GMR, color = cl, label = m, markersize = 9)
end
axislegend(ax, position = :rt)
save(joinpath(OUT, "fig_safe_space_gmr.png"), fig, px_per_unit = 2)

# --- 5. Console report -------------------------------------------------------
println("\n=== Publication BE table (vs T50 = 3.0 h) ===")
show(pub; allrows = true, allcols = true, truncate = 40); println()
println("\n=== VBE pass-rate / error rate (500 trials; sourced CV) ===")
show(select(vbe, [:T50, :metric, :n_subj, :true_GMR, :pass_rate, :risk_type]);
     allrows = true, allcols = true); println()
println("\nWrote CSVs and GMR figure to ", normpath(OUT))
