# =============================================================================
# STEP 7 — IVIVC check: does the PopPK release follow the PK-Sim release?
#
# In : results/perarm/step2_parameters.csv, results/joint/step2_parameters.csv,
#      results/joint/step2_standard_errors.csv
# Out: results/step7_ivivc.csv       per formulation: PK-Sim T50 vs PopPK release T50
#      results/step7_ivivc_fit.csv   log-log regression and the joint-model slope
#      results/fig_step7_ivivc.png
#
# How PopPK relates to PBPK here:
#   PK-Sim T50 = the in vivo Weibull release time we set in each simulation
#   (the "truth"). The PopPK models estimate their own Weibull release, so the
#   comparable number is the PopPK time to 50% released:
#       T50_PopPK = Alpha * ln(2)^(1/Gamma)       (from 1 - exp(-(t/Alpha)^Gamma) = 0.5)
#   per-arm: one T50_PopPK per formulation (8 fits).
#   joint  : T50_PopPK = tvalpha * (T50/3)^slope * ln(2)^(1/Gamma). An IVIVC
#            would give slope ≈ 1 (PopPK release scales with PK-Sim T50);
#            slope ≈ 0 means the PopPK release does not track the formulation.
# A Level A style check: regress log(T50_PopPK) on log(T50_PK-Sim) for the
# per-arm fits; slope 1 and intercept 0 = PopPK recovers the release.
# =============================================================================

using CSV, DataFrames, Statistics, Distributions, CairoMakie
include("settings.jl")

perarm = CSV.read(joinpath(RESULTS, "perarm", "step2_parameters.csv"), DataFrame)
joint  = only(eachrow(CSV.read(joinpath(RESULTS, "joint", "step2_parameters.csv"), DataFrame)))
se     = CSV.read(joinpath(RESULTS, "joint", "step2_standard_errors.csv"), DataFrame)
slope_row = only(eachrow(se[se.parameter .== "slope", :]))

t50_popPK(alpha, gamma) = alpha * log(2)^(1 / gamma)

ivivc = DataFrame(T50_PKSim = parse.(Float64, match.(r"T50 ([0-9.]+) h", perarm.fit) .|> m -> m.captures[1]))
ivivc.T50_PopPK_perarm = t50_popPK.(perarm.tvalpha, perarm.tvgamma)
ivivc.T50_PopPK_joint  = [t50_popPK(joint.tvalpha * (t / REF_T50)^joint.slope, joint.tvgamma) for t in ivivc.T50_PKSim]
ivivc.CL_perarm = perarm.tvcl                    # where the formulation effect ends up instead
CSV.write(joinpath(RESULTS, "step7_ivivc.csv"), ivivc)

# log-log regression, per-arm (ordinary least squares with a 95% CI on the slope)
x, y = log.(ivivc.T50_PKSim), log.(ivivc.T50_PopPK_perarm)
b = cov(x, y) / var(x); a = mean(y) - b * mean(x)
res = y .- (a .+ b .* x); n = length(x)
se_b = sqrt(sum(res .^ 2) / (n - 2) / sum((x .- mean(x)) .^ 2))
tq = quantile(TDist(n - 2), 0.975)
fitrows = DataFrame(
    analysis = ["per-arm: log-log regression slope", "joint model: estimated slope"],
    estimate = [b, joint.slope],
    ci_lower = [b - tq * se_b, slope_row.ci_lower],
    ci_upper = [b + tq * se_b, slope_row.ci_upper],
    ivivc_would_give = [1.0, 1.0])
CSV.write(joinpath(RESULTS, "step7_ivivc_fit.csv"), fitrows)

fig = Figure(size = (650, 520))
ax = Axis(fig[1, 1]; xlabel = "PK-Sim T50 (h, set in the simulation)", ylabel = "PopPK T50 (h, estimated)",
          title = "Metformin: PopPK release vs PK-Sim release", xticks = ivivc.T50_PKSim)
lines!(ax, [0.8, 7.2], [0.8, 7.2]; color = :gray60, linestyle = :dash, label = "identity (perfect IVIVC)")
scatterlines!(ax, ivivc.T50_PKSim, ivivc.T50_PopPK_perarm; color = :royalblue, label = "per-arm fits")
lines!(ax, ivivc.T50_PKSim, ivivc.T50_PopPK_joint; color = :darkorange, label = "joint model")
axislegend(ax; position = :lt)
save(joinpath(RESULTS, "fig_step7_ivivc.png"), fig)

println("PopPK vs PK-Sim release time:"); show(ivivc; allrows = true); println()
println("\nSlope (1 = PopPK release follows PK-Sim):"); show(fitrows; allcols = true); println()
