# =============================================================================
# 03_poppk_metrics.jl  — NCA endpoints from the PopPK individual predictions.
#
# Reads the dense IPRED profiles (output/poppk_ipred_profiles.csv) and computes
# the three BE endpoints per subject per formulation, with the SAME column
# schema as the PBPK side (data/pbpk_pk_parameters_individual.csv):
#       T50, IndividualId, Cmax, AUC_inf, AUC_tEnd
# so that 04 can run identical BE/VBE math on both sources.
#
# AUC_tEnd integrates over 0-24 h (the PBPK profile window); AUC_inf adds the
# terminal-slope extrapolation — same convention as the PK-Sim NCA, so the
# flip-flop / extrapolation behaviour is comparable between the two methods.
#
# Input : output/poppk_ipred_profiles.csv   (from 02)
# Output: output/poppk_pk_parameters_individual.csv
# =============================================================================

using CSV, DataFrames, Statistics

const HERE = @__DIR__
const ROOT = normpath(joinpath(HERE, ".."))
const OUT  = joinpath(ROOT, "output")

include(joinpath(ROOT, "src", "nca_lite.jl"))      # auc_trap, lambda_z, nca_metrics

prof = CSV.read(joinpath(OUT, "poppk_ipred_profiles.csv"), DataFrame)
sort!(prof, [:T50, :IndividualId, :time])

rows = DataFrame()
for g in groupby(prof, [:T50, :IndividualId])
    m = nca_metrics(g.time, g.ipred)
    push!(rows, (T50 = g.T50[1], IndividualId = g.IndividualId[1],
                 Cmax = m.Cmax, AUC_inf = m.AUC_inf, AUC_tEnd = m.AUC_tEnd,
                 Tmax = m.Tmax, pct_extrap = m.pct_extrap, lz_r2adj = m.lz_r2adj);
          cols = :union)
end
sort!(rows, [:T50, :IndividualId])

# Write the BE-pipeline columns (matching PBPK schema) plus diagnostics.
be_cols = rows[:, [:T50, :IndividualId, :Cmax, :AUC_inf, :AUC_tEnd]]
CSV.write(joinpath(OUT, "poppk_pk_parameters_individual.csv"), be_cols)
CSV.write(joinpath(OUT, "poppk_nca_diagnostics.csv"), rows)

gm(x) = (v = collect(skipmissing(x)); v = v[isfinite.(v) .& (v .> 0)]; exp(mean(log.(v))))
summ = combine(groupby(rows, :T50),
               nrow => :n,
               :Cmax     => gm => :gm_Cmax,
               :AUC_tEnd => gm => :gm_AUC_tEnd,
               :AUC_inf  => gm => :gm_AUC_inf,
               :pct_extrap => (x -> round(mean(skipmissing(x)), digits=1)) => :mean_pct_extrap)
println("PopPK NCA geometric-mean endpoints by formulation:")
show(summ; allrows = true, allcols = true); println()
println("\nWrote poppk_pk_parameters_individual.csv (BE schema) + poppk_nca_diagnostics.csv -> $OUT")
