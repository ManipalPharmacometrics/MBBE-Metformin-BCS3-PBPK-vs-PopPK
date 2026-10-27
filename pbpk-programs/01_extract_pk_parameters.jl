# =============================================================================
# 01_extract_pk_parameters.jl   (Julia primary pipeline -- Step 1)
#
# Reads the PK-Sim population exports (one *-PK-Analyses.csv per formulation in
# pksim-output/popsims/) and assembles one tidy per-individual table
# of the BE metrics (Cmax, AUC_inf, AUC_tEnd) across the T50 sweep.
#
# Output: output/pk_parameters_individual.csv
# Run from programs/:  julia --project=. 01_extract_pk_parameters.jl
# =============================================================================

using CSV, DataFrames, Statistics
include(joinpath(@__DIR__, "julia", "be_utils.jl"))
using .BEUtils

const POPSIMS = joinpath(@__DIR__, "..", "pksim-output", "popsims")
const OUT     = joinpath(@__DIR__, "output")
mkpath(OUT)

println("Reading population PK analyses from: ", normpath(POPSIMS))
pk = load_all_formulations(POPSIMS)
sort!(pk, [:T50, :IndividualId])

out_csv = joinpath(OUT, "pk_parameters_individual.csv")
CSV.write(out_csv, pk)

println("\nFormulations (T50, h): ", join(sort(unique(pk.T50)), ", "))
println("Subjects per formulation: ", length(unique(pk.IndividualId)))
println("Reference formulation (T50): ", REF_T50, " h")
println("Wrote ", nrow(pk), " rows -> ", out_csv, "\n")

# Geometric-mean summary per formulation (sanity check vs pksim-output/comp.xlsx).
gm(x) = (v = collect(skipmissing(x)); v = v[isfinite.(v) .& (v .> 0)];
         exp(mean(log.(v))))
summ = combine(groupby(pk, :T50),
               nrow => :n,
               :Cmax     => gm => :gm_Cmax,
               :AUC_inf  => gm => :gm_AUC_inf,
               :AUC_tEnd => gm => :gm_AUC_tEnd)
println("Geometric-mean PK metrics by formulation:")
show(summ; allrows = true); println()
