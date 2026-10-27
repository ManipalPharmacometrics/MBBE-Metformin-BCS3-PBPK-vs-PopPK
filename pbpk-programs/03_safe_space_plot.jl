# =============================================================================
# 03_safe_space_plot.jl   (Julia primary pipeline -- Step 3b, ER scenario)
#
# Simcyp-style SAFE-SPACE plots per endpoint (Cmax, AUC_inf, AUC_0-t, joint),
# each with study-size panels n = 48 and n = 86, 500 trials. Reference T50=3.0h.
# Within-subject CV = metformin EXTENDED-release (BEUtils.WSV_CV; FDA NDA 212595).
#
# Shared logic lives in julia/safespace.jl. Outputs to output/:
#   fig_safespace_Cmax.png, _AUCinf.png, _AUC0t.png, _joint.png
#   safe_space_boundaries.csv
#
# Run from programs/:  julia --project=. 03_safe_space_plot.jl
# =============================================================================

using CSV, DataFrames, Statistics, CairoMakie
include(joinpath(@__DIR__, "julia", "be_utils.jl"))
using .BEUtils
include(joinpath(@__DIR__, "julia", "safespace.jl"))

const OUT = joinpath(@__DIR__, "output")
pk = CSV.read(joinpath(OUT, "pk_parameters_individual.csv"), DataFrame)

bnd = build_safespace(pk, OUT; cv_map = WSV_CV, tag = "", scenario = "ER CV")
CSV.write(joinpath(OUT, "safe_space_boundaries.csv"), bnd)

println("ER safe-space boundaries (power ≥ 0.80):")
show(bnd; allrows = true, allcols = true); println()
println("Wrote 4 ER endpoint figures + safe_space_boundaries.csv to ", normpath(OUT))
