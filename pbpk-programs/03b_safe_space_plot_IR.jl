# =============================================================================
# 03b_safe_space_plot_IR.jl   (Julia primary pipeline -- Step 3b, IR scenario)
#
# Same Simcyp-style safe-space plots as 03, but using the HIGH-variability
# metformin IMMEDIATE-release within-subject CV (BEUtils.WSV_CV_IR:
#   Cmax 23.7%, AUC 14.7%) — a sensitivity bound showing how the safe space
# narrows under higher intra-subject variability.
#
# Outputs to output/:
#   fig_safespace_IR_Cmax.png, _IR_AUCinf.png, _IR_AUC0t.png, _IR_joint.png
#   safe_space_boundaries_IR.csv
#
# Run from programs/:  julia --project=. 03b_safe_space_plot_IR.jl
# =============================================================================

using CSV, DataFrames, Statistics, CairoMakie
include(joinpath(@__DIR__, "julia", "be_utils.jl"))
using .BEUtils
include(joinpath(@__DIR__, "julia", "safespace.jl"))

const OUT = joinpath(@__DIR__, "output")
pk = CSV.read(joinpath(OUT, "pk_parameters_individual.csv"), DataFrame)

bnd = build_safespace(pk, OUT; cv_map = WSV_CV_IR, tag = "IR", scenario = "IR CV")
CSV.write(joinpath(OUT, "safe_space_boundaries_IR.csv"), bnd)

println("IR (high-variability) safe-space boundaries (power ≥ 0.80):")
show(bnd; allrows = true, allcols = true); println()
println("Wrote 4 IR endpoint figures + safe_space_boundaries_IR.csv to ", normpath(OUT))
