# =============================================================================
# 02b_fit_weibull_nopumas.jl  — Pumas-FREE substitute for step 02 (JOINT model).
#
# Open-Julia analog of the baseline joint Pumas fit: for each physiology
# (IndividualId) the 8 formulation arms are fit TOGETHER with shared,
# FORMULATION-INDEPENDENT disposition (CL, Vc, gamma) and one release scale
# `alpha` per arm. This keeps between-subject variability while forbidding the
# formulation effect from leaking into disposition — so `alpha` tracks the PBPK
# T50 (the IVIVC). Writes the same files as step 02, so 03/04 are identical.
#
# Inputs : output/sparse_be_dataset.csv       (from 01)
#          src/weibull_absorption.jl
# Outputs: output/poppk_param_estimates.csv    per-arm summary (gm across physiologies)
#          output/poppk_icoef.csv              per-subject Alpha/CL/Vc
#          output/poppk_ipred_profiles.csv     dense individual predictions
#
#   julia --project=nopumas_env scripts/02b_fit_weibull_nopumas.jl
# =============================================================================

using CSV, DataFrames, Statistics, Optim

const HERE = @__DIR__
const ROOT = normpath(joinpath(HERE, ".."))
const OUT  = joinpath(ROOT, "output");  mkpath(OUT)

include(joinpath(ROOT, "src", "weibull_absorption.jl"))

const DENSE = collect(0.0:0.05:24.0)
const LAG   = 0.25

data = CSV.read(joinpath(OUT, "sparse_be_dataset.csv"), DataFrame)
dose = maximum(skipmissing(data.amt))
ids  = sort(unique(data.IndividualId))
t50s = sort(unique(data.T50))

icoef_rows = DataFrame(); ipred_rows = DataFrame()

for id in ids
    @info "Joint fit (no-Pumas) physiology IndividualId = $id  (shared CL/Vc/γ across $(length(t50s)) arms)"
    arms = NamedTuple[]
    for t50 in t50s
        obs = sort(filter(r -> r.IndividualId == id && r.T50 == t50 && r.evid == 0, data), :time)
        push!(arms, (t50 = t50, tobs = obs.time, cobs = collect(skipmissing(obs.dv))))
    end
    f = fit_physiology_joint(arms, dose; lag = LAG)
    for t50 in t50s
        a = f.alpha[t50]
        push!(icoef_rows, (T50 = t50, IndividualId = id, SUBJID = "T50_$(t50)_ID_$(id)",
                           Alpha = a, CL = f.CL, Vc = f.Vc, gamma = f.gamma, lag = LAG); cols = :union)
        pred = simulate_weibull(DENSE, a, f.gamma, f.CL, f.Vc, LAG, dose)
        for (tk, ck) in zip(DENSE, pred)
            push!(ipred_rows, (T50 = t50, IndividualId = id, time = tk, ipred = ck); cols = :union)
        end
    end
end

gm(x) = exp(mean(log.(x)))
param_rows = combine(groupby(icoef_rows, :T50),
                     :Alpha => gm => :gm_alpha, :CL => gm => :gm_cl,
                     :Vc => gm => :gm_vc, :gamma => first => :gamma)
param_rows.ratio_alpha_over_T50 = param_rows.gm_alpha ./ param_rows.T50

CSV.write(joinpath(OUT, "poppk_param_estimates.csv"), param_rows)
CSV.write(joinpath(OUT, "poppk_icoef.csv"),           icoef_rows)
CSV.write(joinpath(OUT, "poppk_ipred_profiles.csv"),  ipred_rows)

println("\nJoint no-Pumas fit — per-arm geometric means (disposition shared per physiology):")
show(sort(param_rows, :T50); allrows = true, allcols = true); println()
println("\n(ratio_alpha_over_T50 ≈ constant ⇒ Alpha tracks T50 = IVIVC restored.)")
println("Wrote poppk_param_estimates.csv, poppk_icoef.csv, poppk_ipred_profiles.csv -> $OUT")
