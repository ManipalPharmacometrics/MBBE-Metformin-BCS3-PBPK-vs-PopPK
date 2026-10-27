# =============================================================================
# 02b_fit_weibull_nopumas.jl  — Pumas-FREE substitute for step 02.
#
# Fits the SAME Weibull-release 1-cmt absorption model as src/poppk_model.jl,
# but per subject (two-stage) with open-Julia nonlinear least squares (Optim),
# so the entire PBPK->PopPK safe-space pipeline can be run and cross-checked
# WITHOUT a Pumas license. Writes exactly the files step 02 writes, so steps 03
# and 04 are byte-for-byte identical regardless of which fit produced them.
#
# Inputs : output/sparse_be_dataset.csv        (from 01)
#          src/weibull_absorption.jl            (simulator + per-subject fit)
# Outputs: output/poppk_param_estimates.csv     population summary per arm
#          output/poppk_icoef.csv               per-subject Alpha/CL/Vc
#          output/poppk_ipred_profiles.csv       dense individual predictions
#
# Also prints the KEY validation for this project: does the fitted in-vivo
# Weibull release scale (geometric-mean Alpha) recover the PBPK input T50?
#
#   julia --project=<open env> scripts/02b_fit_weibull_nopumas.jl
# =============================================================================

using CSV, DataFrames, Statistics, Optim

const HERE = @__DIR__
const ROOT = normpath(joinpath(HERE, ".."))
const OUT  = get(ENV, "POPPK_OUT", joinpath(ROOT, "output"));  mkpath(OUT)

include(joinpath(ROOT, "src", "weibull_absorption.jl"))   # simulate_weibull, fit_subject_weibull

const DENSE = collect(0.0:0.05:24.0)      # match the PBPK profile window (0-24 h)
const LAG   = 0.25                        # fixed release lag (as in the model)

data = CSV.read(joinpath(OUT, "sparse_be_dataset.csv"), DataFrame)
dose = maximum(skipmissing(data.amt))     # 3019.5 µmol (dose rows, evid==1)
t50s = sort(unique(data.T50))

param_rows = DataFrame()
icoef_rows = DataFrame()
ipred_rows = DataFrame()

for t50 in t50s
    arm = filter(:T50 => ==(t50), data)
    ids = sort(unique(arm.IndividualId))
    @info "Fitting (no-Pumas) T50 = $t50 h  (n = $(length(ids)) subjects)"
    as = Float64[]; gs = Float64[]; cls = Float64[]; vcs = Float64[]
    for id in ids
        obs = filter(r -> r.IndividualId == id && r.evid == 0, arm)
        sort!(obs, :time)
        fit = fit_subject_weibull(obs.time, collect(skipmissing(obs.dv)), dose;
                                  lag = LAG, alpha0 = t50)
        push!(as, fit.alpha); push!(gs, fit.gamma); push!(cls, fit.CL); push!(vcs, fit.Vc)
        push!(icoef_rows, (T50 = t50, IndividualId = id, SUBJID = "T50_$(t50)_ID_$(id)",
                           Alpha = fit.alpha, CL = fit.CL, Vc = fit.Vc,
                           gamma = fit.gamma, lag = LAG); cols = :union)
        # dense individual prediction from the fitted subject parameters
        pred = simulate_weibull(DENSE, fit.alpha, fit.gamma, fit.CL, fit.Vc, LAG, dose)
        for (tk, ck) in zip(DENSE, pred)
            push!(ipred_rows, (T50 = t50, IndividualId = id, time = tk, ipred = ck); cols = :union)
        end
    end
    gmean(x) = exp(mean(log.(x)))
    logvar(x) = var(log.(x))
    push!(param_rows, (T50 = t50, tvalpha = gmean(as), tvgamma = gmean(gs),
                       tvcl = gmean(cls), tvvc = gmean(vcs), tvlag = LAG,
                       omega_alpha = logvar(as), omega_cl = logvar(cls),
                       omega_vc = logvar(vcs), sigma_p = missing, OFV = missing);
          cols = :union)
end

CSV.write(joinpath(OUT, "poppk_param_estimates.csv"), param_rows)
CSV.write(joinpath(OUT, "poppk_icoef.csv"),           icoef_rows)
CSV.write(joinpath(OUT, "poppk_ipred_profiles.csv"),  ipred_rows)

println("\nPopulation parameter estimates per formulation (no-Pumas two-stage):")
show(param_rows; allrows = true, allcols = true); println()

# ---- KEY VALIDATION: does fitted Alpha recover the PBPK input T50? ----------
rec = select(param_rows, :T50, :tvalpha => :fitted_alpha)
rec.ratio_alpha_over_T50 = rec.fitted_alpha ./ rec.T50
println("\nIVIVC check — fitted in-vivo Weibull scale (Alpha) vs PBPK input T50:")
show(rec; allrows = true, allcols = true); println()
println("\nWrote poppk_param_estimates.csv, poppk_icoef.csv, poppk_ipred_profiles.csv -> $OUT")
