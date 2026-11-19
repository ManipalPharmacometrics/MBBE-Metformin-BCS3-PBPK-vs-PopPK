# =============================================================================
# STEP 2 — fit the two PopPK models (needs Pumas).
#
# In : results/step1_dataset.csv, model_perarm.jl, model_joint.jl
# Out: results/<model>/step2_parameters.csv       population estimates and -2 log-likelihood
#      results/<model>/step2_standard_errors.csv  estimates with SE and 95% CI
#      results/<model>/step2_ipred.csv            individual predictions, 0–24 h
#                                                  every 0.05 h (used in step 3)
#   <model> = perarm | joint
#
# perarm: model_perarm.jl fitted 8 times, once per formulation (24 subjects).
# joint : model_joint.jl fitted once to all 8 formulations (192 subject-arms).
# Both: FOCE, default initial values from the model, lag fixed at its initial
# value (1e-3 h).
# =============================================================================

using Pumas, CSV, DataFrames
include("settings.jl")
include("model_perarm.jl")
include("model_joint.jl")

const PRED_TIMES = collect(0.0:0.05:24.0)

data = CSV.read(joinpath(RESULTS, "step1_dataset.csv"), DataFrame)

to_population(df) = read_pumas(df; id = :SUBJID, time = :time, amt = :amt, evid = :evid,
                               cmt = :cmt, observations = [:dv],
                               covariates = [:T50, :IndividualId])

"Fit `model` to `df`; return one row of estimates, the SE table and the predictions."
function fit_and_collect(model, df, label)
    pop = to_population(df)
    fitted = fit(model, pop, init_params(model), FOCE(); constantcoef = (:tvlag,))

    est = coef(fitted)
    params = DataFrame(fit = label, minus2LL = -2 * loglikelihood(fitted),
                       tvalpha = est.tvalpha, tvgamma = est.tvgamma, tvcl = est.tvcl,
                       tvvc = est.tvvc, tvlag = est.tvlag, omega_alpha = est.Ω[1, 1],
                       omega_cl = est.Ω[2, 2], omega_vc = est.Ω[3, 3], sigma_p = est.σ_p)
    haskey(est, :slope) && insertcols!(params, :tvalpha, :slope => est.slope; after = true)

    se = DataFrame(coeftable(infer(fitted)))
    insertcols!(se, 1, :fit => label)

    ipred = DataFrame()
    for (subj, pred) in zip(pop, predict(fitted; obstimes = PRED_TIMES))
        append!(ipred, DataFrame(T50 = subj.covariates(0.0).T50,
                                 IndividualId = Int(subj.covariates(0.0).IndividualId),
                                 time = pred.time, ipred = pred.ipred.dv))
    end
    params, se, ipred
end

for model_name in MODELS
    out = joinpath(RESULTS, model_name); mkpath(out)
    params, se, ipred = DataFrame(), DataFrame(), DataFrame()
    if model_name == "perarm"
        for t50 in sort(unique(data.T50))
            @info "per-arm fit, T50 = $t50 h"
            p, s, i = fit_and_collect(poppk_perarm, data[data.T50 .== t50, :], "T50 $t50 h")
            append!(params, p; cols = :union); append!(se, s); append!(ipred, i)
        end
    else
        @info "joint fit, all formulations"
        params, se, ipred = fit_and_collect(poppk_joint, data, "joint")
    end
    CSV.write(joinpath(out, "step2_parameters.csv"), params)
    CSV.write(joinpath(out, "step2_standard_errors.csv"), se)
    CSV.write(joinpath(out, "step2_ipred.csv"), ipred)
    println("\n$model_name — population estimates:"); show(params; allcols = true); println()
end
