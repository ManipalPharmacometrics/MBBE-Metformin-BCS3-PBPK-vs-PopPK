# =============================================================================
# 02_fit_poppk.jl  — fit the baseline PopPK model to each formulation arm.
#
# Each of the 8 formulations (T50 = 1..7 h) is an independent BE arm of 24
# paired subjects. We fit the same 1-compartment WEIBULL-RELEASE absorption
# model to each arm separately (robust, easy to reason about, and mirrors "each
# formulation characterised independently"). Because the Weibull release IS the
# absorption input, the fitted release scale `Alpha` is the in-vivo analog of
# that arm's PBPK T50 (a direct, checkable IVIVC) and the fitted `Gamma` the
# analog of beta. A joint model with T50 as a covariate on Alpha is a reasonable
# alternative — see README.
#
# NB this step needs a licensed Pumas. For a Pumas-FREE run of the whole
# pipeline, use scripts/02b_fit_weibull_nopumas.jl instead (open-Julia NLS fit);
# it writes the same output files so 03/04 are identical either way.
#
# Inputs : output/sparse_be_dataset.csv   (from 01)
#          src/poppk_model.jl             (the @model)
# Outputs: output/poppk_param_estimates.csv   population params + OFV per arm
#          output/poppk_icoef.csv             per-subject Ka/CL/Vc per arm
#          output/poppk_ipred_profiles.csv    dense individual predictions
#
# Run on a licensed Pumas / JuliaHub environment:
#   julia --project=. scripts/02_fit_poppk.jl
# =============================================================================

using Pumas, CSV, DataFrames, Statistics

const HERE = @__DIR__
const ROOT = normpath(joinpath(HERE, ".."))
const OUT  = get(ENV, "POPPK_OUT", joinpath(ROOT, "output"));  mkpath(OUT)

include(joinpath(ROOT, "src", "poppk_model.jl"))   # defines `poppk_model`

# Dense grid for the individual predictions (0-24 h to MATCH the PBPK profile
# window — same truncation on both sides => fair AUC comparison downstream).
const DENSE = collect(0.0:0.05:24.0)

data = CSV.read(joinpath(OUT, "sparse_be_dataset.csv"), DataFrame)
t50s = sort(unique(data.T50))

param_rows   = DataFrame()
icoef_rows   = DataFrame()
ipred_rows   = DataFrame()
se_rows      = DataFrame()
shrink_rows  = DataFrame()

for t50 in t50s
    arm = filter(:T50 => ==(t50), data)
    pop = read_pumas(arm;
                     id          = :SUBJID,
                     time        = :time,
                     amt         = :amt,
                     evid        = :evid,
                     cmt         = :cmt,
                     observations = [:dv],
                     covariates   = [:T50, :IndividualId])

    @info "Fitting T50 = $t50 h  (n = $(length(pop)) subjects)"
    # tvlag held constant (~0): the lag is barely identifiable from the sparse BE
    # schedule and its discontinuity destabilises the gradient; the Weibull shape
    # absorbs the absorption onset. (See RUN_ON_PUMAS.md, version-sensitive spots.)
    fpm = fit(poppk_model, pop, init_params(poppk_model), FOCE();
              constantcoef = (:tvlag,))

    c   = coef(fpm)
    ofv = try; deviance(fpm); catch; missing; end
    push!(param_rows, (T50 = t50, tvalpha = c.tvalpha, tvgamma = c.tvgamma,
                       tvcl = c.tvcl, tvvc = c.tvvc, tvlag = c.tvlag,
                       omega_alpha = c.Ω[1,1], omega_cl = c.Ω[2,2],
                       omega_vc = c.Ω[3,3], sigma_p = c.σ_p, OFV = ofv);
          cols = :union)

    # ---- standard errors (asymptotic, sandwich) and shrinkage ---------------
    try
        ct = DataFrame(coeftable(infer(fpm)))
        ct[!, :fit] .= "per-arm"
        insertcols!(ct, 1, :T50 => t50)
        append!(se_rows, ct; cols = :union)
    catch err
        @warn "infer() failed" t50 err
    end
    ηs = ηshrinkage(fpm).η
    ϵs = try; only(values(ϵshrinkage(fpm))); catch; missing; end
    push!(shrink_rows, (fit = "per-arm", T50 = t50, eta_shrink_alpha = ηs[1],
                        eta_shrink_cl = ηs[2], eta_shrink_vc = ηs[3], eps_shrink = ϵs); cols = :union)

    # ---- per-subject individual coefficients (Alpha, CL, Vc) via EBEs -------
    ic = icoef(fpm)                       # one entry per subject
    for (subj, ico) in zip(pop, ic)
        icv = ico(0.0)                    # icoef entries are callable interpolants
        push!(icoef_rows, (T50 = t50,
                           IndividualId = Int(subj.covariates(0.0).IndividualId),
                           SUBJID = subj.id,
                           Alpha = icv.Alpha, CL = icv.CL, Vc = icv.Vc,
                           gamma = c.tvgamma, lag = c.tvlag);
              cols = :union)
    end

    # ---- dense individual predictions (IPRED uses the EBEs) -----------------
    preds = predict(fpm; obstimes = DENSE)
    for (subj, pr) in zip(pop, preds)
        id   = Int(subj.covariates(0.0).IndividualId)
        ipr  = pr.ipred.dv               # individual prediction of dv (= cp)
        tt   = pr.time
        for (tk, ck) in zip(tt, ipr)
            push!(ipred_rows, (T50 = t50, IndividualId = id, time = tk, ipred = ck);
                  cols = :union)
        end
    end
end

CSV.write(joinpath(OUT, "poppk_param_estimates.csv"),  param_rows)
CSV.write(joinpath(OUT, "poppk_icoef.csv"),            icoef_rows)
CSV.write(joinpath(OUT, "poppk_ipred_profiles.csv"),   ipred_rows)

CSV.write(joinpath(OUT, "poppk_param_se.csv"),         se_rows)
CSV.write(joinpath(OUT, "poppk_shrinkage.csv"),        shrink_rows)

println("\nPopulation parameter estimates per formulation:")
show(param_rows; allrows = true, allcols = true); println()
println("\nWrote: poppk_param_estimates.csv, poppk_icoef.csv, poppk_ipred_profiles.csv -> $OUT")
