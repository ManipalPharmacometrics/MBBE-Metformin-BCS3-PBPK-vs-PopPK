# =============================================================================
# 02_fit_poppk.jl  — JOINT fit of the Weibull-release PopPK model (baseline).
#
# ONE population = all 8 formulation arms (192 subject-arms). Disposition
# (CL, Vc), Weibull shape (gamma) and lag are shared; the release scale varies
# with formulation through the covariate model
#       Alpha(T50) = tvalpha * (T50/REF)^slope * exp(η).
# This holds disposition formulation-independent (physiologically correct) and
# resolves the flip-flop ambiguity, so `slope`/`Alpha` are an identifiable
# in-vivo IVIVC. Everything downstream (03/04) is unchanged.
#
# Inputs : output/sparse_be_dataset.csv   (from 01)
#          src/poppk_model.jl             (the @model, joint covariate form)
# Outputs: output/poppk_param_estimates.csv   population params + slope + OFV
#          output/poppk_icoef.csv             per-subject Alpha/CL/Vc per arm
#          output/poppk_ipred_profiles.csv    dense individual predictions
#
# Run on a licensed Pumas / JuliaHub environment:
#   julia --project=. scripts/02_fit_poppk.jl
#
# This is the JOINT (IVIVC) fit. The per-arm baseline lives in the other export.
# =============================================================================

using Pumas, CSV, DataFrames, Statistics

const HERE = @__DIR__
const ROOT = normpath(joinpath(HERE, ".."))
const OUT  = get(ENV, "POPPK_OUT", joinpath(ROOT, "output"));  mkpath(OUT)

include(joinpath(ROOT, "src", "poppk_model.jl"))   # defines `poppk_model`, REF_T50_MODEL

# Dense grid for the individual predictions (0-24 h to MATCH the PBPK profile
# window — same truncation on both sides => fair AUC comparison downstream).
const DENSE = collect(0.0:0.05:24.0)

data = CSV.read(joinpath(OUT, "sparse_be_dataset.csv"), DataFrame)

# ONE population over all arms; T50 and IndividualId carried as covariates.
pop = read_pumas(data;
                 id           = :SUBJID,
                 time         = :time,
                 amt          = :amt,
                 evid         = :evid,
                 cmt          = :cmt,
                 observations = [:dv],
                 covariates   = [:T50, :IndividualId])

@info "Joint fit over $(length(pop)) subject-arms (8 formulations, shared disposition)"
fpm = fit(poppk_model, pop, init_params(poppk_model), FOCE();
          constantcoef = (:tvlag,))            # lag pinned; γ estimated (shared in-vivo input shape)

c   = coef(fpm)
ofv = try; deviance(fpm); catch; missing; end

param_rows = DataFrame([(tvalpha = c.tvalpha, slope = c.slope, tvgamma = c.tvgamma,
                        tvcl = c.tvcl, tvvc = c.tvvc, tvlag = c.tvlag,
                        omega_alpha = c.Ω[1,1], omega_cl = c.Ω[2,2],
                        omega_vc = c.Ω[3,3], sigma_p = c.σ_p, OFV = ofv)])

se_rows = DataFrame(); shrink_rows = DataFrame()

# ---- standard errors (asymptotic, sandwich) and shrinkage ---------------
try
    ct = DataFrame(coeftable(infer(fpm)))
    ct[!, :fit] .= "joint"
    insertcols!(ct, 1, :T50 => missing)
    append!(se_rows, ct; cols = :union)
catch err
    @warn "infer() failed" err
end
ηs = ηshrinkage(fpm).η
ϵs = try; only(values(ϵshrinkage(fpm))); catch; missing; end
push!(shrink_rows, (fit = "joint", T50 = missing, eta_shrink_alpha = ηs[1],
                    eta_shrink_cl = ηs[2], eta_shrink_vc = ηs[3], eps_shrink = ϵs); cols = :union)

# ---- per-subject individual coefficients (Alpha, CL, Vc) via EBEs -----------
# icoef returns callable interpolants here; read at t=0.
ic = icoef(fpm)
icoef_rows = DataFrame()
for (subj, ico) in zip(pop, ic)
    a  = ico(0.0)
    push!(icoef_rows, (T50 = subj.covariates(0.0).T50,
                       IndividualId = Int(subj.covariates(0.0).IndividualId),
                       SUBJID = subj.id,
                       Alpha = a.Alpha, CL = a.CL, Vc = a.Vc,
                       gamma = c.tvgamma, lag = c.tvlag); cols = :union)
end

# ---- dense individual predictions (IPRED uses the EBEs) ---------------------
preds = predict(fpm; obstimes = DENSE)
ipred_rows = DataFrame()
for (subj, pr) in zip(pop, preds)
    id   = Int(subj.covariates(0.0).IndividualId)
    t50  = subj.covariates(0.0).T50
    ipr  = pr.ipred.dv
    tt   = pr.time
    for (tk, ck) in zip(tt, ipr)
        push!(ipred_rows, (T50 = t50, IndividualId = id, time = tk, ipred = ck); cols = :union)
    end
end

CSV.write(joinpath(OUT, "poppk_param_estimates.csv"),  param_rows)
CSV.write(joinpath(OUT, "poppk_icoef.csv"),            icoef_rows)
CSV.write(joinpath(OUT, "poppk_ipred_profiles.csv"),   ipred_rows)

CSV.write(joinpath(OUT, "poppk_param_se.csv"),         se_rows)
CSV.write(joinpath(OUT, "poppk_shrinkage.csv"),        shrink_rows)

println("\nJoint population parameter estimates (shared disposition + IVIVC slope):")
show(param_rows; allrows = true, allcols = true); println()

# ---- IVIVC check: geometric-mean fitted Alpha vs PBPK input T50 -------------
gm(x) = exp(mean(log.(x)))
rec = combine(groupby(icoef_rows, :T50), :Alpha => gm => :gm_alpha)
rec.ratio = rec.gm_alpha ./ rec.T50
println("\nIVIVC — fitted in-vivo release scale (Alpha) vs PBPK T50 (slope = $(round(c.slope, digits=3))):")
show(sort(rec, :T50); allrows = true, allcols = true); println()
println("\nWrote poppk_param_estimates.csv, poppk_icoef.csv, poppk_ipred_profiles.csv -> $OUT")
