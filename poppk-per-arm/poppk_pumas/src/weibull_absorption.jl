# =============================================================================
# weibull_absorption.jl  — OPEN-JULIA (no Pumas) simulator + per-subject fit of
# the Weibull-release 1-compartment absorption model. This mirrors exactly the
# structural model in src/poppk_model.jl (the Pumas @model), but uses only open
# packages so the whole PBPK->PopPK safe-space pipeline can be regenerated and
# cross-checked WITHOUT a Pumas license.
#
# Model (same as poppk_model.jl):
#   tau  = t - lag
#   h(t) = (gamma/alpha) * (tau/alpha)^(gamma-1)                 (tau > 0, else 0)
#   Depot'   = -h(t) * Depot         Depot(0) = Dose             (Weibull release)
#   Central' =  h(t) * Depot - (CL/Vc)*Central
#   C(t)     =  Central / Vc
# => Depot(t) = Dose * exp(-((t-lag)/alpha)^gamma), the Weibull release law used
#    on the PBPK (PK-Sim) side; alpha ~ PBPK T50, gamma ~ PBPK beta.
#
# The Pumas step 02 does a proper population (NLME) fit. This fallback does a
# per-subject (two-stage) nonlinear least-squares fit on the log scale (matching
# the proportional residual). It is an APPROXIMATION of the NLME fit, adequate
# for regenerating the decision-equivalence overlay and for validating that the
# fitted release scale recovers the PBPK T50.
#
# `include`d (not a module) after `using Statistics, Optim`.
# =============================================================================

"Weibull release hazard h(t) (time-varying first-order release rate, 1/h)."
@inline function wb_hazard(t, alpha, gamma, lag)
    tau = t - lag
    tau <= 0 && return 0.0
    (gamma / alpha) * (max(tau, 1e-8) / alpha)^(gamma - 1)
end

"""
    simulate_weibull(tout, alpha, gamma, cl, vc, lag, dose; nsub=2) -> Vector

Concentration (µmol/L) at the sorted output times `tout` (h) for the Weibull-
release 1-cmt model, by fixed-step RK4. `nsub` RK4 substeps are taken between
consecutive output times (2 is ample here: ke≈0.1/h, smooth release).
"""
function simulate_weibull(tout::AbstractVector, alpha, gamma, cl, vc, lag, dose; nsub::Int = 2)
    ke = cl / vc
    f = function (u, t)
        h = wb_hazard(t, alpha, gamma, lag)
        (-h * u[1], h * u[1] - ke * u[2])          # (Depot', Central')
    end
    u1 = dose; u2 = 0.0                            # Depot, Central
    out = Vector{Float64}(undef, length(tout))
    t = 0.0
    # value at the first output time (assumed >= 0)
    if tout[1] <= 0
        out[1] = u2 / vc; istart = 2
    else
        istart = 1
    end
    tprev = 0.0
    for i in istart:length(tout)
        tnext = tout[i]
        h = (tnext - tprev) / nsub
        for _ in 1:nsub
            d1 = f((u1, u2), t)
            d2 = f((u1 + 0.5h*d1[1], u2 + 0.5h*d1[2]), t + 0.5h)
            d3 = f((u1 + 0.5h*d2[1], u2 + 0.5h*d2[2]), t + 0.5h)
            d4 = f((u1 + h*d3[1],    u2 + h*d3[2]),    t + h)
            u1 += (h/6)*(d1[1] + 2d2[1] + 2d3[1] + d4[1])
            u2 += (h/6)*(d1[2] + 2d2[2] + 2d3[2] + d4[2])
            t  += h
        end
        out[i] = u2 / vc
        tprev = tnext
    end
    out
end

"Trapezoidal AUC (for the CL starting guess)."
function _trapz(t, c)
    a = 0.0
    @inbounds for i in 1:length(t)-1
        a += (t[i+1]-t[i]) * (c[i]+c[i+1]) / 2
    end
    a
end

"""
    fit_subject_weibull(tobs, cobs, dose; lag=0.25, restarts=3) -> NamedTuple

Per-subject NLS fit of (alpha, gamma, CL, Vc) on the log scale (log-least-
squares == proportional residual). `lag` is fixed (barely identifiable from a BE
schedule starting at 0.5 h). Returns fitted params + the SSE and a converged flag.
"""
function fit_subject_weibull(tobs::AbstractVector, cobs::AbstractVector, dose;
                             lag = 0.25, alpha0 = 3.0, restarts = 3)
    keep = (cobs .> 0) .& isfinite.(cobs)
    t = Float64.(tobs[keep]); y = Float64.(cobs[keep])
    ly = log.(y)
    # data-driven starting values
    area = max(_trapz(t, y), 1e-6)
    cl0  = clamp(dose / area, 1.0, 100.0)          # AUC ~ Dose/CL (complete input)
    vc0  = clamp(cl0 / 0.10, 20.0, 600.0)          # ke ~ 0.1/h typical for metformin
    g0   = 1.5

    obj = function (θ)                              # θ = log([alpha,gamma,cl,vc])
        a, g, cl, vc = exp.(θ)
        (a < 1e-3 || g < 0.1 || cl <= 0 || vc <= 0) && return 1e12
        pred = simulate_weibull(t, a, g, cl, vc, lag, dose)
        any(p -> p <= 0 || !isfinite(p), pred) && return 1e12
        sum(abs2, log.(pred) .- ly)
    end

    inits = [log.([alpha0, g0, cl0, vc0]),
             log.([alpha0*0.6, 1.0, cl0*1.3, vc0*0.7]),
             log.([alpha0*1.6, 2.2, cl0*0.7, vc0*1.4])]
    best = nothing; bestf = Inf
    for k in 1:min(restarts, length(inits))
        r = optimize(obj, inits[k], NelderMead(),
                     Optim.Options(iterations = 2000, g_tol = 1e-8))
        if Optim.minimum(r) < bestf
            bestf = Optim.minimum(r); best = r
        end
    end
    a, g, cl, vc = exp.(Optim.minimizer(best))
    (alpha = a, gamma = g, CL = cl, Vc = vc, lag = lag,
     sse = bestf, converged = Optim.converged(best))
end
