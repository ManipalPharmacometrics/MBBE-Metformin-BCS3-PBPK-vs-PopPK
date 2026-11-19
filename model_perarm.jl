# =============================================================================
# model_perarm.jl — PopPK model 1: Weibull release, 1-compartment, fitted to
# each formulation (T50 arm) separately.
#
#   Release from the tablet follows a Weibull curve:
#       fraction released by time t = 1 - exp(-(t/Alpha)^Gamma)
#   written as a release rate (hazard) on a Depot compartment:
#       h(t)     = (Gamma/Alpha) * (t/Alpha)^(Gamma-1)
#       Depot'   = -h(t) * Depot                       Depot(0) = dose
#       Central' =  h(t) * Depot - (CL/Vc) * Central
#       Cp       =  Central / Vc
#   The released drug goes straight into the blood (metformin XR is
#   release-limited, so there is no separate absorption rate ka).
#
#   Parameters: Alpha (release time scale, h), Gamma (shape), CL (L/h), Vc (L),
#   lag (h; fixed at ~0 in the fit).
#   Between-subject variability: Alpha, CL, Vc (log-normal, diagonal Omega).
#   Residual error: proportional.
#   Units: time h, amount µmol, concentration µmol/L.
# =============================================================================

poppk_perarm = @model begin
    @param begin
        tvalpha ∈ RealDomain(lower = 0.0, init = 3.0)     # release time scale (h)
        tvgamma ∈ RealDomain(lower = 0.1, init = 1.5)     # Weibull shape
        tvcl    ∈ RealDomain(lower = 0.0, init = 11.0)    # clearance (L/h)
        tvvc    ∈ RealDomain(lower = 0.0, init = 120.0)   # volume (L)
        tvlag   ∈ RealDomain(lower = 0.0, init = 1e-3)    # lag (h), held fixed
        Ω       ∈ PDiagDomain(init = [0.09, 0.09, 0.09])  # IIV: Alpha, CL, Vc
        σ_p     ∈ RealDomain(lower = 0.0, init = 0.15)    # proportional error SD
    end

    @random begin
        η ~ MvNormal(Ω)
    end

    @pre begin
        Alpha = tvalpha * exp(η[1])
        CL    = tvcl * exp(η[2])
        Vc    = tvvc * exp(η[3])
        Gamma = tvgamma
        Tlag  = tvlag
    end

    @dynamics begin
        # Weibull release rate; max(., 1e-8) keeps the power finite at t -> 0
        Depot'   = -ifelse(t > Tlag, (Gamma/Alpha) * (max(t - Tlag, 1e-8)/Alpha)^(Gamma - 1), 0.0) * Depot
        Central' =  ifelse(t > Tlag, (Gamma/Alpha) * (max(t - Tlag, 1e-8)/Alpha)^(Gamma - 1), 0.0) * Depot - (CL / Vc) * Central
    end

    @derived begin
        cp := @. Central / Vc
        dv ~ @. Normal(cp, abs(cp) * σ_p + 1e-6)          # proportional error
    end
end
