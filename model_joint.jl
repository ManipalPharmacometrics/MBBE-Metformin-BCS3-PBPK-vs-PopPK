# =============================================================================
# model_joint.jl — PopPK model 2: the same Weibull-release 1-compartment model,
# fitted to ALL formulations at once (one population of 8 x 24 subject-arms).
#
# Differences from model_perarm.jl:
#   - CL, Vc and the Weibull shape Gamma are shared by every formulation.
#   - The release time scale depends on the formulation's T50 (covariate):
#         Alpha = tvalpha * (T50 / 3.0)^slope * exp(eta)
#     3.0 h = the reference T50. `slope` is estimated (an in vivo IVIVC).
# Equations, units and error model are otherwise identical.
# =============================================================================

poppk_joint = @model begin
    @param begin
        tvalpha ∈ RealDomain(lower = 0.0, init = 3.0)     # release time scale at T50 = 3 h (h)
        slope   ∈ RealDomain(init = 1.0)                  # Alpha ∝ (T50/3)^slope
        tvgamma ∈ RealDomain(lower = 0.1, init = 1.2)     # Weibull shape (shared)
        tvcl    ∈ RealDomain(lower = 0.0, init = 11.0)    # clearance (L/h)
        tvvc    ∈ RealDomain(lower = 0.0, init = 120.0)   # volume (L)
        tvlag   ∈ RealDomain(lower = 0.0, init = 1e-3)    # lag (h), held fixed
        Ω       ∈ PDiagDomain(init = [0.09, 0.09, 0.09])  # IIV: Alpha, CL, Vc
        σ_p     ∈ RealDomain(lower = 0.0, init = 0.15)    # proportional error SD
    end

    @random begin
        η ~ MvNormal(Ω)
    end

    @covariates T50

    @pre begin
        Alpha = tvalpha * (T50 / 3.0)^slope * exp(η[1])   # 3.0 = reference T50 (h)
        CL    = tvcl * exp(η[2])
        Vc    = tvvc * exp(η[3])
        Gamma = tvgamma
        Tlag  = tvlag
    end

    @dynamics begin
        Depot'   = -ifelse(t > Tlag, (Gamma/Alpha) * (max(t - Tlag, 1e-8)/Alpha)^(Gamma - 1), 0.0) * Depot
        Central' =  ifelse(t > Tlag, (Gamma/Alpha) * (max(t - Tlag, 1e-8)/Alpha)^(Gamma - 1), 0.0) * Depot - (CL / Vc) * Central
    end

    @derived begin
        cp := @. Central / Vc
        dv ~ @. Normal(cp, abs(cp) * σ_p + 1e-6)
    end
end
