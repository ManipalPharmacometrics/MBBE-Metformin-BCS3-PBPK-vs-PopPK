# =============================================================================
# model_joint_frel.jl — the joint model (model_joint.jl) plus a T50 effect on
# bioavailability, used in step 8 (held-out test).
#
#   Alpha = tvalpha * (T50/3)^slope * exp(eta)     release time scale, as before
#   F     = (T50/3)^theta_f                         relative bioavailability,
#                                                   1 at the reference T50 (3 h)
# Why: in PK-Sim a slower release means less metformin absorbed (AUC 0.54x from
# T50 1 to 7 h). Without F the model puts that loss into CL (step 7).
# Everything else is identical to model_joint.jl.
# =============================================================================

poppk_joint_frel = @model begin
    @param begin
        tvalpha ∈ RealDomain(lower = 0.0, init = 3.0)
        slope   ∈ RealDomain(init = 1.0)                  # Alpha ∝ (T50/3)^slope
        theta_f ∈ RealDomain(init = -0.3)                 # F = (T50/3)^theta_f
        tvgamma ∈ RealDomain(lower = 0.1, init = 1.2)
        tvcl    ∈ RealDomain(lower = 0.0, init = 11.0)
        tvvc    ∈ RealDomain(lower = 0.0, init = 120.0)
        tvlag   ∈ RealDomain(lower = 0.0, init = 1e-3)    # held fixed
        Ω       ∈ PDiagDomain(init = [0.09, 0.09, 0.09])  # IIV: Alpha, CL, Vc
        σ_p     ∈ RealDomain(lower = 0.0, init = 0.15)
    end

    @random begin
        η ~ MvNormal(Ω)
    end

    @covariates T50

    @pre begin
        Alpha = tvalpha * (T50 / 3.0)^slope * exp(η[1])
        CL    = tvcl * exp(η[2])
        Vc    = tvvc * exp(η[3])
        Gamma = tvgamma
        Tlag  = tvlag
    end

    @dosecontrol begin
        bioav = (; Depot = (T50 / 3.0)^theta_f,)          # F; written here, @dosecontrol cannot see @pre
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
