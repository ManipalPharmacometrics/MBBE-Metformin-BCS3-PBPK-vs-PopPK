# =============================================================================
# poppk_model.jl  — baseline PopPK structural model for the PBPK -> PopPK
# decision-equivalence comparison.
#
# BASELINE (original plan, 2026-07-22):  WEIBULL-RELEASE ABSORPTION.
#   1-compartment disposition whose INPUT is a Weibull in-vivo release process,
#   mirroring the PK-Sim mechanism directly. On the PBPK side the metformin XR
#   tablet dissolves per a Weibull profile
#         fraction released F(t) = 1 - exp(-((t-lag)/T50)^beta)
#   (T50 = Weibull scale / release time, beta = shape). Metformin XR is
#   dissolution/release rate-limited (flip-flop kinetics), so we let that same
#   Weibull release BE the absorption input — no separate first-order ka is
#   needed or identifiable. The population fit therefore estimates the in-vivo
#   Weibull release parameters (alpha ~ T50 analog, gamma ~ beta analog) plus
#   disposition (CL, Vc), which is the closest PopPK mirror of the PK-Sim input.
#
# Mechanism (time-varying Weibull hazard on a releasing depot):
#   tau = t - tlag
#   h(t) = (gamma/alpha) * (tau/alpha)^(gamma-1)          # Weibull hazard, tau>0
#   Depot'   = -h(t) * Depot           Depot(0) = Dose    # amount not yet released
#   Central' =  h(t) * Depot - (CL/Vc) * Central          # 1-cmt disposition
# so Depot(t) = Dose * exp(-((t-tlag)/alpha)^gamma) reproduces the Weibull
# release exactly, and alpha/gamma are directly comparable to the PBPK T50/beta.
#
# Units (kept consistent with the PK-Sim exports):
#   time     hours
#   amount   µmol            (dose = 390 mg free base / 129.16 g/mol -> µmol)
#   alpha    hours           (in-vivo Weibull release scale; PBPK-T50 analog)
#   gamma    dimensionless   (Weibull shape; PBPK-beta analog)
#   Vc       L               CL  L/h        Central/Vc -> µmol/L
#
# IIV on alpha, CL, Vc (diagonal Omega). Proportional residual error.
# gamma (shape) and lag are fixed effects (no IIV) for estimation stability from
# the sparse BE design; add IIV later if the fit supports it.
#
# Documented alternative (see README / walkthrough §11): the simpler 1-cmt
# FIRST-ORDER absorption + lag model, where T50 is absorbed implicitly into an
# apparent ka. Kept as a sensitivity model, not the baseline.
#
# `include`d after `using Pumas`.
# =============================================================================

poppk_model = @model begin

    @param begin
        tvalpha ∈ RealDomain(lower = 0.0,  init = 3.0)   # typical Weibull release scale (h) ~ T50
        tvgamma ∈ RealDomain(lower = 0.1,  init = 1.5)   # typical Weibull shape       ~ beta
        tvcl    ∈ RealDomain(lower = 0.0,  init = 11.0)  # typical clearance (L/h)
        tvvc    ∈ RealDomain(lower = 0.0,  init = 120.0) # typical central volume (L)
        tvlag   ∈ RealDomain(lower = 0.0,  init = 1e-3)  # release lag time (h); held ~0 (Weibull shape absorbs onset)
        Ω       ∈ PDiagDomain(init = [0.09, 0.09, 0.09]) # IIV variances: alpha, CL, Vc
        σ_p     ∈ RealDomain(lower = 0.0,  init = 0.15)  # proportional residual SD
    end

    @random begin
        η ~ MvNormal(Ω)
    end

    @pre begin
        Alpha = tvalpha * exp(η[1])       # in-vivo Weibull release scale (h)
        CL    = tvcl    * exp(η[2])
        Vc    = tvvc    * exp(η[3])
        Gamma = tvgamma                   # shape: population fixed effect
        Tlag  = tvlag                     # release lag: population fixed effect
    end

    @dynamics begin
        # Weibull hazard release: tau = t - Tlag; release rate constant
        #   h(t) = (Gamma/Alpha) * (tau/Alpha)^(Gamma-1)   for tau > 0, else 0.
        # The (max(.,1e-8)) guard keeps the power finite at tau -> 0 when Gamma<1.
        Depot'   = -ifelse(t > Tlag, (Gamma/Alpha) * (max(t - Tlag, 1e-8)/Alpha)^(Gamma - 1), 0.0) * Depot
        Central' =  ifelse(t > Tlag, (Gamma/Alpha) * (max(t - Tlag, 1e-8)/Alpha)^(Gamma - 1), 0.0) * Depot - (CL / Vc) * Central
    end

    @derived begin
        cp := @. Central / Vc
        # Proportional error; tiny floor keeps the SD strictly positive if the
        # optimiser transiently pushes a predicted concentration to ~0.
        dv ~ @. Normal(cp, abs(cp) * σ_p + 1e-6)
    end
end
