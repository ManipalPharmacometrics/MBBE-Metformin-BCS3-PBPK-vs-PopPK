# =============================================================================
# poppk_model.jl  — baseline PopPK structural model for the PBPK -> PopPK
# decision-equivalence comparison.
#
# BASELINE (2026-07-22):  WEIBULL-RELEASE absorption, JOINT fit across all arms
# with FORMULATION-INDEPENDENT DISPOSITION.
#
# Why joint / shared disposition: metformin XR is release-rate-limited, so in a
# 1-compartment model the release and the elimination are NOT separately
# identifiable from a single arm (the flip-flop ambiguity) — fitting each arm its
# own CL/Vc lets the formulation effect leak into disposition and leaves the
# release scale flat. Clearance and volume are properties of the drug+subject,
# NOT of the tablet's release rate, so we hold them formulation-independent
# (shared across all 8 arms) and let ONLY the release vary with formulation:
#
#     Alpha(T50) = tvalpha * (T50 / REF_T50)^slope * exp(η)
#
# `slope` is the estimated in-vivo IVIVC exponent (release scale vs T50). REF_T50
# = 3.0 h. Disposition (CL, Vc) is shared/typical; the Weibull shape `tvgamma` is
# a single SHARED ESTIMATED value across arms, and lag is held constant. NB gamma
# here is the IN-VIVO INPUT shape (in-vitro dissolution folded with gut absorption
# and transit) — it is NOT the PBPK dissolution shape β=1.0 and should not be
# fixed to it (doing so is too rigid and inflates the safe space). Only the
# release SCALE (Alpha) varies with formulation, giving an identifiable IVIVC.
#
# Mechanism (time-varying Weibull hazard on a releasing depot):
#   tau  = t - lag
#   h(t) = (gamma/alpha) * (tau/alpha)^(gamma-1)                 (tau > 0, else 0)
#   Depot'   = -h(t) * Depot         Depot(0) = Dose             (Weibull release)
#   Central' =  h(t) * Depot - (CL/Vc) * Central
#   C(t)     =  Central / Vc
# so Depot(t) = Dose * exp(-((t-lag)/alpha)^gamma), the same Weibull release the
# PBPK (PK-Sim) tablet uses.
#
# Units: time h, amount µmol, alpha h, Vc L, CL L/h, Central/Vc -> µmol/L.
#
# Pumas-compatibility notes (from the licensed-host run, 2026-07):
#   * @dynamics uses ifelse(...) (this Pumas rejects the `cond ? a : b` form).
#   * tvlag sits on a log-transformed lower bound, so init = 1e-3 (not 0) and it
#     is held constant in the fit via constantcoef=(:tvlag,); the Weibull shape
#     absorbs the onset.
#   * icoef returns callable interpolants: read EBEs as `ico(0.0).Alpha` (see
#     scripts/02_fit_poppk.jl).
#
# Sensitivity model (see README): the same model fit PER ARM with slope held
# constant (Alpha = tvalpha per arm) — the earlier independent-arm variant.
#
# `include`d after `using Pumas`.
# =============================================================================

const REF_T50_MODEL = 3.0    # reference formulation T50 (h) for the IVIVC scaling

poppk_model = @model begin

    @param begin
        tvalpha ∈ RealDomain(lower = 0.0,  init = 3.0)   # Weibull release scale (h) at REF T50
        slope   ∈ RealDomain(init = 1.0)                 # IVIVC exponent: Alpha ∝ (T50/REF)^slope
        tvgamma ∈ RealDomain(lower = 0.1,  init = 1.2)   # in-vivo INPUT shape (shared, estimated)
        tvcl    ∈ RealDomain(lower = 0.0,  init = 11.0)  # clearance (L/h), formulation-independent
        tvvc    ∈ RealDomain(lower = 0.0,  init = 120.0) # central volume (L), formulation-independent
        tvlag   ∈ RealDomain(lower = 0.0,  init = 1e-3)  # release lag (h); held constant in the fit
        Ω       ∈ PDiagDomain(init = [0.09, 0.09, 0.09]) # IIV variances: Alpha, CL, Vc
        σ_p     ∈ RealDomain(lower = 0.0,  init = 0.15)  # proportional residual SD
    end

    @random begin
        η ~ MvNormal(Ω)
    end

    @covariates T50

    @pre begin
        Alpha = tvalpha * (T50 / REF_T50_MODEL)^slope * exp(η[1])   # in-vivo release scale
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
