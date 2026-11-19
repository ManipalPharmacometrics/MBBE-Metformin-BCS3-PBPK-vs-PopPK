# =============================================================================
# model_window_ge.jl — model_window.jl with a GASTRIC EMPTYING step (step 9).
#
# Same as model_window.jl (known Weibull release, absorption window, shared
# 1-compartment disposition), except that drug released in the stomach is not
# absorbed there: it must first empty into the small intestine (first order,
# kge). However fast a tablet releases, drug cannot reach the absorption site
# faster than the stomach empties, which should damp the Cmax effect of fast
# release.
#
#   Tablet'  = -h(t) * Tablet                              Tablet(0) = dose
#   Stomach' =  h(t) * Tablet - kge * Stomach              released, still in the stomach
#   Gut'     =  kge * Stomach - ka*exp(-t/tau) * Gut       in the small intestine
#   Central' =  ka*exp(-t/tau) * Gut - (CL/Vc) * Central
#
# kge -> infinity gives model_window.jl back, so the two -2LL values compare.
# ka is fixed in step 9 (see model_window.jl). IIV: tau, CL, Vc.
# Units: h, µmol, L, µmol/L.
# =============================================================================

poppk_window_ge = @model begin
    @param begin
        tvka  ∈ RealDomain(lower = 0.01, init = 0.3)     # absorption rate at t = 0 (1/h)
        tvkge ∈ RealDomain(lower = 0.05, init = 1.5)     # gastric emptying rate (1/h)
        tvtau ∈ RealDomain(lower = 0.1,  init = 5.0)     # absorption-window time constant (h)
        tvcl  ∈ RealDomain(lower = 0.0,  init = 25.0)    # clearance (L/h)
        tvvc  ∈ RealDomain(lower = 0.0,  init = 50.0)    # volume (L)
        Ω     ∈ PDiagDomain(init = [0.09, 0.09, 0.09])   # IIV: tau, CL, Vc
        σ_p   ∈ RealDomain(lower = 0.0,  init = 0.15)    # proportional error SD
    end

    @random begin
        η ~ MvNormal(Ω)
    end

    @covariates T50

    @pre begin
        Ka   = tvka
        Kge  = tvkge
        Tau  = tvtau * exp(η[1])
        CL   = tvcl * exp(η[2])
        Vc   = tvvc * exp(η[3])
        Rel  = log(2) * 0.92 / T50^0.92                 # known release: h(t) = Rel * t^(0.92-1)
    end

    @dynamics begin
        Tablet'  = -Rel * max(t, 1e-8)^(0.92 - 1) * Tablet
        Stomach' =  Rel * max(t, 1e-8)^(0.92 - 1) * Tablet - Kge * Stomach
        Gut'     =  Kge * Stomach - Ka * exp(-t / Tau) * Gut
        Central' =  Ka * exp(-t / Tau) * Gut - (CL / Vc) * Central
    end

    @derived begin
        cp := @. Central / Vc
        dv ~ @. Normal(cp, abs(cp) * σ_p + 1e-6)
    end
end
