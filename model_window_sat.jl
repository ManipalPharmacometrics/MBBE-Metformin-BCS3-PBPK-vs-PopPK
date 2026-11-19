# =============================================================================
# model_window_sat.jl — model_window.jl with SATURABLE uptake (step 9).
#
# Same as model_window.jl (known Weibull release, absorption window, shared
# 1-compartment disposition), except that uptake from the gut saturates.
# Metformin is taken up by transporters (OCT3, PMAT, THTR2), so when fast
# release puts a lot of drug in the gut at once, uptake cannot keep up:
#
#   uptake = ka * Km * exp(-t/tau) * Gut / (Km + Gut)        (µmol/h)
#     Gut << Km : ka * exp(-t/tau) * Gut      = model_window.jl (linear)
#     Gut >> Km : ka * Km * exp(-t/tau)       = capped at Vmax = ka*Km
#
#   Tablet' = -h(t) * Tablet                              Tablet(0) = dose
#   Gut'    =  h(t) * Tablet - uptake
#   Central'=  uptake - (CL/Vc) * Central
#
# Km = gut amount (µmol) at half-maximal uptake; dose = 3019.5 µmol.
# ka is fixed in step 9 (see model_window.jl). IIV: tau, CL, Vc.
# Units: h, µmol, L, µmol/L.
# =============================================================================

poppk_window_sat = @model begin
    @param begin
        tvka  ∈ RealDomain(lower = 0.01, init = 0.3)     # uptake rate constant at low gut amount (1/h)
        tvkm  ∈ RealDomain(lower = 1.0,  init = 500.0)   # saturation amount (µmol)
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
        Km   = tvkm
        Tau  = tvtau * exp(η[1])
        CL   = tvcl * exp(η[2])
        Vc   = tvvc * exp(η[3])
        Rel  = log(2) * 0.92 / T50^0.92                 # known release: h(t) = Rel * t^(0.92-1)
    end

    @dynamics begin
        Tablet'  = -Rel * max(t, 1e-8)^(0.92 - 1) * Tablet
        Gut'     =  Rel * max(t, 1e-8)^(0.92 - 1) * Tablet - Ka * Km * exp(-t / Tau) * Gut / (Km + Gut)
        Central' =  Ka * Km * exp(-t / Tau) * Gut / (Km + Gut) - (CL / Vc) * Central
    end

    @derived begin
        cp := @. Central / Vc
        dv ~ @. Normal(cp, abs(cp) * σ_p + 1e-6)
    end
end
