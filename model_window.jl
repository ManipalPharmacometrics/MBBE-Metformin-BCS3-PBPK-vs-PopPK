# =============================================================================
# model_window.jl — "known release + absorption window" PopPK model (step 9).
#
# The release is NOT estimated: it is the formulation's own dissolution curve,
# given as data (T50 covariate; shape 0.92 and no lag, as set in PK-Sim):
#     fraction released = 1 - exp(-ln2 * (t/T50)^0.92)
#     release rate h(t) = ln2 * 0.92 * t^(0.92-1) / T50^0.92
# What the model estimates is what happens after release, which is the same
# for every formulation:
#     absorption rate = ka * exp(-t / tau)
#         metformin is absorbed mainly high in the small intestine; as the
#         tablet moves down, absorptive capacity falls with time since dosing
#         (tau, h). Drug released late is therefore partly never absorbed,
#         so slower release lowers F, and Cmax more than AUC.
#     1-compartment disposition (CL, Vc), shared by all formulations.
#
#   Tablet' = -h(t) * Tablet                              Tablet(0) = dose
#   Gut'    =  h(t) * Tablet - ka*exp(-t/tau) * Gut       dissolved, not yet absorbed
#   Central'=  ka*exp(-t/tau) * Gut - (CL/Vc) * Central
#
# IIV: tau, CL, Vc. Proportional error. Units: h, µmol, L, µmol/L.
# To predict an untested formulation: set its T50, keep everything else.
# =============================================================================

poppk_window = @model begin
    @param begin
        tvka  ∈ RealDomain(lower = 0.01, init = 1.0)     # absorption rate at t = 0 (1/h)
        tvtau ∈ RealDomain(lower = 0.1,  init = 4.0)     # absorption-window time constant (h)
        tvcl  ∈ RealDomain(lower = 0.0,  init = 40.0)    # clearance (L/h)
        tvvc  ∈ RealDomain(lower = 0.0,  init = 150.0)   # volume (L)
        Ω     ∈ PDiagDomain(init = [0.09, 0.09, 0.09])   # IIV: tau, CL, Vc
        σ_p   ∈ RealDomain(lower = 0.0,  init = 0.15)    # proportional error SD
    end

    @random begin
        η ~ MvNormal(Ω)
    end

    @covariates T50

    @pre begin
        Ka   = tvka
        Tau  = tvtau * exp(η[1])
        CL   = tvcl * exp(η[2])
        Vc   = tvvc * exp(η[3])
        Rel  = log(2) * 0.92 / T50^0.92                 # known release: h(t) = Rel * t^(0.92-1)
    end

    @dynamics begin
        Tablet'  = -Rel * max(t, 1e-8)^(0.92 - 1) * Tablet
        Gut'     =  Rel * max(t, 1e-8)^(0.92 - 1) * Tablet - Ka * exp(-t / Tau) * Gut
        Central' =  Ka * exp(-t / Tau) * Gut - (CL / Vc) * Central
    end

    @derived begin
        cp := @. Central / Vc
        dv ~ @. Normal(cp, abs(cp) * σ_p + 1e-6)
    end
end
