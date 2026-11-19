# =============================================================================
# settings.jl — every number the analysis uses, in one place.
# Each step script starts with include("settings.jl").
# =============================================================================

const HERE    = @__DIR__
const DATA    = joinpath(HERE, "data")
const RESULTS = joinpath(HERE, "results")

# ---- the drug and the study -------------------------------------------------
# 500 mg metformin HCl XR tablet = 390 mg free base (PK-Sim doses the base).
const DOSE_MG   = 390.0
const MW        = 129.16                         # metformin base, g/mol
const DOSE_UMOL = DOSE_MG / MW * 1000            # = 3019.5 µmol (model amount unit)

# Rich sampling schedule (h). The dose is at 0 h.
const SCHEDULE = [0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 16.0, 24.0]
const PROP_CV  = 0.15                            # 15% proportional assay error
const SEED_DATASET = 20260610                    # random seed for the assay error

# ---- formulations -----------------------------------------------------------
# PK-Sim simulated 8 XR formulations that differ only in the in vivo release
# time T50 (h). T50 = 3.0 h is the reference product.
const REF_T50 = 3.0

# ---- bioequivalence ---------------------------------------------------------
const BE_LOWER = 0.80
const BE_UPPER = 1.25
const ENDPOINTS = ["Cmax", "AUC_tEnd"]            # AUC_tEnd = AUC 0-24 h

# Within-subject CV used to simulate BE trials (FDA NDA 212595 clinical
# pharmacology review, metformin XR): Cmax 12.2%, AUC0-t 11.2%.
const WSV_CV = Dict("Cmax" => 0.122, "AUC_tEnd" => 0.112)

const N_SUBJ = [48, 86]                          # BE study sizes simulated
const SEED_TRIALS = 1234                         # random seed for the virtual BE trials

# Step 5 (power and safe space)
const N_TRIAL_POWER = 500                        # virtual trials per T50 value
const T50_STEP_POWER = 0.05                      # T50 grid step (h)
const POWER_TARGET  = 0.80                       # safe space = power >= 80%

# Step 6 (type I error)
const N_TRIAL_TYPE1 = 20_000                     # SE of a 5% rate ~0.15 percentage points
const T50_STEP_TYPE1 = 0.005                     # finer grid to locate the boundary
const N_TRIAL_OC    = 2_000                      # trials per point of the operating-characteristic figure

# ---- the two PopPK models ---------------------------------------------------
const MODELS = ["perarm", "joint"]               # model_perarm.jl, model_joint.jl
