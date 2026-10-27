# PopPK-based MBBE analysis — full walkthrough (from PBPK output to decision-equivalence overlay)

**Scope of this document.** Everything that happens *after* PK-Sim produces the PBPK
population output, in complete detail: how that output is consumed, whether a model is
re-fit and re-simulated (yes — both), exactly which between-subject and residual
variabilities (BSVs) are added and where, and how the per-formulation results are
regenerated in the PopPK model and compared against the PBPK model from PK-Sim.

**Code documented:** `programs/poppk_pumas/` (Julia/Pumas). This is *bridge #1* from the
project discussion (project notes (not included), §6):

> *PK-Sim → simulate large virtual population → fit popPK model in Pumas → run MBBE
> statistics in Pumas. "PBPK-informed popPK."*

The alternative (bridge #2 — run TOST directly on the PBPK profiles, no PopPK fit) is what
Morris et al. did and is the PBPK "side" we compare against. The headline deliverable of
the PopPK pipeline is the **PBPK-vs-PopPK decision-equivalence overlay**: *do the two
methods place the bioequivalence safe-space boundaries in the same place and reach the same
BE conclusions?*

> **UPDATE 2026-07-22 — baseline absorption model is now WEIBULL-RELEASE.**
> The baseline PopPK absorption model is a **Weibull-release** model that uses the *same*
> release law the PBPK tablet uses: the depot empties by a Weibull hazard, giving
> `Depot(t) = Dose·exp(−((t−lag)/α)^γ)`. So the fitted in-vivo release scale **α is the
> direct analog of the PBPK T50** and **γ of β**, and the model can reproduce the PBPK
> profile shape. The first-order-absorption description that appears in §2/§4/§7 below is
> **retained only as a sensitivity model** — wherever it says the formulation effect is
> absorbed into an *apparent `ka`*, the Weibull baseline instead estimates **α (~T50) and
> γ (~β) directly**. The current model code is `src/poppk_model.jl`; a plain-language
> methodology is in `../METHODOLOGY_metformin_MBBE.md`; how to run the fit on Pumas is in
> `RUN_ON_PUMAS.md`. Everything downstream of the fit (steps 03/04, the BE/VBE/safe-space
> math) is unchanged.

---

## 0. A note on terminology — "TD50" vs **T50**

There is no `TD50` in this project. The critical quality attribute (CQA) being swept is
**T50** — the **Weibull median dissolution (in-vitro release) time** of the metformin XR
formulation, in hours:

```
% released(t) = 100 · (1 − exp(−(t / T50)^β))
```

T50 is the single scalar formulation knob (analogous to particle radius in the Morris LAI
paper). Larger T50 = slower release. **Important:** T50 is *not* a PK parameter that gets
"regenerated" or estimated by the PopPK model. It is a **formulation label / covariate**.
See §7 for exactly how it is handled on the PK-Sim side vs the PopPK side — that is the
heart of your question.

---

## 1. What "the PBPK output" actually is (the starting point)

PK-Sim (Phase 4 of the master plan) produced, for **8 formulations** (T50 = 1.0, 2.0, 3.0,
3.5, 4.0, 5.0, 6.0, 7.0 h), a **population simulation** of **24 virtual subjects** each. Two
kinds of files were exported and both live under `programs/poppk_pumas/data/`:

| File pattern | What it holds | Used by |
|---|---|---|
| `popsims/Metformin_T50_<x>h_Simpop-Results.csv` | **Dense concentration–time profiles.** One row per subject per time point on a **3-minute grid, 0–24 h** (≈3097 rows = 24 subj × ~129 times). Key column = peripheral-venous plasma metformin, **µmol/L**. | PopPK pipeline (steps 01–03) |
| `pbpk_pk_parameters_individual.csv` | **PBPK per-subject BE endpoints** already NCA'd in PK-Sim: `T50, IndividualId, Cmax, AUC_inf, AUC_tEnd` (192 rows = 8×24). | Overlay step 04 (the PBPK comparator) |

Critical structural fact that makes the whole comparison clean: **PK-Sim simulates the same
24 virtual physiologies (IndividualId 0–23) under every formulation.** So test-vs-reference
is a **fully paired 2×2 crossover** with no period/sequence effect — the textbook
average-BE setting. This pairing is exploited identically on both the PBPK and PopPK sides.

The columns the pipeline reads from the dense Results.csv:

- `IndividualId`
- `Time [min]`
- `Organism|PeripheralVenousBlood|Metformin|Plasma (Peripheral Venous Blood) [µmol/l]`

---

## 2. The pipeline at a glance

`programs/poppk_pumas/` runs four scripts in order (each in a fresh Julia process via
`run_all.jl`, so top-level `const`s never clash and any failure halts the chain):

```
data/popsims/*-Results.csv   (dense PBPK truth, 24 subj, 0–24 h, 3-min grid)
        │
   01_build_dataset.jl   ── sparse-sample at BE schedule  +  add 15% proportional residual error
        │                    → output/sparse_be_dataset.csv   (Pumas NONMEM-style)
        ▼
   02_fit_poppk.jl       ── FIT a 1-cmt PopPK model per formulation (FOCE)         ← THE REFIT
        │                    → poppk_param_estimates.csv  (pop params + OFV)
        │                    → poppk_icoef.csv            (per-subject EBEs Ka/CL/Vc)
        │                    → poppk_ipred_profiles.csv   (dense IPRED 0–24 h)       ← THE RESIM
        ▼
   03_poppk_metrics.jl   ── NCA on the dense IPRED → per-subject Cmax/AUC_tEnd/AUC_inf
        │                    → poppk_pk_parameters_individual.csv  (SAME schema as PBPK)
        ▼
   04_be_safespace_overlay.jl ── identical TOST/VBE/safe-space on PBPK *and* PopPK, overlay
                             → be_overlay_gmr.csv, safe_space_overlay_boundaries.csv
                             → fig_overlay_gmr.png, fig_overlay_safespace_*.png
```

The BE/VBE/safe-space statistics (`src/be_utils.jl`) are **byte-for-byte the same helpers**
used on the PBPK side, so the overlay is a true apples-to-apples comparison — the only thing
that differs between the two curves is *how the per-subject Cmax/AUC numbers were produced*
(PK-Sim NCA vs PopPK fit → NCA).

**Baseline design choices, locked 2026-06-10:**

| Choice | Value | Set in |
|---|---|---|
| Sampling schedule | rich BE: 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 16, 24 h (13 obs; predose = dose only) | `01` `RICH_SCHEDULE` |
| Residual (assay) error | proportional, **15% CV** | `01` `PROP_CV` |
| Structural model | 1-compartment, **Weibull-release absorption** (α = in-vivo release scale ~ T50, γ = shape ~ β); first-order+lag retained as sensitivity model | `src/poppk_model.jl` |
| IIV (BSV) | diagonal Ω on **α, CL, Vc** | `src/poppk_model.jl` |
| Reference formulation | **T50 = 3.0 h** | `src/be_utils.jl` `REF_T50` |
| Dose | 390 mg free base = **3019.5 µmol** (≙ 500 mg HCl tablet) | `01` |
| Fit | **FOCE**, each formulation arm fit independently | `02` |
| VBE | 500 trials; n = 48 and 86; endpoint-specific within-subject CV | `src/be_utils.jl` |

---

## 3. Step 01 — build the simulated BE study (`01_build_dataset.jl`)

The dense PBPK profiles are treated as the underlying **"truth."** We turn them into a
**realistic, observable BE study** — what you would actually measure in a clinic — by doing
two things:

### 3.1 Sparse-sampling
Each subject's dense profile is sampled at the **rich BE schedule** (13 nominal times):
`0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 16, 24 h`. Every nominal time lands exactly on the
PK-Sim 3-minute grid (verified by converting hours→minutes and `findfirst`; an absent time
is a hard error). Predose (t=0) is the **dose event only**, not an observation.

A deferred truly-sparse alternative (`SPARSE_SCHEDULE = [1, 3, 6, 12, 24]`, ~6 points) is
coded but not active — it would stress-test the fit.

### 3.2 Adding residual error (the *first* added variability)
To each sampled "true" concentration `c_true`, a **proportional residual error of 15% CV**
is applied:

```julia
c_obs = c_true * (1 + 0.15 * randn())      # proportional (assay / intra-occasion) noise
c_obs = max(c_obs, 1e-6)                    # guard the rare negative draw
```

This is the "simulated analytical measurement error" the master plan (Phase 4.2) and Morris
both insist on: **without it, zero-noise simulated profiles have unrealistically narrow
Cmax/AUC distributions and the type-I-error / power analysis is meaningless** — every trial
would be a clean replicate. Seed = `20260610` for reproducibility.

### 3.3 Output — a Pumas NONMEM-style dataset
`output/sparse_be_dataset.csv`, one dose row + 13 observation rows per subject per
formulation:

| col | meaning |
|---|---|
| `SUBJID` | `T50_<t50>_ID_<id>` (unique per arm×subject) |
| `T50` | formulation label (carried through as a covariate) |
| `IndividualId` | 0–23, the **paired** physiology id (carried through to preserve crossover pairing) |
| `time` | hours |
| `evid` | 1 = dose, 0 = observation |
| `amt` | 3019.5 µmol at t=0 into depot (cmt 1), else 0 |
| `cmt` | 1 = depot (dose), 2 = central (observation) |
| `dv` | observed conc (µmol/L) = the noisy sample |

**Dose detail:** clinical 500 mg metformin **HCl** tablet = 390 mg free base (factor 0.78);
PK-Sim doses in free base; converted mg→µmol with MW 129.16 → **3019.5 µmol**.

---

## 4. Step 02 — the REFIT: fit a PopPK model per formulation (`02_fit_poppk.jl`)

**Yes — the model is re-fit.** This is the core "PBPK-informed PopPK" step. We do **not**
reuse any PBPK parameter. We hand the sparse, noisy BE dataset to Pumas and **estimate a
population PK model from scratch**, exactly as you would with real clinical BE data — that
is the whole point (it asks: *if you only had a clinical-style sparse study, would PopPK MBBE
reach the same BE decision as the underlying PBPK?*).

### 4.1 The structural model (`src/poppk_model.jl`) — WEIBULL-RELEASE (baseline)
A **1-compartment disposition fed by a Weibull in-vivo release**, mirroring the PK-Sim
tablet dissolution directly. Metformin XR is dissolution/release rate-limited (the slow
"flip-flop" tail is set by how slowly the tablet releases), so the Weibull release *is*
the absorption input — no separate first-order `ka` is needed or identifiable. The depot
empties by a Weibull hazard, so `Depot(t) = Dose·exp(−((t−lag)/α)^γ)` reproduces the same
release curve the PBPK side uses; **α ≈ PBPK T50, γ ≈ PBPK β**. Units match the PK-Sim
exports (h, µmol, L, µmol/L).

```julia
@param begin
    tvalpha ∈ RealDomain(lower=0, init=3.0)    # typical Weibull release scale (h) ~ T50
    tvgamma ∈ RealDomain(lower=0.1, init=1.5)  # typical Weibull shape           ~ β
    tvcl    ∈ RealDomain(lower=0, init=11.0)   # typical clearance (L/h)
    tvvc    ∈ RealDomain(lower=0, init=120.0)  # typical central volume (L)
    tvlag   ∈ RealDomain(lower=0, init=0.25)   # release lag (h)
    Ω       ∈ PDiagDomain(init=[0.09,0.09,0.09])# IIV variances: α, CL, Vc
    σ_p     ∈ RealDomain(lower=0, init=0.15)   # proportional residual SD
end
@random η ~ MvNormal(Ω)
@pre begin
    Alpha = tvalpha*exp(η[1]);  CL = tvcl*exp(η[2]);  Vc = tvvc*exp(η[3])
    Gamma = tvgamma;  Tlag = tvlag           # shape & lag: fixed effects
end
@dynamics begin
    # Weibull hazard release h(t) = (γ/α)((t−lag)/α)^(γ−1), tau>0 else 0
    Depot'   = -(t>Tlag ? (Gamma/Alpha)*(max(t-Tlag,1e-8)/Alpha)^(Gamma-1) : 0.0)*Depot
    Central' =  (t>Tlag ? (Gamma/Alpha)*(max(t-Tlag,1e-8)/Alpha)^(Gamma-1) : 0.0)*Depot - (CL/Vc)*Central
end
@derived begin
    cp := @. Central/Vc
    dv ~ @. Normal(cp, abs(cp)*σ_p + 1e-6)   # proportional residual (tiny floor keeps SD>0)
end
```

*(The prior first-order-absorption model — `tvka` with `@dosecontrol lags`, T50 folded
into an apparent `ka` — is retained only as a sensitivity comparison.)*

### 4.2 What BSVs / variabilities are in the model (the *second* added variability)
This is the explicit answer to "what BSVs were added":

- **Between-subject variability (IIV / η, log-normal) on three parameters:** `α`, `CL`,
  `Vc`, via `Alpha = tvalpha·exp(η₁)`, `CL = tvcl·exp(η₂)`, `Vc = tvvc·exp(η₃)`.
- **Ω is diagonal** (`PDiagDomain`) — no off-diagonal (no IIV correlation modelled),
  initialised at variance **0.09 each** (≈30% CV) for α, CL, Vc respectively.
- **Weibull shape γ and lag time are fixed effects — NO IIV**, deliberately, for
  estimation stability from the sparse BE design (add IIV later if the fit supports it).
- **Residual error: proportional, `σ_p`** (init 0.15) — this is the model's estimate of the
  same 15% proportional noise that step 01 injected.

So there are **two distinct, separately-purposed variabilities** in the pipeline:
1. the **15% proportional residual** *added to the data* in step 01 (the assay noise the
   model then has to estimate back via `σ_p`), and
2. the **IIV (Ω) on ka/CL/Vc** *estimated by the model* in step 02 (genuine between-subject
   physiological spread).
A third, separate within-subject CV is injected later for the VBE trials — see §6.2.

### 4.3 Estimation
- **Each of the 8 formulation arms is read into its own `Pumas` population and fit
  independently** with **`FOCE()`** (`init_params(poppk_model)` as starting values). The
  arms are *not* pooled — this mirrors "each formulation characterised independently." (A
  joint model with formulation as a categorical covariate on `ka` is noted as a reasonable
  alternative in the README, not used in the baseline.)
- `read_pumas` maps: `id=:SUBJID, time=:time, amt=:amt, evid=:evid, cmt=:cmt,
  observations=[:dv], covariates=[:T50, :IndividualId]`. Note **T50 and IndividualId are
  carried as covariates** — they are *not* estimated; they are bookkeeping (formulation
  label + crossover-pairing id).

### 4.4 The RE-SIMULATION: dense individual predictions
After each fit we **regenerate dense profiles**:

- `icoef(fpm)` → per-subject **empirical Bayes estimates (EBEs)** of `Ka, CL, Vc`
  → `poppk_icoef.csv`.
- `predict(fpm; obstimes = 0:0.05:24)` → **dense individual predictions (IPRED)** on a
  0.05-h grid spanning **0–24 h to match the PBPK profile window** (same truncation on both
  sides ⇒ fair AUC comparison) → `poppk_ipred_profiles.csv` (reads `pr.ipred.dv`, `pr.time`).

So "did we re-simulate?" — **yes:** the fitted PopPK model is used to reconstruct smooth
0–24 h concentration–time curves for every subject. (This is an individual-prediction
reconstruction from the EBEs, *not* a fresh stochastic simulation of brand-new subjects —
the 24 paired physiologies are preserved so the crossover pairing survives into the BE step.)

> **Robustness note (README):** because the disposition is a textbook 1-cmt first-order
> model, the dense IPRED can be regenerated **analytically** from `poppk_icoef.csv` without
> re-running Pumas:
> `C(t) = (Dose·Ka)/(Vc·(Ka−ke))·(e^(−ke·(t−lag)) − e^(−Ka·(t−lag)))`, `t>lag`,
> `ke = CL/Vc` (with the `Ka≈ke` limit handled). This is a fallback if `icoef`/`predict`
> field names differ across Pumas versions.

**Outputs of step 02:** `poppk_param_estimates.csv` (tvka, tvcl, tvvc, tvlag, the three Ω
diagonals, σ_p, and OFV via `deviance`, wrapped in try/catch → `missing` if unavailable),
`poppk_icoef.csv`, `poppk_ipred_profiles.csv`.

> This is the **only Pumas-license-sensitive step.** It has *not* been executed on the
> authoring host (the open-registry `Pumas` there is an empty stub); it must run on a
> licensed Pumas / JuliaHub session. Steps 01, 03, 04 use only open packages.

---

## 5. Step 03 — regenerate the BE endpoints from the PopPK profiles (`03_poppk_metrics.jl`)

Now the PopPK dense IPRED profiles are turned back into the **same three BE endpoints the
PBPK side reports**, using a minimal, transparent NCA (`src/nca_lite.jl`) so the comparison
is on identical metric definitions:

- **Cmax** — observed peak of the dense profile.
- **AUC_tEnd** (= AUC₀–t, the 0–24 h window) — **linear-up / log-down trapezoid**
  (`auc_trap`): log-trapezoid on falling segments, linear on rising/flat/zero. The 0.05-h
  grid makes this effectively exact.
- **AUC_inf** = AUC_tEnd + Clast/λz, where **λz** is the terminal slope from a log-linear
  regression that scans the last 3–8 post-peak points and keeps the best adjusted R²
  (`lambda_z`) — standard NCA practice.
- Diagnostics also kept: Tmax, % extrapolated, λz adjusted R².

The output `output/poppk_pk_parameters_individual.csv` is written with **exactly the PBPK
schema**: `T50, IndividualId, Cmax, AUC_inf, AUC_tEnd`. That schema match is what lets step
04 run *identical* BE math on both sources. (`poppk_nca_diagnostics.csv` keeps the extra
diagnostic columns.) The script also prints per-formulation **geometric-mean** Cmax /
AUC_tEnd / AUC_inf and mean % extrapolated.

> The AUC convention is deliberately matched to the PK-Sim NCA: AUC_tEnd over 0–24 h,
> AUC_inf with terminal extrapolation, so flip-flop / extrapolation behaviour is comparable
> between the two methods. AUC_0-t is the preferred BE endpoint for XR / high-extrapolation
> products; AUC_inf is reported alongside.

---

## 6. Step 04 — the headline deliverable: PBPK vs PopPK decision-equivalence overlay (`04_be_safespace_overlay.jl`)

Loads **both** per-subject endpoint tables —
`data/pbpk_pk_parameters_individual.csv` (PBPK) and
`output/poppk_pk_parameters_individual.csv` (PopPK) — and runs the **same** TOST / VBE /
safe-space machinery (`src/be_utils.jl`, the exact PBPK-side code) on each, then overlays.

### 6.1 GMR(T50) with 90% CI — both methods
For every T50 arm and every endpoint, a **paired (crossover) average-BE** test vs the
reference (T50 = 3.0 h) using within-subject log-differences (`be_paired`):

- `d = log(test) − log(ref)` over subjects with finite, positive paired values;
- `GMR = exp(mean d)`; 90% CI from a t-distribution (TOST, two one-sided tests at α=0.05);
- **BE pass** if the whole 90% CI ⊂ **[0.80, 1.25]** (FDA/EMA average-BE limits).

→ `be_overlay_gmr.csv` and `fig_overlay_gmr.png` (GMR vs T50, PBPK firebrick vs PopPK
royalblue, dodged so CIs don't overlap, with the 0.8–1.25 band shaded).

### 6.2 Virtual BE (VBE) power / safe-space curves — both methods
For each endpoint, a **power curve over a fine T50 grid** (`anchors[1]..anchors[end]`, step
**0.05 h**), at **two study sizes n = 48 and n = 86**:

- The 8 anchor T50s give each subject a per-formulation metric; values **between anchors are
  log-interpolated** (`interp_log`) so power can be evaluated on the fine grid.
- At each grid T50, `vbe_replicate` runs **500 bootstrap trials**: resample `n_subj`
  subjects **with replacement** from the 24 paired physiologies (n>24 ⇒ bootstrap), inject a
  **within-subject residual CV** (the *third* added variability), run the paired TOST, and
  record the **pass fraction = power**.
- **Joint endpoint** (`vbe_replicate_joint`): a trial passes only if **both** Cmax and
  AUC₀–t pass.

**The within-subject CVs injected for VBE** (`WSV_CV`, from FDA NDA 212595 ClinPharm review,
metformin ER): **Cmax 0.122, AUC_inf 0.114, AUC_tEnd 0.112**. A high-variability sensitivity
set (`WSV_CV_IR`, immediate-release upper bound) is also available: **Cmax 0.237, AUC 0.147**.
The CV is converted to a log-scale SD (`sw = sqrt(log(cv²+1))`) and applied multiplicatively
to both test and reference draws each trial.

**Safe-space boundary** = the contiguous T50 range around the reference where **power ≥ 0.80**
(`safe_bounds` walks outward from the ref T50 while power stays ≥ target).

→ `fig_overlay_safespace_<endpoint>.png` for Cmax, AUC_tEnd, AUC_inf, and the joint endpoint
(power curves PBPK vs PopPK with their ≥80%-power boundaries marked as dotted verticals), and
the boundary comparison table `safe_space_overlay_boundaries.csv`
(`endpoint, source, n_subj, ss_lo, ss_hi`).

### 6.3 What "decision equivalence" means here
The question the overlay answers: **do PBPK and PopPK put the safe-space boundaries in the
same place and reach the same BE conclusions?** Same anchors, same reference, same TOST, same
VBE CVs, same NCA conventions — the *only* difference between the two curves is the source of
the per-subject Cmax/AUC (direct PK-Sim NCA vs fit-then-NCA through a PopPK model). Any
divergence is therefore attributable to the PopPK abstraction (sparse sampling + structural
model + estimation), which is exactly what the study is measuring.

---

## 7. How T50 is handled in the PopPK model vs the PBPK model (your specific question)

This is the subtle and important part, so it gets its own section.

### On the PK-Sim / PBPK side
T50 is a **mechanistic input**: it is the Weibull median dissolution time of the formulation.
You *set* it, PK-Sim dissolves the tablet accordingly, and the physiology turns that release
profile into plasma concentrations. Each T50 value is a separate population simulation. The
BE endpoints (Cmax/AUC) come straight from NCA of those mechanistic profiles. T50 is a
**cause** in the model.

### On the PopPK side
T50 is **never estimated and never appears as a structural parameter.** It is carried purely
as a **label/covariate** (`covariates = [:T50, :IndividualId]`). Its role is bookkeeping:

1. **It defines the arms.** Step 02 fits the 1-cmt model **independently per T50 arm**. So
   the formulation's effect on release is absorbed entirely into that arm's **estimated
   parameters** — chiefly the **apparent absorption rate `ka`** (and to a lesser extent
   lag). The model comment states this explicitly: *the in-vitro release T50 maps inversely
   onto the apparent ka.* A slower-releasing (larger-T50) formulation is reconstructed by the
   PopPK fit as a smaller apparent `ka`. **T50 → ka is implicit, via independent fitting —
   not a fitted regression and not an estimated covariate effect.**
2. **It is the x-axis** of every downstream comparison (GMR vs T50, power vs T50,
   safe-space boundaries in T50 units), so the PopPK results can be laid directly over the
   PBPK results on the same formulation axis.

### So how are the per-formulation (T50) results "regenerated" in the PopPK model vs PK-Sim?
They are regenerated by the **fit → re-simulate → re-NCA chain**, not by re-deriving T50:

| Stage | PK-Sim / PBPK | PopPK (Pumas) |
|---|---|---|
| Formulation effect | T50 set mechanistically (Weibull dissolution) | absorbed into the arm's estimated `ka`/lag (T50 only a label) |
| Per-subject profile | mechanistic PBPK ODE, dense 0–24 h | re-simulated dense IPRED 0–24 h from the fitted EBEs (`Ka, CL, Vc`) |
| Per-subject Cmax/AUC | PK-Sim NCA | `nca_lite` NCA on the IPRED (identical definitions) |
| GMR(T50), power(T50), safe space | TOST/VBE on PBPK endpoints | **same** TOST/VBE on PopPK endpoints |

In other words, the PopPK model regenerates the **same observable BE quantities at each T50**
(the per-subject Cmax/AUC, and from them GMR and the safe-space boundaries), but it does so
by *estimating a PK model from a sparse noisy study and predicting from it*, whereas PK-Sim
produces them mechanistically. The overlay then asks whether those two regeneration routes
agree.

> If you ever want an *explicit* T50→ka relationship (rather than the implicit one from
> independent per-arm fits), the README's deferred variant — a joint model with T50 as a
> categorical (or continuous) covariate on `ka` — would estimate it directly. The baseline
> deliberately keeps the arms independent.

---

## 8. Exactly what was added vs reused (quick-reference for the question "did we refit / what BSVs")

| Question | Answer |
|---|---|
| Did we re-fit a model? | **Yes.** Step 02 fits a 1-cmt first-order-absorption-with-lag PopPK model **per formulation arm** with **FOCE**, from the sparse noisy BE dataset. No PBPK parameter is reused. |
| Did we re-simulate? | **Yes.** Step 02 reconstructs dense 0–24 h individual predictions (IPRED) from the fitted per-subject EBEs (`predict(fpm; obstimes=0:0.05:24)`). |
| What residual error was added to the data? | **Proportional 15% CV** assay/intra-occasion noise, added to each sparse sample in step 01 (seed 20260610). |
| What BSV/IIV is in the PopPK model? | **Log-normal IIV (η) on ka, CL, Vc**, diagonal Ω (init 0.09 each ≈ 30% CV). **No IIV on lag** (fixed effect, for stability). Residual estimated as proportional `σ_p` (init 0.15). |
| What within-subject CV is added for VBE? | **Cmax 0.122, AUC_tEnd 0.112, AUC_inf 0.114** (FDA NDA 212595 metformin ER); IR sensitivity set 0.237 / 0.147 available. Injected as multiplicative log-normal noise in the 500-trial bootstrap. |
| How is T50 regenerated in the PK model vs PBPK? | T50 is **not** regenerated/estimated in PopPK — it is a label. The formulation effect is absorbed into each arm's estimated **apparent ka**. Per-T50 BE results are regenerated via fit→IPRED→NCA→TOST and overlaid on the mechanistic PBPK results. See §7. |
| Are the BE statistics the same on both sides? | **Yes — identical code** (`src/be_utils.jl`): same TOST, same 0.80–1.25 limits, same 90% CI, same VBE (500 trials, n=48/86), same within-subject CVs, reference T50 = 3.0 h. |

---

## 9. Output inventory (`programs/poppk_pumas/output/`)

| File | Produced by | Contents |
|---|---|---|
| `sparse_be_dataset.csv` | 01 | simulated BE study (Pumas NONMEM format) |
| `poppk_param_estimates.csv` | 02 | population params (tvka/tvcl/tvvc/tvlag, Ω diagonals, σ_p) + OFV per arm |
| `poppk_icoef.csv` | 02 | per-subject EBEs Ka/CL/Vc per arm |
| `poppk_ipred_profiles.csv` | 02 | dense individual predictions, 0–24 h |
| `poppk_pk_parameters_individual.csv` | 03 | PopPK per-subject Cmax/AUC_inf/AUC_tEnd (**PBPK schema**) |
| `poppk_nca_diagnostics.csv` | 03 | + Tmax, % extrapolated, λz adj-R² |
| `be_overlay_gmr.csv` / `fig_overlay_gmr.png` | 04 | GMR(T50) ± 90% CI, both methods |
| `safe_space_overlay_boundaries.csv` / `fig_overlay_safespace_*.png` | 04 | safe-space boundary table + overlay figures |

## 10. How to run

```bash
# from programs/poppk_pumas/ on a licensed Pumas / JuliaHub environment
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # first time only
julia --project=. run_all.jl                            # runs 01 → 04
# or step by step: scripts/01_build_dataset.jl … 04_be_safespace_overlay.jl
```

Only step 02 needs the licensed Pumas; 01/03/04 use open packages. If `icoef`/`predict`
field names differ on your Pumas version, adjust step 02 (see README "If something breaks"),
or use the closed-form 1-cmt fallback (§4.4) to regenerate the profiles analytically.

---

## 11. Deferred variants (coded hooks, not in the baseline)

1. **Truly-sparse schedule (~6 pts):** set `SCHEDULE = SPARSE_SCHEDULE` in `01` — stresses
   the fit, demonstrates robustness under limited sampling.
2. **First-order absorption + lag (sensitivity model):** the earlier `@dynamics` form with
   `tvka` and `@dosecontrol lags`, where T50 is folded into an apparent `ka`. Kept as a
   comparison against the Weibull baseline.
3. **Joint T50-covariate model:** fit all arms together with T50 as a covariate on `α` to
   get an explicit T50→release relationship instead of independent per-arm fits.
4. **Pumas-free run:** `scripts/02b_fit_weibull_nopumas.jl` fits the same Weibull model per
   subject with open packages (no license) and writes the same outputs — a cross-check of
   the Pumas NLME fit.

---

*Sources: `programs/poppk_pumas/{README.md, run_all.jl, src/*.jl, scripts/*.jl}`,
project notes (not included), and the data files under
`programs/poppk_pumas/data/`. Authored 2026-06-25.*
