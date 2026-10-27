# Metformin XR virtual bioequivalence — methodology in plain language

**What this project does, in one sentence:** we use two different **modelling methodologies**
for metformin extended-release (XR) tablets to answer the same regulatory question —
*"how much can the tablet's release rate drift before it would fail a bioequivalence
(BE) study?"* — and we check whether the two methodologies give the **same answer**.

The two methodologies are:

1. **PBPK** — physiology-based modelling, done in **PK-Sim**, and
2. **PopPK** — population-PK modelling, done in **Pumas**.

The single formulation "knob" we turn is **T50** — the Weibull *release time* of the
tablet (hours). Bigger T50 = slower release. Everything below explains, in simple
terms, exactly what was done on each side and how they are compared.

> **A note on the plan.** The intended PopPK absorption model is the **Weibull-release
> model** described in Part 2 — the same release law the PBPK tablet uses. This is the
> baseline throughout; a simpler first-order-absorption form is retained only as a
> secondary sensitivity model.

---

## 0. The concepts you need first (30-second versions)

- **Bioequivalence (BE).** Two products are "bioequivalent" if the body is exposed to
  the drug almost identically. Practically: the **90% confidence interval** of the
  test/reference ratio (the **GMR**) for peak concentration (**Cmax**) and total
  exposure (**AUC**) must sit entirely inside **0.80–1.25**. The test is the
  **two one-sided tests (TOST)** procedure.
- **Virtual / model-based BE (VBE / MBBE).** Run that BE study on a *simulated*
  population instead of in real people.
- **Safe space.** The range of the formulation knob (here T50) where the product not
  only *passes* BE but a realistic-sized study has **≥ 80% power** to pass. The safe
  space is narrower than the raw BE range, because power falls off before the ratio
  actually leaves 0.80–1.25.
- **Why two methodologies?** PBPK and PopPK are the two standard MBBE routes. If they place the
  safe-space boundaries in the same place, that is strong evidence the boundaries are
  real and not an artifact of one modelling choice. That agreement is our headline
  deliverable — **"decision equivalence."**

We are reproducing the *method* of **Morris et al., AAPS J 2026** (a long-acting
injectable safe-space + type-I-error paper), but applied to an **oral MR** drug.

---

## 1. The PBPK side (PK-Sim) — "build the physiology, dissolve the tablet"

**Goal:** a mechanistic model where the drug moves through actual organs, and the only
thing we change between formulations is how fast the tablet releases the drug.

All modelling and simulation on this side is done entirely in **PK-Sim** (the small-molecule
PBPK builder). No other tool is used — no MoBi, no scripting layer — everything is built and
run in the PK-Sim GUI.

### 1.1 What the PBPK methodology is (and what PK-Sim gives us)
Unlike a compartmental model (a few abstract boxes fit to data), a PBPK model represents the
body as its **actual organs** — gut, liver, kidney, muscle, fat, brain, etc. — each with a
real blood flow, volume, and tissue-to-plasma partitioning. Drug-specific properties
(molecular weight, lipophilicity, solubility, permeability, protein binding, and
transporter/enzyme kinetics) drive how the drug moves through that physiology. PK-Sim holds
the physiology database and assembles the model from a small set of building blocks:
**Compound** (the drug), **Formulation** (how it is released), **Administration Protocol**
(dose/route), **Individual** or **Population** (the physiology), and **Simulation** (which
ties them together).

### 1.2 The starting model — a published, validated metformin model
We did **not** build metformin from scratch. We loaded the **published open-source metformin
PBPK model** (a PK-Sim project snapshot) and ran it in PK-Sim:

- **Model source:** https://github.com/Open-Systems-Pharmacology/Metformin-Model
- **Model paper:** Hanke et al., *Clin Pharmacokinet* 2020 — https://doi.org/10.1007/s40262-020-00896-w
- It was validated against **39 clinical studies** (0.001–2550 mg, oral + IV). We confirmed
  the install by reproducing the model's own **immediate-release validation** (the
  Tucker 1981 clinical overlay) before changing anything. The model ships covering **oral
  immediate-release and IV only** — the disposition (distribution + clearance) backbone
  comes for free; the **extended-release layer is what we add** (§1.4).

### 1.3 The drug — how metformin behaves in the body
The concrete metformin properties that matter for this study (all carried in the PK-Sim
Compound):

- **Small, hydrophilic base**, MW ≈ 129.16 g/mol (free base); high aqueous solubility but
  **low passive permeability** → a BCS class III drug (dissolves easily, crosses membranes
  poorly).
- **Absorption is transporter-mediated and occurs high in the gut** — there is an effective
  **absorption window**: drug that has not been released before it passes the upper small
  intestine is absorbed less well. This matters for slow-release tablets (see limitations).
- **Cleared almost entirely by the kidney** (renal secretion, essentially no hepatic
  metabolism), via organic-cation transporters. The model carries the canonical Hanke
  transporter set — **OCT1, OCT2, MATE1, PMAT** (renal/hepatic/intestinal cation
  transport) — which is what makes metformin's disposition mechanistic rather than a
  fitted clearance.

Because clearance is renal and transporter-mediated and does **not** depend on how the
tablet releases, **clearance and volume are the same for every formulation** — a fact we
rely on later when interpreting the PopPK fit (§2.3).

### 1.4 Add the extended-release (XR) tablet as a Weibull dissolution
Metformin XR (e.g. Glucophage XR) is a slow-release matrix tablet. In PK-Sim we attach a
**Weibull dissolution Formulation** to the compound — the fraction of drug released by
time *t*:

```
F(t) = 1 − exp( −((t − lag)/T50)^β )
```

- **T50** = the Weibull **scale** = the release time. **This is the single formulation knob
  (the CQA) we vary.**
- **β** = the Weibull **shape** = how S-shaped the release is. **We hold β fixed at 1.0 for
  every tablet** — so across the study only the *scale* changes, not the shape.
- Reference T50 ≈ **3 h** (representative of Glucophage XR in-vitro release; cf. Timmins
  et al., *Clin Pharmacokinet* 2005; Stepensky et al., *Biopharm Drug Dispos* 2002).

**Dose convention (used consistently for every simulation):** a clinical 500 mg metformin
**HCl** tablet = **390 mg free base** (factor 0.78; MW base 129.16 vs salt 165.62). PK-Sim
doses in free base, single dose, **fasted**, 24 h simulation. Keeping this identical across
all tablets is what makes the test-vs-reference comparison valid.

### 1.5 The formulation sweep
We built **8 tablets** by changing only T50: **1.0, 2.0, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0 h**
(3.0 h = reference). β = 1.0, dose, population and everything else are held fixed — so any
difference in the blood profile is attributable to release rate alone.

### 1.6 The virtual population
For each of the 8 tablets we ran a **population simulation of 24 virtual subjects** built
from PK-Sim's population module:

- Template: **European healthy adults**, ages **18–55**, BMI **18.5–30**, **50/50 male/female**.
- PK-Sim varies each subject's physiology (organ sizes, blood flows, transporter expression)
  in an internally consistent way, so the 24 subjects are genuinely different bodies.
- **Crucially, PK-Sim uses the same 24 physiologies under every tablet.** So test-vs-reference
  is a **fully paired 2×2 crossover** (the same person effectively "takes" both tablets) —
  the cleanest possible average-BE setting, with no between-subject noise from mismatched
  arms.

Output per tablet: dense **peripheral-venous plasma** concentration–time curves (a 3-minute
grid, 0–24 h) for all 24 subjects, exported to CSV — the input to everything downstream.

### 1.7 Known limitations / scope (stated honestly)
- **This is a methodology platform, not an absolute-accuracy exercise.** The population
  over-predicts absolute metformin exposure by a few-fold vs clinic (reference-formulation
  geometric-mean Cmax ≈ 2.8 µg/mL and AUC₀–24 ≈ 31 µg·h/mL, vs roughly 0.6–1 µg/mL and
  ~7–11 µg·h/mL clinically — about 3–4×). That is deliberately accepted: bioequivalence is a
  *relative* test/reference comparison, and because every formulation uses the same
  physiologies and dose the over-prediction is a **shared multiplicative offset that cancels
  exactly** in the GMR/TOST — it does not move the safe-space boundaries. (It would matter
  only for an *absolute* exposure/dose claim, which this study does not make.)
- **Slow tablets show flip-flop kinetics.** When release is slower than elimination, the
  terminal decline reflects the *release* rate, not true elimination. This inflates the
  extrapolated **AUC_inf** for the slowest tablets, so **AUC₀–t** is the primary AUC endpoint
  (AUC_inf is reported alongside but treated cautiously).
- **The absorption window caps the slowest tablets.** With very slow release, some drug
  passes the absorption window before dissolving (fraction absorbed falls from ~0.8 toward
  ~0.5), which is a real mechanistic feature of metformin, not a modelling error.

### 1.8 BE endpoints straight from the PBPK profiles
For each subject we compute, by non-compartmental analysis (NCA):
- **Cmax** (peak),
- **AUC₀–t** (exposure over 0–24 h; the preferred endpoint for XR), and
- **AUC_inf** (extrapolated to infinity).

These per-subject numbers are the PBPK "truth" we later compare against.

*Software:* **PK-Sim** — the small-molecule PBPK tool; docs
https://docs.open-systems-pharmacology.org/working-with-pk-sim/pk-sim-tutorial. The
metformin model was loaded as a PK-Sim project snapshot from the source in §1.2.

---

## 2. The PopPK side (Pumas) — "pretend you only have a clinical study, and fit a model"

**Goal:** answer the honest real-world question: *if all you had was an ordinary,
sparse clinical BE study, would the PopPK methodology reach the same BE decision the
mechanistic PBPK methodology does?*

So we deliberately **throw away** the PBPK detail and rebuild a PK model from a
realistic, noisy, sparsely-sampled dataset.

### 2.1 Turn the PBPK "truth" into a realistic clinical study
1. **Sparse-sample** each subject's dense PBPK curve at a normal clinical BE schedule:
   0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 16, 24 h (13 samples).
2. **Add measurement noise:** a **15% proportional** (assay) error on each sample.
   Without noise the simulated study is unrealistically clean and the power/type-I
   analysis is meaningless — the same point Morris and the FDA guidance stress.
3. Write it out in standard NONMEM/Pumas format (dose row + observation rows per
   subject), keeping the T50 label and the subject id (to preserve crossover pairing).

### 2.2 The PopPK model — a **Weibull-release absorption** model (the original plan)
This is the key modelling choice. On the PBPK side the tablet dissolves by a Weibull
law; **we use that same Weibull law as the absorption input on the PopPK side**, so the
two methodologies speak the same mechanistic language.

Metformin XR is **dissolution/release rate-limited** (its slow tail — "flip-flop"
kinetics — is set by how slowly the tablet releases, not by elimination). So we let the
Weibull release *be* the absorption, with no separate absorption-rate constant (there is
nothing extra to identify). The model is a **1-compartment disposition fed by a Weibull
release**:

```
release rate constant   h(t) = (γ/α) · ((t − lag)/α)^(γ − 1)      (for t > lag)
depot (undissolved drug) Depot'   = −h(t)·Depot,   Depot(0) = Dose
central (in body)        Central' =  h(t)·Depot − (CL/Vc)·Central
concentration            C(t)     =  Central / Vc
```

This makes `Depot(t) = Dose · exp(−((t−lag)/α)^γ)` — **exactly the Weibull release curve
the PBPK side uses.** So:

- **α** (the fitted in-vivo release scale) is the direct analog of the PBPK **T50** — it
  varies with formulation;
- **γ** is the in-vivo **input shape** — the in-vitro dissolution (PBPK β = 1.0) folded
  together with gut absorption and transit; it is estimated, *not* fixed to β = 1 (in the
  per-arm baseline it comes out ≈ 1.7; in the disposition-constrained IVIVC variant it is a
  single shared value ≈ 0.8 — see §4). It does not vary by formulation (neither does the
  shape on the PBPK side);
- **CL, Vc** are clearance and volume (disposition).

Because it uses the same release law, the PopPK model *can reproduce the PBPK profile
shape* at each formulation — which the older first-order model could not.

**On the 1-compartment disposition.** We use one central compartment. In the published
metformin population-PK literature both structures appear: most models are **two-compartment**
(metformin shows a mild biexponential decline and red-blood-cell uptake), but a **one-compartment
model with first-order absorption and lag is documented as adequate for immediate-release
metformin** (e.g. the review by Wang et al. 2022, and the Korean FDC analysis, Chung/Kim et al.
2018). For this *relative* BE comparison the one-compartment disposition is a reasonable,
transparent choice; a two-compartment disposition is the natural sensitivity/extension if the
terminal phase needs it.

**Variability added in the model:**
- **Between-subject variability (IIV):** log-normal random effects on **α, CL, Vc**
  (diagonal Ω). γ and lag are kept as fixed effects for stability with sparse data.
- **Residual error:** proportional (the model's estimate of the 15% assay noise).

### 2.3 Fit the model — this is a real re-fit
Each formulation arm is fit by population maximum-likelihood estimation. We reuse **no**
PBPK parameter — the model is estimated from the sparse noisy data alone, exactly as it
would be with a real clinical study.

**One structural caveat — release vs disposition (flip-flop).** Metformin XR is
release-rate-limited, so its terminal decline is set by the slower of two processes
(release vs elimination). In a one-compartment model these two are **not separately
identifiable from a single arm** — "slow release + fast elimination" and "fast release +
slow elimination" produce the same curve. If every arm is allowed its own clearance and
volume, the fit can absorb the formulation effect into disposition and leave the release
scale α nearly flat. The physiologically correct constraint removes the ambiguity:
**clearance (CL) and volume (Vc) are properties of the drug and subject, not of the tablet's
release rate, so they are held formulation-independent** (shared across arms). Only the
release **scale α** then varies by formulation (γ stays shared), and α becomes an
identifiable, monotonic in-vivo analog of T50 (an in-silico IVIVC). Constraining
disposition also smooths the
endpoint-vs-formulation relationship, so it can shift the safe-space boundaries relative to
an unconstrained per-arm fit — the release-vs-disposition split and the safe space are
coupled. The safe-space comparison should therefore be reported from a single, consistent
fit, stating which one.

### 2.4 Re-simulate and re-derive the same endpoints
From the fit we reconstruct each subject's **dense predicted profile (0–24 h)**, then run
the **same NCA** as the PBPK side to get per-subject **Cmax / AUC₀–t / AUC_inf** — same
column layout, same definitions. Now the two methods produce endpoints that can be
compared apples-to-apples.

---

## 3. Comparing the two — the safe space and "decision equivalence"

Both sides now have per-subject Cmax/AUC at each T50. We run the **identical** statistics
code on each:

1. **GMR(T50) with 90% CI** — paired TOST vs the reference (T50 = 3 h). BE passes if the
   whole 90% CI is inside 0.80–1.25.
2. **Virtual BE power curves** — for a fine grid of T50 values and two realistic study
   sizes (**n = 48** and **n = 86**), resample crossover trials (bootstrap of the 24
   paired subjects), inject the **within-subject variability** measured for metformin ER
   (Cmax CV ≈ 12.2%, AUC ≈ 11.2%; source: FDA NDA 212595 clinical-pharmacology review),
   run TOST, and record the fraction of trials that pass = **power**.
3. **Safe space** = the T50 range around the reference where power ≥ 80%.
4. A **joint** endpoint (must pass both Cmax and AUC₀–t) is also computed.
5. **Simulated BE boundary + Type I error** — invert GMR(T50) to find the T50 where the
   GMR is exactly 0.80 / 1.25, then run replicate trials *at* that T50; the acceptance
   rate there is the Type I error (nominal 5%). See §4.2.

Because the *only* difference between the PBPK and PopPK curves is **how the per-subject
Cmax/AUC were produced** (mechanistic NCA vs fit-then-NCA), any gap between them is
exactly the effect of the PopPK methodology (sparse fit-then-NCA) instead of the full
physiology. That gap is what "decision equivalence" measures.

---

## 4. What we found

**Decision equivalence (the primary result) — from the licensed Pumas fit.** The PopPK
model was fit to each formulation arm and its safe-space boundaries put next to the PBPK
truth. They **nearly coincide**: the sparse-data PopPK methodology reaches the same BE
conclusions as the mechanistic PBPK methodology.

**Safe-space boundaries (power ≥ 80%), PBPK vs PopPK — licensed Pumas run:**

| Endpoint | n | PBPK [lo, hi] h | PopPK [lo, hi] h |
|---|---|---|---|
| Cmax | 48 | [2.00, 4.15] | [1.95, 4.05] |
| Cmax | 86 | [1.80, 4.25] | [1.75, 4.15] |
| AUC₀–t | 48 | [1.00, 4.60] | [1.00, 4.55] |
| AUC₀–t | 86 | [1.00, 4.70] | [1.00, 4.65] |
| AUC_inf | 48 | [1.00, 4.95] | [1.00, 4.95] |
| **Joint (Cmax & AUC₀–t)** | 48 | **[1.95, 4.10]** | **[1.90, 4.05]** |
| **Joint (Cmax & AUC₀–t)** | 86 | **[1.85, 4.25]** | **[1.75, 4.15]** |

**Reading it:** Cmax is the binding endpoint (AUC stays equivalent over a much wider T50
range). The joint safe space is **T50 ≈ 1.9–4.1 h**, and the PBPK and PopPK methods agree to
within ~0.05–0.1 h — i.e. **decision equivalence holds**.

The safe space above is confirmed by a **second, independent Pumas fit** — a joint fit over
all arms with shared typical disposition (Cmax n=48 [1.90, 4.10], Joint [1.85, 4.10]). So
decision equivalence is robust to the fitting strategy. (Both runs side by side:
`PERARM_vs_JOINT_COMPARISON.md`.)

**Visual confirmation:** the mean concentration–time profiles from the two methodologies
overlap closely at every formulation — see `fig_pbpk_vs_poppk_profiles.png` (PBPK solid vs
PopPK dashed, one panel per T50; peak concentrations agree to within ~3%). The only visible
gap is a slight PopPK under-peak / heavier tail at the slowest tablets, the expected
signature of fitting a one-compartment model to flip-flop profiles.

### 4.1 In-vivo release scale (IVIVC) — the honest finding: it is not identifiable here
Ideally the fitted release scale α would recover each tablet's T50 (an in-silico IVIVC). It
does **not** — and this held across **two** licensed Pumas fits:

| PBPK T50 (h) | 1.0 | 2.0 | 3.0 | 3.5 | 4.0 | 5.0 | 6.0 | 7.0 |
|---|---|---|---|---|---|---|---|---|
| α, per-arm fit | 1.42 | 1.60 | 1.69 | 1.65 | 1.68 | 1.65 | 1.69 | 1.71 |
| α, joint fit (fitted IVIVC exponent ≈ 0.10) | 1.43 | 1.57 | 1.68 | 1.67 | 1.72 | 1.71 | 1.70 | 1.76 |

α is essentially **flat** in both. This is the **one-compartment flip-flop
non-identifiability**: from a sparse single read-out the model cannot separate release rate
from elimination rate, so the formulation's effect is absorbed into **apparent disposition**,
not the release scale. The joint fit was built specifically to prevent this by sharing the
*typical* disposition — but the per-subject clearance and volume still carry random effects
that drift with formulation (fitted **CL 10.6→14.7 L/h, Vc 86→190 L across T50 = 1→7 h**),
so the leak persists and the IVIVC exponent comes out ≈ 0.

**Interpretation (and it is a useful one):** PopPK MBBE reproduces the BE *decision* but not
the *mechanism* (the in-vitro→in-vivo release mapping) from sparse data. Recovering the
release mapping requires a **mechanistic model (PBPK)** — which is exactly where PBPK adds
value over PopPK. (An earlier license-free cross-check appeared to recover a monotonic α only
because it forced identical disposition on every arm of a subject — i.e. it imposed the
answer; the licensed NLME fits above are authoritative and supersede it.)

*Note on γ (release shape):* the licensed joint fit estimated a single shared γ ≈ 1.7. A
license-free exploration of fixing γ = 1.0 (the PBPK dissolution β) vs estimating it is in
`GAMMA_COMPARISON_SUMMARY.md`, but those numbers are from the over-constrained cross-check and
are illustrative only — the licensed fits estimate γ and are authoritative.

### 4.2 Type I error at the simulated BE boundary

**What this adds.** The safe space above answers *"where does power fall below 80%?"*.
Morris Part 2's other question is *"where exactly is the BE boundary, and what is the Type I
error there?"* — and it is the sharper question, because Morris flags that an 8% error in the
simulated limits produces a ~400% error in the Type I error. We ran it on the same data, with
the same TOST/VBE machinery (`BEUtils`), in three steps per method and endpoint:

1. **Invert** the CQA → GMR relationship: walk outward from the reference on a 0.005 h T50
   grid and find where the deterministic paired GMR first leaves [0.80, 1.25]. That T50 is the
   *simulated BE boundary*.
2. **Type I error** = the Monte-Carlo acceptance rate at that boundary, where the true GMR sits
   exactly on the limit. Replicates were raised from the 500 used for the power curves to
   **20 000** — at 500 the standard error is ±1.0 percentage point, which cannot resolve a 5%
   tail.
3. **Cross-evaluate**: the PBPK truth's acceptance rate at the boundary *PopPK* declared. This
   is the number a sponsor actually cares about.

**Step 1 — the simulated BE boundaries (T50, h):**

| Endpoint | Side | Limit | PBPK | PopPK per-arm | PopPK joint |
|---|---|---|---|---|---|
| Cmax | lower | 1.25 | 1.457 | 1.359 | 1.250 |
| Cmax | upper | 0.80 | 4.623 | 4.574 | 4.708 |
| AUC₀–t | lower | 1.25 | *not reached* | *not reached* | *not reached* |
| AUC₀–t | upper | 0.80 | 5.184 | 5.221 | 5.526 |
| AUC_inf | upper | — | 5.357 **(hits 1.25)** | 6.158 (hits 0.80) | 6.546 (hits 0.80) |
| **Joint** | both | Cmax binds | **1.457 / 4.623** | **1.359 / 4.574** | **1.250 / 4.708** |

"Not reached" means the GMR never leaves the acceptance corridor within T50 ∈ [1, 7] h — the
boundary is genuinely not invertible there, not merely uncertain. Consistent with §4, **Cmax is
the only endpoint that binds on both sides**, so it also sets the joint boundary.

**Step 2 — both methods return 5%, and that is arithmetic, not validation.** At its own
boundary each method gives 4.8–5.4% across every endpoint and both study sizes; the
scaled-reference control (reference × 0.80) gives 5.15% at n=48 and 5.39% at n=86. This is
guaranteed by TOST: at a true GMR of exactly 0.80 you pass only if the estimate lands ~1.68
standard errors above the truth, a 1-in-20 event *regardless of the variance or the study size*.
So the agreement here is **not** evidence that either model is right — both would return 5% even
if one were badly wrong. **The like-for-like Type I error comparison is necessarily a null
result**, in the same way §4's decision equivalence is a near-null result.

The one exception is **PBPK AUC_inf at 7.65% (n=48) / 6.69% (n=86)** — above nominal because
flip-flop kinetics make its GMR non-monotonic, so the inversion is unstable. Independent support
for the §1.7 decision to treat AUC_inf cautiously.

**Step 3 — the divergence is entirely in boundary placement.** Plotted against *true GMR* the
two operating characteristics lie on top of one another. What differs is which T50 each method
believes sits at the limit — and because the acceptance curve is near-vertical there, a small
placement error is expensive. "T1 reported" is what PopPK computes at its own boundary; "T1
true" is the PBPK truth's acceptance rate at that same declared T50:

| Fit | Endpoint | Side | n | bnd PBPK | bnd PopPK | CQA err | true GMR | T1 reported | T1 true | ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| per-arm | Cmax | upper | 48 | 4.623 | 4.574 | −1.1% | 0.805 | 4.9% | **8.1%** | 1.66× |
| per-arm | Cmax | upper | 86 | 4.623 | 4.574 | −1.1% | 0.805 | 5.2% | **9.6%** | 1.85× |
| per-arm | Cmax | lower | 48 | 1.457 | 1.359 | −6.8% | 1.269 | 5.0% | 1.4% | 0.28× |
| per-arm | Cmax | lower | 86 | 1.457 | 1.359 | −6.8% | 1.269 | 5.1% | 0.7% | 0.14× |
| per-arm | AUC₀–t | upper | 48 | 5.184 | 5.221 | +0.7% | 0.797 | 5.0% | 3.6% | 0.72× |
| per-arm | Joint | upper | 86 | 4.623 | 4.574 | −1.1% | 0.805 | 4.9% | **8.9%** | 1.81× |
| joint | Cmax | upper | 48 | 4.623 | 4.708 | +1.8% | 0.792 | 4.9% | 2.1% | 0.44× |
| joint | Cmax | upper | 86 | 4.623 | 4.708 | +1.8% | 0.792 | 5.2% | 1.6% | 0.30× |
| joint | Cmax | lower | 48 | 1.457 | 1.250 | −14.2% | 1.293 | 5.0% | 0.2% | 0.03× |
| joint | AUC₀–t | upper | 48 | 5.184 | 5.526 | +6.6% | 0.774 | 4.9% | 0.2% | 0.03× |

**Reading it:** a **1.1% error in the declared T50 moves the true Type I error from 5% to 9.6%**
— an 87% relative error from a one-percent CQA error. Morris's warning reproduced, on an oral MR
product. Note that the 8.1 / 9.6 are *percentages*; the 1.66× / 1.85× are the ratios against
nominal.

**What the risk is, and is not.** At GMR 0.805 the product is genuinely bioequivalent, so nothing
unsafe is waved through. What is wrong is the **reported operating characteristic**: a submission
would state "Type I error at my boundary is 5%" for a limit that really carries 8–10%. On the
*fast-release* side the error runs the other way and is the more consequential one — both fits put
PopPK's lower boundary below PBPK's (1.359 h and 1.250 h against 1.457 h, i.e. 7% and 14% too
permissive), into a region where the true GMR is 1.27–1.29 and the product is genuinely **not**
bioequivalent. The practical consequence is a failed study rather than an unsafe product, since
true acceptance there is only 0.2–1.4%.

**A larger study is more fragile, not less.** This is counter-intuitive and worth stating
explicitly: going from n = 48 to n = 86 *raises* the cost of the same boundary error (8.1% →
9.6% on per-arm Cmax). The acceptance curve is **pinned at 5% at the limit for every study size**
— more subjects cannot help you there — and what more subjects do is **steepen** it, because the
CI narrows. A dense scan through the limit (`08_export_zoom_oc.jl`, 20 000 trials × 5 seeds per
point) shows the two curves crossing exactly at the limit:

| T50 (h) | true GMR | where | n = 48 | n = 86 | effect of larger n |
|---|---|---|---|---|---|
| 4.550 | 0.8074 | inside | 10.0% | 12.4% | higher |
| 4.590 | 0.8033 | inside | 6.9% | 7.7% | higher |
| 4.610 | 0.8013 | inside | 5.7% | 5.9% | higher |
| **4.623** | **0.8000** | **on the limit** | **5.0%** | **4.9%** | **none — pinned** |
| 4.630 | 0.7993 | outside | 4.6% | 4.5% | lower |
| 4.650 | 0.7973 | outside | 3.8% | 3.4% | lower |
| 4.710 | 0.7915 | outside | 1.9% | 1.3% | lower |

The n = 48 / n = 86 ordering **flips exactly across the limit**. So sample size cannot buy
robustness against a misplaced limit — which means the accuracy of the CQA → GMR inversion
matters *more* as studies get larger, not less. **That is the strongest argument in this project
for letting the mechanistic model place the boundary even though the sparse-data model reproduces
the BE decision correctly.**

**Limitations of this sub-analysis, stated honestly:**

- **The direction of the error does not replicate across fitting strategies.** On the upper side
  the per-arm fit places PopPK's boundary *inside* PBPK's (inflation, 1.66×) while the joint fit
  places it *outside* (deflation, 0.44×). Only the **magnitude** of the boundary error survives a
  change of fit; report magnitude and treat the sign as fit-dependent.
- **The finding reduces to the boundary error.** Because the operating characteristics coincide
  in GMR space, "T1 true" is fully determined by the true GMR at the declared boundary. The
  amplification is real and is the right unit to report to a regulator, but it is not evidence
  independent of the boundary discrepancy itself.
- **AUC_inf is not comparable at all.** PBPK's boundary hits the *upper* limit while both PopPK
  fits hit the *lower* one, so those rows are flagged `comparable = false` in the output and must
  not be reported as a Type I error ratio.
- **PBPK remains the truth by construction** (as throughout §2–4), and the within-subject CVs are
  healthy-volunteer values from FDA NDA 212595. Every acceptance rate on this page scales with
  them; they are not estimated here.

*Software / reproduction:* `scripts/05_type1_overlay.jl` (inversion, Type I error,
cross-evaluation) → `scripts/06_type1_summary.jl` (reporting table) → `scripts/07_export_type1_json.jl`
and `scripts/08_export_zoom_oc.jl` (curve exports). Outputs: `type1_overlay_boundaries.csv`,
`type1_overlay.csv`, `type1_calibration_control.csv`, `type1_summary_both_fits.csv/.md`,
`type1_zoom_oc.csv`, `fig_overlay_type1_*.png`. Run in both export folders (per-arm and joint).
A plain-language walkthrough of what happens at the boundary — written for review rather than
for the record — is in `TYPE1_ERROR_EXPLAINED.md`, and the same results as a standalone
read-only page are in `TYPE1_ERROR_RESULTS.html`.

---

## 5. References and links

**Models & software**

| What | Link |
|---|---|
| OSP metformin PBPK model (the base model) | https://github.com/Open-Systems-Pharmacology/Metformin-Model |
| Hanke et al. 2020 — metformin PBPK model paper | https://doi.org/10.1007/s40262-020-00896-w |
| PK-Sim / OSP Suite | https://www.open-systems-pharmacology.org/ · docs https://docs.open-systems-pharmacology.org/ |
| Pumas (PopPK fitting) | https://pumas.ai · docs https://docs.pumas.ai/ |

**Methodology papers**

| What | Link |
|---|---|
| Morris et al. AAPS J 2026 — Part 1 (model development) | https://doi.org/10.1208/s12248-026-01213-2 |
| Morris et al. AAPS J 2026 — Part 2 (type-I error + safe space) | https://doi.org/10.1208/s12248-026-01215-0 |
| Bois & Brochot 2024 — Bayesian VBE | https://doi.org/10.3389/fphar.2024.1404619 |
| Chen et al. 2024 — model-integrated BE | https://doi.org/10.1002/psp4.13225 |
| Gajjar et al. 2022 — popPK VBE for LAIs | https://doi.org/10.1016/j.ejps.2022.106296 |

**Clinical / formulation data for metformin XR**

- Timmins et al., *Clin Pharmacokinet* 2005 — Glucophage XR PK.
- Stepensky et al., *Biopharm Drug Dispos* 2002 — dissolution + PK.
- FDA NDA 212595 clinical-pharmacology review — within-subject CVs used in the VBE.

**Regulatory**

| What | Link |
|---|---|
| FDA PBPK biopharmaceutics guidance (2020) | https://www.fda.gov/media/142500/download |
| FDA Statistical Approaches to BE (2001) | https://www.fda.gov/media/70958/download |

---

*Companion documents:* the full technical walkthrough of the PopPK pipeline is in
`POPPK_ANALYSIS_WALKTHROUGH.md`; the master project plan (all 9 phases, PK-Sim build
through write-up) is in project notes (not included); a plain-language
walkthrough of what happens at the BE boundary is in `TYPE1_ERROR_EXPLAINED.md`.
Authored 2026-07-22; §4.2 (Type I error) added 2026-09-09.*
