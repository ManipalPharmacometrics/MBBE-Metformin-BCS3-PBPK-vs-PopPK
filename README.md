# Metformin XR: PopPK QC

This branch holds only the PopPK part of the metformin MBBE work, rewritten as six short steps
for QC. It runs on the corrected PK-Sim simulations of 30 Oct 2026. The PK-Sim model, decks and
papers stay on `main`.

**Question.** Do a PopPK model and the PBPK model (PK-Sim) put the bioequivalence (BE) safe space
in the same place?

**Approach.** PK-Sim simulated 8 XR formulations that differ only in the in vivo release time
T50 (1–7 h; the reference is 3.0 h), each in the same 24 virtual subjects. We sample those
profiles like a BE study, fit PopPK, then run identical BE code on the PK-Sim endpoints and on the
PopPK endpoints.

## Files: read in this order

| # | File | What it does | Output (in `results/`) |
|---|---|---|---|
| 0 | `settings.jl` | Every number used: dose, sampling times, error, CVs, study sizes, seeds | — |
| 0 | `model_perarm.jl`, `model_joint.jl` | The two PopPK models (equations in the header) | — |
| 0 | `be_functions.jl` | Fast paired BE test for the simulated trials (steps 5–6), interpolation between formulations | — |
| 1 | `step1_build_dataset.jl` | Samples the PK-Sim profiles at 13 times, adds 15% error | `step1_dataset.csv` |
| 2 | `step2_fit_poppk.jl` | Fits both models (FOCE) | `<model>/step2_parameters.csv`, `_standard_errors.csv`, `_ipred.csv` |
| 3 | `step3_endpoints.jl` | Cmax and AUC 0–24 h per subject from the PopPK predictions (Pumas `NCA`); checks them against PK-Sim | `<model>/step3_poppk_endpoints.csv`, `step3_check_vs_pbpk.csv` |
| 4 | `step4_gmr.jl` | Test/reference GMR with 90% CI per formulation, PBPK vs PopPK (Pumas `pumas_be`) | `<model>/step4_gmr.csv`, `fig_step4_gmr.png` |
| 5 | `step5_power_safe_space.jl` | Power vs T50 and the safe space (power ≥ 80%) | `<model>/step5_power.csv`, `step5_safe_space.csv`, `fig_step5_power.png` |
| 6 | `step6_type1_error.jl` | BE boundary and the type I error there | `<model>/step6_boundaries.csv`, `step6_type1_error.csv`, `step6_tost_control.csv`, figures |
| 6b | `step6b_check_vs_pumas_be.jl` | Cross-check: our fast BE test vs `pumas_be` on the same simulated trials, at every boundary | `step6b_check_vs_pumas_be.csv` |
| 7 | `step7_ivivc.jl` | IVIVC check: PopPK release time vs the PK-Sim T50 (branch `ivivc-qc`) | `step7_ivivc.csv`, `step7_ivivc_fit.csv`, `fig_step7_ivivc.png` |
| 8 | `step8_heldout_test.jl` | Held-out test: train on T50 1, 3, 7 h, predict 2, 3.5, 4, 5, 6 h (uses `model_joint_frel.jl`) | `heldout/step8_*.csv`, `heldout/fig_step8_gmr.png` |
| 9 | `step9_one_test_formulation.jl` | Fit to the reference + one test formulation (5 h, or 2 h as a check), predict the other 6; ka fixed at 0.3, 1, 3; saturable and gastric-emptying variants tried (uses `model_window.jl`, `model_window_sat.jl`, `model_window_ge.jl`) | `one_test/step9_*.csv`, `one_test/fig_step9_gmr.png` |

`<model>` is `perarm` or `joint`. Inputs are in `data/`:
- `pksim_profiles/`: PK-Sim plasma profiles, 0–24 h, every 3 min.
- `pbpk_endpoints.csv`: PK-Sim Cmax (µmol/L) and AUC 0–24 h (µmol·min/L) per subject.

## The two PopPK models

Both models use a Weibull release feeding a one-compartment model, with variability on release, CL
and Vc and a proportional error.

- **perarm** is fitted separately to each formulation (8 fits of 24 subjects).
- **joint** is fitted once to all 8 formulations. CL, Vc and the shape are shared, and the release
  time scales with T50 as `Alpha = tvalpha·(T50/3)^slope`.

## How to run

You need Julia 1.11 with Pumas. Step 2 uses `Pumas`, step 3 `NCA`, steps 4 and 6b `Bioequivalence`
(all part of the Pumas distribution). From this folder:

```bash
export JULIA="julia --project=."     # or your Pumas command, e.g. with a Pumas sysimage
./run_all.sh                         # all six steps, ~10 min; logs in results/logs/
STEPS="step4_gmr" ./run_all.sh       # one step only
```

## QC check: rerun vs the saved results

On 19 Nov 2026 all six steps were rerun and compared with the saved 30 Oct results
(`learning/metformin_poppk_*_export_edited/poppk_pumas/output_fixed/`). Every value is identical:

| Compared | Per-arm | Joint |
|---|---|---|
| Dataset (step 1) | byte-identical | byte-identical |
| Population estimates (step 2) | identical | identical |
| Cmax and AUC per subject (step 3) | identical | identical |
| GMR, 90% CI, pass/fail (step 4) | identical | identical |
| Safe-space edges (step 5) | identical | identical |
| Boundaries, type I error, TOST control (step 6) | identical | identical |

### Standard Pumas functions (19 Nov 2026)

The endpoints and the BE test now use the standard Pumas packages wherever one exists:

- **Step 3: `NCA.cmax` and `NCA.auc`** (linear-up/log-down) replace our own trapezoid code. Cmax is
  identical and AUC differs by at most 3×10⁻¹² (relative), so steps 5–6 did not need rerunning.
- **Step 4: `pumas_be`** (2x2 crossover ANOVA) replaces our paired t-test. GMRs are identical
  (≤ 10⁻¹³). The 90% CI limits move by at most 0.24%, because the ANOVA also fits a period and a
  sequence effect, and no pass/fail decision changes. The reference row (reference vs itself) is
  dropped, since it has no CI.
- **Steps 5–6 keep the fast paired test (`be_test`)** for their hundreds of thousands of simulated
  trials. **Step 6b** checks this test against `pumas_be`. At every BE boundary, where the true
  GMR is exactly 0.80 or 1.25 and decisions are closest, 2,000 trials ran through both tests:

| | Result |
|---|---|
| Same pass/fail decision | 99.75–99.95% of trials (18 boundary × study-size cases) |
| Largest difference in pass rate | 0.15 percentage points |

So the type I error and power from the fast test match what `pumas_be` gives.

## Key results (rich sampling)

| | PBPK | PopPK per-arm | PopPK joint |
|---|---|---|---|
| Safe space, joint Cmax + AUC, n = 48 (T50, h) | 2.05–4.00 | 2.10–3.90 | 2.10–3.95 |
| Safe space, n = 86 (T50, h) | 1.90–4.10 | 1.95–4.05 | 2.05–4.10 |
| Upper BE boundary, joint (T50, h) | 4.459 | 4.377 | 4.484 |
| PBPK (true) pass rate at that boundary, n = 86 | 4.5% | 14.4% | 2.9% |

The last row is how often a real BE study passes at the boundary each method would declare
(PBPK is the truth). A well-placed boundary gives 5%.
- **Above 5%:** the method's boundary sits inside the true one, which is too strict.
- **Below 5%:** it sits outside the true one, which is too lenient. Formulations between the two
  boundaries would be called BE by the model but are not.

## What to check

1. `step1`: dose 390 mg base = 3019.5 µmol; error applied as `true × (1 + 0.15·z)`.
2. `model_*.jl`: Weibull hazard and units (h, µmol, L, µmol/L).
3. `step2`: lag fixed; FOCE; all fits converge (`results/logs/step2_fit_poppk.log`).
4. `step3`: PopPK/PBPK geometric-mean ratios close to 1 (`step3_check_vs_pbpk.csv`).
5. `be_functions.jl`: the trial simulation and the log-log interpolation; `step6b` shows its
   paired TOST agrees with `pumas_be`.
6. `step6`: the TOST control passes about 5% (reference × 0.80).

## Step 7: how PopPK relates to PBPK (IVIVC check)

PopPK is linked to PBPK in two ways:
- **Through the data.** The PopPK fits use concentrations sampled from the PK-Sim profiles.
- **Through the decisions.** Steps 4–6 run the same BE code on both sets of endpoints.

Step 7 checks a third possible link: does the PopPK release follow the release set in PK-Sim (an
in vivo IVIVC)? It compares the PopPK time to 50% released, `Alpha·ln(2)^(1/Gamma)`, with the
PK-Sim T50.

| | Estimate | 95% CI | An IVIVC would give |
|---|---|---|---|
| Per-arm fits: slope of log(PopPK T50) on log(PK-Sim T50) | −0.012 | −0.045 to 0.021 | 1 |
| Joint model: estimated `slope` | 0.029 | −0.029 to 0.087 | 1 |

**Result: no IVIVC.** The PopPK release time stays at 2.4–2.6 h for every formulation, while the
PK-Sim T50 runs from 1 to 7 h. The formulation effect is absorbed into the apparent clearance
instead (per-arm CL rises from 35 to 59 L/h, `step7_ivivc.csv`). The PopPK model therefore
reproduces each formulation's Cmax and AUC well enough for BE (steps 4–6), but its release
parameters are not the PK-Sim release. They should not be read as an in vitro or in vivo
release-rate estimate.

## Step 8: can PopPK predict a formulation it was not fitted to?

Train on T50 = 1, 3 and 7 h only, then predict the 5 held-out formulations. Each subject keeps the
random effects from its reference (3 h) arm, and only T50 changes. Three methods are compared:
- **joint**: `model_joint.jl`, with T50 acting on the release time only.
- **joint_frel**: `model_joint_frel.jl`, which also lets bioavailability depend on T50:
  F = (T50/3)^θ. The estimate is θ = −0.37.
- **interp3**: no model; per-arm PopPK endpoints at 1, 3 and 7 h, interpolated between them.

Cmax GMR at the held-out formulations:

| T50 (h) | PBPK (truth) | joint | joint_frel | interp3 |
|---|---|---|---|---|
| 2.0 | 1.178 | 1.004 | 1.165 | 1.126 |
| 3.5 | 0.924 | 0.998 | 0.944 | 0.899 |
| 4.0 | 0.857 | 0.997 | 0.898 | 0.820 |
| 5.0 | 0.744 | 0.995 | 0.825 | 0.703 |
| 6.0 | 0.658 | 0.993 | 0.771 | 0.620 |

Safe space, joint Cmax + AUC, with power ≥ 80% (T50, h):

| n | PBPK | joint | joint_frel | interp3 |
|---|---|---|---|---|
| 48 | 2.05–4.00 | 1.0–7.0 (everything) | 2.00–4.45 | 1.75–3.75 |
| 86 | 1.90–4.10 | 1.0–7.0 (everything) | 1.90–4.70 | 1.70–3.85 |

**Result: no method predicts the held-out formulations well enough for BE.**
- **joint** predicts no formulation effect at all (GMR ≈ 1 everywhere), so every T50 looks
  bioequivalent.
- **joint_frel** gets the direction right but underpredicts the effect away from the reference.
  Its safe space runs to 4.45–4.70 h instead of 4.0–4.1 h, so it would call slower
  formulations BE that are not. A single power law on F cannot fit the data: PK-Sim's Cmax falls
  faster than its AUC, but a change in F moves both equally. The fit is off even at the trained
  1 h and 7 h formulations.
- **interp3** shifts the safe space towards faster release (1.75–3.75 h).

So PopPK, even with T50 as a covariate on release and on bioavailability, cannot stand in for PBPK
when predicting untested release rates for metformin. The per-arm results (steps 4–6) hold only for
formulations that were actually simulated. Frel was also tested on all 8 formulations (31 Oct–2 Nov
2026, not part of this branch): it improved the fit (ΔAIC −95) but worsened the BE decisions, and
the base joint model was kept.

## Step 9: one tested formulation → BE for all the others (19 Nov 2026)

Steps 7–8 show why the PopPK models cannot extrapolate: they estimate the release from plasma
data, where it cannot be told apart from absorption. Step 9 changes what is estimated:

- **Release is given, not estimated.** Each formulation's own dissolution curve, the PK-Sim Weibull
  with its T50, shape 0.92 and no lag, is the model input (`model_window.jl`).
- **Absorption window is estimated.** Absorption rate = `ka·exp(−t/τ)`: drug released late is
  partly never absorbed. **Disposition (CL, Vc) is shared** by all formulations.
- **Test:** fit to the reference (3 h) plus ONE test formulation, then predict the other 6 from
  their dissolution curves. The tested formulation is 5 h (main case) or 2 h (check).
- **ka is fixed.** Fitted freely, ka ran to its 0.01/h bound and CL (1.1 L/h) and Vc (1.9 L) lost
  their meaning: from oral data, the absorption scale cannot be separated from CL and Vc. So ka
  is fixed at 0.3, 1 and 3 /h, to show how much the predictions depend on that choice.

Fits (`results/one_test/step9_parameters.csv`):

| Tested | ka (1/h) | −2LL | τ (h) | CL (L/h) | Vc (L) |
|---|---|---|---|---|---|
| 5 h | free → 0.01 (bound) | 1833 | 4.9 | 1.1 | 1.9 |
| 5 h | **0.3 (fixed)** | **1863** | **5.3** | **23** | **50** |
| 5 h | 1 (fixed) | 1989 | 5.5 | 39 | 121 |
| 5 h | 3 (fixed) | 2164 | 2.0 | 27 | 162 |
| 2 h | free → 0.01 (bound) | 2051 | 4.5 | 1.0 | 2.1 |
| 2 h | **0.3 (fixed)** | **2092** | **5.0** | **23** | **55** |
| 2 h | 1 (fixed) | 2171 | → ∞ (no window) | 43 | 159 |
| 2 h | 3 (fixed) | 2435 | 2.2 | 31 | 168 |

**ka = 0.3 /h is the chosen model.** It is the best-fitting fixed value, its CL and Vc are
physiological, and it gives nearly the same answer whichever formulation was tested (5 h or 2 h).
With ka = 1 and 2 h tested, the window collapsed (τ → ∞), so that fit is degenerate.

Cmax and AUC GMR, ka = 0.3 (truth = PBPK; "all 8" = per-arm PopPK fitted to every formulation):

| T50 (h) | PBPK Cmax | tested 5 h | tested 2 h | all 8 | PBPK AUC | tested 5 h | tested 2 h | all 8 |
|---|---|---|---|---|---|---|---|---|
| 1.0 | 1.363 | 1.636 | 1.642 | 1.379 | 1.195 | 1.328 | 1.346 | 1.188 |
| 2.0 | 1.178 | 1.236 | 1.238 | 1.176 | 1.120 | 1.145 | 1.151 | 1.109 |
| 3.5 | 0.924 | 0.914 | 0.914 | 0.924 | 0.941 | 0.940 | 0.938 | 0.939 |
| 4.0 | 0.857 | 0.843 | 0.842 | 0.850 | 0.887 | 0.886 | 0.883 | 0.886 |
| 5.0 | 0.744 | 0.731 | 0.730 | 0.731 | 0.792 | 0.796 | 0.791 | 0.795 |
| 6.0 | 0.658 | 0.646 | 0.645 | 0.625 | 0.714 | 0.722 | 0.716 | 0.718 |
| 7.0 | 0.589 | 0.580 | 0.578 | 0.558 | 0.651 | 0.661 | 0.654 | 0.648 |

Safe space, joint Cmax + AUC, power ≥ 80% (T50, h), and the BE boundaries:

| | PBPK | ka 0.3, tested 5 h | ka 0.3, tested 2 h | all 8 |
|---|---|---|---|---|
| Safe space, n = 48 | 2.05–4.00 | 2.25–3.90 | 2.25–3.90 | 2.10–3.90 |
| Safe space, n = 86 | 1.90–4.10 | 2.15–4.00 | 2.15–4.00 | 1.95–4.05 |
| Lower boundary (true pass rate there, n = 86) | 1.51 (5%) | 1.95 (87%) | 1.95 (88%) | 1.54 (8%) |
| Upper boundary (true pass rate there, n = 86) | 4.46 (5%) | 4.34 (22%) | 4.33 (24%) | 4.38 (14%) |

**What this shows:**
- **Slower release (3.5–7 h) is predicted well** from two arms: Cmax and AUC GMRs within about 2%
  of PBPK, as good as fitting all 8 formulations. This holds even when the tested formulation
  was the faster 2 h one.
- **Faster release is overpredicted**: Cmax GMR 1.64 vs 1.36 at 1 h. This is the case **even at
  2 h when 2 h was the tested formulation** (1.24 vs 1.18). So the bias is not about
  extrapolation or about ka. The absorption-window structure itself makes Cmax too sensitive to
  fast release.
- **The safe space is therefore too strict, never too lenient.** The lower edge is 2.15–2.25 h
  against 1.90–2.05 h. The true pass rate at both predicted boundaries is above 5% (22–88%), so
  no formulation that is not BE gets called BE.
- The fixed ka value matters for AUC at the extremes (ka 1 or 3 move the 7 h AUC GMR by ±10%), so
  ka = 0.3 should be read as the best of the tested values, not as a measured absorption rate.

### Saturable uptake: tried, did not help

Metformin uptake is transporter-mediated, so a saturable version was fitted
(`model_window_sat.jl`): uptake = `ka·Km·exp(−t/τ)·Gut/(Km + Gut)`, with ka fixed at 0.3 and 1. It
contains the linear model (Km → ∞), so −2LL compares directly.

| Tested | ka (1/h) | Km (µmol) | −2LL saturable | −2LL linear | Result |
|---|---|---|---|---|---|
| 5 h | 0.3 | 4.5×10⁹ (→ ∞) | 1863 | 1863 | identical to linear |
| 5 h | 1 | 6.5×10⁶ (→ ∞) | 1989 | 1989 | identical to linear |
| 2 h | 0.3 | 2.0×10¹⁰ (→ ∞) | 2092 | 2092 | identical to linear |
| 2 h | 1 | 1,405 | 2127 | 2171 | window lost (τ → ∞) |

- With the chosen ka = 0.3, **the data push Km to infinity**: no saturation is supported, and the
  predictions are the linear model's.
- The one fit that did saturate lost its absorption window instead. Fast release improved
  (Cmax GMR at 1 h: 1.45 vs PBPK 1.36), but slow release broke (AUC GMR at 7 h: 0.85 vs 0.65), and
  its upper boundary is too lenient (true pass rate there 0.1%, against 5%).
- **Conclusion:** saturable uptake does not explain the fast-release bias. The chosen step 9 model
  stays the linear window model with ka = 0.3.

### Gastric emptying: much better fit, same BE predictions

`model_window_ge.jl` adds a stomach step: released drug must empty into the small intestine (first
order, kge) before the absorption window. kge → ∞ gives the linear model back.

| Tested | ka (1/h) | kge (1/h) | −2LL with stomach | −2LL linear | CL (L/h) | Vc (L) | Cmax GMR at 1 h (PBPK 1.363) |
|---|---|---|---|---|---|---|---|
| 5 h | 0.3 | 1.09 | 1475 | 1863 | 19 | 18 | 1.684 |
| 5 h | 1 | 0.81 | 1449 | 1989 | 34 | 43 | 1.652 |
| 2 h | 0.3 | 0.99 | 1616 | 2092 | 18 | 19 | 1.686 |
| 2 h | 1 | 0.60 | 1631 | 2171 | 31 | 39 | 1.647 |

- **The profiles are described much better:** −2LL falls by 390–540 points for one extra parameter.
  kge settles at 0.6–1.1 /h, a gastric half-time of 0.6–1.2 h. ka is no longer pushed to its lower
  bound: with 5 h tested, ka = 1 now fits better than 0.3.
- **The BE predictions do not change:** Cmax GMR at 1 h is still 1.65–1.69 against 1.36, and the
  safe space is the same or slightly stricter (n = 86: 2.15–2.20 to 3.95–4.00 h, against PBPK
  1.90–4.10 h). The true pass rate at every predicted boundary stays above 5% (14–92%), so it is
  still never too lenient.
- With ka = 1 the parameters are physiological (CL 31–34 L/h, Vc 39–43 L) but slow-release AUC is
  overpredicted by up to 9% (7 h: 0.71 vs 0.65); with ka = 0.3 AUC is right but Vc is low (18 L).

**Conclusion from steps 9a–c.** Three absorption mechanisms were tried (fixed ka, saturable uptake,
gastric emptying). None removes the fast-release Cmax bias, although gastric emptying fits the
profiles far better. A compartmental PopPK model fitted to two arms predicts slower release well,
but extrapolating to faster release seems to need what PK-Sim has and these models lack, such as
absorption that varies by intestinal segment. For faster release, PBPK stays the tool.

## Next steps

1. **Repeat step 9 with sparse sampling.**
2. **Felodipine:** the same idea, with gut-wall CYP3A4 extraction depending on where the drug is
   released in place of the absorption window.

Not in this branch: sparse sampling, AUC to infinity, and the PK-Sim model build.
