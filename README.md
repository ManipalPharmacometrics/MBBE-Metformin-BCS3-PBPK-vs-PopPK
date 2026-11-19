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
| 0 | `be_functions.jl` | BE test, simulated BE trials, interpolation between formulations | — |
| 1 | `step1_build_dataset.jl` | Samples the PK-Sim profiles at 13 times, adds 15% error | `step1_dataset.csv` |
| 2 | `step2_fit_poppk.jl` | Fits both models (FOCE) | `<model>/step2_parameters.csv`, `_standard_errors.csv`, `_ipred.csv` |
| 3 | `step3_endpoints.jl` | Cmax and AUC 0–24 h per subject from the PopPK predictions; checks them against PK-Sim | `<model>/step3_poppk_endpoints.csv`, `step3_check_vs_pbpk.csv` |
| 4 | `step4_gmr.jl` | Test/reference GMR with 90% CI per formulation, PBPK vs PopPK | `<model>/step4_gmr.csv`, `fig_step4_gmr.png` |
| 5 | `step5_power_safe_space.jl` | Power vs T50 and the safe space (power ≥ 80%) | `<model>/step5_power.csv`, `step5_safe_space.csv`, `fig_step5_power.png` |
| 6 | `step6_type1_error.jl` | BE boundary and the type I error there | `<model>/step6_boundaries.csv`, `step6_type1_error.csv`, `step6_tost_control.csv`, figures |
| 7 | `step7_ivivc.jl` | IVIVC check: PopPK release time vs the PK-Sim T50 (branch `ivivc-qc`) | `step7_ivivc.csv`, `step7_ivivc_fit.csv`, `fig_step7_ivivc.png` |
| 8 | `step8_heldout_test.jl` | Held-out test: train on T50 1, 3, 7 h, predict 2, 3.5, 4, 5, 6 h (uses `model_joint_frel.jl`) | `heldout/step8_*.csv`, `heldout/fig_step8_gmr.png` |
| 9 | `step9_one_test_formulation.jl` | Fit to the reference + one test formulation (5 h), predict the other 6 (uses `model_window.jl`) | `one_test/step9_*.csv`, `one_test/fig_step9_gmr.png` |

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

You need Julia 1.11 with Pumas; only step 2 uses Pumas. From this folder:

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
5. `be_functions.jl`: the 90% CI TOST, the trial simulation and the log-log interpolation.
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

## Step 9: one tested formulation → BE for all the others (first run, 19 Nov 2026)

Steps 7–8 show why the PopPK models cannot extrapolate: they estimate the release from plasma
data, where it cannot be told apart from absorption. Step 9 changes what is estimated:

- **Release is given, not estimated.** Each formulation's own dissolution curve, the PK-Sim Weibull
  with its T50, shape 0.92 and no lag, is the model input (`model_window.jl`).
- **Absorption window is estimated.** Absorption rate = `ka·exp(−t/τ)`: drug released late is
  partly never absorbed. **Disposition (CL, Vc) is shared** by all formulations.
- **Test:** fit to the reference (3 h) and one test formulation (5 h) only, then predict the other
  6 from their dissolution curves (`step9_one_test_formulation.jl`, `results/one_test/`).

Cmax GMR (truth = PBPK; "all 8" = per-arm PopPK fitted to every formulation, for comparison):

| T50 (h) | | PBPK | Step 9 (fitted to 3 h + 5 h) | Per-arm, all 8 fitted |
|---|---|---|---|---|
| 1.0 | predicted | 1.363 | 1.662 | 1.379 |
| 2.0 | predicted | 1.178 | 1.245 | 1.176 |
| 3.5 | predicted | 0.924 | 0.912 | 0.924 |
| 4.0 | predicted | 0.857 | 0.839 | 0.850 |
| 5.0 | fitted | 0.744 | 0.724 | 0.731 |
| 6.0 | predicted | 0.658 | 0.639 | 0.625 |
| 7.0 | predicted | 0.589 | 0.572 | 0.558 |

Safe space, joint Cmax + AUC, power ≥ 80% (T50, h):

| n | PBPK | Step 9 | Per-arm, all 8 fitted |
|---|---|---|---|
| 48 | 2.05–4.00 | 2.25–3.85 | 2.10–3.90 |
| 86 | 1.90–4.10 | 2.20–4.00 | 1.95–4.05 |

**What this shows:**
- **On the side that was tested (slower release, 3.5–7 h)** the model predicts the PBPK GMRs within
  about 3%. That is as good as fitting every formulation, from two arms only.
- **On the untested side (faster release, 1–2 h)** it overpredicts the effect: Cmax GMR 1.66 vs
  1.36 at 1 h.
- **The safe space comes out too strict on both sides,** never too lenient: 2.25–3.85 h against
  2.05–4.00 h. The PBPK pass rate at the step 9 boundaries is 20–91%, all above 5%. So no formulation
  that is truly not BE gets called BE, but some that are BE get excluded.

**Caveat (must fix before relying on it):**
- `ka` is on its lower bound (0.01/h), and CL = 1.06 L/h and Vc = 1.9 L are not physiological. With
  oral data only, F, CL and Vc cannot be separated. The fit drifted to a corner where very little is
  absorbed and CL and Vc are scaled down to match.
- The predicted concentrations are unaffected: CL/F ≈ 40 L/h, as in the per-arm fits. The
  individual parameter values are meaningless, though, and the release-time effect may depend on
  that corner.

## Next steps

1. **Fix the scale problem in `model_window.jl`.** Either fix ka at a physiological value and
   estimate τ, CL/F and Vc/F, or write the model in CL/F and Vc/F form. Then rerun step 9.
2. **Repeat with 2 h as the tested formulation**, which tests extrapolation in the other direction.
3. **Repeat with sparse sampling.**
4. **Felodipine:** the same idea, with gut-wall CYP3A4 extraction depending on where the drug is
   released in place of the absorption window.

Not in this branch: sparse sampling, AUC to infinity, and the PK-Sim model build.
