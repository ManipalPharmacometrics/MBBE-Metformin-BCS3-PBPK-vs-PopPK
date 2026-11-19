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
| True type I error at that boundary, n = 86 | 4.5% | 14.4% | 2.9% |

The last row is the PBPK pass rate at the boundary each method would declare. PBPK is the truth,
so it is the real type I error of a safe space set by that method; 5% is correct.

## What to check

1. `step1`: dose 390 mg base = 3019.5 µmol; error applied as `true × (1 + 0.15·z)`.
2. `model_*.jl`: Weibull hazard and units (h, µmol, L, µmol/L).
3. `step2`: lag fixed; FOCE; all fits converge (`results/logs/step2_fit_poppk.log`).
4. `step3`: PopPK/PBPK geometric-mean ratios close to 1 (`step3_check_vs_pbpk.csv`).
5. `be_functions.jl`: the 90% CI TOST, the trial simulation and the log-log interpolation.
6. `step6`: the TOST control passes about 5% (reference × 0.80).

Not in this branch: sparse sampling, AUC to infinity, and the PK-Sim model build.
