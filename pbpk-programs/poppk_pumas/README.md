# PopPK ⇄ PBPK decision-equivalence (Pumas)

Self-contained Julia/Pumas pipeline for the metformin XR MBBE project. It takes
the dense PK-Sim **PBPK** population profiles, builds a simulated BE study from
them (sparse sampling + residual error), fits a **PopPK** model in Pumas,
re-derives the BE endpoints, and overlays the **PopPK vs PBPK** virtual-BE
safe-space boundaries — the core "do the two methods reach the same decision?"
comparison.

> **Run this on a licensed Pumas / JuliaHub environment.** It was authored on a
> host where Pumas can't run (the open-registry `Pumas` package there is an
> empty stub — the real proprietary source loads only in a license-authenticated
> JuliaHub session), so steps 02 (the Pumas fit) have **not been executed/tested
> here**. Steps 01, 03, 04 use only open packages. See "If something breaks".

## Layout

```
poppk_pumas/
  Project.toml                  deps (Pumas + open pkgs)
  run_all.jl                    runs scripts 01-04 in order (fresh process each)
  src/
    poppk_model.jl              the @model (1-cmt, WEIBULL-release absorption, prop. error)
    weibull_absorption.jl       open-Julia Weibull simulator + per-subject NLS fit (no Pumas)
    nca_lite.jl                 minimal NCA (Cmax, AUC_tEnd, AUC_inf)
    be_utils.jl                 TOST / VBE / safe-space — SAME code as the PBPK side
  scripts/
    01_build_dataset.jl         sparse-sample PBPK profiles + 15% prop. residual error
    02_fit_poppk.jl             Pumas fit per formulation -> params, icoef, dense IPRED
    02b_fit_weibull_nopumas.jl  Pumas-FREE substitute for 02 (same outputs, open pkgs only)
    03_poppk_metrics.jl         NCA on IPRED -> per-subject Cmax/AUC (PBPK schema)
    04_be_safespace_overlay.jl  PBPK vs PopPK GMR + safe-space overlay figures
  data/
    popsims/*-Results.csv       8 dense PBPK population profiles (24 subj, 0-24 h)
    pbpk_pk_parameters_individual.csv   PBPK per-subject BE endpoints (for the overlay)
  output/                       all generated CSVs + figures land here
```

## How to run

```bash
# from this folder, on the licensed environment
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # first time only
julia --project=. run_all.jl
```

Or step by step:

```bash
julia --project=. scripts/01_build_dataset.jl
julia --project=. scripts/02_fit_poppk.jl
julia --project=. scripts/03_poppk_metrics.jl
julia --project=. scripts/04_be_safespace_overlay.jl
```

## Design choices (baseline, locked 2026-06-10)

| Choice | Value |
|---|---|
| Sampling schedule | rich BE: 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 16, 24 h (predose = dose only) |
| Residual error | proportional, 15% CV |
| Structural model | 1-compartment, **Weibull-release absorption** (α = in-vivo release scale ~ T50, γ = shape ~ β) |
| IIV | diagonal Ω on α, CL, Vc |
| Reference formulation | T50 = 3.0 h |
| Dose | 390 mg free base = 3019.5 µmol (500 mg HCl tablet equivalent) |
| Fit | FOCE, each formulation arm fit independently |
| VBE | 500 trials, n = 48 and 86, endpoint-specific within-subject CV (Cmax 0.122, AUC 0.112) |

The BE/VBE/safe-space statistics (`src/be_utils.jl`) are byte-for-byte the same
helpers used to produce the PBPK-side results, so the overlay is a true
apples-to-apples comparison.

## Outputs (in `output/`)

- `sparse_be_dataset.csv` — the simulated BE study (Pumas NONMEM format)
- `poppk_param_estimates.csv` — population parameters + OFV per formulation
- `poppk_icoef.csv` — per-subject Ka/CL/Vc (EBEs)
- `poppk_ipred_profiles.csv` — dense individual predictions
- `poppk_pk_parameters_individual.csv` — PopPK per-subject Cmax/AUC (PBPK schema)
- `poppk_nca_diagnostics.csv` — + Tmax, %extrapolated, λz R²
- `be_overlay_gmr.csv`, `fig_overlay_gmr.png` — GMR(T50) ± 90% CI, both methods
- `safe_space_overlay_boundaries.csv`, `fig_overlay_safespace_*.png` — safe-space overlay

## Running without a Pumas license (cross-check)

Step 02 is the only step that needs Pumas. `scripts/02b_fit_weibull_nopumas.jl`
reproduces it with open packages only — it fits the **same** Weibull-release model
per subject by nonlinear least squares (Optim.jl) and writes the **same** output
files, so 03/04 run identically. Use it to run/verify the full pipeline anywhere:

```bash
julia --project=nopumas_env scripts/01_build_dataset.jl
julia --project=nopumas_env scripts/02b_fit_weibull_nopumas.jl
julia --project=nopumas_env scripts/03_poppk_metrics.jl
julia --project=nopumas_env scripts/04_be_safespace_overlay.jl
```

The Pumas NLME fit (`02`) is authoritative; the open-Julia two-stage fit (`02b`) is a
slightly conservative directional cross-check. To run the real fit and share results,
see `RUN_ON_PUMAS.md`.

## Alternative / sensitivity models

1. **Truly-sparse schedule (~6 pts)** — in `01_build_dataset.jl` set
   `const SCHEDULE = SPARSE_SCHEDULE`. Stresses the fit; shows robustness under
   limited sampling.
2. **First-order absorption + lag** — the earlier structural form (T50 absorbed
   implicitly into an apparent `ka`). Kept only as a sensitivity comparison against
   the Weibull baseline.
3. **Joint T50-covariate model** — fit all arms together with T50 as a covariate on
   `α` to estimate an explicit T50→release relationship.

## If something breaks

The only Pumas-version-sensitive code is in `02_fit_poppk.jl`:

- **`icoef(fpm)` field names** — the per-subject coefficients are read as
  `ico.Alpha`, `ico.CL`, `ico.Vc`. If your Pumas version names them differently,
  adjust there.
- **`predict(fpm; obstimes=...)` access** — dense predictions are read as
  `pr.ipred.dv` and `pr.time`. If the prediction object differs, adjust the loop.
- **Time-varying `@dynamics`** — the Weibull release hazard references `t`. If the
  solver struggles at the lag discontinuity, set `tvlag` `init = 0.0` in
  `src/poppk_model.jl` (the BE schedule barely identifies a lag) — the Weibull shape
  absorbs the onset.
- **Pumas-free fallback** — if step 02 can't run at all, `02b_fit_weibull_nopumas.jl`
  regenerates `poppk_ipred_profiles.csv` from the same model without Pumas.

`deviance(fpm)` is wrapped in `try`/`catch` (named `OFV`); if it errors on your
version it just records `missing`.
