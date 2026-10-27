# START HERE — metformin XR Weibull PopPK ⇄ PBPK safe space

This folder is self-contained. Copy the whole `poppk-per-arm/` folder to
your **licensed Pumas / JuliaHub** environment and run one command.

## ▶ The file to run

```bash
cd poppk_pumas
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # first time only
julia --project=. run_all.jl                          # <-- RUN THIS
```

`run_all.jl` runs the four steps in order, each in a fresh Julia process:

1. `scripts/01_build_dataset.jl` — sparse-sample the 8 PBPK profiles + add 15% assay error
2. `scripts/02_fit_poppk.jl` — **the Pumas fit** (Weibull-release model, per formulation) ← needs Pumas
3. `scripts/03_poppk_metrics.jl` — NCA on the fitted profiles
4. `scripts/04_be_safespace_overlay.jl` — PBPK vs PopPK GMR + safe-space overlay figures

Only **step 02** needs a Pumas license; the rest use open packages.

## ▶ What to send back

After it runs, everything lands in `poppk_pumas/output/`. Send back that folder — or at
minimum these two files and I'll rebuild the rest:

- `poppk_pumas/output/poppk_ipred_profiles.csv`  (most important — profiles)
- `poppk_pumas/output/poppk_param_estimates.csv` (population parameters per arm)

Full list and details are in `poppk_pumas/RUN_ON_PUMAS.md`.

## ▶ If Pumas errors

See the "If something breaks" section of `poppk_pumas/README.md` — the only
version-sensitive spots are in `scripts/02_fit_poppk.jl` (`icoef`/`predict` field names,
and the time-varying `@dynamics`; set `tvlag init=0.0` if the solver struggles at the lag).

## ▶ Files in this folder

| Path | What it is |
|---|---|
| `poppk_pumas/run_all.jl` | **the entry point — run this** |
| `poppk_pumas/src/poppk_model.jl` | the Weibull-release PopPK model (baseline) |
| `poppk_pumas/scripts/02_fit_poppk.jl` | the Pumas fit |
| `poppk_pumas/data/` | the 8 PBPK population profiles + PBPK endpoints (inputs) |
| `poppk_pumas/RUN_ON_PUMAS.md` | detailed run + share-back instructions |
| `poppk_pumas/README.md` | pipeline reference |
| `METHODOLOGY_metformin_MBBE.md` | plain-language scientific methodology (PBPK + PopPK) |

## ▶ Note on the optional no-Pumas script

`scripts/02b_fit_weibull_nopumas.jl` (+ `src/weibull_absorption.jl`) is an **optional**
license-free cross-check that fits the same Weibull model with open packages. It is **not**
part of `run_all.jl` and is **not** needed on the Pumas box — ignore it there. (To use it
you would add `Optim` to the environment.)
