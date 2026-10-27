# Run the Weibull PopPK fit on a licensed Pumas environment

You run steps **01 → 02 → 03 → 04** on a Pumas/JuliaHub session, then paste the output
CSVs back here. Only **step 02** actually needs Pumas; 01/03/04 use open packages, but
running the whole chain in one place is simplest.

## 1. Copy this folder to the Pumas machine
Copy the entire `programs/poppk_pumas/` directory. The pieces the fit needs:

```
poppk_pumas/
  Project.toml                     # Pumas + open deps
  run_all.jl                       # runs 01→04 in order
  src/poppk_model.jl               # the Weibull-release @model  (BASELINE)
  scripts/01_build_dataset.jl      # sparse-sample + 15% residual error
  scripts/02_fit_poppk.jl          # THE PUMAS NLME FIT (Weibull)
  scripts/03_poppk_metrics.jl      # NCA on predictions
  scripts/04_be_safespace_overlay.jl
  data/popsims/*-Results.csv       # 8 PBPK profiles (already here)
  data/pbpk_pk_parameters_individual.csv
```

## 2. Instantiate once, then run
```bash
cd poppk_pumas
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # first time only
julia --project=. run_all.jl                           # 01 → 02 → 03 → 04
```

If you prefer step-by-step:
```bash
julia --project=. scripts/01_build_dataset.jl
julia --project=. scripts/02_fit_poppk.jl        # <- the fit; needs Pumas
julia --project=. scripts/03_poppk_metrics.jl
julia --project=. scripts/04_be_safespace_overlay.jl
```

## 3. What the fit does (so you know what to expect)
`02_fit_poppk.jl` fits the Weibull-release 1-cmt model to **each of the 8 T50 arms
independently**. Watch the printed population table: `tvalpha` (in-vivo release scale)
should increase with T50 and roughly track it; `tvgamma` ≈ 1; `tvcl`/`tvvc` stable.

## 4. Share these files back
Paste (or upload) the contents of `output/` — at minimum these CSVs:

| File | Why I need it |
|---|---|
| `output/poppk_param_estimates.csv` | population params per arm (α, γ, CL, Vc, Ω, σ) |
| `output/poppk_icoef.csv` | per-subject α/CL/Vc |
| `output/poppk_ipred_profiles.csv` | dense predicted profiles — **most important**; 03/04 rebuild everything from this |
| `output/poppk_pk_parameters_individual.csv` | per-subject Cmax/AUC (from 03) |
| `output/be_overlay_gmr.csv` | GMR(T50) both methods (from 04) |
| `output/safe_space_overlay_boundaries.csv` | the safe-space boundary table (from 04) |

If you'd rather run only the fit on Pumas, share just `poppk_ipred_profiles.csv` and
`poppk_param_estimates.csv` — I can run 03 and 04 here from those.

## 5. Known version-sensitive spots (only in step 02)
If Pumas errors, it's almost certainly one of these — all isolated to `02_fit_poppk.jl`:

- **`icoef(fpm)` field names.** Read as `ico.Alpha`, `ico.CL`, `ico.Vc`. If your Pumas
  version names the `@pre` outputs differently, adjust there.
- **`predict(fpm; obstimes=…)` access.** Dense predictions read as `pr.ipred.dv` and
  `pr.time`. Adjust the loop if the prediction object differs.
- **Time-varying `@dynamics`.** The model references `t` in the release hazard. If the
  solver struggles at the lag discontinuity, set `tvlag` `init = 0.0` in
  `src/poppk_model.jl` (the BE schedule barely identifies a lag anyway) — the Weibull
  shape absorbs the onset.
- **`deviance(fpm)`** is wrapped in try/catch (records `missing` if unavailable).

Then I'll overlay your Pumas results against the PBPK truth and finalize the numbers in
the methodology.
