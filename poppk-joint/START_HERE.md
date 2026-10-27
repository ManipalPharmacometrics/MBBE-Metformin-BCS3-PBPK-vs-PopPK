# START HERE — JOINT (IVIVC) PopPK run

This folder runs the **joint, disposition-constrained** PopPK fit — the variant that gives
an **IVIVC** (the fitted in-vivo release scale α tracking T50). It is meant to be run on your
**licensed Pumas / JuliaHub** environment and the outputs sent back, so we can compare it
against the per-arm run you already did.

## ▶ The file to run

```bash
cd poppk_pumas
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # first time only
julia --project=. run_all.jl                          # <-- RUN THIS  (01 → 02 → 03 → 04)
```

## ▶ What makes this the JOINT run (vs your per-arm run)

Step 2 here is **one fit over all 8 formulation arms at once** (not eight separate fits):

- **CL, Vc are shared** across all tablets — formulation-independent (a tablet's release rate
  doesn't change how the body clears the drug). This removes the 1-cmt flip-flop ambiguity
  that made α flat in the per-arm run.
- **γ (release shape) is one shared, estimated value.**
- **Only the release scale α varies with formulation**, through a T50 covariate:
  `Alpha(T50) = tvalpha · (T50/3)^slope · exp(η)`. The fitted **`slope`** is the IVIVC exponent.

So the population table (`poppk_param_estimates.csv`) will have a single row with `tvalpha`,
**`slope`**, `tvgamma`, `tvcl`, `tvvc`; and the printout ends with an **IVIVC check** —
geometric-mean fitted α vs T50 (α should now increase with T50).

All the compatibility fixes from your per-arm run are already baked in (`ifelse` dynamics,
`tvlag` `init=1e-3` held via `constantcoef=(:tvlag,)`, callable `icoef(0.0).Alpha`).

## ▶ What to send back

Everything lands in `poppk_pumas/output/`. Send that folder back — or at minimum:

| File | Why |
|---|---|
| `poppk_pumas/output/poppk_param_estimates.csv` | the `slope` + shared CL/Vc/γ (the IVIVC) |
| `poppk_pumas/output/poppk_icoef.csv` | per-subject α/CL/Vc (α-vs-T50) |
| `poppk_pumas/output/poppk_ipred_profiles.csv` | profiles — 03/04 rebuild everything from this |
| `poppk_pumas/output/safe_space_overlay_boundaries.csv` | joint safe space vs PBPK |
| `poppk_pumas/output/be_overlay_gmr.csv` | GMR(T50) both methods |

Expect the joint safe space to differ somewhat from your per-arm run — that difference (and
the recovered IVIVC) is exactly what we want to compare. I'll overlay the two here.

## ▶ If Pumas errors (only step 02 is version-sensitive)

- `@covariates T50` must be recognized (used in `@pre` for the release scale).
- `icoef(fpm)` entries read as callables: `ico(0.0).Alpha`.
- `predict(fpm; obstimes=…)` read as `pr.ipred.dv`, `pr.time`.
- If the lag discontinuity destabilises the gradient, `tvlag` is already pinned at 1e-3 via
  `constantcoef`; leave it.

## ▶ Files here

| Path | What |
|---|---|
| `poppk_pumas/run_all.jl` | **entry point — run this** |
| `poppk_pumas/src/poppk_model.jl` | the **joint** Weibull model (shared CL/Vc/γ, T50 covariate on α) |
| `poppk_pumas/scripts/02_fit_poppk.jl` | the joint Pumas fit |
| `poppk_pumas/scripts/{01,03,04}_*.jl` | dataset build, NCA, PBPK-vs-PopPK overlay |
| `poppk_pumas/data/` | the 8 PBPK profiles + PBPK endpoints (inputs) |
