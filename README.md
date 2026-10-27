# Model-Based Bioequivalence: Metformin XR (BCS III), PBPK vs PopPK

Do a PBPK model and a PopPK model place the bioequivalence (BE) safe space in the same place?
This repository compares the two routes for metformin extended release (XR), using virtual BE
trials across a range of in vivo release rates.

- **PBPK (the reference).** A PK-Sim model of metformin with a Weibull release function, swept over
  8 release times (T50 = 1, 2, 3, 3.5, 4, 5, 6, 7 h).
- **PopPK (the method under test).** A Weibull-input PopPK model in Pumas, fitted to clinical-style
  data sampled from the PBPK profiles with 15% assay error. Two variants: a per-arm fit and a joint
  (IVIVC-style) fit across formulations.
- **Output.** Test/reference GMRs for Cmax and AUC, the T50 range that passes BE (the safe space),
  and the type I error at the boundary each method declares.

## Summary of results

With rich sampling (13 samples per subject) the two routes agree on the safe space to within 0.15 h of
T50, although the PopPK model does not recover the release parameters themselves. With sparse sampling
(5 samples) they disagree: between-subject variability in release becomes non-estimable, the per-arm
fit moves the Cmax boundaries inward, and the joint fit overstates the safe space (1.50–4.75 h against
1.95–4.10 h from PBPK).

## Known limitation

The PK-Sim metformin model overpredicts absolute exposure. Its total clearance is close to GFR only
(about 150 mL/min against a clinical renal clearance near 500 mL/min) and bioavailability is 0.84
(clinical about 0.5–0.6), so Cmax is roughly 3-fold higher than clinical data. Active renal secretion
(OCT2/MATE1) is the likely missing process. The PBPK-vs-PopPK comparison is internally consistent,
because the PopPK data are simulated from the same PBPK profiles, but absolute concentrations should
not be read as clinical predictions.

## Layout

```
pksim-output/      PK-Sim exports: one workbook per T50, plus population runs in popsims/
pbpk-programs/     PBPK side: NCA, BE safe space, VBE trials, reports (Julia and R)
  poppk_pumas/     first PopPK pipeline and the PBPK-vs-PopPK profile figure
poppk-per-arm/     PopPK fitted per formulation
  poppk_pumas/     scripts/, src/, data/, run_all.jl
                   output/ (rich sampling), output_rich_se/ (SEs, shrinkage), output_sparse/
poppk-joint/       PopPK joint fit across formulations, same structure as poppk-per-arm/
docs/              methodology, analysis walkthrough, type I error explainer and results page, slides
environment/       Project.toml of the Julia environment the PopPK fits ran in
```

Start with `docs/METHODOLOGY_metformin_MBBE.md`, then `poppk-per-arm/START_HERE.md`.

## Running it

- Julia 1.11 and Quarto 1.8 were used.
- The PopPK fit step (`scripts/02_fit_poppk.jl`) needs a licensed Pumas installation. Every other step
  uses open-source Julia packages.
- Each `poppk_pumas/` folder has its own `Project.toml`. From inside it:

  ```bash
  julia --project=. -e 'using Pkg; Pkg.instantiate()'   # first time only
  julia --project=. run_all.jl
  ```

  For the sparse-sampling run, set `POPPK_SCHEDULE=sparse` and point `POPPK_OUT` at
  `output_sparse/`.
- Scripts locate each other relative to their own folder, so the repository can be cloned anywhere.
- All data in this repository are simulated. No clinical or patient-level data are included.

## Related work

- [MBBE-Felodipine-ER-BCS2-PBPK-vs-PopPK](https://github.com/ManipalPharmacometrics/MBBE-Felodipine-ER-BCS2-PBPK-vs-PopPK):
  the same comparison for a dissolution-limited (BCS II) extended-release product.
- [Uncertainty Propagation in MBBE](https://github.com/ManipalPharmacometrics/Uncertainity-Propogation-in-MBBE):
  how parameter uncertainty is carried into virtual BE trials.
- The method follows Morris et al., AAPS J 2026, Parts 1 and 2
  ([doi:10.1208/s12248-026-01213-2](https://doi.org/10.1208/s12248-026-01213-2),
  [doi:10.1208/s12248-026-01215-0](https://doi.org/10.1208/s12248-026-01215-0)).

## Author

Nikita Ramwani
