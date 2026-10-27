# Per-arm vs joint PopPK — both licensed-Pumas runs compared

Two PopPK fits were run on the licensed Pumas host and compared against the PBPK truth:

- **Per-arm** — each of the 8 formulations fit independently (own CL, Vc, γ, α).
- **Joint (IVIVC variant)** — one fit over all 192 subject-arms; typical CL, Vc, γ **shared**;
  release scale via a T50 covariate `Alpha = tvalpha·(T50/3)^slope·exp(η)`. Converged, OFV ≈ 6093.6,
  fitted **slope = 0.10**, shared CL ≈ 12.1 L/h, Vc ≈ 129 L, γ ≈ 1.71.

## 1. Headline — decision equivalence holds in BOTH runs

Safe-space boundaries (power ≥ 80%), T50 in hours:

| Endpoint | n | PBPK | PopPK **per-arm** | PopPK **joint** |
|---|---|---|---|---|
| Cmax | 48 | 2.00–4.15 | 1.95–4.05 | 1.90–4.10 |
| Cmax | 86 | 1.80–4.25 | 1.75–4.15 | 1.70–4.25 |
| AUC₀–t | 48 | 1.00–4.60 | 1.00–4.55 | 1.00–4.75 |
| AUC₀–t | 86 | 1.00–4.70 | 1.00–4.65 | 1.00–4.90 |
| AUC_inf | 48 | 1.00–4.95 | 1.00–4.95 | 1.00–5.45 |
| AUC_inf | 86 | 1.00–5.05 | 1.00–5.20 | 1.00–5.80 |
| **Joint (Cmax & AUC₀–t)** | 48 | **1.95–4.10** | **1.90–4.05** | **1.85–4.10** |
| **Joint (Cmax & AUC₀–t)** | 86 | **1.85–4.25** | **1.75–4.15** | **1.70–4.25** |

Cmax (the binding endpoint) and the joint endpoint agree with PBPK to ~0.05–0.15 h in both
fits. AUC_inf runs a little wider than PBPK (more so in the joint fit) — expected: AUC_inf
extrapolation is unstable for the slow, flip-flop arms. **Decision equivalence is robust to
the fitting strategy.**

## 2. IVIVC — the key finding: it does NOT appear (either fit)

Geometric-mean fitted release scale α vs PBPK T50:

| PBPK T50 (h) | 1.0 | 2.0 | 3.0 | 3.5 | 4.0 | 5.0 | 6.0 | 7.0 |
|---|---|---|---|---|---|---|---|---|
| α, per-arm fit | 1.42 | 1.60 | 1.69 | 1.65 | 1.68 | 1.65 | 1.69 | 1.71 |
| α, joint fit (slope 0.10) | 1.43 | 1.57 | 1.68 | 1.67 | 1.72 | 1.71 | 1.70 | 1.76 |

**In both fits α is essentially flat** — it does *not* track T50. The joint model was built
specifically to force the IVIVC (by sharing the *typical* disposition), and it still didn't:
the fitted IVIVC exponent came out ≈ 0.

**Why.** In the joint fit the per-subject clearance and volume still carry random effects
(IIV), and those drift systematically with formulation — the fitted per-subject values rise
**CL 10.6 → 14.7 L/h and Vc 86 → 190 L across T50 = 1 → 7 h**. So even with the *typical*
disposition shared, the formulation's effect was absorbed into the *random* disposition,
leaving α flat. This is the one-compartment **flip-flop non-identifiability**: from a sparse
single-read-out BE study the model cannot tell "slow release + fast disposition" from "fast
release + slow disposition," and it defaults to explaining the change through disposition.

The only way we made α track T50 was the earlier license-free cross-check, which **forced the
same disposition on every arm of a subject** (no arm-to-arm disposition freedom at all) — i.e.
it effectively told the model the answer. The data themselves do not identify it. *(That
cross-check therefore overstated the IVIVC; the licensed NLME fits above are authoritative.)*

## 3. What this means

- **PopPK MBBE reproduces the BE *decision*** (safe space) of the mechanistic PBPK model —
  that is the deliverable, and it is confirmed twice.
- **PopPK MBBE does *not* recover the *mechanism*** (the in-vitro→in-vivo release mapping) from
  a sparse clinical-style study — the formulation effect is statistically absorbed into
  apparent disposition. Recovering the release mapping requires a **mechanistic model (PBPK)**,
  or external constraints on disposition that sparse BE data cannot themselves justify.
- This is a clean, honest contrast and arguably a *stronger* result than a forced IVIVC: it
  delineates exactly what PopPK can and cannot do relative to PBPK.

## 4. Recommendation

Report **decision equivalence** as the headline (per-arm baseline; joint as confirmation).
Report the **flat α / flip-flop non-identifiability** as a finding, not a failure — it is the
reason a mechanistic PBPK model adds value. No further IVIVC-chasing is warranted from this
1-compartment sparse-data setting.

*Sources: licensed-Pumas outputs in `poppk-per-arm/` (per-arm) and
`poppk-joint/` (joint). Authored 2026-07-22.*
