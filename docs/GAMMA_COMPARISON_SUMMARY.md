# γ (Weibull shape) — fixed vs estimated: comparison summary

> **Superseded / illustrative.** These numbers are from the license-free cross-check, which
> forces identical disposition on every arm of a subject and therefore *overstates* the
> release–T50 tracking. The subsequent **licensed Pumas joint fit gave a nearly flat release
> scale (IVIVC exponent ≈ 0.10)** and estimated a single shared γ ≈ 1.7 — see
> `PERARM_vs_JOINT_COMPARISON.md`, which is authoritative. Keep this file only as an
> illustration of *why the shape must be free*, not as a result.

**Question.** In the PopPK Weibull-release model, the release shape γ can either be **fixed
at 1.0** (to match the PBPK dissolution shape β = 1.0) or **estimated as one shared value**
across all formulations. What does each choice do to the IVIVC and to the safe space?

**Setup (both modes identical except γ).** Joint fit; disposition (CL, Vc) held
formulation-independent (shared per subject, estimated — *not* hard-fixed to numbers); only
the release scale α varies with T50. Reference T50 = 3.0 h.

> These numbers are from the license-free cross-check (open-Julia two-stage joint fit). The
> authoritative version is the Pumas NLME fit; run both γ modes there to finalize. The
> qualitative conclusion (below) is expected to carry over.

## 1. IVIVC — fitted release scale α vs PBPK T50

| PBPK T50 (h) | 1.0 | 2.0 | 3.0 | 3.5 | 4.0 | 5.0 | 6.0 | 7.0 |
|---|---|---|---|---|---|---|---|---|
| α, **γ estimated** (≈0.82) | 3.2 | 4.1 | 5.4 | 5.7 | 6.5 | 7.6 | 9.1 | 11.0 |
| α, **γ = 1.0 fixed** | 0.9 | 1.2 | 1.6 | 1.6 | 1.8 | 2.1 | 2.5 | 3.1 |

Both modes give a **monotonic** α vs T50 (an IVIVC either way). The absolute α differs
because α and γ trade off (a smaller shape is compensated by a different scale).

## 2. Safe space (power ≥ 80%), PBPK vs both modes

| Endpoint | n | PBPK | PopPK **γ estimated** | PopPK **γ = 1.0** |
|---|---|---|---|---|
| Cmax | 48 | 2.00–4.15 | **1.70–4.60** | 1.00–6.40 |
| Cmax | 86 | 1.80–4.25 | **1.50–4.80** | 1.00–6.60 |
| AUC₀–t | 48 | 1.00–4.60 | **1.00–6.15** | 1.00–7.00 |
| Joint (Cmax & AUC₀–t) | 48 | 1.95–4.10 | **1.65–4.55** | 1.00–6.35 |
| Joint (Cmax & AUC₀–t) | 86 | 1.85–4.25 | **1.50–4.80** | 1.00–6.60 |

Figures: `programs/poppk_pumas/output/fig_gamma_comparison_*.png`
(black = PBPK, blue = γ estimated, red = γ = 1.0).

## 3. Conclusion

- **γ estimated (shared) is the better model** — it reproduces the PBPK safe space
  reasonably closely (Joint 1.65–4.55 vs PBPK 1.95–4.10).
- **γ = 1.0 is too rigid** — it cannot bend the concentration profile enough as T50 changes,
  so the test/reference ratio moves too slowly and the safe space balloons (Joint 1.0–6.35),
  a poor match to PBPK.
- **Why:** PBPK β = 1.0 is the *in-vitro dissolution* shape. The PopPK γ is the *in-vivo
  input* shape — dissolution folded with gut absorption and transit — which is not 1.
  Forcing them equal imposes a wrong structural assumption.
- **Both are kept** (`02_fit_poppk.jl` = γ estimated baseline; set γ constant to reproduce
  the γ = 1.0 arm) and compared here. **Recommendation: report γ estimated as the baseline**,
  and show γ = 1.0 as the sensitivity that demonstrates why the shape must be free.

*Reproduce:* `julia --project=nopumas_env scripts/05_gamma_comparison.jl` (license-free), or
run the Pumas fit both ways. Authored 2026-07-22.
