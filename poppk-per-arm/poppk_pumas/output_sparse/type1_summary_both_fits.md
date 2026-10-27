# Type I error — PBPK vs PopPK (companion to the safe-space table)

Nominal alpha = 0.05. `bnd` = simulated BE boundary (T50, h) from inverting the CQA -> GMR relationship. `T1 reported` = the Type I error PopPK computes at its own boundary; `T1 true` = the PBPK truth's acceptance rate at that same boundary; `ratio` = true / reported.


## PopPK fit: per-arm

| Endpoint | Side | Limit PBPK | Limit PopPK | n | bnd PBPK | bnd PopPK | CQA err | true GMR at PopPK bnd | T1 reported | T1 true | ratio | comparable |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| AUC_inf | high | 1.25 | 0.8 | 48 | 5.36 | 6.43 | +20.0% | 1.406 | 4.9% | 1.8% | 0.37x | NO — different limits |
| AUC_inf | high | 1.25 | 0.8 | 86 | 5.36 | 6.43 | +20.0% | 1.406 | 5.0% | 0.5% | 0.09x | NO — different limits |
| AUC_tEnd | high | 0.8 | 0.8 | 48 | 5.18 | 5.60 | +8.1% | 0.769 | 5.0% | 0.1% | 0.02x | yes |
| AUC_tEnd | high | 0.8 | 0.8 | 86 | 5.18 | 5.60 | +8.1% | 0.769 | 5.4% | 0.0% | 0.00x | yes |
| Cmax | high | 0.8 | 0.8 | 48 | 4.62 | 4.50 | -2.7% | 0.813 | 5.1% | 15.4% | 3.02x | yes |
| Cmax | high | 0.8 | 0.8 | 86 | 4.62 | 4.50 | -2.7% | 0.813 | 5.3% | 21.4% | 4.01x | yes |
| Cmax | low | 1.25 | 1.25 | 48 | 1.46 | 1.64 | +12.4% | 1.218 | 4.8% | 24.2% | 5.05x | yes |
| Cmax | low | 1.25 | 1.25 | 86 | 1.46 | 1.64 | +12.4% | 1.218 | 4.9% | 36.9% | 7.56x | yes |
| Joint (Cmax & AUC_0-t) | high | 0.8 | 0.8 | 48 | 4.62 | 4.50 | -2.7% | 0.813 | 4.8% | 12.8% | 2.65x | yes |
| Joint (Cmax & AUC_0-t) | high | 0.8 | 0.8 | 86 | 4.62 | 4.50 | -2.7% | 0.813 | 5.0% | 20.6% | 4.10x | yes |
| Joint (Cmax & AUC_0-t) | low | 1.25 | 1.25 | 48 | 1.46 | 1.64 | +12.4% | 1.218 | 4.8% | 25.1% | 5.21x | yes |
| Joint (Cmax & AUC_0-t) | low | 1.25 | 1.25 | 86 | 1.46 | 1.64 | +12.4% | 1.218 | 5.0% | 36.9% | 7.46x | yes |

## PopPK fit: joint

| Endpoint | Side | Limit PBPK | Limit PopPK | n | bnd PBPK | bnd PopPK | CQA err | true GMR at PopPK bnd | T1 reported | T1 true | ratio | comparable |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Cmax | high | 0.8 | 0.8 | 48 | 4.62 | 5.78 | +25.0% | 0.699 | 4.9% | 0.0% | 0.00x | yes |
| Cmax | high | 0.8 | 0.8 | 86 | 4.62 | 5.78 | +25.0% | 0.699 | 5.2% | 0.0% | 0.00x | yes |
| Joint (Cmax & AUC_0-t) | high | 0.8 | 0.8 | 48 | 4.62 | 5.78 | +25.0% | 0.699 | 4.7% | 0.0% | 0.00x | yes |
| Joint (Cmax & AUC_0-t) | high | 0.8 | 0.8 | 86 | 4.62 | 5.78 | +25.0% | 0.699 | 5.1% | 0.0% | 0.00x | yes |
