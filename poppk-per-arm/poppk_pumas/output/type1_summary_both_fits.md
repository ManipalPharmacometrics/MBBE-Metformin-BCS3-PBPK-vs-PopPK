# Type I error — PBPK vs PopPK (companion to the safe-space table)

Nominal alpha = 0.05. `bnd` = simulated BE boundary (T50, h) from inverting the CQA -> GMR relationship. `T1 reported` = the Type I error PopPK computes at its own boundary; `T1 true` = the PBPK truth's acceptance rate at that same boundary; `ratio` = true / reported.


## PopPK fit: per-arm

| Endpoint | Side | Limit PBPK | Limit PopPK | n | bnd PBPK | bnd PopPK | CQA err | true GMR at PopPK bnd | T1 reported | T1 true | ratio | comparable |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| AUC_inf | high | 1.25 | 0.8 | 48 | 5.36 | 6.16 | +14.9% | 1.245 | 5.0% | 9.7% | 1.94x | NO — different limits |
| AUC_inf | high | 1.25 | 0.8 | 86 | 5.36 | 6.16 | +14.9% | 1.245 | 5.0% | 8.6% | 1.73x | NO — different limits |
| AUC_tEnd | high | 0.8 | 0.8 | 48 | 5.18 | 5.22 | +0.7% | 0.797 | 5.0% | 3.6% | 0.72x | yes |
| AUC_tEnd | high | 0.8 | 0.8 | 86 | 5.18 | 5.22 | +0.7% | 0.797 | 5.2% | 3.5% | 0.67x | yes |
| Cmax | high | 0.8 | 0.8 | 48 | 4.62 | 4.57 | -1.1% | 0.805 | 4.9% | 8.1% | 1.66x | yes |
| Cmax | high | 0.8 | 0.8 | 86 | 4.62 | 4.57 | -1.1% | 0.805 | 5.2% | 9.6% | 1.85x | yes |
| Cmax | low | 1.25 | 1.25 | 48 | 1.46 | 1.36 | -6.8% | 1.269 | 5.0% | 1.4% | 0.28x | yes |
| Cmax | low | 1.25 | 1.25 | 86 | 1.46 | 1.36 | -6.8% | 1.269 | 5.1% | 0.7% | 0.14x | yes |
| Joint (Cmax & AUC_0-t) | high | 0.8 | 0.8 | 48 | 4.62 | 4.57 | -1.1% | 0.805 | 3.8% | 6.1% | 1.60x | yes |
| Joint (Cmax & AUC_0-t) | high | 0.8 | 0.8 | 86 | 4.62 | 4.57 | -1.1% | 0.805 | 4.9% | 8.9% | 1.81x | yes |
| Joint (Cmax & AUC_0-t) | low | 1.25 | 1.25 | 48 | 1.46 | 1.36 | -6.8% | 1.269 | 5.2% | 1.3% | 0.26x | yes |
| Joint (Cmax & AUC_0-t) | low | 1.25 | 1.25 | 86 | 1.46 | 1.36 | -6.8% | 1.269 | 5.2% | 0.8% | 0.15x | yes |

## PopPK fit: joint

| Endpoint | Side | Limit PBPK | Limit PopPK | n | bnd PBPK | bnd PopPK | CQA err | true GMR at PopPK bnd | T1 reported | T1 true | ratio | comparable |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| AUC_inf | high | 1.25 | 0.8 | 48 | 5.36 | 6.55 | +22.2% | 1.481 | 4.8% | 0.9% | 0.19x | NO — different limits |
| AUC_inf | high | 1.25 | 0.8 | 86 | 5.36 | 6.55 | +22.2% | 1.481 | 5.0% | 0.1% | 0.03x | NO — different limits |
| AUC_tEnd | high | 0.8 | 0.8 | 48 | 5.18 | 5.53 | +6.6% | 0.774 | 4.9% | 0.2% | 0.03x | yes |
| AUC_tEnd | high | 0.8 | 0.8 | 86 | 5.18 | 5.53 | +6.6% | 0.774 | 5.2% | 0.0% | 0.00x | yes |
| Cmax | high | 0.8 | 0.8 | 48 | 4.62 | 4.71 | +1.8% | 0.792 | 4.9% | 2.1% | 0.44x | yes |
| Cmax | high | 0.8 | 0.8 | 86 | 4.62 | 4.71 | +1.8% | 0.792 | 5.2% | 1.6% | 0.30x | yes |
| Cmax | low | 1.25 | 1.25 | 48 | 1.46 | 1.25 | -14.2% | 1.293 | 5.0% | 0.2% | 0.03x | yes |
| Cmax | low | 1.25 | 1.25 | 86 | 1.46 | 1.25 | -14.2% | 1.293 | 5.2% | 0.0% | 0.01x | yes |
| Joint (Cmax & AUC_0-t) | high | 0.8 | 0.8 | 48 | 4.62 | 4.71 | +1.8% | 0.792 | 4.0% | 1.2% | 0.31x | yes |
| Joint (Cmax & AUC_0-t) | high | 0.8 | 0.8 | 86 | 4.62 | 4.71 | +1.8% | 0.792 | 4.9% | 1.2% | 0.24x | yes |
| Joint (Cmax & AUC_0-t) | low | 1.25 | 1.25 | 48 | 1.46 | 1.25 | -14.2% | 1.293 | 5.2% | 0.2% | 0.04x | yes |
| Joint (Cmax & AUC_0-t) | low | 1.25 | 1.25 | 86 | 1.46 | 1.25 | -14.2% | 1.293 | 5.2% | 0.0% | 0.01x | yes |
