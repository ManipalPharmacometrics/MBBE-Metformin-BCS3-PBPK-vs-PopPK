# How to read the Type I error results — plain-language guide

*Companion to `METHODOLOGY_metformin_MBBE.md` (safe space) and the Type I error
overlay produced by `scripts/05_type1_overlay.jl` / `06_type1_summary.jl`.
Results page: `TYPE1_ERROR_RESULTS.html` in this folder.*

**First, a naming clarification.** In the results table the numbers **8.1 and 9.6
are percentages** (the true type I error). The **ratios** are the 1.66x and 1.85x
sitting next to them. These are two different columns and it is easy to conflate them.

---

## 1. What "the boundary" is

Slower release (bigger T50) means lower Cmax, so as T50 increases the GMR falls.
Somewhere there is a T50 where the true GMR is **exactly 0.80** — sitting dead on
the BE limit. That T50 is the *simulated BE boundary*. PBPK puts it at **4.623 h**.

A tablet at exactly that T50 is a knife-edge product: not quite failing, not quite
passing.

## 2. Why 5% shows up there on its own

Type I error is the chance a study concludes "bioequivalent" when the truth sits
exactly on the limit.

For TOST this is **5% by construction**. At a true GMR of 0.80 you pass only if the
lower confidence bound climbs above 0.80, which requires the estimate to land about
1.68 standard errors above the truth. That is a 1-in-20 fluke — regardless of how
variable the drug is, and regardless of how large the study is.

This is why both methods return 4.8–5.4% at their own boundary, and why that
agreement **is not evidence that either model is correct**. It is arithmetic. Both
would return 5% even if one of the models were badly wrong. The scaled-reference
control (reference x 0.80) confirms the machinery: 5.15% at n=48, 5.39% at n=86.

## 3. What the 8.1% and 9.6% actually are

PopPK draws its line at **T50 = 4.574 h** and reports "type I error here is 4.9%".

But the PBPK truth says that at 4.574 h the true GMR is **0.8049**, not 0.8000 —
a hair *inside* the acceptance range. So a tablet there is, truthfully, just barely
bioequivalent.

The real chance a study at that T50 passes is therefore not 5%. It is **8.1%** at
n=48 and **9.6%** at n=86.

> **8.1% / 9.6% = the actual pass rate of a tablet sitting exactly where PopPK drew
> its line.** The ratio 1.66x is simply 8.1 / 4.9.

### What the risk actually is (and is not)

The product at that T50 is genuinely bioequivalent (GMR 0.805 is inside
0.80–1.25), so nothing unsafe is being waved through. What is wrong is the
**reported operating characteristic**: you would hand a regulator "5%" for a limit
that really carries 8–10%. Morris Part 2's point is that this number is fragile,
which makes it a shaky quantity to put in a submission.

## 4. Why the ratio is larger at n = 86

Picture the acceptance curve as a **cliff face**:

- The cliff is **pinned at 5% exactly at the limit**, for every study size. More
  subjects cannot help you there — you still need the 1-in-20 fluke.
- More subjects make the cliff **steeper**, because the confidence interval narrows.

The 5% anchor does not move, but the wall beside it gets steeper. Take one small
step inland from the edge — the 0.005 GMR offset — and a steeper cliff carries you
higher up:

| | SD of the log-ratio (SE) | pass rate at true GMR 0.8049 |
|---|---|---|
| n = 48 | 0.0175 | 8.1% |
| n = 86 | 0.0131 | 9.6% |

Same offset, narrower CI, more of it converts into a pass. With more subjects,
"genuinely just inside the limit" becomes "reliably demonstrable".

### The check that this is really about steepness, not about n making passing easier

When a boundary errs the *other* way — landing outside the limit, where the product
is truly not bioequivalent — a larger n pushes the rate **down**, not up:

| Boundary | true GMR | n = 48 | n = 86 | direction |
|---|---|---|---|---|
| PopPK per-arm, upper (inside limit) | 0.805 | 8.1% | **9.6%** | up |
| PopPK joint, upper (outside) | 0.792 | 2.1% | **1.6%** | down |
| PopPK per-arm, lower (outside) | 1.269 | 1.4% | **0.7%** | down |
| PopPK joint, lower (outside) | 1.293 | 0.2% | **0.04%** | down |
| PopPK per-arm, AUC0-t upper (outside) | 0.797 | 3.6% | **3.5%** | down |

Every row moves **further from 5%** as n grows — upward if the boundary landed
inside the limit, downward if outside. That is exactly what a steeper cliff pinned
at the same point does. If a larger n simply "made passing easier", the lower four
rows would have risen too.

### Direct verification of the pin and the flip

Scanned with `scripts/08_export_zoom_oc.jl` — 20 000 trials x 5 seeds per point,
PBPK Cmax, so the crossing is not Monte-Carlo noise:

| T50 (h) | true GMR | where | n = 48 | n = 86 | effect of larger n |
|---|---|---|---|---|---|
| 4.550 | 0.8074 | inside | 10.0% | 12.4% | higher |
| 4.570 | 0.8054 | inside | 8.3% | 9.9% | higher |
| 4.590 | 0.8033 | inside | 6.9% | 7.7% | higher |
| 4.610 | 0.8013 | inside | 5.7% | 5.9% | higher |
| **4.623** | **0.8000** | **on the limit** | **5.0%** | **4.9%** | **none — pinned** |
| 4.630 | 0.7993 | outside | 4.6% | 4.5% | lower |
| 4.650 | 0.7973 | outside | 3.8% | 3.4% | lower |
| 4.710 | 0.7915 | outside | 1.9% | 1.3% | lower |
| 4.750 | 0.7876 | outside | 1.2% | 0.7% | lower |

The n=48 / n=86 ordering **flips exactly across the limit**, and both curves pass
through 5% at GMR 0.8000. That crossing is the whole answer in one line: n does
not move the 5% anchor, it only steepens the approach to it.

Raw scan: `output/type1_zoom_oc.csv`.

### Arithmetic check (for review)

Within-subject CV for Cmax is 12.2%, so the log-scale SD is
`sw = sqrt(log(0.122^2 + 1)) = 0.1216`.

- n=48: `SE = 0.1216/sqrt(48) = 0.01754`, `t(0.95, 47) = 1.678`
- n=86: `SE = 0.1216/sqrt(86) = 0.01311`, `t(0.95, 85) = 1.663`

The offset is `delta = log(0.8049/0.80) = 0.00611`. Approximate pass rate is
`P(Z >= t - delta/SE)`:

- n=48: `P(Z >= 1.678 - 0.348) = P(Z >= 1.330) = 9.2%`
- n=86: `P(Z >= 1.663 - 0.466) = P(Z >= 1.197) = 11.6%`

These analytic values (9.2%, 11.6%) run a little above the simulated ones
(8.1%, 9.6%) because the simulation bootstraps the 24 PK-Sim physiologies with
replacement, which adds variance the closed form ignores. The **mechanism and the
direction match**; treat the closed form as the explanation and the simulation as
the reported result.

## 5. The takeaway

A bigger, more expensive study offers **no protection** against a misplaced
boundary — it makes the consequence worse. So the accuracy of the CQA -> GMR
inversion matters **more** as studies get larger, not less.

That is the argument for letting the mechanistic model place the limit even in the
case where the sparse-data model reproduces the BE decision correctly — which is
exactly what the safe-space comparison found it does.

---

*Written 2026-09-09. Numbers from the licensed per-arm and joint Pumas fits;
20 000 replicate trials per type I error estimate.*
