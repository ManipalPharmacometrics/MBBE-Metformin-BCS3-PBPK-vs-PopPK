# =============================================================================
# be_functions.jl — the bioequivalence (BE) calculations used in steps 4–6.
#
#   be_test             one paired crossover comparison: GMR, 90% CI, pass/fail
#   pass_rate           share of simulated BE trials that pass (one endpoint)
#   pass_rate_joint     share that pass on Cmax AND AUC together
#   values_at           per-subject endpoint at any T50 between the simulated ones
#   first_crossing      T50 where the GMR first leaves 0.80–1.25
#
# Design: PK-Sim simulates the SAME virtual subjects under every formulation,
# so test vs reference is a paired 2x2 crossover. BE is tested on the
# within-subject log differences (two one-sided tests = 90% CI in 0.80–1.25).
# =============================================================================

using Distributions, Random, Statistics

"Paired average BE for one endpoint: GMR, 90% CI and pass (CI inside 0.80–1.25)."
function be_test(test, ref)
    d  = log.(test) .- log.(ref)                 # within-subject log difference
    n  = length(d)
    se = std(d) / sqrt(n)
    tq = quantile(TDist(n - 1), 0.95)            # 90% two-sided CI
    gmr = exp(mean(d))
    lo  = exp(mean(d) - tq * se)
    hi  = exp(mean(d) + tq * se)
    (n = n, GMR = gmr, CI_lower = lo, CI_upper = hi, pass = lo >= BE_LOWER && hi <= BE_UPPER)
end

"""
    pass_rate(test, ref; n_subj, n_trial, cv)

Simulate `n_trial` BE trials. Each trial draws `n_subj` subjects (with
replacement) from the virtual population, adds within-subject noise with CV
`cv` to the test and the reference value, and runs `be_test`.
Returns the fraction of trials that pass: the power if the test is truly BE,
the type I error if the true GMR sits on 0.80 or 1.25.
"""
function pass_rate(test, ref; n_subj, n_trial, cv)
    rng = MersenneTwister(SEED_TRIALS)
    sw  = sqrt(log(cv^2 + 1))                    # CV -> log-scale SD
    N   = length(ref)
    npass = 0
    for _ in 1:n_trial
        idx = rand(rng, 1:N, n_subj)
        tv  = test[idx] .* exp.(randn(rng, n_subj) .* sw)
        rv  = ref[idx]  .* exp.(randn(rng, n_subj) .* sw)
        npass += be_test(tv, rv).pass
    end
    npass / n_trial
end

"As `pass_rate`, but a trial passes only if Cmax AND AUC both pass (same subjects)."
function pass_rate_joint(testC, refC, testA, refA; n_subj, n_trial, cvC, cvA)
    rng = MersenneTwister(SEED_TRIALS)
    swC = sqrt(log(cvC^2 + 1)); swA = sqrt(log(cvA^2 + 1))
    N   = length(refC)
    npass = 0
    for _ in 1:n_trial
        idx = rand(rng, 1:N, n_subj)
        tc = testC[idx] .* exp.(randn(rng, n_subj) .* swC)
        rc = refC[idx]  .* exp.(randn(rng, n_subj) .* swC)
        ta = testA[idx] .* exp.(randn(rng, n_subj) .* swA)
        ra = refA[idx]  .* exp.(randn(rng, n_subj) .* swA)
        npass += be_test(tc, rc).pass && be_test(ta, ra).pass
    end
    npass / n_trial
end

"""
    values_at(pk, endpoint, t50)

Per-subject `endpoint` at any T50. Between two simulated T50 values each
subject's value is interpolated linearly in log(endpoint) vs log(T50).
`pk` has columns T50, IndividualId, Cmax, AUC_tEnd. Returns a vector ordered
by IndividualId.
"""
function values_at(pk, endpoint, t50)
    anchors = sort(unique(pk.T50))
    col(t) = sort(pk[pk.T50 .== t, :], :IndividualId)[!, endpoint]
    t50 = clamp(t50, anchors[1], anchors[end])
    t50 in anchors && return col(t50)            # a simulated formulation: exact values
    k = searchsortedlast(anchors, t50)
    w = (log(t50) - log(anchors[k])) / (log(anchors[k+1]) - log(anchors[k]))
    exp.((1 - w) .* log.(col(anchors[k])) .+ w .* log.(col(anchors[k+1])))
end

"""
    first_crossing(grid, gmr, side)

Walk from the reference T50 down (`side = :low`) or up (`:high`) the grid and
return the first T50 where the GMR crosses 0.80 or 1.25 (linear interpolation
between grid points), and which limit it crossed. `missing` if it never does.
"""
function first_crossing(grid, gmr, side)
    i = argmin(abs.(grid .- REF_T50))
    step = side === :low ? -1 : 1
    while 1 <= i + step <= length(grid)
        j = i + step
        for (name, limit) in (("0.80", BE_LOWER), ("1.25", BE_UPPER))
            a, b = gmr[i], gmr[j]
            if (a - limit) * (b - limit) <= 0 && a != b
                w = (limit - a) / (b - a)
                return (t50 = grid[i] + w * (grid[j] - grid[i]), limit = name)
            end
        end
        i = j
    end
    (t50 = missing, limit = "not reached")
end

"""
    pass_rate_at(pk, endpoint, t50, n_subj, n_trial)

Pass rate of simulated BE trials for a test formulation with release time
`t50` against the reference. `endpoint` = "Cmax", "AUC_tEnd" or "Joint"
(both must pass). Within-subject CVs from settings.jl.
"""
function pass_rate_at(pk, endpoint, t50, n_subj, n_trial)
    if endpoint == "Joint"
        pass_rate_joint(values_at(pk, "Cmax", t50), values_at(pk, "Cmax", REF_T50),
                        values_at(pk, "AUC_tEnd", t50), values_at(pk, "AUC_tEnd", REF_T50);
                        n_subj, n_trial, cvC = WSV_CV["Cmax"], cvA = WSV_CV["AUC_tEnd"])
    else
        pass_rate(values_at(pk, endpoint, t50), values_at(pk, endpoint, REF_T50);
                  n_subj, n_trial, cv = WSV_CV[endpoint])
    end
end

"True (noise-free) GMR of `endpoint` at release time `t50`."
gmr_at(pk, endpoint, t50) = be_test(values_at(pk, endpoint, t50), values_at(pk, endpoint, REF_T50)).GMR
