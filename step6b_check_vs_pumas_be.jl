# =============================================================================
# STEP 6b — cross-check: does our fast BE test give the same pass/fail as the
# standard Pumas test, `pumas_be`?
#
# In : data/pbpk_endpoints.csv, results/<model>/step3_poppk_endpoints.csv,
#      results/<model>/step6_boundaries.csv
# Out: results/step6b_check_vs_pumas_be.csv
#
# Why: steps 5–6 simulate hundreds of thousands of BE trials, so they use the
# fast paired test `be_test` (be_functions.jl). Step 4 uses `pumas_be`
# (2x2 crossover ANOVA). This step shows the two reach the same decisions.
#
# Where: at each BE boundary found in step 6. There the true GMR sits on
# 0.80 or 1.25, so trials are as close to the pass/fail line as they get, and
# any difference between the two tests would show up most.
#
# How: the first N_TRIAL_CHECK trials of `pass_rate` (same seed, same draws),
# each one run through BOTH tests. Odd subjects get sequence RT, even TR.
# PBPK rows are the same for both PopPK models, so they are listed once.
# =============================================================================

using CSV, DataFrames, Random, Statistics, Bioequivalence
include("settings.jl")
include("be_functions.jl")

pbpk = CSV.read(joinpath(DATA, "pbpk_endpoints.csv"), DataFrame)
out  = DataFrame()

for model_name in MODELS
    poppk = CSV.read(joinpath(RESULTS, model_name, "step3_poppk_endpoints.csv"), DataFrame)
    bnd   = CSV.read(joinpath(RESULTS, model_name, "step6_boundaries.csv"), DataFrame)
    for row in eachrow(bnd), n in N_SUBJ
        (row.endpoint in ENDPOINTS && !ismissing(row.T50)) || continue
        row.source == "PBPK" && model_name != first(MODELS) && continue
        pk = row.source == "PBPK" ? pbpk : poppk

        test, ref = values_at(pk, row.endpoint, row.T50), values_at(pk, row.endpoint, REF_T50)
        sw  = sqrt(log(WSV_CV[row.endpoint]^2 + 1))
        rng = MersenneTwister(SEED_TRIALS)
        seq = [isodd(i) ? "RT" : "TR" for i in 1:n]
        crossover = DataFrame(id = repeat(1:n, inner = 2), sequence = repeat(seq, inner = 2),
                              period = repeat([1, 2], n))
        pass_own, pass_pumas = falses(N_TRIAL_CHECK), falses(N_TRIAL_CHECK)
        for k in 1:N_TRIAL_CHECK
            idx = rand(rng, 1:length(ref), n)                       # as in pass_rate
            tv  = test[idx] .* exp.(randn(rng, n) .* sw)
            rv  = ref[idx]  .* exp.(randn(rng, n) .* sw)
            pass_own[k] = be_test(tv, rv).pass
            crossover.value = [s[p] == 'T' ? tv[i] : rv[i]
                               for (i, s, p) in zip(crossover.id, crossover.sequence, crossover.period)]
            r = pumas_be(crossover; endpoint = :value, sigdigits = 16).result[1, :]
            pass_pumas[k] = r.LB >= BE_LOWER && r.UB <= BE_UPPER
        end
        push!(out, (source = row.source, model = row.source == "PBPK" ? "both" : model_name,
                    endpoint = row.endpoint, side = row.side, T50 = row.T50, n_subj = n,
                    n_trial = N_TRIAL_CHECK, pass_rate_be_test = mean(pass_own),
                    pass_rate_pumas_be = mean(pass_pumas), same_decision = mean(pass_own .== pass_pumas)))
    end
end

CSV.write(joinpath(RESULTS, "step6b_check_vs_pumas_be.csv"), out)
println("Pass rate at the BE boundary, our be_test vs pumas_be (same trials):")
show(out; allrows = true, allcols = true); println()
println("\nLowest agreement: ", round(100 * minimum(out.same_decision); digits = 2), "% of trials")
