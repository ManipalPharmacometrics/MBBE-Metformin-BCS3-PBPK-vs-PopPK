# =============================================================================
# STEP 4 — test/reference GMR with 90% CI at each formulation, PBPK vs PopPK.
#
# In : data/pbpk_endpoints.csv, results/<model>/step3_poppk_endpoints.csv
# Out: results/<model>/step4_gmr.csv
#      results/<model>/fig_step4_gmr.png
#
# Each formulation is compared with the reference (T50 3.0 h) in the same 24
# virtual subjects, with no within-subject noise added: this is the "true" GMR
# each method sees. The test is the standard Pumas one, `pumas_be`
# (Bioequivalence package): a 2x2 crossover ANOVA, 90% CI of the GMR.
# Odd-numbered subjects get sequence RT, even-numbered TR. The reference row
# (reference vs itself) is GMR 1 by definition and has no CI.
# =============================================================================

using CSV, DataFrames, CairoMakie, Bioequivalence
include("settings.jl")
include("be_functions.jl")

pbpk = CSV.read(joinpath(DATA, "pbpk_endpoints.csv"), DataFrame)
COLOR = Dict("PBPK" => :firebrick, "PopPK" => :royalblue)

for model_name in MODELS
    dir   = joinpath(RESULTS, model_name)
    poppk = CSV.read(joinpath(dir, "step3_poppk_endpoints.csv"), DataFrame)

    gmr = DataFrame()
    for (source, pk) in ("PBPK" => pbpk, "PopPK" => poppk), t50 in sort(unique(pk.T50)), ep in ENDPOINTS
        t50 == REF_T50 && continue
        test, ref = values_at(pk, ep, t50), values_at(pk, ep, REF_T50)
        n   = length(ref)
        seq = [isodd(i) ? "RT" : "TR" for i in 1:n]
        crossover = DataFrame(id = repeat(1:n, inner = 2), sequence = repeat(seq, inner = 2),
                              period = repeat([1, 2], n))
        crossover.value = [s[p] == 'T' ? test[i] : ref[i]
                           for (i, s, p) in zip(crossover.id, crossover.sequence, crossover.period)]
        r = pumas_be(crossover; endpoint = :value, sigdigits = 16).result[1, :]
        push!(gmr, (source = source, T50 = t50, endpoint = ep, n = n, GMR = r.GMR,
                    CI_lower = r.LB, CI_upper = r.UB, pass = r.LB >= BE_LOWER && r.UB <= BE_UPPER))
    end
    CSV.write(joinpath(dir, "step4_gmr.csv"), gmr)

    fig = Figure(size = (1100, 450))
    for (k, ep) in enumerate(ENDPOINTS)
        ax = Axis(fig[1, k]; title = "$ep GMR ($model_name PopPK)", xlabel = "T50 (h)",
                  ylabel = "GMR vs T50 3.0 h", xticks = sort(unique(pbpk.T50)))
        hlines!(ax, [BE_LOWER, BE_UPPER]; color = :gray60, linestyle = :dash)
        for source in ("PBPK", "PopPK")
            d = gmr[(gmr.source .== source) .& (gmr.endpoint .== ep), :]
            scatterlines!(ax, d.T50, d.GMR; color = COLOR[source], label = source)
            rangebars!(ax, d.T50, d.CI_lower, d.CI_upper; color = COLOR[source])
        end
        k == 1 && axislegend(ax; position = :rt)
    end
    save(joinpath(dir, "fig_step4_gmr.png"), fig)

    println("\n$model_name — GMR by formulation:")
    show(unstack(gmr[:, [:source, :T50, :endpoint, :GMR]], [:T50, :endpoint], :source, :GMR); allrows = true)
    println()
end
