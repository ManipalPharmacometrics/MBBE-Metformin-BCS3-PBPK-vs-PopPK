# =============================================================================
# STEP 3 — BE endpoints (Cmax, AUC 0–24 h) per subject from the PopPK
# individual predictions, and a side-by-side check against PK-Sim.
#
# In : results/<model>/step2_ipred.csv, data/pbpk_endpoints.csv
# Out: results/<model>/step3_poppk_endpoints.csv   T50, IndividualId, Cmax, AUC_tEnd
#      results/<model>/step3_check_vs_pbpk.csv     geometric means by formulation,
#                                                  PopPK / PBPK ratio
#
#   Cmax     = NCA.cmax: highest predicted concentration (0.05 h grid)
#   AUC_tEnd = NCA.auc (auctype :last, method :linuplogdown): AUC 0–24 h,
#              linear trapezoid while rising, log trapezoid while falling
#
# Units: Cmax µmol/L. AUC: PopPK µmol·h/L; PK-Sim (data/pbpk_endpoints.csv)
# µmol·min/L. Units cancel in every test/reference ratio, so steps 4–6 use
# both files as they are; the check table below converts PK-Sim AUC to hours.
# =============================================================================

using CSV, DataFrames, Statistics, NCA
include("settings.jl")

geomean(x) = exp(mean(log.(x)))
pbpk = CSV.read(joinpath(DATA, "pbpk_endpoints.csv"), DataFrame)

for model_name in MODELS
    dir   = joinpath(RESULTS, model_name)
    ipred = CSV.read(joinpath(dir, "step2_ipred.csv"), DataFrame)

    # Pumas NCA package: one NCA subject per (T50, IndividualId)
    pop   = read_nca(ipred; id = :IndividualId, time = :time, observations = :ipred,
                     group = [:T50], verbose = false)
    nca   = innerjoin(NCA.cmax(pop), NCA.auc(pop; auctype = :last, method = :linuplogdown),
                      on = [:id, :T50])
    poppk = sort(DataFrame(T50 = nca.T50, IndividualId = parse.(Int, string.(nca.id)),
                           Cmax = nca.cmax, AUC_tEnd = nca.auc), [:T50, :IndividualId])
    CSV.write(joinpath(dir, "step3_poppk_endpoints.csv"), poppk)

    a = combine(groupby(poppk, :T50), :Cmax => geomean => :PopPK_Cmax, :AUC_tEnd => geomean => :PopPK_AUC)
    b = combine(groupby(pbpk, :T50), :Cmax => geomean => :PBPK_Cmax,
                :AUC_tEnd => (x -> geomean(x) / 60) => :PBPK_AUC)     # min -> h
    check = innerjoin(b, a, on = :T50)
    check.ratio_Cmax = check.PopPK_Cmax ./ check.PBPK_Cmax
    check.ratio_AUC  = check.PopPK_AUC ./ check.PBPK_AUC
    CSV.write(joinpath(dir, "step3_check_vs_pbpk.csv"), check)
    println("\n$model_name — geometric means (Cmax µmol/L, AUC 0–24 h µmol·h/L):")
    show(check; allcols = true); println()
end
