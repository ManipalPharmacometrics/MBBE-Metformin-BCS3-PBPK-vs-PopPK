# =============================================================================
# STEP 9 — one tested formulation: predict BE for all the others.
#
# In : results/step1_dataset.csv, data/pbpk_endpoints.csv,
#      results/perarm/step3_poppk_endpoints.csv, model_window.jl
# Out: results/one_test/step9_parameters.csv, step9_predicted_endpoints.csv,
#      step9_gmr.csv, step9_safe_space.csv, step9_boundaries.csv, fig_step9_gmr.png
#
# Scenario: a BE study has been run with the reference (T50 3 h) and ONE test
# formulation. Can a PopPK model fitted to those two arms predict the other 6
# formulations from their dissolution curves alone?
#   Tested formulation: 5 h (main case) and 2 h (check: does the answer depend
#   on which side was tested?).
# Model: model_window.jl — release = the known dissolution curve (not
# estimated); estimated = absorption window (tau) and disposition (CL, Vc).
# ka is FIXED: from oral data alone, the absorption scale (ka) cannot be told
# apart from CL and Vc (fitted freely, ka ran to its 0.01 bound and CL, Vc lost
# their meaning). Fixing ka is the standard remedy; it is fixed at 0.3, 1 and
# 3 1/h to show the predictions do not depend on the value chosen. The free-ka
# fit is kept for comparison. Method names: "test5_ka1" = tested 5 h, ka 1.
# Saturable uptake (model_window_sat.jl, one extra parameter Km) is fitted with
# ka fixed at 0.3 and 1: "test5_sat_ka0.3". It contains the linear model
# (Km -> infinity), so the two -2LL values can be compared directly.
# Gastric emptying (model_window_ge.jl, one extra parameter kge): released drug
# must leave the stomach before the window, ka fixed at 0.3 and 1:
# "test5_ge_ka0.3". kge -> infinity gives the linear model back.
# Prediction: each subject keeps the random effects of its reference arm; only
# T50 (the dissolution curve) changes. Same 24 subjects as PK-Sim (paired BE).
# Compared with PBPK (truth) and with the per-arm PopPK fits, which used data
# from all 8 formulations (the best PopPK can do when every arm is tested).
# =============================================================================

using Pumas, CSV, DataFrames, CairoMakie, NCA
include("settings.jl")
include("be_functions.jl")
include("model_window.jl")
include("model_window_sat.jl")
include("model_window_ge.jl")

const TESTED   = [5.0, 2.0]                  # the one tested formulation (+ the reference)
const KA_FIXED = [0.3, 1.0, 3.0]             # fixed ka values (1/h); "free" = estimated
const KA_FIXED_SAT = [0.3, 1.0]              # fixed ka values for the saturable model
const KA_FIXED_GE  = [0.3, 1.0]              # fixed ka values for the gastric-emptying model
const OUT   = joinpath(RESULTS, "one_test"); mkpath(OUT)
const PRED_TIMES = collect(0.0:0.05:24.0)

data = CSV.read(joinpath(RESULTS, "step1_dataset.csv"), DataFrame)
pbpk = CSV.read(joinpath(DATA, "pbpk_endpoints.csv"), DataFrame)
T50S = sort(unique(data.T50))
to_population(df) = read_pumas(df; id = :SUBJID, time = :time, amt = :amt, evid = :evid, cmt = :cmt,
                               observations = [:dv], covariates = [:T50, :IndividualId])

# ---- 1. fit on the reference + one tested formulation, predict all 8 ---------------
refrows = data[data.T50 .== REF_T50, :]
params = DataFrame(); predicted = Dict{String,DataFrame}()
for tested in TESTED, (kind, model, kas) in (("", poppk_window, ["free"; KA_FIXED]),
                                             ("sat_", poppk_window_sat, KA_FIXED_SAT),
                                             ("ge_", poppk_window_ge, KA_FIXED_GE)), ka in kas
    name  = "test$(Int(tested))_$(kind)ka$(ka isa Number ? (isinteger(ka) ? Int(ka) : ka) : ka)"
    train = to_population(data[in.(data.T50, Ref([REF_T50, tested])), :])
    @info "fitting $name"
    fitted = ka == "free" ? fit(model, train, init_params(model), FOCE()) :
             fit(model, train, merge(init_params(model), (tvka = ka,)), FOCE(); constantcoef = (:tvka,))
    est = coef(fitted)
    push!(params, merge((model = name, tested_T50 = tested, minus2LL = -2 * loglikelihood(fitted)),
                        (; (k => (v isa Number ? v : missing) for (k, v) in pairs(est))...),
                        (omega_tau = est.Ω[1, 1], omega_cl = est.Ω[2, 2], omega_vc = est.Ω[3, 3]));
          cols = :union)

    # random effects of each subject's reference arm
    ebe = Dict(subj.id => η for (subj, η) in zip(train, empirical_bayes(fitted)))
    prof = DataFrame()
    for t in T50S
        rows = copy(refrows); rows.T50 .= t                    # same subject, new formulation
        pop  = to_population(rows)
        for (subj, pred) in zip(pop, predict(model, pop, est, [ebe[s.id] for s in pop]; obstimes = PRED_TIMES))
            append!(prof, DataFrame(T50 = t, IndividualId = Int(subj.covariates(0.0).IndividualId),
                                    time = pred.time, ipred = pred.ipred.dv))
        end
    end
    # Cmax and AUC 0–24 h with Pumas NCA, as in step 3
    pop_nca = read_nca(prof; id = :IndividualId, time = :time, observations = :ipred, group = [:T50], verbose = false)
    nca = innerjoin(NCA.cmax(pop_nca), NCA.auc(pop_nca; auctype = :last, method = :linuplogdown), on = [:id, :T50])
    predicted[name] = DataFrame(T50 = nca.T50, IndividualId = parse.(Int, string.(nca.id)),
                                Cmax = nca.cmax, AUC_tEnd = nca.auc)
end
perarm = CSV.read(joinpath(RESULTS, "perarm", "step3_poppk_endpoints.csv"), DataFrame)
predicted["perarm_all8"] = perarm                       # benchmark: every formulation tested
CSV.write(joinpath(OUT, "step9_parameters.csv"), params)
CSV.write(joinpath(OUT, "step9_predicted_endpoints.csv"),
          reduce(vcat, [insertcols(d, 1, :method => m) for (m, d) in predicted]))

sources = ["PBPK" => pbpk; [m => predicted[m] for m in params.model]; "perarm_all8" => perarm]

# ---- 2. GMR at every T50: predicted vs PBPK ----------------------------------------
gmr = DataFrame()
for t in T50S, ep in ENDPOINTS
    row = (T50 = t, endpoint = ep)
    vals = [Symbol("GMR_", m) => gmr_at(pk, ep, t) for (m, pk) in sources]
    push!(gmr, merge(row, (; vals...)))
end
CSV.write(joinpath(OUT, "step9_gmr.csv"), gmr)

# ---- 3. safe space (joint ABE power >= 80%) ------------------------------------------
grid = collect(range(1.0, 7.0; step = T50_STEP_POWER))
function edges(p)
    lo = hi = argmin(abs.(grid .- REF_T50))
    while lo > 1 && p[lo-1] >= POWER_TARGET; lo -= 1; end
    while hi < length(grid) && p[hi+1] >= POWER_TARGET; hi += 1; end
    (grid[lo], grid[hi])
end
safe = DataFrame()
for n in N_SUBJ, (m, pk) in sources
    lo, hi = edges([pass_rate_at(pk, "Joint", t, n, N_TRIAL_POWER) for t in grid])
    push!(safe, (n_subj = n, method = m, T50_low = lo, T50_high = hi))
end
CSV.write(joinpath(OUT, "step9_safe_space.csv"), safe)

# ---- 4. boundaries and the TRUE type I error there ----------------------------------
fine = collect(range(1.0, 7.0; step = T50_STEP_TYPE1))
bnd = DataFrame()
for (m, pk) in sources, side in (:low, :high)
    hits = [first_crossing(fine, [gmr_at(pk, ep, t) for t in fine], side) for ep in ENDPOINTS]
    ok = [h for h in hits if !ismissing(h.t50)]
    b = isempty(ok) ? missing : ok[argmin([abs(h.t50 - REF_T50) for h in ok])].t50   # joint = first reached
    push!(bnd, (method = m, side = String(side), boundary_T50 = b,
                true_type1_n48 = ismissing(b) ? missing : pass_rate_at(pbpk, "Joint", b, 48, N_TRIAL_TYPE1),
                true_type1_n86 = ismissing(b) ? missing : pass_rate_at(pbpk, "Joint", b, 86, N_TRIAL_TYPE1));
          promote = true)
end
CSV.write(joinpath(OUT, "step9_boundaries.csv"), bnd)

# ---- figure: GMR along T50, one row per tested formulation ---------------------------
fig = Figure(size = (1100, 850))
for (i, tested) in enumerate(TESTED), (k, ep) in enumerate(ENDPOINTS)
    ax = Axis(fig[i, k]; title = "$ep GMR, fitted to T50 3 h + $(Int(tested)) h only", xlabel = "T50 (h)",
              ylabel = "GMR vs T50 3 h", xticks = T50S)
    hlines!(ax, [BE_LOWER, BE_UPPER]; color = :gray70, linestyle = :dash)
    vlines!(ax, setdiff(T50S, [REF_T50, tested]); color = (:gray, 0.15), linewidth = 12)   # predicted
    d = gmr[gmr.endpoint .== ep, :]
    scatterlines!(ax, d.T50, d.GMR_PBPK; color = :black, linewidth = 3, label = "PBPK (truth)")
    for (m, c) in zip(params.model[params.tested_T50 .== tested], [:orange, :royalblue, :seagreen, :purple, :red, :magenta, :brown, :teal])
        lines!(ax, d.T50, d[!, Symbol("GMR_", m)]; color = c, label = m)
    end
    k == 1 && axislegend(ax; position = :rt)
end
save(joinpath(OUT, "fig_step9_gmr.png"), fig)

println("Trained models:"); show(params; allcols = true); println()
println("\nGMR, predicted vs PBPK:"); show(gmr; allrows = true, allcols = true); println()
println("\nSafe space (joint):"); show(safe; allrows = true); println()
println("\nBoundaries and true (PBPK) type I error there:"); show(bnd; allrows = true, allcols = true); println()
