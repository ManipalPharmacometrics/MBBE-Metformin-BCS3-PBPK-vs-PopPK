# =============================================================================
# STEP 8 — held-out test: can a PopPK model trained on 3 formulations predict
# the other 5?
#
# In : results/step1_dataset.csv, data/pbpk_endpoints.csv,
#      results/perarm/step3_poppk_endpoints.csv, model_joint.jl, model_joint_frel.jl
# Out: results/heldout/step8_parameters.csv        estimates of the trained models
#      results/heldout/step8_predicted_endpoints.csv  Cmax, AUC per subject at all 8 T50s
#      results/heldout/step8_gmr.csv               predicted vs PBPK GMR per T50
#      results/heldout/step8_safe_space.csv        safe space (joint, power >= 80%)
#      results/heldout/step8_boundaries.csv        BE boundaries and true type I error there
#      results/heldout/fig_step8_gmr.png
#
# Train on T50 = 1, 3, 7 h (the ends and the reference). Hold out 2, 3.5, 4, 5, 6 h.
# Methods compared:
#   joint      model_joint.jl (T50 on the release time only)
#   joint_frel model_joint_frel.jl (T50 on release time AND bioavailability)
#   interp3    no model: the per-arm PopPK endpoints at 1, 3, 7 h, interpolated
#              between them (what per-arm fitting gives with 3 formulations)
# Prediction for a held-out T50: each subject keeps the random effects estimated
# in its reference (3 h) arm; only T50 changes. Same 24 subjects as PK-Sim, so
# BE stays paired.
# Pass: the predicted GMRs, safe space and boundaries match PBPK, and the TRUE
# (PBPK) type I error at the predicted boundary stays near 5%.
# =============================================================================

using Pumas, CSV, DataFrames, CairoMakie
include("settings.jl")
include("be_functions.jl")
include("model_joint.jl")
include("model_joint_frel.jl")

const TRAIN = [1.0, 3.0, 7.0]
const OUT   = joinpath(RESULTS, "heldout"); mkpath(OUT)
const PRED_TIMES = collect(0.0:0.05:24.0)

data = CSV.read(joinpath(RESULTS, "step1_dataset.csv"), DataFrame)
pbpk = CSV.read(joinpath(DATA, "pbpk_endpoints.csv"), DataFrame)
T50S = sort(unique(data.T50))
to_population(df) = read_pumas(df; id = :SUBJID, time = :time, amt = :amt, evid = :evid, cmt = :cmt,
                               observations = [:dv], covariates = [:T50, :IndividualId])

"AUC by linear-up / log-down trapezoids (as step 3)."
function auc_lin_up_log_down(t, c)
    auc = 0.0
    for i in 1:length(t)-1
        dt, c1, c2 = t[i+1] - t[i], c[i], c[i+1]
        auc += (c1 > 0 && c2 > 0 && c2 < c1) ? dt * (c1 - c2) / (log(c1) - log(c2)) : dt * (c1 + c2) / 2
    end
    auc
end

# ---- 1. fit both models on the 3 training formulations, predict all 8 -------------
train = to_population(data[in.(data.T50, Ref(TRAIN)), :])
refrows = data[data.T50 .== REF_T50, :]
params = DataFrame(); predicted = Dict{String,DataFrame}()
for (name, model) in ("joint" => poppk_joint, "joint_frel" => poppk_joint_frel)
    @info "training $name on T50 = $TRAIN"
    fitted = fit(model, train, init_params(model), FOCE(); constantcoef = (:tvlag,))
    est = coef(fitted)
    push!(params, merge((model = name, minus2LL = -2 * loglikelihood(fitted)),
                        (; (k => (v isa Number ? v : missing) for (k, v) in pairs(est))...),
                        (omega_alpha = est.Ω[1, 1], omega_cl = est.Ω[2, 2], omega_vc = est.Ω[3, 3]));
          cols = :union)

    # random effects of each subject's reference arm
    ebe = Dict(subj.id => η for (subj, η) in zip(train, empirical_bayes(fitted)))
    ends = DataFrame()
    for t in T50S
        rows = copy(refrows); rows.T50 .= t                    # same subject, new formulation
        pop  = to_population(rows)
        for (subj, pred) in zip(pop, predict(model, pop, est, [ebe[s.id] for s in pop]; obstimes = PRED_TIMES))
            push!(ends, (T50 = t, IndividualId = Int(subj.covariates(0.0).IndividualId),
                         Cmax = maximum(pred.ipred.dv), AUC_tEnd = auc_lin_up_log_down(pred.time, pred.ipred.dv)))
        end
    end
    predicted[name] = ends
end
perarm = CSV.read(joinpath(RESULTS, "perarm", "step3_poppk_endpoints.csv"), DataFrame)
predicted["interp3"] = perarm[in.(perarm.T50, Ref(TRAIN)), :]
CSV.write(joinpath(OUT, "step8_parameters.csv"), params)
CSV.write(joinpath(OUT, "step8_predicted_endpoints.csv"),
          reduce(vcat, [insertcols(d, 1, :method => m) for (m, d) in predicted]))

sources = ["PBPK" => pbpk; [m => predicted[m] for m in ("joint", "joint_frel", "interp3")]]

# ---- 2. GMR at every T50: predicted vs PBPK ----------------------------------------
gmr = DataFrame()
for t in T50S, ep in ENDPOINTS
    row = (T50 = t, role = t in TRAIN ? "trained" : "held out", endpoint = ep)
    vals = [Symbol("GMR_", m) => gmr_at(pk, ep, t) for (m, pk) in sources]
    push!(gmr, merge(row, (; vals...)))
end
CSV.write(joinpath(OUT, "step8_gmr.csv"), gmr)

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
CSV.write(joinpath(OUT, "step8_safe_space.csv"), safe)

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
CSV.write(joinpath(OUT, "step8_boundaries.csv"), bnd)

# ---- figure -------------------------------------------------------------------------
COLOR = Dict("PBPK" => :black, "joint" => :darkorange, "joint_frel" => :royalblue, "interp3" => :gray60)
fig = Figure(size = (1100, 450))
for (k, ep) in enumerate(ENDPOINTS)
    ax = Axis(fig[1, k]; title = "$ep GMR, trained on T50 1, 3, 7 h", xlabel = "T50 (h)",
              ylabel = "GMR vs T50 3 h", xticks = T50S)
    hlines!(ax, [BE_LOWER, BE_UPPER]; color = :gray70, linestyle = :dash)
    vlines!(ax, setdiff(T50S, TRAIN); color = (:gray, 0.15), linewidth = 12)   # held-out formulations
    d = gmr[gmr.endpoint .== ep, :]
    for (m, _) in sources
        scatterlines!(ax, d.T50, d[!, Symbol("GMR_", m)]; color = COLOR[m], label = m)
    end
    k == 1 && axislegend(ax; position = :rt)
end
save(joinpath(OUT, "fig_step8_gmr.png"), fig)

println("Trained models:"); show(params; allcols = true); println()
println("\nGMR, predicted vs PBPK:"); show(gmr; allrows = true, allcols = true); println()
println("\nSafe space (joint):"); show(safe; allrows = true); println()
println("\nBoundaries and true (PBPK) type I error there:"); show(bnd; allrows = true, allcols = true); println()
