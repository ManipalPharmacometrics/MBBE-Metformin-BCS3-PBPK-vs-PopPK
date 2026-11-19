# =============================================================================
# STEP 5 — power and the BE safe space, PBPK vs PopPK.
#
# In : data/pbpk_endpoints.csv, results/<model>/step3_poppk_endpoints.csv
# Out: results/<model>/step5_power.csv        power at every T50 on a 0.05 h grid
#      results/<model>/step5_safe_space.csv   T50 range where power >= 80%
#      results/<model>/fig_step5_power.png
#
# Power at a given T50 = share of 500 simulated BE trials (2x2 crossover,
# n = 48 or 86, within-subject CV from settings.jl) that pass.
# Between the 8 simulated formulations, each subject's Cmax/AUC is
# interpolated in log-log (values_at in be_functions.jl).
# Safe space = the run of T50 values around the reference where power >= 80%.
# =============================================================================

using CSV, DataFrames, CairoMakie
include("settings.jl")
include("be_functions.jl")

pbpk = CSV.read(joinpath(DATA, "pbpk_endpoints.csv"), DataFrame)
const GRID = collect(range(minimum(pbpk.T50), maximum(pbpk.T50); step = T50_STEP_POWER))
COLOR = Dict("PBPK" => :firebrick, "PopPK" => :royalblue)

"Edges of the run of grid points around the reference where power >= target."
function safe_space(grid, power)
    lo = hi = argmin(abs.(grid .- REF_T50))
    while lo > 1 && power[lo-1] >= POWER_TARGET; lo -= 1; end
    while hi < length(grid) && power[hi+1] >= POWER_TARGET; hi += 1; end
    (grid[lo], grid[hi])
end

for model_name in MODELS
    dir   = joinpath(RESULTS, model_name)
    poppk = CSV.read(joinpath(dir, "step3_poppk_endpoints.csv"), DataFrame)

    power = DataFrame(); safe = DataFrame()
    for ep in [ENDPOINTS; "Joint"], n in N_SUBJ, (source, pk) in ("PBPK" => pbpk, "PopPK" => poppk)
        p = [pass_rate_at(pk, ep, t, n, N_TRIAL_POWER) for t in GRID]
        append!(power, DataFrame(endpoint = ep, n_subj = n, source = source, T50 = GRID, power = p))
        lo, hi = safe_space(GRID, p)
        push!(safe, (endpoint = ep, n_subj = n, source = source, T50_low = lo, T50_high = hi))
    end
    CSV.write(joinpath(dir, "step5_power.csv"), power)
    CSV.write(joinpath(dir, "step5_safe_space.csv"), safe)

    fig = Figure(size = (1300, 800))
    for (i, ep) in enumerate([ENDPOINTS; "Joint"]), (j, n) in enumerate(N_SUBJ)
        ax = Axis(fig[i, j]; title = "$ep, n = $n ($model_name PopPK)", xlabel = "T50 (h)",
                  ylabel = "power", xticks = sort(unique(pbpk.T50)))
        hlines!(ax, [POWER_TARGET]; color = :gray60, linestyle = :dash)
        for source in ("PBPK", "PopPK")
            d = power[(power.endpoint .== ep) .& (power.n_subj .== n) .& (power.source .== source), :]
            lines!(ax, d.T50, d.power; color = COLOR[source], label = source)
        end
        i == 1 && j == 1 && axislegend(ax; position = :rt)
    end
    save(joinpath(dir, "fig_step5_power.png"), fig)

    println("\n$model_name — safe space (T50 h where power >= 80%):")
    show(safe; allrows = true); println()
end
