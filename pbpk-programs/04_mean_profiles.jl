# =============================================================================
# 04_mean_profiles.jl   (Julia primary pipeline -- Step 4: mean PK profiles)
#
# Builds the population MEAN plasma concentration-time profile for each
# formulation (T50 sweep) from the PK-Sim population Results exports, plus the
# mean fraction-of-dose-absorbed profile. These support the mechanistic report
# (flip-flop kinetics, AUC_inf extrapolation artefact, Tmax plateau, absorption
# window) -> see pk_profiles_report.qmd.
#
# Reads pksim-output/popsims/*-Results.csv. Each file:
#   col1 IndividualId | col2 Time [min] | col3 Peripheral venous plasma [µmol/l]
#   col4 venous-blood container conc | col5 Fraction of dose absorbed into mucosa
# Same 24 subjects, same 3-min time grid (0-24 h) across all formulations.
#
# Outputs under output/:
#   mean_profiles.csv          T50, time_h, conc_mean (µmol/L), fabs_mean
#   profile_diagnostics.csv    per-T50 mean-profile Cmax/Tmax/Clast/Fabs@24h
#   fig_mean_profiles.png      overlaid mean profiles (linear + semi-log)
#   fig_fraction_absorbed.png  overlaid mean fraction-absorbed vs time
#
# Run from programs/:  julia --project=. 04_mean_profiles.jl
# =============================================================================

using CSV, DataFrames, Statistics, CairoMakie
include(joinpath(@__DIR__, "julia", "be_utils.jl"))
using .BEUtils

const POPSIMS = joinpath(@__DIR__, "..", "pksim-output", "popsims")
const OUT     = joinpath(@__DIR__, "output")
mkpath(OUT)

parse_t50(p) = parse(Float64, match(r"T50_([0-9.]+)h", basename(p)).captures[1])

# --- Build mean profiles per formulation -------------------------------------
files = sort(filter(f -> endswith(f, "-Results.csv"), readdir(POPSIMS; join = true)),
             by = parse_t50)
prof = DataFrame()
for f in files
    df = CSV.read(f, DataFrame)
    rename!(df, names(df)[1] => :id, names(df)[2] => :tmin,
            names(df)[3] => :cpvb, names(df)[5] => :fabs)
    g = combine(groupby(df, :tmin),
                :cpvb => (x -> mean(skipmissing(x))) => :conc,
                :fabs => (x -> mean(skipmissing(x))) => :fabs)
    g.T50 .= parse_t50(f)
    g.time_h = g.tmin ./ 60
    append!(prof, select(g, [:T50, :time_h, :conc, :fabs]))
end
sort!(prof, [:T50, :time_h])
CSV.write(joinpath(OUT, "mean_profiles.csv"), prof)

# --- Per-formulation diagnostics from the MEAN profile -----------------------
diag = DataFrame()
for t50 in sort(unique(prof.T50))
    d = sort(filter(:T50 => ==(t50), prof), :time_h)
    imax = argmax(d.conc)
    push!(diag, (T50 = t50, Cmax = d.conc[imax], Tmax = d.time_h[imax],
                 Clast_24h = d.conc[end], Fabs_24h = d.fabs[end]))
end
CSV.write(joinpath(OUT, "profile_diagnostics.csv"), diag)

# --- Plots -------------------------------------------------------------------
t50s = sort(unique(prof.T50))
cols = cgrad(:viridis, length(t50s); categorical = true)
colof = Dict(t => cols[i] for (i, t) in enumerate(t50s))

fig = Figure(size = (1150, 560), fontsize = 16)
Label(fig[0, 1:2], "Metformin XR — population mean plasma profiles by dissolution T50";
      fontsize = 18, font = :bold)
axL = Axis(fig[1, 1], xlabel = "Time (h)", ylabel = "Mean plasma conc (µmol/L)",
           title = "Linear", xlabelsize = 17, ylabelsize = 17,
           xticks = 0:4:24)
axS = Axis(fig[1, 2], xlabel = "Time (h)", ylabel = "Mean plasma conc (µmol/L)",
           title = "Semi-log (terminal slope = flip-flop)", yscale = log10,
           xlabelsize = 17, ylabelsize = 17, xticks = 0:4:24)
for t in t50s
    d = sort(filter(:T50 => ==(t), prof), :time_h)
    lbl = t == REF_T50 ? "T50 $(t) h (ref)" : "T50 $(t) h"
    lw  = t == REF_T50 ? 3.5 : 2
    lines!(axL, d.time_h, d.conc; color = colof[t], linewidth = lw, label = lbl)
    pos = d.conc .> 0
    lines!(axS, d.time_h[pos], d.conc[pos]; color = colof[t], linewidth = lw)
end
ylims!(axS, 1, 40)   # focus the log axis on the informative range
Legend(fig[1, 3], axL, "Formulation"; framevisible = false)
save(joinpath(OUT, "fig_mean_profiles.png"), fig, px_per_unit = 2)

# Fraction absorbed over time.
fig2 = Figure(size = (760, 560), fontsize = 16)
ax2 = Axis(fig2[1, 1], xlabel = "Time (h)", ylabel = "Mean fraction of dose absorbed",
           title = "Metformin XR — fraction absorbed vs time by T50",
           xlabelsize = 17, ylabelsize = 17, xticks = 0:4:24)
for t in t50s
    d = sort(filter(:T50 => ==(t), prof), :time_h)
    lines!(ax2, d.time_h, d.fabs; color = colof[t], linewidth = (t == REF_T50 ? 3.5 : 2),
           label = t == REF_T50 ? "T50 $(t) h (ref)" : "T50 $(t) h")
end
hlines!(ax2, [1.0]; linestyle = :dash, color = :gray50)
ylims!(ax2, 0, 1.02)
Legend(fig2[1, 2], ax2, "Formulation"; framevisible = false)
save(joinpath(OUT, "fig_fraction_absorbed.png"), fig2, px_per_unit = 2)

println("Mean-profile diagnostics:")
show(diag; allrows = true, allcols = true); println()
println("Wrote mean_profiles.csv, profile_diagnostics.csv, and 2 figures to ", normpath(OUT))
