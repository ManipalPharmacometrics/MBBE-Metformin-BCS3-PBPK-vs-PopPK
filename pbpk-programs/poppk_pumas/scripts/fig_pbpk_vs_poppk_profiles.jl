# =============================================================================
# fig_pbpk_vs_poppk_profiles.jl  — overlay mean concentration-time profiles,
# PBPK (PK-Sim truth) vs PopPK (Pumas IPRED), one panel per formulation (T50).
#
# Uses the REAL licensed-Pumas per-arm IPRED as the PopPK profiles.
# License-free (CSV/DataFrames/Statistics/CairoMakie). Run:
#   julia --project=nopumas_env scripts/fig_pbpk_vs_poppk_profiles.jl
# =============================================================================

using CSV, DataFrames, Statistics, CairoMakie

const ROOT     = normpath(joinpath(@__DIR__, ".."))
const POPSIMS  = joinpath(ROOT, "data", "popsims")
const OUT      = joinpath(ROOT, "output"); mkpath(OUT)
# real licensed-Pumas per-arm IPRED (authoritative PopPK profiles)
const POPPK_IPRED = normpath(joinpath(@__DIR__, "..", "..", "..", "poppk-per-arm", "poppk_pumas", "output", "poppk_ipred_profiles.csv"))

const CONC = "Organism|PeripheralVenousBlood|Metformin|Plasma (Peripheral Venous Blood) [µmol/l]"
const TMIN = "Time [min]"
parse_t50(p) = parse(Float64, match(r"T50_([0-9.]+)h", basename(p)).captures[1])

# ---- PBPK mean profiles (per T50) -------------------------------------------
pbpk = Dict{Float64,Tuple{Vector{Float64},Vector{Float64}}}()
for f in sort(filter(x->endswith(x,"-Results.csv"), readdir(POPSIMS; join=true)), by=parse_t50)
    df = CSV.read(f, DataFrame)
    g  = combine(groupby(df, TMIN), CONC => mean => :c)
    sort!(g, TMIN)
    pbpk[parse_t50(f)] = (Float64.(g[!,TMIN]) ./ 60, Float64.(g.c))   # min -> h
end

# ---- PopPK mean profiles (per T50) from real Pumas IPRED --------------------
ip = CSV.read(POPPK_IPRED, DataFrame)
poppk = Dict{Float64,Tuple{Vector{Float64},Vector{Float64}}}()
for gk in groupby(ip, :T50)
    g = combine(groupby(gk, :time), :ipred => mean => :c); sort!(g, :time)
    poppk[gk.T50[1]] = (Float64.(g.time), Float64.(g.c))
end

# ---- figure: 2x4 grid, one panel per T50 ------------------------------------
t50s = sort(collect(keys(pbpk)))
fig = Figure(size = (1500, 720), fontsize = 15)
Label(fig[0, 1:4], "Mean plasma metformin — PBPK (PK-Sim) vs PopPK (Pumas IPRED), by formulation T50";
      fontsize = 19, font = :bold)
for (i, t50) in enumerate(t50s)
    r = (i-1) ÷ 4 + 1; c = (i-1) % 4 + 1
    ax = Axis(fig[r, c]; title = "T50 = $(t50) h",
              xlabel = r == 2 ? "Time (h)" : "", ylabel = c == 1 ? "Conc (µmol/L)" : "")
    tb, cb = pbpk[t50]; tp, cp = poppk[t50]
    lines!(ax, tb, cb; color = :black,     linewidth = 2.5, label = "PBPK (truth)")
    lines!(ax, tp, cp; color = :royalblue, linewidth = 2.5, linestyle = :dash, label = "PopPK (fit)")
    xlims!(ax, 0, 24)
    i == 1 && axislegend(ax; position = :rt, framevisible = false, labelsize = 12)
end
save(joinpath(OUT, "fig_pbpk_vs_poppk_profiles.png"), fig, px_per_unit = 2)
println("Wrote fig_pbpk_vs_poppk_profiles.png -> $OUT")

# quick numeric check: peak of mean profile each method
println("\nT50   PBPK_Cmax(mean)  PopPK_Cmax(mean)  ratio")
for t50 in t50s
    pc = maximum(pbpk[t50][2]); qc = maximum(poppk[t50][2])
    println(rpad(t50,6), rpad(round(pc,digits=2),17), rpad(round(qc,digits=2),18), round(qc/pc,digits=3))
end
