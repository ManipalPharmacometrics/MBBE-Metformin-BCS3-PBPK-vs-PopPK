# =============================================================================
# 01_build_dataset.jl  — assemble the simulated BE-study dataset.
#
# Takes the 8 dense PK-Sim population profiles (data/popsims/*-Results.csv,
# 24 paired subjects each, 3-min grid 0-24 h) as the underlying "truth", then:
#   1. sparse-samples each subject at a realistic BE schedule,
#   2. adds 15% proportional residual (assay/intra-occasion) error,
#   3. writes a Pumas-ready NONMEM-style dataset (one dose row + obs rows per
#      subject per formulation), with T50 and the original IndividualId kept.
#
# Output: output/sparse_be_dataset.csv
#
# Design choices (2026-06-10):
#   schedule = rich (14 nominal times incl. predose);
#   residual = proportional 15% CV;
#   reference formulation = T50 3.0 h.
# Deferred: a truly-sparse (~6 pt) schedule — see SPARSE_SCHEDULE below.
# =============================================================================

using CSV, DataFrames, Random, Statistics

const HERE     = @__DIR__
const ROOT     = normpath(joinpath(HERE, ".."))
const POPSIMS  = joinpath(ROOT, "data", "popsims")
const OUT      = joinpath(ROOT, "output");  mkpath(OUT)

# Dose: clinical 500 mg metformin HCl tablet = 390 mg free base (factor 0.78).
# PK-Sim doses in free base; convert mg -> µmol (MW metformin base = 129.16).
const DOSE_MG_BASE = 390.0
const MW_BASE      = 129.16
const DOSE_UMOL    = DOSE_MG_BASE / MW_BASE * 1000        # ≈ 3019.5 µmol

# Rich BE sampling schedule (hours). Predose (0 h) is the dose event only — its
# zero concentration is not an observation. Every nominal time below lands
# exactly on the PK-Sim 3-min grid (verified).
const RICH_SCHEDULE   = [0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 16.0, 24.0]
# Deferred truly-sparse alternative (set SCHEDULE = SPARSE_SCHEDULE to try):
const SPARSE_SCHEDULE = [1.0, 3.0, 6.0, 12.0, 24.0]
const SCHEDULE        = RICH_SCHEDULE

const PROP_CV = 0.15          # proportional residual error (CV)
const SEED    = 20260610      # reproducibility

# Column in the Results.csv holding peripheral venous plasma (µmol/L).
const CONC_COL = "Organism|PeripheralVenousBlood|Metformin|Plasma (Peripheral Venous Blood) [µmol/l]"
const TIME_COL = "Time [min]"
const ID_COL   = "IndividualId"

parse_t50(path) = parse(Float64, match(r"T50_([0-9.]+)h", basename(path)).captures[1])

files = sort(filter(f -> endswith(f, "-Results.csv"),
                    readdir(POPSIMS; join = true)), by = parse_t50)
isempty(files) && error("No *-Results.csv under $POPSIMS")

rng  = MersenneTwister(SEED)
rows = DataFrame(SUBJID = String[], T50 = Float64[], IndividualId = Int[],
                 time = Float64[], evid = Int[], amt = Union{Missing,Float64}[],
                 cmt = Int[], dv = Union{Missing,Float64}[])

for f in files
    t50 = parse_t50(f)
    df  = CSV.read(f, DataFrame)                       # keep original header text
    tcol = df[!, TIME_COL]; ccol = df[!, CONC_COL]; icol = df[!, ID_COL]
    sched_min = round.(Int, SCHEDULE .* 60)
    for id in sort(unique(icol))
        sub  = view(df, icol .== id, :)
        smin = sub[!, TIME_COL]; sc = sub[!, CONC_COL]
        sid  = "T50_$(t50)_ID_$(id)"
        # dose event at t = 0 into the depot (cmt 1)
        push!(rows, (sid, t50, id, 0.0, 1, DOSE_UMOL, 1, missing))
        for (th, tm) in zip(SCHEDULE, sched_min)
            k = findfirst(==(tm), smin)
            k === nothing && error("time $tm min absent for $sid in $(basename(f))")
            ctrue = Float64(sc[k])
            cobs  = ctrue * (1 + PROP_CV * randn(rng))   # proportional error
            cobs  = max(cobs, 1e-6)                      # guard the rare negative draw
            push!(rows, (sid, t50, id, th, 0, 0.0, 2, cobs))
        end
    end
end

sort!(rows, [:T50, :IndividualId, :time])
out_csv = joinpath(OUT, "sparse_be_dataset.csv")
CSV.write(out_csv, rows)

println("Dose: $(round(DOSE_UMOL, digits=1)) µmol  ($(DOSE_MG_BASE) mg free base)")
println("Schedule ($(length(SCHEDULE)) obs/subject, h): ", join(SCHEDULE, ", "))
println("Residual error: proportional $(round(Int, 100PROP_CV))% CV   seed=$SEED")
println("Formulations (T50 h): ", join(sort(unique(rows.T50)), ", "))
println("Subjects/formulation: ", length(unique(rows.IndividualId)))
println("Wrote $(nrow(rows)) rows -> $out_csv")
