# =============================================================================
# STEP 1 — build the PopPK dataset from the PK-Sim profiles.
#
# In : data/pksim_profiles/Metformin_T50_<x>h_Simpop-Results.csv
#      (8 formulations x 24 virtual subjects, plasma conc. every 3 min, 0–24 h)
# Out: results/step1_dataset.csv   (NONMEM-style: one dose row + 13 samples per
#                                    subject per formulation)
#
# What it does, per formulation and subject:
#   1. take the PK-Sim concentration at each time in SCHEDULE
#   2. add 15% proportional assay error:  obs = true * (1 + 0.15 * N(0,1))
#   3. write a dose row (390 mg base = 3019.5 µmol at t = 0) and the samples
# =============================================================================

using CSV, DataFrames, Random
include("settings.jl")
mkpath(RESULTS)

const CONC = "Organism|PeripheralVenousBlood|Metformin|Plasma (Peripheral Venous Blood) [µmol/l]"
t50_of(file) = parse(Float64, match(r"T50_([0-9.]+)h", file).captures[1])

files = sort(readdir(joinpath(DATA, "pksim_profiles")), by = t50_of)
rng   = MersenneTwister(SEED_DATASET)
rows  = DataFrame()

for file in files
    t50 = t50_of(file)
    sim = CSV.read(joinpath(DATA, "pksim_profiles", file), DataFrame)
    for id in sort(unique(sim.IndividualId))
        subj = sim[sim.IndividualId .== id, :]
        sid  = "T50_$(t50)_ID_$(id)"                       # one ID per subject per formulation
        push!(rows, (SUBJID = sid, T50 = t50, IndividualId = id, time = 0.0,
                     evid = 1, amt = DOSE_UMOL, cmt = 1, dv = missing); promote = true)
        for t in SCHEDULE
            k = findfirst(==(round(Int, t * 60)), subj[!, "Time [min]"])   # PK-Sim time is in minutes
            ctrue = subj[k, CONC]
            cobs  = max(ctrue * (1 + PROP_CV * randn(rng)), 1e-6)       # proportional error
            push!(rows, (SUBJID = sid, T50 = t50, IndividualId = id, time = t,
                         evid = 0, amt = 0.0, cmt = 2, dv = cobs); promote = true)
        end
    end
end

CSV.write(joinpath(RESULTS, "step1_dataset.csv"), rows)
println("Formulations (T50, h): ", join(t50_of.(files), ", "))
println("Subjects per formulation: ", length(unique(rows.IndividualId)))
println("Samples per subject: ", length(SCHEDULE), "  (", join(SCHEDULE, ", "), " h)")
println("Wrote results/step1_dataset.csv  (", nrow(rows), " rows)")
