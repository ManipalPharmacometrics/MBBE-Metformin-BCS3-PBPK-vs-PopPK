# =============================================================================
# be_utils.jl
# Shared helpers for the Metformin XR PBPK virtual-bioequivalence (VBE) and
# safe-space analysis. PRIMARY implementation (Julia / Pumas ecosystem).
#
# Project: PBPK vs PopPK for BE safe space + Type I error. Reference formulation: T50 = 3.0 h (as specified).
#
# Design note -- PAIRED crossover:
#   PK-Sim simulates the SAME virtual individuals (IndividualId 0..N-1) under
#   every formulation, so test-vs-reference is a fully paired 2x2 crossover with
#   no period/sequence effect. We analyse within-subject log-differences, the
#   textbook average-BE (TOST) estimator.
#
# BE package note:
#   Pumas `Bioequivalence`/`NCA` (licensed, in PumasRegistry) can be dropped in
#   once the env is instantiated with a JuliaHub token. The open-registry
#   `Bioequivalence` UUID is broken (does not define its module), so the TOST is
#   implemented here directly. The R scripts under programs/R provide an
#   independent cross-check of these numbers.
# =============================================================================

module BEUtils

using CSV, DataFrames, Distributions, Random, Statistics

export REF_T50, BE_LOWER, BE_UPPER, BE_CI_LEVEL, BE_METRICS,
       WSV_CV, WSV_CV_IR, N_TRIAL_DEFAULT, N_SUBJ_LIST,
       load_all_formulations, be_paired, be_all_metrics, vbe_replicate,
       vbe_replicate_joint

# ---- Configuration ----------------------------------------------------------
const REF_T50      = 3.0          # reference formulation T50 (h)
const BE_LOWER     = 0.80         # FDA/EMA average-BE limits (log scale)
const BE_UPPER     = 1.25
const BE_CI_LEVEL  = 0.90         # 90% CI <=> two one-sided tests at alpha=0.05

# PK metric label => PK-Sim `Parameter` name.
#   AUC_tEnd (AUC_0-t) is preferred for XR / high-extrapolation products;
#   AUC_inf reported alongside.
const BE_METRICS = ["Cmax" => "C_max", "AUC_inf" => "AUC_inf",
                    "AUC_tEnd" => "AUC_tEnd"]

# Endpoint-specific within-subject (intra-subject) CVs for the VBE simulation.
# Metformin extended-release, sourced:
#   FDA NDA 212595 ClinPharm review: Cmax 12.2%, AUC_0-t 11.2%, AUC_inf 11.4%
#   (corroborated Cmax 14.3% / AUC 10.3%, Chung et al. 2018 TCP 26(1):6-9).
const WSV_CV = Dict("Cmax" => 0.122, "AUC_inf" => 0.114, "AUC_tEnd" => 0.112)

# High-variability sensitivity scenario — metformin IMMEDIATE-release / upper
# bound: Cmax 23.7% (IR tablets, lit. range 22.7-25.7%; observed up to 24.8%),
# AUC 14.7% (pooled 142-study mean; metformin IR study reported 13.8%).
const WSV_CV_IR = Dict("Cmax" => 0.237, "AUC_inf" => 0.147, "AUC_tEnd" => 0.147)

# VBE simulation defaults (as specified): 500 replicate trials, study sizes
# of 48 and 86 subjects. NB the PBPK population has 24 unique virtual subjects,
# so n>24 trials bootstrap-with-replacement from those 24 physiologies.
const N_TRIAL_DEFAULT = 500
const N_SUBJ_LIST     = [48, 86]

# ---- Data loading ------------------------------------------------------------

"Extract T50 (h) from a PK-Sim filename like `..._T50_3.5h_Simpop-...csv`."
function parse_t50(path::AbstractString)
    m = match(r"T50_([0-9.]+)h", basename(path))
    m === nothing && error("No T50 token in $(basename(path))")
    parse(Float64, m.captures[1])
end

"""
    read_pk_analyses(path) -> DataFrame

Read one PK-Sim `*-PK-Analyses.csv` (long: IndividualId, QuantityPath,
Parameter, Value, Unit) and return a wide per-individual table with one column
per metric in `BE_METRICS`, tagged with the formulation `T50`.
"""
function read_pk_analyses(path::AbstractString)
    raw = CSV.read(path, DataFrame; normalizenames = true)
    # normalizenames turns "IndividualId" etc. into symbols; Parameter has no spaces.
    long = filter(r -> r.Parameter in last.(BE_METRICS), raw)
    ids  = sort(unique(long.IndividualId))
    out  = DataFrame(T50 = fill(parse_t50(path), length(ids)), IndividualId = ids)
    for (label, pcol) in BE_METRICS
        sub = filter(r -> r.Parameter == pcol, long)
        d   = Dict(sub.IndividualId .=> Float64.(sub.Value))
        out[!, label] = [get(d, id, missing) for id in ids]
    end
    out
end

"Read & stack every `*-PK-Analyses.csv` under `dir`."
function load_all_formulations(dir::AbstractString)
    files = filter(f -> endswith(f, "-PK-Analyses.csv"), readdir(dir; join = true))
    isempty(files) && error("No *-PK-Analyses.csv under $dir")
    reduce(vcat, read_pk_analyses.(files))
end

# ---- Bioequivalence statistics ----------------------------------------------

"""
    be_paired(test, ref; ...) -> NamedTuple

Paired (crossover) average-BE for ONE metric: within-subject log-differences,
GMR, 90% CI (t-based, TOST), and BE pass/fail. Subjects with non-finite or
non-positive values (e.g. AUC_inf extrapolation failures) are dropped.
"""
function be_paired(test, ref; ci_level = BE_CI_LEVEL,
                   lower = BE_LOWER, upper = BE_UPPER)
    ok = @. !ismissing(test) & !ismissing(ref) &
            isfinite(test) & isfinite(ref) & (test > 0) & (ref > 0)
    d  = log.(test[ok]) .- log.(ref[ok])
    n  = length(d)
    md = mean(d)
    se = std(d) / sqrt(n)
    tq = quantile(TDist(n - 1), 1 - (1 - ci_level) / 2)
    lo = exp(md - tq * se)
    hi = exp(md + tq * se)
    (n = n, GMR = exp(md), CI_lower = lo, CI_upper = hi,
     BE = (lo >= lower) && (hi <= upper))
end

"Run `be_paired` for every metric of one test formulation vs the reference."
function be_all_metrics(test_df::DataFrame, ref_df::DataFrame)
    t = sort(test_df, :IndividualId)
    r = sort(ref_df,  :IndividualId)
    @assert t.IndividualId == r.IndividualId "subject IDs must align (paired design)"
    rows = map(first.(BE_METRICS)) do m
        res = be_paired(t[!, m], r[!, m])
        merge((metric = m,), res)
    end
    DataFrame(rows)
end

# ---- Virtual BE replicate trials (power / Type I error) ----------------------

"""
    vbe_replicate(test, ref; n_subj, n_trial, wsv_cv, seed) -> NamedTuple

Resample `n_subj`-subject crossover trials (bootstrap, with optional
within-subject residual CV `wsv_cv`) and return the BE pass fraction:
  * test != reference -> POWER (probability of correctly concluding BE)
  * test == reference -> TYPE I ERROR, but only meaningful with wsv_cv > 0
    (identical paired PBPK values give a degenerate self-comparison).
`wsv_cv` is an ASSUMPTION; set it from metformin within-subject CV literature.
"""
function vbe_replicate(test, ref; n_subj = nothing, n_trial = 2000,
                       wsv_cv = 0.0, seed = 1234,
                       ci_level = BE_CI_LEVEL, lower = BE_LOWER, upper = BE_UPPER)
    rng  = MersenneTwister(seed)
    keep = @. !ismissing(test) & !ismissing(ref) &
              isfinite(test) & isfinite(ref) & (test > 0) & (ref > 0)
    tv0  = Float64.(test[keep]); rv0 = Float64.(ref[keep])
    N    = length(rv0)
    nss  = n_subj === nothing ? N : n_subj
    sw   = sqrt(log(wsv_cv^2 + 1))          # log-scale SD from CV (0 -> 0)
    npass = 0; gmrs = Vector{Float64}(undef, n_trial)
    for t in 1:n_trial
        idx = rand(rng, 1:N, nss)
        tv  = tv0[idx]; rv = rv0[idx]
        if sw > 0
            tv = tv .* exp.(randn(rng, nss) .* sw)
            rv = rv .* exp.(randn(rng, nss) .* sw)
        end
        res = be_paired(tv, rv; ci_level, lower, upper)
        npass += res.BE
        gmrs[t] = res.GMR
    end
    (n_trial = n_trial, n_subj = nss, wsv_cv = wsv_cv,
     pass_rate = npass / n_trial, GMR_med = median(gmrs))
end

"""
    vbe_replicate_joint(testC, refC, testA, refA; cv_C, cv_A, ...) -> NamedTuple

Joint VBE: a trial concludes BE only if BOTH Cmax and AUC pass. Returns the
joint pass-rate. `cv_C`, `cv_A` are the metric-specific within-subject CVs.

Within-subject CV defaults (metformin ER, intra-subject):
  Cmax = 0.122, AUC0-t = 0.112  — FDA NDA 212595 ClinPharm review
  (corroborated: Cmax 0.143 / AUC 0.103, Chung et al. 2018 TCP 26(1):6-9).
"""
function vbe_replicate_joint(testC, refC, testA, refA;
                             n_subj = nothing, n_trial = 2000,
                             cv_C = 0.122, cv_A = 0.112, seed = 1234,
                             ci_level = BE_CI_LEVEL, lower = BE_LOWER, upper = BE_UPPER)
    rng  = MersenneTwister(seed)
    keep = @. !ismissing(testC) & !ismissing(refC) & !ismissing(testA) & !ismissing(refA) &
              isfinite(testC) & isfinite(refC) & isfinite(testA) & isfinite(refA) &
              (testC > 0) & (refC > 0) & (testA > 0) & (refA > 0)
    tc = Float64.(testC[keep]); rc = Float64.(refC[keep])
    ta = Float64.(testA[keep]); ra = Float64.(refA[keep])
    N   = length(tc)
    nss = n_subj === nothing ? N : n_subj
    swC = sqrt(log(cv_C^2 + 1)); swA = sqrt(log(cv_A^2 + 1))
    npass = 0
    for _ in 1:n_trial
        idx = rand(rng, 1:N, nss)
        tcv = tc[idx]; rcv = rc[idx]; tav = ta[idx]; rav = ra[idx]
        if swC > 0
            tcv = tcv .* exp.(randn(rng, nss) .* swC); rcv = rcv .* exp.(randn(rng, nss) .* swC)
        end
        if swA > 0
            tav = tav .* exp.(randn(rng, nss) .* swA); rav = rav .* exp.(randn(rng, nss) .* swA)
        end
        rC = be_paired(tcv, rcv; ci_level, lower, upper)
        rA = be_paired(tav, rav; ci_level, lower, upper)
        npass += (rC.BE && rA.BE)
    end
    (n_trial = n_trial, n_subj = nss, cv_C = cv_C, cv_A = cv_A,
     pass_rate = npass / n_trial)
end

end # module
