# =============================================================================
# 05_type1_overlay.jl  — PBPK vs PopPK Type I error, the companion to 04.
#
# Same design, same machinery (BEUtils), same two sources as the safe-space
# overlay — but instead of asking "where does power fall below 80%?" it asks
# Morris et al.'s question: "where is the SIMULATED BE BOUNDARY, and what is the
# Type I error there?"
#
# Three steps, per source (PBPK, PopPK) and per endpoint:
#   1. INVERT the CQA -> GMR relationship: find the T50 nearest the reference at
#      which the deterministic GMR first leaves [0.80, 1.25]. That T50 is the
#      "simulated BE boundary" (Morris Part 2). Done on the same fine T50 grid
#      and the same per-subject log-interpolation used for the power curves.
#   2. TYPE I ERROR = the Monte-Carlo BE acceptance rate AT that boundary T50
#      (true GMR sits exactly on 0.80 / 1.25). Nominal target = 5%.
#   3. CROSS-EVALUATION (the headline): evaluate the PBPK truth's acceptance
#      rate at the boundary POPPK would have declared, and vice versa. Morris's
#      warning is that the acceptance curve is very steep near the limit, so a
#      small error in the boundary produces a large error in the Type I error
#      estimate. This quantifies that for the PBPK-vs-PopPK choice.
#
# Also emits a pure TOST calibration control (reference scaled by exactly 0.80),
# which isolates the test's own calibration from any CQA-inversion error.
#
# Inputs : data/pbpk_pk_parameters_individual.csv     (PBPK, from the PK-Sim pipeline)
#          output/poppk_pk_parameters_individual.csv  (PopPK, from 03)
# Outputs: output/type1_overlay_boundaries.csv    simulated BE boundaries
#          output/type1_overlay.csv               Type I error, own + cross
#          output/type1_calibration_control.csv   scaled-reference TOST control
#          output/fig_overlay_type1_oc_<ep>.png   acceptance rate vs true GMR
#          output/fig_overlay_type1_bars.png      Type I error, own vs cross
# =============================================================================

using CSV, DataFrames, Statistics, Distributions, CairoMakie

const HERE = @__DIR__
const ROOT = normpath(joinpath(HERE, ".."))
const OUT  = get(ENV, "POPPK_OUT", joinpath(ROOT, "output"));  mkpath(OUT)

include(joinpath(ROOT, "src", "be_utils.jl"));  using .BEUtils

pbpk  = CSV.read(joinpath(ROOT, "data", "pbpk_pk_parameters_individual.csv"), DataFrame)
poppk = CSV.read(joinpath(OUT, "poppk_pk_parameters_individual.csv"), DataFrame)
SOURCES = ["PBPK" => pbpk, "PopPK" => poppk]
SRCMAP  = Dict(SOURCES)
COL = Dict("PBPK" => :firebrick, "PopPK" => :royalblue)

# Type I error is a ~5% tail probability, so it needs far more replicates than a
# power curve: 500 trials gives SE ~1.0 percentage points, 20 000 gives ~0.15.
const N_TRIAL_T1 = 20_000
const N_TRIAL_OC = 2_000        # for the operating-characteristic curves
const GRID_DX    = 0.005        # fine grid for the boundary inversion
const NOMINAL_A  = 0.05

# ---- per-subject interpolation across formulations (identical to 04) --------
function metric_matrix(pk, metric, anchors, ids)
    M = Matrix{Union{Missing,Float64}}(missing, length(ids), length(anchors))
    for (j, t) in enumerate(anchors)
        d = sort(filter(:T50 => ==(t), pk), :IndividualId)
        for (i, id) in enumerate(ids)
            r = findfirst(==(id), d.IndividualId)
            r !== nothing && (M[i, j] = d[r, metric])
        end
    end
    M
end
function interp_log(vals, anchors, t)
    ok = .!ismissing.(vals); sum(ok) < 2 && return missing
    x = anchors[ok]; v = Float64.(vals[ok])
    t <= x[1]   && return v[1]
    t >= x[end] && return v[end]
    k = searchsortedlast(x, t); (t == x[k]) && return v[k]
    w = (log(t) - log(x[k])) / (log(x[k+1]) - log(x[k]))
    exp((1 - w) * log(v[k]) + w * log(v[k+1]))
end

"Per-subject endpoint vector at an arbitrary T50, plus the reference vector."
function test_ref_at(pk, metric, t, anchors, ids)
    M  = metric_matrix(pk, metric, anchors, ids)
    rj = findfirst(==(REF_T50), anchors)
    ([interp_log(M[i, :], anchors, t) for i in eachindex(ids)], M[:, rj])
end

# ---- 1. invert CQA -> GMR to locate the simulated BE boundary ---------------
"Deterministic paired GMR(T50) over `grid`."
function gmr_curve(pk, metric, grid, anchors, ids)
    M  = metric_matrix(pk, metric, anchors, ids)
    rj = findfirst(==(REF_T50), anchors); rv = M[:, rj]
    [BEUtils.be_paired([interp_log(M[i, :], anchors, t) for i in eachindex(ids)], rv).GMR
     for t in grid]
end

"""
    first_exit(grid, gmr, side) -> (t, limit, target)

Walk outward from the reference and return the first T50 at which GMR crosses
either BE limit. `side = :low` walks down in T50, `:high` walks up. Linear
interpolation between the bracketing grid points. Returns `t = missing` if the
GMR never leaves [0.80, 1.25] within the simulated CQA range.
"""
function first_exit(grid, gmr, side)
    ri   = argmin(abs.(grid .- REF_T50))
    step = side === :low ? -1 : 1
    i    = ri
    while (side === :low ? i > 1 : i < length(grid))
        j = i + step
        for (nm, tgt) in (("0.80", BE_LOWER), ("1.25", BE_UPPER))
            a, b = gmr[i], gmr[j]
            if (a - tgt) * (b - tgt) <= 0 && a != b
                w = (tgt - a) / (b - a)
                return (t = grid[i] + w * (grid[j] - grid[i]), limit = nm, target = tgt)
            end
        end
        i = j
    end
    (t = missing, limit = "not reached", target = NaN)
end

# ---- 2. Monte-Carlo BE acceptance rate at a given T50 ----------------------
function accept_single(pk, metric, cv, t, anchors, ids, nss; n_trial = N_TRIAL_T1)
    tv, rv = test_ref_at(pk, metric, t, anchors, ids)
    BEUtils.vbe_replicate(tv, rv; n_subj = nss, n_trial = n_trial, wsv_cv = cv).pass_rate
end
function accept_joint(pk, cvC, cvA, t, anchors, ids, nss; n_trial = N_TRIAL_T1)
    tC, rC = test_ref_at(pk, "Cmax", t, anchors, ids)
    tA, rA = test_ref_at(pk, "AUC_tEnd", t, anchors, ids)
    BEUtils.vbe_replicate_joint(tC, rC, tA, rA; n_subj = nss, n_trial = n_trial,
                                cv_C = cvC, cv_A = cvA).pass_rate
end

"Wilson 95% CI for a Monte-Carlo proportion (so the reader can see MC noise)."
function wilson(p, n; z = 1.96)
    d = 1 + z^2 / n
    c = (p + z^2 / (2n)) / d
    h = z * sqrt(p * (1 - p) / n + z^2 / (4n^2)) / d
    (max(0.0, c - h), min(1.0, c + h))
end

# ---- endpoint list — identical to 04 ---------------------------------------
endpoints = [("Cmax",                   :single, WSV_CV["Cmax"]),
             ("AUC_tEnd",               :single, WSV_CV["AUC_tEnd"]),
             ("AUC_inf",                :single, WSV_CV["AUC_inf"]),
             ("Joint (Cmax & AUC_0-t)", :joint,  (WSV_CV["Cmax"], WSV_CV["AUC_tEnd"]))]
metric_of(label) = label == "AUC_tEnd" ? "AUC_tEnd" :
                   label == "AUC_inf"  ? "AUC_inf"  : "Cmax"

const ANCHORS = sort(unique(pbpk.T50))
const GRID    = collect(range(ANCHORS[1], ANCHORS[end]; step = GRID_DX))
ids_of(pk)    = sort(unique(pk.IndividualId))

# =============================================================================
#  Boundary inversion
# =============================================================================
bnd = DataFrame()
for (label, kind, _) in endpoints, (src, pk) in SOURCES
    ids = ids_of(pk)
    if kind === :joint
        # the joint boundary is whichever endpoint leaves the limits FIRST
        gC = gmr_curve(pk, "Cmax",     GRID, ANCHORS, ids)
        gA = gmr_curve(pk, "AUC_tEnd", GRID, ANCHORS, ids)
        for side in (:low, :high)
            cands = [(first_exit(GRID, gC, side), "Cmax"),
                     (first_exit(GRID, gA, side), "AUC_0-t")]
            cands = filter(c -> !ismissing(c[1].t), cands)
            if isempty(cands)
                push!(bnd, (endpoint = label, source = src, side = String(side),
                            bnd_T50 = missing, limit = "not reached",
                            binding = "none"); cols = :union)
            else
                (e, who) = argmin(c -> abs(c[1].t - REF_T50), cands)
                push!(bnd, (endpoint = label, source = src, side = String(side),
                            bnd_T50 = e.t, limit = e.limit, binding = who); cols = :union)
            end
        end
    else
        g = gmr_curve(pk, metric_of(label), GRID, ANCHORS, ids)
        for side in (:low, :high)
            e = first_exit(GRID, g, side)
            push!(bnd, (endpoint = label, source = src, side = String(side),
                        bnd_T50 = e.t, limit = e.limit,
                        binding = metric_of(label)); cols = :union)
        end
    end
end
CSV.write(joinpath(OUT, "type1_overlay_boundaries.csv"), bnd)

println("\n" * "="^78)
println("  SIMULATED BE BOUNDARIES (CQA inversion: T50 where GMR hits 0.80 / 1.25)")
println("="^78)
show(bnd; allrows = true, allcols = true); println()

# =============================================================================
#  Type I error at own boundary + cross-evaluation
# =============================================================================
t1 = DataFrame()
for (label, kind, cv) in endpoints, side in ("low", "high"), nss in N_SUBJ_LIST
    rows = filter(r -> r.endpoint == label && r.side == side, bnd)
    for (src, pk) in SOURCES
        me = only(filter(r -> r.source == src, rows))
        ismissing(me.bnd_T50) && continue
        other = src == "PBPK" ? "PopPK" : "PBPK"
        oth   = only(filter(r -> r.source == other, rows))

        # Type I error at THIS source's own boundary, judged by THIS source
        own = kind === :joint ?
              accept_joint(pk, cv[1], cv[2], me.bnd_T50, ANCHORS, ids_of(pk), nss) :
              accept_single(pk, metric_of(label), cv, me.bnd_T50, ANCHORS, ids_of(pk), nss)

        # ... and what THIS source says at the OTHER method's boundary.
        # With src = PBPK this is the number that matters: the truth's acceptance
        # rate at the limit PopPK would have set.
        cross, cross_gmr = if ismissing(oth.bnd_T50)
            (missing, missing)
        else
            cr = kind === :joint ?
                 accept_joint(pk, cv[1], cv[2], oth.bnd_T50, ANCHORS, ids_of(pk), nss) :
                 accept_single(pk, metric_of(label), cv, oth.bnd_T50, ANCHORS, ids_of(pk), nss)
            gm = BEUtils.be_paired(test_ref_at(pk, me.binding, oth.bnd_T50, ANCHORS, ids_of(pk))...).GMR
            (cr, gm)
        end

        lo, hi = wilson(own, N_TRIAL_T1)
        push!(t1, (endpoint = label, side = side, n_subj = nss, judged_by = src,
                   own_bnd_T50 = me.bnd_T50, limit = me.limit, binding = me.binding,
                   T1_own = own, T1_own_lo = lo, T1_own_hi = hi,
                   other_bnd_T50 = oth.bnd_T50,
                   dT50 = ismissing(oth.bnd_T50) ? missing : oth.bnd_T50 - me.bnd_T50,
                   T1_at_other_bnd = cross, GMR_at_other_bnd = cross_gmr,
                   n_trial = N_TRIAL_T1); cols = :union)
    end
end
CSV.write(joinpath(OUT, "type1_overlay.csv"), t1)

println("\n" * "="^78)
println("  TYPE I ERROR (acceptance rate at the simulated BE boundary, $(N_TRIAL_T1) trials)")
println("  nominal = $(NOMINAL_A)")
println("="^78)
show(select(t1, :endpoint, :side, :n_subj, :judged_by, :own_bnd_T50, :limit,
            :T1_own, :other_bnd_T50, :dT50, :T1_at_other_bnd, :GMR_at_other_bnd);
     allrows = true, allcols = true); println()

# =============================================================================
#  Pure TOST calibration control — reference scaled to exactly GMR = 0.80.
#  Isolates the test's own calibration from any CQA-inversion error.
# =============================================================================
ctl = DataFrame()
for (label, kind, cv) in endpoints, (src, pk) in SOURCES, nss in N_SUBJ_LIST
    ids = ids_of(pk)
    if kind === :joint
        _, rC = test_ref_at(pk, "Cmax", REF_T50, ANCHORS, ids)
        _, rA = test_ref_at(pk, "AUC_tEnd", REF_T50, ANCHORS, ids)
        p = BEUtils.vbe_replicate_joint(rC .* BE_LOWER, rC, rA .* BE_LOWER, rA;
                                        n_subj = nss, n_trial = N_TRIAL_T1,
                                        cv_C = cv[1], cv_A = cv[2]).pass_rate
    else
        _, rv = test_ref_at(pk, metric_of(label), REF_T50, ANCHORS, ids)
        p = BEUtils.vbe_replicate(rv .* BE_LOWER, rv; n_subj = nss,
                                  n_trial = N_TRIAL_T1, wsv_cv = cv).pass_rate
    end
    lo, hi = wilson(p, N_TRIAL_T1)
    push!(ctl, (endpoint = label, source = src, n_subj = nss,
                T1_scaled_ref = p, ci_lo = lo, ci_hi = hi); cols = :union)
end
CSV.write(joinpath(OUT, "type1_calibration_control.csv"), ctl)
println("\nTOST calibration control (reference x 0.80 -> true GMR exactly 0.80):")
show(ctl; allrows = true, allcols = true); println()

# =============================================================================
#  Figures
# =============================================================================
# (a) operating-characteristic curves: acceptance rate vs TRUE GMR
for (label, kind, cv) in endpoints
    token = replace(label, r"[^A-Za-z0-9]" => "")
    fig = Figure(size = (1200, 520), fontsize = 16)
    Label(fig[0, 1:2],
          "Type I error / operating characteristic — $label — PBPK vs PopPK";
          fontsize = 18, font = :bold)
    for (pcol, nss) in enumerate(N_SUBJ_LIST)
        ax = Axis(fig[1, pcol]; xlabel = "true GMR (test/ref)",
                  ylabel = pcol == 1 ? "P(conclude BE)" : "", title = "n = $nss")
        ylims!(ax, -0.02, 1.02)
        vlines!(ax, [BE_LOWER, BE_UPPER]; linestyle = :dash, color = :gray40)
        hlines!(ax, [NOMINAL_A]; linestyle = :dot, color = :gray50)
        for (src, pk) in SOURCES
            ids  = ids_of(pk)
            mkey = kind === :joint ? "Cmax" : metric_of(label)
            g    = gmr_curve(pk, mkey, GRID, ANCHORS, ids)
            # thin the grid for the OC curve, keep only 0.70 <= GMR <= 1.35
            sel  = [i for i in 1:8:length(GRID) if 0.70 <= g[i] <= 1.35]
            xs   = g[sel]
            ys   = [kind === :joint ?
                    accept_joint(pk, cv[1], cv[2], GRID[i], ANCHORS, ids, nss;
                                 n_trial = N_TRIAL_OC) :
                    accept_single(pk, metric_of(label), cv, GRID[i], ANCHORS, ids, nss;
                                  n_trial = N_TRIAL_OC) for i in sel]
            o = sortperm(xs)
            lines!(ax, xs[o], ys[o]; color = COL[src], linewidth = 2.5, label = src)
            for side in ("low", "high")
                r = filter(x -> x.endpoint == label && x.side == side && x.source == src, bnd)
                nrow(r) == 1 && !ismissing(r[1, :bnd_T50]) &&
                    vlines!(ax, [gmr_curve(pk, mkey, [r[1, :bnd_T50]], ANCHORS, ids)[1]];
                            color = COL[src], linestyle = :dot, linewidth = 1.5)
            end
        end
        pcol == 2 && axislegend(ax; position = :rb, framevisible = false)
    end
    save(joinpath(OUT, "fig_overlay_type1_oc_$(token).png"), fig, px_per_unit = 2)
end

# (b) Type I error bars: own boundary vs the other method's boundary
let
    d = filter(r -> r.judged_by == "PBPK" && !ismissing(r.T1_at_other_bnd), t1)
    if nrow(d) > 0
        labs = ["$(r.endpoint)\n$(r.side), n=$(r.n_subj)" for r in eachrow(d)]
        fig = Figure(size = (max(1100, 90 * nrow(d)), 480), fontsize = 15)
        ax = Axis(fig[1, 1];
                  title = "Type I error judged by PBPK truth: at its own boundary vs at the PopPK boundary",
                  ylabel = "P(conclude BE)", xticks = (1:nrow(d), labs),
                  xticklabelrotation = pi / 4)
        barplot!(ax, (1:nrow(d)) .- 0.19, d.T1_own; width = 0.36,
                 color = :firebrick, label = "at PBPK boundary (nominal)")
        barplot!(ax, (1:nrow(d)) .+ 0.19, Float64.(d.T1_at_other_bnd); width = 0.36,
                 color = :royalblue, label = "at PopPK boundary")
        hlines!(ax, [NOMINAL_A]; linestyle = :dash, color = :gray30)
        axislegend(ax; position = :rt, framevisible = false)
        save(joinpath(OUT, "fig_overlay_type1_bars.png"), fig, px_per_unit = 2)
    end
end

println("\nWrote: type1_overlay_boundaries.csv, type1_overlay.csv,")
println("       type1_calibration_control.csv,")
println("       fig_overlay_type1_oc_*.png, fig_overlay_type1_bars.png -> $OUT")
