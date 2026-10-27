# =============================================================================
# 08_export_zoom_oc.jl — dense, high-replicate operating characteristic through
# the neighbourhood of the upper BE limit, for the "how to read these numbers"
# explainer on the results page.
#
# The point being illustrated: the acceptance curve is PINNED at 5% where the
# true GMR equals 0.80, for every study size, and n only changes how STEEPLY it
# rises away from that point. So the n=48 and n=86 curves must cross at 0.80.
#
# 20 000 trials per point, averaged over 5 seeds, so the crossing is not
# Monte-Carlo noise.
# =============================================================================

using CSV, DataFrames, Statistics, Printf

let
    ROOT = normpath(joinpath(@__DIR__, ".."))
    include(joinpath(ROOT, "src", "be_utils.jl")); BU = Main.BEUtils

    pbpk = CSV.read(joinpath(ROOT, "data", "pbpk_pk_parameters_individual.csv"), DataFrame)
    anchors = sort(unique(pbpk.T50)); ids = sort(unique(pbpk.IndividualId))
    cv = BU.WSV_CV["Cmax"]

    M = Matrix{Union{Missing,Float64}}(missing, length(ids), length(anchors))
    for (j, t) in enumerate(anchors)
        d = sort(filter(:T50 => ==(t), pbpk), :IndividualId)
        for (i, id) in enumerate(ids)
            r = findfirst(==(id), d.IndividualId); r !== nothing && (M[i, j] = d[r, "Cmax"])
        end
    end
    function ilog(v, t)
        ok = .!ismissing.(v); x = anchors[ok]; vv = Float64.(v[ok])
        t <= x[1] && return vv[1]; t >= x[end] && return vv[end]
        k = searchsortedlast(x, t); (t == x[k]) && return vv[k]
        w = (log(t) - log(x[k])) / (log(x[k+1]) - log(x[k]))
        exp((1 - w) * log(vv[k]) + w * log(vv[k+1]))
    end
    rj = findfirst(==(3.0), anchors); rv = M[:, rj]
    tvec(t) = [ilog(M[i, :], t) for i in eachindex(ids)]

    ts = collect(range(4.35, 4.95; length = 31))
    rows = DataFrame()
    for t in ts
        tv = tvec(t)
        g  = BU.be_paired(tv, rv).GMR
        a48 = mean([BU.vbe_replicate(tv, rv; n_subj = 48, n_trial = 20_000,
                                     wsv_cv = cv, seed = s).pass_rate for s in 1:5])
        a86 = mean([BU.vbe_replicate(tv, rv; n_subj = 86, n_trial = 20_000,
                                     wsv_cv = cv, seed = s).pass_rate for s in 1:5])
        push!(rows, (T50 = t, GMR = g, acc48 = a48, acc86 = a86); cols = :union)
    end
    out = joinpath(ROOT, "output", "type1_zoom_oc.csv")
    CSV.write(out, rows)

    println("Wrote $out")
    @printf("%-8s %-9s %-9s %-9s %s\n", "T50", "GMR", "n=48", "n=86", "")
    for r in eachrow(rows)
        mark = r.GMR > 0.80 ? (r.acc86 > r.acc48 ? "86 higher" : "48 higher") :
                              (r.acc86 > r.acc48 ? "86 higher" : "48 higher")
        @printf("%-8.4f %-9.5f %-9.4f %-9.4f %s\n", r.T50, r.GMR, r.acc48, r.acc86, mark)
    end
    cmp(a) = "[" * join([string(round(x, digits = 5)) for x in a], ",") * "]"
    println("\nZ_gmr=", cmp(rows.GMR))
    println("Z_a48=", cmp(rows.acc48))
    println("Z_a86=", cmp(rows.acc86))
end
