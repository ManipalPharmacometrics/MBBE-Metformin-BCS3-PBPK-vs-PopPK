# =============================================================================
# 07_export_type1_json.jl — dump the Type I error curves as JSON for the HTML
# results page, so the page's charts are drawn from the real numbers rather
# than transcribed by hand.
#
# Emits output/type1_page_data.json with:
#   oc[source]      : (gmr, acc48, acc86)  operating characteristic for Cmax
#   gmr_t50[source] : (t50, gmr)           the CQA -> GMR mapping for Cmax
#   bounds          : simulated BE boundaries, both fits
#   summary         : the cross-evaluation table, both fits
# Run after 05 (both export folders) and 06.
# =============================================================================

using CSV, DataFrames, Statistics, JSON3

let
    LEARN   = normpath(joinpath(@__DIR__, "..", "..", ".."))
    PERARM  = joinpath(LEARN, "poppk-per-arm", "poppk_pumas")
    JOINTD  = joinpath(LEARN, "poppk-joint", "poppk_pumas")

    include(joinpath(PERARM, "src", "be_utils.jl"))
    BU = Main.BEUtils

    pbpk    = CSV.read(joinpath(PERARM, "data", "pbpk_pk_parameters_individual.csv"), DataFrame)
    pop_pa  = CSV.read(joinpath(PERARM, "output", "poppk_pk_parameters_individual.csv"), DataFrame)
    pop_jt  = CSV.read(joinpath(JOINTD, "output", "poppk_pk_parameters_individual.csv"), DataFrame)

    anchors = sort(unique(pbpk.T50))
    refT50  = BU.REF_T50
    cv      = BU.WSV_CV["Cmax"]

    function mmat(pk, metric)
        ids = sort(unique(pk.IndividualId))
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
    function ilog(vals, t)
        ok = .!ismissing.(vals); sum(ok) < 2 && return missing
        x = anchors[ok]; v = Float64.(vals[ok])
        t <= x[1] && return v[1]; t >= x[end] && return v[end]
        k = searchsortedlast(x, t); (t == x[k]) && return v[k]
        w = (log(t) - log(x[k])) / (log(x[k+1]) - log(x[k]))
        exp((1 - w) * log(v[k]) + w * log(v[k+1]))
    end

    "operating characteristic + CQA->GMR mapping for one source, Cmax"
    function curves(pk)
        M  = mmat(pk, "Cmax")
        rj = findfirst(==(refT50), anchors); rv = M[:, rj]
        ts = collect(range(anchors[1], anchors[end]; length = 70))
        gmr = Float64[]; a48 = Float64[]; a86 = Float64[]
        for t in ts
            tv = [ilog(M[i, :], t) for i in 1:size(M, 1)]
            push!(gmr, BU.be_paired(tv, rv).GMR)
            push!(a48, BU.vbe_replicate(tv, rv; n_subj = 48, n_trial = 4000, wsv_cv = cv).pass_rate)
            push!(a86, BU.vbe_replicate(tv, rv; n_subj = 86, n_trial = 4000, wsv_cv = cv).pass_rate)
        end
        (t50 = ts, gmr = gmr, acc48 = a48, acc86 = a86)
    end

    payload = Dict(
        "oc" => Dict("PBPK"  => curves(pbpk),
                     "PopPK" => curves(pop_pa)),
        "bounds"  => CSV.read(joinpath(PERARM, "output", "type1_overlay_boundaries.csv"), DataFrame) |>
                     df -> [Dict(string(n) => (ismissing(r[n]) ? nothing : r[n]) for n in names(df))
                            for r in eachrow(df)],
        "bounds_joint" => CSV.read(joinpath(JOINTD, "output", "type1_overlay_boundaries.csv"), DataFrame) |>
                     df -> [Dict(string(n) => (ismissing(r[n]) ? nothing : r[n]) for n in names(df))
                            for r in eachrow(df)],
        "summary" => CSV.read(joinpath(PERARM, "output", "type1_summary_both_fits.csv"), DataFrame) |>
                     df -> [Dict(string(n) => (ismissing(r[n]) ? nothing : r[n]) for n in names(df))
                            for r in eachrow(df)],
        "control" => CSV.read(joinpath(PERARM, "output", "type1_calibration_control.csv"), DataFrame) |>
                     df -> [Dict(string(n) => r[n] for n in names(df)) for r in eachrow(df)],
        "t1_perarm" => CSV.read(joinpath(PERARM, "output", "type1_overlay.csv"), DataFrame) |>
                     df -> [Dict(string(n) => (ismissing(r[n]) ? nothing : r[n]) for n in names(df))
                            for r in eachrow(df)],
    )

    out = joinpath(PERARM, "output", "type1_page_data.json")
    open(out, "w") do io
        JSON3.pretty(io, payload)
    end
    println("Wrote $out  ($(round(filesize(out)/1024, digits=1)) KB)")
    println("OC points per source: ", length(payload["oc"]["PBPK"].gmr))
    println("PBPK  Cmax GMR range: ", round(minimum(payload["oc"]["PBPK"].gmr), digits=3),
            " .. ", round(maximum(payload["oc"]["PBPK"].gmr), digits=3))
    println("PopPK Cmax GMR range: ", round(minimum(payload["oc"]["PopPK"].gmr), digits=3),
            " .. ", round(maximum(payload["oc"]["PopPK"].gmr), digits=3))
end
