# =============================================================================
# 06_type1_summary.jl — consolidate the Type I error results from BOTH licensed
# PopPK fits (per-arm and joint) into one reporting table, in the same layout as
# the safe-space table in METHODOLOGY_metformin_MBBE.md §4.
#
# The reported quantity is the sponsor's-eye view:
#   * PopPK declares a simulated BE boundary and reports the Type I error it
#     computes there  ->  "reported T1".
#   * The PBPK truth's acceptance rate at that same declared boundary
#                        ->  "true T1".
#   * ratio = true / reported = the Type I error mis-statement caused by letting
#     the sparse-data PopPK model place the boundary.
#
# Run AFTER 05_type1_overlay.jl has been executed in both export folders.
# =============================================================================

using CSV, DataFrames, Printf

const LEARN = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SUB   = get(ENV, "POPPK_SUBDIR", "output")   # "output" (rich) or e.g. "output_sparse"
const FITS  = ["per-arm" => joinpath(LEARN, "poppk-per-arm",
                                     "poppk_pumas", SUB),
               "joint"   => joinpath(LEARN, "poppk-joint",
                                     "poppk_pumas", SUB)]
const OUT   = joinpath(LEARN, "poppk-per-arm",
                       "poppk_pumas", SUB)

summ = DataFrame()
for (fit, dir) in FITS
    t1  = CSV.read(joinpath(dir, "type1_overlay.csv"), DataFrame)
    # PopPK's own boundary + the T1 it reports there
    pop = filter(r -> r.judged_by == "PopPK", t1)
    # PBPK truth's acceptance rate at the PopPK boundary lives on the PBPK rows
    pbp = filter(r -> r.judged_by == "PBPK", t1)
    for r in eachrow(pop)
        m = filter(x -> x.endpoint == r.endpoint && x.side == r.side &&
                        x.n_subj == r.n_subj, pbp)
        nrow(m) == 1 || continue
        truth = m[1, :]
        ismissing(truth.T1_at_other_bnd) && continue
        push!(summ, (fit = fit, endpoint = r.endpoint, side = r.side,
                     n_subj = r.n_subj,
                     limit_PBPK = truth.limit, limit_PopPK = r.limit,
                     # the two methods must hit the SAME BE limit for the paired
                     # comparison to mean anything; AUC_inf does not (flip-flop
                     # drives its PBPK GMR up through 1.25 while the PopPK GMR
                     # falls through 0.80), so those rows are flagged.
                     comparable = truth.limit == r.limit,
                     bnd_PBPK = truth.own_bnd_T50,
                     bnd_PopPK = r.own_bnd_T50,
                     dT50 = r.own_bnd_T50 - truth.own_bnd_T50,
                     pct_CQA_err = 100 * (r.own_bnd_T50 - truth.own_bnd_T50) /
                                   truth.own_bnd_T50,
                     T1_reported = r.T1_own,
                     T1_true     = truth.T1_at_other_bnd,
                     GMR_true_at_bnd = truth.GMR_at_other_bnd,
                     ratio = truth.T1_at_other_bnd / r.T1_own); cols = :union)
    end
end
sort!(summ, [:fit, :endpoint, :side, :n_subj])
CSV.write(joinpath(OUT, "type1_summary_both_fits.csv"), summ)

# ---- markdown table, same shape as the safe-space table --------------------
open(joinpath(OUT, "type1_summary_both_fits.md"), "w") do io
    println(io, "# Type I error — PBPK vs PopPK (companion to the safe-space table)\n")
    println(io, "Nominal alpha = 0.05. `bnd` = simulated BE boundary (T50, h) from ",
                "inverting the CQA -> GMR relationship. `T1 reported` = the Type I ",
                "error PopPK computes at its own boundary; `T1 true` = the PBPK ",
                "truth's acceptance rate at that same boundary; `ratio` = true / ",
                "reported.\n")
    for (fit, _) in FITS
        d = filter(r -> r.fit == fit, summ)
        println(io, "\n## PopPK fit: $fit\n")
        println(io, "| Endpoint | Side | Limit PBPK | Limit PopPK | n | bnd PBPK | ",
                    "bnd PopPK | CQA err | true GMR at PopPK bnd | T1 reported | ",
                    "T1 true | ratio | comparable |")
        println(io, "|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for r in eachrow(d)
            @printf(io, "| %s | %s | %s | %s | %d | %.2f | %.2f | %+.1f%% | %.3f | %.1f%% | %.1f%% | %.2fx | %s |\n",
                    r.endpoint, r.side, r.limit_PBPK, r.limit_PopPK, r.n_subj,
                    r.bnd_PBPK, r.bnd_PopPK, r.pct_CQA_err, r.GMR_true_at_bnd,
                    100r.T1_reported, 100r.T1_true, r.ratio,
                    r.comparable ? "yes" : "NO — different limits")
        end
    end
end

println("\n" * "="^100)
println("  TYPE I ERROR SUMMARY — both licensed PopPK fits")
println("="^100)
for (fit, _) in FITS
    d = filter(r -> r.fit == fit, summ)
    println("\n--- PopPK fit: $fit ---")
    @printf("%-24s %-5s %-11s %3s  %8s %8s %8s   %8s %8s %8s %7s %s\n",
            "endpoint", "side", "limitPB/Pop", "n", "bndPBPK", "bndPop", "CQAerr",
            "GMRtrue", "T1rep", "T1true", "ratio", "cmp")
    for r in eachrow(d)
        @printf("%-24s %-5s %4s/%-6s %3d  %8.3f %8.3f %+7.1f%%   %8.3f %7.2f%% %7.2f%% %6.2fx %s\n",
                r.endpoint, r.side, r.limit_PBPK, r.limit_PopPK, r.n_subj,
                r.bnd_PBPK, r.bnd_PopPK, r.pct_CQA_err, r.GMR_true_at_bnd,
                100r.T1_reported, 100r.T1_true, r.ratio, r.comparable ? "" : "<< DIFFERENT LIMITS")
    end
end
println("\nWrote type1_summary_both_fits.csv / .md -> $OUT")
