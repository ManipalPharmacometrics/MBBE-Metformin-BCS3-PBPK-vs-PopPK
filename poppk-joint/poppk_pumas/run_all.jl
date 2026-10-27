# =============================================================================
# run_all.jl  — run the full PBPK -> PopPK decision-equivalence pipeline.
#
# Each step runs in a FRESH Julia process (so top-level `const`s never clash and
# a failure stops the pipeline). Run from the poppk_pumas/ folder on a licensed
# Pumas / JuliaHub environment:
#
#   julia --project=. run_all.jl
#
# (First time only: julia --project=. -e 'using Pkg; Pkg.instantiate()')
# =============================================================================

const HERE    = @__DIR__
const SCRIPTS = ["01_build_dataset.jl",      # sparse-sample + add residual error
                 "02_fit_poppk.jl",          # JOINT Pumas fit (all arms, IVIVC)  (needs Pumas)
                 "03_poppk_metrics.jl",      # NCA endpoints from predictions
                 "04_be_safespace_overlay.jl"] # PBPK vs PopPK overlay

for s in SCRIPTS
    path = joinpath(HERE, "scripts", s)
    println("\n" * "="^78 * "\n  RUN  $s\n" * "="^78)
    cmd = `$(Base.julia_cmd()) --project=$(HERE) $(path)`
    run(cmd)
end
println("\nAll steps complete. Figures and CSVs are in output/.")
