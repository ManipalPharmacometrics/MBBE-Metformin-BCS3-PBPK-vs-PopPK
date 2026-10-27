# =============================================================================
# 01_extract_pk_parameters.R
# Step 1 of the BE/safe-space pipeline.
#
# Reads the PK-Sim population exports (one *-PK-Analyses.csv per formulation in
# pksim-output/popsims/) and assembles a single tidy per-individual
# dataset of the BE metrics (Cmax, AUC_inf, AUC_tEnd) across the T50 sweep.
#
# Output: output/pk_parameters_individual.csv
# Run from the programs/ directory:  Rscript 01_extract_pk_parameters.R
# =============================================================================

this_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
source(file.path("R", "be_utils.R"))

# R is the INDEPENDENT CROSS-CHECK; writes to output_r/ so it never collides
# with the primary Julia outputs in output/.
popsims_dir <- file.path("..", "pksim-output", "popsims")
out_dir     <- "output_r"
dir.create(out_dir, showWarnings = FALSE)

cat("Reading population PK analyses from:", normalizePath(popsims_dir), "\n")
pk <- load_all_formulations(popsims_dir)
pk <- pk[order(pk$T50, pk$IndividualId), ]

out_csv <- file.path(out_dir, "pk_parameters_individual.csv")
write.csv(pk, out_csv, row.names = FALSE)

cat("\nFormulations (T50, h):", paste(sort(unique(pk$T50)), collapse = ", "), "\n")
cat("Subjects per formulation:", length(unique(pk$IndividualId)), "\n")
cat("Reference formulation (T50):", REF_T50, "h\n")
cat("Wrote", nrow(pk), "rows ->", out_csv, "\n\n")

# Quick geometric-mean summary per formulation (sanity check vs comp.csv).
gm <- function(x) exp(mean(log(x[is.finite(x) & x > 0])))
summ <- do.call(rbind, lapply(split(pk, pk$T50), function(d) {
  data.frame(T50 = d$T50[1], n = nrow(d),
             gm_Cmax = gm(d$Cmax),
             gm_AUC_inf = gm(d$AUC_inf),
             gm_AUC_tEnd = gm(d$AUC_tEnd))
}))
cat("Geometric-mean PK metrics by formulation:\n")
print(summ, row.names = FALSE, digits = 5)
