# =============================================================================
# 02_be_safe_space.R
# Steps 2-3 of the pipeline: BE assessment + safe space + VBE error rate.
#
# Reference formulation: T50 = 3.0 h (as specified). Every other formulation
# is the TEST and is assessed against it.
#
# Produces, under output/:
#   be_point_estimates.csv   one-trial paired BE (GMR, 90% CI, pass) per metric
#   vbe_trial_results.csv     bootstrap VBE pass-rate (power / Type I error)
#   safe_space_summary.csv    per-formulation BE verdict + risk classification
#   fig_safe_space_gmr.png    GMR +/- 90% CI vs T50, with 0.80-1.25 band
#   fig_vbe_passrate.png      VBE BE pass-rate vs T50 (safe-space / power curve)
#
# Run from programs/:  Rscript 02_be_safe_space.R
#
# ---- On "error rate" (VBE) ---------------------------------------------------
# A single PBPK population is one realisation. We resample n-subject crossover
# trials (bootstrap, with optional within-subject residual error wsv_cv) to get
# the PROBABILITY a real trial concludes BE. Interpreting that pass-rate:
#   * where the TRUE (population) GMR is inside 0.80-1.25 -> pass-rate = POWER
#     (correctly concluding BE; a low value is consumer-protective but failing
#      a truly-equivalent product = producer risk).
#   * where the TRUE GMR is OUTSIDE 0.80-1.25 -> the product is genuinely NOT
#     bioequivalent, so any BE pass is a false positive: pass-rate = TYPE I
#     ERROR (consumer risk). It should sit at/below ~5%.
# The "safe space" is the band of T50 over which a test formulation reliably
# passes BE against the 3.0 h reference.
# =============================================================================

source(file.path("R", "be_utils.R"))

out_dir <- "output_r"   # R cross-check outputs (primary Julia outputs -> output/)
dir.create(out_dir, showWarnings = FALSE)

pk <- read.csv(file.path(out_dir, "pk_parameters_individual.csv"),
               stringsAsFactors = FALSE)
if (!nrow(pk)) stop("Run 01_extract_pk_parameters.R first.")

ref_df  <- pk[pk$T50 == REF_T50, ]
if (!nrow(ref_df)) stop("Reference T50 = ", REF_T50, " not found in data.")
test_t50 <- sort(setdiff(unique(pk$T50), REF_T50))

# Within-subject CV scenarios for the VBE / Type I error simulation.
#   0.00 = pure PBPK pairing (between-subject variability only)
#   0.15 = illustrative within-subject/analytical CV (replace w/ metformin lit.)
WSV_SCENARIOS <- c(0.00, 0.15)
N_TRIAL <- 2000

# --- 1. Point (single-trial) paired BE for every test formulation ------------
point <- do.call(rbind, lapply(test_t50, function(tt) {
  res <- be_all_metrics(pk[pk$T50 == tt, ], ref_df)
  cbind(T50 = tt, ref_T50 = REF_T50, res)
}))
write.csv(point, file.path(out_dir, "be_point_estimates.csv"), row.names = FALSE)

# --- 2. VBE replicate trials: pass-rate (power / Type I error) ---------------
vbe <- do.call(rbind, lapply(WSV_SCENARIOS, function(cv) {
  do.call(rbind, lapply(test_t50, function(tt) {
    tdf <- pk[pk$T50 == tt, ]
    tdf <- tdf[order(tdf$IndividualId), ]
    rdf <- ref_df[order(ref_df$IndividualId), ]
    do.call(rbind, lapply(names(BE_METRICS), function(m) {
      r <- vbe_replicate(tdf[[m]], rdf[[m]], n_trial = N_TRIAL, wsv_cv = cv)
      # True (population) GMR decides whether pass-rate is power or Type I error.
      # Use only subjects with finite, positive paired values (AUC_inf can fail
      # to extrapolate for very slow-release formulations -> dropped, as in TOST).
      ok <- is.finite(tdf[[m]]) & is.finite(rdf[[m]]) & tdf[[m]] > 0 & rdf[[m]] > 0
      true_gmr <- exp(mean(log(tdf[[m]][ok])) - mean(log(rdf[[m]][ok])))
      truly_be <- (true_gmr >= BE_LOWER) & (true_gmr <= BE_UPPER)
      data.frame(T50 = tt, metric = m, wsv_cv = cv,
                 true_GMR = true_gmr, truly_BE = truly_be,
                 pass_rate = r$pass_rate,
                 risk_type = ifelse(truly_be, "power", "type_I_error"))
    }))
  }))
}))
write.csv(vbe, file.path(out_dir, "vbe_trial_results.csv"), row.names = FALSE)

# --- 3. Safe-space summary (overall BE verdict per formulation) --------------
# Overall BE requires ALL regulatory metrics (Cmax and AUC) to pass.
reg_metrics <- c("Cmax", "AUC_tEnd")   # AUC_0-t preferred for XR; see header
summary_tbl <- do.call(rbind, lapply(test_t50, function(tt) {
  p <- point[point$T50 == tt & point$metric %in% reg_metrics, ]
  data.frame(
    T50          = tt,
    Cmax_GMR     = p$GMR[p$metric == "Cmax"],
    Cmax_CI      = sprintf("[%.3f, %.3f]", p$CI_lower[p$metric=="Cmax"], p$CI_upper[p$metric=="Cmax"]),
    AUCt_GMR     = p$GMR[p$metric == "AUC_tEnd"],
    AUCt_CI      = sprintf("[%.3f, %.3f]", p$CI_lower[p$metric=="AUC_tEnd"], p$CI_upper[p$metric=="AUC_tEnd"]),
    overall_BE   = all(p$BE)
  )
}))
write.csv(summary_tbl, file.path(out_dir, "safe_space_summary.csv"), row.names = FALSE)

# --- 4. Plots ----------------------------------------------------------------
metric_cols <- c(Cmax = "#d62728", AUC_inf = "#1f77b4", AUC_tEnd = "#2ca02c")
all_t50 <- sort(unique(pk$T50))

# 4a. GMR +/- 90% CI vs T50, with the 0.80-1.25 acceptance band.
png(file.path(out_dir, "fig_safe_space_gmr.png"), width = 1000, height = 700, res = 120)
op <- par(mar = c(4.5, 4.5, 3, 1))
plot(NA, xlim = range(all_t50), ylim = c(0.6, 1.5),
     xlab = "Test formulation T50 (h)", ylab = "GMR vs reference (T50 = 3.0 h)",
     main = "Metformin XR BE safe space (paired PBPK crossover)")
rect(par("usr")[1], BE_LOWER, par("usr")[2], BE_UPPER,
     col = "#2ca02c22", border = NA)             # acceptance band
abline(h = c(BE_LOWER, 1, BE_UPPER), lty = c(2, 1, 2), col = c("grey40","grey60","grey40"))
abline(v = REF_T50, lty = 3, col = "grey50")
for (m in names(metric_cols)) {
  d <- point[point$metric == m, ]
  d <- d[order(d$T50), ]
  cl <- metric_cols[[m]]
  arrows(d$T50, d$CI_lower, d$T50, d$CI_upper, angle = 90, code = 3,
         length = 0.04, col = cl)
  lines(d$T50, d$GMR, col = cl, lwd = 2)
  points(d$T50, d$GMR, pch = 19, col = cl)
}
text(REF_T50, 1.47, "reference", pos = 4, cex = 0.8, col = "grey40")
legend("topright", legend = names(metric_cols), col = unlist(metric_cols),
       lwd = 2, pch = 19, bg = "white", cex = 0.9)
par(op); dev.off()

# 4b. VBE pass-rate vs T50 (the safe-space / power curve), wsv_cv = 0 scenario.
png(file.path(out_dir, "fig_vbe_passrate.png"), width = 1000, height = 700, res = 120)
op <- par(mar = c(4.5, 4.5, 3, 1))
plot(NA, xlim = range(all_t50), ylim = c(0, 1),
     xlab = "Test formulation T50 (h)", ylab = "P(trial concludes BE)",
     main = sprintf("VBE pass-rate vs reference T50=3.0 h (%d trials, n=24)", N_TRIAL))
abline(h = c(0.05, 0.80, 0.90), lty = 2, col = c("red","grey50","grey50"))
abline(v = REF_T50, lty = 3, col = "grey50")
v0 <- vbe[vbe$wsv_cv == 0, ]
for (m in names(metric_cols)) {
  d <- v0[v0$metric == m, ]; d <- d[order(d$T50), ]
  lines(d$T50, d$pass_rate, col = metric_cols[[m]], lwd = 2)
  points(d$T50, d$pass_rate, pch = 19, col = metric_cols[[m]])
}
text(par("usr")[1], 0.05, "5% (Type I error ceiling)", pos = 4, cex = 0.7, col = "red")
legend("right", legend = names(metric_cols), col = unlist(metric_cols),
       lwd = 2, pch = 19, bg = "white", cex = 0.9)
par(op); dev.off()

# --- 5. Console report -------------------------------------------------------
cat("\n=== BE point estimates (paired crossover vs T50 = 3.0 h) ===\n")
print(point[, c("T50","metric","n","GMR","CI_lower","CI_upper","BE")],
      row.names = FALSE, digits = 4)
cat("\n=== Safe-space summary (Cmax & AUC_0-t must both pass) ===\n")
print(summary_tbl, row.names = FALSE, digits = 4)
cat("\n=== VBE pass-rate / error rate (wsv_cv = 0) ===\n")
print(v0[, c("T50","metric","true_GMR","truly_BE","pass_rate","risk_type")],
      row.names = FALSE, digits = 4)
cat("\nWrote CSVs and figures to", normalizePath(out_dir), "\n")
