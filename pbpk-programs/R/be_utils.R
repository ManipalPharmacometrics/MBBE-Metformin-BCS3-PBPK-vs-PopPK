# =============================================================================
# be_utils.R
# Shared helpers for the Metformin XR PBPK virtual-bioequivalence (VBE) and
# safe-space analysis.
#
# Project: PBPK vs PopPK for BE safe space / Type I error.
# Data:    PK-Sim population simulations exported per formulation (T50 sweep).
#          One value per virtual subject per PK metric per formulation.
# Design:  The SAME virtual individuals (IndividualId 0..N-1) are simulated
#          under every formulation, so test vs reference is a fully PAIRED
#          crossover with no period/sequence effect. We therefore analyse the
#          within-subject log-differences (the textbook 2x2 crossover estimator).
#
# Base R only (no external package dependencies).
# =============================================================================

# ---- Configuration ----------------------------------------------------------

# Reference formulation for ALL comparisons (as specified): T50 = 3.0 h.
REF_T50 <- 3.0

# BE acceptance limits (FDA/EMA average BE on the log scale).
BE_LOWER <- 0.80
BE_UPPER <- 1.25

# 90% CI => two one-sided tests at alpha = 0.05.
BE_CI_LEVEL <- 0.90

# PK metrics carried through the BE analysis. PK-Sim parameter names on the LHS.
#   C_max     : peak concentration
#   AUC_inf   : AUC extrapolated to infinity
#   AUC_tEnd  : AUC to last simulated time (AUC_0-t) -- preferred for XR / high
#               extrapolation products, reported alongside AUC_inf.
BE_METRICS <- c(Cmax = "C_max", AUC_inf = "AUC_inf", AUC_tEnd = "AUC_tEnd")


# ---- Data loading ------------------------------------------------------------

# Pull the T50 (h) value out of a PK-Sim export filename like
# "Metformin_T50_3.5h_Simpop-PK-Analyses.csv".
parse_t50_from_name <- function(path) {
  m <- regmatches(basename(path), regexpr("T50_[0-9.]+h", basename(path)))
  as.numeric(sub("h$", "", sub("^T50_", "", m)))
}

# Read one PK-Sim "*-PK-Analyses.csv" (long format: IndividualId, QuantityPath,
# Parameter, Value, Unit) and return a wide per-individual data.frame with one
# column per requested metric, tagged with the formulation's T50.
read_pk_analyses <- function(path, metrics = BE_METRICS) {
  raw <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                  fileEncoding = "UTF-8-BOM")
  names(raw) <- trimws(names(raw))
  keep <- raw[raw$Parameter %in% metrics, c("IndividualId", "Parameter", "Value")]
  keep$Value <- as.numeric(keep$Value)

  ids <- sort(unique(keep$IndividualId))
  out <- data.frame(T50 = parse_t50_from_name(path), IndividualId = ids)
  for (nm in names(metrics)) {
    pcol <- metrics[[nm]]
    sub <- keep[keep$Parameter == pcol, ]
    out[[nm]] <- sub$Value[match(ids, sub$IndividualId)]
  }
  out
}

# Read every "*-PK-Analyses.csv" under a directory and stack them.
load_all_formulations <- function(dir, metrics = BE_METRICS) {
  files <- list.files(dir, pattern = "-PK-Analyses\\.csv$", full.names = TRUE)
  if (!length(files)) stop("No *-PK-Analyses.csv files found in ", dir)
  parts <- lapply(files, read_pk_analyses, metrics = metrics)
  do.call(rbind, parts)
}


# ---- Bioequivalence statistics ----------------------------------------------

# Paired (crossover) average-BE for ONE metric, ONE test vs the reference.
#   test_vals / ref_vals : numeric vectors aligned subject-by-subject.
# Returns GMR, 90% CI, and BE pass/fail on the log scale via TOST.
be_paired <- function(test_vals, ref_vals, ci_level = BE_CI_LEVEL,
                      lower = BE_LOWER, upper = BE_UPPER) {
  ok <- is.finite(test_vals) & is.finite(ref_vals) &
        test_vals > 0 & ref_vals > 0
  d  <- log(test_vals[ok]) - log(ref_vals[ok])   # within-subject log-difference
  n  <- length(d)
  md <- mean(d)
  se <- sd(d) / sqrt(n)
  # 90% CI -> alpha 0.05 each side -> t quantile at (1 - 0.05).
  tq <- qt(1 - (1 - ci_level) / 2, df = n - 1)
  ci_lo <- exp(md - tq * se)
  ci_hi <- exp(md + tq * se)
  data.frame(
    n        = n,
    GMR      = exp(md),
    CI_lower = ci_lo,
    CI_upper = ci_hi,
    BE       = (ci_lo >= lower) & (ci_hi <= upper)
  )
}

# Convenience: run be_paired for every metric for one test formulation.
be_all_metrics <- function(test_df, ref_df, metrics = names(BE_METRICS)) {
  ord_t <- test_df[order(test_df$IndividualId), ]
  ord_r <- ref_df[order(ref_df$IndividualId), ]
  stopifnot(identical(ord_t$IndividualId, ord_r$IndividualId))
  do.call(rbind, lapply(metrics, function(m) {
    res <- be_paired(ord_t[[m]], ord_r[[m]])
    cbind(metric = m, res)
  }))
}


# ---- Virtual BE replicate trials (power / Type I error) ----------------------
#
# A single PBPK population is one realisation. To estimate the PROBABILITY a
# real n-subject crossover trial would conclude BE, we resample subjects (with
# replacement) into many trials and record the BE pass fraction:
#   * test formulation != reference  -> that fraction is POWER.
#   * test formulation == reference   -> that fraction is the TYPE I ERROR rate,
#     but only after adding within-subject residual variability (wsv_cv): with
#     identical paired PBPK values the self-comparison is degenerate (GMR=1,
#     CI width 0). wsv_cv injects the within-subject/analytical error that a
#     real crossover carries. It is an ASSUMPTION -- set it from literature
#     (e.g. metformin within-subject CV) and report it explicitly.
#
# n_subj  : subjects per simulated trial (default = population size)
# n_trial : number of replicate trials
# wsv_cv  : within-subject CV (e.g. 0.15 = 15%); 0 = none (pure PBPK pairing)
# seed    : RNG seed for reproducibility
vbe_replicate <- function(test_vals, ref_vals, n_subj = NULL, n_trial = 1000,
                          wsv_cv = 0, seed = 1234,
                          ci_level = BE_CI_LEVEL,
                          lower = BE_LOWER, upper = BE_UPPER) {
  set.seed(seed)
  N <- length(ref_vals)
  if (is.null(n_subj)) n_subj <- N
  sw <- sqrt(log(wsv_cv^2 + 1))   # log-scale SD from CV (0 -> 0)

  pass <- logical(n_trial)
  gmrs <- numeric(n_trial)
  for (t in seq_len(n_trial)) {
    idx <- sample.int(N, n_subj, replace = TRUE)
    tv <- test_vals[idx]
    rv <- ref_vals[idx]
    if (sw > 0) {                       # independent within-subject error / period
      tv <- tv * exp(rnorm(n_subj, 0, sw))
      rv <- rv * exp(rnorm(n_subj, 0, sw))
    }
    res <- be_paired(tv, rv, ci_level = ci_level, lower = lower, upper = upper)
    pass[t] <- res$BE
    gmrs[t] <- res$GMR
  }
  data.frame(
    n_trial   = n_trial,
    n_subj    = n_subj,
    wsv_cv    = wsv_cv,
    pass_rate = mean(pass),
    GMR_med   = median(gmrs)
  )
}
