# Run from the LFS_RT repository with Rscript --vanilla validate_examples.R.
# Uses only simulated example data, base R, and survival; writes no files.

local({
  script <- grep("^--file=", commandArgs(), value = TRUE)
  candidates <- c(
    if (length(script)) dirname(sub("^--file=", "", script[[1L]])),
    ".", "LFS_RT"
  )
  found <- candidates[file.exists(file.path(candidates, "example_data.R"))]
  if (!length(found)) stop("Cannot find the example directory.")
  previous_directory <- getwd()
  on.exit(setwd(previous_directory), add = TRUE)
  setwd(found[[1L]])

  assert_close <- function(actual, expected, label, tolerance = 1e-8) {
    if (length(actual) != length(expected) ||
        any(!is.finite(actual)) || any(!is.finite(expected)) ||
        max(c(0, abs(as.numeric(actual) - as.numeric(expected)))) > tolerance) {
      stop(label, call. = FALSE)
    }
  }
  check_risks <- function(result, label) {
    stopifnot(is.data.frame(result), nrow(result) > 0L,
              all(c("time", "cause", "risk") %in% names(result)))
    if (any(!is.finite(result$risk)) ||
        any(result$risk < -1e-10 | result$risk > 1 + 1e-10)) {
      stop(label, ": risks are not probabilities.", call. = FALSE)
    }
    keys <- setdiff(names(result), c("time", "cause", "risk"))
    curve_keys <- c(keys, "cause")
    curves <- if (length(curve_keys)) {
      split(seq_len(nrow(result)),
            interaction(result[curve_keys], drop = TRUE, lex.order = TRUE))
    } else list(seq_len(nrow(result)))
    for (indices in curves) {
      values <- result$risk[indices[order(result$time[indices])]]
      if (any(diff(values) < -1e-9)) {
        stop(label, ": a cumulative risk decreases.", call. = FALSE)
      }
    }
    sums <- aggregate(result$risk, result[c(keys, "time")], sum)
    if (any(sums$x > 1 + 1e-8)) {
      stop(label, ": competing-cause risks exceed one.", call. = FALSE)
    }
  }

  source("AR/analysis_helpers.R", local = TRUE)
  source("example_data.R", local = TRUE)
  source("RR/RR_breast.R", local = TRUE)
  source("RR/RR_STS.R", local = TRUE)
  source("AR/AR_breast.R", local = TRUE)
  source("AR/AR_STS.R", local = TRUE)

  # A closed-form two-jump case verifies the exponential convention,
  # including a time with zero hazard and a cause with no events.
  hazard <- rbind(c(cancer = 0.2, death = 0.1, mastectomy = 0),
                  c(cancer = 0, death = 0, mastectomy = 0),
                  c(cancer = 0, death = 0.4, mastectomy = 0))
  actual <- cif_from_hazards(c(12, 13, 16), hazard)
  expected_cancer <- rep((1 - exp(-0.3)) * 2 / 3, 3L)
  expected_death <- c(rep((1 - exp(-0.3)) / 3, 2L),
                      (1 - exp(-0.3)) / 3 + exp(-0.3) * (1 - exp(-0.4)))
  assert_close(actual$survival, exp(-c(0.3, 0.3, 0.7)),
               "Exponential survival differs from the analytic result.")
  assert_close(actual$cancer, expected_cancer,
               "Cancer CIF differs from the analytic result.")
  assert_close(actual$death, expected_death,
               "Death CIF differs from the analytic result.")
  assert_close(actual$mastectomy, rep(0, 3),
               "An event-free cause has nonzero risk.")
  assert_close(rowSums(actual[-1L]), rep(1, 3),
               "Probability mass is not conserved.")

  # Weighted risk sets, delayed entry, tied events, and the strict >10
  # boundary are checked against directly calculated hazard increments.
  raw <- data.frame(
    id = 1:6, t1 = c(0, 0, 0, 11, 0, 0),
    t2 = c(10, 12, 12, 14, 14, 17),
    status = c("cancer", "cancer", "death", "censor", "death", "cancer"),
    w = c(1, 2, 3, 4, 2, 1)
  )
  causes <- c("cancer", "death", "mastectomy")
  stacked <- do.call(rbind, lapply(causes, function(cause) {
    transform(raw, cause = cause, del = as.integer(status == cause))
  }))
  sf <- survival::survfit(survival::Surv(t1, t2, del) ~ cause,
                         data = stacked, weights = w)
  observed <- cif_from_survfit(sf, causes, start_age = 10, end_age = 17)
  manual_hazard <- rbind(c(cancer = 2 / 12, death = 3 / 12, mastectomy = 0),
                         c(cancer = 0, death = 2 / 7, mastectomy = 0),
                         c(cancer = 1, death = 0, mastectomy = 0))
  totals <- rowSums(manual_hazard)
  before <- exp(-c(0, head(cumsum(totals), -1L)))
  increments <- manual_hazard * (before * (-expm1(-totals)) / totals)
  expected <- apply(increments, 2L, cumsum)
  for (cause in causes) {
    curve <- observed[observed$cause == cause, ]
    values <- curve$risk[match(c(12, 14, 17), curve$time)]
    assert_close(values, expected[, cause],
                 "Weighted stratified CIF differs from direct risk-set arithmetic.")
  }

  # The matched No-RT comparator uses the RT target weights, excludes
  # unsupported target cells, and makes that exclusion visible.
  cells <- data.frame(group = "test", entry_age = 10, ytaN = 0,
    cell_id = 1:5, W = c(100, 1, 4, 6, 1),
    rt = c("No RT", "No RT", "In-Field", "In-Field", "Out-Field"),
    match_key = c("A", "B", "A", "unsupported", "B"))
  curves <- data.frame(cell_id = 1:5, group = "test", entry_age = 10,
                        time = 20, cause = "cancer", risk = c(.1, .7, .2, .5, .3))
  standardized <- standardize_cif(curves, cells, report_support = FALSE)
  comparator <- subset(standardized$curves,
                        yta == "all" & rt == "No RT" & population == "RT")
  assert_close(comparator$risk, (.1 * 4 + .7 * 1) / 5,
               "No-RT standardization does not use the supported RT target weights.")
  diagnostic <- subset(standardized$support, yta == "all")
  stopifnot(diagnostic$unmatched_cells == 1L, diagnostic$unmatched_weight == 6)

  # Decimal attained ages must remain numeric throughout hazard alignment:
  # converting numeric time groups through printed row names loses precision.
  set.seed(902)
  decimal <- data.frame(id = seq_len(80), t1 = runif(80, 0, 8),
                         t2 = runif(80, 11, 55), x = rbinom(80, 1, .5),
                         status = sample(0:2, 80, replace = TRUE), w = rexp(80))
  decimal$t2[2] <- decimal$t2[1]
  decimal_stacked <- rbind(
    transform(decimal, region = "cancer", xC = x, xD = 0, del = status == 1),
    transform(decimal, region = "death", xC = 0, xD = x, del = status == 2))
  decimal_fit <- survival::coxph(
    survival::Surv(t1, t2, del) ~ xC + xD + strata(region),
    data = decimal_stacked, weights = w, id = id, model = TRUE)
  decimal_profile <- data.frame(cell_id = 1, group = "decimal", entry_age = 10,
    region = c("cancer", "death"), xC = c(1, 0), xD = c(0, 1))
  decimal_map <- data.frame(group = "decimal", region = c("cancer", "death"),
                             cause = c("cancer", "death"))
  decimal_cif <- predict_cif(decimal_fit, decimal_profile, decimal_map)
  check_risks(decimal_cif, "Decimal-age prediction")

  # Public RR tables must report the clustered sandwich standard errors.
  for (outcome in c("breast", "sts")) {
    fits <- get(paste0("rr_", outcome, "_fits"))
    result <- get(paste0("rr_", outcome, "_results"))
    stopifnot(all(c("basic", "full") %in% names(fits)),
              all(c("model", "term", "HR", "lower95", "upper95", "se_log", "p")
                  %in% names(result)))
    for (model in c("basic", "full")) {
      fit <- fits[[model]]
      stopifnot(!is.null(fit$naive.var), fit$method == "efron")
      expected_se <- sqrt(diag(stats::vcov(fit)))
      rows <- result[result$model == model, ]
      assert_close(rows$se_log, expected_se[rows$term],
                   "RR table does not use the clustered covariance.")
      assert_close(rows$HR, exp(stats::coef(fit)[rows$term]),
                   "RR table does not match its Cox model.")
    }
    risks <- get(paste0("ar_", outcome, "_results"))
    check_risks(risks$stratified, paste(outcome, "stratified AR"))
    check_risks(risks$cox, paste(outcome, "Cox AR"))
  }
  stopifnot("mastectomy" %in% ar_breast_results$stratified$cause,
            !"mastectomy" %in% ar_sts_results$stratified$cause,
            !"nos" %in% ar_sts_results$cox$group)

  cat("PASS: simulated examples, analytic CIFs, weighted risk sets, matched standardization, clustered RR SEs, and AR probability checks.\n")
})
