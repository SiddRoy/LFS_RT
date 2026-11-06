# Run from the repository root after source("example_data.R"), or supply ar_breast.
library(survival)
source("AR/analysis_helpers.R", local = TRUE)
if (!exists("ar_breast")) stop("First source('example_data.R') or supply prepared ar_breast input.")

run_ar_breast <- function(input, weights = NULL) {
  # 1. PRIMARY AR: fixed RT groups, cancer + death + bilateral mastectomy.
  # All event timing/censoring is already prepared; no time-varying predictors.
  dat <- add_subject_weights(input$stratified, weights)
  groups <- expand.grid(age = c("all", "t1 <= 20", "t1 > 20"),
                        rt = c("all", "No RT", "Out-Field", "In-Field"),
                        stringsAsFactors = FALSE)
  stratified <- lapply(seq_len(nrow(groups)), function(j) {
    g <- groups[j, ]
    d <- dat
    if (g$age != "all") d <- d[d$young == as.integer(g$age == "t1 <= 20"), ]
    if (g$rt != "all") d <- d[d$rt == g$rt, ]
    if (!nrow(d)) return(NULL)
    sf <- survfit(Surv(t1, t2, del) ~ cause, data = d, weights = w,
                  se.fit = FALSE, na.action = na.fail)
    ci <- cif_from_survfit(sf, causes = c("cancer", "death", "mastectomy"),
                          start_age = if (g$age == "t1 > 20") 30 else 10)
    data.frame(group = "breast", age = g$age, rt = g$rt, ci)
  })
  stratified <- do.call(rbind, stratified)

  # 2. SECONDARY AR: baseline Cox models, followed by g-formula averaging.
  # Each predictor block is zero on rows for the other causes. As in the
  # legacy breast/mastectomy model, sex is not included in these formulas.
  dat <- add_subject_weights(input$cox, weights)
  fit <- coxph(
    Surv(t1, t2, del) ~
      zI0 + zO0 + chm0 + t1_20 + ytaN + t1_20:zI0 +
      rtD + chmD + ytaND + t1_20D + t1_20D:rtD +
      chmM + zIM + zOM + ytaNM + t1_20M + t1_20M:zIM + strata(region),
    data = dat, id = id, weights = w, ties = "efron", model = TRUE, na.action = na.fail
  )
  cells <- make_reference_cells(input, weights, entry_ages = c(10, 30))
  profiles <- prediction_data(cells, input$region_map)
  individual_cif <- predict_cif(fit, profiles, input$region_map, end_age = 60)
  standardized <- standardize_cif(individual_cif, cells)
  list(stratified = stratified, cox = standardized$curves,
       support = standardized$support, fit = fit, cells = cells)
}

ar_breast_results <- run_ar_breast(ar_breast)
message("Breast AR: inspect ar_breast_results$stratified (primary), $cox, and $support.")

# 3. OPTIONAL: one subject weight shared by all causes AND reference cells.
# ar_breast_ci <- bootstrap_ar(ar_breast, run_ar_breast, B = 1000, seed = 1729)
