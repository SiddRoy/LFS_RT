# Run from the repository root after source("example_data.R"), or supply ar_sts.
library(survival)
source("AR/analysis_helpers.R", local = TRUE)
if (!exists("ar_sts")) stop("First source('example_data.R') or supply prepared ar_sts input.")

run_ar_sts <- function(input, weights = NULL) {
  # 1. PRIMARY AR: fixed RT/age groups within each anatomical RT-region group.
  # Prepared first-event data have STS and death as the competing causes.
  dat <- add_subject_weights(input$stratified, weights)
  groups <- expand.grid(group = names(input$if_sites),
                        age = c("all", "t1 <= 20", "t1 > 20"),
                        rt = c("all", "No RT", "Out-Field", "In-Field"),
                        stringsAsFactors = FALSE)
  stratified <- lapply(seq_len(nrow(groups)), function(j) {
    g <- groups[j, ]
    d <- dat[dat$group == g$group, ]
    if (g$age != "all") d <- d[d$young == as.integer(g$age == "t1 <= 20"), ]
    if (g$rt != "all") d <- d[d$rt == g$rt, ]
    if (!nrow(d)) return(NULL)
    sf <- survfit(Surv(t1, t2, del) ~ cause, data = d, weights = w,
                  se.fit = FALSE, na.action = na.fail)
    ci <- cif_from_survfit(sf, causes = c("cancer", "death"),
                          start_age = if (g$age == "t1 > 20") 30 else 10)
    data.frame(group = g$group, age = g$age, rt = g$rt, ci)
  })
  stratified <- do.call(rbind, stratified)

  # 2. SECONDARY AR: common coefficients; separate cancer/death baselines by
  # region. Original subject IDs are shared across the region groups.
  # NOS cancer rows inform the pooled fit; only named anatomical groups predict.
  dat <- add_subject_weights(input$cox, weights)
  fit <- coxph(
    Surv(t1, t2, del) ~
      sex + zI0 + zO0 + chm0 + t1_20 + ytaN + t1_20:zI0 +
      sexD + rtD + chmD + ytaND + t1_20D + t1_20D:rtD + strata(region),
    data = dat, id = id, weights = w, ties = "efron", model = TRUE, na.action = na.fail
  )
  cells <- make_reference_cells(input, weights, entry_ages = c(10, 30))
  profiles <- prediction_data(cells, input$region_map)
  individual_cif <- predict_cif(fit, profiles, input$region_map, end_age = 60)
  standardized <- standardize_cif(individual_cif, cells)
  list(stratified = stratified, cox = standardized$curves,
       support = standardized$support, fit = fit, cells = cells)
}

ar_sts_results <- run_ar_sts(ar_sts)
message("STS AR: inspect ar_sts_results$stratified (primary), $cox, and $support.")

# 3. OPTIONAL: recomputes both models and empirical reference weights each time.
# ar_sts_ci <- bootstrap_ar(ar_sts, run_ar_sts, B = 1000, seed = 1729)
