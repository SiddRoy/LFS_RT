# Entirely synthetic teaching data: no study records or identifiers.
# Run from the repository root, then source scripts in RR/ or AR/.
# This balanced illustration is not intended to mimic a clinical population.
# AR uses baseline covariates; RR additionally illustrates interval splitting.

example_inputs <- local({
  set.seed(1729)

  # Every age/sex/chemotherapy/YTA cell contains unexposed reference subjects.
  people <- expand.grid(young = 0:1, sex = 0:1, chm = 0:1, ytaN = 0:1,
                        rt_pattern = 0:3, repeat_id = 1:15)
  n <- nrow(people)
  reference <- data.frame(
    id = sprintf("synthetic_%04d", seq_len(n)),
    t10 = ifelse(people$young == 1L, 10, 30),
    sex = people$sex, chm = people$chm, ytaN = people$ytaN,
    ch = as.integer(people$rt_pattern %in% c(1L, 3L)),
    hd = as.integer(people$rt_pattern %in% c(2L, 3L)),
    pel = 0L, abs = 0L, leg = 0L, arm = 0L
  )
  sites <- c("ch", "hd", "pel", "abs", "leg", "arm")
  any_rt <- as.integer(rowSums(reference[sites]) > 0)

  # Event times are rounded to half-years to include ties. Follow-up extends
  # through age 75 for some older entrants, and includes administrative censoring.
  half_year <- function(x) pmax(0.5, round(x * 2) / 2)
  death_time <- half_year(rexp(n, rate = 0.007 * exp(
    0.25 * any_rt + 0.3 * reference$chm + 0.35 * (1 - people$young) +
      0.15 * reference$sex + 0.1 * reference$ytaN)))

  make_group <- function(group, cancer_region, site, mastectomy = FALSE) {
    in_field <- reference[[site]]
    out_field <- as.integer(any_rt == 1L & in_field == 0L)
    cancer_time <- half_year(rexp(n, rate = 0.009 * exp(
      0.75 * in_field + 0.25 * out_field + 0.2 * reference$chm +
        0.15 * reference$sex + 0.2 * reference$ytaN +
        0.2 * (1 - people$young))))
    mast_time <- if (mastectomy) half_year(rexp(n, rate = 0.009 * exp(
      0.25 * in_field + 0.15 * reference$chm + 0.2 * reference$ytaN))) else rep(Inf, n)
    event_times <- cbind(cancer_time, death_time, mast_time, censor = 45)
    event <- max.col(-event_times, ties.method = "first")
    exit <- event_times[cbind(seq_len(n), event)]
    cause <- c("cancer", "death", "mastectomy", "censored")[event]
    data.frame(reference, group = group, cancer_region = cancer_region,
               zI0 = in_field, zO0 = out_field,
               t2 = reference$t10 + exit, cause = cause)
  }

  breast <- make_group("breast", "8", "ch", mastectomy = TRUE)
  sts <- rbind(make_group("ch", "ch-8", "ch"),
               make_group("hd", "hd-10", "hd"))

  make_rr <- function(cohort) {
    # New RT or prior tumors can change the RR covariates after entry. AR
    # retains the original baseline values regardless of these later changes.
    pieces <- lapply(seq_len(nrow(cohort)), function(i) {
      person <- cohort[i, ]
      index <- match(person$id, reference$id)
      cutpoints <- sort(unique(c(person$t10, person$t2,
        person$t10 + c(4, 5, 8)[person$t10 + c(4, 5, 8) < person$t2])))
      start <- head(cutpoints, -1L)
      stop <- tail(cutpoints, -1L)
      delayed_rt <- any_rt[index] == 0L & index %% 4L == 0L &
        start >= person$t10 + 5
      # Synthetic later RT is to the chest; it is out of field for the head.
      in_field <- as.integer(person$zI0 == 1L |
        (delayed_rt & person$group != "hd"))
      out_field <- as.integer(in_field == 0L &
        (person$zO0 == 1L | delayed_rt))
      data.frame(id = person$id, region = person$cancer_region,
        t1 = start, t2 = stop,
        del = as.integer(stop == person$t2 & person$cause == "cancer"),
        zI = in_field, zO = out_field, chm = person$chm,
        sex = person$sex, ytaNA = person$ytaN,
        ben = as.integer(index %% 3L == 0L & start >= person$t10 + 4),
        mlg = as.integer(index %% 5L == 0L & start >= person$t10 + 8))
    })
    result <- do.call(rbind, pieces)
    rownames(result) <- NULL
    result
  }

  make_ar <- function(cohort, if_sites, mastectomy = FALSE) {
    causes <- c("cancer", "death", if (mastectomy) "mastectomy")
    cancer_covariates <- c("sex", "zI0", "zO0", "chm0", "t1_20", "ytaN")
    death_covariates <- c("sexD", "rtD", "chmD", "t1_20D", "ytaND")
    mast_covariates <- c("zIM", "zOM", "chmM", "t1_20M", "ytaNM")
    design <- c(cancer_covariates, death_covariates,
                if (mastectomy) mast_covariates)
    rows <- lapply(causes, function(cause) {
      region <- if (cause == "cancer") cohort$cancer_region else if (cause == "death") {
        if (mastectomy) rep("0", nrow(cohort)) else paste0(cohort$group, "-0")
      } else rep("-1", nrow(cohort))
      result <- data.frame(id = cohort$id, group = cohort$group, region = region,
                           t1 = cohort$t10, t2 = cohort$t2,
                           del = as.integer(cohort$cause == cause))
      for (name in design) result[[name]] <- 0L
      if (cause == "cancer") {
        result$sex <- cohort$sex
        result$zI0 <- cohort$zI0
        result$zO0 <- cohort$zO0
        result$chm0 <- cohort$chm
        result$t1_20 <- as.integer(cohort$t10 <= 20)
        result$ytaN <- cohort$ytaN
      } else if (cause == "death") {
        result$sexD <- cohort$sex
        result$rtD <- as.integer(rowSums(cohort[sites]) > 0)
        result$chmD <- cohort$chm
        result$t1_20D <- as.integer(cohort$t10 <= 20)
        result$ytaND <- cohort$ytaN
      } else {
        result$zIM <- cohort$zI0
        result$zOM <- cohort$zO0
        result$chmM <- cohort$chm
        result$t1_20M <- as.integer(cohort$t10 <= 20)
        result$ytaNM <- cohort$ytaN
      }
      stratified <- data.frame(id = cohort$id, group = cohort$group,
        t1 = cohort$t10, t2 = cohort$t2, del = result$del, cause = cause,
        rt = ifelse(cohort$zI0 == 1L, "In-Field",
                    ifelse(cohort$zO0 == 1L, "Out-Field", "No RT")),
        young = as.integer(cohort$t10 <= 20))
      mapping <- unique(data.frame(group = cohort$group, region = region,
                                    cause = cause))
      for (site in sites) mapping[[site]] <- 0L
      if (cause == "cancer") {
        for (group in names(if_sites)) {
          mapping[mapping$group == group, if_sites[[group]]] <- 1L
        }
      }
      # Preserve legacy prediction mapping: mastectomy has no IF site, so
      # all RT is OF in predicted profiles. Its fitted design above still
      # distinguishes actual baseline IF/OF; the stratified AR is unaffected.
      list(cox = result, stratified = stratified, map = mapping)
    })
    list(cox = do.call(rbind, lapply(rows, `[[`, "cox")),
         stratified = do.call(rbind, lapply(rows, `[[`, "stratified")),
         reference = reference,
         region_map = do.call(rbind, lapply(rows, `[[`, "map")),
         if_sites = if_sites)
  }

  list(rr_breast = make_rr(breast), rr_sts = make_rr(sts),
       ar_breast = make_ar(breast, list(breast = "ch"), mastectomy = TRUE),
       ar_sts = make_ar(sts, list(ch = "ch", hd = "hd")))
})

rr_breast <- example_inputs$rr_breast
rr_sts <- example_inputs$rr_sts
ar_breast <- example_inputs$ar_breast
ar_sts <- example_inputs$ar_sts
rm(example_inputs)
