# Relative risk example: soft-tissue sarcoma (Cox hazard ratios).
# From the repository root, source("example_data.R") or supply rr_sts.
# Each person can contribute several region processes. Within each RT-site
# group, retain time to first cancer; region labels identify baseline strata.
# The same subject id must be shared across all that person's intervals/groups.

library(survival)
if (!exists("rr_sts")) {
  stop("Supply rr_sts first: source('example_data.R') or your private dat_setup.R.")
}

# Common IF/OF effects are pooled over region strata. Subject-level sandwich
# variance handles the dependence between each person's regional processes.
rr_sts_fits <- list(
  basic = coxph(Surv(t1, t2, del) ~ zI + zO + chm + sex + strata(region),
                data = rr_sts, cluster = id, ties = "efron", na.action = na.fail),
  full = coxph(Surv(t1, t2, del) ~ zI + zO + chm + sex + ytaNA + ben + mlg +
                 strata(region),
               data = rr_sts, cluster = id, ties = "efron", na.action = na.fail)
)

# Wald confidence intervals and p-values use the subject-robust variance.
rr_sts_results <- do.call(rbind, lapply(names(rr_sts_fits), function(model) {
  fit <- summary(rr_sts_fits[[model]])
  data.frame(model = model, term = rownames(fit$coefficients),
             HR = fit$conf.int[, "exp(coef)"],
             lower95 = fit$conf.int[, "lower .95"],
             upper95 = fit$conf.int[, "upper .95"],
             se_log = fit$coefficients[, "robust se"],
             p = fit$coefficients[, "Pr(>|z|)"])
}))
rownames(rr_sts_results) <- NULL
print(rr_sts_results, digits = 3, row.names = FALSE)
