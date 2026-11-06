# Relative risk example: breast cancer (Cox hazard ratios).
# From the repository root, source("example_data.R") or supply rr_breast.
# Each row is an attained-age interval (t1, t2] for one subject and region.
# zI/zO: current in-field/out-of-field RT; joint zero is the reference.
# ytaNA is fixed; ben and mlg are current prior-tumor indicators.

library(survival)
if (!exists("rr_breast")) {
  stop("Supply rr_breast first: source('example_data.R') or your private dat_setup.R.")
}

# Region strata allow separate baseline hazards and a common exposure effect.
# Subject clustering accounts for repeated intervals belonging to one person.
rr_breast_fits <- list(
  basic = coxph(Surv(t1, t2, del) ~ zI + zO + chm + sex + strata(region),
                data = rr_breast, cluster = id, ties = "efron", na.action = na.fail),
  full = coxph(Surv(t1, t2, del) ~ zI + zO + chm + sex + ytaNA + ben + mlg +
                 strata(region),
               data = rr_breast, cluster = id, ties = "efron", na.action = na.fail)
)

# Wald confidence intervals and p-values use the subject-robust variance.
rr_breast_results <- do.call(rbind, lapply(names(rr_breast_fits), function(model) {
  fit <- summary(rr_breast_fits[[model]])
  data.frame(model = model, term = rownames(fit$coefficients),
             HR = fit$conf.int[, "exp(coef)"],
             lower95 = fit$conf.int[, "lower .95"],
             upper95 = fit$conf.int[, "upper .95"],
             se_log = fit$coefficients[, "robust se"],
             p = fit$coefficients[, "Pr(>|z|)"])
}))
rownames(rr_breast_results) <- NULL
print(rr_breast_results, digits = 3, row.names = FALSE)
