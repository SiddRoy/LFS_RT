# Numerical helpers shared by the two AR examples. Requires only base R + survival.
# The model formulas and the analysis loops remain in AR_breast.R / AR_STS.R.

rt_sites <- c("ch", "hd", "pel", "abs", "leg", "arm")

add_subject_weights <- function(dat, weights = NULL) {
  dat$w <- if (is.null(weights)) 1 else unname(weights[as.character(dat$id)])
  if (anyNA(dat$w) || any(!is.finite(dat$w) | dat$w <= 0))
    stop("Supply one finite, positive weight for every subject ID.")
  dat
}

# Legacy CIF_CR(stype = 2): exponential update, with no parametric model.
# Each row of hazard contains the nonnegative hazard jumps at a single time.
cif_from_hazards <- function(time, hazard) {
  hazard <- as.matrix(hazard)
  stopifnot(length(time) == nrow(hazard), !is.null(colnames(hazard)))
  if (anyNA(hazard) || any(!is.finite(hazard) | hazard < 0))
    stop("Hazard increments must be finite and nonnegative.")
  survival <- 1
  risk <- setNames(rep(0, ncol(hazard)), colnames(hazard))
  out <- matrix(0, nrow(hazard), ncol(hazard) + 1L)
  for (j in seq_along(time)) {
    total <- sum(hazard[j, ])
    if (total > 0) {
      risk <- risk + survival * (-expm1(-total)) * hazard[j, ] / total
      survival <- survival * exp(-total)
    }
    out[j, ] <- c(survival, risk)
  }
  colnames(out) <- c("survival", colnames(hazard))
  data.frame(time = time, out, check.names = FALSE)
}

long_cif <- function(x, causes) {
  do.call(rbind, lapply(causes, function(cause)
    data.frame(time = x$time, cause = cause, risk = x[[cause]])))
}

# survfit(~ cause) supplies cumulative cause-specific hazards. Take differences
# BEFORE restricting time; legacy curves exclude jumps exactly at start_age.
cif_from_survfit <- function(sf, causes, start_age = 10, end_age = 60) {
  lens <- if (is.null(sf$strata)) length(sf$time) else unname(sf$strata)
  labels <- if (is.null(sf$strata)) causes[1L] else sub("^.*=", "", names(sf$strata))
  idx <- split(seq_along(sf$time), rep(seq_along(lens), lens))
  jumps <- lapply(seq_along(idx), function(j) {
    k <- idx[[j]]
    data.frame(time = sf$time[k], cause = labels[j],
               hazard = c(sf$cumhaz[k[1L]], diff(sf$cumhaz[k])))
  })
  jumps <- do.call(rbind, jumps)
  times <- sort(unique(jumps$time[jumps$time > start_age & jumps$time <= end_age]))
  h <- matrix(0, length(times), length(causes), dimnames = list(NULL, causes))
  for (cause in causes) {
    z <- jumps[jumps$cause == cause & jumps$time %in% times, , drop = FALSE]
    if (nrow(z)) h[match(z$time, times), cause] <- z$hazard
  }
  # A zero row anchors every output at the specified starting age.
  long_cif(cif_from_hazards(c(start_age, times), rbind(0, h)), causes)
}

# Individual reference records are retained so bootstrap weights also update
# the empirical standardization population, rather than only the fitted model.
make_reference_cells <- function(input, weights = NULL, entry_ages = c(10, 30)) {
  ref <- add_subject_weights(input$reference, weights)
  keys <- c("sex", "chm", "ytaN", rt_sites)
  cells <- list()
  for (group in names(input$if_sites)) {
    eligible <- unique(input$cox$id[input$cox$group == group])
    for (entry_age in entry_ages) {
      d <- ref[ref$id %in% eligible & ref$t10 > entry_age - 10 &
                 ref$t10 <= entry_age + 10, , drop = FALSE]
      if (!nrow(d)) next
      z <- aggregate(d$w, d[keys], sum)
      names(z)[ncol(z)] <- "W"
      z <- z[order(-z$W), , drop = FALSE]
      z$group <- group
      z$entry_age <- entry_age
      z$rt <- ifelse(rowSums(z[rt_sites]) == 0, "No RT", "Out-Field")
      z$rt[rowSums(z[input$if_sites[[group]]]) > 0] <- "In-Field"
      cells[[length(cells) + 1L]] <- z
    }
  }
  if (!length(cells)) stop("No eligible reference subjects in the entry-age windows.")
  out <- do.call(rbind, cells)
  rownames(out) <- NULL
  out$cell_id <- seq_len(nrow(out))
  out$match_key <- do.call(paste, c(out[c("group", "entry_age", "sex", "chm", "ytaN")], sep = ":"))
  out
}

# Expand empirical cells to fitted cancer/death/mastectomy strata. Site mapping
# gives IF priority over OF, matching the original ar_cnt_fun().
prediction_data <- function(cells, region_map) {
  if (anyDuplicated(region_map$region))
    stop("Use globally unique region labels, e.g. 'ch-8' and 'hd-8'.")
  out <- list()
  design <- c("sex", "zI0", "zO0", "chm0", "t1_20", "ytaN",
              "sexD", "rtD", "chmD", "t1_20D", "ytaND",
              "zIM", "zOM", "chmM", "t1_20M", "ytaNM")
  for (group in unique(cells$group)) {
    c0 <- cells[cells$group == group, , drop = FALSE]
    r0 <- region_map[region_map$group == group, , drop = FALSE]
    if (!nrow(r0)) stop("Missing prediction strata for group ", group)
    ci <- rep(seq_len(nrow(c0)), each = nrow(r0))
    ri <- rep(seq_len(nrow(r0)), times = nrow(c0))
    c1 <- c0[ci, , drop = FALSE]
    r1 <- r0[ri, , drop = FALSE]
    any_rt <- as.integer(rowSums(c1[rt_sites]) > 0)
    in_field <- as.integer(rowSums(as.matrix(c1[rt_sites]) * as.matrix(r1[rt_sites])) > 0)
    out_field <- as.integer(any_rt == 1 & in_field == 0)
    young <- as.integer(c1$entry_age <= 20)
    z <- data.frame(cell_id = c1$cell_id, group = group,
                    entry_age = c1$entry_age, region = as.character(r1$region))
    for (v in design) z[[v]] <- 0
    cancer <- r1$cause == "cancer"
    death <- r1$cause == "death"
    mast <- r1$cause == "mastectomy"
    z$sex[cancer] <- c1$sex[cancer]
    z$zI0[cancer] <- in_field[cancer]
    z$zO0[cancer] <- out_field[cancer]
    z$chm0[cancer] <- c1$chm[cancer]
    z$t1_20[cancer] <- young[cancer]
    z$ytaN[cancer] <- c1$ytaN[cancer]
    z$sexD[death] <- c1$sex[death]
    z$rtD[death] <- any_rt[death]
    z$chmD[death] <- c1$chm[death]
    z$t1_20D[death] <- young[death]
    z$ytaND[death] <- c1$ytaN[death]
    z$zIM[mast] <- in_field[mast]
    z$zOM[mast] <- out_field[mast]
    z$chmM[mast] <- c1$chm[mast]
    z$t1_20M[mast] <- young[mast]
    z$ytaNM[mast] <- c1$ytaN[mast]
    out[[length(out) + 1L]] <- z
  }
  do.call(rbind, out)
}

# Standard survival::survfit predicts one cumulative hazard curve per profile.
# No `id` here: these rows represent simultaneous competing causes, not a
# covariate path. Combine their increments on a common attained-age timeline.
predict_cif <- function(fit, newdata, region_map, end_age = 60) {
  sf <- survival::survfit(fit, newdata = newdata, se.fit = FALSE, conf.type = "none")
  lens <- if (is.null(sf$strata)) length(sf$time) else unname(sf$strata)
  if (length(lens) != nrow(newdata) || is.matrix(sf$cumhaz))
    stop("Expected one cumulative hazard curve per prediction row.")
  idx <- split(seq_along(sf$time), rep(seq_along(lens), lens))
  map_key <- paste(region_map$group, region_map$region, sep = ":")
  cause <- region_map$cause[match(paste(newdata$group, newdata$region, sep = ":"), map_key)]
  if (anyNA(cause)) stop("A prediction stratum is missing from region_map.")
  all_jumps <- do.call(rbind, lapply(seq_along(idx), function(j) {
    k <- idx[[j]]
    data.frame(cell_id = newdata$cell_id[j], time = sf$time[k], cause = cause[j],
               hazard = c(sf$cumhaz[k[1L]], diff(sf$cumhaz[k])))
  }))
  out <- lapply(split(seq_len(nrow(newdata)), newdata$cell_id), function(rows) {
    profile <- newdata[rows[1L], , drop = FALSE]
    jumps <- all_jumps[all_jumps$cell_id == profile$cell_id &
                         all_jumps$time > profile$entry_age & all_jumps$time <= end_age, ]
    causes <- unique(cause[rows])
    times <- sort(unique(jumps$time))
    h <- matrix(0, length(times), length(causes), dimnames = list(NULL, causes))
    for (cc in causes) {
      z <- jumps[jumps$cause == cc, , drop = FALSE]
      if (nrow(z)) {
        # Integer positions avoid losing precision when decimal ages are
        # converted into rowsum's character row names.
        sums <- rowsum(z$hazard, match(z$time, times), reorder = FALSE)
        h[as.integer(rownames(sums)), cc] <- sums[, 1L]
      }
    }
    z <- long_cif(cif_from_hazards(c(profile$entry_age, times), rbind(0, h)), causes)
    cbind(profile[c("cell_id", "group", "entry_age")], z, row.names = NULL)
  })
  do.call(rbind, out)
}

# Preserve legacy targets: observed No-RT population; IF and OF populations;
# and No-RT cells matched to ALL RT cells on age, sex, chemo and YTA.
standardize_cif <- function(curves, cells, report_support = TRUE) {
  out <- support <- list()
  for (group in unique(cells$group)) for (entry_age in unique(cells$entry_age)) {
    c0 <- cells[cells$group == group & cells$entry_age == entry_age, , drop = FALSE]
    if (!nrow(c0)) next
    for (yta in c("all", "B+", "NA")) {
      c1 <- c0
      if (yta != "all") c1 <- c1[c1$ytaN == as.integer(yta == "NA"), , drop = FALSE]
      if (!nrow(c1)) next
      target <- c1[c1$rt != "No RT", , drop = FALSE]
      controls <- c1[c1$rt == "No RT", , drop = FALSE]
      control_id <- controls$cell_id[match(target$match_key, controls$match_key)]
      missing <- is.na(control_id)
      support[[length(support) + 1L]] <- data.frame(group, entry_age, yta,
        rt_cells = nrow(target), unmatched_cells = sum(missing),
        rt_weight = sum(target$W), unmatched_weight = sum(target$W[missing]))
      targets <- list(
        list(d = controls, ids = controls$cell_id, rt = "No RT", population = "No RT"),
        list(d = c1[c1$rt == "In-Field", ], ids = c1$cell_id[c1$rt == "In-Field"],
             rt = "In-Field", population = "RT"),
        list(d = c1[c1$rt == "Out-Field", ], ids = c1$cell_id[c1$rt == "Out-Field"],
             rt = "Out-Field", population = "RT"),
        list(d = target[!missing, ], ids = control_id[!missing], rt = "No RT", population = "RT")
      )
      for (spec in targets) {
        if (!nrow(spec$d)) next
        z <- merge(data.frame(cell_id = spec$ids, W = spec$d$W), curves, by = "cell_id")
        z$weighted <- z$risk * z$W
        a <- aggregate(z[c("weighted", "W")], z[c("time", "cause")], sum)
        a$risk <- a$weighted / a$W
        out[[length(out) + 1L]] <- data.frame(group, entry_age, yta,
          rt = spec$rt, population = spec$population, a[c("time", "cause", "risk")])
      }
    }
  }
  diagnostic <- do.call(rbind, support)
  if (report_support && any(diagnostic$unmatched_cells > 0))
    message("Some RT reference cells lack a matched No-RT cell; see results$support. ",
            "Their weights are omitted from that comparator, as in the legacy analysis.")
  list(curves = do.call(rbind, out), support = diagnostic)
}

# Subject-level exponential multiplier bootstrap, as in make_cW0(boot=TRUE).
# `analysis` reruns BOTH estimators and recomputes reference-cell weights.
bootstrap_summary <- function(observed, B, read_replicate) {
  estimates <- lapply(c("stratified", "cox"), function(nm) {
    x <- observed[[nm]]
    list(data = x, keys = setdiff(names(x), "risk"), values = matrix(NA_real_, nrow(x), B))
  })
  names(estimates) <- c("stratified", "cox")
  failed <- character()
  for (b in seq_len(B)) {
    result <- read_replicate(b)
    if (inherits(result, "error")) {
      failed <- c(failed, paste0(b, ": ", conditionMessage(result)))
      next
    }
    for (nm in names(estimates)) {
      e <- estimates[[nm]]
      key <- function(d) do.call(paste, c(d[e$keys], sep = "\r"))
      value <- result[[nm]]$risk[match(key(e$data), key(result[[nm]]))]
      estimates[[nm]]$values[, b] <- value
    }
  }
  for (nm in names(estimates)) {
    e <- estimates[[nm]]
    bounds <- t(apply(e$values, 1L, function(v) {
      if (all(is.na(v))) return(c(NA_real_, NA_real_))
      quantile(v, c(.025, .975), na.rm = TRUE, names = FALSE)
    }))
    estimates[[nm]] <- cbind(e$data, lower = bounds[, 1L], upper = bounds[, 2L],
                            successful_replicates = rowSums(!is.na(e$values)))
  }
  if (length(failed)) warning(length(failed), " bootstrap replicates failed; see $failures.")
  c(estimates, list(failures = failed, B = B))
}

bootstrap_ar <- function(input, analysis, B = 0L, seed = 1729L, output_dir = NULL) {
  stopifnot(length(B) == 1L, is.finite(B), B >= 0, B == as.integer(B))
  if (!B) return(NULL)
  ids <- unique(as.character(input$reference$id))
  array_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", "")
  replicates <- seq_len(B)
  if (nzchar(array_id)) {
    count <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT", "1"))
    first <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_MIN", array_id))
    step <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_STEP", "1"))
    last <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_MAX", as.character(first + (count - 1L) * step)))
    rank <- (as.integer(array_id) - first) / step
    if (anyNA(c(count, rank, step, last)) || count < 1 || step < 1 ||
        last != first + (count - 1L) * step ||
        rank < 0 || rank >= count || rank != as.integer(rank))
      stop("Use a contiguous or regularly stepped SLURM array with valid index/count/min/step.")
    replicates <- replicates[(replicates - 1L) %% count == rank]
    if (is.null(output_dir))
      stop("For SLURM, supply output_dir, e.g. 'results/bootstrap_breast'.")
  }
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(output_dir)) stop("Cannot create bootstrap output directory.")
    paths <- file.path(output_dir, sprintf("replicate_%05d.rds", replicates))
    if (any(file.exists(paths))) stop("Bootstrap files already exist; choose a new output_dir.")
  }
  one_replicate <- function(b) {
    # Same replicate seed in a serial run and in any SLURM array partition.
    set.seed(seed + b)
    weights <- setNames(rexp(length(ids)), ids)
    result <- tryCatch(suppressMessages(analysis(input, weights = weights)), error = identity)
    if (!inherits(result, "error")) result <- result[c("stratified", "cox")]
    if (!is.null(output_dir)) saveRDS(list(replicate = b, B = B, seed = seed,
      result = result), file.path(output_dir, sprintf("replicate_%05d.rds", b)))
    result
  }
  if (nzchar(array_id)) {
    for (b in replicates) one_replicate(b)
    return(invisible(list(replicates = replicates, files = paths, B = B, seed = seed)))
  }
  observed <- analysis(input)
  c(bootstrap_summary(observed, B, one_replicate), list(seed = seed))
}

# After every array task completes, combine the per-replicate files into CIs.
collect_bootstrap_ar <- function(input, analysis, output_dir) {
  files <- list.files(output_dir, pattern = "^replicate_[0-9]+[.]rds$", full.names = TRUE)
  if (!length(files)) stop("No bootstrap replicate files found.")
  header <- readRDS(files[1L])
  expected <- file.path(output_dir, sprintf("replicate_%05d.rds", seq_len(header$B)))
  if (!all(file.exists(expected))) stop("Some bootstrap replicates are still missing.")
  read_replicate <- function(b) {
    x <- readRDS(expected[b])
    if (!identical(x$B, header$B) || !identical(x$seed, header$seed) || x$replicate != b)
      stop("Bootstrap files mix different runs.")
    x$result
  }
  c(bootstrap_summary(analysis(input), header$B, read_replicate), list(seed = header$seed))
}
