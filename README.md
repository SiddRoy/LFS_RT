# Subsequent Neoplasm Risk in Individuals with Germline TP53 Variants

Main analysis code 

Title: Risk of Therapy-Related Subsequent Neoplasms in Individuals with Germline TP53 Variants, Nature Communication 2026

**Authors**: Anita Villani M.D., M.Sc.<sup>&ast;</sup>, Payal P. Khincha, M.B.B.S., M.S.H.S.<sup>&ast;</sup>, Siddharth Roy, Ph.D., Jessica N. Hatton, M.S., C.G.C., Sara J. Schonfeld, Ph.D., Catherine Goudie, M.D., M.Sc., Leatrisse Oba, M.D., Megan N. Frone, M.S., C.G.C., Kishan Pithadia, M.D., Birte Sänger, B.A., Anh N. Le, D.O., Caitlin Orr, Ph.D., Anne Naumer, C.G.C., Jennie Vagher, C.G.C., Wendy K. Kohlmann, M.S., C.G.C., Orli Michaeli, M.D., Ron Rabinowicz, M.D., Jacob J. K. Carson, M.D., Jamie L. Maciaszek, Ph.D., Victoria E. Goodman, M.P.H., Elizabeth Santana dos Santos, M.D., Lara Reichman, C.G.C., Paul S. Albert, Ph.D., Suzanne P. MacFarland, M.D., Kim E. Nichols, M.D., Derek S. Tsang, M.D., M.Sc., Raymond H. Kim, M.D., Ph.D., Joshua Schiffman, M.D., Kara N. Maxwell, M.D., Ph.D., Christian Kratz, M.D., Maria Isabel Achatz, M.D., Ph.D., Sharon A. Savage, M.D., David Malkin, M.D., Lindsay M. Morton, Ph.D.

<sup>&ast;</sup> These authors contributed equally.

## Relative and absolute risk examples in R

Four short scripts demonstrate breast cancer and soft-tissue sarcoma (STS) analyses. `RR/RR_breast.R` and `RR/RR_STS.R` fit basic and full Cox models. Each `AR/AR_*.R` script shows the primary stratified analysis at the top and the secondary Cox/g-formula analysis below it. The numerical AR operations are in `AR/analysis_helpers.R` so the model formulas remain easy to read.

The public examples need only **base R and `survival`**, which is included in standard R installations. `example_data.R` creates entirely fictional data in memory; it contains no study records. No example writes files automatically.


The analysis files are organized by method; simulated data and validation stay at the repository root:

```text
LFS_RT/
├── RR/
│   ├── RR_breast.R
│   └── RR_STS.R
├── AR/
│   ├── AR_breast.R
│   ├── AR_STS.R
│   └── analysis_helpers.R
├── example_data.R
├── validate_examples.R
├── README.md
└── .gitignore
```

### Run the public examples

Start R with the `LFS_RT` repository root as the working directory, and keep it there when sourcing the scripts:

```r
source("example_data.R")
source("RR/RR_breast.R")
source("RR/RR_STS.R")
source("AR/AR_breast.R")
source("AR/AR_STS.R")

rr_breast_results
rr_sts_results
str(ar_breast_results, max.level = 1)
str(ar_sts_results, max.level = 1)
```

The RR tables contain `model`, `term`, `HR`, `lower95`, `upper95`, `se_log` and `p`. The fitted models are retained in `rr_breast_fits` and `rr_sts_fits`, each with `basic` and `full` entries. RR is shorthand here: Cox estimates are **hazard ratios**, rather than ratios of cumulative risks.

Each AR result is a list containing `stratified` (primary risk curves), `cox` (standardized secondary curves), `support` (matching diagnostics), `fit` (the Cox model), and `cells` (reference-population cells). Curves have `time`, `cause` and `risk` columns; time is attained age in years and risk is a probability.

For example, inspect cancer risk at the last available time at or before age 45, for predictions starting at age 30:

```r
x <- subset(ar_breast_results$cox,
            cause == "cancer" & entry_age == 30 & yta == "all" & time <= 45)
x <- x[order(x$time), ]
x[!duplicated(x[c("rt", "population")], fromLast = TRUE),
  c("rt", "population", "time", "risk")]
```

Run the public numerical checks from a terminal in this folder:

```sh
Rscript --vanilla validate_examples.R
```

They check analytic CIF calculations, weighted risk sets, delayed entry, ties, boundary handling, clustered RR standard errors, and valid AR probabilities. They write no files.

### Supply your own prepared data

Replace the `source("example_data.R")` step with your own preparation or explicit data loading. For example, a reader who has prepared these R objects can use:

```r
rr_breast <- readRDS("private/rr_breast.rds")
rr_sts <- readRDS("private/rr_sts.rds")
ar_breast <- readRDS("private/ar_breast.rds")
ar_sts <- readRDS("private/ar_sts.rds")
```

The analysis scripts use the objects already in the session. They do not find or load private data themselves. Inspect the synthetic input structures with `str(rr_breast)`, `str(ar_breast, max.level = 2)` and their STS equivalents.

#### RR data frames: `rr_breast` and `rr_sts`

| Columns | Type and meaning |
| --- | --- |
| `id`, `region` | Character subject ID and baseline-hazard stratum. Reuse the same ID across all records from one person. |
| `t1`, `t2` | Numeric attained ages defining a positive-length interval `(t1, t2]`; delayed entry is allowed. |
| `del` | Numeric/integer 0 or 1; a cancer event at the interval end. Death and other relevant endpoints censor this cancer process. |
| `zI`, `zO` | 0/1 in-field and out-of-field RT at that interval. They are mutually exclusive; both zero is the reference. |
| `chm`, `sex` | 0/1 chemotherapy and sex covariates. Use the same coding consistently throughout the inputs. |
| `ytaNA` | Fixed 0/1 YTA group; 1 denotes null/A, 0 the other categories. |
| `ben`, `mlg` | 0/1 prior benign and malignant tumor indicators, updated over intervals. |

Create interval splits before running RR whenever a predictor changes. Within each subject/process, intervals must not overlap or continue after its first cancer event. A subject can contribute several STS processes; their region labels distinguish the baseline strata.

Both basic models use `zI + zO + chm + sex + strata(region)`; full models additionally use `ytaNA + ben + mlg`. They pool coefficients across region strata and use Efron ties with subject-clustered sandwich standard errors.

#### AR input lists: `ar_breast` and `ar_sts`

Each input has these five named components. All predictors used for AR are **fixed at baseline**; do not supply changing predictor paths.

| Component | Required contents |
| --- | --- |
| `cox` | Cause-stacked data frame with `id`, `group`, `region`, `t1`, `t2`, `del` and the predictor blocks below. Region labels identify distinct cause-specific baseline hazards within each group. |
| `stratified` | Cause-stacked data frame with character `id`, `group`, `cause`, `rt`; numeric `t1`, `t2`; 0/1 `del`, `young`. `cause` is `cancer`, `death`, or breast-only `mastectomy`; `rt` is `No RT`, `Out-Field`, or `In-Field`; `young` indicates baseline age ≤20. |
| `reference` | One row per subject: character `id`; numeric `t10` (baseline age); 0/1 `sex`, `chm`, `ytaN`, `ch`, `hd`, `pel`, `abs`, `leg`, `arm`. The last six fields indicate baseline RT to chest, head, pelvis, abdomen, leg and arm. |
| `region_map` | One row per `group`/`region`: those two character fields, character `cause`, and six 0/1 site-membership fields named as above. Death and breast mastectomy rows have all memberships zero; see the compatibility convention below. |
| `if_sites` | Named list giving the in-field site column(s) for each group, e.g. `list(breast = "ch")` or `list(ch = "ch", hd = "hd")`. These groups receive AR outputs. |

The `cox` predictor blocks are:

- Cancer: `sex`, `zI0`, `zO0`, `chm0`, `t1_20`, `ytaN`.
- Death: `sexD`, `rtD`, `chmD`, `t1_20D`, `ytaND`.
- Breast mastectomy: `zIM`, `zOM`, `chmM`, `t1_20M`, `ytaNM`.

Each block is zero on rows for other causes. IF/OF indicators are exclusive, `rtD` denotes any RT, `t1_20*` indicates baseline age ≤20, and YTA indicators use the same 0/1 coding as `ytaNA`. The breast Cox formula follows the legacy model and omits sex. Shared original subject IDs must remain consistent across every component. Numeric `del` is 0/1 for the row's cause; a competing event is not a cancer event.

Region labels must be globally unique across groups, for example `ch-8` and `hd-8`. Reusing `8` for both groups would merge their baseline hazards. Input fields used in a model must contain no missing values; the displayed fits fail rather than silently omit records.

### What the AR estimators calculate

The primary analysis stratifies baseline RT and baseline age groups, estimates separate cause-specific cumulative hazards with `survfit`, and combines their hazard jumps. Breast includes cancer, death and bilateral mastectomy; STS includes cancer and death within each RT-site group. Baseline age groups are all, ≤20 and >20. Curves include times strictly **greater than 10** for all/younger groups or **greater than 30** for older groups, through age 60.

Both AR approaches preserve the legacy **exponential CIF** convention: at each time, survival is multiplied by `exp(-sum(hazard_jumps))`; the resulting probability decrease is apportioned to the causes by their hazard shares. This differs from an Aalen–Johansen product-limit update, especially with tied events. The scripts deliberately use this convention for reproducibility.

The secondary analysis fits the displayed baseline Cox model, predicts cause-specific hazards, combines them into individual CIFs, and averages over empirical reference cells. Predictions start at age 10 or 30 and end at 60. Their reference populations include baseline ages in `(0, 20]` and `(20, 40]`, respectively; weights retain the observed sex, chemotherapy, YTA and RT-site patterns. Results are shown for YTA `all`, `B+` (`ytaN = 0`), and `NA` (`ytaN = 1`, the null/A label, not a missing value).

For numerical compatibility, breast mastectomy predictions preserve the legacy convention that every RT pattern is assigned out of field: predicted `zIM = 0` and `zOM = any RT`. The fitted mastectomy model still uses actual baseline IF/OF indicators. This difference between fitting and prediction is intentional here to reproduce the legacy calculation; it does not affect the primary stratified analysis.

The four standardized targets are observed No-RT, observed in-field RT, observed out-of-field RT, and No-RT standardized to the combined RT population. The last matches RT cells to observed No-RT cells on group, entry age, sex, chemotherapy and YTA. Unmatched RT cells are omitted from that comparator; `$support` reports their number and weight. Review it before comparing curves.

#### Optional uncertainty estimates

Bootstrap is off by default. After an AR script has run:

```r
ar_breast_ci <- bootstrap_ar(ar_breast, run_ar_breast, B = 1000, seed = 1729)
```

This preserves the legacy exponential multiplier bootstrap: replicate `b` uses seed `1729 + b`, draws one `Exp(1)` weight per subject, and reuses it across all records, competing causes, groups and empirical reference cells. Both AR estimators and reference weights are recomputed. Results contain percentile `lower`/`upper` bounds, counts of successful replicates, and a `failures` log. A full run can take substantial time; the quickstart does not run it.

For a SLURM array, run the same setup and AR script in each task, then call:

```r
bootstrap_ar(ar_breast, run_ar_breast, B = 1000, seed = 1729,
             output_dir = "results/bootstrap_breast")
```

The helper reads `SLURM_ARRAY_TASK_ID`, `SLURM_ARRAY_TASK_COUNT`, `SLURM_ARRAY_TASK_MIN`, `SLURM_ARRAY_TASK_MAX` and `SLURM_ARRAY_TASK_STEP`, assigns replicates across tasks, and saves one `replicate_00001.rds`-style file per replicate. Contiguous and regularly stepped arrays are supported; irregular sparse lists are rejected. Seeds are identical to serial runs. An array task returns its file list; after **all** tasks complete, combine the results in a separate R session:

```r
ar_breast_ci <- collect_bootstrap_ar(ar_breast, run_ar_breast,
                                      "results/bootstrap_breast")
```

Use a distinct directory such as `results/bootstrap_sts` for STS. Choose a new directory for another run: existing replicate files are never overwritten. Serial runs can also supply `output_dir` to retain individual replicates; the default `NULL` writes nothing.

### Private preparation and publication

The private `dat_setup.R` is kept outside this public repository. It reproduces the existing study preparation using the private `lfsfun` package and read-only legacy helpers. In the original local checkout, start R in the parent `github_new_examples/` folder and run:

```r
source("dat_setup.R")
setwd("LFS_RT")
# Now source RR/RR_breast.R, RR/RR_STS.R, AR/AR_breast.R or AR/AR_STS.R.
```

If the package is in a separate R library, set `R_LIBS_USER` to that library directory before starting R. Private setup also needs the packages loaded by the legacy helpers; the public scripts do not. Public users should prepare the documented input objects or use `example_data.R`.

Private preparation retains IF/OF region-specific analyses: breast uses its first event, while STS takes the first event separately within each configured RT-site aggregate, then combines the processes. This is the legacy teaching example, distinct from the newer pooled-all-regions/any-RT analysis. It freezes AR exposures at entry and preserves legacy cancer/death censoring. Breast preparation includes the legacy bilateral-mastectomy endpoint and its `mast_age + 1` timing rule. Those source-specific fields are not required public inputs.

Publish only `RR/RR_breast.R`, `RR/RR_STS.R`, `AR/AR_breast.R`, `AR/AR_STS.R`, `AR/analysis_helpers.R`, `example_data.R`, `validate_examples.R`, `README.md` and `.gitignore`. Do not publish `dat_setup.R`, `validation_private.R`, data, or generated results. The private validation script remains outside this repository and checks numerical agreement with the legacy analysis in the original local checkout.

To save results explicitly, use a local ignored directory inside this folder, for example `dir.create("results", showWarnings = FALSE)` followed by `saveRDS(ar_breast_results, "results/ar_breast_results.rds")`.
