#!/usr/bin/env Rscript
# =============================================================================
# 39_fit_B_prime_L2_and_specificity_v6.R
# -----------------------------------------------------------------------------
# 処理内容 (v5 -> v6):
#   - [path update only] Option II cascade により Step 05 が v4 -> v5 に更新された
#     ことに伴う path 切替。統計モデル・カラム生成ロジックは v5 と完全一致。
#     * PATTERN="patB" or "patC" を環境変数で受け取る (必須)
#     * 入力: 05_wgs_sample_burden/output_v5_{PATTERN}/
#               sample_burden_L2_and_specificity_v5_{PATTERN}.tsv
#     * 出力: 06_wgs_primary_L2/output_v6_{PATTERN}/
#               B_prime_L2_classes_results_v6_{PATTERN}.tsv
#               B_prime_specificity_groups_results_v6_{PATTERN}.tsv
#               B_prime_run_summary_v6_{PATTERN}.tsv
#     * manifest JSON (paths_v1_manifest.json) は L2_CLASSES 整合性チェック用に
#       引き続き参照 (paths_v1.py は未修正)
#   - 統計モデル・ロジック (10 L2 classes, 5 specificity groups, 3 exposure,
#     3 comparison, 2 SV type, B' covariate 構成, FDR) は v5 と完全一致
# 使い方:
#   export PATTERN="patB"   # or "patC"
#   Rscript 39_fit_B_prime_L2_and_specificity_v6.R \
#     [--burden <sample_burden_L2_and_specificity_v5_patX.tsv>] \
#     [--outdir <outdir>] \
#     [--min-cell-count 5]
# 出力:
#   B_prime_L2_classes_results_v6_{PATTERN}.tsv         (解析 i)
#   B_prime_specificity_groups_results_v6_{PATTERN}.tsv (解析 ii)
#   B_prime_run_summary_v6_{PATTERN}.tsv                (fit 成否サマリ)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(jsonlite)
})

# -------------------------------------------------------------
# utility
# -------------------------------------------------------------
log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              msg))
  flush.console()
}

# -------------------------------------------------------------
# minimal argparse (base R)
# -------------------------------------------------------------
parse_cli_args <- function(raw) {
  out <- list()
  i <- 1
  while (i <= length(raw)) {
    tok <- raw[i]
    if (!startsWith(tok, "--")) {
      stop(sprintf("Unexpected positional argument: '%s'", tok))
    }
    key <- sub("^--", "", tok)
    if (grepl("=", key, fixed = TRUE)) {
      kv <- strsplit(key, "=", fixed = TRUE)[[1]]
      out[[kv[1]]] <- paste(kv[-1], collapse = "=")
      i <- i + 1
    } else {
      if (i + 1 > length(raw)) {
        stop(sprintf("Argument --%s requires a value.", key))
      }
      out[[key]] <- raw[i + 1]
      i <- i + 2
    }
  }
  out
}

# -------------------------------------------------------------
# v6 (2026-04-22): PATTERN env var で patB / patC を切替
# DEFAULT_BURDEN / DEFAULT_OUTDIR は v6_{PATTERN} 用にハードコード override
# (paths_v1.py を修正しないため)
# Option II cascade: Step 05 v5 (silent-sample-loss fix) 出力を入力とする
# -------------------------------------------------------------
PATTERN <- Sys.getenv("PATTERN", unset = "")
if (!(PATTERN %in% c("patB", "patC"))) {
  stop(sprintf(
    "v6 requires PATTERN env var set to 'patB' or 'patC'. Got: '%s'",
    PATTERN))
}
V6_SUFFIX <- sprintf("v6_%s", PATTERN)

PIPELINE_ROOT  <- "/lustre12/home/kushima-pg/tad04212026"
MANIFEST_PATH  <- file.path(PIPELINE_ROOT, "common", "paths_v1_manifest.json")
if (!file.exists(MANIFEST_PATH)) {
  stop(sprintf(
    "paths_v1_manifest.json not found at %s. Run 'python3 -c \"from common.paths_v1 import export_r_manifest; export_r_manifest()\"' first.",
    MANIFEST_PATH))
}
manifest <- jsonlite::fromJSON(MANIFEST_PATH, simplifyVector = TRUE)

# v6 hardcoded override (Pattern B/C: 05_v5_{PATTERN} burden -> 06_v6_{PATTERN} outdir)
DEFAULT_BURDEN <- file.path(PIPELINE_ROOT,
                            "05_wgs_sample_burden",
                            sprintf("output_v5_%s", PATTERN),
                            sprintf("sample_burden_L2_and_specificity_v5_%s.tsv", PATTERN))
DEFAULT_OUTDIR <- file.path(PIPELINE_ROOT,
                            "06_wgs_primary_L2",
                            sprintf("output_v6_%s", PATTERN))

raw_args <- commandArgs(trailingOnly = TRUE)
opt <- if (length(raw_args) > 0) parse_cli_args(raw_args) else list()
burden_path <- ifelse(is.null(opt$burden),  DEFAULT_BURDEN, opt$burden)
outdir      <- ifelse(is.null(opt$outdir),  DEFAULT_OUTDIR, opt$outdir)
min_cell    <- as.integer(ifelse(is.null(opt[["min-cell-count"]]),
                                 5L, opt[["min-cell-count"]]))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------
# constants
# -------------------------------------------------------------
L2_CLASSES_HARDCODE <- c(
  "HPC_Exc-CA", "HPC_Exc-DG", "HPC_Exc-ENT",
  "HPC_Inh-CGE", "HPC_Inh-MGE",
  "PFC_Astro", "PFC_Exc-DL", "PFC_Exc-UL",
  "PFC_Inh-CGE", "PFC_Inh-MGE"
)
stopifnot(
  "Manifest L2_CLASSES mismatch with R hardcoded list" =
    identical(sort(L2_CLASSES_HARDCODE), sort(as.character(manifest$L2_CLASSES)))
)
L2_CLASSES <- as.character(manifest$L2_CLASSES)

GROUP_COLUMNS <- c("group_primary", "group_s2", "group_s3", "group_s4", "group_s5")
SV_TYPES      <- c("DEL", "DUP")
EXPOSURES     <- c("n_boundary", "n_events", "carrier_boundary")

COMPARISONS <- list(
  list(name = "ASD_vs_HC", case = "ASD",     control = "Healthy"),
  list(name = "SZ_vs_HC",  case = "SZ",      control = "Healthy"),
  list(name = "ASD_vs_SZ", case = "ASD",     control = "SZ")
)

FIXED_COVS <- c(
  "Sex_numeric",
  "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"
)

t0 <- Sys.time()
log_msg(strrep("=", 60))
log_msg(sprintf("Start 39_fit_B_prime_L2_and_specificity_v6.R (PATTERN=%s)", PATTERN))
log_msg(sprintf("  burden: %s", burden_path))
log_msg(sprintf("  outdir: %s", outdir))
log_msg(sprintf("  min_cell_count: %d", min_cell))

# -------------------------------------------------------------
# load burden
# -------------------------------------------------------------
log_msg("Loading burden table ...")
bur <- read_tsv(burden_path, show_col_types = FALSE, guess_max = 50000)
log_msg(sprintf("  dim = (%d, %d)", nrow(bur), ncol(bur)))

required_cov <- c("sample_id", "Diagnosis",
                  FIXED_COVS,
                  "log1p_total_del_bases", "log1p_total_dup_bases",
                  "log1p_total_gene_DEL",  "log1p_total_gene_DUP")
missing_cov <- setdiff(required_cov, colnames(bur))
if (length(missing_cov) > 0) {
  stop(sprintf("Burden missing required covariates: %s",
               paste(missing_cov, collapse = ", ")))
}

dx_tab <- table(bur$Diagnosis, useNA = "ifany")
log_msg("Diagnosis distribution:")
for (d in names(dx_tab)) {
  log_msg(sprintf("  %s : %d", d, dx_tab[[d]]))
}

# -------------------------------------------------------------
# helpers
# -------------------------------------------------------------
make_subset <- function(df, comp) {
  sub <- df %>% filter(Diagnosis %in% c(comp$case, comp$control))
  sub$case <- as.integer(sub$Diagnosis == comp$case)
  sub
}

covariate_formula_terms <- function(sv_type) {
  base <- FIXED_COVS
  if (sv_type == "DEL") {
    c(base, "log1p_total_del_bases", "log1p_total_gene_DEL")
  } else if (sv_type == "DUP") {
    c(base, "log1p_total_dup_bases", "log1p_total_gene_DUP")
  } else {
    stop(sprintf("Unknown sv_type: %s", sv_type))
  }
}

fit_logit <- function(sub, exposure_col, sv_type, min_cell = 5) {
  covs <- covariate_formula_terms(sv_type)
  needed <- c("case", exposure_col, covs)
  if (!all(needed %in% colnames(sub))) {
    return(tibble(
      n_case = NA_integer_, n_ctrl = NA_integer_,
      carrier_case = NA_integer_, carrier_ctrl = NA_integer_,
      mean_exp_case = NA_real_, mean_exp_ctrl = NA_real_,
      beta = NA_real_, se = NA_real_, z = NA_real_,
      p_value = NA_real_, or = NA_real_,
      or_lo95 = NA_real_, or_hi95 = NA_real_,
      fit_status = "missing_columns"
    ))
  }
  d <- sub[, needed, drop = FALSE]
  d <- d[complete.cases(d), , drop = FALSE]
  if (nrow(d) == 0) {
    return(tibble(
      n_case = 0L, n_ctrl = 0L,
      carrier_case = 0L, carrier_ctrl = 0L,
      mean_exp_case = NA_real_, mean_exp_ctrl = NA_real_,
      beta = NA_real_, se = NA_real_, z = NA_real_,
      p_value = NA_real_, or = NA_real_,
      or_lo95 = NA_real_, or_hi95 = NA_real_,
      fit_status = "no_rows_complete"
    ))
  }

  n_case <- sum(d$case == 1L)
  n_ctrl <- sum(d$case == 0L)
  carrier_case <- sum(d$case == 1L & d[[exposure_col]] >= 1)
  carrier_ctrl <- sum(d$case == 0L & d[[exposure_col]] >= 1)
  mean_exp_case <- mean(d[[exposure_col]][d$case == 1L])
  mean_exp_ctrl <- mean(d[[exposure_col]][d$case == 0L])

  if (length(unique(d[[exposure_col]])) < 2) {
    return(tibble(
      n_case = n_case, n_ctrl = n_ctrl,
      carrier_case = carrier_case, carrier_ctrl = carrier_ctrl,
      mean_exp_case = mean_exp_case, mean_exp_ctrl = mean_exp_ctrl,
      beta = NA_real_, se = NA_real_, z = NA_real_,
      p_value = NA_real_, or = NA_real_,
      or_lo95 = NA_real_, or_hi95 = NA_real_,
      fit_status = "exposure_no_variance"
    ))
  }
  if ((carrier_case + carrier_ctrl) < min_cell) {
    return(tibble(
      n_case = n_case, n_ctrl = n_ctrl,
      carrier_case = carrier_case, carrier_ctrl = carrier_ctrl,
      mean_exp_case = mean_exp_case, mean_exp_ctrl = mean_exp_ctrl,
      beta = NA_real_, se = NA_real_, z = NA_real_,
      p_value = NA_real_, or = NA_real_,
      or_lo95 = NA_real_, or_hi95 = NA_real_,
      fit_status = "insufficient_carriers"
    ))
  }

  rhs <- paste(c(sprintf("`%s`", exposure_col), covs), collapse = " + ")
  form <- as.formula(sprintf("case ~ %s", rhs))

  fit <- tryCatch(
    suppressWarnings(glm(form, data = d, family = binomial(link = "logit"))),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(tibble(
      n_case = n_case, n_ctrl = n_ctrl,
      carrier_case = carrier_case, carrier_ctrl = carrier_ctrl,
      mean_exp_case = mean_exp_case, mean_exp_ctrl = mean_exp_ctrl,
      beta = NA_real_, se = NA_real_, z = NA_real_,
      p_value = NA_real_, or = NA_real_,
      or_lo95 = NA_real_, or_hi95 = NA_real_,
      fit_status = sprintf("glm_error: %s", conditionMessage(fit))
    ))
  }

  coefs <- summary(fit)$coefficients
  target_name <- sprintf("`%s`", exposure_col)
  row_idx <- which(rownames(coefs) %in% c(target_name, exposure_col))
  if (length(row_idx) == 0) {
    return(tibble(
      n_case = n_case, n_ctrl = n_ctrl,
      carrier_case = carrier_case, carrier_ctrl = carrier_ctrl,
      mean_exp_case = mean_exp_case, mean_exp_ctrl = mean_exp_ctrl,
      beta = NA_real_, se = NA_real_, z = NA_real_,
      p_value = NA_real_, or = NA_real_,
      or_lo95 = NA_real_, or_hi95 = NA_real_,
      fit_status = "no_exposure_term"
    ))
  }
  row <- coefs[row_idx[1], , drop = TRUE]
  beta <- unname(row["Estimate"])
  se   <- unname(row["Std. Error"])
  z    <- unname(row["z value"])
  p    <- unname(row["Pr(>|z|)"])
  or   <- exp(beta)
  or_lo <- exp(beta - 1.959964 * se)
  or_hi <- exp(beta + 1.959964 * se)

  status <- "ok"
  if (!fit$converged) status <- "not_converged"
  if (any(is.na(row))) status <- paste(status, "na_in_coef", sep = ";")

  tibble(
    n_case = n_case, n_ctrl = n_ctrl,
    carrier_case = carrier_case, carrier_ctrl = carrier_ctrl,
    mean_exp_case = mean_exp_case, mean_exp_ctrl = mean_exp_ctrl,
    beta = beta, se = se, z = z,
    p_value = p, or = or,
    or_lo95 = or_lo, or_hi95 = or_hi,
    fit_status = status
  )
}

# -------------------------------------------------------------
# 解析 (i): 10 L2 classes
# -------------------------------------------------------------
log_msg("---- Analysis (i): 10 L2 classes -------------------")

l2_rows <- list()
for (l2 in L2_CLASSES) {
  for (sv in SV_TYPES) {
    for (exp_kind in EXPOSURES) {
      col <- sprintf("%s_%s_%s", exp_kind, l2, sv)
      for (cmp in COMPARISONS) {
        sub <- make_subset(bur, cmp)
        res <- fit_logit(sub, col, sv, min_cell = min_cell)
        l2_rows[[length(l2_rows) + 1]] <- tibble(
          analysis       = "L2_class",
          L2_class       = l2,
          sv_type        = sv,
          exposure       = exp_kind,
          exposure_col   = col,
          comparison     = cmp$name,
          case_label     = cmp$case,
          control_label  = cmp$control,
          !!!res
        )
      }
    }
  }
  log_msg(sprintf("  [L2] done: %s", l2))
}
l2_tbl <- bind_rows(l2_rows)

l2_tbl <- l2_tbl %>%
  group_by(sv_type, exposure, comparison) %>%
  mutate(p_fdr_within_stratum = p.adjust(p_value, method = "BH")) %>%
  ungroup()

out_l2 <- file.path(outdir,
                    sprintf("B_prime_L2_classes_results_%s.tsv", V6_SUFFIX))
log_msg(sprintf("Writing: %s", out_l2))
write_tsv(l2_tbl, out_l2)

# -------------------------------------------------------------
# 解析 (ii): 5 group definitions
# -------------------------------------------------------------
log_msg("---- Analysis (ii): specificity groups -------------")

grp_rows <- list()
for (gcol in GROUP_COLUMNS) {
  pattern <- sprintf("^carrier_boundary_%s__", gcol)
  lbl_cols <- grep(pattern, colnames(bur), value = TRUE)
  labels <- lbl_cols %>%
    str_replace(pattern, "") %>%
    str_replace("_(DEL|DUP)$", "") %>%
    unique()
  log_msg(sprintf("  [group] %s -> labels: %s",
                  gcol, paste(labels, collapse = ", ")))

  for (lbl in labels) {
    for (sv in SV_TYPES) {
      for (exp_kind in EXPOSURES) {
        col <- sprintf("%s_%s__%s_%s", exp_kind, gcol, lbl, sv)
        if (!col %in% colnames(bur)) next
        for (cmp in COMPARISONS) {
          sub <- make_subset(bur, cmp)
          res <- fit_logit(sub, col, sv, min_cell = min_cell)
          grp_rows[[length(grp_rows) + 1]] <- tibble(
            analysis       = "specificity_group",
            group_scheme   = gcol,
            group_label    = lbl,
            sv_type        = sv,
            exposure       = exp_kind,
            exposure_col   = col,
            comparison     = cmp$name,
            case_label     = cmp$case,
            control_label  = cmp$control,
            !!!res
          )
        }
      }
    }
  }
}
grp_tbl <- bind_rows(grp_rows)

grp_tbl <- grp_tbl %>%
  group_by(group_scheme, sv_type, exposure, comparison) %>%
  mutate(p_fdr_within_stratum = p.adjust(p_value, method = "BH")) %>%
  ungroup()

out_grp <- file.path(outdir,
                     sprintf("B_prime_specificity_groups_results_%s.tsv", V6_SUFFIX))
log_msg(sprintf("Writing: %s", out_grp))
write_tsv(grp_tbl, out_grp)

# -------------------------------------------------------------
# 実行サマリ
# -------------------------------------------------------------
status_summary <- bind_rows(
  l2_tbl  %>% mutate(block = "L2_class") %>% select(block, fit_status),
  grp_tbl %>% mutate(block = "specificity_group") %>% select(block, fit_status)
) %>%
  count(block, fit_status)

out_sum <- file.path(outdir,
                     sprintf("B_prime_run_summary_%s.tsv", V6_SUFFIX))
log_msg(sprintf("Writing: %s", out_sum))
write_tsv(status_summary, out_sum)

# -------------------------------------------------------------
# preview
# -------------------------------------------------------------
log_msg("Top hits per block (by p_value):")

log_msg("  -- L2 classes, carrier, (ASD/SZ vs HC), DEL --")
print(
  l2_tbl %>%
    filter(exposure == "carrier_boundary",
           comparison %in% c("ASD_vs_HC", "SZ_vs_HC"),
           sv_type == "DEL") %>%
    arrange(p_value) %>%
    head(12) %>%
    select(L2_class, comparison, n_case, n_ctrl,
           carrier_case, carrier_ctrl,
           or, or_lo95, or_hi95, p_value, fit_status)
)

log_msg("  -- group_primary, carrier, ASD_vs_HC | SZ_vs_HC, DEL --")
print(
  grp_tbl %>%
    filter(group_scheme == "group_primary",
           exposure == "carrier_boundary",
           comparison %in% c("ASD_vs_HC", "SZ_vs_HC"),
           sv_type == "DEL") %>%
    arrange(comparison, group_label) %>%
    select(group_label, comparison, n_case, n_ctrl,
           carrier_case, carrier_ctrl,
           or, or_lo95, or_hi95, p_value, fit_status)
)

elapsed <- difftime(Sys.time(), t0, units = "secs")
log_msg(sprintf("Done. PATTERN=%s. Elapsed: %.2f s = %.2f min",
                PATTERN,
                as.numeric(elapsed), as.numeric(elapsed) / 60))
log_msg(strrep("=", 60))
