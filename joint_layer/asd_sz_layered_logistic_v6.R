# asd_sz_layered_logistic_v5.R
#
# 処理概要 (v5 + critical bug fix):
# - ASD + SZ unified joint-layer logistic regression
# - 6 layer (TAD/TRE/PRS/PTV/CNV/NAHR-CNV) progressive nested model
# - Added-last LRT (conditional contribution test)
# - Multinomial logistic for formal ASD-vs-SZ TAD coefficient contrast
# - Static_DEL sensitivity, sparse layer warnings, full diagnostics
#
# v5 → v6 の critical fix:
# - read_tsv の guess_max を Inf に変更 (Matched_SZ column の class 推定 bug 修正)
# - v5 では default guess_max=1000 のため、Matched_SZ column の最初の non-empty 値が
#   row 1153 (file の最後の方) であり、最初 1000 行が全 NA → column class が logical
#   と推定され、後続の non-empty 値も NA に変換されていた
# - 結果: SZ × CNV layer が 0 carriers と誤計算 (実際は 39 carriers)
# - 修正: cnv_del と他の TSV 読み込みで guess_max = Inf を指定
#
# v4 → v5 の変更 (継承):
# - NIG R 環境に nnet パッケージが無いため、multinomial logistic regression を
#   base R + optim() で manual 実装 (依存パッケージ削除)
# - manual 実装は 3-level multinomial (HC=0 / ASD=1 / SZ=2) を MLE で fit
# - Hessian の inverse から vcov を取得し、β_TAD,ASD - β_TAD,SZ の Wald test を実装
#
# 入出力 path は v4 と同一。出力 file 名のみ v4 → v5 に変更。
#
# 出力 (/home/kushima-pg/tre_prs_04282026/):
# - asd_sz_layered_v6.tsv
# - asd_sz_layered_v6.added_last_lrt.tsv
# - asd_sz_layered_v6.multinomial.tsv
# - asd_sz_layered_v6.static_sensitivity.tsv
# - asd_sz_layered_v6.carrier_counts.tsv
# - asd_sz_layered_v6.cor_matrix.tsv
# - asd_sz_layered_v6.diagnostics.tsv
# - asd_sz_layered_v6.summary.txt
# - asd_sz_layered_v6.cohort.<disorder>.<caller>.tsv

start_time <- Sys.time()
options(warn = 1)
cat(sprintf("[%s] Start asd_sz_layered_logistic_v5.R\n", format(start_time)))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  # nnet 削除: optim() で manual multinomial 実装
})

# ---- Paths (v4 と同一) ----
TAD_PATH         <- "/home/kushima-pg/tad04212026/05_wgs_sample_burden/output_v3/sample_burden_L2_and_specificity_v3.tsv"
EHDN_PATH        <- "/home/kushima-pg/str_12282025/analysis_results_novel/outlier_burden_rare_crossfit_v19.per_sample.tsv"
STRLING_PATH     <- "/home/kushima-pg/str_12282025/analysis_results_strling/strling_outlier_burden_rare_crossfit_inbounds_v9.per_sample.tsv"
ASD_PRS_PATH     <- "/lustre12/home/kushima-pg/PRS/GRIFIN-NCBN_ASD_10082025/scores/merged_data_nodup.QC.rs.autosomes_PRS_asd_noMHC.sscore"
SZ_PRS_PATH      <- "/lustre12/home/kushima-pg/PRS/GRIFIN-NCBN_SCZ_12062025/scores/merged_data_nodup.QC.rs.autosomes_PRS_noMHC.sscore"
PTV_CASES        <- "/lustre12/home/kushima-pg/all_riskVariant_01012026/results_v17/rare_ptv_cases_v17.loeuf0.6.txt"
PTV_HEALTHY      <- "/lustre12/home/kushima-pg/all_riskVariant_01012026/results_v17/rare_ptv_healthy_v17.loeuf0.6.txt"
GENESET_PATH     <- "/lustre12/home/kushima-pg/candidate_genes/candidateGenes_asd_sz_sfari_ndd_12312025.txt"
CNV_DEL_PATH     <- "/lustre12/home/kushima-pg/cnv_01012026/analysis/annotsv_results/annotsv_exonic_riskGenes_DEL_v19.tsv"
NAHR_PATH        <- "/lustre12/home/kushima-pg/cnv_01012026/analysis/annotsv_results/annotsv_genomicDisorder_NAHR_v8.tsv"
CURATED_GD       <- "/lustre12/home/kushima-pg/cnv_01012026/curated_genomic_disorder_cnv_loci_v3.txt"
OUT_DIR          <- "/home/kushima-pg/tre_prs_04282026"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Helpers ----
auc_base <- function(y, p) {
  ranks <- rank(p, ties.method = "average")
  n_pos <- sum(y == 1); n_neg <- sum(y == 0)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  (sum(ranks[y == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}
is_yes <- function(x) tolower(trimws(as.character(x))) %in% c("yes", "y", "true", "t", "1")
fmt_or <- function(or_vec) {
  if (length(or_vec) < 4 || is.na(or_vec[1])) return(c("-", "-"))
  c(sprintf("%.3f (%.3f-%.3f)", or_vec[1], or_vec[2], or_vec[3]),
    sprintf("%.3e", or_vec[4]))
}
extract_or <- function(fit, coef_name) {
  cm <- coef(summary(fit))
  if (!coef_name %in% rownames(cm)) return(c(NA, NA, NA, NA))
  beta <- cm[coef_name, "Estimate"]
  se   <- cm[coef_name, "Std. Error"]
  p    <- cm[coef_name, "Pr(>|z|)"]
  c(exp(beta), exp(beta - 1.96 * se), exp(beta + 1.96 * se), p)
}
assert_unique_id <- function(df, id_col, label) {
  dup <- df[[id_col]][duplicated(df[[id_col]])]
  if (length(dup) > 0) {
    stop(sprintf("[%s] duplicated %s: %d rows; examples: %s",
                 label, id_col, length(dup), paste(head(unique(dup), 10), collapse = ",")))
  }
}

# ---- (v5 NEW) Manual multinomial logistic regression via optim() ----
# 3-level outcome: y ∈ {0, 1, 2} where 0 = baseline (HC), 1 = ASD, 2 = SZ
# Likelihood:
#   P(y=0) = 1/D, P(y=1) = exp(X β_1)/D, P(y=2) = exp(X β_2)/D
#   D = 1 + exp(X β_1) + exp(X β_2)
# Parameters: β_1 (length p), β_2 (length p) where p = ncol(X)
multinom_manual <- function(y, X, par_init = NULL, maxit = 5000) {
  p <- ncol(X)
  n <- nrow(X)
  if (is.null(par_init)) par_init <- rep(0, 2 * p)

  neg_loglik <- function(par) {
    beta1 <- par[1:p]
    beta2 <- par[(p + 1):(2 * p)]
    eta1 <- as.vector(X %*% beta1)
    eta2 <- as.vector(X %*% beta2)
    # numerically stable log-sum-exp for log(1 + exp(eta1) + exp(eta2))
    m <- pmax(0, eta1, eta2)
    log_D <- m + log(exp(0 - m) + exp(eta1 - m) + exp(eta2 - m))
    ll <- sum((y == 1) * eta1 + (y == 2) * eta2 - log_D)
    -ll  # optim minimizes
  }

  # Gradient (analytical for speed and stability)
  grad_neg_loglik <- function(par) {
    beta1 <- par[1:p]
    beta2 <- par[(p + 1):(2 * p)]
    eta1 <- as.vector(X %*% beta1)
    eta2 <- as.vector(X %*% beta2)
    m <- pmax(0, eta1, eta2)
    e0 <- exp(0 - m); e1 <- exp(eta1 - m); e2 <- exp(eta2 - m)
    D  <- e0 + e1 + e2
    P1 <- e1 / D
    P2 <- e2 / D
    grad1 <- crossprod(X, P1 - (y == 1))
    grad2 <- crossprod(X, P2 - (y == 2))
    as.vector(c(grad1, grad2))
  }

  fit <- optim(par_init, fn = neg_loglik, gr = grad_neg_loglik,
               method = "BFGS", hessian = TRUE,
               control = list(maxit = maxit, reltol = 1e-10))

  vcov_mat <- tryCatch(solve(fit$hessian), error = function(e) NULL)

  list(
    par = fit$par,
    coef_ASD = fit$par[1:p],
    coef_SZ  = fit$par[(p + 1):(2 * p)],
    var_names = colnames(X),
    vcov = vcov_mat,
    convergence = fit$convergence,
    loglik = -fit$value,
    n = n, p = p
  )
}

# Wald test on β_var,ASD - β_var,SZ
multinom_contrast <- function(fit, var_name) {
  k <- which(fit$var_names == var_name)
  if (length(k) == 0) return(list(diff_beta = NA, SE_diff = NA, z = NA, p_two = NA, p_one = NA,
                                   beta_asd = NA, beta_sz = NA,
                                   se_asd = NA, se_sz = NA))

  p <- fit$p
  beta_asd <- fit$coef_ASD[k]
  beta_sz  <- fit$coef_SZ[k]

  if (is.null(fit$vcov)) {
    return(list(diff_beta = beta_asd - beta_sz, SE_diff = NA, z = NA, p_two = NA, p_one = NA,
                beta_asd = beta_asd, beta_sz = beta_sz, se_asd = NA, se_sz = NA))
  }

  var_asd <- fit$vcov[k, k]
  var_sz  <- fit$vcov[p + k, p + k]
  cov_as  <- fit$vcov[k, p + k]
  se_asd  <- sqrt(var_asd)
  se_sz   <- sqrt(var_sz)
  var_diff <- var_asd + var_sz - 2 * cov_as
  if (var_diff < 0) {
    return(list(diff_beta = beta_asd - beta_sz, SE_diff = NA, z = NA, p_two = NA, p_one = NA,
                beta_asd = beta_asd, beta_sz = beta_sz, se_asd = se_asd, se_sz = se_sz))
  }
  se_diff <- sqrt(var_diff)
  diff_beta <- beta_asd - beta_sz
  z <- diff_beta / se_diff
  p_two <- 2 * pnorm(-abs(z))
  p_one <- pnorm(-z)
  list(diff_beta = diff_beta, SE_diff = se_diff, z = z,
       p_two = p_two, p_one = p_one,
       beta_asd = beta_asd, beta_sz = beta_sz,
       se_asd = se_asd, se_sz = se_sz)
}

# ---- 1. Load shared data ----
cat("\n=== 1. Loading shared data sources ===\n")
gene_sets <- read.table(GENESET_PATH, header = TRUE, stringsAsFactors = FALSE)
ptv_cases <- read_tsv(PTV_CASES, show_col_types = FALSE, name_repair = "minimal") %>%
  select(SampleID = Pid, SYMBOL)
ptv_hc    <- read_tsv(PTV_HEALTHY, show_col_types = FALSE, name_repair = "minimal") %>%
  select(SampleID = Pid, SYMBOL)
ptv_all   <- bind_rows(ptv_cases, ptv_hc)
cnv_del <- read_tsv(CNV_DEL_PATH, show_col_types = FALSE, name_repair = "minimal", guess_max = Inf)
gd_curated <- read_tsv(CURATED_GD, show_col_types = FALSE)
nahr_per_sv <- read_tsv(NAHR_PATH, show_col_types = FALSE)

tad_raw <- read_tsv(TAD_PATH, show_col_types = FALSE)
required_tad_cols <- c("sample_id", "Diagnosis", "Sex_numeric", paste0("PC", 1:10),
                       "log1p_total_del_bases", "log1p_total_gene_DEL",
                       "n_boundary_group_primary__Diff_specific_n1_DEL",
                       "n_boundary_group_primary__Diff_shared_n2plus_DEL",
                       "n_boundary_group_primary__Static_DEL")
miss <- setdiff(required_tad_cols, names(tad_raw))
if (length(miss) > 0) stop(sprintf("TAD missing: %s", paste(miss, collapse = ", ")))

tad <- tad_raw %>%
  mutate(SampleID = sample_id,
         Diff_any_DEL = `n_boundary_group_primary__Diff_specific_n1_DEL` +
                        `n_boundary_group_primary__Diff_shared_n2plus_DEL`,
         Static_DEL   = `n_boundary_group_primary__Static_DEL`) %>%
  select(SampleID, Diagnosis, Sex_numeric,
         PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
         log1p_total_del_bases, log1p_total_gene_DEL, Diff_any_DEL, Static_DEL)
cat(sprintf("TAD selected: %d rows; group counts: ", nrow(tad))); print(table(tad$Diagnosis))

ehdn    <- read_tsv(EHDN_PATH, show_col_types = FALSE) %>%
  select(SampleID, Group, rare_any_EHdn = rare_any, Depth)
strling <- read_tsv(STRLING_PATH, show_col_types = FALSE) %>%
  select(SampleID, rare_any_STRling = rare_any)

assert_unique_id(tad,     "SampleID", "TAD")
assert_unique_id(ehdn,    "SampleID", "EHdn")
assert_unique_id(strling, "SampleID", "STRling")

build_disorder_layers <- function(disorder, ptv_geneset_col, cnv_matched_col, nahr_assoc_col) {
  cat(sprintf("\n=== Building layers for %s ===\n", disorder))
  geneset_yes <- unique(gene_sets$gene[is_yes(gene_sets[[ptv_geneset_col]])])
  ptv_carriers <- ptv_all %>% filter(SYMBOL %in% geneset_yes) %>% pull(SampleID) %>% unique()
  cnv_carriers <- cnv_del %>%
    filter(!is.na(.data[[cnv_matched_col]]) & trimws(.data[[cnv_matched_col]]) != "") %>%
    pull(Samples_ID) %>% unique()
  nahr_loci <- gd_curated %>% filter(is_yes(nahr), is_yes(.data[[nahr_assoc_col]])) %>% pull(gd_id)
  nahr_carriers <- nahr_per_sv %>%
    rowwise() %>%
    mutate(has_assoc = any(trimws(unlist(strsplit(as.character(gd_id), ","))) %in% nahr_loci)) %>%
    ungroup() %>% filter(has_assoc) %>% pull(Samples_ID) %>% unique()
  cat(sprintf("  PTV: %d genes, %d carriers; CNV: %d carriers; NAHR: %d loci, %d carriers\n",
              length(geneset_yes), length(ptv_carriers), length(cnv_carriers),
              length(nahr_loci), length(nahr_carriers)))
  list(ptv = ptv_carriers, cnv = cnv_carriers, nahr = nahr_carriers,
       n_geneset_genes = length(geneset_yes), n_nahr_loci = length(nahr_loci))
}

asd_layers <- build_disorder_layers("ASD", "asd_fu_fdr005", "Matched_ASD", "asd_associated")
sz_layers  <- build_disorder_layers("SZ",  "scz_chick_fdr005","Matched_SZ", "sz_associated")

load_prs <- function(prs_path, label) {
  prs_raw <- read_tsv(prs_path, show_col_types = FALSE, comment = "")
  names(prs_raw) <- gsub("^#", "", names(prs_raw))
  out <- prs_raw %>% select(SampleID = IID, PRS_raw = SCORE1_AVG)
  assert_unique_id(out, "SampleID", label)
  out
}
asd_prs <- load_prs(ASD_PRS_PATH, "ASD_PRS")
sz_prs  <- load_prs(SZ_PRS_PATH,  "SZ_PRS")

build_joint <- function(disorder, caller_name, prs_df, layers) {
  cat(sprintf("\n--- %s × %s ---\n", disorder, caller_name))
  df <- tad %>% inner_join(ehdn, by = "SampleID") %>% rename(Group_TRE = Group)
  if (caller_name == "EHdn") df$rare_any <- df$rare_any_EHdn
  else if (caller_name == "STRling") {
    df <- df %>% inner_join(strling, by = "SampleID")
    df$rare_any <- df$rare_any_STRling
  }
  df <- df %>% inner_join(prs_df, by = "SampleID") %>%
    mutate(PTV  = as.integer(SampleID %in% layers$ptv),
           CNV  = as.integer(SampleID %in% layers$cnv),
           NAHR = as.integer(SampleID %in% layers$nahr)) %>%
    filter(Diagnosis %in% c(disorder, "Healthy"))
  df$outcome <- as.integer(df$Diagnosis == disorder)
  hc_mean <- mean(df$PRS_raw[df$outcome == 0], na.rm = TRUE)
  hc_sd   <- sd(df$PRS_raw[df$outcome == 0], na.rm = TRUE)
  df$PRS_Z <- (df$PRS_raw - hc_mean) / hc_sd
  required_cols <- c("Sex_numeric", "Depth", paste0("PC", 1:10),
                     "log1p_total_del_bases", "log1p_total_gene_DEL",
                     "Diff_any_DEL", "Static_DEL", "rare_any", "PRS_Z",
                     "PTV", "CNV", "NAHR", "outcome")
  df <- df %>% drop_na(all_of(required_cols))
  cat(sprintf("  N=%d (case=%d, HC=%d)\n", nrow(df), sum(df$outcome == 1), sum(df$outcome == 0)))
  write_tsv(df %>% select(SampleID, Diagnosis, outcome,
                          Diff_any_DEL, Static_DEL, rare_any, PRS_Z, PRS_raw,
                          PTV, CNV, NAHR, Sex_numeric, Depth,
                          all_of(paste0("PC", 1:10)),
                          log1p_total_del_bases, log1p_total_gene_DEL),
            file.path(OUT_DIR, sprintf("asd_sz_layered_v6.cohort.%s.%s.tsv", disorder, caller_name)))
  df
}

detect_sparse_warning <- function(df, var_name) {
  vec <- df[[var_name]]
  if (length(unique(vec)) <= 1) return("constant_predictor")
  if (is.numeric(vec) && all(vec %in% c(0, 1))) {
    cc <- sum(vec == 1 & df$outcome == 1)
    co <- sum(vec == 1 & df$outcome == 0)
    if (cc == 0) return("no_case_carrier")
    if (co == 0) return("no_control_carrier")
    if (cc < 3 || co < 3) return("low_carrier_count")
  }
  ""
}

glm_diagnostics <- function(fit, disorder, caller, model_name) {
  beta_vec <- coef(fit)
  beta_no_int <- beta_vec[setdiff(names(beta_vec), "(Intercept)")]
  max_abs_beta <- if (length(beta_no_int) > 0) max(abs(beta_no_int), na.rm = TRUE) else NA
  data.frame(
    Disorder = disorder, Caller = caller, Model = model_name,
    Converged = fit$converged, N_iter = fit$iter,
    Max_abs_beta = round(max_abs_beta, 4),
    Aliased_coef_count = sum(is.na(beta_vec)),
    Deviance = round(fit$deviance, 3),
    Null_deviance = round(fit$null.deviance, 3),
    DF_residual = fit$df.residual,
    stringsAsFactors = FALSE
  )
}

run_panel <- function(disorder, caller_name, prs_df, layers) {
  df <- build_joint(disorder, caller_name, prs_df, layers)
  cat(sprintf("\n========== %s × %s panel ==========\n", disorder, caller_name))

  base_formula <- outcome ~ Sex_numeric + Depth +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
    log1p_total_del_bases + log1p_total_gene_DEL

  fits <- list(
    M0    = glm(base_formula, family = binomial(), data = df),
    M1    = glm(update(base_formula, . ~ . + Diff_any_DEL),
                family = binomial(), data = df),
    M2    = glm(update(base_formula, . ~ . + Diff_any_DEL + rare_any),
                family = binomial(), data = df),
    M3    = glm(update(base_formula, . ~ . + Diff_any_DEL + rare_any + PRS_Z),
                family = binomial(), data = df),
    M4    = glm(update(base_formula, . ~ . + Diff_any_DEL + rare_any + PRS_Z + PTV),
                family = binomial(), data = df),
    M5    = glm(update(base_formula, . ~ . + Diff_any_DEL + rare_any + PRS_Z + PTV + CNV),
                family = binomial(), data = df),
    Mfull = glm(update(base_formula, . ~ . + Diff_any_DEL + rare_any + PRS_Z + PTV + CNV + NAHR),
                family = binomial(), data = df)
  )

  fit_static_sens <- tryCatch(
    glm(update(base_formula, . ~ . + Diff_any_DEL + Static_DEL + rare_any + PRS_Z + PTV + CNV + NAHR),
        family = binomial(), data = df),
    error = function(e) NULL)

  diag_rows <- list()
  for (m in names(fits)) diag_rows[[m]] <- glm_diagnostics(fits[[m]], disorder, caller_name, m)
  if (!is.null(fit_static_sens)) {
    diag_rows[["Mfull_static"]] <- glm_diagnostics(fit_static_sens, disorder, caller_name, "Mfull_static")
  }
  diag_df <- do.call(rbind, diag_rows)

  layer_keys <- c(TAD = "Diff_any_DEL", TRE = "rare_any", PRS = "PRS_Z",
                  PTV = "PTV", CNV = "CNV", NAHR = "NAHR")
  added_last_lrt <- list()
  for (lk in names(layer_keys)) {
    target_var <- layer_keys[[lk]]
    full_layers <- c("Diff_any_DEL", "rare_any", "PRS_Z", "PTV", "CNV", "NAHR")
    kept <- setdiff(full_layers, target_var)
    f_minus <- as.formula(paste("~ . +", paste(kept, collapse = " + ")))
    sparse_note <- detect_sparse_warning(df, target_var)
    fit_minus <- tryCatch(glm(update(base_formula, f_minus), family = binomial(), data = df),
                          error = function(e) NULL)
    full_aliased <- sum(is.na(coef(fits$Mfull))) > 0
    full_max_abs <- max(abs(coef(fits$Mfull)[setdiff(names(coef(fits$Mfull)), "(Intercept)")]),
                        na.rm = TRUE)
    if (full_aliased) sparse_note <- paste(sparse_note, "aliased_in_full", sep = ";")
    if (full_max_abs > 10) sparse_note <- paste(sparse_note, "possible_separation", sep = ";")
    sparse_note <- gsub("^;|;$", "", sparse_note)

    if (is.null(fit_minus) || !fit_minus$converged) {
      added_last_lrt[[lk]] <- list(p = NA, df_diff = NA, llr_diff = NA,
                                    note = paste(c(sparse_note, "fit_failed_or_unconverged"), collapse = ";"))
    } else {
      a <- anova(fit_minus, fits$Mfull, test = "LRT")
      added_last_lrt[[lk]] <- list(p = a$`Pr(>Chi)`[2], df_diff = a$Df[2],
                                    llr_diff = a$Deviance[2], note = sparse_note)
    }
  }

  fit_null <- glm(outcome ~ 1, family = binomial(), data = df)
  nagelkerke_r2 <- function(fit) {
    n <- nobs(fit); ll_full <- as.numeric(logLik(fit)); ll_null <- as.numeric(logLik(fit_null))
    r2_cs <- 1 - exp((2/n)*(ll_null - ll_full)); r2_max <- 1 - exp((2/n)*ll_null)
    r2_cs / r2_max
  }
  auc_compute <- function(fit) auc_base(df$outcome, predict(fit, type = "response"))

  results <- list()
  for (m in names(fits)) {
    fit <- fits[[m]]
    res_m <- list()
    for (lk in names(layer_keys)) res_m[[lk]] <- extract_or(fit, layer_keys[[lk]])
    res_m$r2 <- nagelkerke_r2(fit); res_m$auc <- auc_compute(fit)
    results[[m]] <- res_m
  }

  static_sens_result <- NULL
  if (!is.null(fit_static_sens) && fit_static_sens$converged) {
    static_sens_result <- list(
      Diff_any = extract_or(fit_static_sens, "Diff_any_DEL"),
      Static = extract_or(fit_static_sens, "Static_DEL"),
      TRE = extract_or(fit_static_sens, "rare_any"),
      r2 = nagelkerke_r2(fit_static_sens),
      auc = auc_compute(fit_static_sens)
    )
  }

  lrt_seq <- list(
    M1_vs_M0=anova(fits$M0,fits$M1,test="LRT")$`Pr(>Chi)`[2],
    M2_vs_M1=anova(fits$M1,fits$M2,test="LRT")$`Pr(>Chi)`[2],
    M3_vs_M2=anova(fits$M2,fits$M3,test="LRT")$`Pr(>Chi)`[2],
    M4_vs_M3=anova(fits$M3,fits$M4,test="LRT")$`Pr(>Chi)`[2],
    M5_vs_M4=anova(fits$M4,fits$M5,test="LRT")$`Pr(>Chi)`[2],
    Mfull_vs_M5=anova(fits$M5,fits$Mfull,test="LRT")$`Pr(>Chi)`[2]
  )

  ctrl <- df %>% filter(outcome == 0)
  layer_vars <- c("Diff_any_DEL", "Static_DEL", "rare_any", "PRS_Z", "PTV", "CNV", "NAHR")
  cor_mat <- tryCatch(cor(ctrl[, layer_vars], use = "pairwise.complete.obs", method = "pearson"),
                      error = function(e) NULL)

  iqr_safe <- function(x) {
    if (length(unique(x)) <= 1) return(c(NA, NA, NA, NA))
    q <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
    c(mean(x, na.rm = TRUE), q[1], q[2], q[3])
  }
  cc_rows <- list()
  N_case <- sum(df$outcome == 1); N_ctrl <- sum(df$outcome == 0)
  for (lk in c("Diff_any_DEL", "Static_DEL", "PRS_Z")) {
    ca <- iqr_safe(df[[lk]][df$outcome == 1]); co <- iqr_safe(df[[lk]][df$outcome == 0])
    cc_rows[[lk]] <- data.frame(
      Disorder=disorder, Caller=caller_name, Layer=lk, Metric_type="continuous_mean",
      N_case=N_case, N_ctrl=N_ctrl, Case_count=NA, Ctrl_count=NA,
      Case_rate=NA, Ctrl_rate=NA,
      Case_mean=round(ca[1],4), Ctrl_mean=round(co[1],4),
      Case_median=round(ca[3],4), Ctrl_median=round(co[3],4),
      Case_Q1=round(ca[2],4), Case_Q3=round(ca[4],4),
      Ctrl_Q1=round(co[2],4), Ctrl_Q3=round(co[4],4),
      stringsAsFactors=FALSE)
  }
  for (lk in c("rare_any", "PTV", "CNV", "NAHR")) {
    cc_count <- sum(df[[lk]] == 1 & df$outcome == 1)
    cn_count <- sum(df[[lk]] == 1 & df$outcome == 0)
    cc_rows[[lk]] <- data.frame(
      Disorder=disorder, Caller=caller_name, Layer=lk, Metric_type="binary_count",
      N_case=N_case, N_ctrl=N_ctrl, Case_count=cc_count, Ctrl_count=cn_count,
      Case_rate=round(cc_count/N_case,4), Ctrl_rate=round(cn_count/N_ctrl,4),
      Case_mean=NA, Ctrl_mean=NA, Case_median=NA, Ctrl_median=NA,
      Case_Q1=NA, Case_Q3=NA, Ctrl_Q1=NA, Ctrl_Q3=NA,
      stringsAsFactors=FALSE)
  }
  carrier_counts <- do.call(rbind, cc_rows)

  list(disorder=disorder, caller=caller_name,
       n_total=nrow(df), n_case=sum(df$outcome==1), n_ctrl=sum(df$outcome==0),
       models=results, lrt_seq=lrt_seq, added_last_lrt=added_last_lrt,
       static_sens=static_sens_result,
       cor_mat=cor_mat, carrier_counts=carrier_counts, diagnostics=diag_df)
}

panels <- list(
  asd_ehdn    = run_panel("ASD", "EHdn",    asd_prs, asd_layers),
  asd_strling = run_panel("ASD", "STRling", asd_prs, asd_layers),
  sz_ehdn     = run_panel("SZ",  "EHdn",    sz_prs,  sz_layers),
  sz_strling  = run_panel("SZ",  "STRling", sz_prs,  sz_layers)
)

# ---- 7. Multinomial logistic (manual via optim) ----
cat("\n=== 7. Manual multinomial logistic (5 rows) ===\n")

build_multi_indep <- function() {
  df <- tad %>% filter(Diagnosis %in% c("Healthy", "ASD", "SZ"))
  cols <- c("Sex_numeric", paste0("PC", 1:10),
            "log1p_total_del_bases", "log1p_total_gene_DEL", "Diff_any_DEL")
  df <- df %>% drop_na(all_of(cols))
  df$y <- ifelse(df$Diagnosis == "Healthy", 0,
                 ifelse(df$Diagnosis == "ASD", 1, 2))
  df
}

build_multi_caller <- function(caller_name) {
  df <- tad %>% inner_join(ehdn, by = "SampleID") %>% rename(Group_TRE = Group)
  if (caller_name == "EHdn") df$rare_any <- df$rare_any_EHdn
  else if (caller_name == "STRling") {
    df <- df %>% inner_join(strling, by = "SampleID")
    df$rare_any <- df$rare_any_STRling
  }
  df <- df %>% filter(Diagnosis %in% c("Healthy", "ASD", "SZ"))
  cols <- c("Sex_numeric", "Depth", paste0("PC", 1:10),
            "log1p_total_del_bases", "log1p_total_gene_DEL", "Diff_any_DEL", "rare_any")
  df <- df %>% drop_na(all_of(cols))
  df$y <- ifelse(df$Diagnosis == "Healthy", 0,
                 ifelse(df$Diagnosis == "ASD", 1, 2))
  df
}

run_multinomial <- function(df_multi, model_label, include_TRE, include_Depth) {
  cat(sprintf("\nMultinomial: %s, N=%d (HC=%d, ASD=%d, SZ=%d)\n",
              model_label, nrow(df_multi),
              sum(df_multi$y == 0), sum(df_multi$y == 1), sum(df_multi$y == 2)))

  formula_terms <- c("Sex_numeric", paste0("PC", 1:10),
                     "log1p_total_del_bases", "log1p_total_gene_DEL", "Diff_any_DEL")
  if (include_Depth) formula_terms <- c(formula_terms, "Depth")
  if (include_TRE)   formula_terms <- c(formula_terms, "rare_any")

  # Build design matrix with intercept
  mf_formula <- as.formula(paste("~", paste(formula_terms, collapse = " + ")))
  X <- model.matrix(mf_formula, data = df_multi)
  y <- df_multi$y

  fit <- multinom_manual(y, X)
  contrast <- multinom_contrast(fit, "Diff_any_DEL")

  data.frame(
    Model_label = model_label,
    Include_TRE = include_TRE, Include_Depth = include_Depth,
    N_total = nrow(df_multi),
    N_HC = sum(y == 0), N_ASD = sum(y == 1), N_SZ = sum(y == 2),
    Convergence = fit$convergence,
    LogLik = sprintf("%.3f", fit$loglik),
    OR_TAD_ASD = sprintf("%.3f (%.3f-%.3f)",
                          exp(contrast$beta_asd),
                          exp(contrast$beta_asd - 1.96 * contrast$se_asd),
                          exp(contrast$beta_asd + 1.96 * contrast$se_asd)),
    OR_TAD_SZ  = sprintf("%.3f (%.3f-%.3f)",
                          exp(contrast$beta_sz),
                          exp(contrast$beta_sz  - 1.96 * contrast$se_sz),
                          exp(contrast$beta_sz  + 1.96 * contrast$se_sz)),
    beta_TAD_ASD = sprintf("%.3f", contrast$beta_asd),
    beta_TAD_SZ  = sprintf("%.3f", contrast$beta_sz),
    diff_beta = ifelse(is.na(contrast$diff_beta), "-", sprintf("%.3f", contrast$diff_beta)),
    SE_diff = ifelse(is.na(contrast$SE_diff), "-", sprintf("%.3f", contrast$SE_diff)),
    z_diff = ifelse(is.na(contrast$z), "-", sprintf("%.3f", contrast$z)),
    P_two_sided = ifelse(is.na(contrast$p_two), "-", sprintf("%.3e", contrast$p_two)),
    P_one_sided_ASD_gt_SZ = ifelse(is.na(contrast$p_one), "-", sprintf("%.3e", contrast$p_one)),
    stringsAsFactors = FALSE
  )
}

df_indep <- build_multi_indep()
df_ehdn  <- build_multi_caller("EHdn")
df_strl  <- build_multi_caller("STRling")

multinomial_panel <- rbind(
  run_multinomial(df_indep, "TAD-only (caller-independent, no Depth)",          FALSE, FALSE),
  run_multinomial(df_ehdn,  "TAD-only (EHdn complete-case subset, +Depth)",     FALSE, TRUE),
  run_multinomial(df_ehdn,  "TAD+TRE (EHdn complete-case subset, +Depth)",      TRUE,  TRUE),
  run_multinomial(df_strl,  "TAD-only (STRling complete-case subset, +Depth)",  FALSE, TRUE),
  run_multinomial(df_strl,  "TAD+TRE (STRling complete-case subset, +Depth)",   TRUE,  TRUE)
)
write_tsv(multinomial_panel, file.path(OUT_DIR, "asd_sz_layered_v6.multinomial.tsv"))

# ---- 8. Format output panels ----
desc_map <- c(M0 = "Base (Sex+Depth+PCs+log1p[del_bases]+log1p[gene_DEL])",
              M1 = "Base + TAD",
              M2 = "Base + TAD + TRE",
              M3 = "Base + TAD + TRE + PRS",
              M4 = "Base + TAD + TRE + PRS + PTV",
              M5 = "Base + TAD + TRE + PRS + PTV + CNV",
              Mfull = "All layers (+ NAHR_CNV)")

format_main_panel <- function(res) {
  rows <- list()
  for (m in names(res$models)) {
    mdl <- res$models[[m]]
    tad_str <- fmt_or(mdl$TAD); tre_str <- fmt_or(mdl$TRE); prs_str <- fmt_or(mdl$PRS)
    ptv_str <- fmt_or(mdl$PTV); cnv_str <- fmt_or(mdl$CNV); nahr_str <- fmt_or(mdl$NAHR)
    rows[[m]] <- data.frame(
      Disorder = res$disorder, Caller = res$caller, Model = m,
      Description = desc_map[[m]],
      N_total = res$n_total, N_case = res$n_case, N_ctrl = res$n_ctrl,
      TAD_OR_95CI = tad_str[1], TAD_P = tad_str[2],
      TRE_OR_95CI = tre_str[1], TRE_P = tre_str[2],
      PRS_OR_95CI = prs_str[1], PRS_P = prs_str[2],
      PTV_OR_95CI = ptv_str[1], PTV_P = ptv_str[2],
      CNV_OR_95CI = cnv_str[1], CNV_P = cnv_str[2],
      NAHR_OR_95CI = nahr_str[1], NAHR_P = nahr_str[2],
      Nagelkerke_R2 = sprintf("%.4f", mdl$r2), AUC = sprintf("%.4f", mdl$auc),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}
main_panel <- do.call(rbind, lapply(panels, format_main_panel))
write_tsv(main_panel, file.path(OUT_DIR, "asd_sz_layered_v6.tsv"))

format_added_last <- function(res) {
  rows <- list()
  for (lk in names(res$added_last_lrt)) {
    al <- res$added_last_lrt[[lk]]
    rows[[lk]] <- data.frame(
      Disorder = res$disorder, Caller = res$caller, Layer = lk,
      DF_diff = ifelse(is.null(al$df_diff)||is.na(al$df_diff), NA, al$df_diff),
      LL_ratio_diff = ifelse(is.null(al$llr_diff)||is.na(al$llr_diff), NA, sprintf("%.3f", al$llr_diff)),
      Conditional_LRT_P = ifelse(is.null(al$p)||is.na(al$p), NA, sprintf("%.3e", al$p)),
      Note = ifelse(is.null(al$note)||al$note=="", "", al$note),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}
added_last_panel <- do.call(rbind, lapply(panels, format_added_last))
write_tsv(added_last_panel, file.path(OUT_DIR, "asd_sz_layered_v6.added_last_lrt.tsv"))

format_static_sens <- function(res) {
  if (is.null(res$static_sens)) {
    return(data.frame(Disorder=res$disorder, Caller=res$caller,
                      Diff_any_OR="-", Diff_any_P="-", Static_OR="-", Static_P="-",
                      TRE_OR="-", TRE_P="-", Nagelkerke_R2="-", AUC="-",
                      Note="fit_failed", stringsAsFactors=FALSE))
  }
  ss <- res$static_sens
  da <- fmt_or(ss$Diff_any); st <- fmt_or(ss$Static); tr <- fmt_or(ss$TRE)
  data.frame(Disorder=res$disorder, Caller=res$caller,
             Diff_any_OR=da[1], Diff_any_P=da[2],
             Static_OR=st[1], Static_P=st[2],
             TRE_OR=tr[1], TRE_P=tr[2],
             Nagelkerke_R2=sprintf("%.4f", ss$r2), AUC=sprintf("%.4f", ss$auc),
             Note="", stringsAsFactors=FALSE)
}
static_sens_panel <- do.call(rbind, lapply(panels, format_static_sens))
write_tsv(static_sens_panel, file.path(OUT_DIR, "asd_sz_layered_v6.static_sensitivity.tsv"))

carrier_count_panel <- do.call(rbind, lapply(panels, function(p) p$carrier_counts))
write_tsv(carrier_count_panel, file.path(OUT_DIR, "asd_sz_layered_v6.carrier_counts.tsv"))

diag_panel <- do.call(rbind, lapply(panels, function(p) p$diagnostics))
write_tsv(diag_panel, file.path(OUT_DIR, "asd_sz_layered_v6.diagnostics.tsv"))

cor_long <- list()
for (panel_id in names(panels)) {
  res <- panels[[panel_id]]
  cm <- res$cor_mat
  if (is.null(cm)) next
  layer_vars <- rownames(cm)
  for (i in seq_along(layer_vars)) for (j in seq_along(layer_vars)) {
    if (j > i) {
      cor_long[[length(cor_long)+1]] <- data.frame(
        Disorder=res$disorder, Caller=res$caller,
        Layer1=layer_vars[i], Layer2=layer_vars[j],
        Pearson_r=sprintf("%.4f", cm[i,j]), stringsAsFactors=FALSE)
    }
  }
}
cor_panel <- do.call(rbind, cor_long)
write_tsv(cor_panel, file.path(OUT_DIR, "asd_sz_layered_v6.cor_matrix.tsv"))

# ---- Text summary ----
sink(file.path(OUT_DIR, "asd_sz_layered_v6.summary.txt"))
cat("=== ASD + SZ Unified Layered Logistic Regression v5 ===\n")
cat(sprintf("Date: %s\n", format(Sys.time())))
cat("\nv4 → v5 change: nnet dependency removed; multinomial logistic implemented\n")
cat("manually via base R + optim() (BFGS with analytical gradient + Hessian)\n")
cat("All other functionality identical to v4.\n")
cat("\n=== Multinomial TAD contrast (5 rows; manual MLE) ===\n")
print(multinomial_panel)
cat("\n=== Static sensitivity ===\n"); print(static_sens_panel)
cat("\n=== Main panel ===\n"); print(main_panel)
cat("\n=== Added-last LRT ===\n"); print(added_last_panel)
cat("\n=== Carrier counts ===\n"); print(carrier_count_panel)
cat("\n=== HC layer correlation ===\n"); print(cor_panel)
cat("\n=== Diagnostics ===\n"); print(diag_panel)
cat("\n=== sessionInfo() ===\n"); print(sessionInfo())
sink()

cat("\n========== MAIN PANEL ==========\n");           print(main_panel)
cat("\n========== MULTINOMIAL ==========\n");          print(multinomial_panel)
cat("\n========== STATIC SENSITIVITY ==========\n");   print(static_sens_panel)
cat("\n========== ADDED-LAST LRT ==========\n");       print(added_last_panel)

end_time <- Sys.time()
elapsed  <- difftime(end_time, start_time, units = "secs")
cat(sprintf("\n[%s] End. Elapsed: %.1f seconds\n", format(end_time), as.numeric(elapsed)))
