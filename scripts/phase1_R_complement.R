#
# Phase 1 R Complement: Pseudobulk DEG + Variance Decomposition
# Requires: R 4.5+ with variancePartition, dream, DESeq2, lme4, lmerTest
# Input: pseudobulk_donor_counts.csv.gz + donor_metadata.csv (from Python step)
#

library(variancePartition)
library(DESeq2)
library(lme4)
library(lmerTest)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)

RES_DIR <- "D:/openclaw_workspace/atherosclerosis_atlas/results"

# -------------------------------------------------------------------
# Load pseudobulk data from Python step
# -------------------------------------------------------------------
message("Loading pseudobulk data...")
pseudobulk <- as.matrix(read.csv(
  file.path(RES_DIR, "pseudobulk_donor_counts.csv.gz"),
  row.names = 1, check.names = FALSE
))
donor_meta <- read.csv(
  file.path(RES_DIR, "donor_metadata.csv"),
  row.names = 1
)

message(sprintf("Pseudobulk: %d donors x %d genes", nrow(pseudobulk), ncol(pseudobulk)))
message(sprintf("Donors per bed: %s", paste(table(donor_meta$plaque_location), collapse = ", ")))

# Align
common_donors <- intersect(rownames(pseudobulk), rownames(donor_meta))
pseudobulk <- pseudobulk[common_donors, ]
donor_meta <- donor_meta[common_donors, ]

# -------------------------------------------------------------------
# 1. Variance Decomposition (P1)
# -------------------------------------------------------------------
message("\n=== Variance Decomposition ===")

# VST normalization for variancePartition
dds <- DESeqDataSetFromMatrix(
  countData = t(pseudobulk),
  colData = donor_meta,
  design = ~ plaque_location + sex
)

# Filter lowly expressed genes
keep <- rowSums(counts(dds) >= 10) >= ncol(dds) * 0.3
dds <- dds[keep, ]
message(sprintf("Genes retained: %d / %d", sum(keep), length(keep)))

vsd <- vst(dds, blind = FALSE)
vst_counts <- assay(vsd)

# variancePartition
form <- ~ (1|plaque_location) + (1|dataset) + (1|sex) + age
vp <- fitExtractVarPartModel(vst_counts, form, donor_meta)

# Plot
pdf(file.path(RES_DIR, "variance_partition.pdf"), width = 10, height = 6)
plotVarPart(vp, label.angle = 45)
dev.off()

# Summary (base R apply workaround for colMedians compatibility)
vp_median <- apply(as.matrix(vp), 2, median)
vp_summary <- data.frame(
  median_variance_explained = vp_median,
  row.names = names(vp_median)
)
vp_summary <- vp_summary[order(vp_summary$median_variance_explained, decreasing = TRUE), , drop = FALSE]
print(vp_summary)

write.csv(vp_summary, file.path(RES_DIR, "variance_partition_summary.csv"))

# -------------------------------------------------------------------
# 2. Differential Expression with dream (P1)
# -------------------------------------------------------------------
message("\n=== Differential Expression (dream) ===")

# Focus on highly variable genes
rv <- rowVars(vst_counts)
top_variable <- order(rv, decreasing = TRUE)[1:2000]
expr_subset <- vst_counts[top_variable, ]

# dream analysis — run each contrast via factor releveling
# Falls back to simpler models if coefficient not estimable
run_dream_contrast <- function(expr_data, meta, ref_level, comp_level, form_de) {
  meta$plaque_location <- relevel(as.factor(meta$plaque_location), ref = ref_level)

  # Align meta to expression data columns
  common_samples <- intersect(colnames(expr_data), rownames(meta))
  expr_aligned <- expr_data[, common_samples, drop = FALSE]
  meta_aligned <- meta[common_samples, , drop = FALSE]

  coef_name <- paste0("plaque_location", comp_level)

  # Strategy 1: dream with full formula
  fit <- tryCatch(
    dream(expr_aligned, form_de, meta_aligned),
    error = function(e) NULL
  )
  if (!is.null(fit) && coef_name %in% colnames(coef(fit))) {
    fit <- eBayes(fit)
    res <- topTable(fit, coef = coef_name, number = Inf, sort.by = "p")
    if (!all(is.na(res$logFC))) {
      return(list(res = res, method = "dream"))
    }
  }

  # Strategy 2: lmFit with full model
  message("  -> trying lmFit ~ plaque_location + sex + age")
  meta_cc <- meta_aligned[complete.cases(meta_aligned[, c("plaque_location", "sex", "age")]), , drop = FALSE]
  expr_cc <- expr_aligned[, rownames(meta_cc), drop = FALSE]
  design <- model.matrix(~ plaque_location + sex + age, data = meta_cc)
  if (coef_name %in% colnames(design)) {
    fit_lm <- lmFit(expr_cc, design)
    fit_lm <- eBayes(fit_lm)
    res <- topTable(fit_lm, coef = coef_name, number = Inf, sort.by = "p")
    if (!all(is.na(res$logFC))) {
      return(list(res = res, method = "lmFit"))
    }
  }

  # Strategy 3: lmFit with plaque_location only
  message("  -> trying lmFit ~ plaque_location")
  design_simple <- model.matrix(~ plaque_location, data = meta_cc)
  if (coef_name %in% colnames(design_simple)) {
    fit_simple <- lmFit(expr_cc, design_simple)
    fit_simple <- eBayes(fit_simple)
    res <- topTable(fit_simple, coef = coef_name, number = Inf, sort.by = "p")
    if (!all(is.na(res$logFC))) {
      return(list(res = res, method = "lmFit_simple"))
    }
  }

  # Not estimable in any model
  message(sprintf("  -> SKIP: coefficient %s not estimable in any model", coef_name))
  return(NULL)
}

form_de <- ~ plaque_location + sex + age + (1|dataset)

run_and_collect <- function(label, ref, comp, expr_data, meta, form_de) {
  message(sprintf("  Running %s...", label))
  result <- run_dream_contrast(expr_data, meta, ref, comp, form_de)
  if (is.null(result)) {
    message(sprintf("  %s: SKIPPED (not estimable)", label))
    return(NULL)
  }
  res <- result$res
  res$gene <- rownames(res)
  res$contrast <- label
  res$method <- result$method
  n_deg <- sum(res$adj.P.Val < 0.05 & abs(res$logFC) > 0.5)
  message(sprintf("  %s: %d DEGs (|logFC|>0.5, adj.P.Val<0.05) [%s]", label, n_deg, result$method))
  res
}

res_cc <- run_and_collect("Coronary_vs_Carotid", "carotid", "coronary", expr_subset, donor_meta, form_de)
res_fc <- run_and_collect("Femoral_vs_Carotid",   "carotid",  "femoral",  expr_subset, donor_meta, form_de)
res_cf <- run_and_collect("Coronary_vs_Femoral",  "femoral",  "coronary", expr_subset, donor_meta, form_de)

# Collect results
de_results <- list()
if (!is.null(res_cc)) de_results[["Coronary_vs_Carotid"]] <- res_cc
if (!is.null(res_fc)) de_results[["Femoral_vs_Carotid"]]   <- res_fc
if (!is.null(res_cf)) de_results[["Coronary_vs_Femoral"]]  <- res_cf

# Combine and save
de_all <- do.call(rbind, de_results)
write.csv(de_all, file.path(RES_DIR, "de_dream_results.csv"), row.names = FALSE)

# Volcano plots
n_contrasts <- length(de_results)
if (n_contrasts > 0) {
  pdf(file.path(RES_DIR, "de_volcano_dream.pdf"), width = 4.5 * n_contrasts, height = 5)
  par(mfrow = c(1, n_contrasts))
  for (i in seq_along(de_results)) {
    res <- de_results[[i]]
    res$sig <- with(res, adj.P.Val < 0.05 & abs(logFC) > 0.5)
    n_sig <- sum(res$sig)
    with(res, plot(logFC, -log10(P.Value),
      pch = 20, cex = 0.4, col = ifelse(sig, "#D55E00", "#999999"),
      main = sprintf("%s\n(%d DEGs)", names(de_results)[i], n_sig),
      xlab = "log2 Fold Change", ylab = "-log10(p)"
    ))
    abline(v = c(-0.5, 0.5), lty = 2, col = "#0072B2")
  }
  dev.off()
}

# -------------------------------------------------------------------
# 3. Linear Mixed Models for module scores (P0)
# -------------------------------------------------------------------
message("\n=== LMM for module scores ===")

donor_scores <- read.csv(file.path(RES_DIR, "donor_level_scores.csv"), row.names = 1)

# donor_scores already has plaque_location, dataset, sex — only need age from donor_meta
donor_df <- donor_scores
donor_df$donor_id <- rownames(donor_df)
meta_add <- donor_meta[, "age", drop = FALSE]
meta_add$donor_id <- rownames(meta_add)
donor_df <- merge(donor_df, meta_add, by = "donor_id")
rownames(donor_df) <- donor_df$donor_id

score_vars <- c("TI_composite", "TI_pca", "Inflammatory_Mac_score",
  "Resident_Mac_score", "Foamy_Mac_score", "Glycolysis_score",
  "FAO_score", "OXPHOS_score"
)
score_vars <- intersect(score_vars, names(donor_df))

lmm_results <- data.frame()
for (v in score_vars) {
  formula_str <- sprintf("%s ~ plaque_location + sex + age + (1|dataset)", v)
  tryCatch({
    model <- lmer(as.formula(formula_str), data = donor_df)
    summary_fixed <- summary(model)$coefficients
    anova_result <- anova(model)
    lmm_results <- rbind(lmm_results, data.frame(
      variable = v,
      plaque_location_F = anova_result["plaque_location", "F value"],
      plaque_location_p = anova_result["plaque_location", "Pr(>F)"],
      stringsAsFactors = FALSE
    ))
  }, error = function(e) {
    message(sprintf("  LMM failed for %s: %s", v, e$message))
  })
}

print(lmm_results)
write.csv(lmm_results, file.path(RES_DIR, "lmm_results.csv"), row.names = FALSE)

# -------------------------------------------------------------------
# 4. Lesion stage trend analysis (P2 temporal proxy)
# -------------------------------------------------------------------
message("\n=== Lesion stage analysis ===")

if ("lesion_stage" %in% names(donor_df)) {
  for (bed in c("carotid", "coronary", "femoral")) {
    bed_data <- donor_df[donor_df$plaque_location == bed, ]
    if (nrow(bed_data) >= 10) {
      model <- lmer(TI_composite ~ lesion_stage + sex + age + (1|dataset),
        data = bed_data
      )
      message(sprintf("  %s stage effect: F=%.2f, p=%.3e",
        bed,
        anova(model)["lesion_stage", "F value"],
        anova(model)["lesion_stage", "Pr(>F)"]
      ))
    }
  }
}

# -------------------------------------------------------------------
# 5. Healthy vs Disease PVAT correlation (P2)
# -------------------------------------------------------------------
message("\n=== Healthy vs Disease PVAT ===")

if (all(c("Healthy_PVAT_score", "Disease_PVAT_score") %in% names(donor_df))) {
  corr_test <- cor.test(donor_df$Healthy_PVAT_score, donor_df$Disease_PVAT_score,
    method = "spearman"
  )
  message(sprintf("  Spearman rho = %.3f, p = %.3e", corr_test$estimate, corr_test$p.value))

  pdf(file.path(RES_DIR, "healthy_vs_disease_pvat.pdf"), width = 7, height = 6)
  ggplot(donor_df, aes(x = Healthy_PVAT_score, y = Disease_PVAT_score, color = plaque_location)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
    labs(
      title = "Healthy vs Disease PVAT Score",
      subtitle = sprintf("Spearman rho = %.3f, p = %.3e", corr_test$estimate, corr_test$p.value)
    ) +
    theme_minimal(base_size = 14)
  dev.off()
}

message("\n=== R complement complete ===")
sink(file.path(RES_DIR, "sessionInfo.txt"))
sessionInfo()
sink()
