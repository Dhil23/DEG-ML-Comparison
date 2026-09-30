###############################################################################
#  MULTI-PIPELINE COMPARISON OF DIFFERENTIAL EXPRESSION ALGORITHMS
#  AND MACHINE LEARNING FOR BREAST CANCER BIOMARKER IDENTIFICATION
#
#  Dataset : GSE47462 (RNA-Seq 3SEQ, FFPE) — Normal / EN / DCIS / IDC
#            72 samples, 25 patients
#  Pipelines: DESeq2 x {RF, SVM} and EdgeR x {RF, SVM}
#  Validation: grouped leave-one-patient-out cross-validation
#
#  Mohammad Fadhil Ihsan (10622010) — SITH ITB
#  Supervisor: Popi Septiani, S.Si., M.Si., Ph.D.
#
#  Section order:
#    0  Libraries & configuration      8  Training 4 pipelines
#    1  Data import & metadata         9  Evaluation metrics
#    2  DESeq2 branch                 10  Significance testing
#    3  EdgeR branch                  11  Feature importance
#    4  Between-branch comparison     12  Core biomarkers
#    5  Artefact gene flagging        13  Figures
#    6  DEG sensitivity analyses      14  Reproducibility & summary
#    7  Feature sets & ML matrices
###############################################################################


# === SECTION 0: LIBRARIES & CONFIGURATION ====================================

# BiocManager::install(c("DESeq2","edgeR","EnhancedVolcano"))
# install.packages(c("caret","e1071","randomForest","kernlab","pROC","ggplot2",
#                    "pheatmap","reshape2","ggVennDiagram","patchwork","dplyr",
#                    "tidyr","scales"))

suppressPackageStartupMessages({
  library(DESeq2); library(edgeR); library(EnhancedVolcano); library(pheatmap)
  library(ggplot2); library(ggVennDiagram); library(reshape2); library(scales)
  library(dplyr); library(tidyr); library(caret); library(e1071)
  library(randomForest); library(pROC)
})

# Namespace clash: randomForest exports margin() and is attached AFTER
# ggplot2, so a bare margin(10, 26, 10, 26) dispatches to
# randomForest:::margin.default and fails with "26 is not a factor". Every
# margin() call below is therefore namespaced as ggplot2::margin. dplyr verbs
# masked by Bioconductor (select, filter, count, slice_*) are likewise prefixed.

DATA_PATH <- "GSE47462 Raw Counts Refseq Genes.txt"      # Adjust as needed

PADJ_CUTOFF  <- 0.01
FC_CUTOFF    <- 2.0
TOP_N        <- 20
N_TREE       <- 500
TUNE_L       <- 3
SEED         <- 42
B_BOOT       <- 1000     # bootstrap resamples, drawn at the patient level
DOM_CUTOFF   <- 0.50     # dominance threshold for artefact flagging
VOLC_TRIM    <- 0.999    # x-axis trimming quantile for volcano plots
JALAN_PAIRED <- TRUE     # run the ~ patient + condition sensitivity check

COND_LEVELS  <- c("Normal", "EN", "DCIS", "IDC")
PIPE_COLORS  <- c("DESeq2 + RF"  = "#B03A3A", "DESeq2 + SVM" = "#E08585",
                  "EdgeR + RF"   = "#0F3D35", "EdgeR + SVM"  = "#2AA181")
PIPE_LEVELS  <- names(PIPE_COLORS)
COND_COLORS  <- c(Normal = "#2196F3", EN = "#4CAF50",
                  DCIS = "#FF9800", IDC = "#F44336")
CLASS_COLORS <- c("#378ADD", "#3F9E52", "#EF9F27", "#E24B4A")

set.seed(SEED)

# Base R subsetting is used instead of dplyr::filter because dplyr drops row
# names, which would silently discard the gene identifiers.
get_deg <- function(res_table, p_col, fc_col,
                      p_cut = PADJ_CUTOFF, fc_cut = FC_CUTOFF) {
  d <- as.data.frame(res_table)
  passed <- !is.na(d[[p_col]]) & d[[p_col]] < p_cut &
           !is.na(d[[fc_col]]) & abs(d[[fc_col]]) > fc_cut
  d <- d[passed, , drop = FALSE]
  d <- d[order(d[[p_col]]), , drop = FALSE]
  d$Gene <- rownames(d)
  d
}


# === SECTION 1: DATA IMPORT & METADATA =======================================

DataRaw <- read.table(DATA_PATH, header = TRUE, row.names = 1,
                      sep = "\t", check.names = FALSE)
all_samples <- colnames(DataRaw)

normal_samples <- all_samples[grepl("normal$", all_samples)]
en_samples     <- all_samples[grepl("_EN$",    all_samples)]
dcis_samples   <- all_samples[grepl("DCIS$",   all_samples)]
idc_samples    <- all_samples[grepl("IDC$",    all_samples)]

selected_samples <- c(normal_samples, en_samples, dcis_samples, idc_samples)
stopifnot(!any(duplicated(selected_samples)))
counts_raw <- DataRaw[, selected_samples]

coldata <- data.frame(
  sample    = selected_samples,
  condition = factor(c(rep("Normal", length(normal_samples)),
                       rep("EN",     length(en_samples)),
                       rep("DCIS",   length(dcis_samples)),
                       rep("IDC",    length(idc_samples))),
                     levels = COND_LEVELS),
  stringsAsFactors = FALSE)
rownames(coldata) <- selected_samples
coldata$patient <- sapply(strsplit(selected_samples, "_"), `[`, 3)

y                 <- coldata$condition
classes           <- levels(y)
N_SAMPEL          <- nrow(coldata)
smallestGroupSize <- min(table(coldata$condition))

cat("Samples:", N_SAMPEL, "| Patients:", length(unique(coldata$patient)),
    "| Raw genes:", nrow(counts_raw), "\n")
print(table(coldata$condition))


# === SECTION 2: DESeq2 BRANCH — FILTER, NORMALIZATION, DEG, VST ==============

dds <- DESeqDataSetFromMatrix(countData = counts_raw, colData = coldata,
                              design = ~ condition)
keep_deseq <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep_deseq, ]
dds <- DESeq(dds)

cat("\n[DESeq2] genes passing filter:", nrow(dds), "of", nrow(counts_raw), "\n")

res_EN_DESeq2       <- results(dds, contrast = c("condition","EN","Normal"),   alpha = 0.05)
res_DCIS_DESeq2     <- results(dds, contrast = c("condition","DCIS","Normal"), alpha = 0.05)
res_IDC_DESeq2      <- results(dds, contrast = c("condition","IDC","Normal"),  alpha = 0.05)
res_IDCvDCIS_DESeq2 <- results(dds, contrast = c("condition","IDC","DCIS"),    alpha = 0.05)

DEG_EN_DESeq2   <- get_deg(res_EN_DESeq2,   "padj", "log2FoldChange")
DEG_DCIS_DESeq2 <- get_deg(res_DCIS_DESeq2, "padj", "log2FoldChange")
DEG_IDC_DESeq2  <- get_deg(res_IDC_DESeq2,  "padj", "log2FoldChange")

cat("[DESeq2] DEG — EN:", nrow(DEG_EN_DESeq2),
    "| DCIS:", nrow(DEG_DCIS_DESeq2), "| IDC:", nrow(DEG_IDC_DESeq2), "\n")

write.csv(DEG_EN_DESeq2,   "P_DESeq2_DEG_EN_vs_Normal.csv",   row.names = FALSE)
write.csv(DEG_DCIS_DESeq2, "P_DESeq2_DEG_DCIS_vs_Normal.csv", row.names = FALSE)
write.csv(DEG_IDC_DESeq2,  "P_DESeq2_DEG_IDC_vs_Normal.csv",  row.names = FALSE)
write.csv(as.data.frame(res_EN_DESeq2),   "P_DESeq2_full_results_EN.csv")
write.csv(as.data.frame(res_DCIS_DESeq2), "P_DESeq2_full_results_DCIS.csv")
write.csv(as.data.frame(res_IDC_DESeq2),  "P_DESeq2_full_results_IDC.csv")

vsd <- if (nrow(dds) >= 1000) vst(dds, blind = FALSE) else
       varianceStabilizingTransformation(dds, blind = FALSE)
mat_deseq2 <- assay(vsd)


# === SECTION 3: EdgeR BRANCH — FILTER, TMM, DEG, logCPM ======================

dge <- DGEList(counts = counts_raw, group = coldata$condition)
keep_edger <- filterByExpr(dge, group = coldata$condition)
dge <- dge[keep_edger, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")

cat("\n[EdgeR] genes passing filter:", nrow(dge), "of", nrow(counts_raw), "\n")

design_matrix <- model.matrix(~ 0 + condition, data = coldata)
colnames(design_matrix) <- levels(coldata$condition)
dge <- estimateDisp(dge, design_matrix)
fit <- glmQLFit(dge, design_matrix)

contrasts_matrix <- makeContrasts(
  EN_vs_Normal   = EN   - Normal,
  DCIS_vs_Normal = DCIS - Normal,
  IDC_vs_Normal  = IDC  - Normal,
  IDC_vs_DCIS    = IDC  - DCIS,
  levels = design_matrix)

res_EN_EdgeR   <- topTags(glmQLFTest(fit, contrast = contrasts_matrix[,"EN_vs_Normal"]),   n = Inf)$table
res_DCIS_EdgeR <- topTags(glmQLFTest(fit, contrast = contrasts_matrix[,"DCIS_vs_Normal"]), n = Inf)$table
res_IDC_EdgeR  <- topTags(glmQLFTest(fit, contrast = contrasts_matrix[,"IDC_vs_Normal"]),  n = Inf)$table

DEG_EN_EdgeR   <- get_deg(res_EN_EdgeR,   "FDR", "logFC")
DEG_DCIS_EdgeR <- get_deg(res_DCIS_EdgeR, "FDR", "logFC")
DEG_IDC_EdgeR  <- get_deg(res_IDC_EdgeR,  "FDR", "logFC")

cat("[EdgeR] DEG — EN:", nrow(DEG_EN_EdgeR),
    "| DCIS:", nrow(DEG_DCIS_EdgeR), "| IDC:", nrow(DEG_IDC_EdgeR), "\n")

write.csv(DEG_EN_EdgeR,   "P_EdgeR_DEG_EN_vs_Normal.csv",   row.names = FALSE)
write.csv(DEG_DCIS_EdgeR, "P_EdgeR_DEG_DCIS_vs_Normal.csv", row.names = FALSE)
write.csv(DEG_IDC_EdgeR,  "P_EdgeR_DEG_IDC_vs_Normal.csv",  row.names = FALSE)
write.csv(res_EN_EdgeR,   "P_EdgeR_full_results_EN.csv")
write.csv(res_DCIS_EdgeR, "P_EdgeR_full_results_DCIS.csv")
write.csv(res_IDC_EdgeR,  "P_EdgeR_full_results_IDC.csv")

mat_edger <- cpm(dge, log = TRUE, prior.count = 2)


# === SECTION 4: DESCRIPTIVE COMPARISON: DESeq2 vs EdgeR ======================
# Descriptive only. NOT used for ML feature selection.

contrast_names  <- c("EN", "DCIS", "IDC")
posthoc  <- list(
  DESeq2 = list(EN = DEG_EN_DESeq2$Gene, DCIS = DEG_DCIS_DESeq2$Gene,
                IDC = DEG_IDC_DESeq2$Gene),
  EdgeR  = list(EN = DEG_EN_EdgeR$Gene,  DCIS = DEG_DCIS_EdgeR$Gene,
                IDC = DEG_IDC_EdgeR$Gene))

jaccard_deg <- do.call(rbind, lapply(contrast_names, function(k) {
  a <- posthoc$DESeq2[[k]]; b <- posthoc$EdgeR[[k]]
  data.frame(Contrast = paste(k, "vs Normal"),
             DESeq2 = length(a), EdgeR = length(b),
             Intersection = length(intersect(a, b)), Union = length(union(a, b)),
             DESeq2_only = length(setdiff(a, b)),
             EdgeR_only  = length(setdiff(b, a)),
             Jaccard = round(length(intersect(a, b)) / length(union(a, b)), 4))
}))
cat("\n=== DEG OVERLAP: DESeq2 vs EdgeR ===\n"); print(jaccard_deg, row.names = FALSE)
write.csv(jaccard_deg, "P_DEG_Overlap_Jaccard.csv", row.names = FALSE)

deg_gradient <- data.frame(
  Stage  = contrast_names,
  DESeq2 = sapply(contrast_names, function(k) length(posthoc$DESeq2[[k]])),
  EdgeR  = sapply(contrast_names, function(k) length(posthoc$EdgeR[[k]])),
  row.names = NULL)
write.csv(deg_gradient, "P_DEG_Gradient_per_Stage.csv", row.names = FALSE)

branch_summary <- data.frame(
  Stage  = c("Genes after filter", "Normalization", "Test statistic",
             "Transformasi ML", "DEG EN", "DEG DCIS", "DEG IDC",
             "DEG IDC Up", "DEG IDC Down"),
  DESeq2 = c(nrow(dds), "median-of-ratios", "Wald (NB)", "VST",
             nrow(DEG_EN_DESeq2), nrow(DEG_DCIS_DESeq2), nrow(DEG_IDC_DESeq2),
             sum(DEG_IDC_DESeq2$log2FoldChange > 0),
             sum(DEG_IDC_DESeq2$log2FoldChange < 0)),
  EdgeR  = c(nrow(dge), "TMM", "Quasi-Likelihood F", "logCPM",
             nrow(DEG_EN_EdgeR), nrow(DEG_DCIS_EdgeR), nrow(DEG_IDC_EdgeR),
             sum(DEG_IDC_EdgeR$logFC > 0), sum(DEG_IDC_EdgeR$logFC < 0)),
  stringsAsFactors = FALSE)
print(branch_summary, row.names = FALSE)
write.csv(branch_summary, "P_Branch_Summary.csv", row.names = FALSE)


# === SECTION 5: ARTEFACT GENE FLAGGING =======================================
# Dominance = share of a gene's total counts contributed by a single sample.
# Genes are flagged, not removed.

total_counts   <- rowSums(counts_raw)
dominance  <- ifelse(total_counts > 0, apply(counts_raw, 1, max) / pmax(total_counts, 1), NA_real_)
n_detected <- rowSums(counts_raw >= 5)

all_deg <- unique(unlist(c(posthoc$DESeq2, posthoc$EdgeR)))

artefact_flag <- data.frame(Gene = rownames(counts_raw), Total_count = total_counts,
                           Dominance = round(dominance, 4),
                           N_samples_detected = n_detected, row.names = NULL) %>%
  dplyr::filter(Gene %in% all_deg) %>%
  dplyr::mutate(Suspect = Dominance > DOM_CUTOFF | N_samples_detected < 5) %>%
  dplyr::arrange(dplyr::desc(Dominance))

cat("\n=== ARTEFACT GENES ===\nUnique DEG:", length(all_deg),
    "| flagged as suspect:", sum(artefact_flag$Suspect), "\n")
write.csv(artefact_flag, "R_Artefact_Gene_Flags.csv", row.names = FALSE)


# === SECTION 6: DEG SENSITIVITY ANALYSES =====================================
# Companion analyses, not replacements. Primary results remain sections 2-3.

# --- 6.1 Native effect-size test (threshold enters the null hypothesis) ---
deg_native <- list(
  DESeq2 = setNames(lapply(contrast_names, function(k) {
    r <- results(dds, contrast = c("condition", k, "Normal"),
                 lfcThreshold = FC_CUTOFF, altHypothesis = "greaterAbs",
                 alpha = PADJ_CUTOFF)
    rownames(r)[!is.na(r$padj) & r$padj < PADJ_CUTOFF]
  }), contrast_names),
  EdgeR = setNames(lapply(contrast_names, function(k) {
    tt <- topTags(glmTreat(fit, contrast = contrasts_matrix[, paste0(k, "_vs_Normal")],
                           lfc = FC_CUTOFF), n = Inf)$table
    rownames(tt)[tt$FDR < PADJ_CUTOFF]
  }), contrast_names))

native_vs_posthoc <- do.call(rbind, lapply(names(posthoc), function(cb)
  do.call(rbind, lapply(contrast_names, function(k) {
    a <- posthoc[[cb]][[k]]; b <- deg_native[[cb]][[k]]
    data.frame(Branch = cb, Contrast = k, Post_hoc = length(a), Native = length(b),
               Intersection = length(intersect(a, b)),
               Post_hoc_only = length(setdiff(a, b)),
               Native_only   = length(setdiff(b, a)))
  }))))
cat("\n=== NATIVE vs POST-HOC TEST ===\n"); print(native_vs_posthoc, row.names = FALSE)
write.csv(native_vs_posthoc, "R_Native_vs_PostHoc.csv", row.names = FALSE)

# --- 6.2 Paired design ~ patient + condition ---
if (JALAN_PAIRED) {
  coldat_p <- as.data.frame(SummarizedExperiment::colData(dds))
  coldat_p$patient <- factor(coldat_p$patient)
  mm <- model.matrix(~ patient + condition, data = coldat_p)
  full_rank <- qr(mm)$rank == ncol(mm)
  keep_cols <- seq_len(ncol(dds))

  if (!full_rank) {
    n_per_patient   <- table(coldat_p$patient)
    drop_pat  <- names(n_per_patient)[n_per_patient < 2]
    keep_cols <- which(!(coldat_p$patient %in% drop_pat))
    mm <- model.matrix(~ patient + condition, data = droplevels(coldat_p[keep_cols, ]))
    full_rank <- qr(mm)$rank == ncol(mm)
    cat("\n[paired] dropping single-sample patients:",
        paste(drop_pat, collapse = ", "), "\n")
  }
  cat("[paired] design", nrow(mm), "x", ncol(mm), "| rank =", qr(mm)$rank,
      "|", ifelse(full_rank, "full rank", "rank-deficient"), "\n")

  if (full_rank) {
    dds_p <- dds[, keep_cols]
    SummarizedExperiment::colData(dds_p)$patient <-
      droplevels(factor(SummarizedExperiment::colData(dds_p)$patient))
    design(dds_p) <- ~ patient + condition
    dds_p <- DESeq(dds_p, quiet = TRUE)

    deg_paired <- setNames(lapply(contrast_names, function(k)
      get_deg(results(dds_p, contrast = c("condition", k, "Normal"),
                        alpha = PADJ_CUTOFF), "padj", "log2FoldChange")$Gene),
      contrast_names)

    paired_comparison <- do.call(rbind, lapply(contrast_names, function(k) {
      a <- posthoc$DESeq2[[k]]; b <- deg_paired[[k]]
      data.frame(Contrast = k, Tanpa_pasien = length(a), Dengan_pasien = length(b),
                 Intersection = length(intersect(a, b)))
    }))
    cat("\n=== PAIRED DESIGN (DESeq2) ===\n")
    print(paired_comparison, row.names = FALSE)
    write.csv(paired_comparison, "R_Paired_Design_DESeq2.csv", row.names = FALSE)
  }
}


# === SECTION 7: FEATURE SETS & ML MATRICES ===================================
# Union of DEG from three contrasts against Normal, per branch. No consensus.

genes_deseq2 <- intersect(unique(unlist(posthoc$DESeq2)), rownames(mat_deseq2))
genes_edger  <- intersect(unique(unlist(posthoc$EdgeR)),  rownames(mat_edger))

if (length(genes_deseq2) < 2 || length(genes_edger) < 2)
  stop("Feature set too small. Relax PADJ_CUTOFF / FC_CUTOFF.")

X_deseq2 <- t(mat_deseq2[genes_deseq2, rownames(coldata), drop = FALSE])
X_edger  <- t(mat_edger[genes_edger,   rownames(coldata), drop = FALSE])
stopifnot(identical(rownames(X_deseq2), rownames(coldata)),
          identical(rownames(X_edger),  rownames(coldata)))

cat("\nFitur ML — DESeq2:", ncol(X_deseq2), "| EdgeR:", ncol(X_edger),
    "| shared:", length(intersect(genes_deseq2, genes_edger)), "\n")

write.csv(data.frame(
  Branch = c("DESeq2","EdgeR"), Matrix = c("VST","logCPM (TMM)"),
  Genes_after_filter = c(nrow(dds), nrow(dge)),
  ML_features = c(ncol(X_deseq2), ncol(X_edger))),
  "P_Feature_Composition.csv", row.names = FALSE)


# === SECTION 8: TRAINING 4 PIPELINES — GROUPED LOOCV =========================
# Folds are grouped by PATIENT: every sample from one patient is held out together.

patients      <- unique(coldata$patient)
train_indices <- lapply(patients, function(p) which(coldata$patient != p))
test_indices  <- lapply(patients, function(p) which(coldata$patient == p))
names(train_indices) <- names(test_indices) <- paste0("Patient_", patients)

ctrl_loocv <- trainControl(method = "cv", index = train_indices,
                           indexOut = test_indices, classProbs = TRUE,
                           savePredictions = "final")

train_model <- function(X, y, algorithm, label_name) {
  cat("  -", label_name, "... ")
  w <- system.time({
    set.seed(SEED)
    m <- if (algorithm == "rf")
      train(x = X, y = y, method = "rf", trControl = ctrl_loocv,
            ntree = N_TREE, tuneLength = TUNE_L)
    else
      train(x = X, y = y, method = "svmRadial", trControl = ctrl_loocv,
            preProcess = c("center","scale"), tuneLength = TUNE_L)
  })
  cat("done (", round(w["elapsed"], 1), "s )\n"); m
}

cat("\n=== TRAINING 4 PIPELINES (", length(patients), "patient folds ) ===\n")
models <- list(
  "DESeq2 + RF"  = train_model(X_deseq2, y, "rf",  "DESeq2 + RF"),
  "DESeq2 + SVM" = train_model(X_deseq2, y, "svm", "DESeq2 + SVM"),
  "EdgeR + RF"   = train_model(X_edger,  y, "rf",  "EdgeR + RF"),
  "EdgeR + SVM"  = train_model(X_edger,  y, "svm", "EdgeR + SVM"))
preds <- lapply(models, function(m) m$pred)


# === SECTION 9: EVALUATION METRICS ===========================================

auc_per_class <- function(d) sapply(classes, function(k)
  as.numeric(auc(roc(as.numeric(d$obs == k), d[[k]],
                     quiet = TRUE, direction = "<"))))

cm_list <- lapply(preds, function(p) confusionMatrix(p$pred, p$obs))

metrics_all <- do.call(rbind, lapply(PIPE_LEVELS, function(nm) {
  cm <- cm_list[[nm]]; a <- auc_per_class(preds[[nm]])
  data.frame(Pipeline = nm,
             DE_method = ifelse(grepl("DESeq2", nm), "DESeq2", "EdgeR"),
             ML_model   = ifelse(grepl("RF", nm), "Random Forest", "SVM"),
             N_features = ifelse(grepl("DESeq2", nm), ncol(X_deseq2), ncol(X_edger)),
             Accuracy = as.numeric(cm$overall["Accuracy"]),
             Kappa    = as.numeric(cm$overall["Kappa"]),
             Balanced_Accuracy = mean(cm$byClass[, "Balanced Accuracy"], na.rm = TRUE),
             Macro_F1 = mean(cm$byClass[, "F1"], na.rm = TRUE),
             Macro_AUC = mean(a), stringsAsFactors = FALSE)
}))
rownames(metrics_all) <- NULL

num_cols <- c("Accuracy","Kappa","Balanced_Accuracy","Macro_F1","Macro_AUC")
metrics_print <- metrics_all; metrics_print[num_cols] <- round(metrics_print[num_cols], 3)

cat("\n=== METRICS, 4 PIPELINES ===\n")
print(metrics_print[, c("Pipeline","N_features", num_cols)], row.names = FALSE)
for (nm in PIPE_LEVELS) { cat("\n>>", nm, "\n"); print(cm_list[[nm]]$table) }

auc_detail <- do.call(rbind, lapply(PIPE_LEVELS, function(nm) {
  a <- auc_per_class(preds[[nm]])
  data.frame(Pipeline = nm, t(round(a, 3)), Macro = round(mean(a), 3),
             check.names = FALSE)
}))
rownames(auc_detail) <- NULL
print(auc_detail, row.names = FALSE)

write.csv(metrics_print, "P_Metrics_4Pipelines_LOOCV.csv", row.names = FALSE)
write.csv(auc_detail,   "P_AUC_per_Kelas.csv",          row.names = FALSE)


# === SECTION 10: SIGNIFICANCE TESTING ========================================
# caret orders m$pred by resample, not by sample. Aligning on rowIndex is
# mandatory before paired = TRUE.

pr <- lapply(preds, function(d) d[order(d$rowIndex), , drop = FALSE])
stopifnot(length(unique(lapply(pr, function(d) d$rowIndex))) == 1)
patient_of_row <- coldata$patient[pr[[1]]$rowIndex]

roc_obj <- function(d, k) {
  o <- as.numeric(d$obs == k)
  if (length(unique(o)) < 2) return(NULL)
  roc(o, d[[k]], quiet = TRUE, direction = "<")
}

# --- 10.1 DeLong confidence intervals ---
auc_ci <- do.call(rbind, lapply(PIPE_LEVELS, function(nm)
  do.call(rbind, lapply(classes, function(k) {
    r <- roc_obj(pr[[nm]], k); if (is.null(r)) return(NULL)
    ci <- as.numeric(ci.auc(r, method = "delong"))
    data.frame(Pipeline = nm, Class = k, AUC = ci[2],
               CI_lower = ci[1], CI_upper = ci[3])
  }))))
cat("\n=== AUC + 95% CI DeLong ===\n")
print(as.data.frame(auc_ci %>% dplyr::mutate(
  dplyr::across(dplyr::where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
write.csv(auc_ci, "R_AUC_CI_DeLong.csv", row.names = FALSE)

# --- 10.2 Paired DeLong test between pipelines ---
pair_mat <- t(combn(PIPE_LEVELS, 2))
delong_test <- do.call(rbind, lapply(seq_len(nrow(pair_mat)), function(i)
  do.call(rbind, lapply(classes, function(k) {
    r1 <- roc_obj(pr[[pair_mat[i,1]]], k); r2 <- roc_obj(pr[[pair_mat[i,2]]], k)
    if (is.null(r1) || is.null(r2)) return(NULL)
    tt <- try(roc.test(r1, r2, method = "delong", paired = TRUE), silent = TRUE)
    if (inherits(tt, "try-error")) return(NULL)
    data.frame(Class = k, Pipeline_A = pair_mat[i,1], Pipeline_B = pair_mat[i,2],
               AUC_A = as.numeric(auc(r1)), AUC_B = as.numeric(auc(r2)),
               Difference = as.numeric(auc(r1)) - as.numeric(auc(r2)), p = tt$p.value)
  }))))
delong_test$p_adj <- p.adjust(delong_test$p, method = "BH")
cat("\n=== PAIRED DeLong TEST ===\nPassing p_adj < 0.05:",
    sum(delong_test$p_adj < 0.05), "of", nrow(delong_test), "tests\n")
print(as.data.frame(delong_test %>% dplyr::arrange(p) %>% dplyr::slice_head(n = 8) %>%
  dplyr::mutate(Difference = round(Difference, 3), p = signif(p, 3),
                p_adj = signif(p_adj, 3)) %>%
  dplyr::select(Class, Pipeline_A, Pipeline_B, Difference, p, p_adj)),
  row.names = FALSE)
write.csv(delong_test, "R_RocTest_DeLong.csv", row.names = FALSE)

# --- 10.3 Patient-level bootstrap (no model refitting) ---
metrics_from_idx <- function(d, idx) {
  cm <- confusionMatrix(factor(d$pred[idx], levels = classes),
                        factor(d$obs[idx],  levels = classes))
  mauc <- mean(sapply(classes, function(k) {
    o <- as.numeric(d$obs[idx] == k)
    if (length(unique(o)) < 2) return(NA_real_)
    as.numeric(auc(roc(o, d[[k]][idx], quiet = TRUE, direction = "<")))
  }), na.rm = TRUE)
  c(Accuracy = as.numeric(cm$overall["Accuracy"]),
    Kappa    = as.numeric(cm$overall["Kappa"]),
    Balanced_Accuracy = mean(cm$byClass[,"Balanced Accuracy"], na.rm = TRUE),
    Macro_F1 = mean(cm$byClass[,"F1"], na.rm = TRUE), Macro_AUC = mauc)
}

set.seed(SEED); up <- unique(patient_of_row)
cat("\n=== PATIENT-LEVEL BOOTSTRAP (B =", B_BOOT, ") ===\nRunning ")
bt <- lapply(seq_len(B_BOOT), function(b) {
  if (b %% 100 == 0) cat(".")
  idx <- unlist(lapply(sample(up, length(up), replace = TRUE),
                       function(p) which(patient_of_row == p)))
  sapply(PIPE_LEVELS, function(nm)
    tryCatch(metrics_from_idx(pr[[nm]], idx),
             error = function(e) rep(NA_real_, 5)))
})
cat(" done\n")

boot_ci <- do.call(rbind, lapply(PIPE_LEVELS, function(nm)
  do.call(rbind, lapply(seq_along(num_cols), function(j) {
    v  <- sapply(bt, function(m) m[j, nm])
    ci <- quantile(v, c(0.025, 0.975), na.rm = TRUE)
    data.frame(Pipeline = nm, Metric = num_cols[j],
               Point = as.numeric(metrics_all[metrics_all$Pipeline == nm, num_cols[j]]),
               CI_lower = ci[1], CI_upper = ci[2])
  }))))
print(as.data.frame(boot_ci %>% dplyr::mutate(
  dplyr::across(dplyr::where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
write.csv(boot_ci, "R_Bootstrap_CI_Metrics.csv", row.names = FALSE)

boot_diff <- do.call(rbind, lapply(seq_len(nrow(pair_mat)), function(i)
  do.call(rbind, lapply(seq_along(num_cols), function(j) {
    d  <- sapply(bt, function(m) m[j, pair_mat[i,1]] - m[j, pair_mat[i,2]])
    ci <- quantile(d, c(0.025, 0.975), na.rm = TRUE)
    data.frame(Metric = num_cols[j], Pipeline_A = pair_mat[i,1],
               Pipeline_B = pair_mat[i,2],
               Selisih_median = median(d, na.rm = TRUE),
               CI_lower = ci[1], CI_upper = ci[2],
               Contains_zero = ci[1] <= 0 && ci[2] >= 0)
  }))))
cat("Pipeline differences whose CI excludes zero:",
    sum(!boot_diff$Contains_zero), "of", nrow(boot_diff), "\n")
write.csv(boot_diff, "R_Bootstrap_Pipeline_Differences.csv", row.names = FALSE)


# === SECTION 11: FEATURE IMPORTANCE ==========================================
# varImp for svmRadial is NOT model-based: caret falls back to a univariate
# ROC filter. The Importance_method column prevents misreading the table later.

imp_full <- setNames(lapply(PIPE_LEVELS, function(nm) {
  im <- varImp(models[[nm]], scale = TRUE)$importance
  im$Gene <- rownames(im)
  sc <- setdiff(colnames(im), "Gene")
  im$Mean_Importance <- rowMeans(im[, sc, drop = FALSE], na.rm = TRUE)
  im[order(-im$Mean_Importance), c("Gene","Mean_Importance")]
}), PIPE_LEVELS)

imp_list <- setNames(lapply(PIPE_LEVELS, function(nm)
  head(imp_full[[nm]], TOP_N) %>%
    dplyr::mutate(Pipeline = nm, Rank = dplyr::row_number(),
                  Importance_method = ifelse(grepl("RF", nm),
                    "Mean decrease in impurity (model-based)",
                    "Univariate ROC filter (NOT model-based)"))),
  PIPE_LEVELS)
imp_all <- do.call(rbind, imp_list); rownames(imp_all) <- NULL

cat("\n=== TOP 10 GENES PER PIPELINE ===\n")
for (nm in PIPE_LEVELS) {
  cat("\n>>", nm, "\n")
  print(as.data.frame(imp_list[[nm]][1:10, c("Rank","Gene","Mean_Importance")]),
        row.names = FALSE)
}
write.csv(imp_all, "P_Feature_Importance.csv", row.names = FALSE)

# Sensitivity of the four-way intersection to TOP_N
topn_sensitivity <- do.call(rbind, lapply(c(5,10,15,20,25,30,40,50,75,100), function(n) {
  ir <- Reduce(intersect, lapply(PIPE_LEVELS, function(nm) head(imp_full[[nm]]$Gene, n)))
  data.frame(TOP_N = n, Intersection = length(ir),
             Genes = ifelse(length(ir) == 0, "-", paste(head(ir, 12), collapse = ", ")))
}))
cat("\n=== TOP_N SENSITIVITY ===\n"); print(topn_sensitivity, row.names = FALSE)
write.csv(topn_sensitivity, "R_TOPN_Sensitivity.csv", row.names = FALSE)


# === SECTION 12: CORE BIOMARKERS =============================================
# Read AFTER all four models are trained. This is not feature selection.

venn_sets <- setNames(lapply(imp_list, function(d) d$Gene), PIPE_LEVELS)
core_biomarker <- Reduce(intersect, venn_sets)

cat("\n=== CORE BIOMARKERS (intersection of top-", TOP_N, ") ===\n", sep = "")
if (length(core_biomarker) > 0) {
  cat("Count:", length(core_biomarker), "genes —",
      paste(core_biomarker, collapse = ", "), "\n")
  core_table <- imp_all %>%
    dplyr::filter(Gene %in% core_biomarker) %>%
    dplyr::select(Gene, Pipeline, Rank) %>%
    tidyr::pivot_wider(names_from = Pipeline, values_from = Rank) %>%
    dplyr::mutate(Mean_rank = rowMeans(dplyr::select(., -Gene), na.rm = TRUE)) %>%
    dplyr::arrange(Mean_rank)
  print(as.data.frame(core_table), row.names = FALSE)
  write.csv(core_table, "P_Core_Biomarker.csv", row.names = FALSE)
} else {
  cat("Empty intersection at TOP_N =", TOP_N, "\n")
}

overlap_matrix <- outer(PIPE_LEVELS, PIPE_LEVELS, Vectorize(function(a, b)
  length(intersect(venn_sets[[a]], venn_sets[[b]]))))
dimnames(overlap_matrix) <- list(PIPE_LEVELS, PIPE_LEVELS)
print(overlap_matrix)
write.csv(as.data.frame(overlap_matrix), "P_Overlap_Matrix.csv")


# === SECTION 13: FIGURES =====================================================

save_plot <- function(p, f, w, h) ggsave(f, p, width = w, height = h, dpi = 300,
                                      bg = "white")

# --- 13.1 PCA (no per-sample labels) ---
pca_d  <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
var_d  <- round(100 * attr(pca_d, "percentVar"), 1)
save_plot(ggplot(pca_d, aes(PC1, PC2, colour = condition)) +
  geom_point(size = 3.2, alpha = 0.85) +
  scale_colour_manual(values = COND_COLORS) +
  labs(title = "PCA — DESeq2 branch (VST)", colour = "Stage",
       subtitle = "Sample labels omitted to keep the point distribution readable",
       x = paste0("PC1 (", var_d[1], "% variance)"),
       y = paste0("PC2 (", var_d[2], "% variance)")) +
  theme_bw(base_size = 12) + theme(plot.title = element_text(face = "bold")),
  "P_DESeq2_PCA.png", 8, 6.5)

top_var <- head(order(apply(mat_edger, 1, var), decreasing = TRUE),
                min(500, nrow(mat_edger)))
pca_e   <- prcomp(t(mat_edger[top_var, ]), scale. = FALSE)
var_e   <- round(100 * pca_e$sdev^2 / sum(pca_e$sdev^2), 1)
save_plot(ggplot(data.frame(PC1 = pca_e$x[,1], PC2 = pca_e$x[,2],
                         condition = coldata$condition),
              aes(PC1, PC2, colour = condition)) +
  geom_point(size = 3.2, alpha = 0.85) +
  scale_colour_manual(values = COND_COLORS) +
  labs(title = "PCA — EdgeR branch (logCPM, TMM)", colour = "Stage",
       subtitle = "500 most variable genes",
       x = paste0("PC1 (", var_e[1], "% variance)"),
       y = paste0("PC2 (", var_e[2], "% variance)")) +
  theme_bw(base_size = 12) + theme(plot.title = element_text(face = "bold")),
  "P_EdgeR_PCA.png", 8, 6.5)

# --- 13.2 Volcano (identical axis limits across all six panels) ---
volc_all <- list(list(res_EN_DESeq2,"log2FoldChange","padj"),
                   list(res_DCIS_DESeq2,"log2FoldChange","padj"),
                   list(res_IDC_DESeq2,"log2FoldChange","padj"),
                   list(res_EN_EdgeR,"logFC","FDR"),
                   list(res_DCIS_EdgeR,"logFC","FDR"),
                   list(res_IDC_EdgeR,"logFC","FDR"))
.fc <- unlist(lapply(volc_all, function(v) as.data.frame(v[[1]])[[v[[2]]]]))
.pv <- unlist(lapply(volc_all, function(v) as.data.frame(v[[1]])[[v[[3]]]]))
.fc <- .fc[is.finite(.fc)]; .pv <- .pv[is.finite(.pv) & .pv > 0]
VOLC_XMAX <- ceiling(as.numeric(quantile(abs(.fc), VOLC_TRIM, na.rm = TRUE)))
VOLC_XLIM <- c(-VOLC_XMAX, VOLC_XMAX)
VOLC_YLIM <- c(0, ceiling(max(-log10(.pv), na.rm = TRUE) * 1.08))
cat("\nSumbu volcano:", paste(VOLC_XLIM, collapse = " s/d "),
    "| y 0 to", VOLC_YLIM[2], "|", sum(abs(.fc) > VOLC_XMAX), "points outside\n")

make_volcano <- function(res_df, kfc, kp, plot_title, file) {
  d <- as.data.frame(res_df)
  nu <- sum(d[[kfc]] >  FC_CUTOFF & d[[kp]] < PADJ_CUTOFF, na.rm = TRUE)
  nd <- sum(d[[kfc]] < -FC_CUTOFF & d[[kp]] < PADJ_CUTOFF, na.rm = TRUE)
  png(file, width = 2400, height = 2000, res = 300)
  print(EnhancedVolcano(d, lab = rownames(d), x = kfc, y = kp, title = plot_title,
    subtitle = paste0(kp," < ",PADJ_CUTOFF,"  |  |",kfc,"| > ",FC_CUTOFF,
                      "  |  Up: ",nu,"  |  Down: ",nd),
    xlim = VOLC_XLIM, ylim = VOLC_YLIM, pCutoff = PADJ_CUTOFF,
    FCcutoff = FC_CUTOFF, pointSize = 1.5, labSize = 3,
    col = c("grey70","#2196F3","#FF9800","#F44336"), colAlpha = 0.7,
    legendLabels = c("NS","log2FC","p-adj","p-adj & log2FC"),
    drawConnectors = TRUE))
  dev.off()
}
make_volcano(res_EN_DESeq2,  "log2FoldChange","padj","EN vs Normal (DESeq2)",  "P_DESeq2_Volcano_EN.png")
make_volcano(res_DCIS_DESeq2,"log2FoldChange","padj","DCIS vs Normal (DESeq2)","P_DESeq2_Volcano_DCIS.png")
make_volcano(res_IDC_DESeq2, "log2FoldChange","padj","IDC vs Normal (DESeq2)", "P_DESeq2_Volcano_IDC.png")
make_volcano(res_EN_EdgeR,   "logFC","FDR","EN vs Normal (EdgeR)",   "P_EdgeR_Volcano_EN.png")
make_volcano(res_DCIS_EdgeR, "logFC","FDR","DCIS vs Normal (EdgeR)", "P_EdgeR_Volcano_DCIS.png")
make_volcano(res_IDC_EdgeR,  "logFC","FDR","IDC vs Normal (EdgeR)",  "P_EdgeR_Volcano_IDC.png")

# --- 13.3 Heatmap Top 50 DEG IDC ---
annot_col <- data.frame(Stage = coldata$condition, row.names = rownames(coldata))
make_heatmap <- function(mat, genes_sel, plot_title, file) {
  g <- intersect(genes_sel, rownames(mat)); if (length(g) < 2) return(invisible())
  png(file, width = 2600, height = 3600, res = 300)
  pheatmap(mat[g, ], annotation_col = annot_col,
           annotation_colors = list(Stage = COND_COLORS), scale = "row",
           clustering_method = "ward.D2", show_colnames = FALSE,
           fontsize_row = 7, main = plot_title,
           color = colorRampPalette(c("#2166AC","white","#D6604D"))(100))
  dev.off()
}
make_heatmap(mat_deseq2, head(DEG_IDC_DESeq2$Gene, 50),
             "Top 50 DEG IDC vs Normal — DESeq2 (VST, row z-score)",
             "P_DESeq2_Heatmap_Top50.png")
make_heatmap(mat_edger, head(DEG_IDC_EdgeR$Gene, 50),
             "Top 50 DEG IDC vs Normal — EdgeR (logCPM, row z-score)",
             "P_EdgeR_Heatmap_Top50.png")

# --- 13.4 Venn of DEG per contrast ---
# ggVennDiagram places set labels OUTSIDE the circles and the ggplot panel
# clips them, so "DESeq2" reads as "ESeq2" and "EdgeR" as "dgeR". The fix
# widens both scales so the labels fall inside the panel, plus a plot margin.
# Do not add coord_sf() here: it REPLACES the coordinate system ggVennDiagram
# already installed. If labels still clip, raise mult in expansion() instead.
#
make_venn2 <- function(a, b, plot_title, lo, hi, file) {
  if (length(a) == 0 || length(b) == 0) return(invisible())
  save_plot(ggVennDiagram(list(DESeq2 = a, EdgeR = b), label_alpha = 0,
                       label = "count", set_size = 4.2, label_size = 3.8) +
    scale_fill_gradient(low = lo, high = hi) +
    scale_x_continuous(expand = expansion(mult = 0.30)) +
    scale_y_continuous(expand = expansion(mult = 0.18)) +
    labs(title = plot_title,
         subtitle = "Descriptive only — not used for ML feature selection") +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"),
          plot.margin = ggplot2::margin(10, 26, 10, 26)),
    file, 8, 6.5)
}
make_venn2(DEG_EN_DESeq2$Gene,  DEG_EN_EdgeR$Gene,
           "DEG EN vs Normal — DESeq2 vs EdgeR",  "#E8F5E9","#2E7D32","P_Venn_DEG_EN.png")
make_venn2(DEG_DCIS_DESeq2$Gene,DEG_DCIS_EdgeR$Gene,
           "DEG DCIS vs Normal — DESeq2 vs EdgeR","#FFF3E0","#E65100","P_Venn_DEG_DCIS.png")
make_venn2(DEG_IDC_DESeq2$Gene, DEG_IDC_EdgeR$Gene,
           "DEG IDC vs Normal — DESeq2 vs EdgeR", "#E8F4FD","#1565C0","P_Venn_DEG_IDC.png")

# --- 13.5 ROC ---
macro_roc <- function(d, label) {
  gr <- seq(0, 1, length.out = 100)
  tp <- sapply(classes, function(k) {
    r <- roc_obj(d, k); if (is.null(r)) return(rep(NA_real_, 100))
    f <- 1 - r$specificities; t <- r$sensitivities; s <- order(f)
    approx(f[s], t[s], xout = gr, ties = "mean", rule = 2)$y })
  data.frame(FPR = gr, TPR = rowMeans(tp, na.rm = TRUE), Pipeline = label)
}
roc_all <- do.call(rbind, lapply(PIPE_LEVELS, function(nm) macro_roc(preds[[nm]], nm)))
roc_all$Pipeline <- factor(roc_all$Pipeline, levels = PIPE_LEVELS)

save_plot(ggplot(roc_all, aes(FPR, TPR, colour = Pipeline)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey55") +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = PIPE_COLORS,
    labels = setNames(paste0(metrics_all$Pipeline, "  (AUC = ",
                             sprintf("%.3f", metrics_all$Macro_AUC), ")"),
                      metrics_all$Pipeline)) +
  coord_equal(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
  labs(title = "ROC curves — four pipelines, grouped LOOCV", colour = NULL,
       subtitle = paste0("Macro-average one-vs-rest | ", N_SAMPEL,
                         " samples from ", length(patients), " patients"),
       x = "False positive rate (1 - specificity)",
       y = "True positive rate (sensitivity)") +
  theme_bw(base_size = 12) + theme(plot.title = element_text(face = "bold")),
  "P_ROC_Combined.png", 10, 6.5)

roc_class <- do.call(rbind, lapply(PIPE_LEVELS, function(nm)
  do.call(rbind, lapply(classes, function(k) {
    r <- roc_obj(preds[[nm]], k); if (is.null(r)) return(NULL)
    data.frame(FPR = 1 - r$specificities, TPR = r$sensitivities, Class = k,
               AUC = as.numeric(auc(r)), Pipeline = nm) }))))
roc_class <- roc_class[order(roc_class$Pipeline, roc_class$Class, roc_class$FPR), ]
roc_class$Pipeline <- factor(roc_class$Pipeline, levels = PIPE_LEVELS)
roc_class$Class    <- factor(roc_class$Class,    levels = classes)

auc_annot <- roc_class %>% dplyr::distinct(Pipeline, Class, AUC) %>%
  dplyr::arrange(Pipeline, Class) %>% dplyr::group_by(Pipeline) %>%
  dplyr::mutate(ypos = 0.30 - 0.075 * (dplyr::row_number() - 1)) %>% dplyr::ungroup()

save_plot(ggplot(roc_class, aes(FPR, TPR, colour = Class)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey60") +
  geom_line(linewidth = 0.95) +
  geom_text(data = auc_annot, inherit.aes = FALSE,
            aes(0.98, ypos, colour = Class,
                label = sprintf("%s: AUC = %.3f", Class, AUC)),
            hjust = 1, size = 3, show.legend = FALSE) +
  scale_colour_manual(values = setNames(CLASS_COLORS, classes)) +
  coord_equal(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
  facet_wrap(~ Pipeline, nrow = 2) +
  labs(title = "One-vs-rest ROC curves by stage", colour = NULL,
       subtitle = "AUC annotated inside each panel",
       x = "False positive rate", y = "True positive rate") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"), legend.position = "top"),
  "P_ROC_per_Class.png", 10, 9.5)

# --- 13.6 Metrics: barplot, heatmap, bootstrap CI ---
# Column names are mapped to display labels so the figures do not carry
# variable-style underscores. CSV files keep the original machine-readable names.
MET_LABEL <- c(Accuracy = "Accuracy", Kappa = "Cohen's kappa",
               Balanced_Accuracy = "Balanced accuracy",
               Macro_F1 = "Macro F1", Macro_AUC = "Macro AUC")

met_long <- metrics_all %>%
  dplyr::select(Pipeline, dplyr::all_of(num_cols)) %>%
  tidyr::pivot_longer(-Pipeline, names_to = "Metric", values_to = "Value") %>%
  dplyr::mutate(Pipeline = factor(Pipeline, levels = PIPE_LEVELS),
                Metric   = factor(MET_LABEL[Metric], levels = MET_LABEL))

save_plot(ggplot(met_long, aes(Metric, Value, fill = Pipeline)) +
  geom_col(position = position_dodge(0.8), width = 0.72) +
  geom_text(aes(label = sprintf("%.3f", Value)), position = position_dodge(0.8),
            vjust = -0.4, size = 2.7) +
  scale_fill_manual(values = PIPE_COLORS) +
  scale_y_continuous(limits = c(0, 1.08), expand = c(0,0)) +
  labs(title = "Evaluation metrics across four pipelines", x = NULL,
       y = "Value", fill = "Pipeline",
       subtitle = paste0("Grouped LOOCV over ", length(patients), " patients")) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top",
        panel.grid.major.x = element_blank()),
  "P_Metrics_Barplot.png", 11, 6)

save_plot(met_long %>% dplyr::group_by(Metric) %>%
  dplyr::mutate(rt = max(Value) - min(Value),
                ns = ifelse(rt == 0, 0.5, (Value - min(Value)) / rt)) %>%
  dplyr::ungroup() %>%
  ggplot(aes(Metric, Pipeline, fill = ns)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = sprintf("%.3f", Value), colour = ns > 0.6),
            size = 4, fontface = "bold", show.legend = FALSE) +
  scale_fill_gradient(low = "#EAF2F8", high = "#1A5C8A") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15")) +
  scale_y_discrete(limits = rev(PIPE_LEVELS)) +
  labs(title = "Performance heatmap", x = NULL, y = NULL,
       fill = "Relative\n(per metric)",
       subtitle = "Colour scaled within each column; printed values are raw") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), panel.grid = element_blank(),
        axis.text.y = element_text(face = "bold")),
  "P_Metrics_Heatmap.png", 9, 5)

save_plot(boot_ci %>%
  dplyr::mutate(Pipeline = factor(Pipeline, levels = rev(PIPE_LEVELS)),
                Metric   = factor(MET_LABEL[Metric], levels = MET_LABEL)) %>%
  ggplot(aes(Point, Pipeline, colour = Pipeline)) +
  geom_linerange(aes(xmin = CI_lower, xmax = CI_upper), linewidth = 0.8) +
  geom_point(size = 2.7) +
  geom_text(aes(label = sprintf("%.3f [%.3f, %.3f]", Point, CI_lower, CI_upper)),
            vjust = -1.05, size = 2.3, colour = "grey20") +
  scale_colour_manual(values = PIPE_COLORS, guide = "none") +
  facet_wrap(~ Metric, scales = "free_x", nrow = 2) +
  labs(title = "Metrics with 95% patient-level bootstrap CI", x = NULL, y = NULL,
       subtitle = paste0(B_BOOT, " resamples drawn at the patient level | ",
                         "overlapping intervals indicate pipelines that are ",
                         "not distinguishable")) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold")),
  "R_Bootstrap_CI_Metrics.png", 12, 7)

# --- 13.7 Confusion matrix & per-class sensitivity ---
cm_df <- do.call(rbind, lapply(PIPE_LEVELS, function(nm) {
  d <- as.data.frame(cm_list[[nm]]$table)
  names(d) <- c("Predicted","Actual","Freq"); d$Pipeline <- nm; d })) %>%
  dplyr::group_by(Pipeline, Actual) %>%
  dplyr::mutate(N_actual = sum(Freq),
                Prop = ifelse(N_actual > 0, Freq / N_actual, 0)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Pipeline = factor(Pipeline, levels = PIPE_LEVELS),
                Actual   = factor(as.character(Actual),   levels = COND_LEVELS),
                Predicted = factor(as.character(Predicted), levels = COND_LEVELS),
                Diagonal = as.character(Predicted) == as.character(Actual))
write.csv(cm_df, "P_Confusion_Matrix.csv", row.names = FALSE)

save_plot(ggplot(cm_df, aes(Actual, Predicted, fill = Prop)) +
  geom_tile(colour = "white", linewidth = 1.1) +
  geom_tile(data = subset(cm_df, Diagonal), fill = NA, colour = "grey15",
            linewidth = 0.9) +
  geom_text(aes(label = ifelse(Freq == 0, "-",
                               sprintf("%d\n%.0f%%", Freq, 100 * Prop)),
                colour = Prop > 0.5), size = 3.4, lineheight = 0.95,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_gradient(low = "#EAF2F8", high = "#1A5C8A", limits = c(0,1),
                      labels = scales::percent) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey20")) +
  scale_y_discrete(limits = rev(COND_LEVELS)) +
  facet_wrap(~ Pipeline, nrow = 2, labeller = labeller(Pipeline = setNames(
    paste0(metrics_all$Pipeline, "   (Acc = ",
           sprintf("%.3f", metrics_all$Accuracy), ")"), metrics_all$Pipeline))) +
  coord_fixed() +
  labs(title = "Confusion matrices — four pipelines",
       fill = "Proportion\nwithin column",
       subtitle = paste0("Top figure = sample count, bottom = proportion of the ",
                         "actual class (diagonal = sensitivity)"),
       x = "Actual stage", y = "Predicted stage") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"), panel.grid = element_blank()),
  "P_Confusion_Matrix.png", 10, 9.5)

sens_df <- cm_df %>% dplyr::filter(Diagonal) %>%
  dplyr::select(Pipeline, Class = Actual, Correct = Freq, N = N_actual, Prop)
save_plot(ggplot(sens_df, aes(Class, Prop, fill = Pipeline)) +
  geom_col(position = position_dodge(0.8), width = 0.72) +
  geom_text(aes(label = sprintf("%d/%d", Correct, N)),
            position = position_dodge(0.8), vjust = -0.4, size = 2.6) +
  scale_fill_manual(values = PIPE_COLORS) +
  scale_y_continuous(limits = c(0, 1.1), labels = scales::percent, expand = c(0,0)) +
  labs(title = "Sensitivity by stage", x = NULL, y = "Sensitivity (recall)",
       fill = NULL,
       subtitle = paste0("Read from the confusion-matrix diagonal; ",
                         "labels give correct / total samples per stage")) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top",
        panel.grid.major.x = element_blank()),
  "P_Sensitivity_per_Class.png", 10, 6)
write.csv(sens_df, "P_Sensitivity_per_Class.csv", row.names = FALSE)

# --- 13.8 Importance, four-way Venn, TOP_N sensitivity ---
# A composite key keeps the bar order correct within each facet.
imp_plot <- imp_all %>% dplyr::group_by(Pipeline) %>% dplyr::slice_head(n = 15) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Pipeline = factor(Pipeline, levels = PIPE_LEVELS),
                key = paste(Pipeline, Gene, sep = "___")) %>%
  dplyr::arrange(Pipeline, Mean_Importance) %>%
  dplyr::mutate(key = factor(key, levels = key))

save_plot(ggplot(imp_plot, aes(Mean_Importance, key, fill = Pipeline)) +
  geom_col(width = 0.74, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.1f", Mean_Importance)), hjust = -0.15,
            size = 2.5, colour = "grey25") +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  scale_x_continuous(limits = c(0, 118), expand = c(0,0)) +
  scale_fill_manual(values = PIPE_COLORS, guide = "none") +
  facet_wrap(~ Pipeline, scales = "free_y", nrow = 2) +
  labs(title = "Feature importance — top 15 genes per pipeline",
       x = "Importance (scaled 0-100)", y = NULL,
       subtitle = paste0("RF: mean decrease in impurity (model-based) | ",
                         "SVM: univariate ROC filter (NOT model-based)")) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 7.5)),
  "P_FeatureImportance_Grid.png", 11, 9)

save_plot(ggVennDiagram(venn_sets, label_alpha = 0, label = "count",
                     set_size = 3.6, label_size = 3.2) +
  scale_fill_gradient(low = "#F2F7FB", high = "#1A5C8A") +
  scale_x_continuous(expand = expansion(mult = 0.30)) +
  scale_y_continuous(expand = expansion(mult = 0.18)) +
  labs(title = paste0("Overlap of the top ", TOP_N,
                      " ranked genes across four pipelines"),
       subtitle = "Genes in the central region are selected by every pipeline") +
  theme(plot.title = element_text(face = "bold"), legend.position = "none",
        plot.margin = ggplot2::margin(10, 28, 10, 28)),
  "P_Venn_Biomarker_4Pipelines.png", 10, 8)

save_plot(ggplot(topn_sensitivity, aes(TOP_N, Intersection)) +
  geom_line(linewidth = 0.8, colour = "#1A5C8A") +
  geom_point(size = 2.6, colour = "#1A5C8A") +
  geom_vline(xintercept = TOP_N, linetype = "dashed", colour = "#B03A3A") +
  geom_text(aes(label = Intersection), vjust = -0.9, size = 3) +
  scale_x_continuous(breaks = topn_sensitivity$TOP_N) +
  labs(title = "Core biomarker count as a function of TOP_N",
       subtitle = paste0("Intersection of the top-N gene lists from all four ",
                         "pipelines; monotonically increasing in N. ",
                         "Dashed line marks the value used (", TOP_N, ")."),
       x = "TOP_N", y = "Genes in the four-way intersection") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()),
  "R_TOPN_Sensitivity.png", 9, 6)


# === SECTION 14: REPRODUCIBILITY & SUMMARY ===================================

writeLines(capture.output(print(sessionInfo())), "R_sessionInfo.txt")
paket <- c("DESeq2","edgeR","caret","randomForest","kernlab","e1071","pROC","ggplot2")
pkg_versions <- data.frame(Package = paket,
  Version = sapply(paket, function(p)
    tryCatch(as.character(packageVersion(p)), error = function(e) "-")),
  row.names = NULL)
write.csv(pkg_versions, "R_Package_Versions.csv", row.names = FALSE)

cat("\n=== RINGKASAN ===\n")
cat("Samples       :", N_SAMPEL, "|", length(patients), "patients\n")
cat("Genes kept    : DESeq2", nrow(dds), "| EdgeR", nrow(dge), "\n")
cat("ML features   : DESeq2", ncol(X_deseq2), "| EdgeR", ncol(X_edger), "\n")
cat("Terbaik (AUC) :", metrics_all$Pipeline[which.max(metrics_all$Macro_AUC)],
    "—", round(max(metrics_all$Macro_AUC), 4), "\n")
cat("DeLong test   :", sum(delong_test$p_adj < 0.05), "/", nrow(delong_test),
    "passed p_adj < 0,05\n")
cat("Bootstrap     :", sum(!boot_diff$Contains_zero), "/", nrow(boot_diff),
    "differences whose CI excludes zero\n")
cat("Core biomarker:", length(core_biomarker),
    if (length(core_biomarker) > 0)
      paste0(" (", paste(core_biomarker, collapse = ", "), ")") else "", "\n")
cat("\n", R.version.string, "\n")
cat("Files P_*.csv, P_*.png, R_*.csv, R_*.png written.\n")

###############################################################################
#  Limitations NOT addressed by this script, to be stated in the manuscript:
#  feature selection and hyperparameter tuning sit outside the validation
#  folds (see the separate nested-CV script); sample size and class imbalance;
#  FFPE tissue and the 3SEQ protocol; absence of clinical covariates, patient
#  outcome, external validation and laboratory confirmation; bulk RNA-seq does
#  not separate tumour epithelium from stroma.
###############################################################################
