###############################################################################
#  KOMPARASI ALGORITMA DEG DAN MACHINE LEARNING MULTI-PIPELINE
#  UNTUK IDENTIFIKASI BIOMARKER KANDIDAT KANKER PAYUDARA
#
#  Dataset : GSE47462 (RNA-Seq 3SEQ, FFPE) — Normal / EN / DCIS / IDC
#            72 sampel, 25 pasien
#  Pipeline: DESeq2 x {RF, SVM} dan EdgeR x {RF, SVM}
#  Validasi: Grouped leave-one-patient-out cross-validation
#
#  Mohammad Fadhil Ihsan (10622010) — SITH ITB
#  Pembimbing: Popi Septiani, S.Si., M.Si., Ph.D.
#
#  Urutan bagian:
#    0  Library & konfigurasi          8  Training 4 pipeline
#    1  Import data & metadata         9  Metrik evaluasi
#    2  Cabang DESeq2                 10  Uji signifikansi
#    3  Cabang EdgeR                  11  Feature importance
#    4  Perbandingan antar cabang     12  Core biomarker
#    5  Penandaan gen artefak         13  Visualisasi
#    6  Sensitivitas DEG              14  Reproduktibilitas
#    7  Set fitur & matriks ML
###############################################################################


# === BAGIAN 0: LIBRARY & KONFIGURASI =========================================

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

# Tabrakan nama: randomForest mengekspor margin() dan dimuat SETELAH ggplot2,
# sehingga margin(10, 26, 10, 26) jatuh ke randomForest:::margin.default dan
# gagal dengan "26 is not a factor". Setiap pemanggilan margin() di skrip ini
# karena itu diberi prefiks ggplot2::. Verba dplyr yang rawan tertutup paket
# Bioconductor (select, filter, count, slice_*) juga diberi prefiks dplyr::.

DATA_PATH <- "/Users/fadhilihsan/TA/GSE47462 Raw Counts Refseq Genes.txt"

PADJ_CUTOFF  <- 0.01
FC_CUTOFF    <- 2.0
TOP_N        <- 20
N_TREE       <- 500
TUNE_L       <- 3
SEED         <- 42
B_BOOT       <- 1000     # ulangan bootstrap per pasien
DOM_CUTOFF   <- 0.50     # ambang dominasi penanda gen artefak
VOLC_TRIM    <- 0.999    # pemangkasan sumbu-x volcano
JALAN_PAIRED <- TRUE     # sensitivitas desain ~ patient + condition

COND_LEVELS  <- c("Normal", "EN", "DCIS", "IDC")
PIPE_COLORS  <- c("DESeq2 + RF"  = "#B03A3A", "DESeq2 + SVM" = "#E08585",
                  "EdgeR + RF"   = "#0F3D35", "EdgeR + SVM"  = "#2AA181")
PIPE_LEVELS  <- names(PIPE_COLORS)
COND_COLORS  <- c(Normal = "#2196F3", EN = "#4CAF50",
                  DCIS = "#FF9800", IDC = "#F44336")
CLASS_COLORS <- c("#378ADD", "#3F9E52", "#EF9F27", "#E24B4A")

set.seed(SEED)

# Subsetting base R dipakai, bukan dplyr::filter, karena dplyr menghapus row
# names dan nama gen akan hilang tanpa pesan error.
ambil_deg <- function(res_table, kolom_p, kolom_fc,
                      p_cut = PADJ_CUTOFF, fc_cut = FC_CUTOFF) {
  d <- as.data.frame(res_table)
  lolos <- !is.na(d[[kolom_p]]) & d[[kolom_p]] < p_cut &
           !is.na(d[[kolom_fc]]) & abs(d[[kolom_fc]]) > fc_cut
  d <- d[lolos, , drop = FALSE]
  d <- d[order(d[[kolom_p]]), , drop = FALSE]
  d$Gene <- rownames(d)
  d
}


# === BAGIAN 1: IMPORT DATA & METADATA ========================================

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

cat("Sampel:", N_SAMPEL, "| Pasien:", length(unique(coldata$patient)),
    "| Gen mentah:", nrow(counts_raw), "\n")
print(table(coldata$condition))


# === BAGIAN 2: CABANG DESeq2 — FILTER, NORMALISASI, DEG, VST =================

dds <- DESeqDataSetFromMatrix(countData = counts_raw, colData = coldata,
                              design = ~ condition)
keep_deseq <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep_deseq, ]
dds <- DESeq(dds)

cat("\n[DESeq2] gen lolos filter:", nrow(dds), "dari", nrow(counts_raw), "\n")

res_EN_DESeq2       <- results(dds, contrast = c("condition","EN","Normal"),   alpha = 0.05)
res_DCIS_DESeq2     <- results(dds, contrast = c("condition","DCIS","Normal"), alpha = 0.05)
res_IDC_DESeq2      <- results(dds, contrast = c("condition","IDC","Normal"),  alpha = 0.05)
res_IDCvDCIS_DESeq2 <- results(dds, contrast = c("condition","IDC","DCIS"),    alpha = 0.05)

DEG_EN_DESeq2   <- ambil_deg(res_EN_DESeq2,   "padj", "log2FoldChange")
DEG_DCIS_DESeq2 <- ambil_deg(res_DCIS_DESeq2, "padj", "log2FoldChange")
DEG_IDC_DESeq2  <- ambil_deg(res_IDC_DESeq2,  "padj", "log2FoldChange")

cat("[DESeq2] DEG — EN:", nrow(DEG_EN_DESeq2),
    "| DCIS:", nrow(DEG_DCIS_DESeq2), "| IDC:", nrow(DEG_IDC_DESeq2), "\n")

write.csv(DEG_EN_DESeq2,   "P_DESeq2_DEG_EN_vs_Normal.csv",   row.names = FALSE)
write.csv(DEG_DCIS_DESeq2, "P_DESeq2_DEG_DCIS_vs_Normal.csv", row.names = FALSE)
write.csv(DEG_IDC_DESeq2,  "P_DESeq2_DEG_IDC_vs_Normal.csv",  row.names = FALSE)
write.csv(as.data.frame(res_EN_DESeq2),   "P_DESeq2_hasil_lengkap_EN.csv")
write.csv(as.data.frame(res_DCIS_DESeq2), "P_DESeq2_hasil_lengkap_DCIS.csv")
write.csv(as.data.frame(res_IDC_DESeq2),  "P_DESeq2_hasil_lengkap_IDC.csv")

vsd <- if (nrow(dds) >= 1000) vst(dds, blind = FALSE) else
       varianceStabilizingTransformation(dds, blind = FALSE)
mat_deseq2 <- assay(vsd)


# === BAGIAN 3: CABANG EdgeR — FILTER, TMM, DEG, logCPM =======================

dge <- DGEList(counts = counts_raw, group = coldata$condition)
keep_edger <- filterByExpr(dge, group = coldata$condition)
dge <- dge[keep_edger, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")

cat("\n[EdgeR] gen lolos filter:", nrow(dge), "dari", nrow(counts_raw), "\n")

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

DEG_EN_EdgeR   <- ambil_deg(res_EN_EdgeR,   "FDR", "logFC")
DEG_DCIS_EdgeR <- ambil_deg(res_DCIS_EdgeR, "FDR", "logFC")
DEG_IDC_EdgeR  <- ambil_deg(res_IDC_EdgeR,  "FDR", "logFC")

cat("[EdgeR] DEG — EN:", nrow(DEG_EN_EdgeR),
    "| DCIS:", nrow(DEG_DCIS_EdgeR), "| IDC:", nrow(DEG_IDC_EdgeR), "\n")

write.csv(DEG_EN_EdgeR,   "P_EdgeR_DEG_EN_vs_Normal.csv",   row.names = FALSE)
write.csv(DEG_DCIS_EdgeR, "P_EdgeR_DEG_DCIS_vs_Normal.csv", row.names = FALSE)
write.csv(DEG_IDC_EdgeR,  "P_EdgeR_DEG_IDC_vs_Normal.csv",  row.names = FALSE)
write.csv(res_EN_EdgeR,   "P_EdgeR_hasil_lengkap_EN.csv")
write.csv(res_DCIS_EdgeR, "P_EdgeR_hasil_lengkap_DCIS.csv")
write.csv(res_IDC_EdgeR,  "P_EdgeR_hasil_lengkap_IDC.csv")

mat_edger <- cpm(dge, log = TRUE, prior.count = 2)


# === BAGIAN 4: PERBANDINGAN DESKRIPTIF DESeq2 vs EdgeR =======================
# Deskriptif. TIDAK dipakai untuk seleksi fitur ML.

kontras  <- c("EN", "DCIS", "IDC")
posthoc  <- list(
  DESeq2 = list(EN = DEG_EN_DESeq2$Gene, DCIS = DEG_DCIS_DESeq2$Gene,
                IDC = DEG_IDC_DESeq2$Gene),
  EdgeR  = list(EN = DEG_EN_EdgeR$Gene,  DCIS = DEG_DCIS_EdgeR$Gene,
                IDC = DEG_IDC_EdgeR$Gene))

jac_deg <- do.call(rbind, lapply(kontras, function(k) {
  a <- posthoc$DESeq2[[k]]; b <- posthoc$EdgeR[[k]]
  data.frame(Kontras = paste(k, "vs Normal"),
             DESeq2 = length(a), EdgeR = length(b),
             Irisan = length(intersect(a, b)), Gabungan = length(union(a, b)),
             Khas_DESeq2 = length(setdiff(a, b)),
             Khas_EdgeR  = length(setdiff(b, a)),
             Jaccard = round(length(intersect(a, b)) / length(union(a, b)), 4))
}))
cat("\n=== IRISAN DEG DESeq2 vs EdgeR ===\n"); print(jac_deg, row.names = FALSE)
write.csv(jac_deg, "P_Overlap_DEG_Jaccard.csv", row.names = FALSE)

gradien_deg <- data.frame(
  Tahap  = kontras,
  DESeq2 = sapply(kontras, function(k) length(posthoc$DESeq2[[k]])),
  EdgeR  = sapply(kontras, function(k) length(posthoc$EdgeR[[k]])),
  row.names = NULL)
write.csv(gradien_deg, "P_Gradien_DEG_per_Tahap.csv", row.names = FALSE)

ringkasan_cabang <- data.frame(
  Tahap  = c("Gen setelah filter", "Normalisasi", "Uji statistik",
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
print(ringkasan_cabang, row.names = FALSE)
write.csv(ringkasan_cabang, "P_Ringkasan_Dua_Cabang.csv", row.names = FALSE)


# === BAGIAN 5: PENANDAAN GEN ARTEFAK =========================================
# Dominasi = proporsi total count sebuah gen yang berasal dari satu sampel.
# Gen ditandai, tidak dibuang.

tot   <- rowSums(counts_raw)
domn  <- ifelse(tot > 0, apply(counts_raw, 1, max) / pmax(tot, 1), NA_real_)
n_det <- rowSums(counts_raw >= 5)

deg_semua <- unique(unlist(c(posthoc$DESeq2, posthoc$EdgeR)))

flag_artefak <- data.frame(Gene = rownames(counts_raw), Total_count = tot,
                           Dominasi = round(domn, 4),
                           N_sampel_terdeteksi = n_det, row.names = NULL) %>%
  dplyr::filter(Gene %in% deg_semua) %>%
  dplyr::mutate(Tersangka = Dominasi > DOM_CUTOFF | N_sampel_terdeteksi < 5) %>%
  dplyr::arrange(dplyr::desc(Dominasi))

cat("\n=== GEN ARTEFAK ===\nDEG unik:", length(deg_semua),
    "| ditandai tersangka:", sum(flag_artefak$Tersangka), "\n")
write.csv(flag_artefak, "R_Gen_Artefak_Flag.csv", row.names = FALSE)


# === BAGIAN 6: SENSITIVITAS DEG ==============================================
# Analisis pendamping, bukan pengganti. Hasil utama tetap dari BAGIAN 2-3.

# --- 6.1 Uji efek-size native (ambang masuk ke hipotesis nol) ---
deg_native <- list(
  DESeq2 = setNames(lapply(kontras, function(k) {
    r <- results(dds, contrast = c("condition", k, "Normal"),
                 lfcThreshold = FC_CUTOFF, altHypothesis = "greaterAbs",
                 alpha = PADJ_CUTOFF)
    rownames(r)[!is.na(r$padj) & r$padj < PADJ_CUTOFF]
  }), kontras),
  EdgeR = setNames(lapply(kontras, function(k) {
    tt <- topTags(glmTreat(fit, contrast = contrasts_matrix[, paste0(k, "_vs_Normal")],
                           lfc = FC_CUTOFF), n = Inf)$table
    rownames(tt)[tt$FDR < PADJ_CUTOFF]
  }), kontras))

banding_uji <- do.call(rbind, lapply(names(posthoc), function(cb)
  do.call(rbind, lapply(kontras, function(k) {
    a <- posthoc[[cb]][[k]]; b <- deg_native[[cb]][[k]]
    data.frame(Cabang = cb, Kontras = k, Post_hoc = length(a), Native = length(b),
               Irisan = length(intersect(a, b)),
               Hanya_post_hoc = length(setdiff(a, b)),
               Hanya_native   = length(setdiff(b, a)))
  }))))
cat("\n=== UJI NATIVE vs POST-HOC ===\n"); print(banding_uji, row.names = FALSE)
write.csv(banding_uji, "R_Uji_Native_vs_PostHoc.csv", row.names = FALSE)

# --- 6.2 Desain berpasangan ~ patient + condition ---
if (JALAN_PAIRED) {
  cd <- as.data.frame(SummarizedExperiment::colData(dds))
  cd$patient <- factor(cd$patient)
  mm <- model.matrix(~ patient + condition, data = cd)
  penuh <- qr(mm)$rank == ncol(mm)
  simpan_kol <- seq_len(ncol(dds))

  if (!penuh) {
    n_pp   <- table(cd$patient)
    buang  <- names(n_pp)[n_pp < 2]
    simpan_kol <- which(!(cd$patient %in% buang))
    mm <- model.matrix(~ patient + condition, data = droplevels(cd[simpan_kol, ]))
    penuh <- qr(mm)$rank == ncol(mm)
    cat("\n[paired] membuang pasien bersampel tunggal:",
        paste(buang, collapse = ", "), "\n")
  }
  cat("[paired] desain", nrow(mm), "x", ncol(mm), "| rank =", qr(mm)$rank,
      "|", ifelse(penuh, "full rank", "rank-deficient"), "\n")

  if (penuh) {
    dds_p <- dds[, simpan_kol]
    SummarizedExperiment::colData(dds_p)$patient <-
      droplevels(factor(SummarizedExperiment::colData(dds_p)$patient))
    design(dds_p) <- ~ patient + condition
    dds_p <- DESeq(dds_p, quiet = TRUE)

    deg_paired <- setNames(lapply(kontras, function(k)
      ambil_deg(results(dds_p, contrast = c("condition", k, "Normal"),
                        alpha = PADJ_CUTOFF), "padj", "log2FoldChange")$Gene),
      kontras)

    banding_paired <- do.call(rbind, lapply(kontras, function(k) {
      a <- posthoc$DESeq2[[k]]; b <- deg_paired[[k]]
      data.frame(Kontras = k, Tanpa_pasien = length(a), Dengan_pasien = length(b),
                 Irisan = length(intersect(a, b)))
    }))
    cat("\n=== DESAIN BERPASANGAN (DESeq2) ===\n")
    print(banding_paired, row.names = FALSE)
    write.csv(banding_paired, "R_Desain_Berpasangan_DESeq2.csv", row.names = FALSE)
  }
}


# === BAGIAN 7: SET FITUR & MATRIKS ML ========================================
# Union DEG tiga kontras terhadap Normal, per cabang. Tanpa konsensus.

genes_deseq2 <- intersect(unique(unlist(posthoc$DESeq2)), rownames(mat_deseq2))
genes_edger  <- intersect(unique(unlist(posthoc$EdgeR)),  rownames(mat_edger))

if (length(genes_deseq2) < 2 || length(genes_edger) < 2)
  stop("Set fitur terlalu sedikit. Longgarkan PADJ_CUTOFF / FC_CUTOFF.")

X_deseq2 <- t(mat_deseq2[genes_deseq2, rownames(coldata), drop = FALSE])
X_edger  <- t(mat_edger[genes_edger,   rownames(coldata), drop = FALSE])
stopifnot(identical(rownames(X_deseq2), rownames(coldata)),
          identical(rownames(X_edger),  rownames(coldata)))

cat("\nFitur ML — DESeq2:", ncol(X_deseq2), "| EdgeR:", ncol(X_edger),
    "| irisan:", length(intersect(genes_deseq2, genes_edger)), "\n")

write.csv(data.frame(
  Cabang = c("DESeq2","EdgeR"), Matriks = c("VST","logCPM (TMM)"),
  Gen_Lolos_Filter = c(nrow(dds), nrow(dge)),
  Fitur_ML = c(ncol(X_deseq2), ncol(X_edger))),
  "P_Komposisi_Fitur.csv", row.names = FALSE)


# === BAGIAN 8: TRAINING 4 PIPELINE — GROUPED LOOCV ===========================
# Lipatan dikelompokkan per PASIEN: seluruh sampel satu pasien ditahan bersama.

pasien        <- unique(coldata$patient)
train_indices <- lapply(pasien, function(p) which(coldata$patient != p))
test_indices  <- lapply(pasien, function(p) which(coldata$patient == p))
names(train_indices) <- names(test_indices) <- paste0("Patient_", pasien)

ctrl_loocv <- trainControl(method = "cv", index = train_indices,
                           indexOut = test_indices, classProbs = TRUE,
                           savePredictions = "final")

latih_model <- function(X, y, algoritma, nama) {
  cat("  -", nama, "... ")
  w <- system.time({
    set.seed(SEED)
    m <- if (algoritma == "rf")
      train(x = X, y = y, method = "rf", trControl = ctrl_loocv,
            ntree = N_TREE, tuneLength = TUNE_L)
    else
      train(x = X, y = y, method = "svmRadial", trControl = ctrl_loocv,
            preProcess = c("center","scale"), tuneLength = TUNE_L)
  })
  cat("selesai (", round(w["elapsed"], 1), "detik )\n"); m
}

cat("\n=== TRAINING 4 PIPELINE (", length(pasien), "lipatan pasien ) ===\n")
models <- list(
  "DESeq2 + RF"  = latih_model(X_deseq2, y, "rf",  "DESeq2 + RF"),
  "DESeq2 + SVM" = latih_model(X_deseq2, y, "svm", "DESeq2 + SVM"),
  "EdgeR + RF"   = latih_model(X_edger,  y, "rf",  "EdgeR + RF"),
  "EdgeR + SVM"  = latih_model(X_edger,  y, "svm", "EdgeR + SVM"))
preds <- lapply(models, function(m) m$pred)


# === BAGIAN 9: METRIK EVALUASI ===============================================

auc_per_kelas <- function(d) sapply(classes, function(k)
  as.numeric(auc(roc(as.numeric(d$obs == k), d[[k]],
                     quiet = TRUE, direction = "<"))))

cm_list <- lapply(preds, function(p) confusionMatrix(p$pred, p$obs))

metrik_all <- do.call(rbind, lapply(PIPE_LEVELS, function(nm) {
  cm <- cm_list[[nm]]; a <- auc_per_kelas(preds[[nm]])
  data.frame(Pipeline = nm,
             Metode_DEG = ifelse(grepl("DESeq2", nm), "DESeq2", "EdgeR"),
             Model_ML   = ifelse(grepl("RF", nm), "Random Forest", "SVM"),
             Jumlah_Fitur = ifelse(grepl("DESeq2", nm), ncol(X_deseq2), ncol(X_edger)),
             Accuracy = as.numeric(cm$overall["Accuracy"]),
             Kappa    = as.numeric(cm$overall["Kappa"]),
             Balanced_Accuracy = mean(cm$byClass[, "Balanced Accuracy"], na.rm = TRUE),
             Macro_F1 = mean(cm$byClass[, "F1"], na.rm = TRUE),
             Macro_AUC = mean(a), stringsAsFactors = FALSE)
}))
rownames(metrik_all) <- NULL

num_cols <- c("Accuracy","Kappa","Balanced_Accuracy","Macro_F1","Macro_AUC")
metrik_print <- metrik_all; metrik_print[num_cols] <- round(metrik_print[num_cols], 3)

cat("\n=== METRIK 4 PIPELINE ===\n")
print(metrik_print[, c("Pipeline","Jumlah_Fitur", num_cols)], row.names = FALSE)
for (nm in PIPE_LEVELS) { cat("\n>>", nm, "\n"); print(cm_list[[nm]]$table) }

auc_detail <- do.call(rbind, lapply(PIPE_LEVELS, function(nm) {
  a <- auc_per_kelas(preds[[nm]])
  data.frame(Pipeline = nm, t(round(a, 3)), Macro = round(mean(a), 3),
             check.names = FALSE)
}))
rownames(auc_detail) <- NULL
print(auc_detail, row.names = FALSE)

write.csv(metrik_print, "P_Metrik_4Pipeline_LOOCV.csv", row.names = FALSE)
write.csv(auc_detail,   "P_AUC_per_Kelas.csv",          row.names = FALSE)


# === BAGIAN 10: UJI SIGNIFIKANSI =============================================
# caret mengurutkan m$pred per resample, bukan per sampel. Penyelarasan lewat
# rowIndex wajib sebelum paired = TRUE.

pr <- lapply(preds, function(d) d[order(d$rowIndex), , drop = FALSE])
stopifnot(length(unique(lapply(pr, function(d) d$rowIndex))) == 1)
pat <- coldata$patient[pr[[1]]$rowIndex]

roc_obj <- function(d, k) {
  o <- as.numeric(d$obs == k)
  if (length(unique(o)) < 2) return(NULL)
  roc(o, d[[k]], quiet = TRUE, direction = "<")
}

# --- 10.1 Selang kepercayaan DeLong ---
auc_ci <- do.call(rbind, lapply(PIPE_LEVELS, function(nm)
  do.call(rbind, lapply(classes, function(k) {
    r <- roc_obj(pr[[nm]], k); if (is.null(r)) return(NULL)
    ci <- as.numeric(ci.auc(r, method = "delong"))
    data.frame(Pipeline = nm, Kelas = k, AUC = ci[2],
               CI_bawah = ci[1], CI_atas = ci[3])
  }))))
cat("\n=== AUC + 95% CI DeLong ===\n")
print(as.data.frame(auc_ci %>% dplyr::mutate(
  dplyr::across(dplyr::where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
write.csv(auc_ci, "R_AUC_CI_DeLong.csv", row.names = FALSE)

# --- 10.2 Uji DeLong berpasangan antar pipeline ---
pasangan <- t(combn(PIPE_LEVELS, 2))
uji_del <- do.call(rbind, lapply(seq_len(nrow(pasangan)), function(i)
  do.call(rbind, lapply(classes, function(k) {
    r1 <- roc_obj(pr[[pasangan[i,1]]], k); r2 <- roc_obj(pr[[pasangan[i,2]]], k)
    if (is.null(r1) || is.null(r2)) return(NULL)
    tt <- try(roc.test(r1, r2, method = "delong", paired = TRUE), silent = TRUE)
    if (inherits(tt, "try-error")) return(NULL)
    data.frame(Kelas = k, Pipeline_A = pasangan[i,1], Pipeline_B = pasangan[i,2],
               AUC_A = as.numeric(auc(r1)), AUC_B = as.numeric(auc(r2)),
               Selisih = as.numeric(auc(r1)) - as.numeric(auc(r2)), p = tt$p.value)
  }))))
uji_del$p_adj <- p.adjust(uji_del$p, method = "BH")
cat("\n=== UJI DeLong BERPASANGAN ===\nLolos p_adj < 0,05:",
    sum(uji_del$p_adj < 0.05), "dari", nrow(uji_del), "uji\n")
print(as.data.frame(uji_del %>% dplyr::arrange(p) %>% dplyr::slice_head(n = 8) %>%
  dplyr::mutate(Selisih = round(Selisih, 3), p = signif(p, 3),
                p_adj = signif(p_adj, 3)) %>%
  dplyr::select(Kelas, Pipeline_A, Pipeline_B, Selisih, p, p_adj)),
  row.names = FALSE)
write.csv(uji_del, "R_RocTest_DeLong.csv", row.names = FALSE)

# --- 10.3 Bootstrap per pasien (tanpa melatih ulang) ---
metrik_dari_idx <- function(d, idx) {
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

set.seed(SEED); up <- unique(pat)
cat("\n=== BOOTSTRAP PER PASIEN (B =", B_BOOT, ") ===\nBerjalan ")
bt <- lapply(seq_len(B_BOOT), function(b) {
  if (b %% 100 == 0) cat(".")
  idx <- unlist(lapply(sample(up, length(up), replace = TRUE),
                       function(p) which(pat == p)))
  sapply(PIPE_LEVELS, function(nm)
    tryCatch(metrik_dari_idx(pr[[nm]], idx),
             error = function(e) rep(NA_real_, 5)))
})
cat(" selesai\n")

boot_ci <- do.call(rbind, lapply(PIPE_LEVELS, function(nm)
  do.call(rbind, lapply(seq_along(num_cols), function(j) {
    v  <- sapply(bt, function(m) m[j, nm])
    ci <- quantile(v, c(0.025, 0.975), na.rm = TRUE)
    data.frame(Pipeline = nm, Metrik = num_cols[j],
               Titik = as.numeric(metrik_all[metrik_all$Pipeline == nm, num_cols[j]]),
               CI_bawah = ci[1], CI_atas = ci[2])
  }))))
print(as.data.frame(boot_ci %>% dplyr::mutate(
  dplyr::across(dplyr::where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
write.csv(boot_ci, "R_Bootstrap_CI_Metrik.csv", row.names = FALSE)

boot_sel <- do.call(rbind, lapply(seq_len(nrow(pasangan)), function(i)
  do.call(rbind, lapply(seq_along(num_cols), function(j) {
    d  <- sapply(bt, function(m) m[j, pasangan[i,1]] - m[j, pasangan[i,2]])
    ci <- quantile(d, c(0.025, 0.975), na.rm = TRUE)
    data.frame(Metrik = num_cols[j], Pipeline_A = pasangan[i,1],
               Pipeline_B = pasangan[i,2],
               Selisih_median = median(d, na.rm = TRUE),
               CI_bawah = ci[1], CI_atas = ci[2],
               Memuat_nol = ci[1] <= 0 && ci[2] >= 0)
  }))))
cat("Selisih pipeline dengan CI tidak memuat nol:",
    sum(!boot_sel$Memuat_nol), "dari", nrow(boot_sel), "\n")
write.csv(boot_sel, "R_Bootstrap_Selisih_Pipeline.csv", row.names = FALSE)


# === BAGIAN 11: FEATURE IMPORTANCE ===========================================
# varImp untuk svmRadial BUKAN model-based: caret jatuh ke filter ROC univariat.
# Kolom Metode_Importance mencegah salah tafsir saat tabel dibaca ulang.

imp_penuh <- setNames(lapply(PIPE_LEVELS, function(nm) {
  im <- varImp(models[[nm]], scale = TRUE)$importance
  im$Gene <- rownames(im)
  sc <- setdiff(colnames(im), "Gene")
  im$Mean_Importance <- rowMeans(im[, sc, drop = FALSE], na.rm = TRUE)
  im[order(-im$Mean_Importance), c("Gene","Mean_Importance")]
}), PIPE_LEVELS)

imp_list <- setNames(lapply(PIPE_LEVELS, function(nm)
  head(imp_penuh[[nm]], TOP_N) %>%
    dplyr::mutate(Pipeline = nm, Peringkat = dplyr::row_number(),
                  Metode_Importance = ifelse(grepl("RF", nm),
                    "Mean decrease in impurity (model-based)",
                    "Filter ROC univariat (BUKAN model-based)"))),
  PIPE_LEVELS)
imp_all <- do.call(rbind, imp_list); rownames(imp_all) <- NULL

cat("\n=== TOP 10 GEN PER PIPELINE ===\n")
for (nm in PIPE_LEVELS) {
  cat("\n>>", nm, "\n")
  print(as.data.frame(imp_list[[nm]][1:10, c("Peringkat","Gene","Mean_Importance")]),
        row.names = FALSE)
}
write.csv(imp_all, "P_Feature_Importance.csv", row.names = FALSE)

# Sensitivitas irisan terhadap TOP_N
sens_topn <- do.call(rbind, lapply(c(5,10,15,20,25,30,40,50,75,100), function(n) {
  ir <- Reduce(intersect, lapply(PIPE_LEVELS, function(nm) head(imp_penuh[[nm]]$Gene, n)))
  data.frame(TOP_N = n, Irisan = length(ir),
             Gen = ifelse(length(ir) == 0, "-", paste(head(ir, 12), collapse = ", ")))
}))
cat("\n=== SENSITIVITAS TOP_N ===\n"); print(sens_topn, row.names = FALSE)
write.csv(sens_topn, "R_Sensitivitas_TOP_N.csv", row.names = FALSE)


# === BAGIAN 12: CORE BIOMARKER ===============================================
# Irisan dibaca SETELAH keempat model selesai dilatih. Bukan seleksi fitur.

venn_sets <- setNames(lapply(imp_list, function(d) d$Gene), PIPE_LEVELS)
core_biomarker <- Reduce(intersect, venn_sets)

cat("\n=== CORE BIOMARKER (irisan Top-", TOP_N, ") ===\n", sep = "")
if (length(core_biomarker) > 0) {
  cat("Jumlah:", length(core_biomarker), "gen —",
      paste(core_biomarker, collapse = ", "), "\n")
  core_table <- imp_all %>%
    dplyr::filter(Gene %in% core_biomarker) %>%
    dplyr::select(Gene, Pipeline, Peringkat) %>%
    tidyr::pivot_wider(names_from = Pipeline, values_from = Peringkat) %>%
    dplyr::mutate(Rerata_Peringkat = rowMeans(dplyr::select(., -Gene), na.rm = TRUE)) %>%
    dplyr::arrange(Rerata_Peringkat)
  print(as.data.frame(core_table), row.names = FALSE)
  write.csv(core_table, "P_Core_Biomarker.csv", row.names = FALSE)
} else {
  cat("Irisan kosong pada TOP_N =", TOP_N, "\n")
}

overlap_matrix <- outer(PIPE_LEVELS, PIPE_LEVELS, Vectorize(function(a, b)
  length(intersect(venn_sets[[a]], venn_sets[[b]]))))
dimnames(overlap_matrix) <- list(PIPE_LEVELS, PIPE_LEVELS)
print(overlap_matrix)
write.csv(as.data.frame(overlap_matrix), "P_Overlap_Matrix.csv")


# === BAGIAN 13: VISUALISASI ==================================================

simpan <- function(p, f, w, h) ggsave(f, p, width = w, height = h, dpi = 300,
                                      bg = "white")

# --- 13.1 PCA (tanpa anotasi sampel) ---
pca_d  <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
var_d  <- round(100 * attr(pca_d, "percentVar"), 1)
simpan(ggplot(pca_d, aes(PC1, PC2, colour = condition)) +
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
simpan(ggplot(data.frame(PC1 = pca_e$x[,1], PC2 = pca_e$x[,2],
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

# --- 13.2 Volcano (sumbu identik di keenam panel) ---
volc_semua <- list(list(res_EN_DESeq2,"log2FoldChange","padj"),
                   list(res_DCIS_DESeq2,"log2FoldChange","padj"),
                   list(res_IDC_DESeq2,"log2FoldChange","padj"),
                   list(res_EN_EdgeR,"logFC","FDR"),
                   list(res_DCIS_EdgeR,"logFC","FDR"),
                   list(res_IDC_EdgeR,"logFC","FDR"))
.fc <- unlist(lapply(volc_semua, function(v) as.data.frame(v[[1]])[[v[[2]]]]))
.pv <- unlist(lapply(volc_semua, function(v) as.data.frame(v[[1]])[[v[[3]]]]))
.fc <- .fc[is.finite(.fc)]; .pv <- .pv[is.finite(.pv) & .pv > 0]
VOLC_XMAX <- ceiling(as.numeric(quantile(abs(.fc), VOLC_TRIM, na.rm = TRUE)))
VOLC_XLIM <- c(-VOLC_XMAX, VOLC_XMAX)
VOLC_YLIM <- c(0, ceiling(max(-log10(.pv), na.rm = TRUE) * 1.08))
cat("\nSumbu volcano:", paste(VOLC_XLIM, collapse = " s/d "),
    "| y 0 s/d", VOLC_YLIM[2], "|", sum(abs(.fc) > VOLC_XMAX), "titik di luar\n")

buat_volcano <- function(res_df, kfc, kp, judul, file) {
  d <- as.data.frame(res_df)
  nu <- sum(d[[kfc]] >  FC_CUTOFF & d[[kp]] < PADJ_CUTOFF, na.rm = TRUE)
  nd <- sum(d[[kfc]] < -FC_CUTOFF & d[[kp]] < PADJ_CUTOFF, na.rm = TRUE)
  png(file, width = 2400, height = 2000, res = 300)
  print(EnhancedVolcano(d, lab = rownames(d), x = kfc, y = kp, title = judul,
    subtitle = paste0(kp," < ",PADJ_CUTOFF,"  |  |",kfc,"| > ",FC_CUTOFF,
                      "  |  Up: ",nu,"  |  Down: ",nd),
    xlim = VOLC_XLIM, ylim = VOLC_YLIM, pCutoff = PADJ_CUTOFF,
    FCcutoff = FC_CUTOFF, pointSize = 1.5, labSize = 3,
    col = c("grey70","#2196F3","#FF9800","#F44336"), colAlpha = 0.7,
    legendLabels = c("NS","log2FC","p-adj","p-adj & log2FC"),
    drawConnectors = TRUE))
  dev.off()
}
buat_volcano(res_EN_DESeq2,  "log2FoldChange","padj","EN vs Normal (DESeq2)",  "P_DESeq2_Volcano_EN.png")
buat_volcano(res_DCIS_DESeq2,"log2FoldChange","padj","DCIS vs Normal (DESeq2)","P_DESeq2_Volcano_DCIS.png")
buat_volcano(res_IDC_DESeq2, "log2FoldChange","padj","IDC vs Normal (DESeq2)", "P_DESeq2_Volcano_IDC.png")
buat_volcano(res_EN_EdgeR,   "logFC","FDR","EN vs Normal (EdgeR)",   "P_EdgeR_Volcano_EN.png")
buat_volcano(res_DCIS_EdgeR, "logFC","FDR","DCIS vs Normal (EdgeR)", "P_EdgeR_Volcano_DCIS.png")
buat_volcano(res_IDC_EdgeR,  "logFC","FDR","IDC vs Normal (EdgeR)",  "P_EdgeR_Volcano_IDC.png")

# --- 13.3 Heatmap Top 50 DEG IDC ---
anot_col <- data.frame(Stage = coldata$condition, row.names = rownames(coldata))
buat_heatmap <- function(mat, gen, judul, file) {
  g <- intersect(gen, rownames(mat)); if (length(g) < 2) return(invisible())
  png(file, width = 2600, height = 3600, res = 300)
  pheatmap(mat[g, ], annotation_col = anot_col,
           annotation_colors = list(Stage = COND_COLORS), scale = "row",
           clustering_method = "ward.D2", show_colnames = FALSE,
           fontsize_row = 7, main = judul,
           color = colorRampPalette(c("#2166AC","white","#D6604D"))(100))
  dev.off()
}
buat_heatmap(mat_deseq2, head(DEG_IDC_DESeq2$Gene, 50),
             "Top 50 DEG IDC vs Normal — DESeq2 (VST, row z-score)",
             "P_DESeq2_Heatmap_Top50.png")
buat_heatmap(mat_edger, head(DEG_IDC_EdgeR$Gene, 50),
             "Top 50 DEG IDC vs Normal — EdgeR (logCPM, row z-score)",
             "P_EdgeR_Heatmap_Top50.png")

# --- 13.4 Venn DEG per kontras ---
# ggVennDiagram menaruh label himpunan DI LUAR lingkaran, dan panel ggplot
# memotongnya sehingga "DESeq2" terbaca "ESeq2" dan "EdgeR" terbaca "dgeR".
# Perbaikannya melebarkan sumbu supaya label masuk ke dalam panel, ditambah
# margin plot. Jangan menambahkan coord_sf() di sini: ia akan MENGGANTI sistem
# koordinat yang sudah dipasang ggVennDiagram. Kalau label masih terpotong,
# naikkan mult pada expansion(), bukan memakai coord_sf().
buat_venn2 <- function(a, b, judul, lo, hi, file) {
  if (length(a) == 0 || length(b) == 0) return(invisible())
  simpan(ggVennDiagram(list(DESeq2 = a, EdgeR = b), label_alpha = 0,
                       label = "count", set_size = 4.2, label_size = 3.8) +
    scale_fill_gradient(low = lo, high = hi) +
    scale_x_continuous(expand = expansion(mult = 0.30)) +
    scale_y_continuous(expand = expansion(mult = 0.18)) +
    labs(title = judul,
         subtitle = "Descriptive only — not used for ML feature selection") +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"),
          plot.margin = ggplot2::margin(10, 26, 10, 26)),
    file, 8, 6.5)
}
buat_venn2(DEG_EN_DESeq2$Gene,  DEG_EN_EdgeR$Gene,
           "DEG EN vs Normal — DESeq2 vs EdgeR",  "#E8F5E9","#2E7D32","P_Venn_DEG_EN.png")
buat_venn2(DEG_DCIS_DESeq2$Gene,DEG_DCIS_EdgeR$Gene,
           "DEG DCIS vs Normal — DESeq2 vs EdgeR","#FFF3E0","#E65100","P_Venn_DEG_DCIS.png")
buat_venn2(DEG_IDC_DESeq2$Gene, DEG_IDC_EdgeR$Gene,
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

simpan(ggplot(roc_all, aes(FPR, TPR, colour = Pipeline)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey55") +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = PIPE_COLORS,
    labels = setNames(paste0(metrik_all$Pipeline, "  (AUC = ",
                             sprintf("%.3f", metrik_all$Macro_AUC), ")"),
                      metrik_all$Pipeline)) +
  coord_equal(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
  labs(title = "ROC curves — four pipelines, grouped LOOCV", colour = NULL,
       subtitle = paste0("Macro-average one-vs-rest | ", N_SAMPEL,
                         " samples from ", length(pasien), " patients"),
       x = "False positive rate (1 - specificity)",
       y = "True positive rate (sensitivity)") +
  theme_bw(base_size = 12) + theme(plot.title = element_text(face = "bold")),
  "P_ROC_Gabungan.png", 10, 6.5)

roc_kelas <- do.call(rbind, lapply(PIPE_LEVELS, function(nm)
  do.call(rbind, lapply(classes, function(k) {
    r <- roc_obj(preds[[nm]], k); if (is.null(r)) return(NULL)
    data.frame(FPR = 1 - r$specificities, TPR = r$sensitivities, Kelas = k,
               AUC = as.numeric(auc(r)), Pipeline = nm) }))))
roc_kelas <- roc_kelas[order(roc_kelas$Pipeline, roc_kelas$Kelas, roc_kelas$FPR), ]
roc_kelas$Pipeline <- factor(roc_kelas$Pipeline, levels = PIPE_LEVELS)
roc_kelas$Kelas    <- factor(roc_kelas$Kelas,    levels = classes)

anot_auc <- roc_kelas %>% dplyr::distinct(Pipeline, Kelas, AUC) %>%
  dplyr::arrange(Pipeline, Kelas) %>% dplyr::group_by(Pipeline) %>%
  dplyr::mutate(ypos = 0.30 - 0.075 * (dplyr::row_number() - 1)) %>% dplyr::ungroup()

simpan(ggplot(roc_kelas, aes(FPR, TPR, colour = Kelas)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey60") +
  geom_line(linewidth = 0.95) +
  geom_text(data = anot_auc, inherit.aes = FALSE,
            aes(0.98, ypos, colour = Kelas,
                label = sprintf("%s: AUC = %.3f", Kelas, AUC)),
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
  "P_ROC_per_Kelas.png", 10, 9.5)

# --- 13.6 Metrik: barplot, heatmap, CI bootstrap ---
# Nama kolom dipetakan ke label tampilan agar gambar tidak memuat garis bawah
# gaya nama variabel. CSV tetap memakai nama kolom asli supaya mesin-terbaca.
MET_LABEL <- c(Accuracy = "Accuracy", Kappa = "Cohen's kappa",
               Balanced_Accuracy = "Balanced accuracy",
               Macro_F1 = "Macro F1", Macro_AUC = "Macro AUC")

met_long <- metrik_all %>%
  dplyr::select(Pipeline, dplyr::all_of(num_cols)) %>%
  tidyr::pivot_longer(-Pipeline, names_to = "Metric", values_to = "Value") %>%
  dplyr::mutate(Pipeline = factor(Pipeline, levels = PIPE_LEVELS),
                Metric   = factor(MET_LABEL[Metric], levels = MET_LABEL))

simpan(ggplot(met_long, aes(Metric, Value, fill = Pipeline)) +
  geom_col(position = position_dodge(0.8), width = 0.72) +
  geom_text(aes(label = sprintf("%.3f", Value)), position = position_dodge(0.8),
            vjust = -0.4, size = 2.7) +
  scale_fill_manual(values = PIPE_COLORS) +
  scale_y_continuous(limits = c(0, 1.08), expand = c(0,0)) +
  labs(title = "Evaluation metrics across four pipelines", x = NULL,
       y = "Value", fill = "Pipeline",
       subtitle = paste0("Grouped LOOCV over ", length(pasien), " patients")) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top",
        panel.grid.major.x = element_blank()),
  "P_Barplot_Metrik.png", 11, 6)

simpan(met_long %>% dplyr::group_by(Metric) %>%
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
  "P_Heatmap_Metrik.png", 9, 5)

simpan(boot_ci %>%
  dplyr::mutate(Pipeline = factor(Pipeline, levels = rev(PIPE_LEVELS)),
                Metrik   = factor(MET_LABEL[Metrik], levels = MET_LABEL)) %>%
  ggplot(aes(Titik, Pipeline, colour = Pipeline)) +
  geom_linerange(aes(xmin = CI_bawah, xmax = CI_atas), linewidth = 0.8) +
  geom_point(size = 2.7) +
  geom_text(aes(label = sprintf("%.3f [%.3f, %.3f]", Titik, CI_bawah, CI_atas)),
            vjust = -1.05, size = 2.3, colour = "grey20") +
  scale_colour_manual(values = PIPE_COLORS, guide = "none") +
  facet_wrap(~ Metrik, scales = "free_x", nrow = 2) +
  labs(title = "Metrics with 95% patient-level bootstrap CI", x = NULL, y = NULL,
       subtitle = paste0(B_BOOT, " resamples drawn at the patient level | ",
                         "overlapping intervals indicate pipelines that are ",
                         "not distinguishable")) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold")),
  "R_Bootstrap_CI_Metrik.png", 12, 7)

# --- 13.7 Confusion matrix & sensitivitas per kelas ---
cm_df <- do.call(rbind, lapply(PIPE_LEVELS, function(nm) {
  d <- as.data.frame(cm_list[[nm]]$table)
  names(d) <- c("Prediksi","Aktual","Freq"); d$Pipeline <- nm; d })) %>%
  dplyr::group_by(Pipeline, Aktual) %>%
  dplyr::mutate(N_aktual = sum(Freq),
                Prop = ifelse(N_aktual > 0, Freq / N_aktual, 0)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Pipeline = factor(Pipeline, levels = PIPE_LEVELS),
                Aktual   = factor(as.character(Aktual),   levels = COND_LEVELS),
                Prediksi = factor(as.character(Prediksi), levels = COND_LEVELS),
                Diagonal = as.character(Prediksi) == as.character(Aktual))
write.csv(cm_df, "P_Confusion_Matrix.csv", row.names = FALSE)

simpan(ggplot(cm_df, aes(Aktual, Prediksi, fill = Prop)) +
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
    paste0(metrik_all$Pipeline, "   (Acc = ",
           sprintf("%.3f", metrik_all$Accuracy), ")"), metrik_all$Pipeline))) +
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
  dplyr::select(Pipeline, Kelas = Aktual, Benar = Freq, N = N_aktual, Prop)
simpan(ggplot(sens_df, aes(Kelas, Prop, fill = Pipeline)) +
  geom_col(position = position_dodge(0.8), width = 0.72) +
  geom_text(aes(label = sprintf("%d/%d", Benar, N)),
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
  "P_Sensitivitas_per_Kelas.png", 10, 6)
write.csv(sens_df, "P_Sensitivitas_per_Kelas.csv", row.names = FALSE)

# --- 13.8 Importance, Venn 4 pipeline, sensitivitas TOP_N ---
# Kunci komposit dipakai supaya urutan batang benar di tiap faset.
imp_plot <- imp_all %>% dplyr::group_by(Pipeline) %>% dplyr::slice_head(n = 15) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Pipeline = factor(Pipeline, levels = PIPE_LEVELS),
                key = paste(Pipeline, Gene, sep = "___")) %>%
  dplyr::arrange(Pipeline, Mean_Importance) %>%
  dplyr::mutate(key = factor(key, levels = key))

simpan(ggplot(imp_plot, aes(Mean_Importance, key, fill = Pipeline)) +
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

simpan(ggVennDiagram(venn_sets, label_alpha = 0, label = "count",
                     set_size = 3.6, label_size = 3.2) +
  scale_fill_gradient(low = "#F2F7FB", high = "#1A5C8A") +
  scale_x_continuous(expand = expansion(mult = 0.30)) +
  scale_y_continuous(expand = expansion(mult = 0.18)) +
  labs(title = paste0("Overlap of the top ", TOP_N,
                      " ranked genes across four pipelines"),
       subtitle = "Genes in the central region are selected by every pipeline") +
  theme(plot.title = element_text(face = "bold"), legend.position = "none",
        plot.margin = ggplot2::margin(10, 28, 10, 28)),
  "P_Venn_Biomarker_4Pipeline.png", 10, 8)

simpan(ggplot(sens_topn, aes(TOP_N, Irisan)) +
  geom_line(linewidth = 0.8, colour = "#1A5C8A") +
  geom_point(size = 2.6, colour = "#1A5C8A") +
  geom_vline(xintercept = TOP_N, linetype = "dashed", colour = "#B03A3A") +
  geom_text(aes(label = Irisan), vjust = -0.9, size = 3) +
  scale_x_continuous(breaks = sens_topn$TOP_N) +
  labs(title = "Core biomarker count as a function of TOP_N",
       subtitle = paste0("Intersection of the top-N gene lists from all four ",
                         "pipelines; monotonically increasing in N. ",
                         "Dashed line marks the value used (", TOP_N, ")."),
       x = "TOP_N", y = "Genes in the four-way intersection") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()),
  "R_Sensitivitas_TOP_N.png", 9, 6)


# === BAGIAN 14: REPRODUKTIBILITAS & RINGKASAN ================================

writeLines(capture.output(print(sessionInfo())), "R_sessionInfo.txt")
paket <- c("DESeq2","edgeR","caret","randomForest","kernlab","e1071","pROC","ggplot2")
versi <- data.frame(Paket = paket,
  Versi = sapply(paket, function(p)
    tryCatch(as.character(packageVersion(p)), error = function(e) "-")),
  row.names = NULL)
write.csv(versi, "R_Versi_Paket.csv", row.names = FALSE)

cat("\n=== RINGKASAN ===\n")
cat("Sampel        :", N_SAMPEL, "|", length(pasien), "pasien\n")
cat("Gen lolos     : DESeq2", nrow(dds), "| EdgeR", nrow(dge), "\n")
cat("Fitur ML      : DESeq2", ncol(X_deseq2), "| EdgeR", ncol(X_edger), "\n")
cat("Terbaik (AUC) :", metrik_all$Pipeline[which.max(metrik_all$Macro_AUC)],
    "—", round(max(metrik_all$Macro_AUC), 4), "\n")
cat("Uji DeLong    :", sum(uji_del$p_adj < 0.05), "/", nrow(uji_del),
    "lolos p_adj < 0,05\n")
cat("Bootstrap     :", sum(!boot_sel$Memuat_nol), "/", nrow(boot_sel),
    "selisih dengan CI tidak memuat nol\n")
cat("Core biomarker:", length(core_biomarker),
    if (length(core_biomarker) > 0)
      paste0(" (", paste(core_biomarker, collapse = ", "), ")") else "", "\n")
cat("\n", R.version.string, "\n")
cat("Berkas P_*.csv, P_*.png, R_*.csv, R_*.png tersimpan.\n")

###############################################################################
#  Keterbatasan yang TIDAK ditangani skrip ini dan perlu dinyatakan di naskah:
#  seleksi fitur dan penalaan hyperparameter dilakukan di luar lipatan validasi
#  (lihat skrip nested CV terpisah); ukuran sampel dan ketidakseimbangan kelas;
#  jaringan FFPE dan protokol 3SEQ; ketiadaan data klinis, luaran pasien,
#  validasi eksternal, dan konfirmasi laboratorium; bulk RNA-seq tidak
#  memisahkan sel tumor dari stroma.
###############################################################################
