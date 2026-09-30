# Comparison of Breast Cancer Candidate Biomarker Pipelines (DEG and Machine Learning)

Multi-pipeline comparison of differential expression algorithms and machine
learning classifiers for candidate biomarker identification across breast
cancer progression stages.

**Dataset:** GSE47462 · **Design:** 2 DE methods × 2 ML algorithms = 4 complete
pipelines · **Validation:** grouped leave-one-patient-out cross-validation

---

## Background

Histopathology assesses a lesion as it is at the moment of biopsy, not its
potential to progress. Around 25% of ductal carcinoma in situ (DCIS) lesions
advance to invasive disease, and the ones that will cannot currently be
identified histologically. This drives both overtreatment of indolent lesions
and undertreatment of lesions that already harbour an invasive component.

Transcriptomic biomarkers could inform that decision. Most published biomarker
studies, however, stop at differential expression using a single pipeline, so
how far the resulting gene list depends on the analytical choices is rarely
tested. DESeq2 and edgeR model the same negative binomial distribution but
differ in normalisation and test statistic, and they return different gene
lists on identical input.

This repository runs four complete analytical pipelines side by side on one
dataset, with identical thresholds and one validation scheme, and reports where
they agree and where they do not.

## Study design

Two branches are kept independent from the raw count matrix onward. Nothing is
shared between them except the raw counts and the sample labels.

| Step | DESeq2 branch | EdgeR branch |
|---|---|---|
| Detection filter | count ≥ 10 in ≥ smallest group | `filterByExpr()` |
| Normalisation | median-of-ratios | TMM |
| Test statistic | Wald (negative binomial) | quasi-likelihood F-test |
| Transformation for ML | VST | logCPM |

Each branch feeds a Random Forest and an SVM with a radial basis kernel,
giving four pipelines: DESeq2+RF, DESeq2+SVM, EdgeR+RF, EdgeR+SVM. The
intersection of the two DE gene lists is reported as a descriptive statistic
only and is never used to select features.

## Dataset

[GSE47462](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47462) —
Brunner et al., *Genome Biology* 15:R71 (2014).

- 72 FFPE tissue samples from 25 patients
- 3SEQ (3′-end enriched RNA-Seq), Illumina Genome Analyzer IIx (GPL10999)
- Four progression stages: Normal (n=24), early neoplasia (n=25), DCIS (n=9),
  invasive ductal carcinoma (n=14)
- 22,774 genes in the raw count matrix

Only 3 of the 25 patients contribute all four stages. Because several samples
come from the same patient, cross-validation folds are grouped by patient:
every sample from a held-out patient is withheld together.

## Requirements

R ≥ 4.2 with Bioconductor.

```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "edgeR", "EnhancedVolcano"))
install.packages(c("caret", "e1071", "randomForest", "kernlab", "pROC",
                   "ggplot2", "pheatmap", "reshape2", "ggVennDiagram",
                   "patchwork", "dplyr", "tidyr", "scales"))
```

## Usage

```r
# 1. Set the path to the count matrix in BAGIAN 0
DATA_PATH <- "GSE47462 Raw Counts Refseq Genes.txt"

# 2. Run
source("multipipeline_deg_ml_final.R")
```

Runtime is roughly 15–25 minutes on an Apple M2 with 16 GB RAM, dominated by
the bootstrap resampling in section 10 and model training in section 8. Lower
`B_BOOT` from 1000 to 500 if that is too slow.

### Analysis parameters

| Parameter | Value | Applies to |
|---|---|---|
| `PADJ_CUTOFF` | 0.01 | padj (DESeq2) / FDR (edgeR) |
| `FC_CUTOFF` | 2.0 | \|log2 fold change\| |
| `TOP_N` | 20 | importance list length per pipeline |
| `N_TREE` | 500 | Random Forest |
| `TUNE_L` | 3 | caret `tuneLength` |
| `B_BOOT` | 1000 | patient-level bootstrap resamples |

Thresholds are identical across both branches by design, so that differences
in the results come from the method rather than from the cutoffs.

## Script structure

`multipipeline_deg_ml_final.R` runs end to end in 14 numbered sections:

| Section | Contents |
|---|---|
| 0 | Libraries and configuration |
| 1 | Data import and metadata |
| 2 | DESeq2 branch: filter, normalisation, DEG, VST |
| 3 | EdgeR branch: filter, TMM, DEG, logCPM |
| 4 | Descriptive comparison between branches (Jaccard index) |
| 5 | Flagging of artefact-prone genes (expression dominance) |
| 6 | DEG sensitivity: native effect-size test, paired design |
| 7 | Feature sets and ML matrices |
| 8 | Training of the four pipelines (grouped LOOCV) |
| 9 | Evaluation metrics and confusion matrices |
| 10 | Significance testing: DeLong test, patient-level bootstrap |
| 11 | Feature importance and TOP_N sensitivity |
| 12 | Core biomarkers |
| 13 | Figures |
| 14 | Reproducibility record and summary |

Sections 5, 6, 10 and the TOP_N sensitivity analysis in 11 are additive
robustness checks. They do not alter the primary analysis; they quantify how
much the primary results depend on specific methodological choices.

## Output

Files prefixed `P_` come from the primary analysis, `R_` from the robustness
checks.

**Tables** — differential expression results per contrast and branch, feature
composition, evaluation metrics, per-class AUC, confusion matrices, feature
importance, core biomarkers, DEG overlap with Jaccard indices, DeLong test
results, bootstrap confidence intervals, TOP_N sensitivity, artefact flags,
package versions.

**Figures** — PCA per branch, six volcano plots on a shared axis scale,
heatmaps of the top 50 DEG, Venn diagrams per contrast, ROC curves (macro-average
and per class), metric barplot and heatmap, bootstrap confidence intervals,
confusion matrices, per-class sensitivity, feature importance, four-way Venn,
TOP_N sensitivity curve.

`R_sessionInfo.txt` records the exact package versions of the run.

## Results

**Differential expression.** The number of genes passing threshold grows along
the progression axis in both branches: 12–17 genes for early neoplasia, 255–281
for DCIS, 409–416 for IDC. Jaccard indices between the two branches are 0.706
(EN), 0.691 (DCIS) and 0.797 (IDC). Disagreement concentrates on genes sitting
near the threshold, which is the expected consequence of applying two hard
cutoffs to methods that differ in normalisation and test statistic.

**Classification.** All four pipelines separate the four stages with a
macro-average AUC of 0.95–0.96 under grouped LOOCV. Discrimination is uneven
across stages: IDC is the easiest to recognise (AUC 0.991–0.993), reflecting
its larger and more widespread expression change, while early neoplasia
overlaps substantially with normal tissue.

**Candidate biomarkers.** Five genes appear in the top-20 importance list of
all four pipelines: **LAMC3, INHBA, TNNI2, SYT8, NDRG2**. INHBA is the only one
that rises with progression; the other four decline.

Numerical differences between pipelines are small. Section 10 tests them
formally with a paired DeLong test and a patient-level bootstrap rather than
comparing point estimates.

## Limitations

These constrain how far the results can be read, and are stated here rather
than left to the reader to infer.

**Feature selection sits outside the validation loop.** Differential
expression is computed once on all 72 samples and the resulting gene list is
used in every fold, so the held-out patient contributes to choosing the
features. Hyperparameters are likewise tuned on the same folds whose
predictions are reported. Both inflate the performance estimates
([Ambroise & McLachlan 2002](https://doi.org/10.1073/pnas.102102699);
[Varma & Simon 2006](https://doi.org/10.1186/1471-2105-7-91)). A nested
cross-validation variant is available for comparison.

**Statistical power.** 72 samples from 25 patients against >13,000 tested
genes. DCIS is represented by 9 samples. Confidence intervals on all metrics
are correspondingly wide.

**The `TOP_N = 20` cutoff determines the size of the core gene set.** The
intersection grows monotonically with N; five genes is a consequence of that
choice, not a quantity that emerges from the data. The sensitivity analysis in
section 11 reports the intersection across a range of N.

**SVM feature importance is not model-based.** caret has no model-specific
importance method for `svmRadial` and falls back to a univariate ROC filter, so
the gene ranking from the two SVM pipelines does not reflect the fitted kernel.
The output tables label this explicitly.

**Material and biology.** FFPE tissue and a 3′-end enriched protocol on a
first-generation sequencer. No clinical covariates, molecular subtype, or
patient outcome accompany the dataset, so no confounder can be adjusted for and
the prognostic value of the candidates cannot be assessed. Bulk RNA-seq does not
separate tumour epithelium from stroma. No external validation cohort and no
experimental confirmation (RT-qPCR, immunohistochemistry) were performed; the
direction of change is read from log2 fold change alone.

**Design.** The comparison is cross-sectional. The Normal → EN → DCIS → IDC
gradient compares different lesions at one point in time; it does not observe a
lesion changing.

## Citation

If this code is useful in your work, please cite the original dataset and this
repository.

```
Brunner AL, Li J, Guo X, et al. A shared transcriptional program in early
breast neoplasias despite genetic and clinical distinctions.
Genome Biology. 2014;15(5):R71. doi:10.1186/gb-2014-15-5-r71
```

## Author

Mohammad Fadhil Ihsan (10622010)
School of Life Sciences and Technology, Institut Teknologi Bandung

Supervisor: Popi Septiani, S.Si., M.Si., Ph.D.

Undergraduate thesis project, BI4092.
