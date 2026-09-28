# Genome-Wide DNA Methylation Analysis Pipeline

This repository contains the final analytical report and R workflow for Group 5, focusing on genome-wide CpG DNA methylation profiling using Illumina Infinium microarray data.

## Authors
Alessia Visani, Veronica Boccaletti, Francesca Crispino, Valentina Molinari, Samuele Hamama, Mahsa Samadi, Tommaso Nicolis

## Project Overview
The pipeline processes raw fluorescence intensity files (`.idat`) to evaluate sample quality, perform background normalization, conduct exploratory data analysis (PCA and hierarchical clustering), and perform differential methylation analysis between Control (`CTRL`) and Disease (`DIS`) cohorts.

---

## Pipeline Workflow

1. **Data Loading (`minfi`)**
   - Imports raw signal data using a sample metadata sheet (`targets`).
   - Generates the initial `RGChannelSet` object (`RGset`).

2. **Fluorescence Inspection**
   - Extracts Red and Green channel intensity matrices.
   - Identifies probe chemistry and layout parameters for specific target addresses (e.g., Type I and Type II probes).

3. **Methylation Signal Conversion**
   - Converts the raw dataset into a `MethylSet` object (`MSet.raw`) to map physical coordinates to biological CpG coordinates.
   - Extracts Methylated and Unmethylated signal matrices.

4. **Quality Control (QC)**
   - **QC Plots:** Computes median methylation/unmethylation signals to identify outlier samples.
   - **Negative Controls:** Verifies background intensity thresholds using control probes.
   - **Detection $p$-values:** Computes detection $p$-values per sample to flag low-quality probes ($p > 0.05$).

5. **Beta and M-Value Extraction**
   - Computes continuous **Beta values** (proportion of methylation, ranging from 0 to 1) and **M-values** ($\log_2$ ratio of methylated to unmethylated signals).
   - Compares density distributions grouped by biological condition (`CTRL` vs. `DIS`).

6. **Normalization (`preprocessNoob`)**
   - Applies Normal-exponential out-of-band (Noob) background correction and dye-bias normalization.
   - Evaluates changes across probe chemistries (Type I vs. Type II) via comparative density plots and boxplots.

7. **Exploratory Data Analysis (PCA)**
   - Performs Principal Component Analysis on normalized Beta values.
   - Evaluates sample clustering against biological groups (`CTRL` vs. `DIS`), sex, and technical batch effects (`Slide`/`Sentrix_ID`).

8. **Differential Methylation Analysis**
   - Utilizes a non-parametric Wilcoxon-Mann-Whitney test to identify differentially methylated probes between groups.
   - Applies multiple testing corrections: **Bonferroni** and **Benjamini-Hochberg (BH)** false discovery rate (FDR) adjustments.

9. **Genomic Visualizations**
   - **Volcano Plot:** Visualizes delta beta values against nominal $- \log_{10}(p\text{-values})$.
   - **Manhattan Plot:** Displays chromosomal distribution of raw $p$-values.
   - **Heatmaps:** Generates hierarchical clustering heatmaps using Average, Complete, and Single linkage methods across the top variable probes.

---

## Key Conclusions
- **Batch Effects:** Exploratory analyses (PCA and clustering) revealed that technical variation (specifically experimental `Slide` batch effects) serves as the primary driver of variance within the dataset, overshadowing biological signals.
- **Multiple Testing:** Following Bonferroni and Benjamini-Hochberg corrections, zero probes reached genome-wide significance ($p_{\text{adj}} \le 0.05$), indicating that uncorrected nominal hits were largely driven by technical noise or lack of correction for confounding artifacts.

---

## Requirements & Dependencies
The analysis requires an R environment with Bioconductor packages:
- `minfi`
- `knitr`
- `IlluminaHumanMethylation450kanno.ilmn12.hg19`
- `factoextra`
- `qqman`
- `gplots`
