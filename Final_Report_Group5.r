# finalreport, Group 5, 2026-06-11

#This is the final report for group 5 which authors are: Alessia Visani, Veronica Boccaletti, Francesca Crispino, Valentina Molinari, Samuele Hamama, Mahsa Samadi, Tommaso Nicolis

knitr::opts_chunk$set(echo = TRUE)

# 1. Load raw data with minfi and create an object called 'RGset' storing the RGChannelSet object.

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("minfi")
library(minfi)
library(knitr)

rm(list = ls())
directory = "/Users/aleonmac/Desktop/UNI/DNA_RNA/final_report/Input_Data" # Insert here your directory
setwd(directory)

baseDir <- (directory)
targets <- read.metharray.sheet(baseDir)
# Verify and print the number of imported samples
paste("Number of found samples: ", nrow(targets))

RGset <- read.metharray.exp(targets = targets)
save(RGset, file = "RGset.RData")
# Investigate if the RGChannelSet object was created correctly:
RGset

# 2. Create the dataframes Red and Green to store the red and green fluorescences respectively.

Red <- data.frame(getRed(RGset))
dim(Red)
head(Red)

Green <- data.frame(getGreen(RGset))
dim(Green)
head(Green)

# 3. Identify the Red and Green fluorescences for the assigned address.

df_I  <- data.frame(getProbeInfo(RGset, type = "I"))
head(df_I)

df_II <- data.frame(getProbeInfo(RGset, type = "II"))
head(df_II)

address <- "47683488"
if (address %in% df_I$AddressA | address %in% df_I$AddressB) {
  matched_row <- df_I[df_I$AddressA == address | df_I$AddressB == address, ]
  p_type  <- "Type I"
  p_color <- matched_row$Color
} else if (address %in% df_II$AddressA) {
  p_type  <- "Type II"
  p_color <- "-"
}
address <- "47683488"
matched_row <- df_II[df_II$AddressA == address,]

report_table <- data.frame(
  Sample = colnames(Red),
  `Red fluor` = as.numeric(Red[address, ]),
  `Green fluor` = as.numeric(Green[address, ]),
  Type = p_type,
  Color = p_color,
  check.names = FALSE
)
report_table

# 3.1. Check in the manifest file if the address corresponds to a Type I or a Type II probe and, in case of Type I probe, report its color.

load(file.path(directory, "Illumina450Manifest_clean.RData"))

type_A <- Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID == "47683488", 'Infinium_Design_Type']
type_B <- Illumina450Manifest_clean[Illumina450Manifest_clean$AddressB_ID == "47683488", 'Infinium_Design_Type']
paste("Type from AddressA_ID:",type_A)

paste("Type from AddressB_ID:",type_B)

# 4. Create the object MSet.raw

MSet.raw <- preprocessRaw(RGset)
MSet.raw
save(MSet.raw,file="MSet_raw.RData")

Meth <- as.matrix(getMeth(MSet.raw))
head(Meth)

Unmeth <- as.matrix(getUnmeth(MSet.raw))
head(Unmeth)

cpg_name <- df_II$Name[df_II$AddressA == address]
cpg_name

meth_values   <- as.numeric(Meth[cpg_name, ])
unmeth_values <- as.numeric(Unmeth[cpg_name, ])

red_values    <- as.numeric(Red[address, ])
green_values  <- as.numeric(Green[address, ])

identical(meth_values, green_values)

identical(unmeth_values, red_values)

# 5. Perform different quality checks and provide a brief comment to each step: QCplot; check the intensity of negative controls using minfi; calculate detection p-values for each sample, how many probes have a detection p-value higher than the threshold of 0.05?**
  # QCplot

qc <- getQC(MSet.raw)
qc

plotQC(qc)

  # Check the intensity of negative controls using `minfi`

getProbeInfo(RGset, type = "Control")

df_TypeControl <- data.frame(getProbeInfo(RGset, type = "Control"))
table(df_TypeControl$Type)

controlStripPlot(RGset, controls="NEGATIVE")

neg_addresses <- df_TypeControl$Address[df_TypeControl$Type == "NEGATIVE"]

valid_neg_addresses <- intersect(neg_addresses, rownames(Green))
neg_green_raw       <- Green[valid_neg_addresses, ]
neg_red_raw         <- Red[valid_neg_addresses, ]

neg_green_log2 <- log2(neg_green_raw)
neg_red_log2   <- log2(neg_red_raw)

mean_neg_green <- colMeans(neg_green_log2, na.rm = TRUE)
mean_neg_red   <- colMeans(neg_red_log2, na.rm = TRUE)

neg_control_summary <- data.frame(
  Mean_Log2_Green = round(mean_neg_green, 3),
  Mean_Log2_Red   = round(mean_neg_red, 3))
neg_control_summary

  # Calculate detection p-values for each sample, how many probes have a detection p-value higher than the threshold of 0.05?

threshold <- 0.05
detP <- detectionP(RGset)
head(detP)

failed <- detP>threshold
table(failed)

true_counts <- colSums(failed, na.rm = TRUE)

failed_table <- data.frame(
  Sample = colnames(failed),
  `n Failed Positions` = true_counts,
  `% Failed Probes` = (means_of_columns <- colMeans(failed))*100,
  row.names = NULL,
  check.names = FALSE
)
failed_table

# 6. Calculate raw beta and M values and plot the densities of mean methylation values, dividing the samples in CTRL and DIS (suggestion: subset the beta and M values matrixes in order to retain CTRL or DIS subjects and apply the function mean to the 2 subsets). Do you see any difference between the two groups?**

beta <- getBeta(MSet.raw)
M <- getM(MSet.raw)

summary(beta)

summary(M)

ctrl_dis_df <- data.frame(
  C_D = RGset@colData@listData[["Group"]],
  Samples = RGset@colData@rownames)
ctrl_dis_df

dis_Beta <- beta[, ctrl_dis_df$C_D=="DIS"]
ctrl_Beta <- beta[, ctrl_dis_df$C_D=="CTRL"]

dis_M <- M[, ctrl_dis_df$C_D=="DIS"]
ctrl_M <- M[, ctrl_dis_df$C_D=="CTRL"]

mean_ctrlBeta <- apply(ctrl_Beta, 1, mean, na.rm = TRUE)
mean_disBeta  <- apply(dis_Beta, 1, mean, na.rm = TRUE)
mean_ctrlM    <- apply(ctrl_M, 1, mean, na.rm = TRUE)
mean_disM     <- apply(dis_M, 1, mean, na.rm = TRUE)

par(mfrow = c(1, 2))
# Density of mean Beta values
plot(density(mean_ctrlBeta, na.rm = TRUE),
     main = "Density of Beta Values",
     col = "darkgoldenrod1", lwd = 2.5,
     xlab = "Mean Beta Values", ylab = "Density",
     xlim = c(-0.1,1.2), ylim = c(0,3.5))

lines(density(mean_disBeta, na.rm = TRUE), col = "lightgreen", lwd = 2.5)
legend("topright", legend = c("CTRL", "DIS"), col = c("darkgoldenrod1", "lightgreen"), lwd = 2, cex = 0.8, seg.len = 1, bty = "n")

# Density of mean M values
plot(density(mean_ctrlM, na.rm = TRUE),
     main = "Density of M Values",
     col = "darkgoldenrod1", lwd = 2.5,
     xlab = "Mean M Values", ylab = "Density",
     ylim = c(0,0.30))
lines(density(mean_disM, na.rm = TRUE), col = "lightgreen", lwd = 2.5)
legend("topright", legend = c("CTRL", "DIS"), col = c("darkgoldenrod1", "lightgreen"), lwd = 2, cex = 0.8, seg.len = 1, bty = "n")

par(mfrow = c(1, 1))

# 7. Normalize the data using the function preprocessNoob() and compare raw data and normalized data. Produce a plot with 6 panels in which, for both raw and normalized data, you show the density plots of beta mean values according to the chemistry of the probes, the density plot of beta standard deviation values according to the chemistry of the probes and the boxplot of beta values. Provide a short comment about the changes you observe. Do you think that the normalization approach that you used is appropriate considering this specific dataset? Try to color the boxplots according to the group (CTRL and DIS) and check whether the distribution of methylation values is different between the two groups, before and after normalization.

dfI <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type == "I", ]
dfII <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type == "II", ]
dfI <- droplevels(dfI)
dfII <- droplevels(dfII)

beta_raw <- getBeta(MSet.raw)
beta_I_raw <- beta_raw[rownames(beta_raw) %in% dfI$IlmnID, ]
beta_II_raw <- beta_raw[rownames(beta_raw) %in% dfII$IlmnID, ]

mean_beta_I_raw <- apply(beta_I_raw, 1, mean, na.rm = TRUE)
mean_beta_II_raw <- apply(beta_II_raw, 1, mean, na.rm = TRUE)
sd_beta_I_raw <- apply(beta_I_raw, 1, sd, na.rm = TRUE)
sd_beta_II_raw <- apply(beta_II_raw, 1, sd, na.rm = TRUE)

BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
MSet.norm <- preprocessNoob(RGset)
beta_norm <- getBeta(MSet.norm)  # beta values of the normalized data

beta_I_norm <- beta_norm[rownames(beta_norm) %in% dfI$IlmnID, ]
beta_II_norm <- beta_norm[rownames(beta_norm) %in% dfII$IlmnID, ]

mean_beta_I_norm <- apply(beta_I_norm, 1, mean, na.rm = TRUE)
mean_beta_II_norm <- apply(beta_II_norm, 1, mean, na.rm = TRUE)
sd_beta_I_norm <- apply(beta_I_norm, 1, sd, na.rm = TRUE)
sd_beta_II_norm <- apply(beta_II_norm, 1, sd, na.rm = TRUE)

dis_Beta_norm <- beta_norm[, ctrl_dis_df$C_D=="DIS"]
ctrl_Beta_norm <- beta_norm[, ctrl_dis_df$C_D=="CTRL"]

par(mfrow = c(2, 3))
# RAW DATA (first row)
# 1. Raw mean Beta densities
plot(density(mean_beta_I_raw, na.rm = TRUE), main = "Raw Beta Mean", col = "lightblue", lwd = 2, xlab = "Mean Beta")
lines(density(mean_beta_II_raw, na.rm = TRUE), col = "salmon", lwd = 2)
legend("topright", legend = c("Type I", "Type II"), col = c("lightblue", "salmon"), lwd = 2,  cex = 0.9, seg.len = 1, bty = "n")

# 2. Raw sd Beta densities
plot(density(sd_beta_I_raw, na.rm = TRUE), main = "Raw Beta SD", col = "lightblue", lwd = 2, xlab = "Beta SD")
lines(density(sd_beta_II_raw, na.rm = TRUE), col = "salmon", lwd = 2)
legend("topright", legend = c("Type I", "Type II"), col = c("lightblue", "salmon"), lwd = 2, cex = 0.9, seg.len = 1, bty = "n")

# 3. Raw Beta boxplot (colored by group)
boxplot(beta_raw, outline = FALSE, subset=c(ctrl_Beta, dis_Beta), col=c("darkgoldenrod1", "darkolivegreen2"),
        main = "Raw Beta Boxplot", ylab = "Beta Values", names = rep("", ncol(beta_raw)), ylim = c(0, 1.3))
legend("topright", legend = c("CTRL", "DIS"), fill = c("darkgoldenrod1", "darkolivegreen2"), bty = "n", cex = 0.7)



# NORMALIZED DATA (second row)
# 4. Normalized mean Beta densities
plot(density(mean_beta_I_norm, na.rm = TRUE), main = "Norm. Beta Mean", col = "lightblue", lwd = 2, xlab = "Mean Beta")
lines(density(mean_beta_II_norm, na.rm = TRUE), col = "salmon", lwd = 2)
legend("topright", legend = c("Type I", "Type II"), col = c("lightblue", "salmon"), lwd = 2, cex = 0.9, seg.len = 1, bty = "n")

# 5. Normalized sd Beta densities
plot(density(sd_beta_I_norm, na.rm = TRUE), main = "Norm. Beta SD", col = "lightblue", lwd = 2, xlab = "Beta SD")
lines(density(sd_beta_II_norm, na.rm = TRUE), col = "salmon", lwd = 2)
legend("topright", legend = c("Type I", "Type II"), col = c("lightblue", "salmon"), lwd = 2, cex = 0.9, seg.len = 1, bty = "n")

# 6. Normalized Beta boxplot (colored by group)
boxplot(beta_norm, outline = FALSE, subset=c(ctrl_Beta_norm, dis_Beta_norm), col=c("darkgoldenrod1", "darkolivegreen2"),
        main = "Norm. Beta Boxplot", ylab = "Beta Values", names = rep("", ncol(beta_norm)), ylim = c(0, 1.3)) 
legend("topright", legend = c("CTRL", "DIS"), fill = c("darkgoldenrod1", "darkolivegreen2"), bty = "n", cex = 0.7)
#Reset plotting layout
par(mfrow = c(1, 1))

# 8. Perform a PCA on the matrix of normalized beta values generated in step 7, after normalization. Comment the plot (do the samples divide according to the group? Do they divide according to the sex of the samples? Do they divide according to the batch, that is the column Sentrix_ID?).

pca_results <- prcomp(t(beta_norm))
options(repos = c(CRAN = "https://cloud.r-project.org"))
#install.packages(c("ggplot2", "ggpubr", "cluster", "factoMineR"))
#install.packages("factoextra")
library(factoextra)
fviz_eig(pca_results, addlabels = TRUE, xlab = 'PC number', ylab = '% of variance', barfill = "darkorange", barcolor = "darkgrey")

#plot(pca_results$x[,1], pca_results$x[,2],cex=1,pch=2)
targets$Group <- as.factor(targets$Group)
#levels(targets$Group)
group_colors <- c("darkgoldenrod1", "lightgreen")
plot(pca_results$x[,1], pca_results$x[,2],cex=4,pch=23, col="white",bg = group_colors[as.numeric(targets$Group)],xlab="PC1",ylab="PC2",xlim=c(-75,75),ylim=c(-75,75))
text(pca_results$x[,1], pca_results$x[,2],labels=rownames(pca_results$x),cex=0.5,pos=1)
legend("bottomright",legend=levels(targets$Group),col=c(1:nlevels(targets$Group)),pt.bg=group_colors,pch=23)

targets$Sex <- as.factor(targets$Sex)
#levels(targets$Sex)
sex_colors <- c("lightpink", "lightblue")
plot(pca_results$x[,1], pca_results$x[,2],cex=4,pch=23,col="white",bg = sex_colors[as.numeric(targets$Sex)],xlab="PC1",ylab="PC2",xlim=c(-75,75),ylim=c(-75,75))
text(pca_results$x[,1], pca_results$x[,2],labels=rownames(pca_results$x),cex=0.5,pos=1)
legend("bottomright",legend=levels(targets$Sex),pt.bg = sex_colors,pch=23)

targets$Slide <- as.factor(targets$Slide)
#in order to assign the colors i search for the unique slides
#unique(targets$Slide)
#table(targets$Slide)
#i have 5 unique slides
slide_colors <- c("coral", "brown1", "darkred", "lightblue", "lightgrey")
plot(pca_results$x[,1], pca_results$x[,2],cex=4,pch=23,col="white",bg = slide_colors[as.numeric(targets$Slide)],xlab = "PC1",ylab="PC2",xlim=c(-75,75),ylim=c(-75,75))
text(pca_results$x[,1], pca_results$x[,2],labels=rownames(pca_results$x),cex=0.5,pos=1)
legend("bottomright",legend=levels(targets$Slide),col="white",pt.bg=slide_colors,pch=23)

# 9. Using the matrix of normalized beta values generated in step 7, identify differentially methylated probes between group CTRL and group DIS using the function assigned to each group.

# This function takes a single row (probe) of beta values and compares them by group
  my_mannwhitney_test <- function(x) {
    # We use exact = FALSE because methylation data often contains ties
    wilcox <- wilcox.test(x ~ targets$Group, exact = FALSE)
    return(wilcox$p.value) 
  }
# Apply the test to every row (probe) in your normalized matrix
pValues_wilcox <- apply(beta_norm, 1, my_mannwhitney_test)
# Combine your beta values and the newly calculated p-values
test_df <- data.frame(beta_norm, pValues_wilcox = pValues_wilcox)
#keep only those with p-value <= 0.05
diff_meth_probes <- test_df[test_df$pValues_wilcox <= 0.05, ]
# Check how many you found
nrow(diff_meth_probes)

# 10. Apply multiple test correction and set a significant threshold of 0.05. How many probes do you identify as differentially methylated considering nominal pValues? How many after Bonferroni correction? How many after BH correction?**\

corrected_pValues_BH <- p.adjust(diff_meth_probes$pValues_wilcox,method="BH")
corrected_pValues_Bonf <-p.adjust(diff_meth_probes$pValues_wilcox,method="bonferroni")
sum(corrected_pValues_BH <= 0.05)
sum(corrected_pValues_Bonf <= 0.05)


# 11.Produce a volcano plot and a Manhattan plot of the results of differential methylation analysis

# Calculate means for each group

ctrl_indices <- which(targets$Group == "CTRL")
dis_indices <- which(targets$Group == "DIS")
mean_ctrl <- rowMeans(beta_norm[, ctrl_indices], na.rm = TRUE)
mean_dis <- rowMeans(beta_norm[, dis_indices], na.rm = TRUE)
delta_beta <- mean_dis - mean_ctrl


toVolcPlot <- data.frame(delta = delta_beta, logP = -log10(test_df$pValues_wilcox))
# 3. Plotting

plot(toVolcPlot$delta, toVolcPlot$logP, 
     pch = 16, cex = 0.5, 
     col = rgb(0.5, 0.5, 0.5, 0.3), # Grigio trasparente
     xlab = "Delta Beta (DIS - CTRL)", 
     ylab = "-log10(p-value)",
     main = "Volcano Plot: Differential Methylation")

abline(h = -log10(0.05), col = "salmon", lwd = 2, lty = 2) # Soglia p=0.05


sig_idx <- which(toVolcPlot$logP > -log10(0.05))
points(toVolcPlot$delta[sig_idx], toVolcPlot$logP[sig_idx], 
       pch = 16, cex = 0.7, col = "chocolate2")
# 6. Legenda
legend("topright", legend = c("NS", "Significant"), 
       col = c("gray", "chocolate2"), pch = 16, bty = "n")

# Manhattan plot
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# We need 'chr' and 'pos' (base position)
anno_df <- as.data.frame(anno[, c("chr", "pos")])
cols <- brewer.pal(8, "Set2")
manhattan_data <- data.frame(
  SNP = rownames(test_df),
  P = test_df$pValues_wilcox,
  CHR = anno_df[rownames(test_df), "chr"],
  BP = anno_df[rownames(test_df), "pos"]
)
manhattan_data$CHR <- as.numeric(gsub("chr", "", manhattan_data$CHR))
manhattan_data <- manhattan_data[!is.na(manhattan_data$CHR), ]
manhattan_data <- manhattan_data[!is.na(manhattan_data$P), ] 

library(qqman)
threshold_fdr <- -log10(0.05) 
manhattan(manhattan_data, 
          main = "Manhattan Plot (FDR corrected)",
          col = c("blue4", "chocolate2"),
          suggestiveline = FALSE, 
          genomewideline = threshold_fdr, 
          ylim = c(0, max(-log10(manhattan_data$P), 2) + 0.5), 
          ylab = "-log10(p-adjusted)")

# 12 Heatmap

#install.packages("gplots")
library(gplots)

vars <- apply(beta_norm, 1, var)
top_probes <- order(vars, decreasing = TRUE)[1:100]
input_heatmap <- as.matrix(beta_norm[top_probes, ])
group_colors_map <- c("CTRL" = "lightgrey", "DIS" = "rosybrown")

library("RColorBrewer")
group_color <- as.character(group_colors_map[targets$Group])
heatmap.2(input_heatmap,
          col = brewer.pal(11,"RdBu"),
          Rowv = TRUE,
          Colv = TRUE,
          hclustfun = function(x) hclust(x, method = 'average'),
          dendrogram = "both",
          key = TRUE,
          ColSideColors = group_color,
          density.info = "none",
          trace = "none",
          scale = "row",
          symm = FALSE,
          main = "Average linkage",
          key.xlab = 'beta-val',
          key.title = NA,
          keysize = 1,
          labRow = NA)
legend("topright",
       legend = levels(targets$Group),
       col = c('lightgrey', "rosybrown"),
       pch = 19,
       cex = 0.7)

heatmap.2(input_heatmap,
          col = brewer.pal(11,"RdBu"),
          Rowv = TRUE,
          Colv = TRUE,
          dendrogram = "both",
          key = TRUE,
          ColSideColors = group_color,
          density.info = "none",
          trace = "none",
          scale = "row",
          symm = FALSE,
          main = "Complete linkage",
          key.xlab = 'beta-val',
          key.title = NA,
          keysize = 1,
          labRow = NA)
legend("topright",
       legend = levels(targets$Group),
       col = c('lightgrey', 'rosybrown'),
       pch = 19,
       cex = 0.7)

heatmap.2(input_heatmap,
          col = brewer.pal(11,"RdBu"),
          Rowv = TRUE,
          Colv = TRUE,
          hclustfun = function(x) hclust(x, method = 'single'),
          dendrogram = "both",
          key = TRUE,
          ColSideColors = group_color,
          density.info = "none",
          trace = "none",
          scale = "row",
          symm = FALSE,
          main = "Single linkage",
          key.xlab = 'beta-val',
          key.title = NA,
          keysize = 1,
          labRow = NA)
legend("topright",
       legend = levels(targets$Group),
       col = c('lightgrey', 'rosybrown'),
       pch = 19,
       cex = 0.7)
