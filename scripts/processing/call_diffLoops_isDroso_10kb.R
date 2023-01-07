## Perform differential loops analysis with DESeq2

library(DESeq2)
library(InteractionSet)
library(glue)
library(stringr)
library(purrr)
library(plyranges)
library(mariner)
library(tidyverse)

# Load extracted loops counts ---------------------------------------------
loopCounts <- readRDS("data/processed/hic/hg38/h5/loops/isDroso/isDroso_loopCounts_10kb.rds")

# Extract count matrix for countData -------------------------------------------
m <- assay(loopCounts)

# Construct colData/metadata ----------------------------------------------

## String split the colnames 

colData <- as.data.frame(do.call(rbind, strsplit(colnames(m), "_")),
                         stringsAsFactors = T)
rownames(colData) <- colnames(m)
colnames(colData) <- c("Project", "Cell_Type", "Treatment", "Replicate")
colData <- colData[, c(1:4)]

## Correct replicate number
colData <- colData |> 
  mutate(Replicate = as.factor(rep(c(1:3), 2)))

## Make sure sample names match
all(colnames(m) == rownames(colData))

# Run DESeq2 --------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = m,
                              colData = colData,
                              design = ~ Replicate + Treatment)

## disable DESeq's default normalization 
sizeFactors(dds) <- rep(1, ncol(dds))

## Hypothesis testing with Wald with `betaPrior = F`
dds <- DESeq(dds)

## Plot dispersion estimates
plotDispEsts(dds)

# QC via data visualization before moving on ------------------------------

## plot PCA
plotPCA(vst(dds), intgroup = "Treatment") + ggplot2::theme(aspect.ratio = 1)
plotPCA(vst(dds), intgroup = "Treatment", returnData = TRUE)

## outputs
res <- results(dds)
resultsNames(dds)
summary(res)
res <- lfcShrink(dds, coef="Treatment_sorbitol_vs_control", type= "apeglm")

pdf(file = "plots/diffLoops_isDroso_10kb_MA.pdf")
plotMA(res, ylim=c(-4,4), main = "Differential Loop Analysis",
       ylab = "LFC",
       xlab = "mean of norm. counts")
dev.off()

# Concatenate loopCounts and DESeqResults ---------------------------------
mcols(interactions(loopCounts)) <- assay(loopCounts)
mcols(interactions(loopCounts)) <- cbind(mcols(interactions(loopCounts)), res)
diff_loopCounts <- loopCounts

## save as .rds
saveRDS(diff_loopCounts, file = "data/processed/hic/hg38/diffLoops/isDroso/diffLoops_isDroso_10kb.rds")

sessionInfo()