## Generate loop counts/ call loops for cont/deg HCT116 cells processed

library(data.table)
library(hictoolsr)
library(mariner)
library(glue)
library(dbscan)
library(GenomeInfoDb)
library(InteractionSet)

# Compile loop files of interest ------------------------------------------
## load loop list
loops <- readRDS("data/processed/hic/hg38/diffLoops/HCT/mergedLoopsParentalCTCF.rds")

# Provide .hic files to extract from --------------------------------------

hic_files <- list.files(glue("data/raw/hic/hg38/HCT/parentalCTCF/{cond}"),
                        full.names = T,
                        pattern = "DEGR")

# Merge bedpe files into one and extract counts  --------------------------

h5 <- "data/processed/hic/hg38/h5/loops/HCT/hct_loopCounts_10kb.h5"
rds <- gsub(".h5$", ".rds", h5)
loopCounts_10kb <- loops |>
  mergePairs(radius = 5e3) |>
  binPairs(binSize = 10e3) |>
  pullHicPixels(binSize = 10e3,
                files = hic_files,
                h5File = h5,
                norm = "NONE",
                matrix = "observed")
saveRDS(loopCounts_10kb, file = rds)
# Create loop sets for each treatment -------------------------------------
source("scripts/utils/saveLoopCounts.R")

aux_loops <- "data/raw/hic/hg38/HCT/parentalCTCF/aux/5kbLoops.txt"
aux_hic <- "data/raw/hic/hg38/HCT/parentalCTCF/aux/MOPS_HCT_CTCFparental_5PhIAA_3h_inter_30.hic"
aux_h5 <- "data/processed/hic/hg38/h5/loops/HCT/aux_loopCounts_10kb.h5"

aux_loopCounts_10kb <- saveLoopCounts(aux_loops, aux_hic, aux_h5)

cont_loops <- "data/raw/hic/hg38/HCT/parentalCTCF/cont/5kbLoops.txt"
cont_hic <- "data/raw/hic/hg38/HCT/parentalCTCF/cont/MOPS_HCT_CTCFparental_Control_0h_inter_30.hic"
cont_h5 <- "data/processed/hic/hg38/h5/loops/HCT/cont_loopCounts_10kb.h5"

cont_loopCounts_10kb <- saveLoopCounts(cont_loops, cont_hic, cont_h5)

sessionInfo()
