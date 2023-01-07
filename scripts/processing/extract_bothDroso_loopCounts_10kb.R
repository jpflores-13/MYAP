## Generate loop counts/ call loops for cont/sorb HEK293 cells processed with both `-isDroso` parameters (T & F)

library(data.table)
library(hictoolsr)
library(mariner)
library(glue)
library(dbscan)
library(GenomeInfoDb)
library(InteractionSet)

# Compile loop files of interest ------------------------------------------

## set a vector to pull from files from all conditions
cond <- c("cont", "sorb", "omega")

# "isDroso" loops are loops called with SIP setting `-isDroso true`
isDroso_loops <- list.files(path = glue("data/raw/hic/hg38/sip-loops/isDroso/{cond}"),
                            full.names = T,
                            pattern = "5kbLoops")

## "noDroso" loops are loops called with SIP setting `-isDroso false`
noDroso_loops <- list.files(path = glue("data/raw/hic/hg38/sip-loops/noDroso/{cond}"),
                            full.names = T,
                            pattern = "5kbLoops")

## combining both "isDroso" and "noDroso" loops
bothDroso_loops <- c(isDroso_loops, noDroso_loops)

# Provide .hic files to extract from --------------------------------------

hic_files <- list.files(path = glue("data/raw/hic/hg38/220722_dietJuicerCore/{cond}"), full.names = T)

# Merge bedpe files into one and extract counts  --------------------------

bothDroso_loops <- lapply(bothDroso_loops, fread) |> 
  lapply(as_ginteractions) 

h5 <- "data/processed/hic/hg38/h5/loops/bothDroso/bothDroso_loopCounts_10kb.h5"
rds <- gsub(".h5$", ".rds", h5)
loopCounts_10kb <- bothDroso_loops |>
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

sorb_loops <- "data/raw/hic/hg38/sip-loops/noDroso/sorb/5kbLoops.txt"
sorb_hic <- "data/raw/hic/hg38/220716_dietJuicerMerge_condition/sorb/YAPP_HEK_sorbitol_inter_30.hic"
sorb_h5 <- "data/processed/hic/hg38/h5/loops/bothDroso/bothDroso_sorb_loops_10kb.h5"

sorb_loopCounts_10kb <- saveLoopCounts(sorb_loops, sorb_hic, sorb_h5)

cont_loops <- "data/raw/hic/hg38/sip-loops/noDroso/cont/5kbLoops.txt"
cont_hic <- "data/raw/hic/hg38/220716_dietJuicerMerge_condition/cont/YAPP_HEK_control_inter_30.hic"
cont_h5 <- "data/processed/hic/hg38/h5/loops/bothDroso/bothDroso_cont_loops_10kb.h5"

cont_loopCounts_10kb <- saveLoopCounts(cont_loops, cont_hic, cont_h5)

omega_loops <- "data/raw/hic/hg38/sip-loops/noDroso/omega/5kbLoops.txt"
omega_hic <- "data/raw/hic/hg38/220716_dietJuicerMerge_omega/YAPP_HEK_inter_30.hic"
omega_h5 <- "data/processed/hic/hg38/h5/loops/bothDroso/bothDroso_omega_loops_10kb.h5"

omega_loopCounts_10kb <- saveLoopCounts(omega_loops, omega_hic, omega_h5)

sessionInfo()