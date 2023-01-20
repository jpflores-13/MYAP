## Intersect Hi-C, ChIP, & ATAC data in untreated and sorbitol-treated EGFP-YAP HEK293 cells


# Load packages / data ----------------------------------------------------

library(plotgardener)
library(InteractionSet)
library(tidyverse)
library(glue)
library(dbscan)
library(data.table)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(mariner)
library(plyranges)

## Load loops
diff_loopCounts <- readRDS("data/processed/hic/hg38/diffLoops/noDroso/diffLoops_noDroso_10kb.rds")
diff_loopCounts <- interactions(diff_loopCounts) |> 
  as.data.frame() |> 
  mariner::as_ginteractions()

## Load ChIP & ATAC data
yap_chip <- list.files("data/raw/chip/hg38/YAP",
                       pattern = ".bw",
                       full.names = T)

tead_chip <- list.files("data/raw/chip/hg38/TEAD1",
                        pattern = ".bw",
                        full.names = T)

atac <- list.files("data/raw/atac/output/mergeSignal",
                   full.names = T)

# Setting static and lost loops -------------------------------------------

## gained loops with a p-adj. value of < 0.1 and a (+) log2FC
gained_adj <- 
  diff_loopCounts |> 
  subset(padj < 0.05 & log2FoldChange > 1)

## filter for the 100 best gained loops

## based off lowest padj value
bestGained <- head(gained_adj[order(gained_adj$padj, decreasing = F)], 100)

# top 100 gained loops
## If the below function doesn't work, might need to use swapAchors()

loopRegions_gained <- 
  GRanges(seqnames = as.character(seqnames(anchors(x = bestGained, "first"))),
          ranges = IRanges(start = start(anchors(bestGained, "first")),
                           end = end(anchors(bestGained, "second"))))

## Expand regions by buffer
buffer <- 200e3
loopRegions_gained <- loopRegions_gained + buffer

# Create Survey Plots -----------------------------------------------------

##make pdf
pdf(file = "plots/EGFP-YAP_ChIP_surveyPlot_10kb.pdf",
    width = 5.75,
    height = 9)

## Loop through each region
for(i in seq_along(loopRegions_gained)){
  
  ## Define parameters
  p <- pgParams(assembly = "hg38",
                resolution = 10e3,
                chrom = as.character(seqnames(loopRegions_gained))[i],
                chromstart = start(loopRegions_gained)[i],
                chromend = end(loopRegions_gained)[i],
                zrange = c(0,100),
                norm = "SCALE",
                x = 0.25,
                width = 5,
                height = 2,
                length = 5,
                fill = "#37a7db",
                linecolor = "#37a7db")
  
  ## Read YAP ChIP data
  
  YAP <- 
    lapply(yap_chip, \(file) {
      readBigwig(file = file,
                 params = p)
    })
  
  
  ## Read TEAD1 ChIP data
  TEAD <- 
    lapply(tead_chip, \(file) {
      readBigwig(file = file,
                 params = p)
    })
  
  ## Read ATAC data
  ATAC <- 
    lapply(atac, \(file) {
      readBigwig(file = file,
                 params = p)
    })
  
  
  # Begin Visualization -----------------------------------------------------
  ## Make page
  pageCreate(width = 5.75, height = 9,
             xgrid = 0, ygrid = 0, showGuides = F)
  
  ## Plot middle Hi-C rectangle + SIP `-isDroso = TRUE` & `-isDroso = TRUE` calls 
  
  control <- plotHicRectangle(data = "data/raw/hic/hg38/220716_dietJuicerMerge_condition/cont/YAPP_HEK_control_inter_30.hic",
                              params = p,
                              y = 0.5)
  
  annoHeatmapLegend(control, orientation = "v",
                    fontsize = 8,
                    fontcolor = "black",
                    digits = 2,
                    x = 5.5,
                    y = 0.5,
                    width = 0.1,
                    height = 1.5,
                    just = c("left", "top"),
                    default.units = "inches")
  
  ## Plot bottom Hi-C rectangle + SIP `-isDroso = TRUE` & `-isDroso = TRUE` calls
  
  sorb <- 
    plotHicRectangle(data = "data/raw/hic/hg38/220716_dietJuicerMerge_condition/sorb/YAPP_HEK_sorbitol_inter_30.hic", 
                     params = p,
                     y = 2.6)
  
  annoHeatmapLegend(sorb, orientation = "v",
                    fontsize = 8,
                    fontcolor = "black",
                    digits = 2,
                    x = 5.5,
                    y = 2.6,
                    width = 0.1,
                    height = 1.5,
                    length = 5,
                    just = c("left", "top"),
                    default.units = "inches")
  
  ## Plot control ChIP YAP track
  plotSignal(param = p,
             data = "data/raw/chip/hg38/YAP/YAPP_HEK293_EGFP-YAP_0hr_YAP.bw",
             x = 0.25,
             y = 4.7,
             height = 0.5,
             range = c(0, max(unlist(lapply(YAP, `[[`, 'score')))))
  
  ## Plot sorbitol ChIP YAP track
  
  plotSignal(param = p,
             data = "data/raw/chip/hg38/YAP/YAPP_HEK293_EGFP-YAP_1hr_YAP.bw",
             x = 0.25,
             y = 5.3,
             height = 0.5,
             range = c(0, max(unlist(lapply(YAP, `[[`, 'score')))))
  
  ## Plot control ChIP TEAD1 track
  plotSignal(param = p,
             data = "data/raw/chip/hg38/TEAD1/YAPP_HEK293_EGFP-YAP_0hr_TEAD1.bw",
             x = 0.25,
             y = 5.9,
             linecolor = "#abcc8e",
             fill = "#abcc8e",
             height = 0.5,
             range = c(0, max(unlist(lapply(TEAD, `[[`, 'score')))))
  
  ## Plot sorbitol ChIP TEAD1 track
  plotSignal(param = p,
             data = "data/raw/chip/hg38/TEAD1/YAPP_HEK293_EGFP-YAP_1hr_TEAD1.bw",
             x = 0.25,
             y = 6.5,
             linecolor = "#abcc8e",
             fill = "#abcc8e",
             height = 0.5,
             range = c(0, max(unlist(lapply(TEAD, `[[`, 'score')))))
  
  ## Plot ATAC control
  plotSignal(param = p, 
             data = "data/raw/atac/output/mergeSignal/YAPP_HEK_cont_0h.bw",
             y = 7.1,
             height = 0.5,
             linecolor = "#FFA500",
             fill = "#FFA500",
             range = c(0, max(unlist(lapply(ATAC, `[[`, 'score')))))
  
  ## Plot ATAC sorbitol
  plotSignal(param = p,
             data = "data/raw/atac/output/mergeSignal/YAPP_HEK_sorb_1h.bw",
             y = 7.7,
             height = 0.5,
             linecolor = "#FFA500",
             fill = "#FFA500",
             range = c(0, max(unlist(lapply(ATAC, `[[`, 'score')))))
  
  ## Plot genes
  plotGenes(param = p,
            chrom = p$chrom,
            x = 0.25,
            y = 8.2,
            height = 0.5)
  
  ## Plot genome label
  plotGenomeLabel(params = p,
                  x = 0.25,
                  y = 8.8)
  
  
  # Annotate Hi-C rectangles by treatment ------------------------------------
  
  plotText(label = "untreated",
           x = 0.25,
           y = 0.5,
           just = c("top", "left"))
  
  plotText(label = "sorbitol",
           x = 0.25,
           y = 2.6,
           just = c("top", "left"))
  
  plotText(label = "EGFP-YAP - untreated",
           x = 0.25,
           y = 4.65,
           fontcolor = "#1d91c0",
           fontsize = 8,
           just = c("top", "left"))
  
  plotText(label = "EGFP-YAP - sorbitol",
           x = 0.25,
           y = 5.25,
           fontsize = 8,
           fontcolor = "#1d91c0",
           just = c("top", "left"))
  
  plotText(label = "TEAD1 - untreated",
           x = 0.25,
           y = 5.85,
           fontsize = 8,
           fontcolor = "#abcc8e",
           just = c("top", "left"))
  
  plotText(label = "TEAD1 - sorbitol",
           x = 0.25,
           y = 6.45,
           fontsize = 8,
           fontcolor = "#abcc8e",
           just = c("top", "left"))
  
  plotText(label = "ATAC - untreated",
           x = 0.25,
           y = 7.05,
           fontsize = 8,
           fontcolor = "#FFA500",
           just = c("top", "left"))
  
  plotText(label = "ATAC - sorbitol",
           x = 0.25,
           y = 7.65,
           fontsize = 8,
           fontcolor = "#FFA500",
           just = c("top", "left"))
}
dev.off()