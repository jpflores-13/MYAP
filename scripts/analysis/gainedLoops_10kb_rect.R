## Gained differential loop visualizations (HiC rectangles)

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

diff_loopCounts <- readRDS("data/processed/hic/hg38/diffLoops/bothDroso/diffLoops_bothDroso_10kb.rds") |> 
  interactions() |> 
  as.data.frame() |> 
  as_ginteractions()

# Setting static and lost loops -------------------------------------------

## gained loops with a p-adj. value of < 0.1 and a (+) log2FC
gained_adj <- 
  diff_loopCounts |> 
  subset(padj < 0.1 & log2FoldChange > 0)

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

# load loop lists --------------------------------------------------------
loops <- list.files(glue("data/processed/hic/hg38/diffLoops/bothDroso/"),
                        full.names = T,
                        pattern = ".rds") |> 
  readRDS() |> 
  interactions() |> 
  as.data.frame() |> 
  as_ginteractions()

lostLoops <- loops |> 
  subset(padj < 0.1 & log2FoldChange < 0)

gainedLoops <- loops |> 
  subset(padj < 0.1 & log2FoldChange > 0)

# Create Survey Plots -----------------------------------------------------

##make pdf
pdf(file = "plots/gainedLoops_10kb_rect_HK2chip.pdf",
    width = 5.75,
    height = 8)

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
                length = 5,
                height = 2,
                fill = "#37a7db",
                linecolor = "#37a7db")
  
  
  # Begin Visualization -----------------------------------------------------
  ## Make page
  pageCreate(width = 5.75, height = 8,
             xgrid = 0, ygrid = 0, showGuides = F)
  
  ## Plot top Hi-C rectangle + lost loops
  
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
  
  annoPixels(control,
             data = lostLoops,
             shift = 0.5,
             type = "arrow",
             col = "#005AB5")
  
  annoPixels(control,
             data = gainedLoops,
             shift = 0.5,
             type = "arrow",
             col = "#DC3220")
  
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
                    just = c("left", "top"),
                    default.units = "inches")
  
  annoPixels(sorb,
             data = lostLoops,
             shift = 0.5,
             type = "arrow",
             col = "#005AB5")
  
  annoPixels(sorb,
             data = gainedLoops,
             shift = 0.5,
             type = "arrow",
             col = "#DC3220")

  
  ## Plot untreated HK2 YAP1 ChIP-seq data
  plotSignal(param = p,
             data = "data/raw/chip/hg38/Share_Cai_HK2_YAP.bw",
             x = 0.25,
             y = 4.7,
             linecolor = "#abcc8e",
             fill = "#abcc8e",
             height = 0.5)
  
  ## Plot untreated HK2 TEAD1 ChIP-seq data
  plotSignal(param = p,
             data = "data/raw/chip/hg38/Share_Cai_HK2_TEAD1.bw",
             x = 0.25,
             y = 5.3,
             linecolor = "#abcc8e",
             fill = "#abcc8e",
             height = 0.5)
  
  ## Plot control ATAC track
  plotSignal(param = p,
             data = "data/raw/atac/output/mergeSignal/YAPP_HEK_cont_0h.bw",
             x = 0.25,
             y = 5.9,
             height = 0.5)
  
  ## Plot sorbitol ATAC track
  plotSignal(param = p,
             data = "data/raw/atac/output/mergeSignal/YAPP_HEK_sorb_1h.bw",
             x = 0.25,
             y = 6.5,
             height = 0.5)

  ## Plot genes
  plotGenes(param = p,
            chrom = p$chrom,
            x = 0.25,
            y = 7.1,
            height = 0.5)
  
  
  ## Plot genome label
  plotGenomeLabel(params = p,
                  x = 0.25,
                  y = 7.7)
  
  # Annotate Hi-C rectangles by treatment ------------------------------------
  
  plotText(label = "untreated",
           x = 0.25,
           y = 0.5,
           just = c("top", "left"))
  
  plotText(label = "+ sorbitol",
           x = 0.25,
           y = 2.6,
           just = c("top", "left"))
  
  plotText(label = "untreated - HK2 YAP1",
           x = 0.25,
           y = 4.65,
           fontcolor = "#1d91c0",
           fontsize = 8,
           just = c("top", "left"))
  
  plotText(label = "untreated - HK2 TEAD1",
           x = 0.25,
           y = 5.25,
           fontsize = 8,
           fontcolor = "#1d91c0",
           just = c("top", "left"))
  
  plotText(label = "untreated - ATAC",
           x = 0.25,
           y = 6.45,
           fontsize = 8,
           fontcolor = "#abcc8e",
           just = c("top", "left"))
  
  plotText(label = "+ sorbitol - ATAC",
           x = 0.25,
           y = 5.85,
           fontsize = 8,
           fontcolor = "#abcc8e",
           just = c("top", "left"))
  

}
dev.off()