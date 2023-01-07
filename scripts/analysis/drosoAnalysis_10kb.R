## Observe differences between Droso datasets at 10kb resolution

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

diff_loopCounts <- readRDS("data/processed/hic/hg38/diffLoops/bothDroso/diffLoops_bothDroso_10kb.rds")
diff_loopCounts <- interactions(diff_loopCounts) |> 
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
drosoType <- c("bothDroso", "isDroso", "noDroso")
loopFiles <- list.files(glue("data/processed/hic/hg38/diffLoops/{drosoType}"),
                        full.names = T,
                        pattern = ".rds")

loops <- loopFiles |> 
  lapply(readRDS) |> 
  lapply(interactions) |> 
  lapply(as.data.frame) |> 
  lapply(as_ginteractions)

prefix <- c(drosoType)

names(loops) <- loopFiles |> 
  as.data.frame() |> 
  mutate(loopFiles = str_remove(loopFiles, paste0("data/processed/hic/hg38/diffLoops/", prefix, "/"))) |> 
  pull()

# Create Survey Plots -----------------------------------------------------

##make pdf
pdf(file = "plots/drosoAnalysis_10kb.pdf",
    width = 11,
    height = 4.25)

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
                width = 3,
                length = 3,
                height = 3,
                fill = "#37a7db",
                linecolor = "#37a7db")
  
  
  # Begin Visualization -----------------------------------------------------
  ## Make page
  pageCreate(width = 11, height = 4.25,
             xgrid = 0, ygrid = 0, showGuides = F)
  
  ## Plot top left Hi-C square + SIP `-isDroso = FALSE` calls 
  
  control_noDroso <- plotHicSquare(data = "data/raw/hic/hg38/220716_dietJuicerMerge_condition/cont/YAPP_HEK_control_inter_30.hic",
                           params = p,
                           y = 0.5,
                           half = "top")
  
  lostLoops_noDroso <- loops$diffLoops_noDroso_10kb.rds |> 
    subset(padj < 0.1 & log2FoldChange < 0)
  
  gainedLoops_noDroso <- loops$diffLoops_noDroso_10kb.rds |> 
    subset(padj < 0.1 & log2FoldChange > 0)
  
  annoPixels(control_noDroso,
             data = lostLoops_noDroso,
             shift = 0.5,
             type = "box",
             col = "#005AB5")
  
  annoPixels(control_noDroso,
             data = gainedLoops_noDroso,
             shift = 0.5,
             type = "box",
             col = "#DC3220")
  
  
  ## Plot bottom left Hi-C triangle + SIP `-isDroso = F` calls
  
  sorb_noDroso <- 
    plotHicSquare(data = "data/raw/hic/hg38/220716_dietJuicerMerge_condition/sorb/YAPP_HEK_sorbitol_inter_30.hic", 
                  params = p,
                  y = 0.5,
                  half = "bottom")
  
  annoHeatmapLegend(sorb_noDroso, orientation = "v",
                    fontsize = 8,
                    fontcolor = "black",
                    digits = 2,
                    x = 3.35,
                    y = 1,
                    width = 0.1,
                    height = 1.5,
                    just = c("left", "top"),
                    default.units = "inches")
  
  annoPixels(sorb_noDroso,
             data = lostLoops_noDroso,
             shift = 0.5,
             type = "box",
             col = "#005AB5")
  
  annoPixels(sorb_noDroso,
             data = gainedLoops_noDroso,
             shift = 0.5,
             type = "box",
             col = "#DC3220")
  
  ## Plot genes
  plotGenes(param = p,
            chrom = p$chrom,
            x = 0.25,
            y = 3.5,
            height = 0.5)
  
  
  ## Plot genome label
  plotGenomeLabel(params = p,
                  x = 0.25,
                  y = 4)
  
  # Annotate Hi-C rectangles by treatment ------------------------------------
  
  plotText(label = "untreated",
           x = 0.25,
           y = 0.5,
           just = c("top", "left"),
           fontsize = 10)
  
  plotText(label = "sorbitol-treated",
           x = 2.3,
           y = 3.35,
           just = c("top", "left"),
           fontsize = 10)
  
  
  ## Plot middle top Hi-C square + SIP `-isDroso true` calls 
  
  control_isDroso <- plotHicSquare(data = "data/raw/hic/hg38/220716_dietJuicerMerge_condition/cont/YAPP_HEK_control_inter_30.hic",
                                     params = p,
                                     x = 4,
                                     y = 0.5,
                                     half = "top")
  
  lostLoops_isDroso <- loops$diffLoops_isDroso_10kb.rds |> 
    subset(padj < 0.1 & log2FoldChange < 0)
  
  gainedLoops_isDroso <- loops$diffLoops_isDroso_10kb.rds |> 
    subset(padj < 0.1 & log2FoldChange > 0)
  
  annoPixels(control_isDroso,
             data = lostLoops_isDroso,
             shift = 0.5,
             type = "box",
             col = "#005AB5")
  
  annoPixels(control_isDroso,
             data = gainedLoops_isDroso,
             shift = 0.5,
             type = "box",
             col = "#DC3220")
  
  
  ## Plot middle bottom Hi-C square + SIP `-isDroso true` calls 
  
  sorb_isDroso <- 
    plotHicSquare(data = "data/raw/hic/hg38/220716_dietJuicerMerge_condition/sorb/YAPP_HEK_sorbitol_inter_30.hic", 
                  params = p,
                  x = 4,
                  y = 0.5,
                  half = "bottom")
  
  annoHeatmapLegend(sorb_isDroso, orientation = "v",
                    fontsize = 8,
                    fontcolor = "black",
                    digits = 2,
                    x = 7.1,
                    y = 1,
                    width = 0.1,
                    height = 1.5,
                    just = c("left", "top"),
                    default.units = "inches")
  
  annoPixels(sorb_isDroso,
             data = lostLoops_isDroso,
             shift = 0.5,
             type = "box",
             col = "#005AB5")
  
  annoPixels(sorb_isDroso,
             data = gainedLoops_isDroso,
             shift = 0.5,
             type = "box",
             col = "#DC3220")
  
  ## Plot genes
  plotGenes(param = p,
            chrom = p$chrom,
            x = 4,
            y = 3.5,
            height = 0.5)
  
  
  ## Plot genome label
  plotGenomeLabel(params = p,
                  x = 4,
                  y = 4)
  
  # Annotate Hi-C rectangles by treatment ------------------------------------
  
  plotText(label = "untreated",
           x = 4,
           y = 0.5,
           just = c("top", "left"),
           fontsize = 10)
  
  plotText(label = "sorbitol-treated",
           x = 6.1,
           y = 3.35,
           just = c("top", "left"),
           fontsize = 10)
  
  
  ## Plot right top Hi-C square + SIP `-bothDroso` calls 
  
  control_bothDroso <- plotHicSquare(data = "data/raw/hic/hg38/220716_dietJuicerMerge_condition/cont/YAPP_HEK_control_inter_30.hic",
                           params = p,
                           x = 7.75,
                           y = 0.5,
                           half = "top")
  
  lostLoops_bothDroso <- loops$diffLoops_bothDroso_10kb.rds |> 
    subset(padj < 0.1 & log2FoldChange < 0)
  
  gainedLoops_bothDroso <- loops$diffLoops_bothDroso_10kb.rds |> 
    subset(padj < 0.1 & log2FoldChange > 0)
  
  annoPixels(control_bothDroso,
             data = lostLoops_bothDroso,
             shift = 0.5,
             type = "box",
             col = "#005AB5")
  
  annoPixels(control_bothDroso,
             data = gainedLoops_bothDroso,
             shift = 0.5,
             type = "box",
             col = "#DC3220")
  
  
  ## Plot right bottom Hi-C square + SIP "bothDroso" calls
  
  sorb_bothDroso <- 
    plotHicSquare(data = "data/raw/hic/hg38/220716_dietJuicerMerge_condition/sorb/YAPP_HEK_sorbitol_inter_30.hic", 
                  params = p,
                  x = 7.75,
                  y = 0.5,
                  half = "bottom")
  
  annoHeatmapLegend(sorb_bothDroso, orientation = "v",
                    fontsize = 8,
                    fontcolor = "black",
                    digits = 2,
                    x = 10.85,
                    y = 1,
                    width = 0.1,
                    height = 1.5,
                    just = c("left", "top"),
                    default.units = "inches")
  
  annoPixels(sorb_bothDroso,
             data = lostLoops_bothDroso,
             shift = 0.5,
             type = "box",
             col = "#005AB5")
  
  annoPixels(sorb_bothDroso,
             data = gainedLoops_bothDroso,
             shift = 0.5,
             type = "box",
             col = "#DC3220")
  
  ## Plot genes
  plotGenes(param = p,
            chrom = p$chrom,
            x = 7.75,
            y = 3.5,
            height = 0.5)
  
  
  ## Plot genome label
  plotGenomeLabel(params = p,
                  x = 7.75,
                  y = 4)
  
  # Annotate Hi-C rectangles by treatment ------------------------------------
  
  plotText(label = "untreated",
           x = 7.75,
           y = 0.5,
           just = c("top", "left"),
           fontsize = 10)
  
  plotText(label = "sorbitol-treated",
           x = 9.85,
           y = 3.35,
           just = c("top", "left"),
           fontsize = 10)
  
  plotText(label = "-isDroso false",
           x = 1.75,
           y = 0.35)
  
  plotText(label = "-isDroso true",
           x = 5.5,
           y = 0.35)
  
  plotText(label = "bothDroso",
           x = 9.25,
           y = 0.35)
  
}
dev.off()