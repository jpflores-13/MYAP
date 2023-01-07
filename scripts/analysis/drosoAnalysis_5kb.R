## Observe differences between Droso datasets at 5kb resolution

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
  as_ginteractions() |> 
  binPairs(binSize = 5e3,
           pos1 = "center", pos2 = "center")

# Filter for gained loops -------------------------------------------------

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
buffer <- 100e3
loopRegions_gained <- loopRegions_gained + buffer

# load loop lists --------------------------------------------------------
## set up glue
drosoType <- c("bothDroso", "isDroso", "noDroso")

## get all loop files
loopFiles <- list.files(glue("data/processed/hic/hg38/diffLoops/{drosoType}"),
                        full.names = T,
                        pattern = ".rds")

## read as GInteractions & bin to 5kb res
loops <- loopFiles |> 
  lapply(readRDS) |> 
  lapply(interactions) |> 
  lapply(as.data.frame) |> 
  lapply(as_ginteractions) |> 
  lapply(binPairs, binSize = 5e3,
                    pos1 = "center", pos2 = "center")

## set up naming loop files
prefix <- c(drosoType)

## assign loop lists to name to keep track
names(loops) <- loopFiles |> 
  as.data.frame() |> 
  mutate(loopFiles = str_remove(loopFiles, paste0("data/processed/hic/hg38/diffLoops/", prefix, "/"))) |> 
  pull()

# Create Survey Plots -----------------------------------------------------

##make pdf
pdf(file = "plots/drosoAnalysis_5kb.pdf",
    width = 11,
    height = 5)

## Loop through each region
for(i in seq_along(loopRegions_gained)){
  
  ## Define parameters
  p <- pgParams(assembly = "hg38",
                resolution = 5e3,
                chrom = as.character(seqnames(loopRegions_gained))[i],
                chromstart = start(loopRegions_gained)[i],
                chromend = end(loopRegions_gained)[i],
                zrange = c(0,50),
                norm = "SCALE",
                x = 0.25,
                width = 3,
                length = 3,
                height = 3,
                fill = "#37a7db",
                linecolor = "#37a7db")
  
  
  # Begin Visualization -----------------------------------------------------
  
  ## Make page
  pageCreate(width = 11, height = 5,
             xgrid = 0, ygrid = 0, showGuides = F)
  
  ## Plot top left Hi-C square + SIP `-isDroso = FALSE` calls 
  
  control_noDroso <- plotHicSquare(data = "data/raw/hic/hg38/220716_dietJuicerMerge_condition/cont/YAPP_HEK_control_inter_30.hic",
                                   params = p,
                                   y = 0.5,
                                   half = "top")
  
  ## filter out lost loops from noDroso loop calls
  lostLoops_noDroso <- loops$diffLoops_noDroso_10kb.rds |> 
    subset(padj < 0.1 & log2FoldChange < 0)
  
  ## filter out gained loops from noDroso loop calls
  gainedLoops_noDroso <- loops$diffLoops_noDroso_10kb.rds |> 
    subset(padj < 0.1 & log2FoldChange > 0)
  
  ## annotate lost loops in control .hic file
  annoPixels(control_noDroso,
             data = lostLoops_noDroso,
             shift = 0.5,
             type = "box",
             col = "#005AB5")
  
  ## annotate gained loops in control .hic file
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
  
  ## add heatmap legend to display zrange
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
  
  ## annotate lost loops in sorb .hic file
  annoPixels(sorb_noDroso,
             data = lostLoops_noDroso,
             shift = 0.5,
             type = "box",
             col = "#005AB5")
  
  ## annotate gained loops in sorb .hic file
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
  

# Annotate Hi-C triangles by treatment ------------------------------------
  
  ## label top triangle
  plotText(label = "untreated",
           x = 0.25,
           y = 0.5,
           just = c("top", "left"),
           fontsize = 10)
  
  ## label bottom triangle
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
