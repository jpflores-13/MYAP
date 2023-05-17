library(InteractionSet)
library(plotgardener)
library(mariner)
library(tidyverse)

# Create survey plots to visualize FISH candidates ------------------------

##make pdf
pdf(file = "plots/rao_fishCandidate_chr4.pdf",
    width = 5.75,
    height = 5.6)

## Loop through each region
## Define parameters
p <- pgParams(assembly = "hg38",
              resolution = 10e3,
              chrom = "chr4",
              chromstart = 40800000,
              chromend = 42100000,
              zrange = c(0,100),
              norm = "SCALE",
              x = 0.25,
              width = 5,
              height = 2,
              length = 5,
              fill = "#37a7db",
              linecolor = "#37a7db")


# Begin Visualization -----------------------------------------------------
## Make page
pageCreate(width = 5.75, height = 5.6,
           xgrid = 0, ygrid = 0, showGuides = F)

## Plot middle Hi-C rectangle + SIP `-isDroso = TRUE` & `-isDroso = TRUE` calls 

control <- plotHicRectangle(data = "data/raw/hic/hg38/HCT/parentalCTCF/cont/MOPS_HCT_CTCFparental_Control_0h_inter_30.hic",
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

aux <- 
  plotHicRectangle(data = "data/raw/hic/hg38/HCT/parentalCTCF/aux/MOPS_HCT_CTCFparental_5PhIAA_3h_inter_30.hic", 
                   params = p,
                   y = 2.6) 

annoHeatmapLegend(aux, orientation = "v",
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

## Plot genes
plotGenes(param = p,
          chrom = p$chrom,
          x = 0.25,
          y = 4.65,
          height = 0.5)

## Plot genome label
plotGenomeLabel(params = p,
                x = 0.25,
                y = 5.25)

## Plot Title
plotText(label = "Fig.1, Rao et al 2017 HCT116 Example Region on chr4",
         x = 0.25,
         y = 0.35,
         just = c("left", "top"))

dev.off()
