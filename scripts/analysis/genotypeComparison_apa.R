## create APA plots comparing original and new Hi-C data

## load packages
library(mariner)
library(plotgardener)
library(glue)
library(hictoolsr)
library(RColorBrewer)

# load data ---------------------------------------------------------------
cond <- c("cont", "sorb")
og_hicFiles <- list.files(glue("data/raw/hic/hg38/220716_dietJuicerMerge_condition/{cond}"),
                          full.names = T)

new_hicFiles <- list.files(glue("data/raw/hic/hg38/230105_dietJuicerMerge_condition"),
                           full.names = T)

ogLoops <- readRDS("data/processed/hic/hg38/diffLoops/noDroso/diffLoops_noDroso_10kb.rds")
ogLoops <- ogLoops |> 
  interactions() |> 
  as.data.frame() |> 
  mariner::as_ginteractions() |> 
  subset(padj < 0.1 & log2FoldChange > 0) |> 
  filterBedpe(res = 10e3,
              buffer = 10)

# calculate APA matrices --------------------------------------------------

## OG 
og_contAPA <- ogLoops |> 
  pixelsToMatrices(buffer = 10) |> 
  pullHicMatrices(binSize = 10e3,
                  files = og_hicFiles[1], ##cont
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum) |> 
  as.matrix()

og_sorbAPA <- ogLoops |> 
  pixelsToMatrices(buffer = 10) |> 
  pullHicMatrices(binSize = 10e3,
                  files = og_hicFiles[2], ##cont
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum) |> 
  as.matrix()

## WT genotype
wt_contAPA <- ogLoops |> 
  pixelsToMatrices(buffer = 10) |> 
  pullHicMatrices(binSize = 10e3,
                  files = new_hicFiles[5], ##cont
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum) |> 
  as.matrix()

wt_sorbAPA <- ogLoops |> 
  pixelsToMatrices(buffer = 10) |> 
  pullHicMatrices(binSize = 10e3,
                  files = new_hicFiles[6], ##sorb
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum) |> 
  as.matrix()

## EGFP-YAP genotype
oe_contAPA <- ogLoops |> 
  pixelsToMatrices(buffer = 10) |> 
  pullHicMatrices(binSize = 10e3,
                  files = new_hicFiles[1], ##cont
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum) |> 
  as.matrix()

oe_sorbAPA <- ogLoops |> 
  pixelsToMatrices(buffer = 10) |> 
  pullHicMatrices(binSize = 10e3,
                  files = new_hicFiles[2], ##sorb
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum) |> 
  as.matrix()


## EGFP-YAPdTAD genotype
dTAD_contAPA <- ogLoops |> 
  pixelsToMatrices(buffer = 10) |> 
  pullHicMatrices(binSize = 10e3,
                  files = new_hicFiles[3], ##cont
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum) |> 
  as.matrix()

dTAD_sorbAPA <- ogLoops |> 
  pixelsToMatrices(buffer = 10) |> 
  pullHicMatrices(binSize = 10e3,
                  files = new_hicFiles[4], ##sorb
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum) |> 
  as.matrix()

## Divide each genotype/condition by nLoops
nLoops <- length(ogLoops)

og_contAPA <- (og_contAPA/nLoops)
og_sorbAPA <- (og_sorbAPA/nLoops)
oe_contAPA <- (oe_contAPA/nLoops)
oe_sorbAPA <- (oe_sorbAPA/nLoops)
wt_contAPA <- (wt_contAPA/nLoops)
wt_sorbAPA <- (wt_sorbAPA/nLoops)
dTAD_contAPA <- (dTAD_contAPA/nLoops)
dTAD_sorbAPA <- (dTAD_sorbAPA/nLoops)

## combine APA matrices to pull out the max value for zrange max
og_mats_combined <- c(og_contAPA,
                      og_sorbAPA)

oe_mats_combined <- c(oe_contAPA,
                      oe_sorbAPA)

wt_mats_combined <- c(wt_contAPA,
                      wt_sorbAPA)

dTAD_mats_combined <- c(dTAD_contAPA,
                        dTAD_sorbAPA)


# Visualize APAs ----------------------------------------------------------
pdf(file = "plots/genotypeComparison_apa.pdf",
    width = 10.5, height = 5.5)
pageCreate(width = 10.5, height = 5.5, showGuides = F)

p <- pgParams(assembly = "hg38",
              height = 2,
              width = 2,
              palette = colorRampPalette(brewer.pal(n = 9, "YlGnBu")))

plotApa(og_contAPA,
        x = 0.5,
        y = 0.5,
        zrange = c(0, max(og_mats_combined)),
        params = p) |> 
  annoHeatmapLegend(x = 2.6,
                    y = 1,
                    width = 0.1,
                    height = 1,
                    fontcolor = "black")

plotApa(og_sorbAPA,
        x = 0.5,
        y = 3,
        zrange = c(0, max(og_mats_combined)),
        params = p) |> 
  annoHeatmapLegend(x = 2.6,
                    y = 3.5,
                    width = 0.1,
                    height = 1,
                    fontcolor = "black")

plotApa(oe_contAPA,
        x = 3,
        y = 0.5,
        zrange = c(0, max(oe_mats_combined)),
        params = p)  |> 
  annoHeatmapLegend(x = 5.1,
                    y = 1,
                    width = 0.1,
                    height = 1,
                    fontcolor = "black")


plotApa(oe_sorbAPA,
        x = 3,
        y = 3,
        zrange = c(0, max(oe_mats_combined)),
        params = p) |> 
  annoHeatmapLegend(x = 5.1,
                    y = 3.5,
                    width = 0.1,
                    height = 1,
                    fontcolor = "black")

plotApa(wt_contAPA,
        x = 5.5,
        y = 0.5,
        zrange = c(0, max(wt_mats_combined)),
        params = p) |> 
  annoHeatmapLegend(x = 7.6,
                    y = 1,
                    width = 0.1,
                    height = 1,
                    fontcolor = "black")

plotApa(wt_sorbAPA,
        x = 5.5,
        y = 3,
        zrange = c(0, max(wt_mats_combined)),
        params = p) |> 
  annoHeatmapLegend(x = 7.6,
                    y = 3.5,
                    width = 0.1,
                    height = 1,
                    fontcolor = "black")

plotApa(dTAD_contAPA,
        x = 8,
        y = 0.5,
        zrange = c(0, max(dTAD_mats_combined)),
        params = p) |> 
  annoHeatmapLegend(x = 10.1,
                    y = 1,
                    width = 0.1,
                    height = 1,
                    fontcolor = "black")

plotApa(dTAD_sorbAPA,
        x = 8,
        y = 3,
        zrange = c(0, max(dTAD_mats_combined)),
        params = p) |> 
  annoHeatmapLegend(x = 10.1,
                    y = 3.5,
                    width = 0.1,
                    height = 1,
                    fontcolor = "black")

plotText(label = "Original", fontcolor = "black",
         x = 1.5,
         y = 0.35,
         just = c("center", "top"))

plotText(label = "OE", fontcolor = "black",
         x = 4,
         y = 0.35,
         just = c("center", "top"))

plotText(label = "WT", fontcolor = "black",
         x = 6.5,
         y = 0.35,
         just = c("center", "top"))

plotText(label = "dTAD", fontcolor = "black",
         x = 9,
         y = 0.35,
         just = c("center", "top"))

plotText(label = "untreated", fontcolor = "black",
         x = 0.35,
         y = 1.5,
         rot = 90)

plotText(label = "sorbitol", fontcolor = "black",
         x = 0.35,
         y = 4,
         rot = 90)
dev.off()