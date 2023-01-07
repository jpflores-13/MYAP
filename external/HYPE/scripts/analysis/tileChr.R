## Load libraries
library(plotgardener)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Set up tiles ------------------------------------------------------------

tile <- 
  tileGenome(seqlengths = seqlengths(TxDb.Hsapiens.UCSC.hg38.knownGene),
             tilewidth = 25e5,
             cut.last.tile.in.chrom = T) |> 
  keepStandardChromosomes(species = "Homo_sapiens", pruning.mode = "coarse")

tile <- tile[seqnames(tile) == "chr1"]

# Set-up params -----------------------------------------------------------
pdf(file = "external/HYPE/plots/HYPE_hic_SCALE_tile_chr1.pdf",
    width = 6,
    height = 7.5)

for (i in seq_along(tile)){
  
  p <- pgParams(assembly = "hg38",
                resolution = 10e3,
                chrom = as.character(seqnames(tile)[i]),
                chromstart = start(tile)[i],
                chromend = end(tile)[i],
                zrange = c(0,50),
                norm = "SCALE",
                x = 0.25,
                width = 5,
                length = 5,
                height = 2,
                fill = "#37a7db",
                linecolor = "#37a7db")
  
  
  # Begin visualization -----------------------------------------------------
  
  ## Make page
  pageCreate(width = 6, height = 7.5,
             xgrid = 0, ygrid = 0, showGuides = F)
  
  
  ## Plot Hi-C maps & legends
  omega <- plotHicRectangle(data = "external/HYPE/data/raw/hic/hg38/220718_dietJuicerMerge_omega/HYPE_T47D_inter_30.hic",
                            params = p,
                            y = 0.5) |> 
    annoHeatmapLegend(orientation = "v",
                      fontsize = 8,
                      fontcolor = "black",
                      digits = 2,
                      x = 5.5,
                      y = 0.5,
                      width = 0.1,
                      height = 1.5,
                      just = c("left", "top"),
                      default.units = "inches")
  
  cont <- plotHicRectangle(data = "external/HYPE/data/raw/hic/hg38/220718_dietJuicerMerge_treatment/cont/HYPE_T47D_None_inter_30.hic",
                           params = p,
                           y = 2.5) |> 
    annoHeatmapLegend(orientation = "v",
                      fontsize = 8,
                      fontcolor = "black",
                      digits = 2,
                      x = 5.5,
                      y = 2.5,
                      width = 0.1,
                      height = 1.5,
                      just = c("left", "top"),
                      default.units = "inches")
  
  nacl <- plotHicRectangle(data = "external/HYPE/data/raw/hic/hg38/220718_dietJuicerMerge_treatment/nacl/HYPE_T47D_NaCl_inter_30.hic",
                           params = p,
                           y = 4.5) |> 
    annoHeatmapLegend(orientation = "v",
                      fontsize = 8,
                      fontcolor = "black",
                      digits = 2,
                      x = 5.5,
                      y = 4.5,
                      width = 0.1,
                      height = 1.5,
                      just = c("left", "top"),
                      default.units = "inches")
  
  ## Plot Gene Track
  plotGenes(params = p,
            chrom = p$chrom,
            height = 0.75,
            y = 6.5)
  
  
  ## Plot Gene Region
  plotGenomeLabel(params = p,
                  y = 7.25)
  
  ## Plot Text
  plotText(label = "omega",
           x = 0.25,
           y = 0.5,
           just = c("top", "left"))
  
  plotText(label = "control",
           x = 0.25,
           y = 2.5,
           just = c("top", "left"))
  
  plotText(label = "NaCl",
           x = 0.25,
           y = 4.5,
           just = c("top", "left"))
  
}

dev.off()