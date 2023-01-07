saveLoopCounts <- function(bedpe, hic, h5){
  
  ## convert bedpe files to GInteractions
  loops <- lapply(bedpe, fread) |> 
    lapply(as_ginteractions) 
  
  h5File <- h5
  rds <- gsub(".h5$", ".rds", h5File)
  ## extract counts from .hic files 
  loopCounts_10kb <- loops |>
    mergePairs(radius = 5e3) |>
    binPairs(binSize = 10e3) |>
    pullHicPixels(binSize = 10e3,
                  files = hic,
                  h5File = h5File,
                  norm = "NONE",
                  matrix = "observed")
  
  saveRDS(loopCounts_10kb, file = rds)
}