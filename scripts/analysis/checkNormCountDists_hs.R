## check distributions of counts by normalization

library(data.table)
library(mariner)
library(InteractionSet)
library(strawr)
library(tidyverse)
library(glue)

## load in loops to test
testLoops <- readRDS("data/processed/hic/hg38/diffLoops/bothDroso/diffLoops_bothDroso_10kb.rds") |> 
  interactions() |> 
  as.data.frame() |> 
  as_ginteractions() |> 
  binPairs(binSize = 1e3,
           pos1 = "center", pos2 = "center")

## set up different treatments for glue
cond <- c("sorb", "cont")
hicFiles <- list.files(glue("data/raw/hic/hg38/220716_dietJuicerMerge_condition/{cond}"),
                       full.names = T)

## get normalization methods for just one .hic file
hicNorms <- lapply(hicFiles[1], readHicNormTypes) |> 
  unlist()

## for loop to extract counts from each .hic file for each normalization method
extractedPixels <- list()
for(i in seq_along(hicFiles)){
  norms <- list()
  for(j in seq_along(hicNorms)){
    counts <- testLoops |> 
      pullHicPixels(binSize = 1e3,
                    files = hicFiles[i],
                    half = "both",
                    norm = hicNorms[j],
                    matrix = "observed") |> 
      assay() 
    
    print(counts)
    
    norms[[j]] <- counts
  }
  extractedPixels[[i]] <- norms
}
extractedPixels <- do.call(cbind, extractedPixels) |> 
  as.matrix()

# colnames(extractedPixels) <- 
colnames(extractedPixels)[1:5] <- paste0(hicNorms, "_cont")
colnames(extractedPixels)[6:10] <- paste0(hicNorms, "_sorb")

extractedPixels_df <- extractedPixels |> 
  as.data.frame() 


## Visualize distribution of counts by norm 
preTransform <- extractedPixels_df |> 
  rename(VC.SQRT_cont = VC_SQRT_cont,
         VC.SQRT_sorb = VC_SQRT_sorb) |> 
  pivot_longer(cols = everything(),
               values_to = "value",
               names_to = "norm") |> 
               separate(norm, into = c("norm", "file"), sep = '_')


preTransform |>  
  ggplot(aes(y = value, color = norm)) +
  geom_density() + 
  coord_flip()

## Visualize distribution of counts by norm after log transform
logTransform <- log(extractedPixels_df) |> 
  rename(VC.SQRT_cont = VC_SQRT_cont,
         VC.SQRT_sorb = VC_SQRT_sorb) |> 
  pivot_longer(cols = everything(),
               values_to = "value",
               names_to = "norm") |> 
  separate(norm, into = c("norm", "file"), sep = '_')

logTransform |> 
  ggplot(aes(y = value, color = norm)) +
  geom_density() +
  coord_flip() +
  theme_minimal() +
  theme(panel.grid.major = element_blank()) +
  ggtitle(label = "Extracted Counts by Normalization",
          subtitle = "log transformed")

## square root transform to make data normal
sqrtTransform <- sqrt(extractedPixels_df) |> 
  rename(VC.SQRT_cont = VC_SQRT_cont,
         VC.SQRT_sorb = VC_SQRT_sorb) |> 
  pivot_longer(cols = everything(),
               values_to = "value",
               names_to = "norm") |> 
  separate(norm, into = c("norm", "file"), sep = '_')

sqrtTransform |> 
  ggplot(aes(y = value, color = norm)) +
  geom_density() +
  # ylim(0,1000) +
  coord_flip() +
  theme_minimal() +
  theme(panel.grid.major = element_blank()) +
  ggtitle(label = "Extracted Counts by Normalization",
          subtitle = "sqrt transformed")



