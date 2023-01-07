## check distributions of counts by normalization

library(data.table)
library(mariner)
library(InteractionSet)
library(strawr)
library(tidyverse)

## load in 
testLoops <- fread("data/processed/hic/dm6/loops/Pc_loopsKc_1kb.txt") |> 
  as_ginteractions()

hicNorms <- readHicNormTypes("data/raw/hic/dm/Kc_allcombined.hic")

extractedPixels <- list()
for(i in seq_along(hicNorms)){
  counts <- testLoops |> 
    pullHicPixels(binSize = 1e3,
                  files = "data/raw/hic/dm/Kc_allcombined.hic",
                  half = "both",
                  norm = hicNorms[i],
                  matrix = "observed") |> 
    assay() |> 
    as.matrix()
  extractedPixels[[i]] <- counts
  
}
extractedPixels <- do.call(cbind, extractedPixels) |> 
  as.matrix()
colnames(extractedPixels) <- hicNorms
extractedPixels_df <- extractedPixels |> 
  as.data.frame() 

## Visualize distribution of counts by norm 
preTransform <- extractedPixels_df |> 
  pivot_longer(cols = everything(),
               values_to = "value",
               names_to = "norm")

preTransform |> 
  ggplot(aes(y = value, color = norm)) +
  geom_density() +
  ylim(0,1000) +
  coord_flip() +
  theme_minimal() +
  theme(panel.grid.major = element_blank()) +
  ggtitle(label = "Extracted Counts by Normalization")


## log transform to make data normal
logTransform <- log(extractedPixels_df) |> 
  pivot_longer(cols = everything(),
               values_to = "value",
               names_to = "norm")

logTransform |> 
  ggplot(aes(y = value, color = norm)) +
  geom_density() +
  # ylim(0,1000) +
  coord_flip() +
  theme_minimal() +
  theme(panel.grid.major = element_blank()) +
  ggtitle(label = "Extracted Counts by Normalization",
          subtitle = "log transformed")

## square root transform to make data normal
sqrtTransform <- sqrt(extractedPixels_df) |> 
  pivot_longer(cols = everything(),
               values_to = "value",
               names_to = "norm")

sqrtTransform |> 
  ggplot(aes(y = value, color = norm)) +
  geom_density() +
  # ylim(0,1000) +
  coord_flip() +
  theme_minimal() +
  theme(panel.grid.major = element_blank()) +
  ggtitle(label = "Extracted Counts by Normalization",
          subtitle = "sqrt transformed")
  



