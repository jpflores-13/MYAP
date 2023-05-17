## Designing FISH probes for HCT116 Lost loops cells

library(InteractionSet)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)

## Function for designing FISH probes
# designFISHprobes <- function(){
# }

lost_loopCands <- readRDS("data/processed/fish/candLoops/fish_LostLoopCandidates_hct.rds") |> 
  keepStandardChromosomes(pruning.mode = "coarse")

# Locate loop regions for lost loop candidates ----------------------------

## first anchor
loopCands_lost_first_chr12 <- lost_loopCands |> 
  anchors(type = "first") |> 
  subset(seqnames == "chr12") 

loopCands_lost_first_chr12 <- c(loopCands_lost_first_chr12[start(loopCands_lost_first_chr12) == 114660000])

## second anchor
loopCands_lost_second_chr12 <- lost_loopCands |> 
  anchors(type = "second") |> 
  subset(seqnames == "chr12") 

loopCands_lost_second_chr12 <- c(loopCands_lost_second_chr12[end(loopCands_lost_second_chr12) == 115760000])

loopCands_lost <- c(loopCands_lost_first_chr12, loopCands_lost_second_chr12)

## resize to 50kb anchors
loopCands_lost <- loopCands_lost + 20e3

# get sequences for loop candidates ---------------------------------------
## prep for bed file
loopCands_lostBed <- loopCands_lost |> 
  as.data.frame() |> 
  select(c(1:3))

loopCands_lostBed |> 
  write_tsv(file = "data/processed/fish/loopCands_lostBed_hct.txt",
              col_names = F)

## prep for fasta file
loopCands_lost_seq <- getSeq(Hsapiens, loopCands_lost)

loopCands_lost_seq |> 
  writeXStringSet(filepath = "data/processed/fish/loopCands_lost_hct.fasta")

