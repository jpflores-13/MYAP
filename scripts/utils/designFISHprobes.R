## Function for designing FISH probes

library(InteractionSet)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)

# designFISHprobes <- function(){
#   }

gained_loopCands <- readRDS("data/processed/hic/hg38/diffLoops/noDroso/fish_gainedLoopCandidates.rds") |> 
  keepStandardChromosomes(pruning.mode = "coarse")

lost_loopCands <- readRDS("data/processed/hic/hg38/diffLoops/noDroso/fish_lostLoopCandidates.rds") |> 
  keepStandardChromosomes(pruning.mode = "coarse")

# Locate loop regions of gained loop candidates ---------------------------
## first anchor
loopCands_gained_first_chr3 <- gained_loopCands |> 
  anchors(type = "first") |> 
  subset(seqnames == "chr3")

loopCands_gained_first_chr3 <- loopCands_gained_first_chr3[2]

loopCands_gained_first_chr9 <- gained_loopCands |> 
  anchors(type = "first") |> 
  subset(seqnames == "chr9")

loopCands_gained_first_chr9 <- loopCands_gained_first_chr9[3:5]

loopCands_gained_first_chr17 <- gained_loopCands |> 
  anchors(type = "first") |> 
  subset(seqnames == "chr17")

loopCands_gained_first_chr17 <- loopCands_gained_first_chr17[c(1,7)]

## second anchor
loopCands_gained_second_chr3 <- gained_loopCands |> 
  anchors(type = "second") |> 
  subset(seqnames == "chr3")

loopCands_gained_second_chr3 <- loopCands_gained_second_chr3[2]

loopCands_gained_second_chr9 <- gained_loopCands |> 
  anchors(type = "second") |> 
  subset(seqnames == "chr9")

loopCands_gained_second_chr9 <- loopCands_gained_second_chr9[3:5]

loopCands_gained_second_chr17 <- gained_loopCands |> 
  anchors(type = "second") |> 
  subset(seqnames == "chr17")

loopCands_gained_second_chr17 <- loopCands_gained_second_chr17[c(1,7)]

## concatentate GRanges
loopCands_gained <- c(loopCands_gained_first_chr3, loopCands_gained_first_chr9,
                      loopCands_gained_first_chr17, loopCands_gained_second_chr3,
                      loopCands_gained_second_chr9, loopCands_gained_second_chr17)

## resize on either end 
# loopCands_gained |> 
#   flank(50e3)

# get sequences for loop candidates ---------------------------------------
## prep for bed file
loopCands_gainedBed <- loopCands_gained |> 
  as.data.frame() |> 
  select(c(1:3))

loopCands_gainedBed |> 
  write.table(file = "data/processed/fish/loopCands_gainedBed.txt",
              col.names = F,
              row.names = F, 
              quote = F)

## prep for fasta file
loopCands_gained_seq <- getSeq(Hsapiens, loopCands_gained)

loopCands_gained_seq |> 
  writeXStringSet(filepath = "data/processed/fish/loopCands_gained.fasta")


# Locate loop regions for lost loop candidates ----------------------------

## first anchor
loopCands_lost_first_chr1 <- lost_loopCands |> 
  anchors(type = "first") |> 
  subset(seqnames == "chr1") 

loopCands_lost_first_chr1 <- c(loopCands_lost_first_chr1[start(loopCands_lost_first_chr1) == 183470000],
                               loopCands_lost_first_chr1[start(loopCands_lost_first_chr1) == 203270000],
                               loopCands_lost_first_chr1[start(loopCands_lost_first_chr1) == 47510000])

loopCands_lost_first_chr5 <- lost_loopCands |> 
  anchors(type = "first") |> 
  subset(seqnames == "chr5")

loopCands_lost_first_chr5 <- loopCands_lost_first_chr5[start(loopCands_lost_first_chr5) == 174410000] |> 
  unique()

loopCands_lost_first_chr13 <- lost_loopCands |> 
  anchors(type = "first") |> 
  subset(seqnames == "chr13")

loopCands_lost_first_chr13 <- loopCands_lost_first_chr13[start(loopCands_lost_first_chr13) == 84890000]

loopCands_lost_first_chr16 <- lost_loopCands |> 
  anchors(type = "first") |> 
  subset(seqnames == "chr16")

loopCands_lost_first_chr16 <- loopCands_lost_first_chr16[start(loopCands_lost_first_chr16) == 75210000]

## second anchor
loopCands_lost_second_chr1 <- lost_loopCands |> 
  anchors(type = "second") |> 
  subset(seqnames == "chr1") 

loopCands_lost_second_chr1 <- c(loopCands_lost_second_chr1[end(loopCands_lost_second_chr1) == 183640000],
                                loopCands_lost_second_chr1[end(loopCands_lost_second_chr1) == 203490000],
                                loopCands_lost_second_chr1[end(loopCands_lost_second_chr1) == 47710000]) |> 
  unique()

loopCands_lost_second_chr5 <- lost_loopCands |> 
  anchors(type = "second") |> 
  subset(seqnames == "chr5")

loopCands_lost_second_chr5 <- loopCands_lost_second_chr5[end(loopCands_lost_second_chr5) == 174800000]

loopCands_lost_second_chr13 <- lost_loopCands |> 
  anchors(type = "second") |> 
  subset(seqnames == "chr13")

loopCands_lost_second_chr13 <- loopCands_lost_second_chr13[end(loopCands_lost_second_chr13) == 85810000] |> 
  unique()

loopCands_lost_second_chr16 <- lost_loopCands |> 
  anchors(type = "second") |> 
  subset(seqnames == "chr16")

loopCands_lost_second_chr16 <- loopCands_lost_second_chr16[end(loopCands_lost_second_chr16) == 75390000] |> 
  unique()

loopCands_lost <- c(loopCands_lost_first_chr1,
                    loopCands_lost_first_chr5,
                    loopCands_lost_first_chr13,
                    loopCands_lost_first_chr16,
                    loopCands_lost_second_chr1,
                    loopCands_lost_second_chr5,
                    loopCands_lost_second_chr13,
                    loopCands_lost_second_chr16)

## resize on either end 
# loopCands_lost|> 
#   flank(50e3)

# get sequences for loop candidates ---------------------------------------
## prep for bed file
loopCands_lostBed <- loopCands_lost |> 
  as.data.frame() |> 
  select(c(1:3))

loopCands_lostBed |> 
  write_tsv(file = "data/processed/fish/loopCands_lostBed.txt",
            col_names = F)

## prep for fasta file
loopCands_lost_seq <- getSeq(Hsapiens, loopCands_lost)

loopCands_lost_seq |> 
  writeXStringSet(filepath = "data/processed/fish/loopCands_lost.fasta")



