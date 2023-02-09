## Obj: Look at which genes have promoters at the ends of gained loop anchors

library(InteractionSet)
library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(RMariaDB)
library(AnnotationHub)
library(GO.db)
library(scales)
library(forcats)
library(mariner)
library(memes)

## load in YAPP Hi-C differential loops
noDroso_loops <- readRDS("data/processed/hic/hg38/diffLoops/noDroso/diffLoops_noDroso_10kb.rds") |> 
  interactions() |> 
  as.data.frame() |>
  mariner::as_ginteractions()

## create new metadata columns for loop type & loop size
mcols(noDroso_loops)$loop_type <- case_when(
  mcols(noDroso_loops)$padj < 0.1 & mcols(noDroso_loops)$log2FoldChange > 0 ~ "gained",
  mcols(noDroso_loops)$padj < 0.1 & mcols(noDroso_loops)$log2FoldChange < 0 ~ "lost",
  mcols(noDroso_loops)$padj > 0.1 ~ "static",
  is.character("NA") ~ "other")

mcols(noDroso_loops)$loop_size <- pairdist(noDroso_loops)

## Add seqinfo to object
txdb <- 
  TxDb.Hsapiens.UCSC.hg38.knownGene |>
  keepStandardChromosomes()

seqlevels(noDroso_loops) <- seqlevels(txdb)
seqinfo(noDroso_loops) <- seqinfo(txdb)

## subset gained loops
gainedLoops <- noDroso_loops |> 
  subset(loop_type == "gained")

ctcfLoops <- noDroso_loops |> 
  subset(!loop_type == "gained")

## genes
txdb_UCSC <- makeTxDbFromUCSC(genome="hg38", tablename="knownGene",
                              transcript_ids=NULL,
                              circ_seqs=NULL,
                              url="http://genome.ucsc.edu/cgi-bin/",
                              goldenPath.url=getOption("UCSC.goldenPath.url"),
                              taxonomyId=NA,
                              miRBaseBuild=NA)
genes <- genes(txdb_UCSC)
genes <- data.frame(genes)
genes <- GRanges(seqnames = Rle(genes$seqnames), 
                 ranges = IRanges(start = genes$start, end = genes$end), 
                 strand = genes$strand,
                 UCSC_gene_id = genes$gene_id)

## subset all anchors at starts of gained loops
gainedAnchors_first <- GRanges(seqnames = Rle(seqnames1(gainedLoops)), 
                               ranges = IRanges(start = start(anchors(gainedLoops, type = "first")),
                                                end = end(anchors(gainedLoops, type = "first"))))

## subset all anchors at ends of gained loops
gainedAnchors_second <- GRanges(seqnames = Rle(seqnames2(gainedLoops)), 
                                ranges = IRanges(start = start(anchors(gainedLoops, type = "second")),
                                                 end = end(anchors(gainedLoops, type = "second"))))

## concatenate 
gainedAnchors <- c(gainedAnchors_first, gainedAnchors_second)

## subset all anchors at starts of existing/lost loops
ctcfAnchors_first <- GRanges(seqnames = Rle(seqnames1(ctcfLoops)), 
                             ranges = IRanges(start = start(anchors(ctcfLoops, type = "first")),
                                              end = end(anchors(ctcfLoops, type = "first"))))

## subset all anchors at ends of existing/lost loops
ctcfAnchors_second <- GRanges(seqnames = Rle(seqnames2(ctcfLoops)), 
                              ranges = IRanges(start = start(anchors(ctcfLoops, type = "second")),
                                               end = end(anchors(ctcfLoops, type = "second"))))

## concatenate
ctcfAnchors <- c(ctcfAnchors_first, ctcfAnchors_second)

## modify promoters
promoters <- promoters(genes) |> 
  keepStandardChromosomes(pruning.mode = c("coarse"))

##  overlap promoter regions & gained loop anchors
promoters_gained <- subsetByOverlaps(gainedAnchors, promoters)
promoters_ctcf <- subsetByOverlaps(ctcfAnchors, promoters)

# Prep for MEME ----------------------------------------------------------
human.genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

## get sequences of background peaks
sequence_promoters <- promoters |>
  get_sequence(human.genome)

memes::write_fasta(seq = sequence_promoters)

## get sequences of focal peaks
sequence_focal_gained <- promoters_gained |> 
  get_sequence(human.genome) 

memes::write_fasta(seq = sequence_focal_gained)

sequence_focal_ctcf <- promoters_ctcf |> 
  get_sequence(human.genome) 

memes::write_fasta(seq = sequence_focal_ctcf)

# input sequences on web browser MEME suite -------------------------------

ame_gained <- read.table("tables/loops/ame/gained/ame.tsv",
                         fill = T,
                         header = T)

ame_lost <- read.table("tables/loops/ame/lost/ame.tsv",
                       header = T)


# Visualization -----------------------------------------------------------
ame_gained_top <- ame_gained |> 
  mutate(log10pval = (-log10(p.value))) |>
  mutate(short_motif_id = str_remove(motif_ID, "_HUMAN.H11MO")) |> 
  arrange(rank) |> 
  slice_head(n = 50)

ame_lost_top <- ame_lost |> 
  mutate(log10pval = (-log10(p.value))) |> 
  mutate(short_motif_id = str_remove(motif_ID, "_HUMAN.H11MO")) |> 
  arrange(rank) |> 
  slice_head(n = 50)

top_ame_gained <- ame_gained_top |> 
  ggplot(aes(x = fct_reorder(short_motif_id, -log10pval), y = log10pval)) +
  geom_col(fill = "steelblue") +
  labs(title = "Gained Loop Anchors", x = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust=1))

top_ame_lost <- ame_lost_top |> 
  ggplot(aes(x = fct_reorder(short_motif_id, -log10pval), y = log10pval)) +
  geom_col(fill = "steelblue") +
  labs(title = "Lost Loop Anchors", x = "motif_id") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust=1))

library(plotgardener)

##make pdf
pdf(file = "plots/motifEnrichment_ame.pdf",
    width = 10,
    height = 7)

# Begin Visualization -----------------------------------------------------
## Make page
pageCreate(width = 10, height = 7,
           xgrid = 0, ygrid = 0, showGuides = F)

plotGG(top_ame_gained, 
       x = 0.5,
       y = 0.5,
       width = 9,
       height = 3)

plotGG(top_ame_lost,
       x = 0.5,
       y = 3.5,
       width = 9,
       height = 3)

dev.off()
