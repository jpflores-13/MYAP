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
library(ggplot2)
library(plyranges)

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


# Figure out % of loops overlapping promoters -----------------------------
##  overlap promoter regions & gained loops
promoters_gained <- subsetByOverlaps(gainedLoops, promoters(genes))
promoters_ctcf <- subsetByOverlaps(ctcfLoops, promoters(genes))

## % of loops at promoters
(length(promoters_gained)/length(gainedLoops))
(length(promoters_ctcf)/length(ctcfLoops))

# Stacked Barplot Visualization -------------------------------------------
## get % of gained loop anchors @ promoters

# create a dataset
loop_type <- c("Gained", "Gained", "Lost/Existing", "Lost/Existing")
condition <- c("Gained", "notGained", "Lost/Existing", "notLE")
value <- c((length(promoters_gained)/length(gainedLoops)),
           (length(promoters_ctcf)/length(ctcfLoops)),
           (1 - (length(promoters_gained)/length(gainedLoops))),
           (1 - (length(promoters_ctcf)/length(ctcfLoops))))
data <- data.frame(loop_type, condition, value)

# Stacked + percent
ggplot(data, aes(fill=condition, y=value, x = loop_type)) + 
  geom_bar(position="fill", stat="identity", width = 0.3) +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 1, hjust = 1),
        legend.position = "none") +
  labs(y = "",
       x = "") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#1D91C0", "#33A02C", "lightgrey", "lightgrey"))


# Figure out % of anchors overlapping promoters ---------------------------

## subset all anchors at starts of gained loops
gainedAnchors_first <- GRanges(seqnames = Rle(seqnames1(gainedLoops)), 
                               ranges = IRanges(start = start(anchors(gainedLoops, type = "first")),
                                                end = end(anchors(gainedLoops, type = "first")))) 

## subset all anchors at ends of gained loops
gainedAnchors_second <- GRanges(seqnames = Rle(seqnames2(gainedLoops)), 
                               ranges = IRanges(start = start(anchors(gainedLoops, type = "second")),
                                                end = end(anchors(gainedLoops, type = "second"))))

## concatenate 
gainedAnchors <- bind_ranges(gainedAnchors_first, gainedAnchors_second) |> 
  plyranges::mutate(loop_number = rep(paste0("loop_", c(1:355)), 2)) 

## subset all anchors at starts of existing/lost loops
ctcfAnchors_first <- GRanges(seqnames = Rle(seqnames1(ctcfLoops)), 
                               ranges = IRanges(start = start(anchors(ctcfLoops, type = "first")),
                                                end = end(anchors(ctcfLoops, type = "first"))))

## subset all anchors at ends of existing/lost loops
ctcfAnchors_second <- GRanges(seqnames = Rle(seqnames2(ctcfLoops)), 
                                ranges = IRanges(start = start(anchors(ctcfLoops, type = "second")),
                                                 end = end(anchors(ctcfLoops, type = "second"))))

## concatenate
ctcfAnchors <- bind_ranges(ctcfAnchors_first, ctcfAnchors_second) |> 
  plyranges::mutate(loop_number = rep(paste0("loop_", c(1:34102)), 2)) 
  

##  overlap promoter regions & gained loop anchors
promoters_gained <- subsetByOverlaps(promoters(genes),gainedAnchors)
promoters_ctcf <- subsetByOverlaps(promoters(genes), ctcfAnchors)

## % of anchors at promoters
(length(promoters_gained)/length(gainedAnchors))
(length(promoters_ctcf)/length(ctcfAnchors))

# Stacked Barplot Visualization -------------------------------------------
## get % of gained loop anchors @ promoters
library(ggplot2)

# create a dataset
loop_type <- c("Gained", "Gained", "Lost/Existing", "Lost/Existing")
condition <- c("Gained", "notGained", "Lost/Existing", "notLE")
value <- c((length(promoters_gained)/length(gainedAnchors)),
           (length(promoters_ctcf)/length(ctcfAnchors)),
           (1 - (length(promoters_gained)/length(gainedAnchors))),
           (1 - (length(promoters_ctcf)/length(ctcfAnchors))))
data <- data.frame(loop_type, condition, value)

# Stacked + percent
ggplot(data, aes(fill=condition, y=value, x = loop_type)) + 
  geom_bar(position="fill", stat="identity", width = 0.3) +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 1, hjust = 1),
        legend.position = "none") +
  labs(y = "",
       x = "") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#1D91C0", "#33A02C", "lightgrey", "lightgrey"))



# % of loops connecting 0 promoters ---------------------------------------
gainedLoops |> 
  linkOverlaps(promoters(genes))

ctcfLoops |> 
  linkOverlaps(promoters(genes))

## % of loops connecting 0 promoters
(length(promoters_gained)/length(gainedLoops))
(length(promoters_ctcf)/length(ctcfLoops))

# Stacked Barplot Visualization -------------------------------------------
## % of loops connecting 0 promoters
library(ggplot2)

# create a dataset
loop_type <- c("Gained", "Gained", "Lost/Existing", "Lost/Existing")
condition <- c("Gained", "notGained", "Lost/Existing", "notLE")
value <- c((length(promoters_gained)/length(gainedAnchors)),
           (length(promoters_ctcf)/length(ctcfAnchors)),
           (1 - (length(promoters_gained)/length(gainedAnchors))),
           (1 - (length(promoters_ctcf)/length(ctcfAnchors))))
data <- data.frame(loop_type, condition, value)

# Stacked + percent
ggplot(data, aes(fill=condition, y=value, x = loop_type)) + 
  geom_bar(position="fill", stat="identity", width = 0.3) +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 1, hjust = 1),
        legend.position = "none") +
  labs(y = "",
       x = "") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#1D91C0", "#33A02C", "lightgrey", "lightgrey"))



# % of E-P pairs ----------------------------------------------------------
enh_gained <- gainedAnchors[!gainedAnchors %over% promoters(genes)]
enh_ctcf <- ctcfAnchors[!ctcfAnchors %over% promoters(genes)]

gainedLoops_ep <- gainedLoops |> 
  linkOverlaps(promoters(genes), enh_gained)

ctcfLoops_ep <- ctcfLoops |> 
  linkOverlaps(promoters(genes), enh_ctcf)

## % of loops that connect enhancers to promoters
(nrow(gainedLoops_ep)/length(gainedAnchors))
(nrow(ctcfLoops_ep)/length(ctcfAnchors))

# Stacked Barplot Visualization -------------------------------------------
## % of EP pairs
library(ggplot2)

# create a dataset
loop_type <- c("Gained", "Gained", "Lost/Existing", "Lost/Existing")
condition <- c("Gained", "notGained", "Lost/Existing", "notLE")
value <- c((nrow(gainedLoops_ep)/length(gainedAnchors)),
           (nrow(ctcfLoops_ep)/length(ctcfAnchors)),
           (1 - (nrow(gainedLoops_ep)/length(gainedAnchors))),
           (1 - (nrow(ctcfLoops_ep)/length(ctcfAnchors))))
data <- data.frame(loop_type, condition, value)

# Stacked + percent
ggplot(data, aes(fill=condition, y=value, x = loop_type)) + 
  geom_bar(position="fill", stat="identity", width = 0.3) +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 1, hjust = 1),
        legend.position = "none") +
  labs(y = "",
       x = "") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#1D91C0", "#33A02C", "lightgrey", "lightgrey"))


# % of PP pairs -----------------------------------------------------------
gainedLoops_pp <- gainedLoops |> 
  linkOverlaps(promoters(genes))

ctcfLoops_pp <- ctcfLoops |> 
  linkOverlaps(promoters(genes))

## % of loops that connect promoters to promoters
(nrow(gainedLoops_pp)/length(gainedLoops))
(nrow(ctcfLoops_pp)/length(ctcfLoops))

# Stacked Barplot Visualization -------------------------------------------
## % of PP pairs
library(ggplot2)

# create a dataset
loop_type <- c("Gained", "Gained", "Lost/Existing", "Lost/Existing")
condition <- c("Gained", "notGained", "Lost/Existing", "notLE")
value <- c((nrow(gainedLoops_pp)/length(gainedLoops)),
           ((nrow(ctcfLoops_pp)/length(ctcfLoops))),
           (1 - (nrow(gainedLoops_pp)/length(gainedLoops))),
           (1 - (nrow(ctcfLoops_pp)/length(ctcfLoops))))
data <- data.frame(loop_type, condition, value)

# Stacked + percent
ggplot(data, aes(fill=condition, y=value, x = loop_type)) + 
  geom_bar(position="fill", stat="identity", width = 0.3) +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 1, hjust = 1),
        legend.position = "none") +
  labs(y = "",
       x = "") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#1D91C0", "#33A02C", "lightgrey", "lightgrey"))
