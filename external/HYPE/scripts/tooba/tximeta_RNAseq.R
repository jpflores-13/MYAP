## RNA-seq Differential Gene Expression Analysis for hg19

## Load required libraries
library(data.table)
library(plyranges)
library(tximeta)
library(readr)
library(DESeq2)
library(InteractionSet)
library(dplyr)
library(liftOver)

## Read in sample sheet as coldata
coldata <-
  fread("./data/rnaseq/YAPP_HEK_wt_1_RNApipeSamplesheet.txt") %>%
  as.data.frame()

## Use Just six hour time point 
coldata <- coldata[c(1,2,7,8),]


## Or try combining timepoints as replicates
#coldata$Replicates <- c(1,2, seq(1,6))


## Edit quant paths 
coldata$files <- paste0("data/rnaseq/quant/", coldata$sn, "/quant.sf")
setnames(coldata, "Name", "names")

## Check that quant paths are correct
file.exists(coldata$files)
  
## Load Linked Txome
loadLinkedTxome("GENCODE.v19.salmon_1.4.0.LinkedTxome.json")

## Import data with tximeta & summarize to gene
se <- tximeta(coldata)
gse <- summarizeToGene(se)

## Convert to factors
colData(gse)[] <- lapply(colData(gse), as.factor)

## Build DESeq object
dds <- DESeqDataSet(gse, design = ~Replicates + Treatment)

## Filter out lowly expressed genes (at least 10 counts in at least 4 samples)
keep <- rowSums(counts(dds) >= 10) >= 4
dds <- dds[keep,]

## Fit model
dds <- DESeq(dds)

## Get DE gene results as GRanges
de_genes <- results(dds,
                    name = "Condition_Sorb_vs_Control",
                    format = "GRanges") %>%
  names_to_column("gene_id")

## Add gene symbols to de_genes and
## remove non-standard chromosomes
de_genes <- 
  inner_join(x = as.data.table(de_genes),
             y = as.data.table(rowData(gse)) |>
               dplyr::select(c("gene_id", "symbol", "tx_ids")),
             by = "gene_id") |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE) |>
  keepStandardChromosomes(pruning.mode = "coarse")

## Add seqinfo to object
txdb <- 
  TxDb.Hsapiens.UCSC.hg38.knownGene |>
  keepStandardChromosomes()

seqlevels(de_genes) <- seqlevels(txdb)
seqinfo(de_genes) <- seqinfo(txdb)

## Lift over to hg19 -----------------------------------------------------------

## Create new object and change seqlevelstyles
de_genes_hg19 <- de_genes

## Import chain
ch <-
  system.file("extdata", "hg38ToHg19.over.chain",
              package = "liftOver") |>
  import.chain()

## Lift over and set genome
de_genes_hg19 <- 
  liftOver(de_genes_hg19, ch) |>
  unlist()

## Add seqinfo to object
txdb <- 
  TxDb.Hsapiens.UCSC.hg19.knownGene |>
  keepStandardChromosomes()

seqlevels(de_genes_hg19) <- seqlevels(txdb)
seqinfo(de_genes_hg19) <- seqinfo(txdb)

## Save results to a file
saveRDS(de_genes_hg19, file = "data/YAPsorbDEGHg19.rds")


