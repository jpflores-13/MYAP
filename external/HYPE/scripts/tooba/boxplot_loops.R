library(InteractionSet)

## Overlap blah blah

YAPsorbDEGHg19 <- readRDS("data/YAPsorbDEGHg19.rds")


## Reading in loops
## Read in data & format as GInteractions

loops <- read.table(file = "./data/20210830_YAPP-diffLoopCounts-LRT-p0.05-FC0-formatted.txt", header = TRUE)


loop_gi <- as_ginteractions(loops[,-4])

## separate loops into gained, lost, static

gainedLoops <- loop_gi[loop_gi$padj <= 0.01 &
                         loop_gi$log2FoldChange > 0]

lostLoops <- loop_gi[loop_gi$padj <= 0.01 &
                       loop_gi$log2FoldChange < 0]

staticLoops <- loop_gi[loop_gi$padj > 0.01]


## Find overlap between unique loop anchors

gained_hits <- subsetByOverlaps(YAPsorbDEGHg19, gainedLoops)
lost_hits <- subsetByOverlaps(YAPsorbDEGHg19, lostLoops)
static_hits <- subsetByOverlaps(YAPsorbDEGHg19, staticLoops)


## create dataframes of loops
gained_df <- as.data.frame(gained_hits)
gained_df$type <- "gained" 


lost_df <- as.data.frame(lost_hits)
lost_df$type <- "lost" 

static_df <- as.data.frame(static_hits)
static_df$type <- "static" 

library(dplyr)

combined1 <- bind_rows(gained_df, lost_df)
combined <- bind_rows(combined1, static_df)


## Plot the RNA log2FoldChange values of "gained", "static", "lost" loops as boxplot
library(ggplot2)

ggplot(data = combined,
       aes(x = reorder(type ,desc(type)),
           y = log2FoldChange)) +
  geom_hline(yintercept= 0,color = "gray", linetype = "dashed")+
  geom_boxplot(outlier.color = NA,
               fill = "#2E86C1",
               alpha = .5) +
  coord_cartesian(ylim = c(-1,1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x= "Loop") 

## wilcox test for significance 
wilcox.test(gained_df$log2FoldChange)


#Create List of genes in gained loops and save as .txt
cat(gained_df$symbol, sep = ",\n")

capture.output(cat(static_df$symbol, sep = ",\n"), file = "gainedloopgenes.txt")
