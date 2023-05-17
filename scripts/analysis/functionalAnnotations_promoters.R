## Motif/GO enrichment of promoters & genes at gained loop anchors & existing/lost loops

## get gene symbols
ah <- AnnotationHub()
orgs <- subset(ah, ah$rdataclass == "OrgDb")
orgdb <- query(orgs, "Homo sapiens")[[1]]
genes_gainedLoops <- select(orgdb, keys = promoters_gained$UCSC_gene_id, 
                            columns = c("SYMBOL", "GENENAME")) |> 
  as.data.frame()

genes_ctcfLoops <- select(orgdb, keys = promoters_ctcf$UCSC_gene_id, 
                          columns = c("SYMBOL", "GENENAME")) |> 
  as.data.frame()


nrow(genes_gainedLoops)

## get GO terms
go_gained <- select(orgdb, keys = genes_gainedLoops$SYMBOL, "GO", "SYMBOL")
go2_gained <- select(GO.db, genes_gainedLoops$GO, c("TERM","DEFINITION"), "GOID") 
go2_gained$TERM

go_ctcf <- select(orgdb, keys = genes_ctcfLoops$SYMBOL, "GO", "SYMBOL")
go2_ctcf <- select(GO.db, genes_ctcfLoops$GO, c("TERM","DEFINITION"), "GOID") 
go2_ctcf$TERM

write_csv(genes_ctcfLoops, "tables/genes_atGainedLoopsAnchors.csv")
write_csv(genes_gainedLoops, "tables/genes_atGainedLoopsAnchors.csv")

# ## do any of these genes have multiple gained loops at their promoters?
# loi <- promoters |> 
#   subset(UCSC_gene_id %in% c("7004", "83937", "8463", "163732", "79142"))
# 
# countOverlaps(loi, gainedLoops)
