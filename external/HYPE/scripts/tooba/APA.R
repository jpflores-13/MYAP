library("hictoolsr") 
library("InteractionSet")
library("plotgardener")
library("purrr")
library("RColorBrewer")

## Read in data & format as GInteractions

loops <- read.table(file = "./data/20210830_YAPP-diffLoopCounts-LRT-p0.05-FC0-formatted.txt", header = TRUE)

  
loop_gi <- as_ginteractions(loops[,-4])

## separate loops into significant gained

gainedLoops <- loop_gi[loop_gi$padj <= 0.01 &
                      loop_gi$log2FoldChange > 0]

# calc apa
looplist <- list(gainLoops = gainedLoops)

res <-10e3
buffer <- 10

binned <- binBedpe(bedpe = gainedLoops, 
                   res = 10e3,
                   a1Pos = "center",
                   a2Pos = "center")

looplist2 <- list (binloops = binned)


filteredLoops <- lapply(X = looplist2,
                        FUN = filterBedpe,
                        res = res,
                        buffer = buffer) |>
  `names<-`(value = names(looplist))


lapply(filteredLoops, summary)

# actually calc apa

ctlHicPath <- "./data/YAPP_HEK_control_inter_30.hic"
sorHicPath <- "./data/YAPP_HEK_sorbitol_inter_30.hic"

ctl2HiCPath <- "./data/HYPE_T47D_None_inter_30.hic"
naclHiCPath <- "./data/HYPE_T47D_NaCl_inter_30.hic"

loopApaCtlHic <-
  lapply(X = filteredLoops,
         FUN = calcApa,
         hic = ctlHicPath,
         norm = "KR",
         buffer = buffer)


loopApaSorHic <-
  lapply(X = filteredLoops,
         FUN = calcApa,
         hic = sorHicPath,
         norm = "KR",
         buffer = buffer)

loopApaCtl2Hic <-
  lapply(X = filteredLoops,
         FUN = calcApa,
         hic = ctl2HiCPath,
         norm = "KR",
         buffer = buffer)


loopApaNaClHic <-
  lapply(X = filteredLoops,
         FUN = calcApa,
         hic = naclHiCPath,
         norm = "KR",
         buffer = buffer)



## normalize?
nLoops <- lapply(filteredLoops, length)

loopApaCtlHic <- Map("/", loopApaCtlHic, nLoops)
loopApaSorHic <- Map("/", loopApaSorHic, nLoops)
loopApaCtl2Hic <- Map("/", loopApaCtl2Hic, nLoops)
loopApaNaClHic <- Map("/", loopApaNaClHic, nLoops)



#plotting with plotgardner 
pageCreate(width = 4.25, height = 3, showGuides = FALSE)



p <- pgParams(x = 1,
              y = .5,
              width = 1,
              height = 1,
              space = 0.075,
              zrange = c(0, max(unlist(c(loopApaCtlHic, loopApaSorHic)))),
              palette = colorRampPalette(brewer.pal(9, 'YlGnBu')))

#zrange = c(0, quantile(unlist(c(loopApaCtlHic, loopApaSorHic)),.95))

p2 <- pgParams(x = 1,
              y = .5,
              width = 1,
              height = 1,
              space = 0.075,
              zrange = c(0, (max(unlist(c(loopApaCtl2Hic,loopApaNaClHic)))-1)),
              palette = colorRampPalette(brewer.pal(9, 'YlGnBu')))



xpos <- c(p$x, p$x + p$width + p$space)
ypos <- c(p$y, p$y + p$height + p$space)

wt_plots <- 
  pmap(list(loopApaCtlHic, xpos[1], ypos[1]), \(a, x, y) {
    plotApa(params = p, apa = a, x = x, y = y)
  })

sor_plots <- 
  pmap(list(loopApaSorHic, xpos[1], ypos[2]), \(a, x, y) {
    plotApa(params = p, apa = a, x = x, y = y)
  })

ctl2_plots <- 
  pmap(list(loopApaCtl2Hic, xpos[2], ypos[1]), \(a, x, y) {
    plotApa(params = p2, apa = a, x = x, y = y)
  })

nacl_plots <- 
  pmap(list(loopApaNaClHic, xpos[2], ypos[2]), \(a, x, y) {
    plotApa(params = p2, apa = a, x = x, y = y)
  })



plotText(label = c("HEK", "T47D"),
         x = xpos + p$width / 2,
         y = ypos[1] - p$space,
         just = c('center', 'bottom'))

plotText(label = c("Control", "Treatment"),
         x = xpos[1] - p$space,
         y = ypos[1:2] + p$height / 2,
         rot = 90,
         just = c('center', 'bottom'))

annoHeatmapLegend(plot = wt_plots[[1]],
                  orientation = "h",
                  x = 1,
                  y = 2.7,
                  height = .125,
                  width = 1,
                  fontcolor = 'black')

annoHeatmapLegend(plot = ctl2_plots[[1]],
                  orientation = "h",
                  x = 2.1,
                  y = 2.7,
                  height = .125,
                  width = 1,
                  fontcolor = 'black')

