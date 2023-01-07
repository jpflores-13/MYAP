#load libraries
  library("plotgardener")
  library("org.Hs.eg.db")
  library("TxDb.Hsapiens.UCSC.hg19.knownGene")
  library("tidyverse")

#get file names
  hek_ctl = "./data/YAPP_HEK_control_inter_30.hic"
  hek_sor = "./data/YAPP_HEK_sorbitol_inter_30.hic"
  t47d_ctl = "./data/HYPE_T47D_None_inter_30.hic"
  t47d_nacl = "./data/HYPE_T47D_NaCl_inter_30.hic"

#get loops
  loops = read.table(file = "./data/20210830_YAPP-diffLoopCounts-LRT-p0.05-FC0-formatted.txt", header = TRUE)
  
  gainloops = loops[loops$log2FoldChange > 0,]
  
  gainloops_nochr <-transform(gainloops,seqnames1=gsub("chr", "", seqnames1))
  gainloops_nochr <-transform(gainloops_nochr,seqnames2=gsub("chr", "", seqnames2))

#set variables  
  zquant = .85
  res = 5000
  filename = "surveyplot_res5000_Z85.pdf"
  

#for loop to print plots
  pdf(filename, width = 10.5, height = 6)
  
  for(i in 1:nrow(gainloops_nochr)){
    
    #set window 
    loopsize <- gainloops_nochr$end2[i]-gainloops_nochr$start1[i]
    windowStart <- gainloops_nochr$start1[i] - loopsize*0.3
    windowEnd <- gainloops_nochr$end2[i] + loopsize*0.3
    
    #assign chromosome
    nochr <- gainloops_nochr$seqnames1[i]
    chr <- gainloops$seqnames1[i]
    
    #set parameters 
    sharedParams_nochr<- pgParams(chrom = nochr,
                             assembly = "hg19",
                             norm = "KR",
                             resolution = res,
                             width = 4,
                             height = 1.5,
                             chromstart = windowStart,
                             chromend = windowEnd) 
    
    sharedParams_chr<- pgParams(chrom = chr,
                              assembly = "hg19",
                              norm = "KR",
                              resolution = res,
                              width = 4,
                              height = .5,
                              chromstart = windowStart,
                              chromend = windowEnd) 
    
    #set Z ranges
    hekctl = readHic(file = hek_ctl, 
                     params = sharedParams_nochr)
    Zhek <- quantile(hekctl$counts[hekctl$counts[i] > 0], zquant)
    
    t47dctl = readHic(file = t47d_ctl, 
                      params = sharedParams_nochr)
    Zt47d = quantile(t47dctl$counts[t47dctl$counts[i] > 0], zquant)
    
    #plots
    pageCreate(width = 10.5,
               height = 6,
               showGuides = FALSE)
    
    hekctlplot = plotHicRectangle(data = hek_ctl,
                                params = sharedParams_nochr,
                                x = 1, 
                                y = 1,
                                zrange = c(0,Zhek))
    heksorplot = plotHicRectangle(data = hek_sor,
                                params = sharedParams_nochr,
                                x = 1, 
                                y = 2.6,
                                zrange = c(0,Zhek))
    t47dctlplot = plotHicRectangle(data = t47d_ctl,
                                params = sharedParams_nochr,
                                x = 5.5, 
                                y = 1,
                                zrange = c(0,Zt47d))
    t47dnaclplot = plotHicRectangle(data = t47d_nacl,
                                params = sharedParams_nochr,
                                x = 5.5, 
                                y = 2.6,
                                zrange = c(0,Zt47d))
    #plot text
    plotText(label = "HEK Control",
             x = 1,
             y = 1.1,
             just = "left")
    plotText(label = "HEK Sorbitol",
             x = 1,
             y = 2.7,
             just = "left")
    plotText(label = "T47D Control",
             x = 5.5,
             y = 1.1,
             just = "left")
    plotText(label = "T47D NaCl",
             x = 5.5,
             y = 2.7,
             just = "left")
  
    #annotate loops
    annoPixels(plot = hekctlplot, 
               data = gainloops_nochr[i, -4])
  
    annoPixels(plot = heksorplot, 
               data = gainloops_nochr[i, -4])
    
    annoPixels(plot = t47dctlplot, 
               data = gainloops_nochr[i, -4])
    
    annoPixels(plot = t47dnaclplot, 
               data = gainloops_nochr[i, -4])
    
    #heatmaps
    annoHeatmapLegend(hekctlplot,
                      x= 5.125,
                      y= 1,
                      width = .125,
                      height = 1)
    annoHeatmapLegend(t47dctlplot,
                      x= 9.625,
                      y= 1,
                      width = .125,
                      height = 1)
    
    #genome and gene labels
    plotGenes(params = sharedParams_chr,
              x = 1,
              y = 4.15)
    
    annoGenomeLabel(hekctlplot,
                    x = 1,
                    y = 4.7,
                    width = 4,
                    scale = "Kb")
    plotGenes(params = sharedParams_chr,
              x = 5.5,
              y = 4.15)
    
    annoGenomeLabel(t47dctlplot,
                    x = 5.5,
                    y = 4.7,
                    width = 4,
                    scale = "Kb")
    }
   
  
  dev.off()


