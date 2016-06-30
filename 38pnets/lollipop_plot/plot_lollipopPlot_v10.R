setwd("/home/agoepfert/documents/pnet/38PNETs/Heatmap/20160527_matrix_v3/createPlots_server_v2")

source("https://bioconductor.org/biocLite.R")
#biocLite("rtracklayer")
#biocLite("Gviz")
## IMPORTANT check if the right version of trackViewer is installed -> 1.8.3 is necessary
## If syou have problems with BiocInstaller version 3.2 is too old for R version 3.3
## TRy troubleshooting or download the newest package of bioconductor into your R folder and extract it there
## https://www.bioconductor.org/packages/release/bioc/html/BiocInstaller.html
library("Gviz")
library("rtracklayer")
library("biomaRt")
#https://www.bioconductor.org/packages/release/bioc/vignettes/trackViewer/inst/doc/trackViewer.html#lolliplot
library("trackViewer")

getGene <- function(geneName){
  ## GET GENE STRUCTURE
  ##https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf
  ##############################################################################
  #ensembl_hg38 = useMart("ensembl",dataset="hsapiens_gene_ensembl",host="grch37.ensembl.org",path="/biomart/martservice")
  #ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
  ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
  x <- listAttributes(ensembl)
  
  geneLocations <- getBM(attributes = c("ensembl_gene_id","external_gene_name","exon_chrom_start","exon_chrom_end",
                                        "chromosome_name","start_position","end_position","ensembl_exon_id",
                                        "transcript_start","transcript_end",	"5_utr_start","5_utr_end","3_utr_start","3_utr_end"),
                         # filters = "ensembl_gene_id", values = geneName, mart=ensembl)
                         filters = "external_gene_name", values = geneName, mart=ensembl)
  
  geneLocations <- geneLocations[which(nchar(geneLocations$chromosome_name)<=2),]
  
  geneLocation.Chr <- geneLocations$chromosome_name[1]
  geneLocation.Chr.START <- geneLocations$start_position[1]
  
  geneLocations.EXON <- geneLocations[which(is.na(geneLocations$`5_utr_start`) &  is.na(geneLocations$`3_utr_start`)),]
  geneLocations.EXON$exon_chrom_length <- geneLocations.EXON$exon_chrom_end - geneLocations.EXON$exon_chrom_start
  
  geneLocations.5UTR <- geneLocations[which(geneLocations$`5_utr_start`>0),]
  geneLocations.5UTR$"5_utr_length" <- geneLocations.5UTR$`5_utr_end` - geneLocations.5UTR$`5_utr_start`
  
  geneLocations.3UTR <- geneLocations[which(geneLocations$`3_utr_start`>0),]
  geneLocations.3UTR$"3_utr_length" <- geneLocations.3UTR$`3_utr_end` - geneLocations.3UTR$`3_utr_start`
  
  result <- list(geneLocations.EXON,geneLocations.5UTR,geneLocations.3UTR,geneLocation.Chr)
  return(result)
}

extendInputMatrix <- function(data_matrix,normalSamples){
  colnames(data_matrix) <- gsub("_VSC ","", colnames(data_matrix))
  
  for (i in normalSamples){
    #cat(i)
    if (is.na(x = match(i,colnames(data_matrix)))==TRUE){
      #  if (match(i,colnames(data_matrix))==NA){
      #cat(i,"\n") 
      tmp_matrix <- rep.int(0,nrow(data_matrix))
      data_matrix <- cbind(data_matrix, tmp_matrix)
      nn<-ncol(data_matrix)
      colnames(data_matrix)[ncol(data_matrix)] <- i
    }
  }
  
  variants <- rownames(data_matrix)
  genes <- c()
  mutationtype <- c()
  SNP.start <- c()
  
  position1stUDS <- gregexpr(pattern="_",variants) #UDS = underscore
  for (i in 1:length(position1stUDS)){
    positionUDS <- min(position1stUDS[[i]])
    genes <- c(genes, substr(variants[i],1,positionUDS-1))
    position4thUDS <- position1stUDS[[i]][4]
    position5thUDS <- position1stUDS[[i]][5]
    mutationtype <- c(mutationtype, substr(variants[i],position4thUDS+1,position5thUDS-1))
    SNP.start <- c(SNP.start, substr(variants[i],position1stUDS[[i]][1]+1,position1stUDS[[i]][2]-1))
    
    ##cat(mutationtype)
  }
  timesMutated<-apply(data_matrix,1,function(x) {sum(x!=0)})
  
  #add genenames and mutationtype, mutationBurdon
  data_matrix <- cbind(data_matrix, genes)
  data_matrix <- cbind(data_matrix, mutationtype)
  data_matrix <- cbind(data_matrix, SNP.start)
  data_matrix <- cbind(data_matrix, timesMutated)
  
  
  remove(genes,mutationtype,SNP.start,timesMutated)
  
  return(data_matrix)
}

create_Lollipop_Plot <- function(geneStructure, data,single.gene,fileName,normalSamples){
  geneLocations.EXON <- as.data.frame(geneStructure[1])
  geneLocations.5UTR <- as.data.frame(geneStructure[2])
  geneLocations.3UTR <- as.data.frame(geneStructure[3])
  if (geneStructure[4]=="X"){
    geneStructure[4]=23
  }
  if (geneStructure[4]=="Y"){
    geneStructure[4]=24
  }
  geneLocation.Chr <- as.integer(geneStructure[4])
  #cat("DONE\n")
  #print(geneLocations.EXON)
  
  #print(geneLocations.5UTR)
  
  #print(geneLocations.3UTR)
  
  SNP <- as.double(as.vector(data$SNP.start))
  SNP_noDoubleEntries <- as.double(as.vector(data$SNP.start))
  
  mutation.SNP <- GRanges(geneLocation.Chr, IRanges(SNP, width=1, names=paste0("snp", SNP)))
  #print(mutation.SNP)
  
  #mutation.SNP$label <- as.character(1:length(mutation.SNP))
  mutation.SNP$score <-  as.double(as.vector(data$timesMutated))
  mutation.SNP$label.parameter.rot <- 45
  #cat(as.double(as.vector(data$timesMutated)),"\n")
  #mutation.SNP$score <- sample.int(5, length(mutation.SNP), replace = TRUE)
  
  if(nrow(geneLocations.EXON)>0){
    features <- GRanges(geneLocation.Chr, IRanges(geneLocations.EXON$exon_chrom_start, 
                                                  width=geneLocations.EXON$exon_chrom_length,
                                                  names=replicate(length(geneLocations.EXON$exon_chrom_start),"exon")))
    features$fill <- c("orange")
    
  }else{features <- GRanges(geneLocation.Chr, IRanges(geneLocations.5UTR$start_position[1], 
                                                  width=(geneLocations.5UTR$end_position[1]-geneLocations.5UTR$start_position[1]),
                                                  names="no exon"))
        features$fill <- c("white")
  }

  if(nrow(geneLocations.5UTR)>0){
    features2 <- GRanges(geneLocation.Chr, IRanges(geneLocations.5UTR$X5_utr_start, 
                                                   width=geneLocations.5UTR$X5_utr_length,
                                                   names=replicate(length(geneLocations.5UTR$X5_utr_start),"5utr")))
    features2$fill <- c("blue")
    }else{features2 <- GRanges(geneLocation.Chr, IRanges(geneLocations.EXON$start_position[1], 
                                                      width=(geneLocations.EXON$end_position[1]-geneLocations.EXON$start_position[1]),
                                                      names="no 5UTR"))
          features2$fill <- c("white")
}
  
  if(nrow(geneLocations.5UTR)>0){
  features3 <- GRanges(geneLocation.Chr, IRanges(geneLocations.3UTR$X3_utr_start, 
                                                 width=geneLocations.3UTR$X3_utr_length,
                                                 names=replicate(length(geneLocations.3UTR$X3_utr_start),"3utr")))
  features3$fill <- c("red")
  }else{features3 <- GRanges(geneLocation.Chr, IRanges(geneLocations.EXON$start_position[1], 
                                                       width=(geneLocations.EXON$end_position[1]-geneLocations.EXON$start_position[1]),
                                                       names="no 3UTR"))
        features3$fill <- c("white")
  }
  
  ### PLOT 1
  ##Legend Coloring
  colLegend = c("nonsynonymousSNV" = "blue", 
                "stopgain" = "red", 
                "stoploss" = "purple", 
                "frameshiftsubstitution" = "#008000", 
                "nonframeshiftsubstitution" = "orange", 
                "unkown" = "yellow")
  data$mutationtype <- sapply(data$mutationtype,as.character) 
  if(is.na(match("NA",data$mutationtype))==FALSE | is.na(match(".",data$mutationtype))==FALSE){
    data$mutationtype[data$mutationtype=="NA"| data$mutationtype=="."] <- "unkown"
  }
  
  for (j in 1:length(data$mutationtype)){
    data$colorMutationType[j] <- colLegend[as.character(data$mutationtype[j])]
  }
  
  legend <- unique(data$colorMutationType) ## legend fill color
  names(legend) <- unique(data$mutationtype)
  mutation.SNP$color <- data$colorMutationType  #sample.int(6, length(SNP), replace=TRUE)
  mutation.SNP$value1 <- 0
  mutation.SNP$value2 <- 0
  
  ### PLOT2
  data$countNormals <- rowSums(data[normalSamples] != 0)
  data$countTumors <- as.integer(data$timesMutated) - as.integer(data$countNormals)
  
  normal_tumor.gr <- mutation.SNP
  normal_tumor.gr$color <- rep(list(c('green', 'purple')), length(SNP))
  normal_tumor.gr$value1 <- as.double(as.vector(data$countNormal))
  normal_tumor.gr$value2 <- as.double(as.vector(data$countTumor))
  
  normal_tumor.Legend <- list(labels=c("Normal", "Tumor"), fill=c('green', 'purple'))
  
  lolliplot(mutation.SNP, c(features,features2,features3), legend=legend, type="circle",ylab="Amount of SNP")

  grid.text(paste(single.gene," - Chr.",geneLocation.Chr,sep=""), x=.5, y=.01, just="bottom")
  grid.text(paste(fileName,"\n Mutation/ SNP overview"), x=.5, y=.98, just="top",
            gp=gpar(cex=1, fontface="bold"))
  
  lolliplot(normal_tumor.gr, c(features,features2,features3), type="pie",legend=normal_tumor.Legend,ylab="Amount of Normal/Tumor Samples")

  grid.text(paste(single.gene," - Chr.",geneLocation.Chr,sep=""), x=.5, y=.01, just="bottom")
  grid.text(paste(fileName,"\n Overview tumor & normal samples"), x=.5, y=.98, just="top",
            gp=gpar(cex=1, fontface="bold"))
  
  library(gridExtra)
  x = ncol(data)-7
  y = ncol(data)
  tableShort1 <- data[,which(!apply(data[,1:x],2,FUN = function(x){all(x == 0)}))] # removes all columns/ samples with no SNP
  rowNamesTableShort1 <- rownames(data)
  
  tableShort2 <- data[,(x+1):y]
  
  tableShort1 <- lapply(tableShort1, function(x) (as.numeric(as.character(x))) )
  tableShort1 <- as.data.frame(lapply(tableShort1,round,4))
  rownames(tableShort1) <- rowNamesTableShort1
  
  plot.new()
  
  pushViewport(viewport(layout = grid.layout(3, 2, heights = unit(c(1, 5, 5), "null"))))   
  
  tt <- ttheme_default(base_size = 4,colhead=list(fg_params = list(parse=FALSE)))
  
  grid.table(d=tableShort1, theme=tt,vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))
  grid.table(d=tableShort2, theme=tt,vp = viewport(layout.pos.row = 3, layout.pos.col = 1:2))
  grid.text(paste("Input data set used for plotting gene ",single.gene,"\n(",fileName,")",sep=""), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
}


##INPUT
######################
normalSamples <- c("ch12","r3n","r4n")

##FILES
files <- list.files(path=".",pattern = "*.Rda")
files <- files[ !grepl(".pdf",files) ]

#numberFiles <- length(files)

#4 Include normals TRUE/ FALSE
includeNormals = TRUE

for (y in files){
  y<-"Matrix_gene_matrix_31052016_pathogenic_noMUC_VAF0.05_TD100.Rda" #works
  load(y)
  print(y)
  cat("Input File geladen\n")
  data_matrix <- as.data.frame(extendInputMatrix(data_matrix,normalSamples))
  input_genes <- unique(data_matrix$genes)
  
#  i=0
  pdf(file=paste("lollipopPlot_",y,".pdf",sep=""))#,onefile = TRUE
  
  for (single.gene in input_genes){
    ## gene structure
    ########################
    cat("Getting gene structure information of:",single.gene,"\n")
    
    geneStructure <- getGene(single.gene)
    input <- data_matrix[which(data_matrix$genes == single.gene),]
    cat("Start drawing lollippop plot of gene:",single.gene,"\n")
    create_Lollipop_Plot(geneStructure,input,single.gene,y,normalSamples)
#    i=i+1
    break
#     if(i==28){
#       break
#     }
  }
  dev.off()
  break
}


