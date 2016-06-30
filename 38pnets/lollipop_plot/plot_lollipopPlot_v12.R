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

##get the development versions of trackviewer
#library(BiocInstaller)
# useDevel()
# biocValid()              # checks for out of date packages
# biocLite("trackViewer")               # (optional) updates out of date packages
## get documentation
# library("utils")
# vignette("trackViewer")

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
                                        "transcript_start","transcript_end",	"5_utr_start","5_utr_end","3_utr_start","3_utr_end",
                                        "ensembl_transcript_id","transcript_length"
                                        ),
                         # filters = "ensembl_gene_id", values = geneName, mart=ensembl)
                         filters = "external_gene_name", values = geneName, mart=ensembl)
  
  geneLocations <- geneLocations[which(nchar(geneLocations$chromosome_name)<=2),]
  
  transcriptIDs <- unique(geneLocations$ensembl_transcript_id)
  
  #listCompleteTranscripts<-vector("list")
  listCompleteTranscripts <- data.frame()
  
  #i=1
  for (id in transcriptIDs){
    single.transcript <- geneLocations[which(geneLocations$ensembl_transcript_id==id),]
    cat("single.transcript$`5_utr_start`[1]: ",single.transcript$`5_utr_start`[1],"single.transcript$`3_utr_start`[nrow(single.transcript)]: ",single.transcript$`3_utr_start`[nrow(single.transcript)],"\n")
    if (!is.na(single.transcript$`5_utr_start`[1]) && !is.na(single.transcript$`3_utr_start`[nrow(single.transcript)])){
        single.transcript$"5utr_length" <- single.transcript$`5_utr_end` - single.transcript$`5_utr_start`
        single.transcript$"3utr_length" <- single.transcript$`3_utr_end` - single.transcript$`3_utr_start`
        single.transcript$"exon_chrom_length" <- single.transcript$exon_chrom_end - single.transcript$exon_chrom_start
        listCompleteTranscripts <- rbind(listCompleteTranscripts,single.transcript)
        #listCompleteTranscripts[[i]] <- single.transcript
        #i = i+1
    }
  }
  
  for (i in 1:nrow(listCompleteTranscripts)){
    #lolliplot features - necessary for the plot
    if(!is.na(listCompleteTranscripts$`5_utr_start`[i])){
      listCompleteTranscripts$lolliplot.names[i] <- "5utr"
      listCompleteTranscripts$lolliplot.start[i] <- listCompleteTranscripts$`5_utr_start`[i]
      listCompleteTranscripts$lolliplot.width[i] <- listCompleteTranscripts$`5utr_length`[i]
      listCompleteTranscripts$lolliplot.fill[i] <- "red"
      
    }else if(!is.na(listCompleteTranscripts$`3_utr_start`[i])){
      listCompleteTranscripts$lolliplot.names[i] <- "3utr"
      listCompleteTranscripts$lolliplot.start[i] <- listCompleteTranscripts$`3_utr_start`[i]
      listCompleteTranscripts$lolliplot.width[i] <- listCompleteTranscripts$`3utr_length`[i]
      listCompleteTranscripts$lolliplot.fill[i] <- "orange"
      
    }else{
      listCompleteTranscripts$lolliplot.names[i] <- "exon"
      listCompleteTranscripts$lolliplot.start[i] <- listCompleteTranscripts$exon_chrom_start[i]
      listCompleteTranscripts$lolliplot.width[i] <- listCompleteTranscripts$exon_chrom_length[i]
      listCompleteTranscripts$lolliplot.fill[i] <- "blue"
    }
  }

  geneLocation.Chr <- geneLocations$chromosome_name[1]
  geneLocation.Chr.START <- geneLocations$start_position[1]

  result <- list(listCompleteTranscripts,geneLocation.Chr)
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
  geneStructure.table <- as.data.frame(geneStructure[1])
  if (geneStructure[2]=="X"){
    geneStructure[2]=23
  }
  if (geneStructure[2]=="Y"){
    geneStructure[2]=24
  }
  geneLocation.Chr <- as.integer(geneStructure[2])
  #cat("DONE\n")

  SNP <- as.double(as.vector(data$SNP.start))
  SNP_noDoubleEntries <- as.double(as.vector(data$SNP.start))
  
  mutation.SNP <- GRanges(geneLocation.Chr, IRanges(SNP, width=1, names=paste0("snp", SNP)))
  #print(mutation.SNP)
  
  #mutation.SNP$label <- as.character(1:length(mutation.SNP))
  mutation.SNP$score <-  as.double(as.vector(data$timesMutated))
  mutation.SNP$label.parameter.rot <- 45
  #cat(as.double(as.vector(data$timesMutated)),"\n")
  #mutation.SNP$score <- sample.int(5, length(mutation.SNP), replace = TRUE)
  
  ##Features for multilayer
  features.mul <- GRanges(as.character(geneStructure[2]), IRanges(as.double(as.vector(geneStructure.table$lolliplot.start)), 
                                      width=as.double(as.vector(geneStructure.table$lolliplot.width)),
                                      names=as.vector(geneStructure.table$lolliplot.names)))
  features.mul$featureLayerID <- geneStructure.table$ensembl_transcript_id
  features.mul$fill <- geneStructure.table$lolliplot.fill
  
  

  
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
  
  lolliplot(mutation.SNP, features.mul, legend=legend, type="circle",ylab="Amount of SNP")

  grid.text(paste(single.gene," - Chr.",geneLocation.Chr,sep=""), x=.5, y=.01, just="bottom")
  grid.text(paste(fileName,"\n Mutation/ SNP overview"), x=.5, y=.98, just="top",
            gp=gpar(cex=1, fontface="bold"))
  
  lolliplot(normal_tumor.gr, features.mul, type="pie",legend=normal_tumor.Legend,ylab="Amount of Normal/Tumor Samples")

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
  
  i=0
  pdf(file=paste("lollipopPlot_",y,".pdf",sep=""))#,onefile = TRUE
  
  for (single.gene in input_genes){
    ## gene structure
    ########################
    cat("Getting gene structure information of:",single.gene,"\n")
    
    geneStructure <- getGene(single.gene)
    input <- data_matrix[which(data_matrix$genes == single.gene),]
    cat("Start drawing lollippop plot of gene:",single.gene,"\n")
    create_Lollipop_Plot(geneStructure,input,single.gene,y,normalSamples)
    i=i+1
#    break
     if(i==4){
       break
     }
  }
  dev.off()
  break
}


