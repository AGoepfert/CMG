setwd("/home/agoepfert/documents/pnet/38PNETs/Heatmap/20160527_matrix_v3/createPlots_server_v2")

source("https://bioconductor.org/biocLite.R")
#biocLite("rtracklayer")
#biocLite("Gviz")
## IMPORTANT check if the right version of trackViewer is installed -> 1.8.3 is necessary
## If syou have problems with BiocInstaller version 3.2 is too old for R version 3.3
## TRy troubleshooting or download the newest package of bioconductor into your R folder and extract it there
## https://www.bioconductor.org/packages/release/bioc/html/BiocInstaller.html
library(Gviz)
library(rtracklayer)
library(biomaRt)
# https://www.bioconductor.org/packages/release/bioc/vignettes/trackViewer/inst/doc/trackViewer.html#lolliplot
library(trackViewer)

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

extendInputMatrix <- function(data_matrix){
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
  
  remove(genes)
  return(data_matrix)
}






create_Lollipop_Plot <- function(geneStructure, data){
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
  cat("DONE\n")
  SNP <- as.double(as.vector(data$SNP.start))

  sample.gr <- GRanges(geneLocation.Chr, IRanges(SNP, width=1, names=paste0("snp", SNP)))
  #sample.gr$label <- as.character(1:length(sample.gr))
  sample.gr$score <-  sample.int(data$timesMutated, length(sample.gr), replace=TRUE)

  
  
  features <- GRanges(geneLocation.Chr, IRanges(geneLocations.EXON$exon_chrom_start, 
                                                width=geneLocations.EXON$exon_chrom_length,
                                                names=replicate(length(geneLocations.EXON$exon_chrom_start),"exon")))
  
  features2 <- GRanges(geneLocation.Chr, IRanges(geneLocations.5UTR$X5_utr_start, 
                                                 width=geneLocations.5UTR$X5_utr_length,
                                                 names=replicate(length(geneLocations.5UTR$X5_utr_start),"5utr")))
  
  features3 <- GRanges(geneLocation.Chr, IRanges(geneLocations.3UTR$X3_utr_start, 
                                                 width=geneLocations.3UTR$X3_utr_length,
                                                 names=replicate(length(geneLocations.3UTR$X3_utr_start),"3utr")))
  
  features$fill <- c("orange")
  features2$fill <- c("blue")
  features3$fill <- c("red")
  
  
  lolliplot(sample.gr, c(features,features2,features3))
}


















##INPUT
######################
files <- list.files(path=".",pattern = "*.Rda")
files <- files[ !grepl(".pdf",files) ]

#numberFiles <- length(files)

#4 Include normals TRUE/ FALSE
includeNormals = TRUE

for (y in files){
  y<-"Matrix_gene_matrix_31052016_pathogenic_noMUC_VAF0.05_TD100.Rda"
  load(y)
  print(y)
  cat("Input File geladen\n")
  data_matrix <- as.data.frame(extendInputMatrix(data_matrix))
  input_genes <- unique(data_matrix$genes)
  
  for (single.gene in input_genes){
    ## gene structure
    ########################
    cat("Getting gene structure information of:",single.gene,"\n")
    
    geneStructure <- getGene(single.gene)
    input <- data_matrix[which(data_matrix$genes == single.gene),]
    cat("Start drawing lollippop plot of gene:",single.gene,"\n")
    create_Lollipop_Plot(geneStructure,input)
    
    break
  }
  
  break
}


