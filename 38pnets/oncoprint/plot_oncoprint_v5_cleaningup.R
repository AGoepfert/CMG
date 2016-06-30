#setwd("/home/agoepfert/documents/pnet/38PNETs/Heatmap/20160527_matrix_v3/createPlots_server")

source("https://bioconductor.org/biocLite.R")
#library(circlize)
#biocLite("ComplexHeatmap") #install ComplexHeatmap
#biocLite()
library(gplots)
#complexHeatmap
library(ComplexHeatmap)
library(circlize)
require(lattice)


createOncoprint <- function(normalSamplesMatrix,tm_SamplesMatrix,data_matrix,plotTitle,includeNormals){
  ##3.adjust input data set
  #remove _VSC in colnames
  colnames(data_matrix) <- gsub("_VSC ","", colnames(data_matrix))
  #remove normal samples
  if (includeNormals==FALSE){
    data_matrix <- data_matrix[ , !colnames(data_matrix) %in% c(normalSamplesMatrix[,'normalsamples'])]
  }else{
    tm_SamplesMatrix <- rbind(tm_SamplesMatrix,normalSamplesMatrix)
  }
  
  for (i in tm_SamplesMatrix[,'tm_samples']){
    #cat(i)
    if (is.na(x = match(i,colnames(data_matrix)))==TRUE){
   #  if (match(i,colnames(data_matrix))==NA){
      #cat(i,"\n") 
      tmp_matrix <- rep.int(0,nrow(data_matrix))
      data_matrix <- cbind(data_matrix, tmp_matrix)
      nn<-ncol(data_matrix)
      colnames(data_matrix)[ncol(data_matrix)] <- i
    }else {
      #cat("_______________________")
    }
  }
  
  if (ncol(data_matrix) != nrow(tm_SamplesMatrix)){
    cat("ERROR: ncol data_matrix != 38!","\n")
    cat("colnames: ", colnames(data_matrix),"\n")
    cat("difference tm_samples #38: ", setdiff(colnames(data_matrix),tm_SamplesMatrix[,'tm_samples']))
    stop()
  }
  
  mutationNumberPatients<-apply(data_matrix,1,function(x) {sum(x!=0)})
  #print(ncol(data_matrix))
  timesMutated<-apply(data_matrix,1,function(x) {sum(x!=0)})
  t <- data_matrix
  
  #add additional columns to data_matrix
  variants <- rownames(data_matrix)
  genes <- c()
  mutationtype <- c()
  genesUnique <- unique(genes)
  
  position1stUDS <- gregexpr(pattern="_",variants) #UDS = underscore
  for (i in 1:length(position1stUDS)){
    positionUDS <- min(position1stUDS[[i]])
    genes <- c(genes, substr(variants[i],1,positionUDS-1))
    position4thUDS <- position1stUDS[[i]][4]
    position5thUDS <- position1stUDS[[i]][5]
    mutationtype <- c(mutationtype, substr(variants[i],position4thUDS+1,position5thUDS-1))
  }
  data_matrix <- cbind(data_matrix, genes)
  data_matrix <- cbind(data_matrix, mutationtype)
  
  mutationgroups <- unique(data_matrix[,'mutationtype'])
  
  mutationgroups <- replace(mutationgroups,mutationgroups==".","unkown")
  mutationgroups <- replace(mutationgroups,mutationgroups=="NA","unkown")

  ##
  ## create tabel for oncoprint plot -> matrixOncoPrint
  ##
  matrixOncoPrint <- matrix("",ncol = length(unique(data_matrix[,'genes'])),nrow = ncol(t))
  colnames(matrixOncoPrint) <- unique(data_matrix[,'genes'])
  rownames(matrixOncoPrint) <- colnames(t)
  
  data_matrix2 <- data_matrix[,1:(ncol(data_matrix)-2)]
  
  location_of_variances <- which(data_matrix2 != "0",arr.ind = TRUE)
  
  #                                                                                 row col
  # MEN1_64575505_C_T_nonsynonymousSNV_DOWNSTREAM_.                                  25   1
  # PIF1_65114700_C_G_nonsynonymousSNV_DOWNSTREAM_.                                  36   1
  # MEN1_64572092_C_CG_frameshiftsubstitution_DOWNSTREAM_.                           23   2

  for (i in 1:nrow(location_of_variances)){
    #cat("location_of_variances: col",location_of_variances[i,2]," row:", location_of_variances[i,1],"\n")
    #cat("data_matrix: column name:", colnames(data_matrix)[location_of_variances[i,2]], "gene: ", data_matrix[location_of_variances[i,1],'genes'], "\n")
    sample <- colnames(data_matrix)[location_of_variances[i,2]]
    gene <- data_matrix[location_of_variances[i,1],'genes']
    mutation <- data_matrix[location_of_variances[i,1],'mutationtype']
    if (mutation == "."){
      mutation = "unkown"
    }
    if (mutation == "NA"){
      mutation = "unkown"
    }
    #cat(mutation,"\n")
    if (nchar(matrixOncoPrint[sample,gene]) > 0){
      matrixOncoPrint[sample,gene] = paste(matrixOncoPrint[sample,gene],mutation,sep=";")
    }else{
      matrixOncoPrint[sample,gene]=mutation
    }
  }
  
  pdf(plotTitle)
  
  alter_fun = list(
    background = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
    },
    nonsynonymousSNV = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "blue", col = NA))
    },
    stopgain = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),gp = gpar(fill = "red", col = NA))
    },
    frameshiftsubstitution = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#008000", col = NA))
    },
    unkown = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "yellow", col = NA))
    },
    nonframeshiftsubstitution = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "orange", col = NA))
    },
    stoploss = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "purple", col = NA))
    }
  )
  
  col = c("nonsynonymousSNV" = "blue", 
          "stopgain" = "red", 
          "stoploss" = "purple", 
          "frameshiftsubstitution" = "#008000", 
          "nonframeshiftsubstitution" = "orange", 
          "unkown" = "yellow")
  
  require(lattice)
  
   h = oncoPrint(t(matrixOncoPrint), get_type = function(x) strsplit(x, ";")[[1]],
            alter_fun = alter_fun, 
            col = col, 
            column_title = plotTitle,
            remove_empty_columns = FALSE,
            show_column_names = TRUE,
            pct_gp = gpar(fontsize = 10),
            pct_digits =2,
            column_names_gp = gpar(fontsize = 5),
            heatmap_legend_param = list(title = "Alternations", at = c(mutationgroups), 
                                        labels =c(mutationgroups)))

  print(h)
  dev.off()
}

#1 Input.matrix with all normal files
normalsamples <- c("ch12","r3n","r4n")
normalstage <- c("control","control","control")
normalgrade <- c("control","control","control")
normalSamplesMatrix <- cbind(normalsamples,normalstage,normalgrade)
rm(normalsamples,normalstage,normalgrade)

#2 Input.matrix with all affected samples
tm_samples <- c("b11","b12","b13","b15","b16","b1","b23","b25","b3","b4","b5","b8","ch11","ch13","ch15","ch17","ch19","ch1","ch5","ch7","r10t",
                "r11t","r13t","r14t","r15t","r17t","r18t","r19t","r1t","r20t","r2t","r3t","r4t","r5t","r6t","r7t","r8t","r9t")
tm_stage <- c("tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","metastasis","tumor",
              "metastasis","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","metastasis","tumor","tumor","metastasis",
              "tumor","tumor","tumor","metastasis","metastasis","metastasis","metastasis","tumor","metastasis")
tm_grade <- c("G1","G1","G1","G1","G1","G1","G1","G1","G1","G1","G1","G1","G2","G1","G2","G2","G1","G2","G2","G1","G1","G1","G1","G1","G2",
              "G2","G1","NA","G1","G1","G1","G1","G2","G2","G2","G2","G1","G1")
tm_SamplesMatrix <- cbind(tm_samples,tm_stage,tm_grade)
rm(tm_samples,tm_stage,tm_grade)

#3 Input: data_matrix
files <- list.files(path=".",pattern = "*.Rda")
files <- files[ !grepl(".pdf",files) ]

#numberFiles <- length(files)

#4 Include normals TRUE/ FALSE
includeNormals = TRUE

for (y in files){
  load(y)
  print(y)
  plotTitle <- paste("oncoPrint_",y,".pdf",sep="")
  createOncoprint(normalSamplesMatrix,tm_SamplesMatrix,data_matrix,plotTitle,includeNormals)
}
