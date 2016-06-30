setwd("/home/agoepfert/documents/pnet/38PNETs/Heatmap/20160527_matrix_v3/createPlots")

source("https://bioconductor.org/biocLite.R")
#library(circlize)
#biocLite("ComplexHeatmap") #install ComplexHeatmap
#biocLite()
library(gplots)
#complexHeatmap
library(ComplexHeatmap)
library(circlize)

createHeatmap <- function(normalSamplesMatrix,tm_SamplesMatrix,data_matrix,plotTitle,includeNormals){

colnames(data_matrix) <- gsub("_VSC ","", colnames(data_matrix))
#remove normal samples if includeNormals==FALSE, otherwise add them to the matrix
if (includeNormals==FALSE){
  data_matrix <- data_matrix[ , !colnames(data_matrix) %in% c(normalSamplesMatrix[,'normalsamples'])]
} else{
  tm_SamplesMatrix <- rbind(tm_SamplesMatrix,normalSamplesMatrix)
}

rownames(tm_SamplesMatrix) <- tm_SamplesMatrix[,'tm_samples']
tm_SamplesMatrix <- tm_SamplesMatrix[order(rownames(tm_SamplesMatrix)), ] 


for (i in tm_SamplesMatrix[,'tm_samples']){
  #cat(i)
  if (is.na(x = match(i,colnames(data_matrix)))==TRUE){
    #  if (match(i,colnames(data_matrix))==NA){
    #cat(i,"\n") 
    tmp_matrix <- rep.int(0,nrow(data_matrix))
    data_matrix <- cbind(data_matrix, tmp_matrix)
    nn<-ncol(data_matrix)
    colnames(data_matrix)[ncol(data_matrix)] <- i
  } else {
    #cat("_______________________")
  }
}

if (ncol(data_matrix) != nrow(tm_SamplesMatrix)){
  cat("ERROR: ncol data_matrix != 38!","\n")
  cat("colnames: ", colnames(data_matrix),"\n")
  cat("difference tm_samples #38: ", setdiff(colnames(data_matrix),tm_SamplesMatrix[,'tm_samples']))
  stop()
}

timesMutated<-apply(data_matrix,1,function(x) {sum(x!=0)})
data_matrix <- cbind(data_matrix, timesMutated)


data_matrix <- as.data.frame(data_matrix[timesMutated > 0,])

df = data.frame(stage = tm_SamplesMatrix[,'tm_stage'], grade = tm_SamplesMatrix[,'tm_grade'])

###### STAGE
colfunc <- colorRampPalette(c( "steelblue4","mediumseagreen","orangered3"))
listOfColors <- colfunc(length(unique(tm_SamplesMatrix[,'tm_stage'])))
stagelist  <- listOfColors
names(stagelist) <- unique(tm_SamplesMatrix[,'tm_stage'])

###### STAGE
colfunc <- colorRampPalette(c( "lightsteelblue1","royalblue4","mediumseagreen","grey94"))
listOfColors <- colfunc(length(unique(tm_SamplesMatrix[,'tm_grade'])))
gradelist  <- listOfColors
names(gradelist) <- unique(tm_SamplesMatrix[,'tm_grade'])

ha <- columnAnnotation(df = df,gap = unit(c(1, 1),"mm"), 
                       col = list(stage = stagelist, grade= gradelist),
                       annotation_legend_param = list(
                         stage=list(title="stage",title_gp = gpar(fontsize = 10),
                                    labels_gp = gpar(fontsize = 8),
                                    lables= c(stagelist), 
                                    at= c(names(stagelist))),
                         grade=list(title="grade",title_gp = gpar(fontsize = 10),
                                    labels_gp = gpar(fontsize = 8),
                                    lables= c(gradelist), 
                                    at= c(names(gradelist)))
                       )
)


#add additional columns to data_matrix
#calculate mutationBurdon

variants <- rownames(data_matrix)
genes <- c()
mutationtype <- c()

position1stUDS <- gregexpr(pattern="_",variants) #UDS = underscore
for (i in 1:length(position1stUDS)){
  positionUDS <- min(position1stUDS[[i]])
  genes <- c(genes, substr(variants[i],1,positionUDS-1))
  position4thUDS <- position1stUDS[[i]][4]
  position5thUDS <- position1stUDS[[i]][5]
  mutationtype <- c(mutationtype, substr(variants[i],position4thUDS+1,position5thUDS-1))
  ##cat(mutationtype)
}

#add genenames and mutationtype, mutationBurdon
data_matrix <- cbind(data_matrix, genes)
data_matrix <- cbind(data_matrix, mutationtype)

remove(timesMutated)
remove(genes)

#calculate mutation burdon per gene
#http://stackoverflow.com/questions/16242367/sum-up-specific-rows-based-on-name-in-r-and-create-new-column

data_matrix <- data.frame(data_matrix,stringsAsFactors=TRUE)
names(data_matrix$genes) <-"NULL"
names(data_matrix$timesMutated) <-"NULL"

data_matrix$timesMutated <- as.numeric(as.character(data_matrix$timesMutated))
data_matrix<- transform(data_matrix, mutation_per_gene = ave(timesMutated,genes, FUN = sum))

data_matrix <- data_matrix[order(-data_matrix$mutation_per_gene,data_matrix$genes),]

t <- data_matrix[,c(1:(ncol(data_matrix)-5))]
tolower(colnames(t))

##Order the variants per gene either by their max. allelic fraction or by sum over the whole row
data_matrix$maxAlleleFrequency <- apply(t, 1, max) 
data_matrix <- data_matrix[order(-data_matrix$mutation_per_gene,data_matrix$genes,-data_matrix$maxAlleleFrequency),]

data_matrix$sumAlleleFrequency <- apply(t, 1, sum) 
#data_matrix <- data_matrix[order(-data_matrix$mutation_per_gene,data_matrix$genes,-data_matrix$sumAlleleFrequency),]

t <- data_matrix[,c(1:(ncol(data_matrix)-6))]
t <- t[, order(colnames(t))] 

pdf(plotTitle)


##Add Annotation about samples
dfGENES = data.frame(anno_genes = data_matrix$genes)

colfunc <- colorRampPalette(c("blanchedalmond","burlywood1","burlywood3","lightsalmon2","indianred1","indianred3","orangered3","maroon","mistyrose4","lightslategray","royalblue1","steelblue1","paleturquoise1","paleturquoise3"))
listOfColors <- colfunc(length(unique(data_matrix$genes)))
genelist  <- listOfColors
names(genelist) <- unique(data_matrix$genes)

annotation_genes = rowAnnotation(df = dfGENES, width = unit(1,"cm"),gap = unit(0,"mm"),
                                 col = list(anno_genes=genelist),
                                 annotation_legend_param = list(
                                   anno_genes = list(title ="genes",title_gp = gpar(fontsize = 10),
                                                     labels_gp = gpar(fontsize = 8),
                                                     lables= c(genelist), 
                                                     at= c(names(genelist)))))

# ##Add Annotation for genes
# dfGENES = data.frame(anno_genes = data_matrix$genes)
# annotation_genes = rowAnnotation(df = dfGENES, width = unit(1,"cm"))

heatmap= Heatmap(scale(t, center=FALSE, scale=FALSE),na_col ="lightsteelblue1", name="Allelic fraction",
                 col = colorRamp2(c(0,0.05,0.1,0.4,1), c("aliceblue","deepskyblue4" ,"mediumaquamarine","orange","red")),
                 row_title ="", 
                 row_title_side ="right",
                 row_title_rot = 0,
                 column_title ="samples",column_title_side ="bottom",
                 #column_dend_side ="bottom",
                 column_dend_height = unit(2,"cm"),
                 row_dend_width = unit(2,"cm"),
                 column_names_gp = gpar(fontsize = 5),
                 bottom_annotation = ha,
                 row_names_gp = gpar(fontsize = 1),
                 show_row_names = FALSE,
                 column_title_gp = gpar(fontsize = 10),
                 row_title_gp = gpar(fontsize = 5),
                 #km=kmeans(t, 4),
                 ##split=data_matrix$genes, #WRONG ORDER
                 ## check data_matrix$genes and c(data_matrix$genes)
                 split = factor(data_matrix$genes,levels(data_matrix$genes)[unique(c(data_matrix$genes))]),
                 #combined_name_fun = NULL, # doesn't change font size
                 #heatmap_legend_param = list(title_gp = gpar(fontsize = 10),labels_gp = gpar(fontsize = 8)),
                 cluster_rows = FALSE,
                 #gap = unit(c(10, 5),"mm"),
                 cluster_columns = FALSE,
                 #clustering_distance_columns = ("euclidean"),
                 show_row_dend = FALSE,
                 heatmap_legend_param = list(color_bar ="continuous", legend_direction ="vertical",
                                             legend_width = unit(5,"cm"))
)

draw(heatmap, heatmap_legend_side ="left",annotation_legend_side ="left")
#heatmap#+annotation_genes

for(an in colnames(df)) {
  decorate_annotation(an, {
    ## annotation names on the right
    grid.text(an, unit(1,"npc") + unit(0.5,"mm"), 0.5, default.units ="npc", just ="left",gp=gpar(fontsize=10))
    ## annotation names on the left
    #grid.text(an, unit(0,"npc") - unit(2,"mm"), 0.5, default.units ="npc", just ="right")
  })
}


dev.off()

#print(h)

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
  plotTitle <- paste("heatmap_",y,".pdf",sep="")
  createHeatmap(normalSamplesMatrix,tm_SamplesMatrix,data_matrix,plotTitle,includeNormals)
}
