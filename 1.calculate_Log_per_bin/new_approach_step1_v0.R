library(gdata)

##################################
# 1. Get input files and sort them
##################################

#Path of input files
path_input="/home/agoepfert/documents/mesothelioma_marieke/original_input/CNV_Results_LPFG2_new/"

#find all csv files
myFiles <- list.files(path=path_input, pattern = "*.binned.50000.count.*.csv",recursive=T,include.dirs=FALSE)

#cat("________________FILE: ", myFiles)
#cat("________________FILE: ", myFiles[1])
#cat ("________________INFO: ", 
#     length(myFiles), # number of elements or components 
#     class(myFiles))  # class or type of an object

length_files=length(myFiles)
sampleInformation <- list()
sampleInformationMatrix <- matrix(NA,ncol=2)

## Find all -csv files
## Identify subfolder and filename
## saved in a list, called sampleInformation
for (i in 1:length_files){
  #cat(" _______________loop: ", myFiles[i],sep="\n")
  x <- gregexpr(pattern ='/',myFiles[i])
  #cat(" ____ position /: ", x[[1]], class(x[[1]]))
  subfolder=substring(myFiles[i],1,x[[1]]-1)
  
  if (substring(subfolder,1,1)=="e"){
    #cat(" _______________subfolder: ", subfolder, substring(subfolder,2,))
    subfolder <- substring(subfolder,2,)
    #cat(" _______________subfolder2: ", substring ,gsub("-","", substring),as.numeric(as.character(gsub("-","", substring))))
  }
  
  subfolder <- as.numeric(as.character(gsub("-","",  subfolder)))  
  
  if (subfolder<=9){
    #cat(" _______________subfolderLength: ", subfolder)
    subfolder <- sprintf("%02d", subfolder)
    #cat(" _______________subfolderLength2: ", subfolder)
  }
  
  #cat( subfolder,class(substring))
  filename=substring(myFiles[i],x[[1]]+1)
  #cat( "----", subfolder,filename,sep="\n")
  
  sampleInformation <- c(sampleInformation, list(subfolder,filename))
  x <- matrix(c(subfolder,filename),ncol=2)
  sampleInformationMatrix <- rbind(sampleInformationMatrix,x)
}


#delete first row -> NA/NA
sampleInformationMatrix <- sampleInformationMatrix[-1,]
#cat(sampleInformationMatrix)

require("BioPhysConnectoR")
sampleInformationMatrix2 <- mat.sort(sampleInformationMatrix,1)



#####################################
# 2. Calculate log per bin per sample
####################################

i=1 #19 example################################################################################################################################################################################################### 1
while (i < length(sampleInformationMatrix2)/2){
  #print (sampleInformationMatrix2[i+1,2])
  #cat(class(sampleInformationMatrix2[i+1,2]))
  
  tumor <- sampleInformationMatrix2[i,2]
  control <- sampleInformationMatrix2[i+1,2]

  controlFolderX <- gregexpr(pattern ='.binned',control)
  controlFolder <- substring(control,1,controlFolderX[[1]]-1)
  
  tumorFolderX <- gregexpr(pattern ='.binned',tumor)
  tumorFolder <- substring(tumor,1,tumorFolderX[[1]]-1)
  
  #cat(tumor,tumorFolderX[[1]], controlFolder, "_____",sep = "\n")
  #cat(control,controlFolderX[[1]],tumorFolder,"_____",sep = "\n")
  
  newTableName <- paste(tumorFolder,controlFolder,sep=".vs.")
  cat("Calculation log per bin of ", newTableName, "\n")
  
  tumorPath <-paste(path_input, tumorFolder, sep = "")
  tumorPath <-paste(tumorPath, tumor, sep = "/")
  controlPath <- paste(path_input, controlFolder, sep = "")
  controlPath <- paste(controlPath, control, sep = "/")
  
  #cat(tumorPath,sep = "\n")
  #cat(controlPath,sep = "\n")
  
  tumorFile=read.csv(tumorPath, header = TRUE, sep = ",", quote = "\"",dec = ".", fill = TRUE, comment.char = "",stringsAsFactors=FALSE)
  controlFile=read.csv(controlPath, header = TRUE, sep = ",", quote = "\"",dec = ".", fill = TRUE, comment.char = "",stringsAsFactors=FALSE)

  outputfile <- data.frame(chromosom=character(),
                           bin.start=integer(),
                           bin.end=integer(),
                           log=character(),
                           stringsAsFactors=FALSE)
  
  
  for (j in 1:nrow(tumorFile)){
    #ignore sex chromosoms
    if(tumorFile[j,"CHR"]=="chrX" | tumorFile[j,"CHR"]=="chrY"){
      #print(tumorFile[j,"CHR"])
    } 
    else {
      chromosom <- tumorFile[j,"CHR"]
      startBin <- tumorFile[j,"BIN.START"]
      startBinControl <- controlFile[j,"BIN.START"]
      endBin <- tumorFile[j,"BIN.END"]
      if (startBin!=startBinControl){
        cat("Warning: bin unequal __startBinTumor: ", startBin, " _____startBinControl: ", startBinControl,"\n")
      }
      else{
        #cat("___startBinTumor: ", startBin, " _____startBinControl: ", startBinControl,"\n")
        logValue <- log2(tumorFile[j,14]/controlFile[j,14])
        #cat(chromosom, startBin, endBin, logValue,"\n")
        outputfile[nrow(outputfile)+1, ] <- c(chromosom, startBin, endBin, logValue)
      }
    }
  }  

  outputfile_name <- paste(newTableName,".txt",sep="")
  write.table(outputfile,file=outputfile_name,quote=FALSE,row.names = FALSE)
  i=i+2
}

