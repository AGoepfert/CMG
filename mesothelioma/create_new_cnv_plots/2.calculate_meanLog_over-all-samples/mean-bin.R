library(gdata)

##################################
# 1. Get input files and sort them
##################################
cat("Step 1: Load input files \n")

#Path of input files
path_input <- "/home/agoepfert/documents/mesothelioma_marieke/new_approach/LPWG2/1.calculate_Log_per_bin/"
path_output <- "/home/agoepfert/documents/mesothelioma_marieke/new_approach/LPWG2/2.calculate_meanLog_over-all-samples/"

#find all input files
myFiles <- list.files(path=path_input, pattern = "*.vs.*.txt")

# cat("________________FILE: ", myFiles)
# cat("________________FILE: ", myFiles[1])
# cat ("________________INFO: ", 
#      length(myFiles), # number of elements or components 
#      class(myFiles))  # class or type of an object

length_files=length(myFiles)

#create a matrix - important that the rom number is already given - otherwise no cbind is possible
file <- paste(path_input,myFiles[1],sep="")
default_inputFile <- read.table(file, header = TRUE, sep = "", quote = "\"",dec = ".", fill = TRUE, comment.char = "",stringsAsFactors=FALSE)
default_length <- length(default_inputFile$log)
logMatrix <- matrix(nrow=length(default_inputFile$log), ncol=0) 

cat("Step 2: Combine all log information in one matrix \n")
##############################################
# 2. Combine all log-information in one matrix
##############################################
for (i in 1:length_files){
  cat("  ", myFiles[i],"\n")
  file <- paste(path_input,myFiles[i],sep="")
  inputFile <- read.table(file, header = TRUE, sep = "", quote = "\"",dec = ".", fill = TRUE, comment.char = "",stringsAsFactors=FALSE)
  
  if (default_length != length(inputFile$log)){
    cat("   WARNING: length in files differs from each other (",default_length,"vs",length(inputFile$log),")")
  } 
  else{
    # Check if start.bin in two files is the same
    if (identical(default_inputFile$bin.start,inputFile$bin.start)){
          cat("    -> Length log column:",length(inputFile$log),"\n")
          logMatrix <- cbind.data.frame(logMatrix,inputFile$log)
          colnames(logMatrix)[colnames(logMatrix)=="inputFile$log"] <- myFiles[i]
          
    }
    else{
      cat("   WARNING: bin.start differs between default file and ",myFiles[i],"\n")
    }
  }
}


####################################
# 3. Calculate mean of all the logs
####################################
cat("Step 3. Calculate mean of all the logs \n")


#Problem: NA und Inf -> if 1 occurs in a row, it get that value
#solution: set NA and Inf/-Inf to zero
#logMatrixNACor <- logMatrix
removeRows <- c()
lmn <- as.numeric(ncol(logMatrix))

for (i in 1:nrow(logMatrix)) {
  #message(length(which(is.na(logMatrix[i,]))))
          
  if (length(which(is.na(logMatrix[i,]))) == lmn) {
    #message(paste("Remove row : ",i))
    removeRows <- c(removeRows,i)
  }
}


if (length(removeRows) > 0) {
  logMatrixNACor <- logMatrix[-removeRows,]
} else {
  logMatrixNACor <- logMatrix
}



logMatrixNAInfiniteCor <- as.matrix(logMatrixNACor)
cat("  test - format of matrices: ",class(logMatrixNAInfiniteCor), class(logMatrixNACor), class(logMatrix))
logMatrixNAInfiniteCor[is.infinite(logMatrixNAInfiniteCor)] <- 0


#generating the outputfile
outputfile <- data.frame(chromosom=character(),
                         bin.start=integer(),
                         bin.end=integer(),
                         log.overAllSamples.mean=integer(),
                         stringsAsFactors=FALSE)

#Problem: NA und Inf -> if 1 occurs in a row, it get that value
#solution A: ignore NA and Inf
#solution B: set NA and Inf/-Inf to zero
if (length(removeRows) > 0) {
  default_inputFile <- default_inputFile[-removeRows,]

} 
  

  for (i in 1:nrow(logMatrixNAInfiniteCor)){
    rowLog <- logMatrixNAInfiniteCor[i,]
    #print(rowLog)
    x <- mean(as.numeric(logMatrixNAInfiniteCor[i,1:length_files]),na.rm=T)
    #print(x)
    #print("_____________________________________________________")
    chromosom <- default_inputFile[i,1]
    startBin <- default_inputFile[i,2]
    endBin <- default_inputFile[i,3]
    logValue <- x
    outputfile[nrow(outputfile)+1, ] <- c(chromosom, startBin, endBin, logValue)
  }

filename <- paste(path_output,"MeanOverAllSamples.txt",sep="")
write.table(outputfile,file=filename,quote=FALSE,row.names = FALSE)
write.table(default_inputFile[,1:2],file=paste(filename,'.pos',sep=''),row.names=F,quote=F)
write.table(logMatrixNAInfiniteCor,file=paste(filename,'.all.logR',sep=''),row.names=F,quote=F)
