# GET arguments
# args <- commandArgs(TRUE)
# sample1 <- args[1] #name sample 1
# sample2 <- args[2] #name sample 2
# file1 <- args[3] # csv file tumor 'CHR','BIN.START','BIN.END','BIN.GC.CONTENT','BIN.N.COUNT' -> *.binned.50000.count.filtered.reads.csv
# file2 <- args[4] #csv file normal 'CHR','BIN.START','BIN.END','BIN.GC.CONTENT','BIN.N.COUNT'
# out_file <- args[5]  #output file

#QGP WES
setwd("/home/agoepfert/documents/pnet/cell_lines/WES_new2015/cnv/qgp/")
sample1 <- "qgpresist4" #name sample 1
sample2 <- "qgpsens3" #name sample 2
file1 <- "/home/agoepfert/documents/pnet/cell_lines/WES_new2015/cnv/qgp/qgpresist4.binned.50000.count.filtered.reads.csv" # csv file tumor 'CHR','BIN.START','BIN.END','BIN.GC.CONTENT','BIN.N.COUNT' -> *.binned.50000.count.filtered.reads.csv
file2 <- "/home/agoepfert/documents/pnet/cell_lines/WES_new2015/cnv/qgp/qgpsens3.binned.50000.count.filtered.reads.csv" #csv file normal 'CHR','BIN.START','BIN.END','BIN.GC.CONTENT','BIN.N.COUNT'
out_file <- "CNV_qgpresist4.vs.qgpsens3"  #output file

#BON WES
# setwd("/home/agoepfert/documents/pnet/cell_lines/WES_new2015/cnv/bon/")
# sample1 <- "bonresistent2" #name sample 1
# sample2 <- "bonbonsens1" #name sample 2
# file1 <- "/home/agoepfert/documents/pnet/cell_lines/WES_new2015/cnv/bon/bonresist2.binned.50000.count.filtered.reads.csv" # csv file tumor 'CHR','BIN.START','BIN.END','BIN.GC.CONTENT','BIN.N.COUNT' -> *.binned.50000.count.filtered.reads.csv
# file2 <- "/home/agoepfert/documents/pnet/cell_lines/WES_new2015/cnv/bon/bonsens1.binned.50000.count.filtered.reads.csv" #csv file normal 'CHR','BIN.START','BIN.END','BIN.GC.CONTENT','BIN.N.COUNT'
# out_file <- "CNV_bonresistent2.vs.bonbonsens1"  #output file

# ##QGP WGS
# setwd("/home/agoepfert/documents/pnet/cell_lines/CNV_calling_BON-QGP_WGS")
# sample1 <- "qgp_r" #name sample 1
# sample2 <- "qgp_s" #name sample 2
# file1 <- "/home/agoepfert/documents/pnet/cell_lines/CNV_calling_BON-QGP_WGS/qgp_r.binned.50000.count.filtered.reads.csv" # csv file tumor 'CHR','BIN.START','BIN.END','BIN.GC.CONTENT','BIN.N.COUNT' -> *.binned.50000.count.filtered.reads.csv
# file2 <- "/home/agoepfert/documents/pnet/cell_lines/CNV_calling_BON-QGP_WGS/qgp_s.binned.50000.count.filtered.reads.csv" #csv file normal 'CHR','BIN.START','BIN.END','BIN.GC.CONTENT','BIN.N.COUNT'
# out_file <- "CNV_qgp_evo.vs.qgp_dmso"  #output file

#BON WGS
# setwd("/home/agoepfert/documents/pnet/cell_lines/CNV_calling_BON-QGP_WGS")
# sample1 <- "bon_r" #name sample 1
# sample2 <- "bon_s" #name sample 2
# file1 <- "/home/agoepfert/documents/pnet/cell_lines/CNV_calling_BON-QGP_WGS/bon_r.binned.50000.count.filtered.reads.csv" # csv file tumor 'CHR','BIN.START','BIN.END','BIN.GC.CONTENT','BIN.N.COUNT' -> *.binned.50000.count.filtered.reads.csv
# file2 <- "/home/agoepfert/documents/pnet/cell_lines/CNV_calling_BON-QGP_WGS/bon_s.binned.50000.count.filtered.reads.csv" #csv file normal 'CHR','BIN.START','BIN.END','BIN.GC.CONTENT','BIN.N.COUNT'
# out_file <- "CNV_bon_evo.vs.bon_vehikel"  #output file

library(DNAcopy)
library("zoo")
options(scipen=999)
options(bitmapType="cairo")
data1 <- read.csv(file1,header=TRUE)
data2 <- read.csv(file2,header=TRUE)


blacklist <- '/home/agoepfert/documents/01_general_datasets/common_cnv/CNV_common_NIPT.bed'
genelist <- '/home/agoepfert/documents/pnet/cell_lines/gene_list/mTor_38Pnets_pathcrds.genecards.txt'

# read the CNV blacklist
common_cnv <- read.table(blacklist,header=FALSE)
colnames(common_cnv) <- c('chr','start','stop')
common_cnv$chr <- sub('X',23,common_cnv$chr)
common_cnv$chr <- sub('Y',24,common_cnv$chr)
previousSetDown=0
# read the gene list
interesting_genelist <- read.table(genelist,header=FALSE)
colnames(interesting_genelist) <- c('chr','start','stop','geneName')
interesting_genelist$chr <- sub('X',23,interesting_genelist$chr)
interesting_genelist$chr <- sub('Y',24,interesting_genelist$chr)

# merge and remove zero
data <- merge(data1,data2,all.x=TRUE,by=c('CHR','BIN.START','BIN.END','BIN.GC.CONTENT','BIN.N.COUNT'))
data$CHR <- sub('chr','',data$CHR)

#rm <- which(data[,18] == 0 )
#data <- data[-rm,]
#rm <- which(data[,7] == 0 )
#data <- data[-rm,]
#data[56613,]
# take ratio.
data$LogR <- log2(data[,15]/data[,27])

data[data$CHR == 'X','LogR'] <- log2(data[data$CHR=='X',16]/data[data$CHR=='X',28])
data[data$CHR == 'Y','LogR'] <- log2(data[data$CHR=='Y',17]/data[data$CHR=='Y',29])
data$CHR <- sub('X',23,data$CHR)
data$CHR <- sub('Y',24,data$CHR)
rm <- which(is.na(data$LogR))
data <- data[-rm,]
rm <- which(is.infinite(data$LogR))
data <- data[-rm,]
# normalize by median.
data$LogR <- data$LogR - median(data$LogR)

#segment
CNA.object <- CNA(data$LogR,data$CHR,data$BIN.START,sampleid=paste(sample1,'/',sample2,sep=''),data.type='logratio')
smoothed.CNA.object <- smooth.CNA(CNA.object)
logR.sd = sd(data$LogR)
# round one segment: no undo
round_one <- segment(smoothed.CNA.object,verbose=0)
nr_seg <- nrow(round_one$output)
undo_thresh <- 0.02
undo_type <- 'prune'
if (logR.sd > 1 || nr_seg > 250) {
	message('high logR.SD or many raw segments');
	undo_thresh <- 2
	undo_type <- 'SD'
	prune.CNA.object <- segment(smoothed.CNA.object,undo.splits='sdundo',undo.SD=undo_thresh,verbose=2)
} else {
	prune.CNA.object <- segment(smoothed.CNA.object,undo.splits='prune',undo.prune=undo_thresh,verbose=2)
}
#write full table
calls <- prune.CNA.object$output
write.table(calls,file=paste(out_file,'.full_list',sep=''),quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
# get the relevant onces
# minimal logR deviation
cnvs <- calls[(calls$seg.mean >= 0.25 | calls$seg.mean <= -0.25),]
# minimal nr of probes
cnvs <- cnvs[cnvs$num.mark >= 5,]
## in common?
cnvs$comments = rep('',nrow(cnvs))
if (nrow(cnvs) > 0) {
   for (i in 1:nrow(cnvs)) {
	chr = cnvs[i,'chrom']
	subcommon <- common_cnv[common_cnv$chr == chr,]
	if (nrow(subcommon) == 0) {
		next
	}
	for (j in 1:nrow(subcommon)) {
		## only set to common if completely inside common region.
		if ((cnvs$loc.start[i] >= subcommon$start[j]) && (cnvs$loc.end[i] <= subcommon$stop[j])) {
			cnvs$comments[i] <- 'Common CNV'
			break
		}
		## mention overlap (start inside common)
		if ((cnvs$loc.start[i] >= subcommon$start[j]) && (cnvs$loc.start[i] <= subcommon$stop[j]) && (cnvs$loc.end[i] > subcommon$stop[j])) {
			frac_in_common <- round((subcommon$stop[j] - cnvs$loc.start[i] + 1) / (cnvs$loc.end[i] - cnvs$loc.start[i] + 1),2)*100
			cnvs$comments[i] <- paste(frac_in_common,'% at 5prime end in common CNV region',sep='')
		}
		## mention overlap (stop inside common)
		if ((cnvs$loc.start[i] < subcommon$start[j]) && (cnvs$loc.end[i] <= subcommon$stop[j]) && (cnvs$loc.end[i] >= subcommon$start[j])) {
			frac_in_common <- round((cnvs$loc.end[i] - subcommon$start[j] + 1) / (cnvs$loc.end[i] - cnvs$loc.start[i] + 1),2)*100
			cnvs$comments[i] <- paste(frac_in_common,'% at 3prime end in common CNV region',sep='')
		}
		## mention overlap (extends both sides)
		if ((cnvs$loc.start[i] < subcommon$start[j]) && (cnvs$loc.end[i] > subcommon$stop[j]) ) {
			frac_in_common <- round((subcommon$stop[j] - subcommon$start[j] + 1) / (cnvs$loc.end[i] - cnvs$loc.start[i] + 1),2)*100
			cnvs$comments[i] <- paste(frac_in_common,'% center part in common CNV region',sep='')
		}
	}

   }
}
# write out the cnvs
write.table(cnvs,file=out_file,quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)



#cnvs <- data.frame(chrom=c())
# plot
#i=1
for (i in 1:12) {
 	png(file=paste(sample1,".vs.",sample2,".page_",i,".png",sep=""), bg="white", width=1060, height=800)
  cat(paste("/tmp/",sample1,".vs.",sample2,".page_",i,".png",sep=""),"\n")
	par(mfrow=c(2,1))
	# first chr
	chr = 2*i -1
	message(paste('plotting chromosome: ',chr,sep=''))
	subcnv <- cnvs[cnvs$chrom == chr,]
	subcommon <- common_cnv[common_cnv$chr == chr,]
	subcommon_genelist <- interesting_genelist[interesting_genelist$chr == chr,]
	
	subdata <- data[data$CHR == chr,]
	subdata <- subdata[order(subdata$BIN.START),]
	rows <-nrow(subdata)
	if (rows > 0) {
		plot(subdata$BIN.START/1000000,subdata$LogR, main=paste(sample1,"/",sample2, " : LogR for Chr", chr, sep=" "),xlab="Position (Mb)", ylab="Log2(R)", pch=".",ylim=c(-4,2),xaxs="i", yaxs="i")
		grid(nx=NA,ny=NULL,col='red',lty='solid',lwd=1)
		cat("plot", chr,"\n")
		if (rows > 75) {
			smoothval <- c()
			smoothval <- rollmean(subdata$LogR,51)
			left <- c()
			for (j in 1:25) {
				left[j] <- mean(subdata[0:(j+25),"LogR"],trim=0.1) 
			}
			j <- 1
			startidx <- rows-24
			right <- c()
			for (k in startidx:rows	) {
				right[j] <- mean(subdata[(k-25):rows,"LogR"],trim=0.1)
			j <- j+1
			}
			smoothval <- c(left,smoothval,right)
			lines(subdata$BIN.START/1000000, smoothval, col='blue')	
			#lines(subdata$BIN.START/1000000,smooth(subdata$LogR),col='blue')
		}
	} else {
		message('no data')
	  cat("no data", chr,"\n")
		plot(c(0,1),c(0,0), main=paste(sample1,"/",sample2, " : LogR for Chr", chr, sep=" "),xlab="Position (Mb)", ylab="Log2(R)", pch="",ylim=c(-4,2),xaxs="i", yaxs="i")
		grid(nx=NA,ny=NULL,col='red',lty='solid',lwd=1)	
	}
	# sample cnvs
	if (nrow(subcnv)>0) {
		for (j in 1:nrow(subcnv)) {
			# the box hightlight
			xs <- c(subcnv$loc.start[j]/1000000, subcnv$loc.end[j]/1000000, subcnv$loc.end[j]/1000000, subcnv$loc.start[j]/1000000)
			ys <- c(-4,-4,2,2)
			if (subcnv$seg.mean[j] < -0.25) {
				polygon(xs, ys, col=rgb(160,32,240,100,maxColorValue=255),border=NA)
				abline(v=xs[1],lwd=1,col=rgb(160,32,240,50,maxColorValue=255))
				abline(v=xs[2],lwd=1,col=rgb(160,32,240,50,maxColorValue=255))

			}
			if (subcnv$seg.mean[j] > 0.25) {
				polygon(xs, ys, col=rgb(0,128,0,100,maxColorValue=255),border=NA)
				abline(v=xs[1],lwd=1,col=rgb(0,128,0,50,maxColorValue=255))
				abline(v=xs[2],lwd=1,col=rgb(0,128,0,50,maxColorValue=255))
			}
			# the mean logR
			y <- subcnv$seg.mean[j]
			segments(subcnv$loc.start[j]/1000000,y, subcnv$loc.end[j]/1000000,y,col='red',lwd=2)
			
		} 
	}
	#  plot the common cnvs
	if (nrow(subcommon)>0) {
		for (j in 1:nrow(subcommon)) {
			
			xs <- c(subcommon$start[j]/1000000, subcommon$stop[j]/1000000, subcommon$stop[j]/1000000, subcommon$start[j]/1000000)
			ys <- c(-3.95,-3.95,-3.8,-3.8)
			
			polygon(xs, ys, col=rgb(75,75,75,250,maxColorValue=255),border=NA)
		} 
	}
	
	# plot interesting genes
	if (nrow(subcommon_genelist)>0) {
	  if (subcommon_genelist$chr[1]==16 | subcommon_genelist$chr[1]==17 | subcommon_genelist$chr[1]==18 | subcommon_genelist$chr[1]==19){
	  }
	  else if (subcommon_genelist$chr[1]==20){
	    xx=0.75
	  }
	  else if (subcommon_genelist$chr[1]==21 | subcommon_genelist$chr[1]==22){
	    xx=0.45
	  }
	  else if (subcommon_genelist$chr[1]==1 | subcommon_genelist$chr[1]==2){
	    xx=2.8
	  }
	  else if (subcommon_genelist$chr[1]==8 | subcommon_genelist$chr[1]==9 | subcommon_genelist$chr[1]==10 | subcommon_genelist$chr[1]==11 | subcommon_genelist$chr[1]==12){
	    xx=1.6
	  }
	  else if (subcommon_genelist$chr[1]==13 | subcommon_genelist$chr[1]==14 | subcommon_genelist$chr[1]==15){
	    xx=1.1
	  }
	  else  {
	    xx=2.2
	  }
	  
	  length_genelist=nrow(subcommon_genelist)
	  for (w in 1:(length_genelist)){
	    if (w==length_genelist){
	      subcommon_genelist$geneDistance[w] <- NA
	    }
	    else {
	      subcommon_genelist$geneDistance[w] <- subcommon_genelist$start[w+1] - subcommon_genelist$stop[w] 
	    }
	  }
	  print(subcommon_genelist)
	  previousSetDown=FALSE
	  for (j in 1:nrow(subcommon_genelist)) {
	    print(j)
	    xs  <- c(subcommon_genelist$start[j]/1000000, subcommon_genelist$stop[j]/1000000+0.2, subcommon_genelist$stop[j]/1000000+0.2, subcommon_genelist$start[j]/1000000)
	    ys <- c(-0.95,-0.95,-0.8,-0.8)
	    
	    polygon(xs, ys, col="red",border=NA)
	    
	    if(is.na(subcommon_genelist$geneDistance[j])){
	      yy=-1
	      print("test")
	      previousSetDown=FALSE
	    }else if (subcommon_genelist$chr[1]==20 | subcommon_genelist$chr[1]==21 | subcommon_genelist$chr[1]==22){
	      yy=-1
	      previousSetDown=FALSE
	    }else if (subcommon_genelist$chr[1]==1 | subcommon_genelist$chr[1]==2 | subcommon_genelist$chr[1]==3 | subcommon_genelist$chr[1]==4 | subcommon_genelist$chr[1]==5 | subcommon_genelist$chr[1]==6 | subcommon_genelist$chr[1]==7){
	      if(subcommon_genelist$geneDistance[j] < 1500000){
	        if(previousSetDown==FALSE){
	          yy=-2.1
	          previousSetDown=TRUE
	        }else{
	          yy=-0.2
	          previousSetDown=FALSE
	        }
	      }else{
	        yy=-1
	        previousSetDown=FALSE
	      }
	    }else if(subcommon_genelist$geneDistance[j] < 500000){
	      if(previousSetDown==FALSE){
	        yy=-2.1
	        previousSetDown=TRUE
	      }else{
	        yy=-0.2
	        previousSetDown=FALSE
	      }
	    }else{
	      yy=-1
	      previousSetDown=FALSE
	    }
	    print("hey")
	    #print(subcommon_genelist$geneDistance[j])
	    #print(typeof(subcommon_genelist$geneDistance[j]))
	    #yy=-1
	    text(subcommon_genelist$start[j]/1000000+xx, yy, pos=2, labels = subcommon_genelist$geneName[j],col="red",cex=0.7,srt=90)
	  } 
	}
	
	
	
	
	## add legend
	if (i == 1) {
		labs <- paste("log(S1/S2).SD: ",round(logR.sd,2),"\nS1: ",sample1," median(count): ",median(data1[,7]),"\nS2:",sample2," median(count): ",median(data2[,7]),"\n","undo.segment by ",undo_type,", threshold: ",undo_thresh)
		x <- max(subdata$BIN.START)/1000000
		y <- -3.45 
		text(x,y,labels = labs,pos=2,cex=0.9)
	}
	# second chr
	chr = 2*i
	message(paste('plotting chromosome: ',chr,sep=''))
	subcnv <- cnvs[cnvs$chrom == chr,]
	subdata <- data[data$CHR == chr,]
	subcommon <- common_cnv[common_cnv$chr == chr,]
	subdata <- subdata[order(subdata$BIN.START),]
	subcommon_genelist <- interesting_genelist[interesting_genelist$chr == chr,]
	rows <-nrow(subdata)
	if (rows > 0) {
		plot(subdata$BIN.START/1000000,subdata$LogR, main=paste(sample1,"/",sample2, " : LogR for Chr", chr, sep=" "),xlab="Position (Mb)", ylab="Log2(R)", pch=".",ylim=c(-4,2),xaxs="i", yaxs="i")
		grid(nx=NA,ny=NULL,col='red',lty='solid',lwd=1)
		if (rows > 75) {
			smoothval <- c()
			smoothval <- rollmean(subdata$LogR,51)
			left <- c()
			for (j in 1:25) {
				left[j] <- mean(subdata[0:(j+25),"LogR"],trim=0.1) 
			}
			j <- 1
			startidx <- rows-24
			right <- c()
			for (k in startidx:rows	) {
				right[j] <- mean(subdata[(k-25):rows,"LogR"],trim=0.1)
				j <- j+1
			}
			smoothval <- c(left,smoothval,right)
			lines(subdata$BIN.START/1000000, smoothval, col='blue')	
		}
	 } else {
		plot(c(0,1),c(0,0), main=paste(sample1,"/",sample2, " : LogR for Chr", chr, sep=" "),xlab="Position (Mb)", ylab="Log2(R)", pch="",ylim=c(-4,2),xaxs="i", yaxs="i")
		grid(nx=NA,ny=NULL,col='red',lty='solid',lwd=1)	
	 }
	 rows <- nrow(subcnv)
	 if (rows > 0) {
		for (j in 1:nrow(subcnv)) {
			
			xs <- c(subcnv$loc.start[j]/1000000, subcnv$loc.end[j]/1000000, subcnv$loc.end[j]/1000000, subcnv$loc.start[j]/1000000)
			ys <- c(-4,-4,2,2)
			
			if (subcnv$seg.mean[j] < -0.25) {
				polygon(xs, ys, col=rgb(160,32,240,100,maxColorValue=255),border=NA)
				abline(v=xs[1],lwd=1,col=rgb(160,32,240,50,maxColorValue=255))
				abline(v=xs[2],lwd=1,col=rgb(160,32,240,50,maxColorValue=255))

			}
			if (subcnv$seg.mean[j] > 0.25) {
				polygon(xs, ys, col=rgb(0,128,0,100,maxColorValue=255),border=NA)
				abline(v=xs[1],lwd=1,col=rgb(0,128,0,50,maxColorValue=255))
				abline(v=xs[2],lwd=1,col=rgb(0,128,0,50,maxColorValue=255))
			}
			# the mean logR
			y <- subcnv$seg.mean[j]
			segments(subcnv$loc.start[j]/1000000,y, subcnv$loc.end[j]/1000000,y,col='red',lwd=2)

		} 
	}
	rows <- nrow(subcommon)
	if (rows > 0) {
		for (j in 1:nrow(subcommon)) {
			
			xs <- c(subcommon$start[j]/1000000, subcommon$stop[j]/1000000, subcommon$stop[j]/1000000, subcommon$start[j]/1000000)
			ys <- c(-3.95,-3.95,-3.8,-3.8)
			
			polygon(xs, ys, col=rgb(75,75,75,250,maxColorValue=255),border=NA)
			#abline(v=xs[1],lwd=1,col=rgb(200,200,200,250,maxColorValue=255))
			#abline(v=xs[2],lwd=1,col=rgb(200,200,200,250,maxColorValue=255))
		} 
	}
	
	# plot interesting genes
	if (nrow(subcommon_genelist)>0) {
	  if (subcommon_genelist$chr[1]==16 | subcommon_genelist$chr[1]==17 | subcommon_genelist$chr[1]==18 | subcommon_genelist$chr[1]==19){
	  }
	  else if (subcommon_genelist$chr[1]==20){
	    xx=0.75
	  }
	  else if (subcommon_genelist$chr[1]==21 | subcommon_genelist$chr[1]==22){
	    xx=0.45
	  }
	  else if (subcommon_genelist$chr[1]==1 | subcommon_genelist$chr[1]==2){
	    xx=2.8
	  }
	  else if (subcommon_genelist$chr[1]==8 | subcommon_genelist$chr[1]==9 | subcommon_genelist$chr[1]==10 | subcommon_genelist$chr[1]==11 | subcommon_genelist$chr[1]==12){
	    xx=1.6
	  }
	  else if (subcommon_genelist$chr[1]==13 | subcommon_genelist$chr[1]==14 | subcommon_genelist$chr[1]==15){
	    xx=1.1
	  }
	  else  {
	    xx=2.2
	  }
	  
	  length_genelist=nrow(subcommon_genelist)
	  for (w in 1:(length_genelist)){
	    if (w==length_genelist){
	      subcommon_genelist$geneDistance[w] <- NA
	    }
	    else {
	      subcommon_genelist$geneDistance[w] <- subcommon_genelist$start[w+1] - subcommon_genelist$stop[w] 
	    }
	  }
	  print(subcommon_genelist)
	  
	  for (j in 1:nrow(subcommon_genelist)) {
	    print(j)
	    xs  <- c(subcommon_genelist$start[j]/1000000, subcommon_genelist$stop[j]/1000000+0.2, subcommon_genelist$stop[j]/1000000+0.2, subcommon_genelist$start[j]/1000000)
	    ys <- c(-0.95,-0.95,-0.8,-0.8)
	    
	    polygon(xs, ys, col="red",border=NA)
	    
	    if(is.na(subcommon_genelist$geneDistance[j])){
	      yy=-1
	      print("test")
	      previousSetDown=FALSE
	    }else if (subcommon_genelist$chr[1]==20 | subcommon_genelist$chr[1]==21 | subcommon_genelist$chr[1]==22){
	      yy=-1
	      previousSetDown=FALSE
	    }else if (subcommon_genelist$chr[1]==1 | subcommon_genelist$chr[1]==2 | subcommon_genelist$chr[1]==3 | subcommon_genelist$chr[1]==4 | subcommon_genelist$chr[1]==5 | subcommon_genelist$chr[1]==6 | subcommon_genelist$chr[1]==7){
	      if(subcommon_genelist$geneDistance[j] < 1500000){
	        if(previousSetDown==FALSE){
	          yy=-2.1
	          previousSetDown=TRUE
	        }else{
	          yy=-0.2
	          previousSetDown=FALSE
	        }
	      }else{
	        yy=-1
	        previousSetDown=FALSE
	      }
	    }else if(subcommon_genelist$geneDistance[j] < 500000){
	      if(previousSetDown==FALSE){
	        yy=-2.1
	        previousSetDown=TRUE
	      }else{
	        yy=-0.2
	        previousSetDown=FALSE
	      }
	    }else{
	      yy=-1
	      previousSetDown=FALSE
	    }
	    #print(subcommon_genelist$geneDistance[j])
	    #print(typeof(subcommon_genelist$geneDistance[j]))
	    #yy=-1
	    text(subcommon_genelist$start[j]/1000000+xx, yy, pos=2, labels = subcommon_genelist$geneName[j],col="red",cex=0.7,srt=90)
	  } 
	}
	

	graphics.off()

}
message('finished.')