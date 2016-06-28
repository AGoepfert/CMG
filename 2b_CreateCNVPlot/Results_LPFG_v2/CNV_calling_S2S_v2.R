# GET arguments
args <- commandArgs(TRUE)
file <- args[1]
pos_file <- paste(file,'.pos',sep='')
logR_file <- paste(file,'.all.logR',sep='')
out_file <- args[2]
library(DNAcopy)
library("fields")
library("zoo")
options(scipen=999)
options(bitmapType="cairo")
data <- read.csv(file,header=TRUE,sep=' ')
colnames(data) <- c('CHR','BIN.START','BIN.END','LogR')

blacklist <- '/home/agoepfert/documents/01_general_datasets/common_cnv/CNV_common_NIPT.bed'
genelist <- '/home/agoepfert/documents/mesothelioma_marieke/new_approach/LPWG2/2b_CreateCNVPlot/Results_LPFG_v2/Cancer_census_genes_v2.txt'

# read the CNV blacklist
common_cnv <- read.table(blacklist,header=FALSE)
colnames(common_cnv) <- c('chr','start','stop')
common_cnv$chr <- sub('X',23,common_cnv$chr)
common_cnv$chr <- sub('Y',24,common_cnv$chr)

# read the CNV blacklist
interesting_genelist <- read.table(genelist,header=FALSE)
colnames(interesting_genelist) <- c('chr','start','stop','geneName')
interesting_genelist$chr <- sub('X',23,interesting_genelist$chr)
interesting_genelist$chr <- sub('Y',24,interesting_genelist$chr)

# merge and remove zero
data$CHR <- sub('chr','',data$CHR)


data$CHR <- sub('X',23,data$CHR)
data$CHR <- sub('Y',24,data$CHR)
# normalize by median.
data$LogR <- data$LogR - median(data$LogR)

#segment
CNA.object <- CNA(data$LogR,data$CHR,data$BIN.START,sampleid='Avg.(Tumor/Norma)',data.type='logratio')
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
cnvs <- calls[(calls$seg.mean >= 0.15 | calls$seg.mean <= -0.15),]
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

all.LogR <- read.csv(logR_file,header=TRUE,sep=' ')
all.LogR[all.LogR < -4] <- -3
all.LogR[all.LogR > 4] <- 4
pos_data <- read.csv(pos_file,header=TRUE,sep=' ')
pos_data[,1] <- sub('chr','',pos_data[,1])
sample1 <- 'Tumor'
sample2 <- 'Normal'

pal.1=colorRampPalette(c("red","white","green"),space="rgb")
breaks <- seq(min(all.LogR,na.rm = T),max(all.LogR,na.rm=T),length.out=100)

#cnvs <- data.frame(chrom=c())
# plot
for (i in 1:22) {
 	png(file=paste("Avg.",sample1,".vs.",sample2,".page_",i,".png",sep=""), bg="white", width=1060, height=800)
	par(mfrow=c(2,1))
	# first chr
	#chr = 2*i -1
	chr = i
	message(paste('plotting chromosome: ',chr,sep=''))
	subcnv <- cnvs[cnvs$chrom == chr,]
	subcommon <- common_cnv[common_cnv$chr == chr,]
	subcommon_genelist <- interesting_genelist[interesting_genelist$chr == chr,]
	
	subdata <- data[data$CHR == chr,]
	subdata <- subdata[order(subdata$BIN.START),]
	rows <-nrow(subdata)
	if (rows > 0) {
		plot(subdata$BIN.START/1000000,subdata$LogR, main=paste('Avg.(',sample1,"/",sample2, ") : LogR for Chr", chr, sep=" "),xlab="Position (Mb)", ylab="Log2(R)", pch=".",ylim=c(-3,2),xaxs="i", yaxs="i")
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
			#lines(subdata$BIN.START/1000000,smooth(subdata$LogR),col='blue')
		}
	} else {
		message('no data')
		plot(c(0,1),c(0,0), main=paste(sample1,"/",sample2, " : LogR for Chr", chr, sep=" "),xlab="Position (Mb)", ylab="Log2(R)", pch="",ylim=c(-3,2),xaxs="i", yaxs="i")
		grid(nx=NA,ny=NULL,col='red',lty='solid',lwd=1)	
	}
	
	# sample cnvs
	if (nrow(subcnv)>0) {
		for (j in 1:nrow(subcnv)) {
			# the box hightlight
			xs <- c(subcnv$loc.start[j]/1000000, subcnv$loc.end[j]/1000000, subcnv$loc.end[j]/1000000, subcnv$loc.start[j]/1000000)
			ys <- c(-3,-3,2,2)
			if (subcnv$seg.mean[j] < -0.15) {
				polygon(xs, ys, col=rgb(160,32,240,100,maxColorValue=255),border=NA)
				abline(v=xs[1],lwd=1,col=rgb(160,32,240,50,maxColorValue=255))
				abline(v=xs[2],lwd=1,col=rgb(160,32,240,50,maxColorValue=255))

			}
			if (subcnv$seg.mean[j] > 0.15) {
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
			ys <- c(-2.95,-2.95,-2.8,-2.8)
			
			polygon(xs, ys, col=rgb(75,75,75,250,maxColorValue=255),border=NA)
		} 

	}
	
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
	  print(subcommon_genelist)
	  for (j in 1:nrow(subcommon_genelist)) {
	    print(j)
	    xs  <- c(subcommon_genelist$start[j]/1000000, subcommon_genelist$stop[j]/1000000+0.2, subcommon_genelist$stop[j]/1000000+0.2, subcommon_genelist$start[j]/1000000)
	    ys <- c(-0.95,-0.95,-0.8,-0.8)
	    
	    polygon(xs, ys, col="red",border=NA)
	    
	 
	    #cat(subcommon_genelist$start[j],subcommon_genelist$stop[j],subcommon_genelist$start[j]-subcommon_genelist$stop[j],"test",subcommon_genelist$geneName[j],"\n")
	    text(subcommon_genelist$start[j]/1000000+xx, -1, pos=2, labels = subcommon_genelist$geneName[j],col="red",cex=0.75,srt=90)
	  } 
  }

	## add legend
	if (i == 1) {
		labs <- paste("log(S1/S2).SD: ",round(logR.sd,2),"undo.segment by ",undo_type,", threshold: ",undo_thresh)
		x <- max(subdata$BIN.START)/1000000
		y <- 1.2
		text(x,y,labels = labs,pos=2,cex=0.9)
	}
	needed.rows <- which(pos_data[,1] == i)
	needed.data <- all.LogR[needed.rows,]

	image.plot(pos_data[needed.rows,2],seq(1,ncol(all.LogR)),xlab='',as.matrix(needed.data),col=pal.1(length(breaks)-1),legend.mar=4,breaks=breaks,horizontal=T)
	graphics.off()
}
message('finished.')
