setwd("/home/agoepfert/documents/pnet/38PNETs/Heatmap/20160527_matrix_v3/createPlots_server_v2")

source("https://bioconductor.org/biocLite.R")
#biocLite("rtracklayer")
#biocLite("Gviz")


## 
##https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf
##############################################################################
library("biomaRt")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

x <- listAttributes(ensembl)

geneLocations <- getBM(attributes = c("external_gene_name","exon_chrom_start","exon_chrom_end",
                     "chromosome_name","start_position","end_position","ensembl_exon_id",
                     "transcript_start","transcript_end",	"5_utr_start","5_utr_end","3_utr_start","3_utr_end","cds_start",
                     "cds_end","ensembl_exon_id"),
      filters = "external_gene_name", values = "BRCA2", mart=ensembl)



## IMPORTANT check if the right version of trackViewer is installed -> 1.8.3 is necessary
## If syou have problems with BiocInstaller version 3.2 is too old for R version 3.3
## TRy troubleshooting or download the newest package of bioconductor into your R folder and extract it there
## https://www.bioconductor.org/packages/release/bioc/html/BiocInstaller.html
library(Gviz)
library(rtracklayer)
library(biomaRt)
# https://www.bioconductor.org/packages/release/bioc/vignettes/trackViewer/inst/doc/trackViewer.html#lolliplot
##################################################################################
library(trackViewer)
SNP <- c(10, 100, 105, 108, 400, 410, 420, 600, 700, 805, 840, 1400, 1402)
x <- sample.int(100, length(SNP))
sample.gr <- GRanges("chr2", IRanges(SNP, width=1, names=paste0("snp", SNP)))
features <- GRanges("chr2", IRanges(c(1, 501, 1001), 
                                    width=c(120, 400, 405),
                                    names=paste0("block", 1:3)))
#features$fill <- c("#FF8833", "#51C6E6", "#DFA32D")

lolliplot(sample.gr, features)

#https://www.bioconductor.org/packages/release/bioc/vignettes/trackViewer/inst/doc/trackViewer.html#lolliplot 
#######################################################5.9 Variant Call Format (VCF) data
# library(VariantAnnotation)
# library(trackViewer)
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# library(org.Hs.eg.db)
# fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
# gr <- GRanges("22", IRanges(50968014, 50970514, names="TYMP"))
# tab <- TabixFile(fl)
# vcf <- readVcf(fl, "hg19", param=gr)
# mutation.frequency <- rowRanges(vcf)
# mcols(mutation.frequency) <- cbind(mcols(mutation.frequency), 
#                                    VariantAnnotation::info(vcf))
# mutation.frequency$border <- "gray30"
# mutation.frequency$color <- 
#   ifelse(grepl("^rs", names(mutation.frequency)), "lightcyan", "lavender")
# ## plot Global Allele Frequency based on AC/AN
# mutation.frequency$score <- mutation.frequency$AF*100
# seqlevelsStyle(gr) <- seqlevelsStyle(mutation.frequency) <- "UCSC" 
# trs <- geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg19.knownGene,
#                          org.Hs.eg.db,
#                          gr=gr)
# features <- c(range(trs[[1]]$dat), range(trs[[5]]$dat))
# names(features) <- c(trs[[1]]$name, trs[[5]]$name)
# features$fill <- c("lightblue", "mistyrose")
# features$height <- c(.02, .04)
# lolliplot(mutation.frequency, features, ranges=gr)








#https://bioconductor.riken.jp/packages/3.3/bioc/vignettes/GenVisR/inst/doc/GenVisR_intro.html
################################################################################
##Create input data
##biocLite("GenVisR")
# library(GenVisR)
# 
# data <- brcaMAF[brcaMAF$Hugo_Symbol == "TP53", c("Hugo_Symbol", "amino_acid_change_WU")]
# data <- as.data.frame(cbind(data, "ENST00000269305"))
# colnames(data) <- c("gene", "amino_acid_change", "transcript_name")
# 
# ## Call lolliplot
# lolliplot(data)
