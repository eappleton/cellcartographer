#This script plots results for BC counting analysis

#Load libraries
library(dict)
library(ggplot2)
library(ggpubr)
library(svglite)

### Create plots ###
plotFreq <- function(filename) {
	
	plotMat <- read.table(filename, sep="," , header=TRUE)
	nameBase <- tools::file_path_sans_ext(filename)

	#Make a normalized frequency version of the counts
	total = sum(plotMat[,2])
	plotMat["BC_COUNTS_N"] <- plotMat["BC_COUNTS"]/total
	
	#Histogram of Number of appearances for top 200 genes
	x <- ggbarplot(plotMat, x = "Unnamed..0", y = "BC_COUNTS", fill = "lightgray", xlab = "Barcoded TF", ylab = "Number of Appearances", sort.val = "desc", x.text.angle = 45, lab.size = 2)
	y <- ggbarplot(plotMat, x = "Unnamed..0", y = "BC_COUNTS_N", fill = "lightgray", xlab = "Barcoded TF", ylab = "Relative Appearance", sort.val = "desc", x.text.angle = 45, lab.size = 2)
	if (grepl("ALL", nameBase)) {
		ggsave(paste0(nameBase,".png"), x, width=30, height=5, limitsize = FALSE)
		ggsave(paste0(nameBase,"_ND.png"), y, width=30, height=5, limitsize = FALSE)
	} else {
		ggsave(paste0(nameBase,".png"), x, width=8, height=5, limitsize = FALSE)
		ggsave(paste0(nameBase,".svg"), x, width=8, height=5, limitsize = FALSE)
		ggsave(paste0(nameBase,"_ND.png"), y, width=8, height=5, limitsize = FALSE)
		ggsave(paste0(nameBase,"_ND.svg"), y, width=8, height=5, limitsize = FALSE)
	}	
}

### DESC ###
setwd("barcode/results")
files <- list.files()
barcodeFiles <- files[grepl(".csv", files, ignore.case=TRUE)]
# print(barcodeFiles)
for (f in barcodeFiles) {
	print(f)
	plotFreq(f)
}
