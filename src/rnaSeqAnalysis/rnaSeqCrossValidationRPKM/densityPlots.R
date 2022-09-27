#This R script makes density plots of RPMK values for a given cell type
library(ggplot2)

#Get the CT variable from the command line
CT <- Sys.getenv("CT")
print(CT)

cType <- read.table(paste0(CT,"_TFs.tsv"), sep='\t', header = TRUE)
cType_RPKM <- cType$RPKM + 1
ggsave(filename=paste0(CT, "_log2_RPKM_plus1.png"),ggplot(data.frame(x = cType_RPKM), aes(x = x)) +
	geom_density(fill="lightblue") +
	geom_rug() +
	labs(x='RPKM') + 
	scale_x_continuous(trans='log2'), device='png')
