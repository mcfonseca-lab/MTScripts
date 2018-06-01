#############################################################################
###                                                                       ###
###   Antisense Project - 2018                                            ###
###   Produce Correlation Plots between Replicates                        ###
###   Author: Rui Luís (MCFonsecaLab)                                     ###
###                                                                       ###
#############################################################################

#Packages
source("https://bioconductor.org/biocLite.R")
#biocLite("scater")
library(scater)
library(ggpubr)

countdata <- read.table("total_counts_Nucleoplasm.txt", header=TRUE)
# Remove first five columns (chr, start, end, strand, length)

countdata <- countdata[ ,c(1,7:ncol(countdata))]
# Remove .bam or .sam from filenames

colnames(countdata) <- gsub("\\.[sb]am$", "", colnames(countdata))


PC_genes65 <- read.table("PC.txt", header=F)
AS_genes65 <- read.table("AS.txt", header=F)

countdata <- as.data.frame(countdata)
resdataPC <- countdata[countdata[,1] %in% PC_genes65[,1],]
resdataAS <- countdata[countdata[,1] %in% AS_genes65[,1],]

resdataPC <- cbind(resdataPC,rep("Protein Coding",nrow(resdataPC)))
resdataAS <- cbind(resdataAS,rep("Natural Antisense",nrow(resdataAS)))

colnames(resdataPC) <- c("GeneID","Cod1Replicate1","Cod1Replicate2","Cod2Replicate1","Cod2Replicate2","GeneType")
colnames(resdataAS) <- c("GeneID","Cod1Replicate1","Cod1Replicate2","Cod2Replicate1","Cod2Replicate2","GeneType")

DTTotal <- rbind(resdataPC, resdataAS)

ggscatter(DTTotal, x = "Cod1Replicate1", y = "Cod1Replicate2",
          add = "reg.line",  # Add regressin linel
          xlim = c(0,12500),
          #ylim = c(0,10000),
          xlab = "Replicate 1 (Reads Number)",
          ylab = "Replicate 2 (Reads Number)",
          color = "GeneType",
          palette = c( "red", "blue"),
          repel = TRUE,
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.sep = "\n")
) + ggtitle("siEx3 Nucleoplasm Fraction Replicates Correlation") +
  theme(plot.title = element_text(hjust = 0.5, face="bold",size=16))



ggscatter(DTTotal, x = "Cod2Replicate1", y = "Cod2Replicate2",
          add = "reg.line",  # Add regressin linel
          xlim = c(0,13000),
          #ylim = c(0,10000),
          xlab = "Replicate 1 (Reads Number)",
          ylab = "Replicate 2 (Reads Number)",
          color = "GeneType",
          repel = TRUE,
          palette = c( "red", "blue"),
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.sep = "\n")
) + ggtitle("siLuc Nucleoplasm Fraction Replicates Correlation") +
  theme(plot.title = element_text(hjust = 0.5, face="bold",size=16))


