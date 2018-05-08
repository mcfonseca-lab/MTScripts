#############################################################################
###                                                                       ###
###   Antisense Project - 2018                                            ###
###   Produce PolyA plot (polyA+ and ployA-) and related Statistics.      ###
###   Author: Rui Luís (MCFonsecaLab)                                     ###
###                                                                       ###
#############################################################################

####
#Modules Used
library(ggplot2)
library(ggpubr)

####
#Data Input

PCpA_minus <- read.table(file = 'ProteinCoding_filter_quantification_tab_pA-_rep1.tsv', sep = '\t', header = FALSE)
PCpA_plus <- read.table(file = 'ProteinCoding_filter_quantification_tab_pA+_rep1.tsv', sep = '\t', header = FALSE)

AntipA_minus <- read.table(file = 'Antisense_filter_quantification_tab_pA-_rep1.tsv', sep = '\t', header = FALSE)
AntipA_plus <- read.table(file = 'Antisense_filter_quantification_tab_pA+_rep1.tsv', sep = '\t', header = FALSE)

Hist_minus <- read.table(file = 'histones_Hugo_Ensembl_id_pA-_rep1.txt', sep = '\t', header =FALSE)
Hist_plus <- read.table(file = 'histones_Hugo_Ensembl_id_pA+_rep1.txt', sep = '\t', header =FALSE)

PC_TOTAL_NOT_NAT_minus <- read.table(file = 'ProteinCoding_TOTAL_NOT_NAT_pA-_rep1.tsv', sep = '\t', header =FALSE)
PC_TOTAL_NOT_NAT_plus <- read.table(file = 'ProteinCoding_TOTAL_NOT_NAT_pA+_rep1.tsv', sep = '\t', header =FALSE)
####
#Build R DataFrame for polyA datsets with TPM (Transcripts per Kibobase Million)
PCTPMpA_minus <- PCpA_minus$V9
PCTPMpA_plus <- PCpA_plus$V9
AntiTPMpA_minus <- AntipA_minus$V9
AntiTPMpA_plus <- AntipA_plus$V9
HistTPMpA_minus <- Hist_minus$V9
HistTPMpA_plus <- Hist_plus$V9
PC_TOTAL_TPM_NOT_NAT_minus <- PC_TOTAL_NOT_NAT_minus$V9
PC_TOTAL_TPM_NOT_NAT_plus <- PC_TOTAL_NOT_NAT_plus$V9

Type <- c(rep("PolyA+",65),rep("PolyA-",65),rep("PolyA+",65),rep("PolyA-",65),rep("PolyA+",85),rep("PolyA-",85),rep("PolyA+",2776),rep("PolyA-",2776))
Descrip <- c(rep("Antisense",65),rep("Antisense",65),rep("ProteinCoding",65),rep("ProteinCoding",65),rep("Histones",85),rep("Histones",85),
             rep("ProteinCoding_TOTAL_NOT_NAT",2776),rep("ProteinCoding_TOTAL_NOT_NAT",2776))
TPMvalue <- c(log2(AntiTPMpA_plus),log2(AntiTPMpA_minus),log2(PCTPMpA_plus),log2(PCTPMpA_minus),
              log2(HistTPMpA_plus),log2(HistTPMpA_minus),log2(PC_TOTAL_TPM_NOT_NAT_plus),log2(PC_TOTAL_TPM_NOT_NAT_minus))

PolyATPM <- data.frame(Type,Descrip,TPMvalue)
PolyATPM$Descrip <- factor(PolyATPM$Descrip, c('Antisense', 'ProteinCoding', 'Histones','ProteinCoding_TOTAL_NOT_NAT'))


####
#Plot in ggplot2 boxplot as pair grouped side by side polyA+ and polyA- (y Axix in log2(TPMs))
base_plot <- ggplot(PolyATPM, aes(x=Descrip, y=TPMvalue,fill=Type))
base_plot  + geom_boxplot()   + scale_fill_manual(values=c("#FF0000", "#0000FF")) +
  stat_compare_means(method = "wilcox.test") + labs(x="Genes Type", y = "Log2( TPM )")

####
#Statistics with wilcox test ; also known as Mann-Whitney (two-samples, non-parametric)
wilcox.test(AntiTPMpA_minus,AntiTPMpA_plus)
wilcox.test(PCTPMpA_minus,PCTPMpA_plus)
wilcox.test(HistTPMpA_minus,HistTPMpA_plus)
x<- wilcox.test(PC_TOTAL_TPM_NOT_NAT_minus,PC_TOTAL_TPM_NOT_NAT_plus)

####
#Change data to present as polyA+ divided by polyA- 
Descrip2 <- c(rep("Antisense",65),rep("ProteinCoding",65),rep("Histones",85),rep("ProteinCoding_TOTAL_NOT_NAT",2776))
TPMvalue2 <- c(log2(AntiTPMpA_plus/AntiTPMpA_minus),log2(PCTPMpA_plus/PCTPMpA_minus), log2(HistTPMpA_plus/HistTPMpA_minus),
               log2(PC_TOTAL_TPM_NOT_NAT_plus/PC_TOTAL_TPM_NOT_NAT_minus))
PolyATPM2 <- data.frame(Descrip2,TPMvalue2)

####
#Plot in ggplot2 boxplot as polyA+ divided by polyA-  (y Axix in log2(TPMs))
PolyATPM2$Descrip2 <- factor(PolyATPM2$Descrip2, c('Antisense', 'ProteinCoding', 'Histones','ProteinCoding_TOTAL_NOT_NAT'))

my_comparisons <- list( c("Antisense", "ProteinCoding"), c("Antisense", "Histones"), c("ProteinCoding", "Histones"),c("Antisense","ProteinCoding_TOTAL_NOT_NAT"))
base_plot <- ggplot(PolyATPM2, aes(x=Descrip2, y=TPMvalue2, fill=Descrip2)) 
base_plot  + geom_boxplot() + scale_fill_manual(values=c("#FF0000", "#0000FF", 	"#808080","#808080"))  + stat_compare_means(comparisons = my_comparisons) +
  labs(x="Genes Type", y = "Log2( TPMs pA+ / pA- Nucleoplasm)")




##################################################################################
#Based on :
#Schlackow M, Nojima T, Gomes T, Dhir A, Carmo-Fonseca M, Proudfoot NJ. 
#Distinctive Patterns of Transcription and RNA Processing for Human lincRNAs. 
#Molecular Cell. 2017;65(1):25-38. doi:10.1016/j.molcel.2016.11.029.
#
#Using only fragment numbers instead of using log2(TPM) values.

####
#Data Input
PCpA_minus <- read.table(file = 'ProteinCoding_feature_counts_pA-_rep1.txt', sep = ' ', header = FALSE)
PCpA_plus <- read.table(file = 'ProteinCoding_feature_counts_pA+_rep1.txt', sep = ' ', header = FALSE)

AntipA_minus <- read.table(file = 'Antisense_feature_counts_pA-_rep1.txt', sep = ' ', header = FALSE)
AntipA_plus <- read.table(file = 'Antisense_feature_counts_pA+_rep1.txt', sep = ' ', header = FALSE)


#Build R DataFrame for polyA datsets with Number of Fragments (previously calculated with featureCounts)
FragCounts_PC_minus <- PCpA_minus$V2
FragCounts_PC_plus <- PCpA_plus$V2
FragCounts_Anti_minus <- AntipA_minus$V2
FragCounts_Anti_plus <- AntipA_plus$V2

Type <- c(rep("PolyA+",65),rep("PolyA-",65),rep("PolyA+",65),rep("PolyA-",65))
Descrip <- c(rep("Antisense",65),rep("Antisense",65),rep("ProteinCoding",65),rep("ProteinCoding",65))
TPMvalue <- c(FragCounts_Anti_plus,FragCounts_Anti_minus,FragCounts_PC_plus,FragCounts_PC_minus)
PolyATPM <- data.frame(Type,Descrip,TPMvalue)

base_plot <- ggplot(PolyATPM, aes(x=Descrip, y=TPMvalue,fill=Type))
base_plot  + geom_boxplot()   + scale_fill_manual(values=c("#FF0000", "#0000FF")) +
  stat_compare_means(method = "wilcox.test") + labs(x="Genes Type", y = "# Fragments")


#############################################################################################

Descrip2 <- c(rep("Antisense",65),rep("ProteinCoding",65))

TPMvalue2 <- c(log10(FragCounts_Anti_plus/FragCounts_Anti_minus),log10(FragCounts_PC_plus/FragCounts_PC_minus))

PolyATPM2 <- data.frame(Descrip2,TPMvalue2)
base_plot <- ggplot(PolyATPM2, aes(x=Descrip2, y=TPMvalue2, fill=Descrip2)) 
base_plot  + geom_boxplot() + scale_fill_manual(values=c("#FF0000", "#0000FF"))  + stat_compare_means(method = "wilcox.test")  +
  labs(x="Genes Type", y = "Log2( pA+ / pA- Nucleoplasm)")




