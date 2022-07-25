library(AlphaSimR)
library(Rcpp)
library(ggplot2)
library(dplyr)
setwd("/home/rke27/Documents/")
set.seed(420)

##Reading in CO intervals from US NAM population
NAM <- read.table("NAM_US_COs_v4.txt", header = TRUE)
NAM <- NAM[order(NAM$Chr,NAM$Start),]
NAM$Start <- round(NAM$Start)

##Looking at gene density along chromosomes
#ref <- read.table("Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.1.gff3", header = TRUE)
#ref_genes <- ref[which(ref$chromosome == 'gene'),]
ref <- read.table("referencefile.csv", header = TRUE, sep =",")
ref_genes <- ref[which(ref$chromosome == 'gene'),]
ref_genes1 <- ref_genes[which(ref_genes$Chr1 == 'Chr1'),]
ref_genes2 <- ref_genes[which(ref_genes$Chr1 == 'Chr2'),]
ref_genes3 <- ref_genes[which(ref_genes$Chr1 == 'Chr3'),]
ref_genes4 <- ref_genes[which(ref_genes$Chr1 == 'Chr4'),]
ref_genes5 <- ref_genes[which(ref_genes$Chr1 == 'Chr5'),]
ref_genes6 <- ref_genes[which(ref_genes$Chr1 == 'Chr6'),]
ref_genes7 <- ref_genes[which(ref_genes$Chr1 == 'Chr7'),]
ref_genes8 <- ref_genes[which(ref_genes$Chr1 == 'Chr8'),]
ref_genes9 <- ref_genes[which(ref_genes$Chr1 == 'Chr9'),]
ref_genes10 <- ref_genes[which(ref_genes$Chr1 == 'Chr10'),]

chr1_CO <- NAM[ which(NAM$Chr == "chr1"),]
chr1 <- ggplot(data = chr1_CO) + geom_histogram(aes(x = Start), bins = 300)
chr1_p <- ggplot_build(chr1)
ref1 <- ggplot(data = ref_genes1) + geom_histogram(aes(x = X1), bins = 300)
ref1_p <- ggplot_build(ref1)
chr1_com <- cbind(chr1_p$data[[1]][,c(3,5,2)], ref1_p$data[[1]][,2])
colnames(chr1_com) <- c("Start", "End", "#COs", "#genes")
chr1_com$length <- (chr1_com$End-chr1_com$Start)/1000000
chr1_com$rate <- ((chr1_com$`#COs`/4713)*100)/chr1_com$length
chr1_com$density <- chr1_com$`#genes`/sum(chr1_com$`#genes`)

chr2_CO <- NAM[ which(NAM$Chr == "chr2"),]
chr2 <- ggplot(data = chr2_CO) + geom_histogram(aes(x = Start), bins = 244)
chr2_p <- ggplot_build(chr2)
ref2 <- ggplot(data = ref_genes2) + geom_histogram(aes(x = X1), bins = 244)
ref2_p <- ggplot_build(ref2)
chr2_com <- cbind(chr2_p$data[[1]][,c(3,5,2)], ref2_p$data[[1]][,2])
colnames(chr2_com) <- c("Start", "End", "#COs", "#genes")
chr2_com$length <- (chr2_com$End-chr2_com$Start)/1000000
chr2_com$rate <- ((chr2_com$`#COs`/4713)*100)/chr2_com$length
chr2_com$density <- chr2_com$`#genes`/sum(chr2_com$`#genes`)

chr3_CO <- NAM[ which(NAM$Chr == "chr3"),]
chr3 <- ggplot(data = chr3_CO) + geom_histogram(aes(x = Start), bins = 236)
chr3_p <- ggplot_build(chr3)
ref3 <- ggplot(data = ref_genes3) + geom_histogram(aes(x = X1), bins = 236)
ref3_p <- ggplot_build(ref3)
chr3_com <- cbind(chr3_p$data[[1]][,c(3,5,2)], ref3_p$data[[1]][,2])
colnames(chr3_com) <- c("Start", "End", "#COs", "#genes")
chr3_com$length <- (chr3_com$End-chr3_com$Start)/1000000
chr3_com$rate <- ((chr3_com$`#COs`/4713)*100)/chr3_com$length
chr3_com$density <- chr3_com$`#genes`/sum(chr3_com$`#genes`)

chr4_CO <- NAM[ which(NAM$Chr == "chr4"),]
chr4 <- ggplot(data = chr4_CO) + geom_histogram(aes(x = Start), bins = 246)
chr4_p <- ggplot_build(chr4)
ref4 <- ggplot(data = ref_genes4) + geom_histogram(aes(x = X1), bins = 246)
ref4_p <- ggplot_build(ref4)
chr4_com <- cbind(chr4_p$data[[1]][,c(3,5,2)], ref4_p$data[[1]][,2])
colnames(chr4_com) <- c("Start", "End", "#COs", "#genes")
chr4_com$length <- (chr4_com$End-chr4_com$Start)/1000000
chr4_com$rate <- ((chr4_com$`#COs`/4713)*100)/chr4_com$length
chr4_com$density <- chr4_com$`#genes`/sum(chr4_com$`#genes`)

chr5_CO <- NAM[ which(NAM$Chr == "chr5"),]
chr5 <- ggplot(data = chr5_CO) + geom_histogram(aes(x = Start), bins = 223)
chr5_p <- ggplot_build(chr5)
ref5 <- ggplot(data = ref_genes5) + geom_histogram(aes(x = X1), bins = 223)
ref5_p <- ggplot_build(ref5)
chr5_com <- cbind(chr5_p$data[[1]][,c(3,5,2)], ref5_p$data[[1]][,2])
colnames(chr5_com) <- c("Start", "End", "#COs", "#genes")
chr5_com$length <- (chr5_com$End-chr5_com$Start)/1000000
chr5_com$rate <- ((chr5_com$`#COs`/4713)*100)/chr5_com$length
chr5_com$density <- chr5_com$`#genes`/sum(chr5_com$`#genes`)

chr6_CO <- NAM[ which(NAM$Chr == "chr6"),]
chr6 <- ggplot(data = chr6_CO) + geom_histogram(aes(x = Start), bins = 173)
chr6_p <- ggplot_build(chr6)
ref6 <- ggplot(data = ref_genes6) + geom_histogram(aes(x = X1), bins = 173)
ref6_p <- ggplot_build(ref6)
chr6_com <- cbind(chr6_p$data[[1]][,c(3,5,2)], ref6_p$data[[1]][,2])
colnames(chr6_com) <- c("Start", "End", "#COs", "#genes")
chr6_com$length <- (chr6_com$End-chr6_com$Start)/1000000
chr6_com$rate <- ((chr6_com$`#COs`/4713)*100)/chr6_com$length
chr6_com$density <- chr6_com$`#genes`/sum(chr6_com$`#genes`)

chr7_CO <- NAM[ which(NAM$Chr == "chr7"),]
chr7 <- ggplot(data = chr7_CO) + geom_histogram(aes(x = Start), bins = 182)
chr7_p <- ggplot_build(chr7)
ref7 <- ggplot(data = ref_genes7) + geom_histogram(aes(x = X1), bins = 182)
ref7_p <- ggplot_build(ref7)
chr7_com <- cbind(chr7_p$data[[1]][,c(3,5,2)], ref7_p$data[[1]][,2])
colnames(chr7_com) <- c("Start", "End", "#COs", "#genes")
chr7_com$length <- (chr7_com$End-chr7_com$Start)/1000000
chr7_com$rate <- ((chr7_com$`#COs`/4713)*100)/chr7_com$length
chr7_com$density <- chr7_com$`#genes`/sum(chr7_com$`#genes`)

chr8_CO <- NAM[ which(NAM$Chr == "chr8"),]
chr8 <- ggplot(data = chr8_CO) + geom_histogram(aes(x = Start), bins = 181)
chr8_p <- ggplot_build(chr8)
ref8 <- ggplot(data = ref_genes8) + geom_histogram(aes(x = X1), bins = 181)
ref8_p <- ggplot_build(ref8)
chr8_com <- cbind(chr8_p$data[[1]][,c(3,5,2)], ref8_p$data[[1]][,2])
colnames(chr8_com) <- c("Start", "End", "#COs", "#genes")
chr8_com$length <- (chr8_com$End-chr8_com$Start)/1000000
chr8_com$rate <- ((chr8_com$`#COs`/4713)*100)/chr8_com$length
chr8_com$density <- chr8_com$`#genes`/sum(chr8_com$`#genes`)

chr9_CO <- NAM[ which(NAM$Chr == "chr9"),]
chr9 <- ggplot(data = chr9_CO) + geom_histogram(aes(x = Start), bins = 159)
chr9_p <- ggplot_build(chr9)
ref9 <- ggplot(data = ref_genes9) + geom_histogram(aes(x = X1), bins = 159)
ref9_p <- ggplot_build(ref9)
chr9_com <- cbind(chr9_p$data[[1]][,c(3,5,2)], ref9_p$data[[1]][,2])
colnames(chr9_com) <- c("Start", "End", "#COs", "#genes")
chr9_com$length <- (chr9_com$End-chr9_com$Start)/1000000
chr9_com$rate <- ((chr9_com$`#COs`/4713)*100)/chr9_com$length
chr9_com$density <- chr9_com$`#genes`/sum(chr9_com$`#genes`)

chr10_CO <- NAM[ which(NAM$Chr == "chr10"),]
chr10 <- ggplot(data = chr10_CO) + geom_histogram(aes(x = Start), bins = 150)
chr10_p <- ggplot_build(chr10)
ref10 <- ggplot(data = ref_genes10) + geom_histogram(aes(x = X1), bins = 150)
ref10_p <- ggplot_build(ref10)
chr10_com <- cbind(chr10_p$data[[1]][,c(3,5,2)], ref10_p$data[[1]][,2])
colnames(chr10_com) <- c("Start", "End", "#COs", "#genes")
chr10_com$length <- (chr10_com$End-chr10_com$Start)/1000000
chr10_com$rate <- ((chr10_com$`#COs`/4713)*100)/chr10_com$length
chr10_com$density <- chr10_com$`#genes`/sum(chr10_com$`#genes`)

##calculating recombination rate per bin of CO data
library(ggplot2)

#recombination frequency calc used:
# recomb. freq. = (# of COs/ size of population *100%)/ length of bin in Mb

##Spearmen correlation test to find correlation between gene density & recombination rate
#500kb correlation
cor.test(chr1_com$rate, chr1_com$density, method = "spearman", alternative = "greater")

cor.test(chr2_com$rate, chr2_com$density, method = "spearman", alternative = "greater")

cor.test(chr3_com$rate, chr3_com$density, method = "spearman", alternative = "greater")

cor.test(chr4_com$rate, chr4_com$density, method = "spearman", alternative = "greater")

cor.test(chr5_com$rate, chr5_com$density, method = "spearman", alternative = "greater")

cor.test(chr6_com$rate, chr6_com$density, method = "spearman", alternative = "greater")

cor.test(chr7_com$rate, chr7_com$density, method = "spearman", alternative = "greater")

cor.test(chr8_com$rate, chr8_com$density, method = "spearman", alternative = "greater")

cor.test(chr9_com$rate, chr9_com$density, method = "spearman", alternative = "greater")

cor.test(chr10_com$rate, chr10_com$density,  method = "spearman", alternative = "greater")

genomewide <- rbind(chr1_com, chr2_com, chr3_com, chr4_com, chr5_com, chr6_com, chr7_com, chr8_com, chr9_com, chr10_com)
cor.test(genomewide$rate, genomewide$density,  method = "spearman", alternative = "greater")
