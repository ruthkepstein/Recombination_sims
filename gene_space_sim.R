library(AlphaSimR)
library(Rcpp)
library(ggplot2)
library(dplyr)
setwd("C:/Users/16192/Documents/PNAS_Simulations")
set.seed(420)

##reading in SNPs from B73xMo17 based on v4 B73 ref
final_snps <- read.table("SNP_V4.bed", header = FALSE)
colnames(final_snps) <- c("Chr#", "SNP Start", "SNP End")
#2000 SNPs genome-wide
final_snps <- sample_n(final_snps, 2000)
final_snps <- final_snps[order(final_snps$`Chr#`,final_snps$`SNP Start`),]
write.csv(final_snps, "C:/Users/16192/Documents/PNAS_Simulations/final_snps.csv", row.names = FALSE)

chr1_snp <- final_snps[ which(final_snps$`Chr#` == "chr1"),]
hist(chr1_snp$`SNP Start`, breaks = 100)
chr1_snp$rate <- NA
#making SNPs start at 0
chr1_snp$`SNP End` <- chr1_snp$`SNP End` - min(chr1_snp$`SNP Start`)
chr1_snp$`SNP Start` <- chr1_snp$`SNP Start`- min(chr1_snp$`SNP Start`)

chr2_snp <- final_snps[ which(final_snps$`Chr#` == "chr2"),]
hist(chr2_snp$`SNP Start`, breaks = 100)
chr2_snp$rate <- NA
chr2_snp$`SNP End` <- chr2_snp$`SNP End` - min(chr2_snp$`SNP Start`)
chr2_snp$`SNP Start` <- chr2_snp$`SNP Start`- min(chr2_snp$`SNP Start`)

chr3_snp <- final_snps[ which(final_snps$`Chr#` == "chr3"),]
hist(chr3_snp$`SNP Start`, breaks = 100)
chr3_snp$rate <- NA
chr3_snp$`SNP End` <- chr3_snp$`SNP End` - min(chr3_snp$`SNP Start`)
chr3_snp$`SNP Start` <- chr3_snp$`SNP Start`- min(chr3_snp$`SNP Start`)

chr4_snp <- final_snps[ which(final_snps$`Chr#` == "chr4"),]
hist(chr4_snp$`SNP Start`, breaks = 100)
chr4_snp$rate <- NA
chr4_snp$`SNP End` <- chr4_snp$`SNP End` - min(chr4_snp$`SNP Start`)
chr4_snp$`SNP Start` <- chr4_snp$`SNP Start`- min(chr4_snp$`SNP Start`)

chr5_snp <- final_snps[ which(final_snps$`Chr#` == "chr5"),]
hist(chr5_snp$`SNP Start`, breaks = 100)
chr5_snp$rate <- NA
chr5_snp$`SNP End` <- chr5_snp$`SNP End` - min(chr5_snp$`SNP Start`)
chr5_snp$`SNP Start` <- chr5_snp$`SNP Start`- min(chr5_snp$`SNP Start`)

chr6_snp <- final_snps[ which(final_snps$`Chr#` == "chr6"),]
hist(chr6_snp$`SNP Start`, breaks = 100)
chr6_snp$rate <- NA
chr6_snp$`SNP End` <- chr6_snp$`SNP End` - min(chr6_snp$`SNP Start`)
chr6_snp$`SNP Start` <- chr6_snp$`SNP Start`- min(chr6_snp$`SNP Start`)

chr7_snp <- final_snps[ which(final_snps$`Chr#` == "chr7"),]
hist(chr7_snp$`SNP Start`, breaks = 100)
chr7_snp$rate <- NA
chr7_snp$`SNP End` <- chr7_snp$`SNP End` - min(chr7_snp$`SNP Start`)
chr7_snp$`SNP Start` <- chr7_snp$`SNP Start`- min(chr7_snp$`SNP Start`)

chr8_snp <- final_snps[ which(final_snps$`Chr#` == "chr8"),]
hist(chr8_snp$`SNP Start`, breaks = 100)
chr8_snp$rate <- NA
chr8_snp$`SNP End` <- chr8_snp$`SNP End` - min(chr8_snp$`SNP Start`)
chr8_snp$`SNP Start` <- chr8_snp$`SNP Start`- min(chr8_snp$`SNP Start`)

chr9_snp <- final_snps[ which(final_snps$`Chr#` == "chr9"),]
hist(chr9_snp$`SNP Start`, breaks = 100)
chr9_snp$rate <- NA
chr9_snp$`SNP End` <- chr9_snp$`SNP End` - min(chr9_snp$`SNP Start`)
chr9_snp$`SNP Start` <- chr9_snp$`SNP Start`- min(chr9_snp$`SNP Start`)

chr10_snp <- final_snps[ which(final_snps$`Chr#` == "chr10"),]
hist(chr10_snp$`SNP Start`, breaks = 100)
chr10_snp$rate <- NA
chr10_snp$`SNP End` <- chr10_snp$`SNP End` - min(chr10_snp$`SNP Start`)
chr10_snp$`SNP Start` <- chr10_snp$`SNP Start`- min(chr10_snp$`SNP Start`)

##Reading in CO intervals from US NAM population
NAM <- read.table("NAM_US_COs.txt", header = FALSE)
NAM <- NAM[,-c(2:4)]
colnames(NAM) <- c("Chr", "CO Start", "CO End")
NAM <- NAM[order(NAM$Chr,NAM$`CO Start`),]
#pop_size <- read.table("pop_size_nam.txt", header = TRUE)
#US_pop_size <- pop_size[1:4713,]

chr1_CO <- NAM[ which(NAM$Chr == 1),]
chr1_CO$midpoint <- (chr1_CO$`CO Start`+ chr1_CO$`CO End`)/2
hist(chr1_CO$`CO Start`, breaks = 300)

chr2_CO <- NAM[ which(NAM$Chr == 2),]
hist(chr2_CO$`CO Start`, breaks = 300)
chr2_CO$midpoint <- (chr2_CO$`CO Start`+ chr2_CO$`CO End`)/2

chr3_CO <- NAM[ which(NAM$Chr == 3),]
hist(chr3_CO$`CO Start`, breaks = 300)
chr3_CO$midpoint <- (chr3_CO$`CO Start`+ chr3_CO$`CO End`)/2

chr4_CO <- NAM[ which(NAM$Chr == 4),]
hist(chr4_CO$`CO Start`, breaks = 300)
chr4_CO$midpoint <- (chr4_CO$`CO Start`+ chr4_CO$`CO End`)/2

chr5_CO <- NAM[ which(NAM$Chr == 5),]
hist(chr5_CO$`CO Start`, breaks = 300)
chr5_CO$midpoint <- (chr5_CO$`CO Start`+ chr5_CO$`CO End`)/2

chr6_CO <- NAM[ which(NAM$Chr == 6),]
hist(chr6_CO$`CO Start`, breaks = 300)
chr6_CO$midpoint <- (chr6_CO$`CO Start`+ chr6_CO$`CO End`)/2

chr7_CO <- NAM[ which(NAM$Chr == 7),]
hist(chr7_CO$`CO Start`, breaks = 300)
chr7_CO$midpoint <- (chr7_CO$`CO Start`+ chr7_CO$`CO End`)/2

chr8_CO <- NAM[ which(NAM$Chr == 8),]
hist(chr8_CO$`CO Start`, breaks = 300)
chr8_CO$midpoint <- (chr8_CO$`CO Start`+ chr8_CO$`CO End`)/2

chr9_CO <- NAM[ which(NAM$Chr == 9),]
hist(chr9_CO$`CO Start`, breaks = 300)
chr9_CO$midpoint <- (chr9_CO$`CO Start`+ chr9_CO$`CO End`)/2

chr10_CO <- NAM[ which(NAM$Chr == 10),]
hist(chr10_CO$`CO Start`, breaks = 300)
chr10_CO$midpoint <- (chr10_CO$`CO Start`+ chr10_CO$`CO End`)/2

###using CO rate to infer genetic map distances

##calculating recombination rate per bin of CO data
library(dlookr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(OneR)

#recombination frequency calc used:
# recomb. freq. = (# of COs/ size of population *100%)/ length of bin in Mb

#bin crossovers into 200 uneven bins
chr1_CO <- chr1_CO[order(chr1_CO$`CO Start`),]
chr1_bin <- binning(chr1_CO$midpoint, nbins = 300, type = "kmeans")
chr1_bin <- as.data.frame(summary(chr1_bin))
#transforming data; making bin interval into 2 columns
chr1_bin <- within(chr1_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr1_bin$levels), ',', fixed=TRUE))))
chr1_bin <- do.call(data.frame, chr1_bin)
chr1_bin <- chr1_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr1_bin <- chr1_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr1_bin[1,4] <- 502954
#making intervals start at 0
chr1_bin$foo.X1 <- chr1_bin$foo.X1 - 502954
chr1_bin$foo.X2 <- chr1_bin$foo.X2 - 502954
#expanding last bin to include last SNP site to avoid NAs in future
chr1_bin[300,5] <- max(chr1_snp$`SNP End`)
#adding length of bin as column and making in Mb
chr1_bin$length <- (chr1_bin$foo.X2-chr1_bin$foo.X1)/1000000
chr1_bin$rate <- ((chr1_bin$freq/4713)*100)/chr1_bin$length

chr2_CO <- chr2_CO[order(chr2_CO$`CO Start`),]
chr2_bin <- binning(chr2_CO$midpoint, nbins = 220, type = "kmeans")
chr2_bin <- as.data.frame(summary(chr2_bin))
chr2_bin <- within(chr2_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr2_bin$levels), ',', fixed=TRUE))))
chr2_bin <- do.call(data.frame, chr2_bin)
chr2_bin <- chr2_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr2_bin <- chr2_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr2_bin[1,4] <- 440104
chr2_bin$foo.X1 <- chr2_bin$foo.X1 - 440104
chr2_bin$foo.X2 <- chr2_bin$foo.X2 - 440104
chr2_bin[220,5] <- max(chr2_snp$`SNP End`)
chr2_bin$length <- (chr2_bin$foo.X2-chr2_bin$foo.X1)/1000000
chr2_bin$rate <- ((chr2_bin$freq/4713)*100)/chr2_bin$length

#have not converted rest of chromosomes to what I did in 1 & 2
chr3_CO <- chr3_CO[order(chr3_CO$`CO Start`),]
chr3_bin <- binning(chr3_CO$midpoint, nbins = 200, type = "kmeans")
chr3_bin <- as.data.frame(summary(chr3_bin))
chr3_bin <- within(chr3_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr3_bin$levels), ',', fixed=TRUE))))
chr3_bin <- do.call(data.frame, chr3_bin)
chr3_bin <- chr3_bin %>% mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr3_bin <- chr3_bin %>% mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr3_bin[1,4] <- 865390
chr3_bin$foo.X1 <- chr3_bin$foo.X1 - 865390
chr3_bin$foo.X2 <- chr3_bin$foo.X2 - 865390
chr3_bin[200,5] <- max(chr3_snp$`SNP End`)
chr3_bin$length <- (chr3_bin$foo.X2-chr3_bin$foo.X1)/1000000
chr3_bin$rate <- ((chr3_bin$freq/4713)*100)/chr3_bin$length

chr4_CO <- chr4_CO[order(chr4_CO$`CO Start`),]
chr4_bin <- binning(chr4_CO$midpoint, nbins = 200, type = "kmeans")
chr4_bin <- as.data.frame(summary(chr4_bin))
chr4_bin <- within(chr4_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr4_bin$levels), ',', fixed=TRUE))))
chr4_bin <- do.call(data.frame, chr4_bin)
chr4_bin <- chr4_bin %>% mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr4_bin <- chr4_bin %>% mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr4_bin[1,4] <- 272401
chr4_bin$foo.X1 <- chr4_bin$foo.X1 - 272401
chr4_bin$foo.X2 <- chr4_bin$foo.X2 - 272401
chr4_bin[200,5] <- max(chr4_snp$`SNP End`)
chr4_bin$length <- (chr4_bin$foo.X2-chr4_bin$foo.X1)/1000000
chr4_bin$rate <- ((chr4_bin$freq/4713)*100)/chr4_bin$length

chr5_CO <- chr5_CO[order(chr5_CO$`CO Start`),]
chr5_bin <- binning(chr5_CO$midpoint, nbins = 120, type = "kmeans")
chr5_bin <- as.data.frame(summary(chr5_bin))
chr5_bin <- within(chr5_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr5_bin$levels), ',', fixed=TRUE))))
chr5_bin <- do.call(data.frame, chr5_bin)
chr5_bin <- chr5_bin %>% mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr5_bin <- chr5_bin %>% mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr5_bin[1,4] <- 267335.5
chr5_bin$foo.X1 <- chr5_bin$foo.X1 - 267335.5
chr5_bin$foo.X2 <- chr5_bin$foo.X2 - 267335.5
chr5_bin[120,5] <- max(chr5_snp$`SNP End`)
chr5_bin$length <- (chr5_bin$foo.X2-chr5_bin$foo.X1)/1000000
chr5_bin$rate <- ((chr5_bin$freq/4713)*100)/chr5_bin$length

chr6_CO <- chr6_CO[order(chr6_CO$`CO Start`),]
chr6_bin <- binning(chr6_CO$midpoint, nbins = 80, type = "kmeans")
chr6_bin <- as.data.frame(summary(chr6_bin))
chr6_bin <- within(chr6_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr6_bin$levels), ',', fixed=TRUE))))
chr6_bin <- do.call(data.frame, chr6_bin)
chr6_bin <- chr6_bin %>% mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr6_bin <- chr6_bin %>% mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr6_bin[1,4] <- 197266.5
chr6_bin$foo.X1 <- chr6_bin$foo.X1 - 197266.5
chr6_bin$foo.X2 <- chr6_bin$foo.X2 - 197266.5
chr6_bin[80,5] <- max(chr6_snp$`SNP End`)
chr6_bin$length <- (chr6_bin$foo.X2-chr6_bin$foo.X1)/1000000
chr6_bin$rate <- ((chr6_bin$freq/4713)*100)/chr6_bin$length

chr7_CO <- chr7_CO[order(chr7_CO$`CO Start`),]
chr7_bin <- binning(chr7_CO$midpoint, nbins = 200, type = "kmeans")
chr7_bin <- as.data.frame(summary(chr7_bin))
chr7_bin <- within(chr7_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr7_bin$levels), ',', fixed=TRUE))))
chr7_bin <- do.call(data.frame, chr7_bin)
chr7_bin <- chr7_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr7_bin <- chr7_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr7_bin[1,4] <- 375904
chr7_bin$foo.X1 <- chr7_bin$foo.X1 - 375904
chr7_bin$foo.X2 <- chr7_bin$foo.X2 - 375904
chr7_bin[200,5] <- max(chr7_snp$`SNP End`)
chr7_bin$length <- (chr7_bin$foo.X2-chr7_bin$foo.X1)/1000000
chr7_bin$rate <- ((chr7_bin$freq/4713)*100)/chr7_bin$length

chr8_CO <- chr8_CO[order(chr8_CO$`CO Start`),]
chr8_bin <- binning(chr8_CO$midpoint, nbins = 150, type = "kmeans")
chr8_bin <- as.data.frame(summary(chr8_bin))
chr8_bin <- within(chr8_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr8_bin$levels), ',', fixed=TRUE))))
chr8_bin <- do.call(data.frame, chr8_bin)
chr8_bin <- chr8_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr8_bin <- chr8_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr8_bin[1,4] <- 132181
chr8_bin$foo.X1 <- chr8_bin$foo.X1 - 132181
chr8_bin$foo.X2 <- chr8_bin$foo.X2 - 132181
chr8_bin[150,5] <- max(chr8_snp$`SNP End`)
chr8_bin$length <- (chr8_bin$foo.X2-chr8_bin$foo.X1)/1000000
chr8_bin$rate <- ((chr8_bin$freq/4713)*100)/chr8_bin$length

chr9_CO <- chr9_CO[order(chr9_CO$`CO Start`),]
chr9_bin <- binning(chr9_CO$midpoint, nbins = 120, type = "kmeans")
chr9_bin <- as.data.frame(summary(chr9_bin))
chr9_bin <- within(chr9_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr9_bin$levels), ',', fixed=TRUE))))
chr9_bin <- do.call(data.frame, chr9_bin)
chr9_bin <- chr9_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr9_bin <- chr9_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr9_bin[1,4] <- 317217.5
chr9_bin$foo.X1 <- chr9_bin$foo.X1 - 317217.5
chr9_bin$foo.X2 <- chr9_bin$foo.X2 - 317217.5
chr9_bin[120,5] <- max(chr9_snp$`SNP End`)
chr9_bin$length <- (chr9_bin$foo.X2-chr9_bin$foo.X1)/1000000
chr9_bin$rate <- ((chr9_bin$freq/4713)*100)/chr9_bin$length

chr10_CO <- chr10_CO[order(chr10_CO$`CO Start`),]
chr10_bin <- binning(chr10_CO$midpoint, nbins = 150, type = "kmeans")
chr10_bin <- as.data.frame(summary(chr10_bin))
chr10_bin <- within(chr10_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr10_bin$levels), ',', fixed=TRUE))))
chr10_bin <- do.call(data.frame, chr10_bin)
chr10_bin <- chr10_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr10_bin <- chr10_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr10_bin[1,4] <- 698530
chr10_bin$foo.X1 <- chr10_bin$foo.X1 - 698530
chr10_bin$foo.X2 <- chr10_bin$foo.X2 - 698530
chr10_bin[150,5] <- max(chr10_snp$`SNP End`)
chr10_bin$length <- (chr10_bin$foo.X2-chr10_bin$foo.X1)/1000000
chr10_bin$rate <- ((chr10_bin$freq/4713)*100)/chr10_bin$length
#use plot to look at distribution of k-means

##assigning frequency to SNPs based on recombination frequency in each bin
snp_rate <- function(chr_bin, chr_snp){
  for(i in 1:nrow(chr_snp)){
    for(k in 1:nrow(chr_bin)){
      if(isTRUE((chr_snp$`SNP Start`[i] >= chr_bin$foo.X1[k]) && (chr_snp$`SNP Start`[i] <= chr_bin$foo.X2[k]))){
        chr_snp$rate[i] <- chr_bin$rate[k]
      }
    }
  }
  print(chr_snp)
}
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

genes_bin1 <- as.data.frame(summary(binning(ref_genes1$X1, nbins = 300, type = "kmeans")))
genes_bin1 <- within(genes_bin1, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin1$levels), ',', fixed=TRUE))))
genes_bin1 <- do.call(data.frame, genes_bin1)
genes_bin1 <- genes_bin1 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin1 <- genes_bin1 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin1[1,4]<-44289
genes_bin1$foo.X1 <- genes_bin1$foo.X1 - 44289
genes_bin1$foo.X2 <- genes_bin1$foo.X2 - 44289
genes_bin1[300,5] <- max(chr1_snp$`SNP End`)

genes_bin2 <- as.data.frame(summary(binning(ref_genes2$X1, nbins = 300, type = "kmeans")))
genes_bin2 <- within(genes_bin2, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin2$levels), ',', fixed=TRUE))))
genes_bin2 <- do.call(data.frame, genes_bin2)
genes_bin2 <- genes_bin2 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin2 <- genes_bin2 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin2[1,4]<-40178
genes_bin2$foo.X1 <- genes_bin2$foo.X1 - 40178
genes_bin2$foo.X2 <- genes_bin2$foo.X2 - 40178
genes_bin2[300,5] <- max(chr2_snp$`SNP End`)

genes_bin3 <- as.data.frame(summary(binning(ref_genes3$X1, nbins = 300, type = "kmeans")))
genes_bin3 <- within(genes_bin3, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin3$levels), ',', fixed=TRUE))))
genes_bin3 <- do.call(data.frame, genes_bin3)
genes_bin3 <- genes_bin3 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin3 <- genes_bin3 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin3[1,4]<-191993
genes_bin3$foo.X1 <- genes_bin3$foo.X1 - 191993
genes_bin3$foo.X2 <- genes_bin3$foo.X2 - 191993
genes_bin3[300,5] <- max(chr3_snp$`SNP End`)

genes_bin4 <- as.data.frame(summary(binning(ref_genes4$X1, nbins = 300, type = "kmeans")))
genes_bin4 <- within(genes_bin4, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin4$levels), ',', fixed=TRUE))))
genes_bin4 <- do.call(data.frame, genes_bin4)
genes_bin4 <- genes_bin4 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin4 <- genes_bin4 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin4[1,4]<-50320
genes_bin4$foo.X1 <- genes_bin4$foo.X1 - 50320
genes_bin4$foo.X2 <- genes_bin4$foo.X2 - 50320
genes_bin4[300,5] <- max(chr4_snp$`SNP End`)

genes_bin5 <- as.data.frame(summary(binning(ref_genes5$X1, nbins = 300, type = "kmeans")))
genes_bin5 <- within(genes_bin5, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin5$levels), ',', fixed=TRUE))))
genes_bin5 <- do.call(data.frame, genes_bin5)
genes_bin5 <- genes_bin5 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin5 <- genes_bin5 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin5[1,4]<-32368
genes_bin5$foo.X1 <- genes_bin5$foo.X1 - 32368
genes_bin5$foo.X2 <- genes_bin5$foo.X2 - 32368
genes_bin5[300,5] <- max(chr5_snp$`SNP End`)

genes_bin6 <- as.data.frame(summary(binning(ref_genes6$X1, nbins = 300, type = "kmeans")))
genes_bin6 <- within(genes_bin6, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin6$levels), ',', fixed=TRUE))))
genes_bin6 <- do.call(data.frame, genes_bin6)
genes_bin6 <- genes_bin6 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin6 <- genes_bin6 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin6[1,4]<-168702
genes_bin6$foo.X1 <- genes_bin6$foo.X1 - 168702
genes_bin6$foo.X2 <- genes_bin6$foo.X2 - 168702
genes_bin6[300,5] <- max(chr6_snp$`SNP End`)

genes_bin7 <- as.data.frame(summary(binning(ref_genes7$X1, nbins = 300, type = "kmeans")))
genes_bin7 <- within(genes_bin7, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin7$levels), ',', fixed=TRUE))))
genes_bin7 <- do.call(data.frame, genes_bin7)
genes_bin7 <- genes_bin7 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin7 <- genes_bin7 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin7[1,4]<-107948
genes_bin7$foo.X1 <- genes_bin7$foo.X1 - 107948
genes_bin7$foo.X2 <- genes_bin7$foo.X2 - 107948
genes_bin7[300,5] <- max(chr7_snp$`SNP End`)

genes_bin8 <- as.data.frame(summary(binning(ref_genes8$X1, nbins = 300, type = "kmeans")))
genes_bin8 <- within(genes_bin8, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin8$levels), ',', fixed=TRUE))))
genes_bin8 <- do.call(data.frame, genes_bin8)
genes_bin8 <- genes_bin8 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin8 <- genes_bin8 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin8[1,4]<-77622
genes_bin8$foo.X1 <- genes_bin8$foo.X1 - 77622
genes_bin8$foo.X2 <- genes_bin8$foo.X2 - 77622
genes_bin8[300,5] <- max(chr8_snp$`SNP End`)

genes_bin9 <- as.data.frame(summary(binning(ref_genes9$X1, nbins = 300, type = "kmeans")))
genes_bin9 <- within(genes_bin9, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin9$levels), ',', fixed=TRUE))))
genes_bin9 <- do.call(data.frame, genes_bin9)
genes_bin9 <- genes_bin9 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin9 <- genes_bin9 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin9[1,4]<-59372
genes_bin9$foo.X1 <- genes_bin9$foo.X1 - 59372
genes_bin9$foo.X2 <- genes_bin9$foo.X2 - 59372
genes_bin9[300,5] <- max(chr9_snp$`SNP End`)

genes_bin10 <- binning(ref_genes10$X1, nbins = 300, type = "kmeans")
plot(genes_bin10)
genes_bin10 <- as.data.frame(summary(binning(ref_genes10$X1, nbins = 300, type = "kmeans")))
genes_bin10 <- within(genes_bin10, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin10$levels), ',', fixed=TRUE))))
genes_bin10 <- do.call(data.frame, genes_bin10)
genes_bin10 <- genes_bin10 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin10 <- genes_bin10 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin10[1,4]<-130112
genes_bin10$foo.X1 <- genes_bin10$foo.X1 - 130112
genes_bin10$foo.X2 <- genes_bin10$foo.X2 - 130112
genes_bin10[300,5] <- max(chr10_snp$`SNP End`)


#REASSIGN FREQ loop through snp positions, all snp that fall in bins with low density get 1/4 of mean rate
snp_uniform_rate <- function(genes_bin, chr_snp, chr_snps_mean){
  for(i in 1:nrow(chr_snp)){
    for(k in 1:nrow(genes_bin)){
      if(isTRUE((chr_snp$`SNP Start`[i] >= genes_bin$foo.X1[k]) && (chr_snp$`SNP Start`[i] <= genes_bin$foo.X2[k]))){
        if(isTRUE((genes_bin$Quartile[k]==1))){
          chr_snp$rate[i] <- 0.5*chr_snps_mean
        } 
        else{
          chr_snp$rate[i] <- chr_snps_mean
        }
      }
    }
  }
  return(chr_snp)
}

#using function, converted SNP start to Mb to get cM/Mb for final genetic position
chr1_snp2 <- snp_rate(chr1_bin, chr1_snp)
#reassigning recombination rates, calling uniform function
chr1_snps_mean <- mean(chr1_snp2$rate)
genes_bin1$Quartile<-cut(genes_bin1$freq,quantile(genes_bin1$freq),include.lowest=TRUE,labels=FALSE)
chr1_snp2 <-snp_uniform_rate(genes_bin1, chr1_snp, chr1_snps_mean)
#chr1_snp2$rate[1:100] <- 1.097985
#chr1_snp2$rate[101:199] <- 0.73199
#chr1_snp2$rate[200:300] <- 1.097985
chr1_snp2$`SNP Start`<- chr1_snp2$`SNP Start`/1000000
chr1_snp2 <- chr1_snp2[order(chr1_snp2$`SNP Start`),]
#smoothing the recombination rate so transitions between bins are not so abrupt
chr1_spl <- smooth.spline(chr1_snp2$rate, spar = 0.8)
#creation of genetic positions from smoothed recombination rate
chr1_snp2$pos <- (chr1_snp2$`SNP Start`*chr1_spl$y)
#graph to look at Mb vs. cM along chromosome
plot(chr1_snp2$`SNP Start`, chr1_snp2$pos)
ggplot(chr1_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
#graph to look at Mb vs. cM/Mb to see recombination rate along chromosome
plot(chr1_snp2$`SNP Start`, chr1_snp2$pos/chr1_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 1 Recombination Distribution")
chr1_finalpos <- chr1_snp2[order(chr1_snp2$pos),]
#want False to input into AlphaSimR
is.unsorted(chr1_finalpos$pos)
#plot again to make sure it looks the same
plot(chr1_snp2$`SNP Start`, chr1_finalpos$pos/chr1_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 1 Recombination Distribution")
plot(chr1_finalpos$`SNP Start`, chr1_finalpos$pos)

chr2_snp2 <- snp_rate(chr2_bin, chr2_snp)
chr2_snps_mean <- mean(chr2_snp2$rate)
genes_bin2$Quartile<-cut(genes_bin2$freq,quantile(genes_bin2$freq),include.lowest=TRUE,labels=FALSE)
chr2_snp2 <-snp_uniform_rate(genes_bin2, chr2_snp, chr2_snps_mean)
#chr2_snp2$rate[1:79] <- 1.318753
#chr2_snp2$rate[80:158] <- 0.8791686
#chr2_snp2$rate[159:237] <- 1.318753
chr2_snp2$`SNP Start` <- chr2_snp2$`SNP Start`/1000000
chr2_snp2 <- chr2_snp2[-(196:237),]
chr2_spl <- smooth.spline(chr2_snp2$rate, spar = 0.7)
chr2_snp2$pos <- (chr2_snp2$`SNP Start`*chr2_spl$y)
plot(chr2_snp2$`SNP Start`, chr2_snp2$pos)
ggplot(chr2_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
plot(chr2_snp2$`SNP Start`, chr2_snp2$pos/chr2_snp2$`SNP Start`, type = "l")
chr2_finalpos <- chr2_snp2[order(chr2_snp2$pos),]
is.unsorted(chr2_finalpos$pos)
plot(chr2_snp2$`SNP Start`, chr2_finalpos$pos/chr2_snp2$`SNP Start`, type = "l")

chr3_snp2 <- snp_rate(chr3_bin, chr3_snp)
#chr3_snp2$rate[1:73] <- 0.9275442
#chr3_snp2$rate[74:147] <- 0.6183628
#chr3_snp2$rate[148:219] <- 0.9275442
chr3_snps_mean <- mean(chr3_snp2$rate)
genes_bin3$Quartile<-cut(genes_bin3$freq,quantile(genes_bin3$freq),include.lowest=TRUE,labels=FALSE)
chr3_snp2 <-snp_uniform_rate(genes_bin3, chr3_snp, chr3_snps_mean)
chr3_snp2$`SNP Start` <- chr3_snp2$`SNP Start`/1000000
chr3_snp2 <- chr3_snp2[-(205:219),]
chr3_spl <- smooth.spline(chr3_snp2$rate, spar = 0.8)
chr3_snp2$pos <- (chr3_snp2$`SNP Start`*chr3_spl$y)
plot(chr3_snp2$`SNP Start`, chr3_snp2$pos)
ggplot(chr3_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
plot(chr3_snp2$`SNP Start`, chr3_snp2$pos/chr3_snp2$`SNP Start`, type = "l")
chr3_finalpos <- chr3_snp2[order(chr3_snp2$pos),]
is.unsorted(chr3_finalpos$pos)
plot(chr3_snp2$`SNP Start`, chr3_finalpos$pos/chr3_snp2$`SNP Start`, type = "l")

chr4_snp2 <- snp_rate(chr4_bin, chr4_snp)
#chr4_snp2$rate[1:85] <- 0.7232293
#chr4_snp2$rate[86:170] <- 0.4821528
#chr4_snp2$rate[171:256] <- 0.7232293
chr4_snps_mean <- mean(chr4_snp2$rate)
genes_bin4$Quartile<-cut(genes_bin4$freq,quantile(genes_bin4$freq),include.lowest=TRUE,labels=FALSE)
chr4_snp2 <-snp_uniform_rate(genes_bin4, chr4_snp, chr4_snps_mean)
chr4_snp2$`SNP Start` <- chr4_snp2$`SNP Start`/1000000
chr4_spl <- smooth.spline(chr4_snp2$rate, spar = 0.8)
chr4_snp2$pos <- (chr4_snp2$`SNP Start`*chr4_spl$y)
plot(chr4_snp2$`SNP Start`, chr4_snp2$pos)
ggplot(chr4_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
plot(chr4_snp2$`SNP Start`, chr4_snp2$pos/chr4_snp2$`SNP Start`, type = "l")
chr4_finalpos <- chr4_snp2[order(chr4_snp2$pos),]
is.unsorted(chr4_finalpos$pos)
plot(chr4_snp2$`SNP Start`, chr4_finalpos$pos/chr4_snp2$`SNP Start`, type = "l")

chr5_snp2 <- snp_rate(chr5_bin, chr5_snp)
chr5_snps_mean <- mean(chr5_snp2$rate)
genes_bin5$Quartile<-cut(genes_bin5$freq,quantile(genes_bin5$freq),include.lowest=TRUE,labels=FALSE)
chr5_snp2 <-snp_uniform_rate(genes_bin5, chr5_snp, chr5_snps_mean)
chr5_snp2$`SNP Start` <- chr5_snp2$`SNP Start`/1000000
chr5_spl <- smooth.spline(chr5_snp2$rate, spar = 1)
chr5_snp2$pos <- (chr5_snp2$`SNP Start`*chr5_spl$y)
plot(chr5_snp2$`SNP Start`, chr5_snp2$pos)
ggplot(chr5_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
plot(chr5_snp2$`SNP Start`, chr5_snp2$pos/chr5_snp2$`SNP Start`, type = "l")
chr5_finalpos <- chr5_snp2[order(chr5_snp2$pos),]
is.unsorted(chr5_finalpos$pos)
plot(chr5_snp2$`SNP Start`, chr5_finalpos$pos/chr5_snp2$`SNP Start`, type = "l")

#chr 6 is lowkey fuked up
chr6_snp2 <- snp_rate(chr6_bin, chr6_snp)
chr6_snps_mean <- mean(chr6_snp2$rate)
genes_bin6$Quartile<-cut(genes_bin6$freq,quantile(genes_bin6$freq),include.lowest=TRUE,labels=FALSE)
chr6_snp2 <-snp_uniform_rate(genes_bin6, chr6_snp, chr6_snps_mean)
chr6_snp2$`SNP Start` <- chr6_snp2$`SNP Start`/1000000
chr6_spl <- smooth.spline(chr6_snp2$rate, spar = 1)
chr6_snp2$pos <- (chr6_snp2$`SNP Start`*chr6_spl$y)
plot(chr6_snp2$`SNP Start`, chr6_snp2$pos)
ggplot(chr6_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
plot(chr6_snp2$`SNP Start`, chr6_snp2$pos/chr6_snp2$`SNP Start`, type = "l")
chr6_finalpos <- chr6_snp2[order(chr6_snp2$pos),]
is.unsorted(chr6_finalpos$pos)
plot(chr6_snp2$`SNP Start`, chr6_finalpos$pos/chr6_snp2$`SNP Start`, type = "l")

chr7_snp2 <- snp_rate(chr7_bin, chr7_snp)
chr7_snps_mean <- mean(chr7_snp2$rate)
genes_bin7$Quartile<-cut(genes_bin7$freq,quantile(genes_bin7$freq),include.lowest=TRUE,labels=FALSE)
chr7_snp2 <-snp_uniform_rate(genes_bin7, chr7_snp, chr7_snps_mean)
chr7_snp2$`SNP Start` <- chr7_snp2$`SNP Start`/1000000
chr7_spl <- smooth.spline(chr7_snp2$rate, spar = .9)
chr7_snp2$pos <- (chr7_snp2$`SNP Start`*chr7_spl$y)
plot(chr7_snp2$`SNP Start`, chr7_snp2$pos)
ggplot(chr7_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
plot(chr7_snp2$`SNP Start`, chr7_snp2$pos/chr7_snp2$`SNP Start`, type = "l")
chr7_finalpos <- chr7_snp2[order(chr7_snp2$pos),]
is.unsorted(chr7_finalpos$pos)
plot(chr7_snp2$`SNP Start`, chr7_finalpos$pos/chr7_snp2$`SNP Start`, type = "l")

chr8_snp2 <- snp_rate(chr8_bin, chr8_snp)
chr8_snps_mean <- mean(chr8_snp2$rate)
genes_bin8$Quartile<-cut(genes_bin8$freq,quantile(genes_bin8$freq),include.lowest=TRUE,labels=FALSE)
chr8_snp2 <-snp_uniform_rate(genes_bin8, chr8_snp, chr8_snps_mean)
chr8_snp2$`SNP Start` <- chr8_snp2$`SNP Start`/1000000
chr8_spl <- smooth.spline(chr8_snp2$rate, spar = .9)
chr8_snp2$pos <- (chr8_snp2$`SNP Start`*chr8_spl$y)
plot(chr8_snp2$`SNP Start`, chr8_snp2$pos)
ggplot(chr8_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
plot(chr8_snp2$`SNP Start`, chr8_snp2$pos/chr8_snp2$`SNP Start`, type = "l")
chr8_finalpos <- chr8_snp2[order(chr8_snp2$pos),]
is.unsorted(chr8_finalpos$pos)
plot(chr8_snp2$`SNP Start`, chr8_finalpos$pos/chr8_snp2$`SNP Start`, type = "l")

chr9_snp2 <- snp_rate(chr9_bin, chr9_snp)
chr9_snps_mean <- mean(chr9_snp2$rate)
genes_bin9$Quartile<-cut(genes_bin9$freq,quantile(genes_bin9$freq),include.lowest=TRUE,labels=FALSE)
chr9_snp2 <-snp_uniform_rate(genes_bin9, chr9_snp, chr9_snps_mean)
chr9_snp2$`SNP Start` <- chr9_snp2$`SNP Start`/1000000
chr9_spl <- smooth.spline(chr9_snp2$rate, spar = 1.1)
chr9_snp2$pos <- (chr9_snp2$`SNP Start`*chr9_spl$y)
plot(chr9_snp2$`SNP Start`, chr9_snp2$pos)
ggplot(chr9_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
plot(chr9_snp2$`SNP Start`, chr9_snp2$pos/chr9_snp2$`SNP Start`, type = "l")
chr9_finalpos <- chr9_snp2[order(chr9_snp2$pos),]
is.unsorted(chr9_finalpos$pos)
plot(chr9_snp2$`SNP Start`, chr9_finalpos$pos/chr9_snp2$`SNP Start`, type = "l")

chr10_snp2 <- snp_rate(chr10_bin, chr10_snp)
chr10_snps_mean <- mean(chr10_snp2$rate)
genes_bin10$Quartile<-cut(genes_bin10$freq,quantile(genes_bin10$freq),include.lowest=TRUE,labels=FALSE)
chr10_snp2 <-snp_uniform_rate(genes_bin10, chr10_snp, chr10_snps_mean)
chr10_snp2$`SNP Start` <- chr10_snp2$`SNP Start`/1000000
chr10_spl <- smooth.spline(chr10_snp2$rate, spar = 1.15)
chr10_snp2$pos <- (chr10_snp2$`SNP Start`*chr10_spl$y)
plot(chr10_snp2$`SNP Start`, chr10_snp2$pos)
ggplot(chr10_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
plot(chr10_snp2$`SNP Start`, chr10_snp2$pos/chr10_snp2$`SNP Start`, type = "l")
chr10_finalpos <- chr10_snp2[order(chr10_snp2$pos),]
is.unsorted(chr10_finalpos$pos)
plot(chr10_snp2$`SNP Start`, chr10_finalpos$pos/chr10_snp2$`SNP Start`, type = "l")

#Putting the final genetic map together
chr1 <- chr1_finalpos$pos/100
chr1len <- length(chr1)
dim(chr1) <- c(chr1len,1)
chr1 <- list(chr1)

chr2 <- chr2_finalpos$pos/100
chr2len <- length(chr2)
dim(chr2) <- c(chr2len,1)
chr2 <- list(chr2)

chr3 <- chr3_finalpos$pos/100
chr3len <- length(chr3)
dim(chr3) <- c(chr3len,1)
chr3 <- list(chr3)

chr4 <- chr4_finalpos$pos/100
chr4len <- length(chr4)
dim(chr4) <- c(chr4len,1)
chr4 <- list(chr4)

chr10 <- chr10_finalpos$pos/100
chr10len <- length(chr10)
dim(chr10) <- c(chr10len,1)
chr10 <- list(chr10)

chr5 <- chr5_finalpos$pos/100
chr5len <- length(chr5)
dim(chr5) <- c(chr5len,1)
chr5 <- list(chr5)

chr6 <- chr6_finalpos$pos/100
chr6len <- length(chr6)
dim(chr6) <- c(chr6len,1)
chr6 <- list(chr6)

chr7 <- chr7_finalpos$pos/100
chr7len <- length(chr7)
dim(chr7) <- c(chr7len,1)
chr7 <- list(chr7)

chr8 <- chr8_finalpos$pos/100
chr8len <- length(chr8)
dim(chr8) <- c(chr8len,1)
chr8 <- list(chr8)

chr9 <- chr9_finalpos$pos/100
chr9len <- length(chr9)
dim(chr9) <- c(chr9len,1)
chr9 <- list(chr9)

final_map <- list(chr1[[1]], chr2[[1]], 
                  chr3[[1]], chr4[[1]], chr5[[1]], 
                  chr6[[1]], chr7[[1]], chr8[[1]], 
                  chr9[[1]], chr10[[1]])

#Creating vector of centromere positions
real_centromere <- c(69.19425, 29.31842, 43.34131, 42.1754, 47.11507,
                     43.40838, 30.01601, 43.69602, 47.12911, 21.17685)
real_centromere <- real_centromere/100

#Creating haplotypes with assumption of LD present
is_in_LD <- function(chr_finalpos, LD_snps){
  for(i in 1:nrow(chr_finalpos)){
    if(chr_finalpos$rate[i] <= 1.5){
      LD_snps[i] <- chr_finalpos$`SNP Start`[i]
    }
  }
  return(LD_snps)
}
chr1_LD_snps <- c()
chr1_LD_snps <- is_in_LD(chr1_finalpos, chr1_LD_snps)
chr1_LD_snps <- as.data.frame(chr1_LD_snps)

chr2_LD_snps <- c()
chr2_LD_snps <- is_in_LD(chr2_finalpos, chr2_LD_snps)
chr2_LD_snps <- as.data.frame(chr2_LD_snps)

chr3_LD_snps <- c()
chr3_LD_snps <- is_in_LD(chr3_finalpos, chr3_LD_snps)
chr3_LD_snps <- as.data.frame(chr3_LD_snps)

chr4_LD_snps <- c()
chr4_LD_snps <- is_in_LD(chr4_finalpos, chr4_LD_snps)
chr4_LD_snps <- as.data.frame(chr4_LD_snps)

chr5_LD_snps <- c()
chr5_LD_snps <- is_in_LD(chr5_finalpos, chr5_LD_snps)
chr5_LD_snps <- as.data.frame(chr5_LD_snps)

chr6_LD_snps <- c()
chr6_LD_snps <- is_in_LD(chr6_finalpos, chr6_LD_snps)
chr6_LD_snps <- as.data.frame(chr6_LD_snps)

chr7_LD_snps <- c()
chr7_LD_snps <- is_in_LD(chr7_finalpos, chr7_LD_snps)
chr7_LD_snps <- as.data.frame(chr7_LD_snps)

chr8_LD_snps <- c()
chr8_LD_snps <- is_in_LD(chr8_finalpos, chr8_LD_snps)
chr8_LD_snps <- as.data.frame(chr8_LD_snps)

chr9_LD_snps <- c()
chr9_LD_snps <- is_in_LD(chr9_finalpos, chr9_LD_snps)
chr9_LD_snps <- as.data.frame(chr9_LD_snps)

chr10_LD_snps <- c()
chr10_LD_snps <- is_in_LD(chr10_finalpos, chr10_LD_snps)
chr10_LD_snps <- as.data.frame(chr10_LD_snps)

#change row.names to change # of individuals we want
chr1_haplo <- matrix(data = NA, nrow = 200, ncol = nrow(chr1_finalpos))
row.names(chr1_haplo) <- 1:200
colnames(chr1_haplo) <- chr1_finalpos$`SNP Start`

chr2_haplo <- matrix(data = NA, nrow = 200, ncol = nrow(chr2_finalpos))
row.names(chr2_haplo) <- 1:200
colnames(chr2_haplo) <- chr2_finalpos$`SNP Start`

chr3_haplo <- matrix(data = NA, nrow = 200, ncol = nrow(chr3_finalpos))
row.names(chr3_haplo) <- 1:200
colnames(chr3_haplo) <- chr3_finalpos$`SNP Start`

chr4_haplo <- matrix(data = NA, nrow = 200, ncol = nrow(chr4_finalpos))
row.names(chr4_haplo) <- 1:200
colnames(chr4_haplo) <- chr4_finalpos$`SNP Start`

chr5_haplo <- matrix(data = NA, nrow = 200, ncol = nrow(chr5_finalpos))
row.names(chr5_haplo) <- 1:200
colnames(chr5_haplo) <- chr5_finalpos$`SNP Start`

chr6_haplo <- matrix(data = NA, nrow = 200, ncol = nrow(chr6_finalpos))
row.names(chr6_haplo) <- 1:200
colnames(chr6_haplo) <- chr6_finalpos$`SNP Start`

chr7_haplo <- matrix(data = NA, nrow = 200, ncol = nrow(chr7_finalpos))
row.names(chr7_haplo) <- 1:200
colnames(chr7_haplo) <- chr7_finalpos$`SNP Start`

chr8_haplo <- matrix(data = NA, nrow = 200, ncol = nrow(chr8_finalpos))
row.names(chr8_haplo) <- 1:200
colnames(chr8_haplo) <- chr8_finalpos$`SNP Start`

chr9_haplo <- matrix(data = NA, nrow = 200, ncol = nrow(chr9_finalpos))
row.names(chr9_haplo) <- 1:200
colnames(chr9_haplo) <- chr9_finalpos$`SNP Start`

chr10_haplo <- matrix(data = NA, nrow = 200, ncol = nrow(chr10_finalpos))
row.names(chr10_haplo) <- 1:200
colnames(chr10_haplo) <- chr10_finalpos$`SNP Start`

fill_matrix <- function(chr_LD_snps, chr_haplo){
  for(i in 1:nrow(chr_LD_snps)){
    if(is.na(chr_LD_snps[i,])){
      for(k in 1:nrow(chr_haplo)){
        chr_haplo[k,i] = sample(0:1,1)
      }
    }
    else{
      x <- sample(0:1,1)
      for(k in 1:nrow(chr_haplo)){
        chr_haplo[k,i] = x
      }
    }
  }
  return(chr_haplo)
}
chr1_haplo <- fill_matrix(chr1_LD_snps, chr1_haplo)

chr2_haplo <- fill_matrix(chr2_LD_snps, chr2_haplo)

chr3_haplo <- fill_matrix(chr3_LD_snps, chr3_haplo)

chr4_haplo <- fill_matrix(chr4_LD_snps, chr4_haplo)

chr5_haplo <- fill_matrix(chr5_LD_snps, chr5_haplo)
chr5_haplo[is.na(chr5_haplo)] <- sample(0:1,1)

chr6_haplo <- fill_matrix(chr6_LD_snps, chr6_haplo)
chr6_haplo[is.na(chr6_haplo)] <- sample(0:1,1)

chr7_haplo <- fill_matrix(chr7_LD_snps, chr7_haplo)

chr8_haplo <- fill_matrix(chr8_LD_snps, chr8_haplo)

chr9_haplo <- fill_matrix(chr9_LD_snps, chr9_haplo)

chr10_haplo <- fill_matrix(chr10_LD_snps, chr10_haplo)

final_haplo <- list(chr1_haplo, chr2_haplo, chr3_haplo, chr4_haplo, chr5_haplo,
                    chr6_haplo, chr7_haplo, chr8_haplo, chr9_haplo, chr10_haplo)

#write.table(final_haplo, "C:/Users/sajai/example_haplotypes.csv")
#inal_haplo <- read.table("example_haplotypes.csv")

###Simulating a realistic breeding program in maize

##Setting up program with altered map
#create founder population with gen map & haplotypes
Founder_pop <- newMapPop(genMap = final_map, haplotypes = final_haplo, inbred = TRUE, ploidy = 2L)

#updating with actual centromere positions
Founder_pop@centromere <- real_centromere

#change depending on if polygenic or oligenic trait
nQtlPerChr <- 1
#creating simulation parameters with the founder population
#tracking recombination as well
SP = SimParam$new(Founder_pop)$setTrackRec(TRUE)
#assuming crossover interference
SP$v = 2.6
#number of COs coming from non-interfering pathway
SP$p = 0.2
#adding an additive trait
SP$addTraitA(
  nQtlPerChr,
  mean = 0,
  var = 1,
  corA = NULL,
  gamma = FALSE,
  shape = 1,
  force = FALSE
)
#setting variability for the trait
SP$setVarE(h2=0.5)

##initiation of first population from founders
pop1 <- newPop(Founder_pop, simParam = SP)
#find phenotypes for the first population in order to do phenotypic selection
pop1 <- setPheno(
  pop1,
  h2 = NULL,
  H2 = NULL,
  onlyPheno = FALSE,
  simParam = SP
)
#selecting top individuals from first population
#select top 10 using phenotype for one trait
pop1_sel <- selectInd(pop1, nInd = 10, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)

#put it together to iterate one program 40 times with 20 generations of selection
#after 20 gen of selection, inbred or DH
pop1_pheno <- matrix(nrow=40,ncol=10)
pop2_pheno <- matrix(nrow=40,ncol=10)
pop3_pheno <- matrix(nrow=40,ncol=10)
pop4_pheno <- matrix(nrow=40,ncol=10)
pop5_pheno <- matrix(nrow=40,ncol=10)
pop6_pheno <- matrix(nrow=40,ncol=10)
pop7_pheno <- matrix(nrow=40,ncol=10)
pop8_pheno <- matrix(nrow=40,ncol=10)
pop9_pheno <- matrix(nrow=40,ncol=10)
pop10_pheno <- matrix(nrow=40,ncol=10)
pop11_pheno <- matrix(nrow=40,ncol=10)
pop12_pheno <- matrix(nrow=40,ncol=10)
pop13_pheno <- matrix(nrow=40,ncol=10)
pop14_pheno <- matrix(nrow=40,ncol=10)
pop15_pheno <- matrix(nrow=40,ncol=10)
pop16_pheno <- matrix(nrow=40,ncol=10)
pop17_pheno <- matrix(nrow=40,ncol=10)
pop18_pheno <- matrix(nrow=40,ncol=10)
pop19_pheno <- matrix(nrow=40,ncol=10)
pop20_pheno <- matrix(nrow=40,ncol=10)

  for(i in 1:40){
    founderPop <- newMapPop(genMap = final_map, haplotypes = final_haplo, inbred = TRUE, ploidy = 2L)
    founderPop@centromere <- real_centromere
    SP = SimParam$new(founderPop)
    SP$setTrackRec(TRUE)
    SP$v = 2.6
    SP$p = 0.2
    nQtlPerChr = 1
    SP$addTraitA(nQtlPerChr)
    SP$setVarE(h2=0.9)
    
    pop <- newPop(founderPop, simParam = SP)
    pop <- setPheno(pop = pop, h2 = 0.9, simParam = SP)
    
    pop_F1 <- randCross(pop, nCrosses = 50, nProgeny = 10, simParam = SP)
    pop_DH <- makeDH(pop_F1, nDH = 1, simParam = SP)
    pop_DH <- setPheno(pop_DH, h2 = 0.9, simParam = SP)
    
    pop1_sel <- selectInd(pop_DH, nInd = 50, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
    pop1_sel_cross <- randCross(pop1_sel, 10, nProgeny = 50, simParam = SP)
    pop1_sel_cross <- setPheno(pop1_sel_cross, h2 = 0.9, simParam = SP)
    
    pop1_sel2 <- selectInd(pop1_sel_cross, nInd = 10, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
    pop1_sel2_cross <- randCross(pop1_sel2, 10, nProgeny = 50, simParam = SP)
    pop1_sel2_cross <- setPheno(pop1_sel2_cross, h2 = 0.9, simParam = SP)
    
    pop1_sel3 <- selectInd(pop1_sel2_cross, nInd = 10, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
    pop1_sel3_cross <- randCross(pop1_sel3, 10, nProgeny = 50, simParam = SP)
    pop1_sel3_cross <- setPheno(pop1_sel3_cross, h2 = 0.9, simParam = SP)
    
    pop1_sel4 <- selectInd(pop1_sel3_cross, nInd = 10, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
    pop1_sel4_cross <- randCross(pop1_sel4, 10, nProgeny = 50, simParam = SP)
    pop1_sel4_cross <- setPheno(pop1_sel4_cross, h2 = 0.9, simParam = SP)
    
    pop1_sel5 <- selectInd(pop1_sel4_cross, nInd = 10, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
    pop1_sel5_cross <- randCross(pop1_sel5, 10, nProgeny = 50, simParam = SP)
    pop1_sel5_cross <- setPheno(pop1_sel5_cross, h2 = 0.9, simParam = SP)
    
    pop1_sel6 <- selectInd(pop1_sel5_cross, nInd = 10, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
    pop1_sel6_cross <- randCross(pop1_sel6, 10, nProgeny = 50, simParam = SP)
    pop1_sel6_cross <- setPheno(pop1_sel6_cross, h2 = 0.9, simParam = SP)
    
    pop1_sel7 <- selectInd(pop1_sel6_cross, nInd = 10, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
    pop1_sel7_cross <- randCross(pop1_sel7, 10, nProgeny = 50, simParam = SP)
    pop1_sel7_cross <- setPheno(pop1_sel7_cross, h2 = 0.9, simParam = SP)
    
    pop1_sel8 <- selectInd(pop1_sel7_cross, nInd = 10, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
    pop1_sel8_cross <- randCross(pop1_sel8, 10, nProgeny = 50, simParam = SP)
    pop1_sel8_cross <- setPheno(pop1_sel8_cross, h2 = 0.9, simParam = SP)
    
    pop1_sel9 <- selectInd(pop1_sel8_cross, nInd = 10, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
    pop1_sel9_cross <- randCross(pop1_sel9, 10, nProgeny = 50, simParam = SP)
    pop1_sel9_cross <- setPheno(pop1_sel9_cross, h2 = 0.9, simParam = SP)
    
    pop1_sel10 <- selectInd(pop1_sel9_cross, nInd = 10, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
    pop1_sel10_cross <- randCross(pop1_sel10, 10, nProgeny = 50, simParam = SP)
    pop1_sel10_cross <- setPheno(pop1_sel10_cross, h2 = 0.9, simParam = SP)
    
    pop1_sel11 <- selectInd(pop1_sel10_cross, nInd = 10, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
    pop1_sel11_cross <- randCross(pop1_sel11, 10, nProgeny = 50, simParam = SP)
    pop1_sel11_cross <- setPheno(pop1_sel11_cross, h2 = 0.9, simParam = SP)
    
    pop1_sel12 <- selectInd(pop1_sel11_cross, nInd = 10, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
    pop1_sel12_cross <- randCross(pop1_sel12, 10, nProgeny = 50, simParam = SP)
    pop1_sel12_cross <- setPheno(pop1_sel12_cross, h2 = 0.9, simParam = SP)
    
    pop1_sel13 <- selectInd(pop1_sel12_cross, nInd = 10, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
    pop1_sel13_cross <- randCross(pop1_sel13, 10, nProgeny = 50, simParam = SP)
    pop1_sel13_cross <- setPheno(pop1_sel13_cross, h2 = 0.9, simParam = SP)
    
    pop1_sel14 <- selectInd(pop1_sel13_cross, nInd = 10, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
    pop1_sel14_cross <- randCross(pop1_sel14, 10, nProgeny = 50, simParam = SP)
    pop1_sel14_cross <- setPheno(pop1_sel14_cross, h2 = 0.9, simParam = SP)
    
    pop1_sel15 <- selectInd(pop1_sel14_cross, nInd = 10, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
    pop1_sel15_cross <- randCross(pop1_sel15, 10, nProgeny = 50, simParam = SP)
    pop1_sel15_cross <- setPheno(pop1_sel15_cross, h2 = 0.9, simParam = SP)
    
    pop1_sel16 <- selectInd(pop1_sel15_cross, nInd = 10, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
    pop1_sel16_cross <- randCross(pop1_sel16, 10, nProgeny = 50, simParam = SP)
    pop1_sel16_cross <- setPheno(pop1_sel16_cross, h2 = 0.9, simParam = SP)
    
    pop1_sel17 <- selectInd(pop1_sel16_cross, nInd = 10, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
    pop1_sel17_cross <- randCross(pop1_sel17, 10, nProgeny = 50, simParam = SP)
    pop1_sel17_cross <- setPheno(pop1_sel17_cross, h2 = 0.9, simParam = SP)
    
    pop1_sel18 <- selectInd(pop1_sel17_cross, nInd = 10, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
    pop1_sel18_cross <- randCross(pop1_sel18, 10, nProgeny = 50, simParam = SP)
    pop1_sel18_cross <- setPheno(pop1_sel18_cross, h2 = 0.9, simParam = SP)
    
    pop1_sel19 <- selectInd(pop1_sel18_cross, nInd = 10, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
    pop1_sel19_cross <- randCross(pop1_sel19, 10, nProgeny = 50, simParam = SP)
    pop1_sel19_cross <- setPheno(pop1_sel19_cross, h2 = 0.9, simParam = SP)
    
    pop1_sel20 <- selectInd(pop1_sel19_cross, nInd = 10, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
    pop1_sel20_cross <- randCross(pop1_sel20, 10, nProgeny = 50, simParam = SP)
    pop1_sel20_cross <- setPheno(pop1_sel20_cross, h2 = 0.9, simParam = SP)
    
    final_sel <- selectInd(pop1_sel20_cross, nInd = 10, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
    pop18_pheno[i,]<- gv(final_sel)
  }

mega_pop <- rbind(pop1_pheno, pop2_pheno, pop3_pheno, pop4_pheno, pop5_pheno, pop6_pheno, pop7_pheno, pop8_pheno, pop9_pheno, pop10_pheno, pop11_pheno, pop12_pheno, pop13_pheno, pop14_pheno, pop15_pheno, pop16_pheno, pop17_pheno, pop18_pheno)
mega_pop <- na.omit(mega_pop)

#Creating confidence intervals
pop1_mean <- mean(mega_pop)
pop1_sd <- sd(mega_pop)
pop1_size <- founderPop@nInd
pop1_se <- pop1_sd/sqrt(pop1_size)
alpha = 0.01
degrees.freedom = pop1_size - 1
t.score = qt(p=alpha/2, df=degrees.freedom,lower.tail=F)
margin.error <- t.score * pop1_se
lower.bound <- pop1_mean - margin.error
upper.bound <- pop1_mean + margin.error
print(c(lower.bound,upper.bound))

#Plotting gv on histogram
hist(mega_pop)
library(MASS)
write.matrix(mega_pop,file="uniform_genespace_gv.csv")
#Plotting confidence intervals
library(ggplot2)
pop1_mean
data <- data.frame(x = 1,
                   y = pop1_mean,
                   lower = lower.bound,
                   upper = upper.bound)
p <- ggplot(data, aes(x, y)) +        # ggplot2 plot with confidence intervals
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper))
p + labs(title = "99% Confidence Interval for Population Means", x="population", y="population mean")
