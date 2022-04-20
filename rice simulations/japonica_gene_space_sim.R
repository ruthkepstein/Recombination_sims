library(AlphaSimR)
library(Rcpp)
library(ggplot2)
library(dplyr)
#setwd("C:/Users/16192/Documents/PNAS_Simulations")
set.seed(420)

##reading in SNPs from B73xMo17 based on v4 B73 ref
jap_snps <- read.table("japonica_SNPs.bed", header =FALSE)
colnames(jap_snps) <- c("Chr#", "SNP Start", "SNP End")
#sample SNPs?
jap_snps <- sample_n(jap_snps, 4000)
jap_snps <- jap_snps[order(jap_snps$`Chr#`,jap_snps$`SNP Start`),]

#Japonica
jap_chr1_snp <- jap_snps[ which(jap_snps$`Chr#` == "Chr1"),]
jap_chr1_snp$rate <- NA
jap_chr1_snp$`SNP End` <- jap_chr1_snp$`SNP End` - min(jap_chr1_snp$`SNP Start`)
jap_chr1_snp$`SNP Start` <- jap_chr1_snp$`SNP Start`- min(jap_chr1_snp$`SNP Start`)

jap_chr2_snp <- jap_snps[ which(jap_snps$`Chr#` == "Chr2"),]
jap_chr2_snp$rate <- NA
jap_chr2_snp$`SNP End` <- jap_chr2_snp$`SNP End` - min(jap_chr2_snp$`SNP Start`)
jap_chr2_snp$`SNP Start` <- jap_chr2_snp$`SNP Start`- min(jap_chr2_snp$`SNP Start`)

jap_chr3_snp <- jap_snps[ which(jap_snps$`Chr#` == "Chr3"),]
jap_chr3_snp$rate <- NA
jap_chr3_snp$`SNP End` <- jap_chr3_snp$`SNP End` - min(jap_chr3_snp$`SNP Start`)
jap_chr3_snp$`SNP Start` <- jap_chr3_snp$`SNP Start`- min(jap_chr3_snp$`SNP Start`)

jap_chr4_snp <- jap_snps[ which(jap_snps$`Chr#` == "Chr4"),]
jap_chr4_snp$rate <- NA
jap_chr4_snp$`SNP End` <- jap_chr4_snp$`SNP End` - min(jap_chr4_snp$`SNP Start`)
jap_chr4_snp$`SNP Start` <- jap_chr4_snp$`SNP Start`- min(jap_chr4_snp$`SNP Start`)

jap_chr5_snp <- jap_snps[ which(jap_snps$`Chr#` == "Chr5"),]
jap_chr5_snp$rate <- NA
jap_chr5_snp$`SNP End` <- jap_chr5_snp$`SNP End` - min(jap_chr5_snp$`SNP Start`)
jap_chr5_snp$`SNP Start` <- jap_chr5_snp$`SNP Start`- min(jap_chr5_snp$`SNP Start`)

jap_chr6_snp <- jap_snps[ which(jap_snps$`Chr#` == "Chr6"),]
jap_chr6_snp$rate <- NA
jap_chr6_snp$`SNP End` <- jap_chr6_snp$`SNP End` - min(jap_chr6_snp$`SNP Start`)
jap_chr6_snp$`SNP Start` <- jap_chr6_snp$`SNP Start`- min(jap_chr6_snp$`SNP Start`)

jap_chr7_snp <- jap_snps[ which(jap_snps$`Chr#` == "Chr7"),]
jap_chr7_snp$rate <- NA
jap_chr7_snp$`SNP End` <- jap_chr7_snp$`SNP End` - min(jap_chr7_snp$`SNP Start`)
jap_chr7_snp$`SNP Start` <- jap_chr7_snp$`SNP Start`- min(jap_chr7_snp$`SNP Start`)

jap_chr8_snp <- jap_snps[ which(jap_snps$`Chr#` == "Chr8"),]
jap_chr8_snp$rate <- NA
jap_chr8_snp$`SNP End` <- jap_chr8_snp$`SNP End` - min(jap_chr8_snp$`SNP Start`)
jap_chr8_snp$`SNP Start` <- jap_chr8_snp$`SNP Start`- min(jap_chr8_snp$`SNP Start`)

jap_chr9_snp <- jap_snps[ which(jap_snps$`Chr#` == "Chr9"),]
jap_chr9_snp$rate <- NA
jap_chr9_snp$`SNP End` <- jap_chr9_snp$`SNP End` - min(jap_chr9_snp$`SNP Start`)
jap_chr9_snp$`SNP Start` <- jap_chr9_snp$`SNP Start`- min(jap_chr9_snp$`SNP Start`)

jap_chr10_snp <- jap_snps[ which(jap_snps$`Chr#` == "Chr10"),]
jap_chr10_snp$rate <- NA
jap_chr10_snp$`SNP End` <- jap_chr10_snp$`SNP End` - min(jap_chr10_snp$`SNP Start`)
jap_chr10_snp$`SNP Start` <- jap_chr10_snp$`SNP Start`- min(jap_chr10_snp$`SNP Start`)

jap_chr11_snp <- jap_snps[ which(jap_snps$`Chr#` == "Chr11"),]
jap_chr11_snp$rate <- NA
jap_chr11_snp$`SNP End` <- jap_chr11_snp$`SNP End` - min(jap_chr11_snp$`SNP Start`)
jap_chr11_snp$`SNP Start` <- jap_chr11_snp$`SNP Start`- min(jap_chr11_snp$`SNP Start`)

jap_chr12_snp <- jap_snps[ which(jap_snps$`Chr#` == "Chr12"),]
jap_chr12_snp$rate <- NA
jap_chr12_snp$`SNP End` <- jap_chr12_snp$`SNP End` - min(jap_chr12_snp$`SNP Start`)
jap_chr12_snp$`SNP Start` <- jap_chr12_snp$`SNP Start`- min(jap_chr12_snp$`SNP Start`)


##Reading in recombination rates
jap_CO <- read.table("japonica_rec_rate.bed", header = FALSE)
colnames(jap_CO) <- c("Chr", "CO Start", "CO End", "rate")
jap_CO <- jap_CO[order(jap_CO$Chr,jap_CO$`CO Start`),]

jap_chr1_CO <- jap_CO[ which(jap_CO$Chr == "chr01"),]
jap_chr1_CO$midpoint <- (jap_chr1_CO$`CO Start`+ jap_chr1_CO$`CO End`)/2
jap_chr1_CO <- jap_chr1_CO[order(jap_chr1_CO$`CO Start`),]

jap_chr2_CO <- jap_CO[ which(jap_CO$Chr == "chr02"),]
jap_chr2_CO$midpoint <- (jap_chr2_CO$`CO Start`+ jap_chr2_CO$`CO End`)/2
jap_chr2_CO <- jap_chr2_CO[order(jap_chr2_CO$`CO Start`),]

jap_chr3_CO <- jap_CO[ which(jap_CO$Chr == "chr03"),]
jap_chr3_CO$midpoint <- (jap_chr3_CO$`CO Start`+ jap_chr3_CO$`CO End`)/2
jap_chr3_CO <- jap_chr3_CO[order(jap_chr3_CO$`CO Start`),]

jap_chr4_CO <- jap_CO[ which(jap_CO$Chr == "chr04"),]
jap_chr4_CO$midpoint <- (jap_chr4_CO$`CO Start`+ jap_chr4_CO$`CO End`)/2
jap_chr4_CO <- jap_chr4_CO[order(jap_chr4_CO$`CO Start`),]

jap_chr5_CO <- jap_CO[ which(jap_CO$Chr == "chr05"),]
jap_chr5_CO$midpoint <- (jap_chr5_CO$`CO Start`+ jap_chr5_CO$`CO End`)/2
jap_chr5_CO <- jap_chr5_CO[order(jap_chr5_CO$`CO Start`),]

jap_chr6_CO <- jap_CO[ which(jap_CO$Chr == "chr06"),]
jap_chr6_CO$midpoint <- (jap_chr6_CO$`CO Start`+ jap_chr6_CO$`CO End`)/2
jap_chr6_CO <- jap_chr6_CO[order(jap_chr6_CO$`CO Start`),]

jap_chr7_CO <- jap_CO[ which(jap_CO$Chr == "chr07"),]
jap_chr7_CO$midpoint <- (jap_chr7_CO$`CO Start`+ jap_chr7_CO$`CO End`)/2
jap_chr7_CO <- jap_chr7_CO[order(jap_chr7_CO$`CO Start`),]

jap_chr8_CO <- jap_CO[ which(jap_CO$Chr == "chr08"),]
jap_chr8_CO$midpoint <- (jap_chr8_CO$`CO Start`+ jap_chr8_CO$`CO End`)/2
jap_chr8_CO <- jap_chr8_CO[order(jap_chr8_CO$`CO Start`),]

jap_chr9_CO <- jap_CO[ which(jap_CO$Chr == "chr09"),]
jap_chr9_CO$midpoint <- (jap_chr9_CO$`CO Start`+ jap_chr9_CO$`CO End`)/2
jap_chr9_CO <- jap_chr9_CO[order(jap_chr9_CO$`CO Start`),]

jap_chr10_CO <- jap_CO[ which(jap_CO$Chr == "chr10"),]
jap_chr10_CO$midpoint <- (jap_chr10_CO$`CO Start`+ jap_chr10_CO$`CO End`)/2
jap_chr10_CO <- jap_chr10_CO[order(jap_chr10_CO$`CO Start`),]

jap_chr11_CO <- jap_CO[ which(jap_CO$Chr == "chr11"),]
jap_chr11_CO$midpoint <- (jap_chr11_CO$`CO Start`+ jap_chr11_CO$`CO End`)/2
jap_chr11_CO <- jap_chr11_CO[order(jap_chr11_CO$`CO Start`),]

jap_chr12_CO <- jap_CO[ which(jap_CO$Chr == "chr12"),]
jap_chr12_CO$midpoint <- (jap_chr12_CO$`CO Start`+ jap_chr12_CO$`CO End`)/2
jap_chr12_CO <- jap_chr12_CO[order(jap_chr12_CO$`CO Start`),]
jap_chr1_CO$`CO End` <- jap_chr1_CO$`CO End` - min(jap_chr1_CO$`CO Start`)
jap_chr1_CO$`CO Start` <- jap_chr1_CO$`CO Start` - min(jap_chr1_CO$`CO Start`)

jap_chr2_CO$`CO End` <- jap_chr2_CO$`CO End` - min(jap_chr2_CO$`CO Start`)
jap_chr2_CO$`CO Start` <- jap_chr2_CO$`CO Start` - min(jap_chr2_CO$`CO Start`)

jap_chr3_CO$`CO End` <- jap_chr3_CO$`CO End` - min(jap_chr3_CO$`CO Start`)
jap_chr3_CO$`CO Start` <- jap_chr3_CO$`CO Start` - min(jap_chr3_CO$`CO Start`)

jap_chr4_CO$`CO End` <- jap_chr4_CO$`CO End` - min(jap_chr4_CO$`CO Start`)
jap_chr4_CO$`CO Start` <- jap_chr4_CO$`CO Start` - min(jap_chr4_CO$`CO Start`)

jap_chr5_CO$`CO End` <- jap_chr5_CO$`CO End` - min(jap_chr5_CO$`CO Start`)
jap_chr5_CO$`CO Start` <- jap_chr5_CO$`CO Start` - min(jap_chr5_CO$`CO Start`)

jap_chr6_CO$`CO End` <- jap_chr6_CO$`CO End` - min(jap_chr6_CO$`CO Start`)
jap_chr6_CO$`CO Start` <- jap_chr6_CO$`CO Start` - min(jap_chr6_CO$`CO Start`)

jap_chr7_CO$`CO End` <- jap_chr7_CO$`CO End` - min(jap_chr7_CO$`CO Start`)
jap_chr7_CO$`CO Start` <- jap_chr7_CO$`CO Start` - min(jap_chr7_CO$`CO Start`)

jap_chr8_CO$`CO End` <- jap_chr8_CO$`CO End` - min(jap_chr8_CO$`CO Start`)
jap_chr8_CO$`CO Start` <- jap_chr8_CO$`CO Start` - min(jap_chr8_CO$`CO Start`)

jap_chr9_CO$`CO End` <- jap_chr9_CO$`CO End` - min(jap_chr9_CO$`CO Start`)
jap_chr9_CO$`CO Start` <- jap_chr9_CO$`CO Start` - min(jap_chr9_CO$`CO Start`)

jap_chr10_CO$`CO End` <- jap_chr10_CO$`CO End` - min(jap_chr10_CO$`CO Start`)
jap_chr10_CO$`CO Start` <- jap_chr10_CO$`CO Start` - min(jap_chr10_CO$`CO Start`)

jap_chr11_CO$`CO End` <- jap_chr11_CO$`CO End` - min(jap_chr11_CO$`CO Start`)
jap_chr11_CO$`CO Start` <- jap_chr11_CO$`CO Start` - min(jap_chr11_CO$`CO Start`)

jap_chr12_CO$`CO End` <- jap_chr12_CO$`CO End` - min(jap_chr12_CO$`CO Start`)
jap_chr12_CO$`CO Start` <- jap_chr12_CO$`CO Start` - min(jap_chr12_CO$`CO Start`)

###using CO rate to infer genetic map distances

##calculating recombination rate per bin of CO data
library(dlookr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(OneR)

#recombination frequency calc used:
# recomb. freq. = (# of COs/ size of population *100%)/ length of bin in Mb

#bin rates into ~1 Mb regions and average each region
fill_start<- function(chr_CO){
  l<-0
  for(k in 1:nrow(chr_CO)){
    if(isFALSE(is.na(chr_CO$rates[k]))){
      temp<-chr_CO$`CO Start`[k]
      chr_CO$`CO Start`[k] <-l
      l<-temp
    }
  }
  print(chr_CO)
}

library(zoo)
jap_chr1_CO_2 <- jap_chr1_CO
bins<-as.integer(max(jap_chr1_CO$`CO End`)/100000)
jap_chr1_CO_2$rates<- rollapply(jap_chr1_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
jap_chr1_CO_2<-fill_start(jap_chr1_CO_2)
jap_chr1_CO_2<- jap_chr1_CO_2 %>% drop_na(rates)

jap_chr2_CO_2 <- jap_chr2_CO
bins<-as.integer(max(jap_chr2_CO$`CO End`)/100000)
jap_chr2_CO_2$rates<- rollapply(jap_chr2_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
jap_chr2_CO_2<-fill_start(jap_chr2_CO_2)
jap_chr2_CO_2<- jap_chr2_CO_2 %>% drop_na(rates)

jap_chr3_CO_2 <- jap_chr3_CO
bins<-as.integer(max(jap_chr3_CO$`CO End`)/100000)
jap_chr3_CO_2$rates<- rollapply(jap_chr3_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
jap_chr3_CO_2<-fill_start(jap_chr3_CO_2)
jap_chr3_CO_2<- jap_chr3_CO_2 %>% drop_na(rates)

jap_chr4_CO_2 <- jap_chr4_CO
bins<-as.integer(max(jap_chr4_CO$`CO End`)/100000)
jap_chr4_CO_2$rates<- rollapply(jap_chr4_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
jap_chr4_CO_2<-fill_start(jap_chr4_CO_2)
jap_chr4_CO_2<- jap_chr4_CO_2 %>% drop_na(rates)

jap_chr5_CO_2 <- jap_chr5_CO
bins<-as.integer(max(jap_chr5_CO$`CO End`)/100000)
jap_chr5_CO_2$rates<- rollapply(jap_chr5_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
jap_chr5_CO_2<-fill_start(jap_chr5_CO_2)
jap_chr5_CO_2<- jap_chr5_CO_2 %>% drop_na(rates)

jap_chr6_CO_2 <- jap_chr6_CO
bins<-as.integer(max(jap_chr6_CO$`CO End`)/100000)
jap_chr6_CO_2$rates<- rollapply(jap_chr6_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
jap_chr6_CO_2<-fill_start(jap_chr6_CO_2)
jap_chr6_CO_2<- jap_chr6_CO_2 %>% drop_na(rates)

jap_chr7_CO_2 <- jap_chr7_CO
bins<-as.integer(max(jap_chr7_CO$`CO End`)/100000)
jap_chr7_CO_2$rates<- rollapply(jap_chr7_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
jap_chr7_CO_2<-fill_start(jap_chr7_CO_2)
jap_chr7_CO_2<- jap_chr7_CO_2 %>% drop_na(rates)

jap_chr8_CO_2 <- jap_chr8_CO
bins<-as.integer(max(jap_chr8_CO$`CO End`)/100000)
jap_chr8_CO_2$rates<- rollapply(jap_chr8_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
jap_chr8_CO_2<-fill_start(jap_chr8_CO_2)
jap_chr8_CO_2<- jap_chr8_CO_2 %>% drop_na(rates)

jap_chr9_CO_2 <- jap_chr9_CO
bins<-as.integer(max(jap_chr9_CO$`CO End`)/100000)
jap_chr9_CO_2$rates<- rollapply(jap_chr9_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
jap_chr9_CO_2<-fill_start(jap_chr9_CO_2)
jap_chr9_CO_2<- jap_chr9_CO_2 %>% drop_na(rates)

jap_chr10_CO_2 <- jap_chr10_CO
bins<-as.integer(max(jap_chr10_CO$`CO End`)/100000)
jap_chr10_CO_2$rates<- rollapply(jap_chr10_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
jap_chr10_CO_2<-fill_start(jap_chr10_CO_2)
jap_chr10_CO_2<- jap_chr10_CO_2 %>% drop_na(rates)

jap_chr11_CO_2 <- jap_chr11_CO
bins<-as.integer(max(jap_chr11_CO$`CO End`)/100000)
jap_chr11_CO_2$rates<- rollapply(jap_chr11_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
jap_chr11_CO_2<-fill_start(jap_chr11_CO_2)
jap_chr11_CO_2<- jap_chr11_CO_2 %>% drop_na(rates)

jap_chr12_CO_2 <- jap_chr12_CO
bins<-as.integer(max(jap_chr12_CO$`CO End`)/100000)
jap_chr12_CO_2$rates<- rollapply(jap_chr12_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
jap_chr12_CO_2<-fill_start(jap_chr12_CO_2)
jap_chr12_CO_2<- jap_chr12_CO_2 %>% drop_na(rates)

##assigning frequency to SNPs based on frequency in each bin
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
##Looking at gene density along chromosomes, binned based on gene density
#ref <- read.table("Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.1.gff3", header = TRUE)
#genome annnotation file from Nipponbare (https://rapdb.dna.affrc.go.jp/download/irgsp1.html)
ref_genes <- read.table("locus.csv", header = TRUE, sep =",")
ref_genes1 <- ref_genes[which(ref_genes$Ref == 'chr01'),]
ref_genes2 <- ref_genes[which(ref_genes$Ref == 'chr02'),]
ref_genes3 <- ref_genes[which(ref_genes$Ref == 'chr03'),]
ref_genes4 <- ref_genes[which(ref_genes$Ref == 'chr04'),]
ref_genes5 <- ref_genes[which(ref_genes$Ref == 'chr05'),]
ref_genes6 <- ref_genes[which(ref_genes$Ref == 'chr06'),]
ref_genes7 <- ref_genes[which(ref_genes$Ref == 'chr07'),]
ref_genes8 <- ref_genes[which(ref_genes$Ref == 'chr08'),]
ref_genes9 <- ref_genes[which(ref_genes$Ref == 'chr09'),]
ref_genes10 <- ref_genes[which(ref_genes$Ref == 'chr10'),]
ref_genes11 <- ref_genes[which(ref_genes$Ref == 'chr11'),]
ref_genes12 <- ref_genes[which(ref_genes$Ref == 'chr12'),]

genes_bin1 <- as.data.frame(summary(binning(ref_genes1$X1, nbins = round(max(ref_genes1$X307041717)/100000), type = "kmeans")))
genes_bin1 <- within(genes_bin1, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin1$levels), ',', fixed=TRUE))))
genes_bin1 <- do.call(data.frame, genes_bin1)
genes_bin1 <- genes_bin1 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin1 <- genes_bin1 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin1[1,4]<-2982
genes_bin1$length <- (genes_bin1$foo.X2-genes_bin1$foo.X1)/100000
genes_bin1[round(max(ref_genes1$X307041717)/100000),5] <- max(ref_genes1$X307041717)
genes_bin1$density <- (genes_bin1$freq/genes_bin1$length)
plot(genes_bin1$foo.X1, genes_bin1$density, type = "l")

genes_bin2 <- as.data.frame(summary(binning(ref_genes2$X1, nbins = round(max(ref_genes2$X307041717)/100000), type = "kmeans")))
genes_bin2 <- within(genes_bin2, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin2$levels), ',', fixed=TRUE))))
genes_bin2 <- do.call(data.frame, genes_bin2)
genes_bin2 <- genes_bin2 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin2 <- genes_bin2 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin2[1,4]<-391
genes_bin2$length <- (genes_bin2$foo.X2-genes_bin2$foo.X1)/100000
genes_bin2[round(max(ref_genes2$X307041717)/100000),5] <- max(ref_genes2$X307041717)
genes_bin2$density <- (genes_bin2$freq/genes_bin2$length)

genes_bin3 <- as.data.frame(summary(binning(ref_genes3$X1, nbins = round(max(ref_genes3$X307041717)/100000), type = "kmeans")))
genes_bin3 <- within(genes_bin3, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin3$levels), ',', fixed=TRUE))))
genes_bin3 <- do.call(data.frame, genes_bin3)
genes_bin3 <- genes_bin3 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin3 <- genes_bin3 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin3[1,4]<-4582
genes_bin3$length <- (genes_bin3$foo.X2-genes_bin3$foo.X1)/100000
genes_bin3[round(max(ref_genes4$X307041717)/100000),5] <- max(ref_genes3$X307041717)
genes_bin3$density <- (genes_bin3$freq/genes_bin3$length)

genes_bin4 <- as.data.frame(summary(binning(ref_genes4$X1, nbins = round(max(ref_genes4$X307041717)/100000), type = "kmeans")))
genes_bin4 <- within(genes_bin4, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin4$levels), ',', fixed=TRUE))))
genes_bin4 <- do.call(data.frame, genes_bin4)
genes_bin4 <- genes_bin4 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin4 <- genes_bin4 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin4[1,4]<-58788
genes_bin4$length <- (genes_bin4$foo.X2-genes_bin4$foo.X1)/100000
genes_bin4[round(max(ref_genes4$X307041717)/100000),5] <- max(ref_genes4$X307041717)
genes_bin4$density <- (genes_bin4$freq/genes_bin4$length)

genes_bin5 <- as.data.frame(summary(binning(ref_genes5$X1, nbins = round(max(ref_genes5$X307041717)/100000), type = "kmeans")))
genes_bin5 <- within(genes_bin5, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin5$levels), ',', fixed=TRUE))))
genes_bin5 <- do.call(data.frame, genes_bin5)
genes_bin5 <- genes_bin5 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin5 <- genes_bin5 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin5[1,4]<-10946
genes_bin5$length <- (genes_bin5$foo.X2-genes_bin5$foo.X1)/100000
genes_bin5[round(max(ref_genes5$X307041717)/100000),5] <- max(ref_genes5$X307041717)
genes_bin5$density <- (genes_bin5$freq/genes_bin5$length)

genes_bin6 <- as.data.frame(summary(binning(ref_genes6$X1, nbins = round(max(ref_genes6$X307041717)/100000), type = "kmeans")))
genes_bin6 <- within(genes_bin6, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin6$levels), ',', fixed=TRUE))))
genes_bin6 <- do.call(data.frame, genes_bin6)
genes_bin6 <- genes_bin6 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin6 <- genes_bin6 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin6[1,4]<-38339
genes_bin6$length <- (genes_bin6$foo.X2-genes_bin6$foo.X1)/100000
genes_bin6[round(max(ref_genes6$X307041717)/100000),5] <- max(ref_genes6$X307041717)
genes_bin6$density <- (genes_bin6$freq/genes_bin6$length)

genes_bin7 <- as.data.frame(summary(binning(ref_genes7$X1, nbins = round(max(ref_genes7$X307041717)/100000), type = "kmeans")))
genes_bin7 <- within(genes_bin7, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin7$levels), ',', fixed=TRUE))))
genes_bin7 <- do.call(data.frame, genes_bin7)
genes_bin7 <- genes_bin7 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin7 <- genes_bin7 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin7[1,4]<-11647
genes_bin7$length <- (genes_bin7$foo.X2-genes_bin7$foo.X1)/100000
genes_bin7[round(max(ref_genes7$X307041717)/100000),5] <- max(ref_genes7$X307041717)
genes_bin7$density <- (genes_bin7$freq/genes_bin7$length)

genes_bin8 <- as.data.frame(summary(binning(ref_genes8$X1, nbins = round(max(ref_genes8$X307041717)/100000), type = "kmeans")))
genes_bin8 <- within(genes_bin8, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin8$levels), ',', fixed=TRUE))))
genes_bin8 <- do.call(data.frame, genes_bin8)
genes_bin8 <- genes_bin8 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin8 <- genes_bin8 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin8[1,4]<-17350
genes_bin8$length <- (genes_bin8$foo.X2-genes_bin8$foo.X1)/100000
genes_bin8[round(max(ref_genes8$X307041717)/100000),5] <- max(ref_genes8$X307041717)
genes_bin8$density <- (genes_bin8$freq/genes_bin8$length)

genes_bin9 <- as.data.frame(summary(binning(ref_genes9$X1, nbins = round(max(ref_genes9$X307041717)/100000), type = "kmeans")))
genes_bin9 <- within(genes_bin9, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin9$levels), ',', fixed=TRUE))))
genes_bin9 <- do.call(data.frame, genes_bin9)
genes_bin9 <- genes_bin9 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin9 <- genes_bin9 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin9[1,4]<-145180
genes_bin9$length <- (genes_bin9$foo.X2-genes_bin9$foo.X1)/100000
genes_bin9[round(max(ref_genes9$X307041717)/100000),5] <- max(ref_genes9$X307041717)
genes_bin9$density <- (genes_bin9$freq/genes_bin9$length)

genes_bin10 <- binning(ref_genes10$X1, nbins = round(max(ref_genes10$X307041717)/100000), type = "kmeans")
genes_bin10 <- as.data.frame(summary(binning(ref_genes10$X1, nbins = round(max(ref_genes10$X307041717)/100000), type = "kmeans")))
genes_bin10 <- within(genes_bin10, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin10$levels), ',', fixed=TRUE))))
genes_bin10 <- do.call(data.frame, genes_bin10)
genes_bin10 <- genes_bin10 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin10 <- genes_bin10 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin10[1,4]<-44901
genes_bin10$length <- (genes_bin10$foo.X2-genes_bin10$foo.X1)/100000
genes_bin10[round(max(ref_genes10$X307041717)/100000),5] <- max(ref_genes10$X307041717)
genes_bin10$density <- (genes_bin10$freq/genes_bin10$length)

genes_bin11 <- binning(ref_genes11$X1, nbins = round(max(ref_genes11$X307041717)/100000), type = "kmeans")
genes_bin11 <- as.data.frame(summary(binning(ref_genes11$X1, nbins = round(max(ref_genes11$X307041717)/100000), type = "kmeans")))
genes_bin11 <- within(genes_bin11, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin11$levels), ',', fixed=TRUE))))
genes_bin11 <- do.call(data.frame, genes_bin11)
genes_bin11 <- genes_bin11 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin11 <- genes_bin11 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin11[1,4]<-1179
genes_bin11$length <- (genes_bin11$foo.X2-genes_bin11$foo.X1)/100000
genes_bin11[round(max(ref_genes11$X307041717)/100000),5] <- max(ref_genes11$X307041717)
genes_bin11$density <- (genes_bin11$freq/genes_bin11$length)

genes_bin12 <- binning(ref_genes12$X1, nbins = round(max(ref_genes12$X307041717)/100000), type = "kmeans")
genes_bin12 <- as.data.frame(summary(binning(ref_genes12$X1, nbins = round(max(ref_genes12$X307041717)/100000), type = "kmeans")))
genes_bin12 <- within(genes_bin12, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin12$levels), ',', fixed=TRUE))))
genes_bin12 <- do.call(data.frame, genes_bin12)
genes_bin12 <- genes_bin12 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin12 <- genes_bin12 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin12[1,4]<-2681
genes_bin12$length <- (genes_bin12$foo.X2-genes_bin12$foo.X1)/100000
genes_bin12[round(max(ref_genes12$X307041717)/100000),5] <- max(ref_genes12$X307041717)
genes_bin12$density <- (genes_bin12$freq/genes_bin12$length)


##Spearmen correlation test to find correlation between gene density & recombination rate
chr1_corr <- genes_bin1$density
chr1_corr <- as.data.frame(chr1_corr)
colnames(chr1_corr) <- "density"
chr1_corr$start <- genes_bin1$foo.X1
chr1_corr$end <- genes_bin1$foo.X2
chr1_corr$rate <- NA

chr2_corr <- genes_bin2$density
chr2_corr <- as.data.frame(chr2_corr)
colnames(chr2_corr) <- "density"
chr2_corr$start <- genes_bin2$foo.X1
chr2_corr$end <- genes_bin2$foo.X2
chr2_corr$rate <- NA

chr3_corr <- genes_bin3$density
chr3_corr <- as.data.frame(chr3_corr)
colnames(chr3_corr) <- "density"
chr3_corr$start <- genes_bin3$foo.X1
chr3_corr$end <- genes_bin3$foo.X2
chr3_corr$rate <- NA

chr4_corr <- genes_bin4$density
chr4_corr <- as.data.frame(chr4_corr)
colnames(chr4_corr) <- "density"
chr4_corr$start <- genes_bin4$foo.X1
chr4_corr$end <- genes_bin4$foo.X2
chr4_corr$rate <- NA

chr5_corr <- genes_bin5$density
chr5_corr <- as.data.frame(chr5_corr)
colnames(chr5_corr) <- "density"
chr5_corr$start <- genes_bin5$foo.X1
chr5_corr$end <- genes_bin5$foo.X2
chr5_corr$rate <- NA

chr6_corr <- genes_bin6$density
chr6_corr <- as.data.frame(chr6_corr)
colnames(chr6_corr) <- "density"
chr6_corr$start <- genes_bin6$foo.X1
chr6_corr$end <- genes_bin6$foo.X2
chr6_corr$rate <- NA

chr7_corr <- genes_bin7$density
chr7_corr <- as.data.frame(chr7_corr)
colnames(chr7_corr) <- "density"
chr7_corr$start <- genes_bin7$foo.X1
chr7_corr$end <- genes_bin7$foo.X2
chr7_corr$rate <- NA

chr8_corr <- genes_bin8$density
chr8_corr <- as.data.frame(chr8_corr)
colnames(chr8_corr) <- "density"
chr8_corr$start <- genes_bin8$foo.X1
chr8_corr$end <- genes_bin8$foo.X2
chr8_corr$rate <- NA

chr9_corr <- genes_bin9$density
chr9_corr <- as.data.frame(chr9_corr)
colnames(chr9_corr) <- "density"
chr9_corr$start <- genes_bin9$foo.X1
chr9_corr$end <- genes_bin9$foo.X2
chr9_corr$rate <- NA

chr10_corr <- genes_bin10$density
chr10_corr <- as.data.frame(chr10_corr)
colnames(chr10_corr) <- "density"
chr10_corr$start <- genes_bin10$foo.X1
chr10_corr$end <- genes_bin10$foo.X2
chr10_corr$rate <- NA

chr11_corr <- genes_bin11$density
chr11_corr <- as.data.frame(chr11_corr)
colnames(chr11_corr) <- "density"
chr11_corr$start <- genes_bin11$foo.X1
chr11_corr$end <- genes_bin11$foo.X2
chr11_corr$rate <- NA

chr12_corr <- genes_bin12$density
chr12_corr <- as.data.frame(chr12_corr)
colnames(chr12_corr) <- "density"
chr12_corr$start <- genes_bin12$foo.X1
chr12_corr$end <- genes_bin12$foo.X2
chr12_corr$rate <- NA

assign_rate <- function(chr_bin, chr_corr){
  for(i in 1:nrow(chr_bin)){
    for(k in 1:nrow(chr_corr)){
      if(isTRUE((chr_bin$`CO Start`[i] <= chr_corr$end[k]) && (chr_bin$`CO End`[i] >= chr_corr$end[k]))){
        chr_corr$rate[k] <- chr_bin$rate[i]
      }
    }
  }
  return(chr_corr)
}
#200kb correlation
chr1_corr_rate <- assign_rate(jap_chr1_CO_2, chr1_corr)
chr1_corr_rate <- na.omit(chr1_corr_rate)
cor.test(chr1_corr_rate$rate, chr1_corr_rate$density,  method = "spearman", alternative = "greater")

chr2_corr_rate <- assign_rate(jap_chr2_CO_2, chr2_corr)
chr2_corr_rate <-  na.omit(chr2_corr_rate)
cor.test(chr2_corr_rate$rate, chr2_corr_rate$density,  method = "spearman", alternative = "greater")

chr3_corr_rate <- assign_rate(jap_chr3_CO_2, chr3_corr)
chr3_corr_rate <- na.omit(chr3_corr_rate)
cor.test(chr3_corr_rate$rate, chr3_corr_rate$density,  method = "spearman", alternative = "greater")

chr4_corr_rate <- assign_rate(jap_chr4_CO_2, chr4_corr)
chr4_corr_rate <- na.omit(chr4_corr_rate)
cor.test(chr4_corr_rate$rate, chr4_corr_rate$density,  method = "spearman", alternative = "greater")

chr5_corr_rate <- assign_rate(jap_chr5_CO_2, chr5_corr)
chr5_corr_rate <- na.omit(chr5_corr_rate)
cor.test(chr5_corr_rate$rate, chr5_corr_rate$density,  method = "spearman", alternative = "greater")

chr6_corr_rate <- assign_rate(jap_chr6_CO_2, chr6_corr)
chr6_corr_rate <- na.omit(chr6_corr_rate)
cor.test(chr6_corr_rate$rate, chr6_corr_rate$density,  method = "spearman", alternative = "greater")

chr7_corr_rate <- assign_rate(jap_chr7_CO_2, chr7_corr)
chr7_corr_rate <-na.omit(chr7_corr_rate)
cor.test(chr7_corr_rate$rate, chr7_corr_rate$density,  method = "spearman", alternative = "greater")

chr8_corr_rate <- assign_rate(jap_chr8_CO_2, chr8_corr)
chr8_corr_rate <- na.omit(chr8_corr_rate)
cor.test(chr8_corr_rate$rate, chr8_corr_rate$density,  method = "spearman", alternative = "greater")

chr9_corr_rate <- assign_rate(jap_chr9_CO_2, chr9_corr)
chr9_corr_rate <- na.omit(chr9_corr_rate)
cor.test(chr9_corr_rate$rate, chr9_corr_rate$density,  method = "spearman", alternative = "greater")

chr10_corr_rate <- assign_rate(jap_chr10_CO_2, chr10_corr)
chr10_corr_rate <- na.omit(chr10_corr_rate)
cor.test(chr10_corr_rate$rate, chr10_corr_rate$density,  method = "spearman", alternative = "greater")

chr11_corr_rate <- assign_rate(jap_chr11_CO_2, chr11_corr)
chr11_corr_rate <- na.omit(chr11_corr_rate)
cor.test(chr11_corr_rate$rate, chr11_corr_rate$density,  method = "spearman", alternative = "greater")

chr12_corr_rate <- assign_rate(jap_chr12_CO_2, chr12_corr)
chr12_corr_rate <- na.omit(chr12_corr_rate)
cor.test(chr12_corr_rate$rate, chr12_corr_rate$density,  method = "spearman", alternative = "greater")

genomewide <- rbind(chr1_corr_rate, chr2_corr_rate, chr3_corr_rate, chr4_corr_rate, chr5_corr_rate,
                    chr6_corr_rate, chr7_corr_rate, chr8_corr_rate, chr9_corr_rate, chr10_corr_rate)

cor.test(genomewide$rate, genomewide$density,  method = "spearman", alternative = "greater")

#REASSIGN FREQ loop through snp positions, all snp that fall in bins with low gene density get 1/4 of mean rate
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

snp_uniform_rate_chr2 <- function(genes_bin, chr_snp, chr_snps_mean){
  for(i in 1:nrow(chr_snp)){
    for(k in 1:nrow(genes_bin)){
      if(isTRUE((chr_snp$`SNP Start`[i] >= genes_bin$foo.X1[k]) && (chr_snp$`SNP Start`[i] <= genes_bin$foo.X2[k]))){
        if(isTRUE((genes_bin$Quartile[k]==1 || genes_bin$Quartile[k]==2))){
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
chr2_snp2 <-snp_uniform_rate_chr2(genes_bin2, chr2_snp, chr2_snps_mean)
chr2_snp2$`SNP Start` <- chr2_snp2$`SNP Start`/1000000
chr2_snp2 <- chr2_snp2[-(196:237),]
#chr2_spl <- smooth.spline(chr2_snp2$rate, spar = 0.6)
chr2_snp2$pos <- (chr2_snp2$`SNP Start`*chr2_snp2$rate)
plot(chr2_snp2$`SNP Start`, chr2_snp2$pos)
ggplot(chr2_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
plot(chr2_snp2$`SNP Start`, chr2_snp2$pos/chr2_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 2 Recombination Distribution")
chr2_finalpos <- chr2_snp2[order(chr2_snp2$pos),]
is.unsorted(chr2_finalpos$pos)
plot(chr2_snp2$`SNP Start`, chr2_finalpos$pos/chr2_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 2 Recombination Distribution")

chr3_snp2 <- snp_rate(chr3_bin, chr3_snp)
#chr3_snp2$rate[1:73] <- 0.9275442
#chr3_snp2$rate[74:147] <- 0.6183628
#chr3_snp2$rate[148:219] <- 0.9275442
chr3_snps_mean <- mean(chr3_snp2$rate)
genes_bin3$Quartile<-cut(genes_bin3$freq,quantile(genes_bin3$freq),include.lowest=TRUE,labels=FALSE)
chr3_snp2 <-snp_uniform_rate(genes_bin3, chr3_snp, chr3_snps_mean)
chr3_snp2$`SNP Start` <- chr3_snp2$`SNP Start`/1000000
chr3_snp2 <- chr3_snp2[-(205:219),]
chr3_spl <- smooth.spline(chr3_snp2$rate, spar = 1.2)
chr3_snp2$pos <- (chr3_snp2$`SNP Start`*chr3_spl$y)
plot(chr3_snp2$`SNP Start`, chr3_snp2$pos)
ggplot(chr3_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
plot(chr3_snp2$`SNP Start`, chr3_snp2$pos/chr3_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 3 Recombination Distribution")
chr3_finalpos <- chr3_snp2[order(chr3_snp2$pos),]
is.unsorted(chr3_finalpos$pos)
plot(chr3_snp2$`SNP Start`, chr3_finalpos$pos/chr3_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 3 Recombination Distribution")

chr4_snp2 <- snp_rate(chr4_bin, chr4_snp)
#chr4_snp2$rate[1:85] <- 0.7232293
#chr4_snp2$rate[86:170] <- 0.4821528
#chr4_snp2$rate[171:256] <- 0.7232293
chr4_snps_mean <- mean(chr4_snp2$rate)
genes_bin4$Quartile<-cut(genes_bin4$freq,quantile(genes_bin4$freq),include.lowest=TRUE,labels=FALSE)
chr4_snp2 <-snp_uniform_rate(genes_bin4, chr4_snp, chr4_snps_mean)
chr4_snp2$`SNP Start` <- chr4_snp2$`SNP Start`/1000000
chr4_spl <- smooth.spline(chr4_snp2$rate, spar = 1)
chr4_snp2$pos <- (chr4_snp2$`SNP Start`*chr4_spl$y)
plot(chr4_snp2$`SNP Start`, chr4_snp2$pos)
ggplot(chr4_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
plot(chr4_snp2$`SNP Start`, chr4_snp2$pos/chr4_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 4 Recombination Distribution")
chr4_finalpos <- chr4_snp2[order(chr4_snp2$pos),]
is.unsorted(chr4_finalpos$pos)
plot(chr4_snp2$`SNP Start`, chr4_finalpos$pos/chr4_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 4 Recombination Distribution")

chr5_snp2 <- snp_rate(chr5_bin, chr5_snp)
chr5_snps_mean <- mean(chr5_snp2$rate)
genes_bin5$Quartile<-cut(genes_bin5$freq,quantile(genes_bin5$freq),include.lowest=TRUE,labels=FALSE)
chr5_snp2 <-snp_uniform_rate(genes_bin5, chr5_snp, chr5_snps_mean)
chr5_snp2$`SNP Start` <- chr5_snp2$`SNP Start`/1000000
chr5_spl <- smooth.spline(chr5_snp2$rate, spar = 1)
chr5_snp2$pos <- (chr5_snp2$`SNP Start`*chr5_spl$y)
plot(chr5_snp2$`SNP Start`, chr5_snp2$pos)
ggplot(chr5_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
plot(chr5_snp2$`SNP Start`, chr5_snp2$pos/chr5_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 5 Recombination Distribution")
chr5_finalpos <- chr5_snp2[order(chr5_snp2$pos),]
is.unsorted(chr5_finalpos$pos)
plot(chr5_snp2$`SNP Start`, chr5_finalpos$pos/chr5_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 5 Recombination Distribution")

chr6_snp2 <- snp_rate(chr6_bin, chr6_snp)
chr6_snps_mean <- mean(chr6_snp2$rate)
genes_bin6$Quartile<-cut(genes_bin6$freq,quantile(genes_bin6$freq),include.lowest=TRUE,labels=FALSE)
chr6_snp2 <-snp_uniform_rate(genes_bin6, chr6_snp, chr6_snps_mean)
chr6_snp2$`SNP Start` <- chr6_snp2$`SNP Start`/1000000
chr6_spl <- smooth.spline(chr6_snp2$rate, spar = 1.1)
chr6_snp2$pos <- (chr6_snp2$`SNP Start`*chr6_spl$y)
plot(chr6_snp2$`SNP Start`, chr6_snp2$pos)
ggplot(chr6_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
plot(chr6_snp2$`SNP Start`, chr6_snp2$pos/chr6_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 6 Recombination Distribution")
chr6_finalpos <- chr6_snp2[order(chr6_snp2$pos),]
is.unsorted(chr6_finalpos$pos)
plot(chr6_snp2$`SNP Start`, chr6_finalpos$pos/chr6_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 6 Recombination Distribution")

chr7_snp2 <- snp_rate(chr7_bin, chr7_snp)
chr7_snps_mean <- mean(chr7_snp2$rate)
genes_bin7$Quartile<-cut(genes_bin7$freq,quantile(genes_bin7$freq),include.lowest=TRUE,labels=FALSE)
chr7_snp2 <-snp_uniform_rate(genes_bin7, chr7_snp, chr7_snps_mean)
chr7_snp2$`SNP Start` <- chr7_snp2$`SNP Start`/1000000
chr7_spl <- smooth.spline(chr7_snp2$rate, spar = 1)
chr7_snp2$pos <- (chr7_snp2$`SNP Start`*chr7_spl$y)
plot(chr7_snp2$`SNP Start`, chr7_snp2$pos)
ggplot(chr7_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
plot(chr7_snp2$`SNP Start`, chr7_snp2$pos/chr7_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 7 Recombination Distribution")
chr7_finalpos <- chr7_snp2[order(chr7_snp2$pos),]
is.unsorted(chr7_finalpos$pos)
plot(chr7_snp2$`SNP Start`, chr7_finalpos$pos/chr7_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 7 Recombination Distribution")

chr8_snp2 <- snp_rate(chr8_bin, chr8_snp)
chr8_snps_mean <- mean(chr8_snp2$rate)
genes_bin8$Quartile<-cut(genes_bin8$freq,quantile(genes_bin8$freq),include.lowest=TRUE,labels=FALSE)
chr8_snp2 <-snp_uniform_rate(genes_bin8, chr8_snp, chr8_snps_mean)
chr8_snp2$`SNP Start` <- chr8_snp2$`SNP Start`/1000000
chr8_spl <- smooth.spline(chr8_snp2$rate, spar =.9)
chr8_snp2$pos <- (chr8_snp2$`SNP Start`*chr8_spl$y)
plot(chr8_snp2$`SNP Start`, chr8_snp2$pos)
ggplot(chr8_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
plot(chr8_snp2$`SNP Start`, chr8_snp2$pos/chr8_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 8 Recombination Distribution")
chr8_finalpos <- chr8_snp2[order(chr8_snp2$pos),]
is.unsorted(chr8_finalpos$pos)
plot(chr8_snp2$`SNP Start`, chr8_finalpos$pos/chr8_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 8 Recombination Distribution")

chr9_snp2 <- snp_rate(chr9_bin, chr9_snp)
chr9_snps_mean <- mean(chr9_snp2$rate)
genes_bin9$Quartile<-cut(genes_bin9$freq,quantile(genes_bin9$freq),include.lowest=TRUE,labels=FALSE)
chr9_snp2 <-snp_uniform_rate(genes_bin9, chr9_snp, chr9_snps_mean)
chr9_snp2$`SNP Start` <- chr9_snp2$`SNP Start`/1000000
chr9_spl <- smooth.spline(chr9_snp2$rate, spar = 1.1)
chr9_snp2$pos <- (chr9_snp2$`SNP Start`*chr9_spl$y)
plot(chr9_snp2$`SNP Start`, chr9_snp2$pos)
ggplot(chr9_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
plot(chr9_snp2$`SNP Start`, chr9_snp2$pos/chr9_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 9 Recombination Distribution")
chr9_finalpos <- chr9_snp2[order(chr9_snp2$pos),]
is.unsorted(chr9_finalpos$pos)
plot(chr9_snp2$`SNP Start`, chr9_finalpos$pos/chr9_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 9 Recombination Distribution")

chr10_snp2 <- snp_rate(chr10_bin, chr10_snp)
chr10_snps_mean <- mean(chr10_snp2$rate)
genes_bin10$Quartile<-cut(genes_bin10$freq,quantile(genes_bin10$freq),include.lowest=TRUE,labels=FALSE)
chr10_snp2 <-snp_uniform_rate(genes_bin10, chr10_snp, chr10_snps_mean)
chr10_snp2$`SNP Start` <- chr10_snp2$`SNP Start`/1000000
chr10_spl <- smooth.spline(chr10_snp2$rate, spar = 1)
chr10_snp2$pos <- (chr10_snp2$`SNP Start`*chr10_spl$y)
plot(chr10_snp2$`SNP Start`, chr10_snp2$pos)
ggplot(chr10_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
plot(chr10_snp2$`SNP Start`, chr10_snp2$pos/chr10_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 10 Recombination Distribution")
chr10_finalpos <- chr10_snp2[order(chr10_snp2$pos),]
is.unsorted(chr10_finalpos$pos)
plot(chr10_snp2$`SNP Start`, chr10_finalpos$pos/chr10_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 10 Recombination Distribution")

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
