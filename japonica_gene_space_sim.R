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
jap_snps <- sample_n(jap_snps, 2000)
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
bins<-as.integer(nrow(jap_chr1_CO)/40)
jap_chr1_CO_2$rates<- rollapply(jap_chr1_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
jap_chr1_CO_2<-fill_start(jap_chr1_CO_2)
jap_chr1_CO_2<- jap_chr1_CO_2 %>% drop_na(rates)

jap_chr2_CO_2 <- jap_chr2_CO
bins<-as.integer(nrow(jap_chr2_CO)/40)
jap_chr2_CO_2$rates<- rollapply(jap_chr2_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
jap_chr2_CO_2<-fill_start(jap_chr2_CO_2)
jap_chr2_CO_2<- jap_chr2_CO_2 %>% drop_na(rates)

jap_chr3_CO_2 <- jap_chr3_CO
bins<-as.integer(nrow(jap_chr3_CO)/40)
jap_chr3_CO_2$rates<- rollapply(jap_chr3_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
jap_chr3_CO_2<-fill_start(jap_chr3_CO_2)
jap_chr3_CO_2<- jap_chr3_CO_2 %>% drop_na(rates)

jap_chr4_CO_2 <- jap_chr4_CO
bins<-as.integer(nrow(jap_chr4_CO)/40)
jap_chr4_CO_2$rates<- rollapply(jap_chr4_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
jap_chr4_CO_2<-fill_start(jap_chr4_CO_2)
jap_chr4_CO_2<- jap_chr4_CO_2 %>% drop_na(rates)

jap_chr5_CO_2 <- jap_chr5_CO
bins<-as.integer(nrow(jap_chr5_CO)/40)
jap_chr5_CO_2$rates<- rollapply(jap_chr5_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
jap_chr5_CO_2<-fill_start(jap_chr5_CO_2)
jap_chr5_CO_2<- jap_chr5_CO_2 %>% drop_na(rates)

jap_chr6_CO_2 <- jap_chr6_CO
bins<-as.integer(nrow(jap_chr6_CO)/40)
jap_chr6_CO_2$rates<- rollapply(jap_chr6_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
jap_chr6_CO_2<-fill_start(jap_chr6_CO_2)
jap_chr6_CO_2<- jap_chr6_CO_2 %>% drop_na(rates)

jap_chr7_CO_2 <- jap_chr7_CO
bins<-as.integer(nrow(jap_chr7_CO)/40)
jap_chr7_CO_2$rates<- rollapply(jap_chr7_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
jap_chr7_CO_2<-fill_start(jap_chr7_CO_2)
jap_chr7_CO_2<- jap_chr7_CO_2 %>% drop_na(rates)

jap_chr8_CO_2 <- jap_chr8_CO
bins<-as.integer(nrow(jap_chr8_CO)/40)
jap_chr8_CO_2$rates<- rollapply(jap_chr8_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
jap_chr8_CO_2<-fill_start(jap_chr8_CO_2)
jap_chr8_CO_2<- jap_chr8_CO_2 %>% drop_na(rates)

jap_chr9_CO_2 <- jap_chr9_CO
bins<-as.integer(nrow(jap_chr9_CO)/40)
jap_chr9_CO_2$rates<- rollapply(jap_chr9_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
jap_chr9_CO_2<-fill_start(jap_chr9_CO_2)
jap_chr9_CO_2<- jap_chr9_CO_2 %>% drop_na(rates)

jap_chr10_CO_2 <- jap_chr10_CO
bins<-as.integer(nrow(jap_chr10_CO)/40)
jap_chr10_CO_2$rates<- rollapply(jap_chr10_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
jap_chr10_CO_2<-fill_start(jap_chr10_CO_2)
jap_chr10_CO_2<- jap_chr10_CO_2 %>% drop_na(rates)

jap_chr11_CO_2 <- jap_chr11_CO
bins<-as.integer(nrow(jap_chr11_CO)/40)
jap_chr11_CO_2$rates<- rollapply(jap_chr11_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
jap_chr11_CO_2<-fill_start(jap_chr11_CO_2)
jap_chr11_CO_2<- jap_chr11_CO_2 %>% drop_na(rates)

jap_chr12_CO_2 <- jap_chr12_CO
bins<-as.integer(nrow(jap_chr12_CO)/40)
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
genes_bin1$length <- (genes_bin1$foo.X2)-(genes_bin1$foo.X1)
genes_bin1[round(max(ref_genes1$X307041717)/100000),5] <- max(ref_genes1$X307041717)
genes_bin1$density <- (genes_bin1$freq/(genes_bin1$length/1000000))/1000
plot(genes_bin1$foo.X1/1000000, genes_bin1$density, type = "l")
saveRDS(genes_bin1$foo.X1/1000000,file="chr1_locipos")
saveRDS(genes_bin1$density,file="chr1_geneProb")

genes_bin2 <- as.data.frame(summary(binning(ref_genes2$X1, nbins = round(max(ref_genes2$X307041717)/100000), type = "kmeans")))
genes_bin2 <- within(genes_bin2, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin2$levels), ',', fixed=TRUE))))
genes_bin2 <- do.call(data.frame, genes_bin2)
genes_bin2 <- genes_bin2 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin2 <- genes_bin2 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin2[1,4]<-391
genes_bin2$length <- (genes_bin2$foo.X2-genes_bin2$foo.X1)
genes_bin2[round(max(ref_genes2$X307041717)/100000),5] <- max(ref_genes2$X307041717)
genes_bin2$density <- (genes_bin2$freq/(genes_bin2$length/1000000))/1000
plot(genes_bin2$foo.X1/1000000, genes_bin2$density,type="l")
saveRDS(genes_bin2$foo.X1/1000000,file="chr2_locipos")
saveRDS(genes_bin2$density,file="chr2_geneProb")

genes_bin3 <- as.data.frame(summary(binning(ref_genes3$X1, nbins = round(max(ref_genes3$X307041717)/100000), type = "kmeans")))
genes_bin3 <- within(genes_bin3, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin3$levels), ',', fixed=TRUE))))
genes_bin3 <- do.call(data.frame, genes_bin3)
genes_bin3 <- genes_bin3 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin3 <- genes_bin3 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin3[1,4]<-4582
genes_bin3$length <- (genes_bin3$foo.X2-genes_bin3$foo.X1)
genes_bin3[round(max(ref_genes4$X307041717)/100000),5] <- max(ref_genes3$X307041717)
genes_bin3$density <- (genes_bin3$freq/(genes_bin3$length/1000000))/1000
plot(genes_bin3$foo.X1/1000000, genes_bin3$density,type="l")
saveRDS(genes_bin3$foo.X1/1000000,file="chr3_locipos")
saveRDS(genes_bin3$density,file="chr3_geneProb")

genes_bin4 <- as.data.frame(summary(binning(ref_genes4$X1, nbins = round(max(ref_genes4$X307041717)/100000), type = "kmeans")))
genes_bin4 <- within(genes_bin4, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin4$levels), ',', fixed=TRUE))))
genes_bin4 <- do.call(data.frame, genes_bin4)
genes_bin4 <- genes_bin4 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin4 <- genes_bin4 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin4[1,4]<-58788
genes_bin4$length <- (genes_bin4$foo.X2-genes_bin4$foo.X1)
genes_bin4[round(max(ref_genes4$X307041717)/100000),5] <- max(ref_genes4$X307041717)
genes_bin4$density <- (genes_bin4$freq/(genes_bin4$length/1000000))/1000
plot(genes_bin4$foo.X1/1000000, genes_bin4$density,type="l")
saveRDS(genes_bin4$foo.X1/1000000,file="chr4_locipos")
saveRDS(genes_bin4$density,file="chr4_geneProb")

genes_bin5 <- as.data.frame(summary(binning(ref_genes5$X1, nbins = round(max(ref_genes5$X307041717)/100000), type = "kmeans")))
genes_bin5 <- within(genes_bin5, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin5$levels), ',', fixed=TRUE))))
genes_bin5 <- do.call(data.frame, genes_bin5)
genes_bin5 <- genes_bin5 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin5 <- genes_bin5 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin5[1,4]<-10946
genes_bin5$length <- (genes_bin5$foo.X2-genes_bin5$foo.X1)
genes_bin5[round(max(ref_genes5$X307041717)/100000),5] <- max(ref_genes5$X307041717)
genes_bin5$density <- (genes_bin5$freq/(genes_bin5$length/1000000))/1000
plot(genes_bin5$foo.X1/1000000, genes_bin5$density,type="l")
saveRDS(genes_bin5$foo.X1/1000000,file="chr5_locipos")
saveRDS(genes_bin5$density,file="chr5_geneProb")

genes_bin6 <- as.data.frame(summary(binning(ref_genes6$X1, nbins = round(max(ref_genes6$X307041717)/100000), type = "kmeans")))
genes_bin6 <- within(genes_bin6, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin6$levels), ',', fixed=TRUE))))
genes_bin6 <- do.call(data.frame, genes_bin6)
genes_bin6 <- genes_bin6 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin6 <- genes_bin6 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin6[1,4]<-38339
genes_bin6$length <- (genes_bin6$foo.X2-genes_bin6$foo.X1)
genes_bin6[round(max(ref_genes6$X307041717)/100000),5] <- max(ref_genes6$X307041717)
genes_bin6$density <- (genes_bin6$freq/(genes_bin6$length/1000000))/1000
plot(genes_bin6$foo.X1/1000000, genes_bin6$density,type="l")
saveRDS(genes_bin6$foo.X1/1000000,file="chr6_locipos")
saveRDS(genes_bin6$density,file="chr6_geneProb")

genes_bin7 <- as.data.frame(summary(binning(ref_genes7$X1, nbins = round(max(ref_genes7$X307041717)/100000), type = "kmeans")))
genes_bin7 <- within(genes_bin7, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin7$levels), ',', fixed=TRUE))))
genes_bin7 <- do.call(data.frame, genes_bin7)
genes_bin7 <- genes_bin7 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin7 <- genes_bin7 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin7[1,4]<-11647
genes_bin7$length <- (genes_bin7$foo.X2-genes_bin7$foo.X1)
genes_bin7[round(max(ref_genes7$X307041717)/100000),5] <- max(ref_genes7$X307041717)
genes_bin7$density <- (genes_bin7$freq/(genes_bin7$length/1000000))/1000
plot(genes_bin7$foo.X1/1000000, genes_bin7$density,type="l")
saveRDS(genes_bin7$foo.X1/1000000,file="chr7_locipos")
saveRDS(genes_bin7$density,file="chr7_geneProb")

genes_bin8 <- as.data.frame(summary(binning(ref_genes8$X1, nbins = round(max(ref_genes8$X307041717)/100000), type = "kmeans")))
genes_bin8 <- within(genes_bin8, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin8$levels), ',', fixed=TRUE))))
genes_bin8 <- do.call(data.frame, genes_bin8)
genes_bin8 <- genes_bin8 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin8 <- genes_bin8 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin8[1,4]<-17350
genes_bin8$length <- (genes_bin8$foo.X2-genes_bin8$foo.X1)
genes_bin8[round(max(ref_genes8$X307041717)/100000),5] <- max(ref_genes8$X307041717)
genes_bin8$density <- (genes_bin8$freq/(genes_bin8$length/1000000))/1000
plot(genes_bin8$foo.X1/1000000, genes_bin8$density,type="l")
saveRDS(genes_bin8$foo.X1/1000000,file="chr8_locipos")
saveRDS(genes_bin8$density,file="chr8_geneProb")

genes_bin9 <- as.data.frame(summary(binning(ref_genes9$X1, nbins = round(max(ref_genes9$X307041717)/100000), type = "kmeans")))
genes_bin9 <- within(genes_bin9, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin9$levels), ',', fixed=TRUE))))
genes_bin9 <- do.call(data.frame, genes_bin9)
genes_bin9 <- genes_bin9 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin9 <- genes_bin9 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin9[1,4]<-145180
genes_bin9$length <- (genes_bin9$foo.X2-genes_bin9$foo.X1)
genes_bin9[round(max(ref_genes9$X307041717)/100000),5] <- max(ref_genes9$X307041717)
genes_bin9$density <- (genes_bin9$freq/(genes_bin9$length/1000000))/1000
plot(genes_bin9$foo.X1/1000000, genes_bin9$density,type="l")
saveRDS(genes_bin9$foo.X1/1000000,file="chr9_locipos")
saveRDS(genes_bin9$density,file="chr9_geneProb")

genes_bin10 <- binning(ref_genes10$X1, nbins = round(max(ref_genes10$X307041717)/100000), type = "kmeans")
genes_bin10 <- as.data.frame(summary(binning(ref_genes10$X1, nbins = round(max(ref_genes10$X307041717)/100000), type = "kmeans")))
genes_bin10 <- within(genes_bin10, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin10$levels), ',', fixed=TRUE))))
genes_bin10 <- do.call(data.frame, genes_bin10)
genes_bin10 <- genes_bin10 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin10 <- genes_bin10 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin10[1,4]<-44901
genes_bin10$length <- (genes_bin10$foo.X2-genes_bin10$foo.X1)
genes_bin10[round(max(ref_genes10$X307041717)/100000),5] <- max(ref_genes10$X307041717)
genes_bin10$density <- (genes_bin10$freq/(genes_bin10$length/1000000))/1000
plot(genes_bin10$foo.X1/1000000, genes_bin10$density,type="l")
saveRDS(genes_bin10$foo.X1/1000000,file="chr10_locipos")
saveRDS(genes_bin10$density,file="chr10_geneProb")

genes_bin11 <- binning(ref_genes11$X1, nbins = round(max(ref_genes11$X307041717)/100000), type = "kmeans")
genes_bin11 <- as.data.frame(summary(binning(ref_genes11$X1, nbins = round(max(ref_genes11$X307041717)/100000), type = "kmeans")))
genes_bin11 <- within(genes_bin11, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin11$levels), ',', fixed=TRUE))))
genes_bin11 <- do.call(data.frame, genes_bin11)
genes_bin11 <- genes_bin11 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin11 <- genes_bin11 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin11[1,4]<-1179
genes_bin11$length <- (genes_bin11$foo.X2-genes_bin11$foo.X1)
genes_bin11[round(max(ref_genes11$X307041717)/100000),5] <- max(ref_genes11$X307041717)
genes_bin11$density <- (genes_bin11$freq/(genes_bin11$length/1000000))/1000
plot(genes_bin11$foo.X1/1000000, genes_bin11$density,type="l")
saveRDS(genes_bin11$foo.X1/1000000,file="chr11_locipos")
saveRDS(genes_bin11$density,file="chr11_geneProb")

genes_bin12 <- binning(ref_genes12$X1, nbins = round(max(ref_genes12$X307041717)/100000), type = "kmeans")
genes_bin12 <- as.data.frame(summary(binning(ref_genes12$X1, nbins = round(max(ref_genes12$X307041717)/100000), type = "kmeans")))
genes_bin12 <- within(genes_bin12, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin12$levels), ',', fixed=TRUE))))
genes_bin12 <- do.call(data.frame, genes_bin12)
genes_bin12 <- genes_bin12 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin12 <- genes_bin12 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin12[1,4]<-2681
genes_bin12$length <- (genes_bin12$foo.X2-genes_bin12$foo.X1)
genes_bin12[round(max(ref_genes12$X307041717)/100000),5] <- max(ref_genes12$X307041717)
genes_bin12$density <- (genes_bin12$freq/(genes_bin12$length/1000000))/1000
plot(genes_bin12$foo.X1/1000000, genes_bin12$density,type="l")
saveRDS(genes_bin12$foo.X1/1000000,file="chr12_locipos")
saveRDS(genes_bin12$density,file="chr12_geneProb")


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
