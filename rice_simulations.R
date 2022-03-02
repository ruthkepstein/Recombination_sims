library(AlphaSimR)
library(Rcpp)
library(ggplot2)
library(dplyr)

setwd("C:/Users/16192/Documents/PNAS_Simulations")
set.seed(420)

jap_snps <- read.table("japonica_SNPs.bed", header =FALSE)
ind_snps <- read.table("indica_snps.bed", header =FALSE)
colnames(jap_snps) <- c("Chr#", "SNP Start", "SNP End")
colnames(ind_snps) <- c("Chr#", "SNP Start", "SNP End")
#sample SNPs?
jap_snps <- sample_n(jap_snps, 2000)
ind_snps <- sample_n(ind_snps, 2000)
jap_snps <- jap_snps[order(jap_snps$`Chr#`,jap_snps$`SNP Start`),]
ind_snps <- ind_snps[order(ind_snps$`Chr#`,ind_snps$`SNP Start`),]

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

#ind
ind_chr1_snp <- ind_snps[ which(ind_snps$`Chr#` == "Chr1"),]
ind_chr1_snp$rate <- NA
ind_chr1_snp$`SNP End` <- ind_chr1_snp$`SNP End` - min(ind_chr1_snp$`SNP Start`)
ind_chr1_snp$`SNP Start` <- ind_chr1_snp$`SNP Start`- min(ind_chr1_snp$`SNP Start`)

ind_chr2_snp <- ind_snps[ which(ind_snps$`Chr#` == "Chr2"),]
ind_chr2_snp$rate <- NA
ind_chr2_snp$`SNP End` <- ind_chr2_snp$`SNP End` - min(ind_chr2_snp$`SNP Start`)
ind_chr2_snp$`SNP Start` <- ind_chr2_snp$`SNP Start`- min(ind_chr2_snp$`SNP Start`)

ind_chr3_snp <- ind_snps[ which(ind_snps$`Chr#` == "Chr3"),]
ind_chr3_snp$rate <- NA
ind_chr3_snp$`SNP End` <- ind_chr3_snp$`SNP End` - min(ind_chr3_snp$`SNP Start`)
ind_chr3_snp$`SNP Start` <- ind_chr3_snp$`SNP Start`- min(ind_chr3_snp$`SNP Start`)

ind_chr4_snp <- ind_snps[ which(ind_snps$`Chr#` == "Chr4"),]
ind_chr4_snp$rate <- NA
ind_chr4_snp$`SNP End` <- ind_chr4_snp$`SNP End` - min(ind_chr4_snp$`SNP Start`)
ind_chr4_snp$`SNP Start` <- ind_chr4_snp$`SNP Start`- min(ind_chr4_snp$`SNP Start`)

ind_chr5_snp <- ind_snps[ which(ind_snps$`Chr#` == "Chr5"),]
ind_chr5_snp$rate <- NA
ind_chr5_snp$`SNP End` <- ind_chr5_snp$`SNP End` - min(ind_chr5_snp$`SNP Start`)
ind_chr5_snp$`SNP Start` <- ind_chr5_snp$`SNP Start`- min(ind_chr5_snp$`SNP Start`)

ind_chr6_snp <- ind_snps[ which(ind_snps$`Chr#` == "Chr6"),]
ind_chr6_snp$rate <- NA
ind_chr6_snp$`SNP End` <- ind_chr6_snp$`SNP End` - min(ind_chr6_snp$`SNP Start`)
ind_chr6_snp$`SNP Start` <- ind_chr6_snp$`SNP Start`- min(ind_chr6_snp$`SNP Start`)

ind_chr7_snp <- ind_snps[ which(ind_snps$`Chr#` == "Chr7"),]
ind_chr7_snp$rate <- NA
ind_chr7_snp$`SNP End` <- ind_chr7_snp$`SNP End` - min(ind_chr7_snp$`SNP Start`)
ind_chr7_snp$`SNP Start` <- ind_chr7_snp$`SNP Start`- min(ind_chr7_snp$`SNP Start`)

ind_chr8_snp <- ind_snps[ which(ind_snps$`Chr#` == "Chr8"),]
ind_chr8_snp$rate <- NA
ind_chr8_snp$`SNP End` <- ind_chr8_snp$`SNP End` - min(ind_chr8_snp$`SNP Start`)
ind_chr8_snp$`SNP Start` <- ind_chr8_snp$`SNP Start`- min(ind_chr8_snp$`SNP Start`)

ind_chr9_snp <- ind_snps[ which(ind_snps$`Chr#` == "Chr9"),]
ind_chr9_snp$rate <- NA
ind_chr9_snp$`SNP End` <- ind_chr9_snp$`SNP End` - min(ind_chr9_snp$`SNP Start`)
ind_chr9_snp$`SNP Start` <- ind_chr9_snp$`SNP Start`- min(ind_chr9_snp$`SNP Start`)

ind_chr10_snp <- ind_snps[ which(ind_snps$`Chr#` == "Chr10"),]
ind_chr10_snp$rate <- NA
ind_chr10_snp$`SNP End` <- ind_chr10_snp$`SNP End` - min(ind_chr10_snp$`SNP Start`)
ind_chr10_snp$`SNP Start` <- ind_chr10_snp$`SNP Start`- min(ind_chr10_snp$`SNP Start`)

ind_chr11_snp <- ind_snps[ which(ind_snps$`Chr#` == "Chr11"),]
ind_chr11_snp$rate <- NA
ind_chr11_snp$`SNP End` <- ind_chr11_snp$`SNP End` - min(ind_chr11_snp$`SNP Start`)
ind_chr11_snp$`SNP Start` <- ind_chr11_snp$`SNP Start`- min(ind_chr11_snp$`SNP Start`)

ind_chr12_snp <- ind_snps[ which(ind_snps$`Chr#` == "Chr12"),]
ind_chr12_snp$rate <- NA
ind_chr12_snp$`SNP End` <- ind_chr12_snp$`SNP End` - min(ind_chr12_snp$`SNP Start`)
ind_chr12_snp$`SNP Start` <- ind_chr12_snp$`SNP Start`- min(ind_chr12_snp$`SNP Start`)

#crossover data & midpoint
jap_CO <- read.table("japonica_sequenceLDhotspots.bed", header = FALSE)
jap_CO <- jap_CO[,-c(4:5)]
colnames(jap_CO) <- c("Chr", "CO Start", "CO End")
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

ind_CO <- read.table("indica_sequenceLDhotspots.bed", header = FALSE)
ind_CO <- ind_CO[,-c(4:5)]
colnames(ind_CO) <- c("Chr", "CO Start", "CO End")
ind_CO <- ind_CO[order(ind_CO$Chr,ind_CO$`CO Start`),]

ind_chr1_CO <- ind_CO[ which(ind_CO$Chr == "chr01"),]
ind_chr1_CO$midpoint <- (ind_chr1_CO$`CO Start`+ ind_chr1_CO$`CO End`)/2
ind_chr1_CO <- ind_chr1_CO[order(ind_chr1_CO$`CO Start`),]

ind_chr2_CO <- ind_CO[ which(ind_CO$Chr == "chr02"),]
ind_chr2_CO$midpoint <- (ind_chr2_CO$`CO Start`+ ind_chr2_CO$`CO End`)/2
ind_chr2_CO <- ind_chr2_CO[order(ind_chr2_CO$`CO Start`),]

ind_chr3_CO <- ind_CO[ which(ind_CO$Chr == "chr03"),]
ind_chr3_CO$midpoint <- (ind_chr3_CO$`CO Start`+ ind_chr3_CO$`CO End`)/2
ind_chr3_CO <- ind_chr3_CO[order(ind_chr3_CO$`CO Start`),]

ind_chr4_CO <- ind_CO[ which(ind_CO$Chr == "chr04"),]
ind_chr4_CO$midpoint <- (ind_chr4_CO$`CO Start`+ ind_chr4_CO$`CO End`)/2
ind_chr4_CO <- ind_chr4_CO[order(ind_chr4_CO$`CO Start`),]

ind_chr5_CO <- ind_CO[ which(ind_CO$Chr == "chr05"),]
ind_chr5_CO$midpoint <- (ind_chr5_CO$`CO Start`+ ind_chr5_CO$`CO End`)/2
ind_chr5_CO <- ind_chr5_CO[order(ind_chr5_CO$`CO Start`),]

ind_chr6_CO <- ind_CO[ which(ind_CO$Chr == "chr06"),]
ind_chr6_CO$midpoint <- (ind_chr6_CO$`CO Start`+ ind_chr6_CO$`CO End`)/2
ind_chr6_CO <- ind_chr6_CO[order(ind_chr6_CO$`CO Start`),]

ind_chr7_CO <- ind_CO[ which(ind_CO$Chr == "chr07"),]
ind_chr7_CO$midpoint <- (ind_chr7_CO$`CO Start`+ ind_chr7_CO$`CO End`)/2
ind_chr7_CO <- ind_chr7_CO[order(ind_chr7_CO$`CO Start`),]

ind_chr8_CO <- ind_CO[ which(ind_CO$Chr == "chr08"),]
ind_chr8_CO$midpoint <- (ind_chr8_CO$`CO Start`+ ind_chr8_CO$`CO End`)/2
ind_chr8_CO <- ind_chr8_CO[order(ind_chr8_CO$`CO Start`),]

ind_chr9_CO <- ind_CO[ which(ind_CO$Chr == "chr09"),]
ind_chr9_CO$midpoint <- (ind_chr9_CO$`CO Start`+ ind_chr9_CO$`CO End`)/2
ind_chr9_CO <- ind_chr9_CO[order(ind_chr9_CO$`CO Start`),]

ind_chr10_CO <- ind_CO[ which(ind_CO$Chr == "chr10"),]
ind_chr10_CO$midpoint <- (ind_chr10_CO$`CO Start`+ ind_chr10_CO$`CO End`)/2
ind_chr10_CO <- ind_chr10_CO[order(ind_chr10_CO$`CO Start`),]

ind_chr11_CO <- ind_CO[ which(ind_CO$Chr == "chr11"),]
ind_chr11_CO$midpoint <- (ind_chr11_CO$`CO Start`+ ind_chr11_CO$`CO End`)/2
ind_chr11_CO <- ind_chr11_CO[order(ind_chr11_CO$`CO Start`),]

ind_chr12_CO <- ind_CO[ which(ind_CO$Chr == "chr12"),]
ind_chr12_CO$midpoint <- (ind_chr12_CO$`CO Start`+ ind_chr12_CO$`CO End`)/2
ind_chr12_CO <- ind_chr12_CO[order(ind_chr12_CO$`CO Start`),]

#calculating recombination rate per bin of CO data
library(dlookr)
library(tidyverse)
library(OneR)

#bin crossovers japonica data
#recomb. freq. = (# of COs/ size of population *100%)/ length of bin in Mb
jap_chr1_bin <- binning(jap_chr1_CO$midpoint, nbins = 50, type = "kmeans")
jap_chr1_bin <- as.data.frame(summary(jap_chr1_bin))
#transforming data; making bin interval into 2 columns
jap_chr1_bin <- within(jap_chr1_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(jap_chr1_bin$levels), ',', fixed=TRUE))))
jap_chr1_bin <- do.call(data.frame, jap_chr1_bin)
jap_chr1_bin <- jap_chr1_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
jap_chr1_bin <- jap_chr1_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
jap_chr1_bin[1,4] <- 2597
#making intervals start at 0
jap_chr1_bin$foo.X1 <- jap_chr1_bin$foo.X1 - 2597
jap_chr1_bin$foo.X2 <- jap_chr1_bin$foo.X2 - 2597
#expanding last bin to include last SNP site to avoid NAs in future
jap_chr1_bin[50,5] <- max(japonica_chr1_snp$`SNP End`)
#adding length of bin as column and making in Mb
jap_chr1_bin$length <- (jap_chr1_bin$foo.X2-jap_chr1_bin$foo.X1)/1000000
jap_chr1_bin$rate <- ((jap_chr1_bin$freq/75)*100)/jap_chr1_bin$length

jap_chr2_bin <- binning(jap_chr2_CO$midpoint, nbins = 100, type = "kmeans")
jap_chr2_bin <- as.data.frame(summary(jap_chr2_bin))
jap_chr2_bin <- within(jap_chr2_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(jap_chr2_bin$levels), ',', fixed=TRUE))))
jap_chr2_bin <- do.call(data.frame, jap_chr2_bin)
jap_chr2_bin <- jap_chr2_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
jap_chr2_bin <- jap_chr2_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
jap_chr2_bin[1,4] <- 4236059
jap_chr2_bin$foo.X1 <- jap_chr2_bin$foo.X1 - 4236059
jap_chr2_bin$foo.X2 <- jap_chr2_bin$foo.X2 - 4236059
jap_chr2_bin[50,5] <- max(japonica_chr2_snp$`SNP End`)
jap_chr2_bin$length <- (jap_chr2_bin$foo.X2-jap_chr2_bin$foo.X1)/1000000
jap_chr2_bin$rate <- ((jap_chr2_bin$freq/75)*100)/jap_chr2_bin$length

jap_chr3_bin <- binning(jap_chr3_CO$midpoint, nbins = 50, type = "kmeans")
jap_chr3_bin <- as.data.frame(summary(jap_chr3_bin))
jap_chr3_bin <- within(jap_chr3_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(jap_chr3_bin$levels), ',', fixed=TRUE))))
jap_chr3_bin <- do.call(data.frame, jap_chr3_bin)
jap_chr3_bin <- jap_chr3_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
jap_chr3_bin <- jap_chr3_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
jap_chr3_bin[1,4] <- 4303164
jap_chr3_bin$foo.X1 <- jap_chr3_bin$foo.X1 - 4303164
jap_chr3_bin$foo.X2 <- jap_chr3_bin$foo.X2 - 4303164
jap_chr3_bin[50,5] <- max(japonica_chr3_snp$`SNP End`)
jap_chr3_bin$length <- (jap_chr3_bin$foo.X2-jap_chr3_bin$foo.X1)/1000000
jap_chr3_bin$rate <- ((jap_chr3_bin$freq/75)*100)/jap_chr3_bin$length

jap_chr4_bin <- binning(jap_chr4_CO$midpoint, nbins = 50, type = "kmeans")
jap_chr4_bin <- as.data.frame(summary(jap_chr4_bin))
jap_chr4_bin <- within(jap_chr4_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(jap_chr4_bin$levels), ',', fixed=TRUE))))
jap_chr4_bin <- do.call(data.frame, jap_chr4_bin)
jap_chr4_bin <- jap_chr4_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
jap_chr4_bin <- jap_chr4_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
jap_chr4_bin[1,4] <- 7369
jap_chr4_bin$foo.X1 <- jap_chr4_bin$foo.X1 - 7369
jap_chr4_bin$foo.X2 <- jap_chr4_bin$foo.X2 - 7369
jap_chr4_bin[50,5] <- max(japonica_chr4_snp$`SNP End`)
jap_chr4_bin$length <- (jap_chr4_bin$foo.X2-jap_chr4_bin$foo.X1)/1000000
jap_chr4_bin$rate <- ((jap_chr4_bin$freq/75)*100)/jap_chr4_bin$length

jap_chr5_bin <- binning(jap_chr5_CO$midpoint, nbins = 50, type = "kmeans")
jap_chr5_bin <- as.data.frame(summary(jap_chr5_bin))
jap_chr5_bin <- within(jap_chr5_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(jap_chr5_bin$levels), ',', fixed=TRUE))))
jap_chr5_bin <- do.call(data.frame, jap_chr5_bin)
jap_chr5_bin <- jap_chr5_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
jap_chr5_bin <- jap_chr5_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
jap_chr5_bin[1,4] <- 111888
jap_chr5_bin$foo.X1 <- jap_chr5_bin$foo.X1 - 111888
jap_chr5_bin$foo.X2 <- jap_chr5_bin$foo.X2 - 111888
jap_chr5_bin[50,5] <- max(japonica_chr5_snp$`SNP End`)
jap_chr5_bin$length <- (jap_chr5_bin$foo.X2-jap_chr5_bin$foo.X1)/1000000
jap_chr5_bin$rate <- ((jap_chr5_bin$freq/75)*100)/jap_chr5_bin$length

jap_chr6_bin <- binning(jap_chr6_CO$midpoint, nbins = 50, type = "kmeans")
jap_chr6_bin <- as.data.frame(summary(jap_chr6_bin))
jap_chr6_bin <- within(jap_chr6_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(jap_chr6_bin$levels), ',', fixed=TRUE))))
jap_chr6_bin <- do.call(data.frame, jap_chr6_bin)
jap_chr6_bin <- jap_chr6_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
jap_chr6_bin <- jap_chr6_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
jap_chr6_bin[1,4] <- 1092668
jap_chr6_bin$foo.X1 <- jap_chr6_bin$foo.X1 - 1092668
jap_chr6_bin$foo.X2 <- jap_chr6_bin$foo.X2 - 1092668
jap_chr6_bin[50,5] <- max(japonica_chr6_snp$`SNP End`)
jap_chr6_bin$length <- (jap_chr6_bin$foo.X2-jap_chr6_bin$foo.X1)/1000000
jap_chr6_bin$rate <- ((jap_chr6_bin$freq/75)*100)/jap_chr6_bin$length

jap_chr7_bin <- binning(jap_chr7_CO$midpoint, nbins = 50, type = "kmeans")
jap_chr7_bin <- as.data.frame(summary(jap_chr7_bin))
jap_chr7_bin <- within(jap_chr7_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(jap_chr7_bin$levels), ',', fixed=TRUE))))
jap_chr7_bin <- do.call(data.frame, jap_chr7_bin)
jap_chr7_bin <- jap_chr7_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
jap_chr7_bin <- jap_chr7_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
jap_chr7_bin[1,4] <- 11552
jap_chr7_bin$foo.X1 <- jap_chr7_bin$foo.X1 - 11552
jap_chr7_bin$foo.X2 <- jap_chr7_bin$foo.X2 - 11552
jap_chr7_bin[50,5] <- max(japonica_chr7_snp$`SNP End`)
jap_chr7_bin$length <- (jap_chr7_bin$foo.X2-jap_chr7_bin$foo.X1)/1000000
jap_chr7_bin$rate <- ((jap_chr7_bin$freq/75)*100)/jap_chr7_bin$length

jap_chr8_bin <- binning(jap_chr8_CO$midpoint, nbins = 50, type = "kmeans")
jap_chr8_bin <- as.data.frame(summary(jap_chr8_bin))
jap_chr8_bin <- within(jap_chr8_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(jap_chr8_bin$levels), ',', fixed=TRUE))))
jap_chr8_bin <- do.call(data.frame, jap_chr8_bin)
jap_chr8_bin <- jap_chr8_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
jap_chr8_bin <- jap_chr8_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
jap_chr8_bin[1,4] <- 15685
jap_chr8_bin$foo.X1 <- jap_chr8_bin$foo.X1 - 15685
jap_chr8_bin$foo.X2 <- jap_chr8_bin$foo.X2 - 15685
jap_chr8_bin[50,5] <- max(japonica_chr8_snp$`SNP End`)
jap_chr8_bin$length <- (jap_chr8_bin$foo.X2-jap_chr8_bin$foo.X1)/1000000
jap_chr8_bin$rate <- ((jap_chr8_bin$freq/75)*100)/jap_chr8_bin$length

jap_chr9_bin <- binning(jap_chr9_CO$midpoint, nbins = 50, type = "kmeans")
jap_chr9_bin <- as.data.frame(summary(jap_chr9_bin))
jap_chr9_bin <- within(jap_chr9_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(jap_chr9_bin$levels), ',', fixed=TRUE))))
jap_chr9_bin <- do.call(data.frame, jap_chr9_bin)
jap_chr9_bin <- jap_chr9_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
jap_chr9_bin <- jap_chr9_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
jap_chr9_bin[1,4] <- 38596
jap_chr9_bin$foo.X1 <- jap_chr9_bin$foo.X1 - 38596
jap_chr9_bin$foo.X2 <- jap_chr9_bin$foo.X2 - 38596
jap_chr9_bin[50,5] <- max(japonica_chr9_snp$`SNP End`)
jap_chr9_bin$length <- (jap_chr9_bin$foo.X2-jap_chr9_bin$foo.X1)/1000000
jap_chr9_bin$rate <- ((jap_chr9_bin$freq/75)*100)/jap_chr9_bin$length

jap_chr10_bin <- binning(jap_chr10_CO$midpoint, nbins = 50, type = "kmeans")
jap_chr10_bin <- as.data.frame(summary(jap_chr10_bin))
jap_chr10_bin <- within(jap_chr10_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(jap_chr10_bin$levels), ',', fixed=TRUE))))
jap_chr10_bin <- do.call(data.frame, jap_chr10_bin)
jap_chr10_bin <- jap_chr10_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
jap_chr10_bin <- jap_chr10_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
jap_chr10_bin[1,4] <- 49838
jap_chr10_bin$foo.X1 <- jap_chr10_bin$foo.X1 - 49838
jap_chr10_bin$foo.X2 <- jap_chr10_bin$foo.X2 - 49838
jap_chr10_bin[50,5] <- max(japonica_chr10_snp$`SNP End`)
jap_chr10_bin$length <- (jap_chr10_bin$foo.X2-jap_chr10_bin$foo.X1)/1000000
jap_chr10_bin$rate <- ((jap_chr10_bin$freq/75)*100)/jap_chr10_bin$length

jap_chr11_bin <- binning(jap_chr11_CO$midpoint, nbins = 50, type = "kmeans")
jap_chr11_bin <- as.data.frame(summary(jap_chr11_bin))
jap_chr11_bin <- within(jap_chr11_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(jap_chr11_bin$levels), ',', fixed=TRUE))))
jap_chr11_bin <- do.call(data.frame, jap_chr11_bin)
jap_chr11_bin <- jap_chr11_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
jap_chr11_bin <- jap_chr11_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
jap_chr11_bin[1,4] <- 4252839
jap_chr11_bin$foo.X1 <- jap_chr11_bin$foo.X1 - 4252839
jap_chr11_bin$foo.X2 <- jap_chr11_bin$foo.X2 - 4252839
jap_chr11_bin[50,5] <- max(japonica_chr11_snp$`SNP End`)
jap_chr11_bin$length <- (jap_chr11_bin$foo.X2-jap_chr11_bin$foo.X1)/1000000
jap_chr11_bin$rate <- ((jap_chr11_bin$freq/75)*100)/jap_chr11_bin$length

jap_chr12_bin <- binning(jap_chr12_CO$midpoint, nbins = 50, type = "kmeans")
jap_chr12_bin <- as.data.frame(summary(jap_chr12_bin))
jap_chr12_bin <- within(jap_chr12_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(jap_chr12_bin$levels), ',', fixed=TRUE))))
jap_chr12_bin <- do.call(data.frame, jap_chr12_bin)
jap_chr12_bin <- jap_chr12_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
jap_chr12_bin <- jap_chr12_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
jap_chr12_bin[1,4] <- 114745
jap_chr12_bin$foo.X1 <- jap_chr12_bin$foo.X1 - 114745
jap_chr12_bin$foo.X2 <- jap_chr12_bin$foo.X2 - 114745
jap_chr12_bin[50,5] <- max(japonica_chr12_snp$`SNP End`)
jap_chr12_bin$length <- (jap_chr12_bin$foo.X2-jap_chr12_bin$foo.X1)/1000000
jap_chr12_bin$rate <- ((jap_chr12_bin$freq/75)*100)/jap_chr12_bin$length


#bin crossovers indica data
#recomb. freq. = (# of COs/ size of population *100%)/ length of bin in Mb
ind_chr1_bin <- binning(ind_chr1_CO$midpoint, nbins = 50, type = "kmeans")
ind_chr1_bin <- as.data.frame(summary(ind_chr1_bin))
ind_chr1_bin <- within(ind_chr1_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(ind_chr1_bin$levels), ',', fixed=TRUE))))
ind_chr1_bin <- do.call(data.frame, ind_chr1_bin)
ind_chr1_bin <- ind_chr1_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
ind_chr1_bin <- ind_chr1_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
ind_chr1_bin[1,4] <- 84703
ind_chr1_bin$foo.X1 <- ind_chr1_bin$foo.X1 - 84703
ind_chr1_bin$foo.X2 <- ind_chr1_bin$foo.X2 - 84703
ind_chr1_bin[50,5] <- max(indica_chr1_snp$`SNP End`)
ind_chr1_bin$length <- (ind_chr1_bin$foo.X2-ind_chr1_bin$foo.X1)/1000000
ind_chr1_bin$rate <- ((ind_chr1_bin$freq/75)*100)/ind_chr1_bin$length

ind_chr2_bin <- binning(ind_chr2_CO$midpoint, nbins = 50, type = "kmeans")
ind_chr2_bin <- as.data.frame(summary(ind_chr2_bin))
ind_chr2_bin <- within(ind_chr2_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(ind_chr2_bin$levels), ',', fixed=TRUE))))
ind_chr2_bin <- do.call(data.frame, ind_chr2_bin)
ind_chr2_bin <- ind_chr2_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
ind_chr2_bin <- ind_chr2_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
ind_chr2_bin[1,4] <- 24880
ind_chr2_bin$foo.X1 <- ind_chr2_bin$foo.X1 - 24880
ind_chr2_bin$foo.X2 <- ind_chr2_bin$foo.X2 - 24880
ind_chr2_bin[50,5] <- max(indica_chr2_snp$`SNP End`)
ind_chr2_bin$length <- (ind_chr2_bin$foo.X2-ind_chr2_bin$foo.X1)/1000000
ind_chr2_bin$rate <- ((ind_chr2_bin$freq/75)*100)/ind_chr2_bin$length

ind_chr3_bin <- binning(ind_chr3_CO$midpoint, nbins = 50, type = "kmeans")
ind_chr3_bin <- as.data.frame(summary(ind_chr3_bin))
ind_chr3_bin <- within(ind_chr3_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(ind_chr3_bin$levels), ',', fixed=TRUE))))
ind_chr3_bin <- do.call(data.frame, ind_chr3_bin)
ind_chr3_bin <- ind_chr3_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
ind_chr3_bin <- ind_chr3_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
ind_chr3_bin[1,4] <- 156144
ind_chr3_bin$foo.X1 <- ind_chr3_bin$foo.X1 - 156144
ind_chr3_bin$foo.X2 <- ind_chr3_bin$foo.X2 - 156144
ind_chr3_bin[50,5] <- max(indica_chr3_snp$`SNP End`)
ind_chr3_bin$length <- (ind_chr3_bin$foo.X2-ind_chr3_bin$foo.X1)/1000000
ind_chr3_bin$rate <- ((ind_chr3_bin$freq/75)*100)/ind_chr3_bin$length

ind_chr4_bin <- binning(ind_chr4_CO$midpoint, nbins = 50, type = "kmeans")
ind_chr4_bin <- as.data.frame(summary(ind_chr4_bin))
ind_chr4_bin <- within(ind_chr4_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(ind_chr4_bin$levels), ',', fixed=TRUE))))
ind_chr4_bin <- do.call(data.frame, ind_chr4_bin)
ind_chr4_bin <- ind_chr4_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
ind_chr4_bin <- ind_chr4_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
ind_chr4_bin[1,4] <- 58437
ind_chr4_bin$foo.X1 <- ind_chr4_bin$foo.X1 - 58437
ind_chr4_bin$foo.X2 <- ind_chr4_bin$foo.X2 - 58437
ind_chr4_bin[50,5] <- max(indica_chr4_snp$`SNP End`)
ind_chr4_bin$length <- (ind_chr4_bin$foo.X2-ind_chr4_bin$foo.X1)/1000000
ind_chr4_bin$rate <- ((ind_chr4_bin$freq/75)*100)/ind_chr4_bin$length

ind_chr5_bin <- binning(ind_chr5_CO$midpoint, nbins = 50, type = "kmeans")
ind_chr5_bin <- as.data.frame(summary(ind_chr5_bin))
ind_chr5_bin <- within(ind_chr5_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(ind_chr5_bin$levels), ',', fixed=TRUE))))
ind_chr5_bin <- do.call(data.frame, ind_chr5_bin)
ind_chr5_bin <- ind_chr5_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
ind_chr5_bin <- ind_chr5_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
ind_chr5_bin[1,4] <- 49396
ind_chr5_bin$foo.X1 <- ind_chr5_bin$foo.X1 - 49396
ind_chr5_bin$foo.X2 <- ind_chr5_bin$foo.X2 - 49396
ind_chr5_bin[50,5] <- max(indica_chr5_snp$`SNP End`)
ind_chr5_bin$length <- (ind_chr5_bin$foo.X2-ind_chr5_bin$foo.X1)/1000000
ind_chr5_bin$rate <- ((ind_chr5_bin$freq/75)*100)/ind_chr5_bin$length

ind_chr6_bin <- binning(ind_chr6_CO$midpoint, nbins = 50, type = "kmeans")
ind_chr6_bin <- as.data.frame(summary(ind_chr6_bin))
ind_chr6_bin <- within(ind_chr6_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(ind_chr6_bin$levels), ',', fixed=TRUE))))
ind_chr6_bin <- do.call(data.frame, ind_chr6_bin)
ind_chr6_bin <- ind_chr6_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
ind_chr6_bin <- ind_chr6_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
ind_chr6_bin[1,4] <- 126510
ind_chr6_bin$foo.X1 <- ind_chr6_bin$foo.X1 - 126510
ind_chr6_bin$foo.X2 <- ind_chr6_bin$foo.X2 - 126510
ind_chr6_bin[50,5] <- max(indica_chr6_snp$`SNP End`)
ind_chr6_bin$length <- (ind_chr6_bin$foo.X2-ind_chr6_bin$foo.X1)/1000000
ind_chr6_bin$rate <- ((ind_chr6_bin$freq/75)*100)/ind_chr6_bin$length

ind_chr7_bin <- binning(ind_chr7_CO$midpoint, nbins = 50, type = "kmeans")
ind_chr7_bin <- as.data.frame(summary(ind_chr7_bin))
ind_chr7_bin <- within(ind_chr7_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(ind_chr7_bin$levels), ',', fixed=TRUE))))
ind_chr7_bin <- do.call(data.frame, ind_chr7_bin)
ind_chr7_bin <- ind_chr7_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
ind_chr7_bin <- ind_chr7_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
ind_chr7_bin[1,4] <- 11552
ind_chr7_bin$foo.X1 <- ind_chr7_bin$foo.X1 - 11552
ind_chr7_bin$foo.X2 <- ind_chr7_bin$foo.X2 - 11552
ind_chr7_bin[50,5] <- max(indica_chr7_snp$`SNP End`)
ind_chr7_bin$length <- (ind_chr7_bin$foo.X2-ind_chr7_bin$foo.X1)/1000000
ind_chr7_bin$rate <- ((ind_chr7_bin$freq/75)*100)/ind_chr7_bin$length

ind_chr8_bin <- binning(ind_chr8_CO$midpoint, nbins = 50, type = "kmeans")
ind_chr8_bin <- as.data.frame(summary(ind_chr8_bin))
ind_chr8_bin <- within(ind_chr8_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(ind_chr8_bin$levels), ',', fixed=TRUE))))
ind_chr8_bin <- do.call(data.frame, ind_chr8_bin)
ind_chr8_bin <- ind_chr8_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
ind_chr8_bin <- ind_chr8_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
ind_chr8_bin[1,4] <- 21142
ind_chr8_bin$foo.X1 <- ind_chr8_bin$foo.X1 - 21142
ind_chr8_bin$foo.X2 <- ind_chr8_bin$foo.X2 - 21142
ind_chr8_bin[50,5] <- max(indica_chr8_snp$`SNP End`)
ind_chr8_bin$length <- (ind_chr8_bin$foo.X2-ind_chr8_bin$foo.X1)/1000000
ind_chr8_bin$rate <- ((ind_chr8_bin$freq/75)*100)/ind_chr8_bin$length

ind_chr9_bin <- binning(ind_chr9_CO$midpoint, nbins = 50, type = "kmeans")
ind_chr9_bin <- as.data.frame(summary(ind_chr9_bin))
ind_chr9_bin <- within(ind_chr9_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(ind_chr9_bin$levels), ',', fixed=TRUE))))
ind_chr9_bin <- do.call(data.frame, ind_chr9_bin)
ind_chr9_bin <- ind_chr9_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
ind_chr9_bin <- ind_chr9_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
ind_chr9_bin[1,4] <- 153162
ind_chr9_bin$foo.X1 <- ind_chr9_bin$foo.X1 - 153162
ind_chr9_bin$foo.X2 <- ind_chr9_bin$foo.X2 - 153162
ind_chr9_bin[50,5] <- max(indica_chr9_snp$`SNP End`)
ind_chr9_bin$length <- (ind_chr9_bin$foo.X2-ind_chr9_bin$foo.X1)/1000000
ind_chr9_bin$rate <- ((ind_chr9_bin$freq/75)*100)/ind_chr9_bin$length

ind_chr10_bin <- binning(ind_chr10_CO$midpoint, nbins = 50, type = "kmeans")
ind_chr10_bin <- as.data.frame(summary(ind_chr10_bin))
ind_chr10_bin <- within(ind_chr10_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(ind_chr10_bin$levels), ',', fixed=TRUE))))
ind_chr10_bin <- do.call(data.frame, ind_chr10_bin)
ind_chr10_bin <- ind_chr10_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
ind_chr10_bin <- ind_chr10_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
ind_chr10_bin[1,4] <- 80838
ind_chr10_bin$foo.X1 <- ind_chr10_bin$foo.X1 - 80838
ind_chr10_bin$foo.X2 <- ind_chr10_bin$foo.X2 - 80838
ind_chr10_bin[50,5] <- max(indica_chr10_snp$`SNP End`)
ind_chr10_bin$length <- (ind_chr10_bin$foo.X2-ind_chr10_bin$foo.X1)/1000000
ind_chr10_bin$rate <- ((ind_chr10_bin$freq/75)*100)/ind_chr10_bin$length

ind_chr11_bin <- binning(ind_chr11_CO$midpoint, nbins = 50, type = "kmeans")
ind_chr11_bin <- as.data.frame(summary(ind_chr11_bin))
ind_chr11_bin <- within(ind_chr11_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(ind_chr11_bin$levels), ',', fixed=TRUE))))
ind_chr11_bin <- do.call(data.frame, ind_chr11_bin)
ind_chr11_bin <- ind_chr11_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
ind_chr11_bin <- ind_chr11_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
ind_chr11_bin[1,4] <- 4252839
ind_chr11_bin$foo.X1 <- ind_chr11_bin$foo.X1 - 4252839
ind_chr11_bin$foo.X2 <- ind_chr11_bin$foo.X2 - 4252839
ind_chr11_bin[50,5] <- max(indica_chr11_snp$`SNP End`)
ind_chr11_bin$length <- (ind_chr11_bin$foo.X2-ind_chr11_bin$foo.X1)/1000000
ind_chr11_bin$rate <- ((ind_chr11_bin$freq/75)*100)/ind_chr11_bin$length

ind_chr12_bin <- binning(ind_chr12_CO$midpoint, nbins = 50, type = "kmeans")
ind_chr12_bin <- as.data.frame(summary(ind_chr12_bin))
ind_chr12_bin <- within(ind_chr12_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(ind_chr12_bin$levels), ',', fixed=TRUE))))
ind_chr12_bin <- do.call(data.frame, ind_chr12_bin)
ind_chr12_bin <- ind_chr12_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
ind_chr12_bin <- ind_chr12_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
ind_chr12_bin[1,4] <- 44123
ind_chr12_bin$foo.X1 <- ind_chr12_bin$foo.X1 - 44123
ind_chr12_bin$foo.X2 <- ind_chr12_bin$foo.X2 - 44123
ind_chr12_bin[50,5] <- max(indica_chr12_snp$`SNP End`)
ind_chr12_bin$length <- (ind_chr12_bin$foo.X2-ind_chr12_bin$foo.X1)/1000000
ind_chr12_bin$rate <- ((ind_chr12_bin$freq/75)*100)/ind_chr12_bin$length

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

##JAPONICA
#using function, converted SNP start to Mb to get cM/Mb for final genetic position
jap_chr1_snp2 <- snp_rate(jap_chr1_bin, jap_chr1_snp)
jap_chr1_snp2$`SNP Start`<- jap_chr1_snp2$`SNP Start`/1000000
jap_chr1_snp2 <- jap_chr1_snp2[order(jap_chr1_snp2$`SNP Start`),]
#smoothing the recombination rate so transitions between bins are not so abrupt
jap_chr1_spl <- smooth.spline(jap_chr1_snp2$rate, spar = .9)
#creation of genetic positions from smoothed recombination rate
jap_chr1_snp2$pos <- (jap_chr1_snp2$`SNP Start`*jap_chr1_spl$y)
#graph to look at Mb vs. cM along chromosome
plot(jap_chr1_snp2$`SNP Start`, jap_chr1_snp2$pos)
ggplot(jap_chr1_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
#graph to look at Mb vs. cM/Mb to see recombination rate along chromosome
plot(jap_chr1_snp2$`SNP Start`, jap_chr1_snp2$pos/jap_chr1_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 1 Recombination Distribution")
jap_chr1_finalpos <- jap_chr1_snp2[order(jap_chr1_snp2$pos),]
#want False to input into AlphaSimR
is.unsorted(jap_chr1_finalpos$pos)
#plot again to make sure it looks the same
plot(jap_chr1_snp2$`SNP Start`, jap_chr1_finalpos$pos/jap_chr1_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 1 Recombination Distribution")
plot(jap_chr1_finalpos$`SNP Start`, jap_chr1_finalpos$pos)

jap_chr2_snp2 <- snp_rate(jap_chr2_bin, jap_chr2_snp)
jap_chr2_snp2$`SNP Start` <- jap_chr2_snp2$`SNP Start`/1000000
jap_chr2_spl <- smooth.spline(jap_chr2_snp2$rate, spar = 1.2)
jap_chr2_snp2$pos <- (jap_chr2_snp2$`SNP Start`*jap_chr2_spl$y)
plot(jap_chr2_snp2$`SNP Start`, jap_chr2_snp2$pos)
plot(jap_chr2_snp2$`SNP Start`, jap_chr2_snp2$pos/jap_chr2_snp2$`SNP Start`, type = "l")
jap_chr2_finalpos <- jap_chr2_snp2[order(jap_chr2_snp2$pos),]
is.unsorted(jap_chr2_finalpos$pos)
plot(jap_chr2_snp2$`SNP Start`, jap_chr2_finalpos$pos/jap_chr2_snp2$`SNP Start`, type = "l")

jap_chr3_snp2 <- snp_rate(jap_chr3_bin, jap_chr3_snp)
jap_chr3_snp2$`SNP Start` <- jap_chr3_snp2$`SNP Start`/1000000
jap_chr3_spl <- smooth.spline(jap_chr3_snp2$rate, spar = 1.1)
jap_chr3_snp2$pos <- (jap_chr3_snp2$`SNP Start`*jap_chr3_spl$y)
plot(jap_chr3_snp2$`SNP Start`, jap_chr3_snp2$pos)
plot(jap_chr3_snp2$`SNP Start`, jap_chr3_snp2$pos/jap_chr3_snp2$`SNP Start`, type = "l")
jap_chr3_finalpos <- jap_chr3_snp2[order(jap_chr3_snp2$pos),]
is.unsorted(jap_chr3_finalpos$pos)
plot(jap_chr3_snp2$`SNP Start`, jap_chr3_finalpos$pos/jap_chr3_snp2$`SNP Start`, type = "l")

jap_chr4_snp2 <- snp_rate(jap_chr4_bin, jap_chr4_snp)
jap_chr4_snp2$`SNP Start` <- jap_chr4_snp2$`SNP Start`/1000000
jap_chr4_spl <- smooth.spline(jap_chr4_snp2$rate, spar = 1.15)
jap_chr4_snp2$pos <- (jap_chr4_snp2$`SNP Start`*jap_chr4_spl$y)
plot(jap_chr4_snp2$`SNP Start`, jap_chr4_snp2$pos)
plot(jap_chr4_snp2$`SNP Start`, jap_chr4_snp2$pos/jap_chr4_snp2$`SNP Start`, type = "l")
jap_chr4_finalpos <- jap_chr4_snp2[order(jap_chr4_snp2$pos),]
is.unsorted(jap_chr4_finalpos$pos)
plot(jap_chr4_snp2$`SNP Start`, jap_chr4_finalpos$pos/jap_chr4_snp2$`SNP Start`, type = "l")

jap_chr5_snp2 <- snp_rate(jap_chr5_bin, jap_chr5_snp)
jap_chr5_snp2$`SNP Start` <- jap_chr5_snp2$`SNP Start`/1000000
jap_chr5_spl <- smooth.spline(jap_chr5_snp2$rate, spar = 1.1)
jap_chr5_snp2$pos <- (jap_chr5_snp2$`SNP Start`*jap_chr5_spl$y)
plot(jap_chr5_snp2$`SNP Start`, jap_chr5_snp2$pos)
plot(jap_chr5_snp2$`SNP Start`, jap_chr5_snp2$pos/jap_chr5_snp2$`SNP Start`, type = "l")
jap_chr5_finalpos <- jap_chr5_snp2[order(jap_chr5_snp2$pos),]
is.unsorted(jap_chr5_finalpos$pos)
plot(jap_chr5_snp2$`SNP Start`, jap_chr5_finalpos$pos/jap_chr5_snp2$`SNP Start`, type = "l")

jap_chr6_snp2 <- snp_rate(jap_chr6_bin, jap_chr6_snp)
jap_chr6_snp2$`SNP Start` <- jap_chr6_snp2$`SNP Start`/1000000
jap_chr6_spl <- smooth.spline(jap_chr6_snp2$rate, spar = 1)
jap_chr6_snp2$pos <- (jap_chr6_snp2$`SNP Start`*jap_chr6_spl$y)
plot(jap_chr6_snp2$`SNP Start`, jap_chr6_snp2$pos)
plot(jap_chr6_snp2$`SNP Start`, jap_chr6_snp2$pos/jap_chr6_snp2$`SNP Start`, type = "l")
jap_chr6_finalpos <- jap_chr6_snp2[order(jap_chr6_snp2$pos),]
is.unsorted(jap_chr6_finalpos$pos)
plot(jap_chr6_snp2$`SNP Start`, jap_chr6_finalpos$pos/jap_chr6_snp2$`SNP Start`, type = "l")

jap_chr7_snp2 <- snp_rate(jap_chr7_bin, jap_chr7_snp)
jap_chr7_snp2$`SNP Start` <- jap_chr7_snp2$`SNP Start`/1000000
jap_chr7_spl <- smooth.spline(jap_chr7_snp2$rate, spar = 1.15)
jap_chr7_snp2$pos <- (jap_chr7_snp2$`SNP Start`*jap_chr7_spl$y)
plot(jap_chr7_snp2$`SNP Start`, jap_chr7_snp2$pos)
plot(jap_chr7_snp2$`SNP Start`, jap_chr7_snp2$pos/jap_chr7_snp2$`SNP Start`, type = "l")
jap_chr7_finalpos <- jap_chr7_snp2[order(jap_chr7_snp2$pos),]
is.unsorted(jap_chr7_finalpos$pos)
plot(jap_chr7_snp2$`SNP Start`, jap_chr7_finalpos$pos/jap_chr7_snp2$`SNP Start`, type = "l")

jap_chr8_snp2 <- snp_rate(jap_chr8_bin, jap_chr8_snp)
jap_chr8_snp2$`SNP Start` <- jap_chr8_snp2$`SNP Start`/1000000
jap_chr8_spl <- smooth.spline(jap_chr8_snp2$rate, spar = 1.1)
jap_chr8_snp2$pos <- (jap_chr8_snp2$`SNP Start`*jap_chr8_spl$y)
plot(jap_chr8_snp2$`SNP Start`, jap_chr8_snp2$pos)
plot(jap_chr8_snp2$`SNP Start`, jap_chr8_snp2$pos/jap_chr8_snp2$`SNP Start`, type = "l")
jap_chr8_finalpos <- jap_chr8_snp2[order(jap_chr8_snp2$pos),]
is.unsorted(jap_chr8_finalpos$pos)
plot(jap_chr8_snp2$`SNP Start`, jap_chr8_finalpos$pos/jap_chr8_snp2$`SNP Start`, type = "l")

jap_chr9_snp2 <- snp_rate(jap_chr9_bin, jap_chr9_snp)
jap_chr9_snp2$`SNP Start` <- jap_chr9_snp2$`SNP Start`/1000000
jap_chr9_spl <- smooth.spline(jap_chr9_snp2$rate, spar = 1.1)
jap_chr9_snp2$pos <- (jap_chr9_snp2$`SNP Start`*jap_chr9_spl$y)
plot(jap_chr9_snp2$`SNP Start`, jap_chr9_snp2$pos)
plot(jap_chr9_snp2$`SNP Start`, jap_chr9_snp2$pos/jap_chr9_snp2$`SNP Start`, type = "l")
jap_chr9_finalpos <- jap_chr9_snp2[order(jap_chr9_snp2$pos),]
is.unsorted(jap_chr9_finalpos$pos)
plot(jap_chr9_snp2$`SNP Start`, jap_chr9_finalpos$pos/jap_chr9_snp2$`SNP Start`, type = "l")

jap_chr10_snp2 <- snp_rate(jap_chr10_bin, jap_chr10_snp)
jap_chr10_snp2$`SNP Start` <- jap_chr10_snp2$`SNP Start`/1000000
jap_chr10_spl <- smooth.spline(jap_chr10_snp2$rate, spar = 1.2)
jap_chr10_snp2$pos <- (jap_chr10_snp2$`SNP Start`*jap_chr10_spl$y)
plot(jap_chr10_snp2$`SNP Start`, jap_chr10_snp2$pos)
plot(jap_chr10_snp2$`SNP Start`, jap_chr10_snp2$pos/jap_chr10_snp2$`SNP Start`, type = "l")
jap_chr10_finalpos <- jap_chr10_snp2[order(jap_chr10_snp2$pos),]
is.unsorted(jap_chr10_finalpos$pos)
plot(jap_chr10_snp2$`SNP Start`, jap_chr10_finalpos$pos/jap_chr10_snp2$`SNP Start`, type = "l")

jap_chr11_snp2 <- snp_rate(jap_chr11_bin, jap_chr11_snp)
jap_chr11_snp2$`SNP Start` <- jap_chr11_snp2$`SNP Start`/1000000
jap_chr11_spl <- smooth.spline(jap_chr11_snp2$rate, spar = 1.2)
jap_chr11_snp2$pos <- (jap_chr11_snp2$`SNP Start`*jap_chr11_spl$y)
plot(jap_chr11_snp2$`SNP Start`, jap_chr11_snp2$pos)
plot(jap_chr11_snp2$`SNP Start`, jap_chr11_snp2$pos/jap_chr11_snp2$`SNP Start`, type = "l")
jap_chr11_finalpos <- jap_chr11_snp2[order(jap_chr11_snp2$pos),]
is.unsorted(jap_chr11_finalpos$pos)
plot(jap_chr11_snp2$`SNP Start`, jap_chr11_finalpos$pos/jap_chr11_snp2$`SNP Start`, type = "l")

jap_chr12_snp2 <- snp_rate(jap_chr12_bin, jap_chr12_snp)
jap_chr12_snp2$`SNP Start` <- jap_chr12_snp2$`SNP Start`/1000000
jap_chr12_spl <- smooth.spline(jap_chr12_snp2$rate, spar = 1.2)
jap_chr12_snp2$pos <- (jap_chr12_snp2$`SNP Start`*jap_chr12_spl$y)
plot(jap_chr12_snp2$`SNP Start`, jap_chr12_snp2$pos)
plot(jap_chr12_snp2$`SNP Start`, jap_chr12_snp2$pos/jap_chr12_snp2$`SNP Start`, type = "l")
jap_chr12_finalpos <- jap_chr12_snp2[order(jap_chr12_snp2$pos),]
is.unsorted(jap_chr12_finalpos$pos)
plot(jap_chr12_snp2$`SNP Start`, jap_chr12_finalpos$pos/jap_chr12_snp2$`SNP Start`, type = "l")

##INDICA 
ind_chr1_snp2 <- snp_rate(ind_chr1_bin, ind_chr1_snp)
ind_chr1_snp2$`SNP Start`<- ind_chr1_snp2$`SNP Start`/1000000
ind_chr1_snp2 <- ind_chr1_snp2[order(ind_chr1_snp2$`SNP Start`),]
ind_chr1_spl <- smooth.spline(ind_chr1_snp2$rate, spar = 1)
ind_chr1_snp2$pos <- (ind_chr1_snp2$`SNP Start`*ind_chr1_spl$y)
plot(ind_chr1_snp2$`SNP Start`, ind_chr1_snp2$pos)
ggplot(ind_chr1_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
plot(ind_chr1_snp2$`SNP Start`, ind_chr1_snp2$pos/ind_chr1_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 1 Recombination Distribution")
ind_chr1_finalpos <- ind_chr1_snp2[order(ind_chr1_snp2$pos),]
is.unsorted(ind_chr1_finalpos$pos)
plot(ind_chr1_snp2$`SNP Start`, ind_chr1_finalpos$pos/ind_chr1_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 1 Recombination Distribution")
plot(ind_chr1_finalpos$`SNP Start`, ind_chr1_finalpos$pos)

ind_chr2_snp2 <- snp_rate(ind_chr2_bin, ind_chr2_snp)
ind_chr2_snp2$`SNP Start` <- ind_chr2_snp2$`SNP Start`/1000000
ind_chr2_spl <- smooth.spline(ind_chr2_snp2$rate, spar = 1.2)
ind_chr2_snp2$pos <- (ind_chr2_snp2$`SNP Start`*ind_chr2_spl$y)
plot(ind_chr2_snp2$`SNP Start`, ind_chr2_snp2$pos)
plot(ind_chr2_snp2$`SNP Start`, ind_chr2_snp2$pos/ind_chr2_snp2$`SNP Start`, type = "l")
ind_chr2_finalpos <- ind_chr2_snp2[order(ind_chr2_snp2$pos),]
is.unsorted(ind_chr2_finalpos$pos)
plot(ind_chr2_snp2$`SNP Start`, ind_chr2_finalpos$pos/ind_chr2_snp2$`SNP Start`, type = "l")

ind_chr3_snp2 <- snp_rate(ind_chr3_bin, ind_chr3_snp)
ind_chr3_snp2$`SNP Start` <- ind_chr3_snp2$`SNP Start`/1000000
ind_chr3_spl <- smooth.spline(ind_chr3_snp2$rate, spar = 1.1)
ind_chr3_snp2$pos <- (ind_chr3_snp2$`SNP Start`*ind_chr3_spl$y)
plot(ind_chr3_snp2$`SNP Start`, ind_chr3_snp2$pos)
plot(ind_chr3_snp2$`SNP Start`, ind_chr3_snp2$pos/ind_chr3_snp2$`SNP Start`, type = "l")
ind_chr3_finalpos <- ind_chr3_snp2[order(ind_chr3_snp2$pos),]
is.unsorted(ind_chr3_finalpos$pos)
plot(ind_chr3_snp2$`SNP Start`, ind_chr3_finalpos$pos/ind_chr3_snp2$`SNP Start`, type = "l")

ind_chr4_snp2 <- snp_rate(ind_chr4_bin, ind_chr4_snp)
ind_chr4_snp2$`SNP Start` <- ind_chr4_snp2$`SNP Start`/1000000
ind_chr4_spl <- smooth.spline(ind_chr4_snp2$rate, spar = 1.15)
ind_chr4_snp2$pos <- (ind_chr4_snp2$`SNP Start`*ind_chr4_spl$y)
plot(ind_chr4_snp2$`SNP Start`, ind_chr4_snp2$pos)
plot(ind_chr4_snp2$`SNP Start`, ind_chr4_snp2$pos/ind_chr4_snp2$`SNP Start`, type = "l")
ind_chr4_finalpos <- ind_chr4_snp2[order(ind_chr4_snp2$pos),]
is.unsorted(ind_chr4_finalpos$pos)
plot(ind_chr4_snp2$`SNP Start`, ind_chr4_finalpos$pos/ind_chr4_snp2$`SNP Start`, type = "l")

ind_chr5_snp2 <- snp_rate(ind_chr5_bin, ind_chr5_snp)
ind_chr5_snp2$`SNP Start` <- ind_chr5_snp2$`SNP Start`/1000000
ind_chr5_spl <- smooth.spline(ind_chr5_snp2$rate, spar = 1.1)
ind_chr5_snp2$pos <- (ind_chr5_snp2$`SNP Start`*ind_chr5_spl$y)
plot(ind_chr5_snp2$`SNP Start`, ind_chr5_snp2$pos)
plot(ind_chr5_snp2$`SNP Start`, ind_chr5_snp2$pos/ind_chr5_snp2$`SNP Start`, type = "l")
ind_chr5_finalpos <- ind_chr5_snp2[order(ind_chr5_snp2$pos),]
is.unsorted(ind_chr5_finalpos$pos)
plot(ind_chr5_snp2$`SNP Start`, ind_chr5_finalpos$pos/ind_chr5_snp2$`SNP Start`, type = "l")

ind_chr6_snp2 <- snp_rate(ind_chr6_bin, ind_chr6_snp)
ind_chr6_snp2$`SNP Start` <- ind_chr6_snp2$`SNP Start`/1000000
ind_chr6_spl <- smooth.spline(ind_chr6_snp2$rate, spar = 1)
ind_chr6_snp2$pos <- (ind_chr6_snp2$`SNP Start`*ind_chr6_spl$y)
plot(ind_chr6_snp2$`SNP Start`, ind_chr6_snp2$pos)
plot(ind_chr6_snp2$`SNP Start`, ind_chr6_snp2$pos/ind_chr6_snp2$`SNP Start`, type = "l")
ind_chr6_finalpos <- ind_chr6_snp2[order(ind_chr6_snp2$pos),]
is.unsorted(ind_chr6_finalpos$pos)
plot(ind_chr6_snp2$`SNP Start`, ind_chr6_finalpos$pos/ind_chr6_snp2$`SNP Start`, type = "l")

ind_chr7_snp2 <- snp_rate(ind_chr7_bin, ind_chr7_snp)
ind_chr7_snp2$`SNP Start` <- ind_chr7_snp2$`SNP Start`/1000000
ind_chr7_spl <- smooth.spline(ind_chr7_snp2$rate, spar = 1.15)
ind_chr7_snp2$pos <- (ind_chr7_snp2$`SNP Start`*ind_chr7_spl$y)
plot(ind_chr7_snp2$`SNP Start`, ind_chr7_snp2$pos)
plot(ind_chr7_snp2$`SNP Start`, ind_chr7_snp2$pos/ind_chr7_snp2$`SNP Start`, type = "l")
ind_chr7_finalpos <- ind_chr7_snp2[order(ind_chr7_snp2$pos),]
is.unsorted(ind_chr7_finalpos$pos)
plot(ind_chr7_snp2$`SNP Start`, ind_chr7_finalpos$pos/ind_chr7_snp2$`SNP Start`, type = "l")

ind_chr8_snp2 <- snp_rate(ind_chr8_bin, ind_chr8_snp)
ind_chr8_snp2$`SNP Start` <- ind_chr8_snp2$`SNP Start`/1000000
ind_chr8_spl <- smooth.spline(ind_chr8_snp2$rate, spar = 1.1)
ind_chr8_snp2$pos <- (ind_chr8_snp2$`SNP Start`*ind_chr8_spl$y)
plot(ind_chr8_snp2$`SNP Start`, ind_chr8_snp2$pos)
plot(ind_chr8_snp2$`SNP Start`, ind_chr8_snp2$pos/ind_chr8_snp2$`SNP Start`, type = "l")
ind_chr8_finalpos <- ind_chr8_snp2[order(ind_chr8_snp2$pos),]
is.unsorted(ind_chr8_finalpos$pos)
plot(ind_chr8_snp2$`SNP Start`, ind_chr8_finalpos$pos/ind_chr8_snp2$`SNP Start`, type = "l")

ind_chr9_snp2 <- snp_rate(ind_chr9_bin, ind_chr9_snp)
ind_chr9_snp2$`SNP Start` <- ind_chr9_snp2$`SNP Start`/1000000
ind_chr9_spl <- smooth.spline(ind_chr9_snp2$rate, spar = 1.1)
ind_chr9_snp2$pos <- (ind_chr9_snp2$`SNP Start`*ind_chr9_spl$y)
plot(ind_chr9_snp2$`SNP Start`, ind_chr9_snp2$pos)
plot(ind_chr9_snp2$`SNP Start`, ind_chr9_snp2$pos/ind_chr9_snp2$`SNP Start`, type = "l")
ind_chr9_finalpos <- ind_chr9_snp2[order(ind_chr9_snp2$pos),]
is.unsorted(ind_chr9_finalpos$pos)
plot(ind_chr9_snp2$`SNP Start`, ind_chr9_finalpos$pos/ind_chr9_snp2$`SNP Start`, type = "l")

ind_chr10_snp2 <- snp_rate(ind_chr10_bin, ind_chr10_snp)
ind_chr10_snp2$`SNP Start` <- ind_chr10_snp2$`SNP Start`/1000000
ind_chr10_spl <- smooth.spline(ind_chr10_snp2$rate, spar = 1.2)
ind_chr10_snp2$pos <- (ind_chr10_snp2$`SNP Start`*ind_chr10_spl$y)
plot(ind_chr10_snp2$`SNP Start`, ind_chr10_snp2$pos)
plot(ind_chr10_snp2$`SNP Start`, ind_chr10_snp2$pos/ind_chr10_snp2$`SNP Start`, type = "l")
chr10_finalpos <- ind_chr10_snp2[order(ind_chr10_snp2$pos),]
is.unsorted(ind_chr10_finalpos$pos)
plot(ind_chr10_snp2$`SNP Start`, ind_chr10_finalpos$pos/ind_chr10_snp2$`SNP Start`, type = "l")

ind_chr11_snp2 <- snp_rate(ind_chr11_bin, ind_chr11_snp)
ind_chr11_snp2$`SNP Start` <- ind_chr11_snp2$`SNP Start`/1000000
ind_chr11_spl <- smooth.spline(ind_chr11_snp2$rate, spar = 1.2)
ind_chr11_snp2$pos <- (ind_chr11_snp2$`SNP Start`*ind_chr11_spl$y)
plot(ind_chr11_snp2$`SNP Start`, ind_chr11_snp2$pos)
plot(ind_chr11_snp2$`SNP Start`, ind_chr11_snp2$pos/ind_chr11_snp2$`SNP Start`, type = "l")
ind_chr11_finalpos <- ind_chr11_snp2[order(ind_chr11_snp2$pos),]
is.unsorted(ind_chr11_finalpos$pos)
plot(ind_chr11_snp2$`SNP Start`, ind_chr11_finalpos$pos/ind_chr11_snp2$`SNP Start`, type = "l")

ind_chr12_snp2 <- snp_rate(ind_chr12_bin, ind_chr12_snp)
ind_chr12_snp2$`SNP Start` <- ind_chr12_snp2$`SNP Start`/1000000
ind_chr12_spl <- smooth.spline(ind_chr12_snp2$rate, spar = 1.2)
ind_chr12_snp2$pos <- (ind_chr12_snp2$`SNP Start`*ind_chr12_spl$y)
plot(ind_chr12_snp2$`SNP Start`, ind_chr12_snp2$pos)
plot(ind_chr12_snp2$`SNP Start`, ind_chr12_snp2$pos/ind_chr12_snp2$`SNP Start`, type = "l")
ind_chr12_finalpos <- ind_chr12_snp2[order(ind_chr12_snp2$pos),]
is.unsorted(ind_chr12_finalpos$pos)
plot(ind_chr12_snp2$`SNP Start`, ind_chr12_finalpos$pos/ind_chr12_snp2$`SNP Start`, type = "l")

#Japonica final genetic map
jap_chr1 <- jap_chr1_finalpos$pos/100
jap_chr1len <- length(jap_chr1)
dim(jap_chr1) <- c(jap_chr1len,1)
jap_chr1 <- list(jap_chr1)

jap_chr2 <- jap_chr2_finalpos$pos/100
jap_chr2len <- length(jap_chr2)
dim(jap_chr2) <- c(jap_chr2len,1)
jap_chr2 <- list(jap_chr2)

jap_chr3 <- jap_chr3_finalpos$pos/100
jap_chr3len <- length(jap_chr3)
dim(jap_chr3) <- c(jap_chr3len,1)
jap_chr3 <- list(jap_chr3)

jap_chr4 <- jap_chr4_finalpos$pos/100
jap_chr4len <- length(jap_chr4)
dim(jap_chr4) <- c(jap_chr4len,1)
jap_chr4 <- list(jap_chr4)

jap_chr5 <- jap_chr5_finalpos$pos/100
jap_chr5len <- length(jap_chr5)
dim(jap_chr5) <- c(jap_chr5len,1)
jap_chr5 <- list(jap_chr5)

jap_chr5 <- jap_chr5_finalpos$pos/100
jap_chr5len <- length(jap_chr5)
dim(jap_chr5) <- c(jap_chr5len,1)
jap_chr5 <- list(jap_chr5)

jap_chr7 <- jap_chr7_finalpos$pos/100
jap_chr7len <- length(jap_chr7)
dim(jap_chr7) <- c(jap_chr7len,1)
jap_chr7 <- list(jap_chr7)

jap_chr8 <- jap_chr8_finalpos$pos/100
jap_chr8len <- length(jap_chr8)
dim(jap_chr8) <- c(jap_chr8len,1)
jap_chr8 <- list(jap_chr8)

jap_chr9 <- jap_chr9_finalpos$pos/100
jap_chr9len <- length(jap_chr9)
dim(jap_chr9) <- c(jap_chr9len,1)
jap_chr9 <- list(jap_chr9)

jap_chr10 <- jap_chr10_finalpos$pos/100
jap_chr10len <- length(jap_chr10)
dim(jap_chr10) <- c(jap_chr10len,1)
jap_chr10 <- list(jap_chr10)

jap_chr11 <- jap_chr11_finalpos$pos/100
jap_chr11len <- length(jap_chr11)
dim(jap_chr11) <- c(jap_chr11len,1)
jap_chr11 <- list(jap_chr11)

jap_chr12 <- jap_chr12_finalpos$pos/100
jap_chr12len <- length(jap_chr12)
dim(jap_chr12) <- c(jap_chr12len,1)
jap_chr12 <- list(jap_chr12)

final_map <- list(jap_chr1[[1]], jap_chr2[[1]], 
                  jap_chr3[[1]], jap_chr4[[1]], jap_chr5[[1]], 
                  jap_chr5[[1]], jap_chr7[[1]], jap_chr8[[1]], 
                  jap_chr9[[1]], jap_chr10[[1]],jap_chr11[[1]], jap_chr12[[1]])

#Indica final genetic map 
ind_chr1 <- ind_chr1_finalpos$pos/100
ind_chr1len <- length(ind_chr1)
dim(ind_chr1) <- c(ind_chr1len,1)
ind_chr1 <- list(ind_chr1)

ind_chr2 <- ind_chr2_finalpos$pos/100
ind_chr2len <- length(ind_chr2)
dim(ind_chr2) <- c(ind_chr2len,1)
ind_chr2 <- list(ind_chr2)

ind_chr3 <- ind_chr3_finalpos$pos/100
ind_chr3len <- length(ind_chr3)
dim(ind_chr3) <- c(ind_chr3len,1)
ind_chr3 <- list(ind_chr3)

ind_chr4 <- ind_chr4_finalpos$pos/100
ind_chr4len <- length(ind_chr4)
dim(ind_chr4) <- c(ind_chr4len,1)
ind_chr4 <- list(ind_chr4)

chr10 <- chr10_finalpos$pos/100
chr10len <- length(chr10)
dim(chr10) <- c(chr10len,1)
chr10 <- list(chr10)

ind_chr5 <- ind_chr5_finalpos$pos/100
ind_chr5len <- length(ind_chr5)
dim(ind_chr5) <- c(ind_chr5len,1)
ind_chr5 <- list(ind_chr5)

ind_chr6 <- ind_chr6_finalpos$pos/100
ind_chr6len <- length(ind_chr6)
dim(ind_chr6) <- c(ind_chr6len,1)
ind_chr6 <- list(ind_chr6)

ind_chr7 <- ind_chr7_finalpos$pos/100
ind_chr7len <- length(ind_chr7)
dim(ind_chr7) <- c(ind_chr7len,1)
ind_chr7 <- list(ind_chr7)

ind_chr8 <- ind_chr8_finalpos$pos/100
ind_chr8len <- length(ind_chr8)
dim(ind_chr8) <- c(ind_chr8len,1)
ind_chr8 <- list(ind_chr8)

ind_chr9 <- ind_chr9_finalpos$pos/100
ind_chr9len <- length(ind_chr9)
dim(ind_chr9) <- c(ind_chr9len,1)
ind_chr9 <- list(ind_chr9)

ind_chr10 <- ind_chr10_finalpos$pos/100
ind_chr10len <- length(ind_chr10)
dim(ind_chr10) <- c(ind_chr10len,1)
ind_chr10 <- list(ind_chr10)

ind_chr11 <- ind_chr11_finalpos$pos/100
ind_chr11len <- length(ind_chr11)
dim(ind_chr11) <- c(ind_chr11len,1)
ind_chr11 <- list(ind_chr11)

ind_chr12 <- ind_chr12_finalpos$pos/100
ind_chr12len <- length(ind_chr12)
dim(ind_chr12) <- c(ind_chr12len,1)
ind_chr12 <- list(ind_chr12)

final_map <- list(ind_chr1[[1]], ind_chr2[[1]], 
                  ind_chr3[[1]], ind_chr4[[1]], ind_ind_chr5[[1]], 
                  ind_chr6[[1]], ind_ind_chr7[[1]], ind_chr8[[1]], 
                  ind_chr9[[1]], ind_chr10[[1]], ind_chr11[[1]], ind_chr12[[1]])

