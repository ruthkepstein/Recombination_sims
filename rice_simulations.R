library(AlphaSimR)
library(Rcpp)
library(ggplot2)
library(dplyr)

setwd("C:/Users/16192/Documents/PNAS_Simulations")
set.seed(420)

japonica_snps <- read.table("japonica_SNPs.bed", header =FALSE)
indica_snps <- read.table("indica_snps.bed", header =FALSE)
colnames(japonica_snps) <- c("Chr#", "SNP Start", "SNP End")
colnames(indica_snps) <- c("Chr#", "SNP Start", "SNP End")
#sample SNPs?
#final_snps <- sample_n(final_snps, 2000)
japonica_snps <- japonica_snps[order(japonica_snps$`Chr#`,japonica_snps$`SNP Start`),]
indica_snps <- indica_snps[order(indica_snps$`Chr#`,indica_snps$`SNP Start`),]

#Japonica
japonica_chr1_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr1"),]
japonica_chr1_snp$rate <- NA
japonica_chr1_snp$`SNP End` <- japonica_chr1_snp$`SNP End` - min(japonica_chr1_snp$`SNP Start`)
japonica_chr1_snp$`SNP Start` <- japonica_chr1_snp$`SNP Start`- min(japonica_chr1_snp$`SNP Start`)

japonica_chr2_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr2"),]
japonica_chr2_snp$rate <- NA
japonica_chr2_snp$`SNP End` <- japonica_chr2_snp$`SNP End` - min(japonica_chr2_snp$`SNP Start`)
japonica_chr2_snp$`SNP Start` <- japonica_chr2_snp$`SNP Start`- min(japonica_chr2_snp$`SNP Start`)

japonica_chr3_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr3"),]
japonica_chr3_snp$rate <- NA
japonica_chr3_snp$`SNP End` <- japonica_chr3_snp$`SNP End` - min(japonica_chr3_snp$`SNP Start`)
japonica_chr3_snp$`SNP Start` <- japonica_chr3_snp$`SNP Start`- min(japonica_chr3_snp$`SNP Start`)

japonica_chr4_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr4"),]
japonica_chr4_snp$rate <- NA
japonica_chr4_snp$`SNP End` <- japonica_chr4_snp$`SNP End` - min(japonica_chr4_snp$`SNP Start`)
japonica_chr4_snp$`SNP Start` <- japonica_chr4_snp$`SNP Start`- min(japonica_chr4_snp$`SNP Start`)

japonica_chr5_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr5"),]
japonica_chr5_snp$rate <- NA
japonica_chr5_snp$`SNP End` <- japonica_chr5_snp$`SNP End` - min(japonica_chr5_snp$`SNP Start`)
japonica_chr5_snp$`SNP Start` <- japonica_chr5_snp$`SNP Start`- min(japonica_chr5_snp$`SNP Start`)

japonica_chr6_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr6"),]
japonica_chr6_snp$rate <- NA
japonica_chr6_snp$`SNP End` <- japonica_chr6_snp$`SNP End` - min(japonica_chr6_snp$`SNP Start`)
japonica_chr6_snp$`SNP Start` <- japonica_chr6_snp$`SNP Start`- min(japonica_chr6_snp$`SNP Start`)

japonica_chr7_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr7"),]
japonica_chr7_snp$rate <- NA
japonica_chr7_snp$`SNP End` <- japonica_chr7_snp$`SNP End` - min(japonica_chr7_snp$`SNP Start`)
japonica_chr7_snp$`SNP Start` <- japonica_chr7_snp$`SNP Start`- min(japonica_chr7_snp$`SNP Start`)

japonica_chr8_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr8"),]
japonica_chr8_snp$rate <- NA
japonica_chr8_snp$`SNP End` <- japonica_chr8_snp$`SNP End` - min(japonica_chr8_snp$`SNP Start`)
japonica_chr8_snp$`SNP Start` <- japonica_chr8_snp$`SNP Start`- min(japonica_chr8_snp$`SNP Start`)

japonica_chr9_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr9"),]
japonica_chr9_snp$rate <- NA
japonica_chr9_snp$`SNP End` <- japonica_chr9_snp$`SNP End` - min(japonica_chr9_snp$`SNP Start`)
japonica_chr9_snp$`SNP Start` <- japonica_chr9_snp$`SNP Start`- min(japonica_chr9_snp$`SNP Start`)

japonica_chr10_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr10"),]
japonica_chr10_snp$rate <- NA
japonica_chr10_snp$`SNP End` <- japonica_chr10_snp$`SNP End` - min(japonica_chr10_snp$`SNP Start`)
japonica_chr10_snp$`SNP Start` <- japonica_chr10_snp$`SNP Start`- min(japonica_chr10_snp$`SNP Start`)

japonica_chr11_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr11"),]
japonica_chr11_snp$rate <- NA
japonica_chr11_snp$`SNP End` <- japonica_chr11_snp$`SNP End` - min(japonica_chr11_snp$`SNP Start`)
japonica_chr11_snp$`SNP Start` <- japonica_chr11_snp$`SNP Start`- min(japonica_chr11_snp$`SNP Start`)

japonica_chr12_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr12"),]
japonica_chr12_snp$rate <- NA
japonica_chr12_snp$`SNP End` <- japonica_chr12_snp$`SNP End` - min(japonica_chr12_snp$`SNP Start`)
japonica_chr12_snp$`SNP Start` <- japonica_chr12_snp$`SNP Start`- min(japonica_chr12_snp$`SNP Start`)

#indica
indica_chr1_snp <- indica_snps[ which(indica_snps$`Chr#` == "Chr1"),]
indica_chr1_snp$rate <- NA
indica_chr1_snp$`SNP End` <- indica_chr1_snp$`SNP End` - min(indica_chr1_snp$`SNP Start`)
indica_chr1_snp$`SNP Start` <- indica_chr1_snp$`SNP Start`- min(indica_chr1_snp$`SNP Start`)

indica_chr2_snp <- indica_snps[ which(indica_snps$`Chr#` == "Chr2"),]
indica_chr2_snp$rate <- NA
indica_chr2_snp$`SNP End` <- indica_chr2_snp$`SNP End` - min(indica_chr2_snp$`SNP Start`)
indica_chr2_snp$`SNP Start` <- indica_chr2_snp$`SNP Start`- min(indica_chr2_snp$`SNP Start`)

indica_chr3_snp <- indica_snps[ which(indica_snps$`Chr#` == "Chr3"),]
indica_chr3_snp$rate <- NA
indica_chr3_snp$`SNP End` <- indica_chr3_snp$`SNP End` - min(indica_chr3_snp$`SNP Start`)
indica_chr3_snp$`SNP Start` <- indica_chr3_snp$`SNP Start`- min(indica_chr3_snp$`SNP Start`)

indica_chr4_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr4"),]
indica_chr4_snp$rate <- NA
indica_chr4_snp$`SNP End` <- indica_chr4_snp$`SNP End` - min(indica_chr4_snp$`SNP Start`)
indica_chr4_snp$`SNP Start` <- indica_chr4_snp$`SNP Start`- min(indica_chr4_snp$`SNP Start`)

indica_chr5_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr5"),]
indica_chr5_snp$rate <- NA
indica_chr5_snp$`SNP End` <- indica_chr5_snp$`SNP End` - min(indica_chr5_snp$`SNP Start`)
indica_chr5_snp$`SNP Start` <- indica_chr5_snp$`SNP Start`- min(indica_chr5_snp$`SNP Start`)

indica_chr6_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr6"),]
indica_chr6_snp$rate <- NA
indica_chr6_snp$`SNP End` <- indica_chr6_snp$`SNP End` - min(indica_chr6_snp$`SNP Start`)
indica_chr6_snp$`SNP Start` <- indica_chr6_snp$`SNP Start`- min(indica_chr6_snp$`SNP Start`)

indica_chr7_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr7"),]
indica_chr7_snp$rate <- NA
indica_chr7_snp$`SNP End` <- indica_chr7_snp$`SNP End` - min(indica_chr7_snp$`SNP Start`)
indica_chr7_snp$`SNP Start` <- indica_chr7_snp$`SNP Start`- min(indica_chr7_snp$`SNP Start`)

indica_chr8_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr8"),]
indica_chr8_snp$rate <- NA
indica_chr8_snp$`SNP End` <- indica_chr8_snp$`SNP End` - min(indica_chr8_snp$`SNP Start`)
indica_chr8_snp$`SNP Start` <- indica_chr8_snp$`SNP Start`- min(indica_chr8_snp$`SNP Start`)

indica_chr9_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr9"),]
indica_chr9_snp$rate <- NA
indica_chr9_snp$`SNP End` <- indica_chr9_snp$`SNP End` - min(indica_chr9_snp$`SNP Start`)
indica_chr9_snp$`SNP Start` <- indica_chr9_snp$`SNP Start`- min(indica_chr9_snp$`SNP Start`)

indica_chr10_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr10"),]
indica_chr10_snp$rate <- NA
indica_chr10_snp$`SNP End` <- indica_chr10_snp$`SNP End` - min(indica_chr10_snp$`SNP Start`)
indica_chr10_snp$`SNP Start` <- indica_chr10_snp$`SNP Start`- min(indica_chr10_snp$`SNP Start`)

indica_chr11_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr11"),]
indica_chr11_snp$rate <- NA
indica_chr11_snp$`SNP End` <- indica_chr11_snp$`SNP End` - min(indica_chr11_snp$`SNP Start`)
indica_chr11_snp$`SNP Start` <- indica_chr11_snp$`SNP Start`- min(indica_chr11_snp$`SNP Start`)

indica_chr12_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr12"),]
indica_chr12_snp$rate <- NA
indica_chr12_snp$`SNP End` <- indica_chr12_snp$`SNP End` - min(indica_chr12_snp$`SNP Start`)
indica_chr12_snp$`SNP Start` <- indica_chr12_snp$`SNP Start`- min(indica_chr12_snp$`SNP Start`)

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
jap_chr1_bin$rate <- ((jap_chr1_bin$freq/4713)*50)/jap_chr1_bin$length

jap_chr2_bin <- binning(jap_chr2_CO$midpoint, nbins = 50, type = "kmeans")
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
jap_chr2_bin$rate <- ((jap_chr2_bin$freq/4713)*50)/jap_chr2_bin$length

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
jap_chr3_bin$rate <- ((jap_chr3_bin$freq/4713)*50)/jap_chr3_bin$length

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
jap_chr4_bin$rate <- ((jap_chr4_bin$freq/4713)*50)/jap_chr4_bin$length

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
jap_chr5_bin$rate <- ((jap_chr5_bin$freq/4713)*50)/jap_chr5_bin$length

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
jap_chr6_bin$rate <- ((jap_chr6_bin$freq/4713)*50)/jap_chr6_bin$length

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
jap_chr7_bin$rate <- ((jap_chr7_bin$freq/4713)*50)/jap_chr7_bin$length

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
jap_chr8_bin$rate <- ((jap_chr8_bin$freq/4713)*50)/jap_chr8_bin$length

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
jap_chr9_bin$rate <- ((jap_chr9_bin$freq/4713)*50)/jap_chr9_bin$length

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
jap_chr10_bin$rate <- ((jap_chr10_bin$freq/4713)*50)/jap_chr10_bin$length

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
jap_chr11_bin$rate <- ((jap_chr11_bin$freq/4713)*50)/jap_chr11_bin$length

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
jap_chr12_bin$rate <- ((jap_chr12_bin$freq/4713)*50)/jap_chr12_bin$length


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
ind_chr1_bin$rate <- ((ind_chr1_bin$freq/4713)*50)/ind_chr1_bin$length

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
ind_chr2_bin$rate <- ((ind_chr2_bin$freq/4713)*50)/ind_chr2_bin$length

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
ind_chr3_bin$rate <- ((ind_chr3_bin$freq/4713)*50)/ind_chr3_bin$length

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
ind_chr4_bin$rate <- ((ind_chr4_bin$freq/4713)*50)/ind_chr4_bin$length

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
ind_chr5_bin$rate <- ((ind_chr5_bin$freq/4713)*50)/ind_chr5_bin$length

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
ind_chr6_bin$rate <- ((ind_chr6_bin$freq/4713)*50)/ind_chr6_bin$length

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
ind_chr7_bin$rate <- ((ind_chr7_bin$freq/4713)*50)/ind_chr7_bin$length

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
ind_chr8_bin$rate <- ((ind_chr8_bin$freq/4713)*50)/ind_chr8_bin$length

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
ind_chr9_bin$rate <- ((ind_chr9_bin$freq/4713)*50)/ind_chr9_bin$length

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
ind_chr10_bin$rate <- ((ind_chr10_bin$freq/4713)*50)/ind_chr10_bin$length

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
ind_chr11_bin$rate <- ((ind_chr11_bin$freq/4713)*50)/ind_chr11_bin$length

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
ind_chr12_bin$rate <- ((ind_chr12_bin$freq/4713)*50)/ind_chr12_bin$length
