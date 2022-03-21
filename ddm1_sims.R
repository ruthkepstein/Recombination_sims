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
#write.csv(final_snps, "C:/Users/16192/Documents/PNAS_Simulations/final_snps.csv", row.names = FALSE)

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

#bin crossovers into 300 uneven bins
#normalize CO count!
chr1_CO <- chr1_CO[order(chr1_CO$`CO Start`),]
chr1_bin <- binning(chr1_CO$midpoint, nbins = max(chr1_snp$`SNP End`)/5000000, type = "kmeans")
chr1_bin <- as.data.frame(summary(chr1_bin))
#chr1_bin$freq <- chr1_bin$freq/nrow(chr1_CO)
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
chr1_bin[max(chr1_snp$`SNP End`)/5000000,5] <- max(chr1_snp$`SNP End`)
#adding length of bin as column and making in Mb
chr1_bin$length <- (chr1_bin$foo.X2-chr1_bin$foo.X1)/1000000
chr1_bin$rate <- ((chr1_bin$freq/4713)*100)/chr1_bin$length

chr2_CO <- chr2_CO[order(chr2_CO$`CO Start`),]
chr2_bin <- binning(chr2_CO$midpoint, nbins = max(chr2_snp$`SNP End`)/1000000, type = "kmeans")
chr2_bin <- as.data.frame(summary(chr2_bin))
#chr2_bin$freq <- chr2_bin$freq/nrow(chr2_CO)
chr2_bin <- within(chr2_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr2_bin$levels), ',', fixed=TRUE))))
chr2_bin <- do.call(data.frame, chr2_bin)
chr2_bin <- chr2_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr2_bin <- chr2_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr2_bin[1,4] <- 440104
chr2_bin$foo.X1 <- chr2_bin$foo.X1 - 440104
chr2_bin$foo.X2 <- chr2_bin$foo.X2 - 440104
chr2_bin[max(chr2_snp$`SNP End`)/1000000,5] <- max(chr2_snp$`SNP End`)
chr2_bin$length <- (chr2_bin$foo.X2-chr2_bin$foo.X1)/1000000
chr2_bin$rate <- ((chr2_bin$freq/4713)*100)/chr2_bin$length

#have not converted rest of chromosomes to what I did in 1 & 2
chr3_CO <- chr3_CO[order(chr3_CO$`CO Start`),]
chr3_bin <- binning(chr3_CO$midpoint, nbins = max(chr3_snp$`SNP End`)/1000000, type = "kmeans")
chr3_bin <- as.data.frame(summary(chr3_bin))
#chr3_bin$freq <- chr3_bin$freq/nrow(chr3_CO)
chr3_bin <- within(chr3_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr3_bin$levels), ',', fixed=TRUE))))
chr3_bin <- do.call(data.frame, chr3_bin)
chr3_bin <- chr3_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr3_bin <- chr3_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr3_bin[1,4] <- 865390
chr3_bin$foo.X1 <- chr3_bin$foo.X1 - 865390
chr3_bin$foo.X2 <- chr3_bin$foo.X2 - 865390
chr3_bin[max(chr3_snp$`SNP End`)/1000000,5] <- max(chr3_snp$`SNP End`)
chr3_bin$length <- (chr3_bin$foo.X2-chr3_bin$foo.X1)/1000000
chr3_bin$rate <- ((chr3_bin$freq/4713)*100)/chr3_bin$length

chr4_CO <- chr4_CO[order(chr4_CO$`CO Start`),]
chr4_bin <- binning(chr4_CO$midpoint, nbins = max(chr4_snp$`SNP End`)/1000000, type = "kmeans")
chr4_bin <- as.data.frame(summary(chr4_bin))
#chr4_bin$freq <- chr4_bin$freq/nrow(chr4_CO)
chr4_bin <- within(chr4_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr4_bin$levels), ',', fixed=TRUE))))
chr4_bin <- do.call(data.frame, chr4_bin)
chr4_bin <- chr4_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr4_bin <- chr4_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr4_bin[1,4] <- 272401
chr4_bin$foo.X1 <- chr4_bin$foo.X1 - 272401
chr4_bin$foo.X2 <- chr4_bin$foo.X2 - 272401
chr4_bin[max(chr2_snp$`SNP End`)/1000000,5] <- max(chr4_snp$`SNP End`)
chr4_bin$length <- (chr4_bin$foo.X2-chr4_bin$foo.X1)/1000000
chr4_bin$rate <- ((chr4_bin$freq/4713)*100)/chr4_bin$length

chr5_CO <- chr5_CO[order(chr5_CO$`CO Start`),]
chr5_bin <- binning(chr5_CO$midpoint, nbins = max(chr5_snp$`SNP End`)/1000000, type = "kmeans")
chr5_bin <- as.data.frame(summary(chr5_bin))
#chr5_bin$freq <- chr5_bin$freq/nrow(chr5_CO)
chr5_bin <- within(chr5_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr5_bin$levels), ',', fixed=TRUE))))
chr5_bin <- do.call(data.frame, chr5_bin)
chr5_bin <- chr5_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr5_bin <- chr5_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr5_bin[1,4] <- 267335.5
chr5_bin$foo.X1 <- chr5_bin$foo.X1 - 267335.5
chr5_bin$foo.X2 <- chr5_bin$foo.X2 - 267335.5
chr5_bin[max(chr2_snp$`SNP End`)/1000000,5] <- max(chr5_snp$`SNP End`)
chr5_bin$length <- (chr5_bin$foo.X2-chr5_bin$foo.X1)/1000000
chr5_bin$rate <- ((chr5_bin$freq/4713)*100)/chr5_bin$length

chr6_CO <- chr6_CO[order(chr6_CO$`CO Start`),]
chr6_bin <- binning(chr6_CO$midpoint, nbins = max(chr6_snp$`SNP End`)/1000000, type = "kmeans")
chr6_bin <- as.data.frame(summary(chr6_bin))
#chr6_bin$freq <- chr6_bin$freq/nrow(chr6_CO)
chr6_bin <- within(chr6_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr6_bin$levels), ',', fixed=TRUE))))
chr6_bin <- do.call(data.frame, chr6_bin)
chr6_bin <- chr6_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr6_bin <- chr6_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr6_bin[1,4] <- 197266.5
chr6_bin$foo.X1 <- chr6_bin$foo.X1 - 197266.5
chr6_bin$foo.X2 <- chr6_bin$foo.X2 - 197266.5
chr6_bin[max(chr6_snp$`SNP End`)/1000000,5] <- max(chr6_snp$`SNP End`)
chr6_bin$length <- (chr6_bin$foo.X2-chr6_bin$foo.X1)/1000000
chr6_bin$rate <- ((chr6_bin$freq/4713)*100)/chr6_bin$length

chr7_CO <- chr7_CO[order(chr7_CO$`CO Start`),]
chr7_bin <- binning(chr7_CO$midpoint, nbins = max(chr7_snp$`SNP End`)/1000000, type = "kmeans")
chr7_bin <- as.data.frame(summary(chr7_bin))
#chr7_bin$freq <- chr7_bin$freq/nrow(chr7_CO)
chr7_bin <- within(chr7_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr7_bin$levels), ',', fixed=TRUE))))
chr7_bin <- do.call(data.frame, chr7_bin)
chr7_bin <- chr7_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr7_bin <- chr7_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr7_bin[1,4] <- 375904
chr7_bin$foo.X1 <- chr7_bin$foo.X1 - 375904
chr7_bin$foo.X2 <- chr7_bin$foo.X2 - 375904
chr7_bin[max(chr7_snp$`SNP End`)/1000000,5] <- max(chr7_snp$`SNP End`)
chr7_bin$length <- (chr7_bin$foo.X2-chr7_bin$foo.X1)/1000000
chr7_bin$rate <- ((chr7_bin$freq/4713)*100)/chr7_bin$length

chr8_CO <- chr8_CO[order(chr8_CO$`CO Start`),]
chr8_bin <- binning(chr8_CO$midpoint, nbins = max(chr8_snp$`SNP End`)/1000000, type = "kmeans")
chr8_bin <- as.data.frame(summary(chr8_bin))
#chr8_bin$freq <- chr8_bin$freq/nrow(chr8_CO)
chr8_bin <- within(chr8_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr8_bin$levels), ',', fixed=TRUE))))
chr8_bin <- do.call(data.frame, chr8_bin)
chr8_bin <- chr8_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr8_bin <- chr8_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr8_bin[1,4] <- 132181
chr8_bin$foo.X1 <- chr8_bin$foo.X1 - 132181
chr8_bin$foo.X2 <- chr8_bin$foo.X2 - 132181
chr8_bin[max(chr8_snp$`SNP End`)/1000000,5] <- max(chr8_snp$`SNP End`)
chr8_bin$length <- (chr8_bin$foo.X2-chr8_bin$foo.X1)/1000000
chr8_bin$rate <- ((chr8_bin$freq/4713)*100)/chr8_bin$length

chr9_CO <- chr9_CO[order(chr9_CO$`CO Start`),]
chr9_bin <- binning(chr9_CO$midpoint, nbins = max(chr9_snp$`SNP End`)/1000000, type = "kmeans")
chr9_bin <- as.data.frame(summary(chr9_bin))
#chr9_bin$freq <- chr9_bin$freq/nrow(chr9_CO)
chr9_bin <- within(chr9_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr9_bin$levels), ',', fixed=TRUE))))
chr9_bin <- do.call(data.frame, chr9_bin)
chr9_bin <- chr9_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr9_bin <- chr9_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr9_bin[1,4] <- 317217.5
chr9_bin$foo.X1 <- chr9_bin$foo.X1 - 317217.5
chr9_bin$foo.X2 <- chr9_bin$foo.X2 - 317217.5
chr9_bin[max(chr9_snp$`SNP End`)/1000000,5] <- max(chr9_snp$`SNP End`)
chr9_bin$length <- (chr9_bin$foo.X2-chr9_bin$foo.X1)/1000000
chr9_bin$rate <- ((chr9_bin$freq/4713)*100)/chr9_bin$length

chr10_CO <- chr10_CO[order(chr10_CO$`CO Start`),]
chr10_bin <- binning(chr10_CO$midpoint, nbins = max(chr10_snp$`SNP End`)/1000000, type = "kmeans")
chr10_bin <- as.data.frame(summary(chr10_bin))
#chr10_bin$freq <- chr10_bin$freq/nrow(chr10_CO)
chr10_bin <- within(chr10_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr10_bin$levels), ',', fixed=TRUE))))
chr10_bin <- do.call(data.frame, chr10_bin)
chr10_bin <- chr10_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr10_bin <- chr10_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr10_bin[1,4] <- 698530
chr10_bin$foo.X1 <- chr10_bin$foo.X1 - 698530
chr10_bin$foo.X2 <- chr10_bin$foo.X2 - 698530
chr10_bin[max(chr10_snp$`SNP End`)/1000000,5] <- max(chr10_snp$`SNP End`)
chr10_bin$length <- (chr10_bin$foo.X2-chr10_bin$foo.X1)/1000000
chr10_bin$rate <- ((chr10_bin$freq/4713)*100)/chr10_bin$length
#use plot to look at distribution of k-means

##assigning frequency to SNPs based on recombination frequency in each bin
snp_rate_ddm1 <- function(chr_w_ddm1, chr_snp){
  for(i in 1:nrow(chr_snp)){
    for(k in 1:nrow(chr_w_ddm1)){
      if(isTRUE((chr_snp$`SNP Start`[i] >= chr_w_ddm1$foo.X1[k]) && (chr_snp$`SNP Start`[i] <= chr_w_ddm1$foo.X2[k]))){
        chr_snp$rate[i] <- chr_w_ddm1$final[k]
      }
    }
  }
  print(chr_snp)
}

##Using ddm1 recombination landscape--> 6% increase in COs
ddm1_dist <- read.table("maize_genome_ddm1_zmet2.txt", header = FALSE)
colnames(ddm1_dist) <- c("Chr", "Start", "End", "Female WT", "Male WT", "ddm1_1", "ddm1_2", "zmet2")
#normalizing the data
ddm1_dist$ddm1_1 <- ddm1_dist$ddm1_1*2/39
ddm1_dist$ddm1_2 <- ddm1_dist$ddm1_2*2/40
#ddm1_dist$`Female WT` <- ddm1_dist$`Female WT`*2/122
ddm1_dist$`Male WT`<- ddm1_dist$`Male WT`*2/135
#ddm1_dist$WT <- (ddm1_dist$`Male WT`)/2
ddm1_dist$ddm1 <- (ddm1_dist$ddm1_1+ddm1_dist$ddm1_2)/2

ddm1_dist$diffwt <- 0
#ddm1_dist$diffwt2 <- ddm1_dist$diffwt*1
#ddm1_dist$diffwt2[ddm1_dist$diffwt2 <= 0] <- 0
#ddm1_dist$diffwt2 <- abs(ddm1_dist$diffwt2)
#~160% increase in telomeres

#make generalization about ddm1

chr1_distddm1 <- ddm1_dist[ which(ddm1_dist$Chr == 1),]
chr2_distddm1 <- ddm1_dist[ which(ddm1_dist$Chr == 2),]
chr3_distddm1 <- ddm1_dist[ which(ddm1_dist$Chr == 3),]
chr4_distddm1 <- ddm1_dist[ which(ddm1_dist$Chr == 4),]
chr5_distddm1 <- ddm1_dist[ which(ddm1_dist$Chr == 5),]
chr6_distddm1 <- ddm1_dist[ which(ddm1_dist$Chr == 6),]
chr7_distddm1 <- ddm1_dist[ which(ddm1_dist$Chr == 7),]
chr8_distddm1 <- ddm1_dist[ which(ddm1_dist$Chr == 8),]
chr9_distddm1 <- ddm1_dist[ which(ddm1_dist$Chr == 9),]
chr10_distddm1 <- ddm1_dist[ which(ddm1_dist$Chr == 10),]

chr1_distddm1[1:9,10] <- 1.6
chr1_distddm1[50:60,10] <- 1.6

chr2_distddm1[1:8,10] <- 1.6
chr2_distddm1[39:47,10] <- 1.6

chr3_distddm1[1:7,10] <- 1.6
chr3_distddm1[39:46,10] <- 1.6

chr4_distddm1[1:8,10] <- 1.6
chr4_distddm1[40:48,10] <- 1.6

chr5_distddm1[1:7,10] <- 1.6
chr5_distddm1[36:43,10] <- 1.6

chr6_distddm1[1:5,10] <- 1.6
chr6_distddm1[28:33,10] <- 1.6

chr7_distddm1[1:6,10] <- 1.6
chr7_distddm1[29:35,10] <- 1.6

chr8_distddm1[1:6,10] <- 1.6
chr8_distddm1[29:35,10] <- 1.6

chr9_distddm1[1:5,10] <- 1.6
chr9_distddm1[26:31,10] <- 1.6

chr10_distddm1[1:5,10] <- 1.6
chr10_distddm1[24:29,10] <- 1.6


ddm1_wt <- function(chr_bin, ddm1_dist){
  for(i in 1:nrow(chr_bin)){
    for(k in 1:nrow(ddm1_dist)){
      if(isTRUE(chr_bin$foo.X1[i] >= ddm1_dist$Start[k] && chr_bin$foo.X2 <= ddm1_dist$End[k])){
        chr_bin$final[i] <- chr_bin$rate[i] + (chr_bin$rate[i]*ddm1_dist$diffwt[k])
      }
    }
  }
  return(chr_bin)
}
chr1_w_ddm1 <- ddm1_wt(chr1_bin, chr1_distddm1)
chr2_w_ddm1 <- ddm1_wt(chr2_bin, chr2_distddm1)
chr3_w_ddm1 <- ddm1_wt(chr3_bin, chr3_distddm1)
chr4_w_ddm1 <- ddm1_wt(chr4_bin, chr4_distddm1)
chr5_w_ddm1 <- ddm1_wt(chr5_bin, chr5_distddm1)
chr6_w_ddm1 <- ddm1_wt(chr6_bin, chr6_distddm1)
chr7_w_ddm1 <- ddm1_wt(chr7_bin, chr7_distddm1)
chr8_w_ddm1 <- ddm1_wt(chr8_bin, chr8_distddm1)
chr9_w_ddm1 <- ddm1_wt(chr9_bin, chr9_distddm1)
chr10_w_ddm1 <- ddm1_wt(chr10_bin, chr10_distddm1)

#using function, converted SNP start to Mb to get cM/Mb for final genetic position
chr1_snp2 <- snp_rate_ddm1(chr1_w_ddm1, chr1_snp)
chr1_snp2$`SNP Start`<- chr1_snp2$`SNP Start`/1000000
chr1_snp2 <- chr1_snp2[order(chr1_snp2$`SNP Start`),]
#smoothing the recombination rate so transitions between bins are not so abrupt
chr1_spl2 <- smooth.spline(chr1_snp2$rate, spar = 1.25)
#creation of genetic positions from smoothed recombination rate
chr1_snp2$pos <- (chr1_snp2$`SNP Start`*chr1_spl2$y)
#graph to look at Mb vs. cM along chromosome
plot(chr1_snp2$`SNP Start`, chr1_snp2$pos, type = "l", col = "blue", 
     main = "Chr1. Genetic Maps", xlab = "Physical Positions (Mb)", ylab = "Recombination rate (cM/Mb)")
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


chr2_snp2 <- snp_rate_ddm1(chr2_w_ddm1, chr2_snp)
chr2_snp2$`SNP Start` <- chr2_snp2$`SNP Start`/1000000
chr2_snp2 <- chr2_snp2[-(228:237),]
chr2_spl <- smooth.spline(chr2_snp2$rate, spar = 1.25)
chr2_snp2$pos <- (chr2_snp2$`SNP Start`*chr2_spl$y)
plot(chr2_snp2$`SNP Start`, chr2_snp2$pos)
plot(chr2_snp2$`SNP Start`, chr2_snp2$pos/chr2_snp2$`SNP Start`, type = "l")
chr2_finalpos <- chr2_snp2[order(chr2_snp2$pos),]
is.unsorted(chr2_finalpos$pos)
plot(chr2_snp2$`SNP Start`, chr2_finalpos$pos/chr2_snp2$`SNP Start`, type = "l")

chr3_snp2 <- snp_rate_ddm1(chr3_w_ddm1, chr3_snp)
chr3_snp2$`SNP Start` <- chr3_snp2$`SNP Start`/1000000
chr3_spl <- smooth.spline(chr3_snp2$rate, spar = 1.25)
chr3_snp2$pos <- (chr3_snp2$`SNP Start`*chr3_spl$y)
plot(chr3_snp2$`SNP Start`, chr3_snp2$pos)
plot(chr3_snp2$`SNP Start`, chr3_snp2$pos/chr3_snp2$`SNP Start`, type = "l")
chr3_finalpos <- chr3_snp2[order(chr3_snp2$pos),]
is.unsorted(chr3_finalpos$pos)
plot(chr3_snp2$`SNP Start`, chr3_finalpos$pos/chr3_snp2$`SNP Start`, type = "l")

chr4_snp2 <- snp_rate_ddm1(chr4_w_ddm1, chr4_snp)
chr4_snp2$`SNP Start` <- chr4_snp2$`SNP Start`/1000000
chr4_spl <- smooth.spline(chr4_snp2$rate, spar = 1.25)
chr4_snp2$pos <- (chr4_snp2$`SNP Start`*chr4_spl$y)
plot(chr4_snp2$`SNP Start`, chr4_snp2$pos)
plot(chr4_snp2$`SNP Start`, chr4_snp2$pos/chr4_snp2$`SNP Start`, type = "l")
chr4_finalpos <- chr4_snp2[order(chr4_snp2$pos),]
is.unsorted(chr4_finalpos$pos)
plot(chr4_snp2$`SNP Start`, chr4_finalpos$pos/chr4_snp2$`SNP Start`, type = "l")

chr5_snp2 <- snp_rate_ddm1(chr5_w_ddm1, chr5_snp)
chr5_snp2$`SNP Start` <- chr5_snp2$`SNP Start`/1000000
chr5_spl <- smooth.spline(chr5_snp2$rate, spar = 1.2)
chr5_snp2$pos <- (chr5_snp2$`SNP Start`*chr5_spl$y)
plot(chr5_snp2$`SNP Start`, chr5_snp2$pos)
plot(chr5_snp2$`SNP Start`, chr5_snp2$pos/chr5_snp2$`SNP Start`, type = "l")
chr5_finalpos <- chr5_snp2[order(chr5_snp2$pos),]
is.unsorted(chr5_finalpos$pos)
plot(chr5_snp2$`SNP Start`, chr5_finalpos$pos/chr5_snp2$`SNP Start`, type = "l")

#chr 6 is lowkey fuked up
chr6_snp2 <- snp_rate_ddm1(chr6_w_ddm1, chr6_snp)
chr6_snp2$`SNP Start` <- chr6_snp2$`SNP Start`/1000000
chr6_spl <- smooth.spline(chr6_snp2$rate, spar = 1)
chr6_snp2$pos <- (chr6_snp2$`SNP Start`*chr6_spl$y)
plot(chr6_snp2$`SNP Start`, chr6_snp2$pos)
plot(chr6_snp2$`SNP Start`, chr6_snp2$pos/chr6_snp2$`SNP Start`, type = "l")
chr6_finalpos <- chr6_snp2[order(chr6_snp2$pos),]
is.unsorted(chr6_finalpos$pos)
plot(chr6_snp2$`SNP Start`, chr6_finalpos$pos/chr6_snp2$`SNP Start`, type = "l")

chr7_snp2 <- snp_rate_ddm1(chr7_w_ddm1, chr7_snp)
chr7_snp2$`SNP Start` <- chr7_snp2$`SNP Start`/1000000
chr7_spl <- smooth.spline(chr7_snp2$rate, spar = 1.2)
chr7_snp2$pos <- (chr7_snp2$`SNP Start`*chr7_spl$y)
plot(chr7_snp2$`SNP Start`, chr7_snp2$pos)
plot(chr7_snp2$`SNP Start`, chr7_snp2$pos/chr7_snp2$`SNP Start`, type = "l")
chr7_finalpos <- chr7_snp2[order(chr7_snp2$pos),]
is.unsorted(chr7_finalpos$pos)
plot(chr7_snp2$`SNP Start`, chr7_finalpos$pos/chr7_snp2$`SNP Start`, type = "l")

chr8_snp2 <- snp_rate_ddm1(chr8_w_ddm1, chr8_snp)
chr8_snp2$`SNP Start` <- chr8_snp2$`SNP Start`/1000000
chr8_spl <- smooth.spline(chr8_snp2$rate, spar = 1.2)
chr8_snp2$pos <- (chr8_snp2$`SNP Start`*chr8_spl$y)
plot(chr8_snp2$`SNP Start`, chr8_snp2$pos)
plot(chr8_snp2$`SNP Start`, chr8_snp2$pos/chr8_snp2$`SNP Start`, type = "l")
chr8_finalpos <- chr8_snp2[order(chr8_snp2$pos),]
is.unsorted(chr8_finalpos$pos)
plot(chr8_snp2$`SNP Start`, chr8_finalpos$pos/chr8_snp2$`SNP Start`, type = "l")

chr9_snp2 <- snp_rate_ddm1(chr9_w_ddm1, chr9_snp)
chr9_snp2$`SNP Start` <- chr9_snp2$`SNP Start`/1000000
chr9_spl <- smooth.spline(chr9_snp2$rate, spar = 1.15)
chr9_snp2$pos <- (chr9_snp2$`SNP Start`*chr9_spl$y)
plot(chr9_snp2$`SNP Start`, chr9_snp2$pos)
plot(chr9_snp2$`SNP Start`, chr9_snp2$pos/chr9_snp2$`SNP Start`, type = "l")
chr9_finalpos <- chr9_snp2[order(chr9_snp2$pos),]
is.unsorted(chr9_finalpos$pos)
plot(chr9_snp2$`SNP Start`, chr9_finalpos$pos/chr9_snp2$`SNP Start`, type = "l")

chr10_snp2 <- snp_rate_ddm1(chr10_w_ddm1, chr10_snp)
chr10_snp2$`SNP Start` <- chr10_snp2$`SNP Start`/1000000
chr10_spl <- smooth.spline(chr10_snp2$rate, spar = 1.2)
chr10_snp2$pos <- (chr10_snp2$`SNP Start`*chr10_spl$y)
plot(chr10_snp2$`SNP Start`, chr10_snp2$pos)
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

ddm1_map <- list(chr1[[1]], chr2[[1]], 
                  chr3[[1]], chr4[[1]], chr5[[1]], 
                  chr6[[1]], chr7[[1]], chr8[[1]], 
                  chr9[[1]], chr10[[1]])

#Creating vector of centromere positions
#chnage this to reflect new pos
ddm1_centromere <- c(171.0977, 231.1713, 120.6425, 94.27555, 227.7858,
                     67.5, 79.2089, 119.8078, 137.7477, 84.6446)
ddm1_centromere <- (ddm1_centromere/100)

ddm1_sel2_gv <- matrix(data = NA, nrow = 100, ncol = 20)
ddm1_sel3_gv <- matrix(data = NA, nrow = 100, ncol = 20)
ddm1_sel4_gv <- matrix(data = NA, nrow = 100, ncol = 20)
ddm1_sel5_gv <- matrix(data = NA, nrow = 100, ncol = 20)
ddm1_sel6_gv <- matrix(data = NA, nrow = 100, ncol = 20)
ddm1_sel7_gv <- matrix(data = NA, nrow = 100, ncol = 20)
ddm1var1 <- matrix(data = NA, ncol = 2, nrow = 200)
ddm1var2 <- matrix(data = NA, ncol = 2, nrow = 200)
ddm1var3 <- matrix(data = NA, ncol = 2, nrow = 200)
ddm1var4 <- matrix(data = NA, ncol = 2, nrow = 200)
ddm1var5 <- matrix(data = NA, ncol = 2, nrow = 200)
ddm1var6 <- matrix(data = NA, ncol = 2, nrow = 200)
ddm1var7 <- matrix(data = NA, ncol = 2, nrow = 200)
for(i in 1:100){
  SP$switchGenMap(ddm1_map, centromere = ddm1_centromere)
  pop <- randCross2(goodpop[49,], badpop[155,], nCrosses = 200, nProgeny = 1, simParam = SP)
  pop <- setPheno(pop = pop, h2 = 0.8, simParam = SP)
  ddm1var1[i] <- varA(pop)
  
  pop1_sel <- selectInd(pop, nInd = 80, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_2 <- selectInd(pop1_sel, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_cross <- randCross2(pop1_sel1_2, goodpop[49,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel1_cross <- setPheno(pop1_sel1_cross, h2 = 0.8, simParam = SP)
  ddm1var2[i] <- varA(pop1_sel1_2)
  
  ddm1_sel2_gv[i,] <- gv(pop1_sel1_2)
  pop1_sel2 <- selectInd(pop1_sel1_cross, nInd = 20, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_2 <- selectInd(pop1_sel2, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_cross <- randCross2(pop1_sel2_2, goodpop[49,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel2_cross <- setPheno(pop1_sel2_cross, h2 = 0.8, simParam = SP)
  ddm1var3[i] <- varA(pop1_sel2_2)
  
  ddm1_sel3_gv[i,] <- gv(pop1_sel2_2)
  pop1_sel3 <- selectInd(pop1_sel2_cross, nInd = 20, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_2 <- selectInd(pop1_sel3, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_cross <- randCross2(pop1_sel3_2, goodpop[49,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel3_cross <- setPheno(pop1_sel3_cross, h2 = 0.8, simParam = SP)
  ddm1var4[i] <- varA(pop1_sel3_2)
  
  ddm1_sel4_gv[i,] <- gv(pop1_sel3_2)
  pop1_sel4 <- selectInd(pop1_sel3_cross, nInd = 20, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_2 <- selectInd(pop1_sel4, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_cross <- randCross2(pop1_sel4_2, goodpop[49,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel4_cross <- setPheno(pop1_sel4_cross, h2 = 0.8, simParam = SP)
  ddm1var5[i] <- varA(pop1_sel4_2)
  
  ddm1_sel5_gv[i,] <- gv(pop1_sel4_2)
  pop1_sel5 <- selectInd(pop1_sel4_cross, nInd = 20, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_2 <- selectInd(pop1_sel5, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_cross <- randCross2(pop1_sel5_2, goodpop[49,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel5_cross <- setPheno(pop1_sel5_cross, h2 = 0.8, simParam = SP)
  ddm1var6[i] <- varA(pop1_sel5_2)
  
  ddm1_sel6_gv[i,] <- gv(pop1_sel5_2)
}

wt <- c(pop1_sel2_gv[,1:10], pop1_sel3_gv[,1:10], pop1_sel4_gv[,1:10], pop1_sel5_gv[,1:10], pop1_sel6_gv[,1:10])

ddm1 <- c(ddm1_sel2_gv[,1:10], ddm1_sel3_gv[,1:10], ddm1_sel4_gv[,1:10], ddm1_sel5_gv[,1:10], ddm1_sel6_gv[,1:10])
          
zmet2 <- c(zmet2_sel2_gv[,1:10], zmet2_sel3_gv[,1:10], zmet2_sel4_gv[,1:10], zmet2_sel5_gv[,1:10], zmet2_sel6_gv[,1:10])

recq4 <- c(recq4_sel2_gv[,1:10], recq4_sel3_gv[,1:10], recq4_sel4_gv[,1:10], recq4_sel5_gv[,1:10], recq4_sel6_gv[,1:10])


wt2 <- c(pop1_sel2_gv[,11:20], pop1_sel3_gv[,11:20], pop1_sel4_gv[,11:20], pop1_sel5_gv[,11:20], pop1_sel6_gv[,11:20])

ddm12 <- c(ddm1_sel2_gv[,11:20], ddm1_sel3_gv[,11:20], ddm1_sel4_gv[,11:20], ddm1_sel5_gv[,11:20], ddm1_sel6_gv[,11:20])

zmet22 <- c(zmet2_sel2_gv[,11:20], zmet2_sel3_gv[,11:20], zmet2_sel4_gv[,11:20], zmet2_sel5_gv[,11:20], zmet2_sel6_gv[,11:20])

recq42 <- c(recq4_sel2_gv[,11:20], recq4_sel3_gv[,11:20], recq4_sel4_gv[,11:20], recq4_sel5_gv[,11:20], recq4_sel6_gv[,11:20])


wtvarall <- c(wtvar1, wtvar2, wtvar3, wtvar4, wtvar5,
              wtvar6)

ddm1varall <- c(ddm1var1, ddm1var2, ddm1var3, ddm1var4, ddm1var5,
                ddm1var6)

zmet2varall <- c(zmet2var1, zmet2var2, zmet2var3, zmet2var4, zmet2var5,
                 zmet2var6)

recq4varall <- c(recq4var1, recq4var2, recq4var3, recq4var4, recq4var5,
                 recq4var6)

all <- cbind(wt, ddm1, zmet2, recq4)
all <- as.data.frame(all)
all2 <- as.data.frame(matrix(data = NA, nrow = 20000))
all2$gv <- c(all$wt, all$ddm1, all$zmet2, all$recq4)
all2$generation <- rep(1:5, size = 4, each = 1000)
all2$gen_map <- rep(c("wt","ddm1", "zmet2", "recq4"), size = 4, each = 5000)
gen2 <- all2[which(all2$generation == 5), ]

reduced_introgression <- aov(gv ~ as.factor(gen_map), data = gen2)
full_introgression <- aov(gv ~ as.factor(gen_map)*generation, data = all2)
introgression <- anova(reduced_introgression, full_introgression)
summary(reduced_introgression)
TukeyHSD(full_introgression, conf.level=.95)
TukeyHSD(reduced_introgression, conf.level=.95)

trait2 <- cbind(wt2, ddm12, zmet22, recq42)
trait2 <- as.data.frame(trait2)
trait22 <- as.data.frame(matrix(data = NA, nrow = 20000))
trait22$gv <- c(trait2$wt2, trait2$ddm12, trait2$zmet22, trait2$recq42)
trait22$generation <- rep(1:5, size = 4, each = 1000)
trait22$gen_map <- rep(c("wt","ddm1", "zmet2", "recq4"), size = 4, each = 5000)
intro_lastgen <- trait22[which(trait22$generation == 5), ]

#many extreme values
ggqqplot(residuals(reduced_introgression))

ggplot(all2, aes(x=as.factor(generation), y= gv, fill=gen_map)) + 
  geom_boxplot() + theme_bw() + xlab("Generation") + ylab("GVs") + 
  scale_fill_manual(values=c("#aae4c2", "#fff87a", "#35c7ff", '#ff8b77')) + ggtitle("Introgression Genetic Values through 6 Generations")+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=14), legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14), plot.title = element_text(size=20))

ggplot(trait22, aes(x=as.factor(generation), y= gv, fill=gen_map)) + 
  geom_count(aes(fill=gen_map))+ theme_bw() + xlab("Generation") + ylab("GVs")
  #scale_fill_manual(values=c("#aae4c2", "#fff87a", "#35c7ff", "#ff8b77")) + ggtitle("Introgression Trait 2 GVs through 6 Generations")

allvar <- cbind(wtvarall, ddm1varall, zmet2varall, recq4varall)
allvar <- as.data.frame(allvar)
allvar2 <- as.data.frame(matrix(data = NA, nrow = 9600))
allvar2$var <- c(allvar$wtvarall, allvar$ddm1varall, allvar$zmet2varall, allvar$recq4varall)
allvar2$gen <- rep(1:5, size = 4, each = 480)
allvar2$gen_map <- rep(c("wt","ddm1", "zmet", "recq4"), size = 4, each = 2400)
lastvar <- allvar2[which(allvar2$gen == 5),]

introgressionvar <- aov(var ~ as.factor(gen_map), data = lastvar)
summary(introgressionvar)
TukeyHSD(introgressionvar, conf.level = .95)

ggqqplot(residuals(introgressionvar))

ggplot(allvar2, aes(x=as.factor(gen), y=var, fill=gen_map)) + 
  geom_boxplot() + theme_bw() + xlab("Generation") + ylab("Genetic Variances") + 
  scale_fill_manual(values=c("#aae4c2", "#fff87a", "#35c7ff",'#ff8b77')) + ggtitle("Introgression Genetic Variance through 6 Generations") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=14), legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14), plot.title = element_text(size=20)) + ylim(-0.05,0.2)
