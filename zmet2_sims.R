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
chr3_bin <- chr3_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr3_bin <- chr3_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
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
chr4_bin <- chr4_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr4_bin <- chr4_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
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
chr5_bin <- chr5_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr5_bin <- chr5_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
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
chr6_bin <- chr6_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr6_bin <- chr6_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
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

##Using zmet2 recombination landscape--> 20% increase in COs
set.seed(420)
zmet2_dist <- read.table("maize_genome_ddm1_zmet2.txt", header = FALSE)
colnames(zmet2_dist) <- c("Chr", "Start", "End", "Female WT", "Male WT", "ddm1_1", "ddm1_2", "zmet2")
zmet2_dist$`Female WT` <- zmet2_dist$`Female WT`*2/122
zmet2_dist$`Male WT` <- zmet2_dist$`Male WT`*2/135
zmet2_dist$zmet2 <- zmet2_dist$zmet2*2/96
zmet2_dist$WT <- (zmet2_dist$`Female WT`+zmet2_dist$`Male WT`)/2
zmet2_dist$diff <- zmet2_dist$WT-zmet2_dist$zmet2
zmet2_dist$diff2 <- zmet2_dist$diff*1
zmet2_dist$diff2[zmet2_dist$diff2 >= 0] <- 0
zmet2_dist$diff2 <- abs(zmet2_dist$diff2)

chr1_dist <- zmet2_dist[ which(zmet2_dist$Chr == 1),]
chr2_dist <- zmet2_dist[ which(zmet2_dist$Chr == 2),]
chr3_dist <- zmet2_dist[ which(zmet2_dist$Chr == 3),]
chr4_dist <- zmet2_dist[ which(zmet2_dist$Chr == 4),]
chr5_dist <- zmet2_dist[ which(zmet2_dist$Chr == 5),]
chr6_dist <- zmet2_dist[ which(zmet2_dist$Chr == 6),]
chr7_dist <- zmet2_dist[ which(zmet2_dist$Chr == 7),]
chr8_dist <- zmet2_dist[ which(zmet2_dist$Chr == 8),]
chr9_dist <- zmet2_dist[ which(zmet2_dist$Chr == 9),]
chr10_dist <- zmet2_dist[ which(zmet2_dist$Chr == 10),]

zmet2_map <- function(chr_bin, chr_dist){
  for(i in 1:nrow(chr_bin)){
    for(k in 1:nrow(chr_dist)){
      if(isTRUE(chr_bin$foo.X1[i] >= chr_dist$Start[k] && chr_bin$foo.X2[i] <= chr_dist$End[k])){
        chr_bin$rate[i] <- (chr_bin$rate[i]+(chr_bin$rate[i]*chr_dist$diff2[k]))
      }else
        chr_bin$rate[i] <- chr_bin$rate[i]
    }
  }
  return(chr_bin)
}
chr1_bin2 <- zmet2_map(chr1_bin, chr1_dist)
chr2_bin2 <- zmet2_map(chr2_bin, chr2_dist)
chr3_bin2 <- zmet2_map(chr3_bin, chr3_dist)
chr4_bin2 <- zmet2_map(chr4_bin, chr4_dist)
chr5_bin2 <- zmet2_map(chr5_bin, chr5_dist)
chr6_bin2 <- zmet2_map(chr6_bin, chr6_dist)
chr7_bin2 <- zmet2_map(chr7_bin, chr7_dist)
chr8_bin2 <- zmet2_map(chr8_bin, chr8_dist)
chr9_bin2 <- zmet2_map(chr9_bin, chr9_dist)
chr10_bin2 <- zmet2_map(chr10_bin, chr10_dist)

#using function, converted SNP start to Mb to get cM/Mb for final genetic position
chr1_snp2 <- snp_rate(chr1_bin2, chr1_snp)
chr1_snp2$`SNP Start`<- chr1_snp2$`SNP Start`/1000000
chr1_snp2 <- chr1_snp2[order(chr1_snp2$`SNP Start`),]
#smoothing the recombination rate so transitions between bins are not so abrupt
chr1_spl <- smooth.spline(chr1_snp2$rate, spar = 1.1)
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

chr2_snp2 <- snp_rate(chr2_bin2, chr2_snp)
chr2_snp2$`SNP Start` <- chr2_snp2$`SNP Start`/1000000
chr2_snp2 <- chr2_snp2[-(228:237),]
chr2_spl <- smooth.spline(chr2_snp2$rate, spar = 1.2)
chr2_snp2$pos <- (chr2_snp2$`SNP Start`*chr2_spl$y)
plot(chr2_snp2$`SNP Start`, chr2_snp2$pos)
plot(chr2_snp2$`SNP Start`, chr2_snp2$pos/chr2_snp2$`SNP Start`, type = "l")
chr2_finalpos <- chr2_snp2[order(chr2_snp2$pos),]
is.unsorted(chr2_finalpos$pos)
plot(chr2_snp2$`SNP Start`, chr2_finalpos$pos/chr2_snp2$`SNP Start`, type = "l")

chr3_snp2 <- snp_rate(chr3_bin2, chr3_snp)
chr3_snp2$`SNP Start` <- chr3_snp2$`SNP Start`/1000000
chr3_spl <- smooth.spline(chr3_snp2$rate, spar = 1.2)
chr3_snp2$pos <- (chr3_snp2$`SNP Start`*chr3_spl$y)
plot(chr3_snp2$`SNP Start`, chr3_snp2$pos)
plot(chr3_snp2$`SNP Start`, chr3_snp2$pos/chr3_snp2$`SNP Start`, type = "l")
chr3_finalpos <- chr3_snp2[order(chr3_snp2$pos),]
is.unsorted(chr3_finalpos$pos)
plot(chr3_snp2$`SNP Start`, chr3_finalpos$pos/chr3_snp2$`SNP Start`, type = "l")

chr4_snp2 <- snp_rate(chr4_bin2, chr4_snp)
chr4_snp2$`SNP Start` <- chr4_snp2$`SNP Start`/1000000
chr4_spl <- smooth.spline(chr4_snp2$rate, spar = 1.2)
chr4_snp2$pos <- (chr4_snp2$`SNP Start`*chr4_spl$y)
plot(chr4_snp2$`SNP Start`, chr4_snp2$pos)
plot(chr4_snp2$`SNP Start`, chr4_snp2$pos/chr4_snp2$`SNP Start`, type = "l")
chr4_finalpos <- chr4_snp2[order(chr4_snp2$pos),]
is.unsorted(chr4_finalpos$pos)
plot(chr4_snp2$`SNP Start`, chr4_finalpos$pos/chr4_snp2$`SNP Start`, type = "l")

chr5_snp2 <- snp_rate(chr5_bin2, chr5_snp)
chr5_snp2$`SNP Start` <- chr5_snp2$`SNP Start`/1000000
chr5_spl <- smooth.spline(chr5_snp2$rate, spar = 1.1)
chr5_snp2$pos <- (chr5_snp2$`SNP Start`*chr5_spl$y)
plot(chr5_snp2$`SNP Start`, chr5_snp2$pos)
plot(chr5_snp2$`SNP Start`, chr5_snp2$pos/chr5_snp2$`SNP Start`, type = "l")
chr5_finalpos <- chr5_snp2[order(chr5_snp2$pos),]
is.unsorted(chr5_finalpos$pos)
plot(chr5_snp2$`SNP Start`, chr5_finalpos$pos/chr5_snp2$`SNP Start`, type = "l")

#chr 6 is lowkey fuked up
chr6_snp2 <- snp_rate(chr6_bin2, chr6_snp)
chr6_snp2$`SNP Start` <- chr6_snp2$`SNP Start`/1000000
chr6_spl <- smooth.spline(chr6_snp2$rate, spar = 1)
chr6_snp2$pos <- (chr6_snp2$`SNP Start`*chr6_spl$y)
plot(chr6_snp2$`SNP Start`, chr6_snp2$pos)
plot(chr6_snp2$`SNP Start`, chr6_snp2$pos/chr6_snp2$`SNP Start`, type = "l")
chr6_finalpos <- chr6_snp2[order(chr6_snp2$pos),]
is.unsorted(chr6_finalpos$pos)
plot(chr6_snp2$`SNP Start`, chr6_finalpos$pos/chr6_snp2$`SNP Start`, type = "l")

chr7_snp2 <- snp_rate(chr7_bin2, chr7_snp)
chr7_snp2$`SNP Start` <- chr7_snp2$`SNP Start`/1000000
chr7_spl <- smooth.spline(chr7_snp2$rate, spar = 1.15)
chr7_snp2$pos <- (chr7_snp2$`SNP Start`*chr7_spl$y)
plot(chr7_snp2$`SNP Start`, chr7_snp2$pos)
plot(chr7_snp2$`SNP Start`, chr7_snp2$pos/chr7_snp2$`SNP Start`, type = "l")
chr7_finalpos <- chr7_snp2[order(chr7_snp2$pos),]
is.unsorted(chr7_finalpos$pos)
plot(chr7_snp2$`SNP Start`, chr7_finalpos$pos/chr7_snp2$`SNP Start`, type = "l")

chr8_snp2 <- snp_rate(chr8_bin2, chr8_snp)
chr8_snp2$`SNP Start` <- chr8_snp2$`SNP Start`/1000000
chr8_spl <- smooth.spline(chr8_snp2$rate, spar = 1.15)
chr8_snp2$pos <- (chr8_snp2$`SNP Start`*chr8_spl$y)
plot(chr8_snp2$`SNP Start`, chr8_snp2$pos)
plot(chr8_snp2$`SNP Start`, chr8_snp2$pos/chr8_snp2$`SNP Start`, type = "l")
chr8_finalpos <- chr8_snp2[order(chr8_snp2$pos),]
is.unsorted(chr8_finalpos$pos)
plot(chr8_snp2$`SNP Start`, chr8_finalpos$pos/chr8_snp2$`SNP Start`, type = "l")

chr9_snp2 <- snp_rate(chr9_bin2, chr9_snp)
chr9_snp2$`SNP Start` <- chr9_snp2$`SNP Start`/1000000
chr9_spl <- smooth.spline(chr9_snp2$rate, spar = 1.1)
chr9_snp2$pos <- (chr9_snp2$`SNP Start`*chr9_spl$y)
plot(chr9_snp2$`SNP Start`, chr9_snp2$pos)
plot(chr9_snp2$`SNP Start`, chr9_snp2$pos/chr9_snp2$`SNP Start`, type = "l")
chr9_finalpos <- chr9_snp2[order(chr9_snp2$pos),]
is.unsorted(chr9_finalpos$pos)
plot(chr9_snp2$`SNP Start`, chr9_finalpos$pos/chr9_snp2$`SNP Start`, type = "l")

chr10_snp2 <- snp_rate(chr10_bin2, chr10_snp)
chr10_snp2$`SNP Start` <- chr10_snp2$`SNP Start`/1000000
chr10_spl <- smooth.spline(chr10_snp2$rate, spar = 1.15)
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

final_map <- list(chr1[[1]], chr2[[1]], 
                  chr3[[1]], chr4[[1]], chr5[[1]], 
                  chr6[[1]], chr7[[1]], chr8[[1]], 
                  chr9[[1]], chr10[[1]])

#Creating vector of centromere positions for zmet2
real_centromere <- c(68.7975, 83.0778, 45.34253, 42.34426, 49.09542, 
                     43.64207, 33.1046, 45.02515, 49.45323, 23.55771)
real_centromere <- real_centromere/100

###Simulating a realistic breeding program in maize
founderPop <- quickHaplo(nInd = 200, nChr = 10, segSites = c(nrow(chr1_finalpos), nrow(chr2_finalpos), 
                                                             nrow(chr3_finalpos), nrow(chr4_finalpos), nrow(chr5_finalpos), 
                                                             nrow(chr6_finalpos), nrow(chr7_finalpos), nrow(chr8_finalpos),
                                                             nrow(chr9_finalpos), nrow(chr10_finalpos)))

lociMap <- new("LociMap",nLoci= as.integer(1), lociPerChr=as.integer(1),lociLoc=as.integer(c(1,2,3,4)))
founderPop@genMap <- final_map
founderPop@centromere <- real_centromere
SP = SimParam$new(founderPop)
SP$setTrackRec(TRUE)
SP$v = 2.6
SP$p = 0.2
nQtlPerChr = 1
SP$addTraitAD(nQtlPerChr)
SP$setVarE(h2=0.5)

#Bad population first with low breeding values
pop_bad_sel10 <- vector(mode = "list", length = 20)
for(i in 1:20){
  SP$manAddTrait(lociMap)
  pop_bad <- newPop(founderPop, simParam = SP)
  pop_bad <- setPheno(pop_bad, h2 = 0.5, simParam = SP)
  
  pop_bad1 <- randCross(pop_bad, nCrosses = 100, nProgeny=100, simParam = SP)
  pop_bad1 <- setPheno(pop_bad1, h2 = 0.5, simParam = SP)
  
  pop_bad_sel <- selectInd(pop_bad1, nInd = 50, use = "bv", trait = 1, selectTop = FALSE, returnPop = TRUE, simParam = SP)
  pop_bad2 <- randCross(pop_bad_sel, nCrosses = 50, nProgeny = 100, simParam = SP)
  pop_bad2 <- setPheno(pop_bad2, h2 = 0.5, simParam = SP)
  
  pop_bad_sel2 <- selectInd(pop_bad2, nInd = 50, use = "bv", trait = 1, selectTop = FALSE, returnPop = TRUE, simParam = SP)
  pop_bad3 <- randCross(pop_bad_sel2, nCrosses = 50, nProgeny = 100, simParam = SP)
  pop_bad3 <- setPheno(pop_bad3, h2 = 0.5, simParam = SP)
  
  pop_bad_sel3 <- selectInd(pop_bad3, nInd = 50, use = "bv", trait = 1, selectTop = FALSE, returnPop = TRUE, simParam = SP)
  pop_bad4 <- randCross(pop_bad_sel2, nCrosses = 50, nProgeny = 100, simParam = SP)
  pop_bad4 <- setPheno(pop_bad4, h2 = 0.5, simParam = SP)
  
  pop_bad_sel4 <- selectInd(pop_bad4, nInd = 50, use = "bv", trait = 1, selectTop = FALSE, returnPop = TRUE, simParam = SP)
  pop_bad5 <- randCross(pop_bad_sel4, nCrosses = 50, nProgeny = 100, simParam = SP)
  pop_bad5 <- setPheno(pop_bad5, h2 = 0.5, simParam = SP)
  
  pop_bad_sel5 <- selectInd(pop_bad5, nInd = 50, use = "bv", trait = 1, selectTop = FALSE, returnPop = TRUE, simParam = SP)
  pop_bad6 <- randCross(pop_bad_sel5, nCrosses = 50, nProgeny = 100, simParam = SP)
  pop_bad6 <- setPheno(pop_bad6, h2 = 0.5, simParam = SP)
  
  pop_bad_sel6 <- selectInd(pop_bad6, nInd = 50, use = "bv", trait = 1, selectTop = FALSE, returnPop = TRUE, simParam = SP)
  pop_bad7 <- randCross(pop_bad_sel6, nCrosses = 50, nProgeny = 100, simParam = SP)
  pop_bad7 <- setPheno(pop_bad7, h2 = 0.5, simParam = SP)
  
  pop_bad_sel7 <- selectInd(pop_bad7, nInd = 50, use = "bv", trait = 1, selectTop = FALSE, returnPop = TRUE, simParam = SP)
  pop_bad8 <- randCross(pop_bad_sel7, nCrosses = 50, nProgeny = 100, simParam = SP)
  pop_bad8 <- setPheno(pop_bad8, h2 = 0.5, simParam = SP)
  
  pop_bad_sel8 <- selectInd(pop_bad8, nInd = 50, use = "bv", trait = 1, selectTop = FALSE, returnPop = TRUE, simParam = SP)
  pop_bad9 <- randCross(pop_bad_sel8, nCrosses = 50, nProgeny = 100, simParam = SP)
  pop_bad9 <- setPheno(pop_bad9, h2 = 0.5, simParam = SP)
  
  pop_bad_sel9 <- selectInd(pop_bad9, nInd = 50, use = "bv", trait = 1, selectTop = FALSE, returnPop = TRUE, simParam = SP)
  pop_bad10 <- randCross(pop_bad_sel9, nCrosses = 50, nProgeny = 100, simParam = SP)
  pop_bad10 <- setPheno(pop_bad10, h2 = 0.5, simParam = SP)
  
  pop_bad_sel10[[i]] <- selectInd(pop_bad10, nInd = 20, use = "bv", trait = 1, selectop =FALSE, returnPop = TRUE, simParam = SP)
}

#Creating the "good"/elite pop with high breeding values
pop_good_sel10 <- vector(mode = "list", length = 20)
for(i in 1:20){
  pop_good <- newPop(founderPop, simParam = SP)
  pop_good <- setPheno(pop_good, h2 = 0.5, simParam = SP)
  
  pop_good1 <- randCross(pop_good, nCrosses = 100, nProgeny=100, simParam = SP)
  pop_good1 <- setPheno(pop_good1, h2 = 0.5, simParam = SP)
  
  pop_good_sel <- selectInd(pop_good1, nInd = 50, use = "bv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good2 <- randCross(pop_good_sel, nCrosses = 50, nProgeny = 100, simParam = SP)
  pop_good2 <- setPheno(pop_good2, h2 = 0.5, simParam = SP)
  
  pop_good_sel2 <- selectInd(pop_good2, nInd = 50, use = "bv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good3 <- randCross(pop_good_sel2, nCrosses = 50, nProgeny = 100, simParam = SP)
  pop_good3 <- setPheno(pop_good3, h2 = 0.5, simParam = SP)
  
  pop_good_sel3 <- selectInd(pop_good3, nInd = 50, use = "bv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good4 <- randCross(pop_good_sel2, nCrosses = 50, nProgeny = 100, simParam = SP)
  pop_good4 <- setPheno(pop_good4, h2 = 0.5, simParam = SP)
  
  pop_good_sel4 <- selectInd(pop_good4, nInd = 50, use = "bv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good5 <- randCross(pop_good_sel4, nCrosses = 50, nProgeny = 100, simParam = SP)
  pop_good5 <- setPheno(pop_good5, h2 = 0.5, simParam = SP)
  
  pop_good_sel5 <- selectInd(pop_good5, nInd = 50, use = "bv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good6 <- randCross(pop_good_sel5, nCrosses = 50, nProgeny = 100, simParam = SP)
  pop_good6 <- setPheno(pop_good6, h2 = 0.5, simParam = SP)
  
  pop_good_sel6 <- selectInd(pop_good6, nInd = 50, use = "bv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good7 <- randCross(pop_good_sel6, nCrosses = 50, nProgeny = 100, simParam = SP)
  pop_good7 <- setPheno(pop_good7, h2 = 0.5, simParam = SP)
  
  pop_good_sel7 <- selectInd(pop_good7, nInd = 50, use = "bv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good8 <- randCross(pop_good_sel7, nCrosses = 50, nProgeny = 100, simParam = SP)
  pop_good8 <- setPheno(pop_good8, h2 = 0.5, simParam = SP)
  
  pop_good_sel8 <- selectInd(pop_good8, nInd = 50, use = "bv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good9 <- randCross(pop_good_sel8, nCrosses = 50, nProgeny = 100, simParam = SP)
  pop_good9 <- setPheno(pop_good9, h2 = 0.5, simParam = SP)
  
  pop_good_sel9 <- selectInd(pop_good9, nInd = 50, use = "bv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good10 <- randCross(pop_good_sel9, nCrosses = 50, nProgeny = 100, simParam = SP)
  pop_good10 <- setPheno(pop_good10, h2 = 0.5, simParam = SP)
  
  pop_good_sel10[[i]] <- selectInd(pop_good10, nInd = 20, use = "bv", trait = 1, selectop = TRUE, returnPop = TRUE, simParam = SP)
}

popList = list(pop_good_sel10, pop_bad_sel10)
mergedpops = mergePops(popList)

#put it together to iterate one program 40 times with 20 generations of selection
zmet2_pheno <- matrix(nrow=40,ncol=10)
for(i in 1:40){
  pop <- newPop(mergedpops, simParam = SP)
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
  zmet2[i,]<- pheno(final_sel)
}

#Creating confidence intervals
pop1_mean <- mean(pop1_gv)
pop1_sd <- sd(pop1_gv)
pop1_size <- founderpop@nInd
pop1_se <- pop1_sd/sqrt(pop1_size)
alpha = 0.01
degrees.freedom = pop1_size - 1
t.score = qt(p=alpha/2, df=degrees.freedom,lower.tail=F)
margin.error <- t.score * pop1_se
lower.bound <- pop1_mean - margin.error
upper.bound <- pop1_mean + margin.error
print(c(lower.bound,upper.bound))

#Plotting gv on histogram
hist(pop1_gv)

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