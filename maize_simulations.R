library(AlphaSimR)
library(Rcpp)
library(ggplot2)
library(dplyr)

setwd("/home/rke27/Documents")
set.seed(420)

##reading in SNPs from B73xMo17 based on v4 B73 ref
final_snps <- read.table("SNP_V4.bed", header = FALSE)
colnames(final_snps) <- c("Chr#", "SNP Start", "SNP End")
#2000 SNPs genome-wide
final_snps <- sample_n(final_snps, 2000)
final_snps <- final_snps[order(final_snps$`Chr#`,final_snps$`SNP Start`),]
#write.csv(final_snps, "C:/Users/16192/Documents/PNAS_Simulations/final_snps.csv", row.names = FALSE)

chr1_snp <- final_snps[ which(final_snps$`Chr#` == "chr1"),]
#hist(chr1_snp$`SNP Start`, breaks = 100, xlab = "Physical Pos (bp)", main = "Chr.1 SNP Density", col = "black")
chr1_snp$rate <- NA
#making SNPs start at 0
chr1_snp$`SNP End` <- chr1_snp$`SNP End` - min(chr1_snp$`SNP Start`)
chr1_snp$`SNP Start` <- chr1_snp$`SNP Start`- min(chr1_snp$`SNP Start`)

chr2_snp <- final_snps[ which(final_snps$`Chr#` == "chr2"),]
#hist(chr2_snp$`SNP Start`, breaks = 100)
chr2_snp$rate <- NA
chr2_snp$`SNP End` <- chr2_snp$`SNP End` - min(chr2_snp$`SNP Start`)
chr2_snp$`SNP Start` <- chr2_snp$`SNP Start`- min(chr2_snp$`SNP Start`)

chr3_snp <- final_snps[ which(final_snps$`Chr#` == "chr3"),]
#hist(chr3_snp$`SNP Start`, breaks = 100)
chr3_snp$rate <- NA
chr3_snp$`SNP End` <- chr3_snp$`SNP End` - min(chr3_snp$`SNP Start`)
chr3_snp$`SNP Start` <- chr3_snp$`SNP Start`- min(chr3_snp$`SNP Start`)

chr4_snp <- final_snps[ which(final_snps$`Chr#` == "chr4"),]
#hist(chr4_snp$`SNP Start`, breaks = 100)
chr4_snp$rate <- NA
chr4_snp$`SNP End` <- chr4_snp$`SNP End` - min(chr4_snp$`SNP Start`)
chr4_snp$`SNP Start` <- chr4_snp$`SNP Start`- min(chr4_snp$`SNP Start`)

chr5_snp <- final_snps[ which(final_snps$`Chr#` == "chr5"),]
#hist(chr5_snp$`SNP Start`, breaks = 100)
chr5_snp$rate <- NA
chr5_snp$`SNP End` <- chr5_snp$`SNP End` - min(chr5_snp$`SNP Start`)
chr5_snp$`SNP Start` <- chr5_snp$`SNP Start`- min(chr5_snp$`SNP Start`)

chr6_snp <- final_snps[ which(final_snps$`Chr#` == "chr6"),]
#hist(chr6_snp$`SNP Start`, breaks = 100)
chr6_snp$rate <- NA
chr6_snp$`SNP End` <- chr6_snp$`SNP End` - min(chr6_snp$`SNP Start`)
chr6_snp$`SNP Start` <- chr6_snp$`SNP Start`- min(chr6_snp$`SNP Start`)

chr7_snp <- final_snps[ which(final_snps$`Chr#` == "chr7"),]
#hist(chr7_snp$`SNP Start`, breaks = 100)
chr7_snp$rate <- NA
chr7_snp$`SNP End` <- chr7_snp$`SNP End` - min(chr7_snp$`SNP Start`)
chr7_snp$`SNP Start` <- chr7_snp$`SNP Start`- min(chr7_snp$`SNP Start`)

chr8_snp <- final_snps[ which(final_snps$`Chr#` == "chr8"),]
#hist(chr8_snp$`SNP Start`, breaks = 100)
chr8_snp$rate <- NA
chr8_snp$`SNP End` <- chr8_snp$`SNP End` - min(chr8_snp$`SNP Start`)
chr8_snp$`SNP Start` <- chr8_snp$`SNP Start`- min(chr8_snp$`SNP Start`)

chr9_snp <- final_snps[ which(final_snps$`Chr#` == "chr9"),]
#hist(chr9_snp$`SNP Start`, breaks = 100)
chr9_snp$rate <- NA
chr9_snp$`SNP End` <- chr9_snp$`SNP End` - min(chr9_snp$`SNP Start`)
chr9_snp$`SNP Start` <- chr9_snp$`SNP Start`- min(chr9_snp$`SNP Start`)

chr10_snp <- final_snps[ which(final_snps$`Chr#` == "chr10"),]
#hist(chr10_snp$`SNP Start`, breaks = 100)
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
#hist(chr1_CO$`CO Start`, breaks = 300, xlab = "Physical Pos (bp)", main = "WT Chr.1 CO Density", col = "black")

chr2_CO <- NAM[ which(NAM$Chr == 2),]
#hist(chr2_CO$`CO Start`, breaks = 300)
chr2_CO$midpoint <- (chr2_CO$`CO Start`+ chr2_CO$`CO End`)/2

chr3_CO <- NAM[ which(NAM$Chr == 3),]
#hist(chr3_CO$`CO Start`, breaks = 300)
chr3_CO$midpoint <- (chr3_CO$`CO Start`+ chr3_CO$`CO End`)/2

chr4_CO <- NAM[ which(NAM$Chr == 4),]
#hist(chr4_CO$`CO Start`, breaks = 300)
chr4_CO$midpoint <- (chr4_CO$`CO Start`+ chr4_CO$`CO End`)/2

chr5_CO <- NAM[ which(NAM$Chr == 5),]
#hist(chr5_CO$`CO Start`, breaks = 300)
chr5_CO$midpoint <- (chr5_CO$`CO Start`+ chr5_CO$`CO End`)/2

chr6_CO <- NAM[ which(NAM$Chr == 6),]
#hist(chr6_CO$`CO Start`, breaks = 300)
chr6_CO$midpoint <- (chr6_CO$`CO Start`+ chr6_CO$`CO End`)/2

chr7_CO <- NAM[ which(NAM$Chr == 7),]
#hist(chr7_CO$`CO Start`, breaks = 300)
chr7_CO$midpoint <- (chr7_CO$`CO Start`+ chr7_CO$`CO End`)/2

chr8_CO <- NAM[ which(NAM$Chr == 8),]
#hist(chr8_CO$`CO Start`, breaks = 300)
chr8_CO$midpoint <- (chr8_CO$`CO Start`+ chr8_CO$`CO End`)/2

chr9_CO <- NAM[ which(NAM$Chr == 9),]
#hist(chr9_CO$`CO Start`, breaks = 300)
chr9_CO$midpoint <- (chr9_CO$`CO Start`+ chr9_CO$`CO End`)/2

chr10_CO <- NAM[ which(NAM$Chr == 10),]
#hist(chr10_CO$`CO Start`, breaks = 300)
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

#bin crossovers into ~1Mb uneven bins
chr1_CO <- chr1_CO[order(chr1_CO$`CO Start`),]
chr1_bin <- binning(chr1_CO$midpoint, nbins = max(chr1_snp$`SNP End`)/1000000, type = "kmeans")
chr1_bin <- as.data.frame(summary(chr1_bin))
#chr1_bin$freq <- chr1_bin$freq*2/4713
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
chr1_bin[max(chr1_snp$`SNP End`)/1000000,5] <- max(chr1_snp$`SNP End`)
#adding length of bin as column and making in Mb
chr1_bin$length <- (chr1_bin$foo.X2-chr1_bin$foo.X1)/1000000
chr1_bin$rate <- ((chr1_bin$freq/4713)*100)/chr1_bin$length

chr2_CO <- chr2_CO[order(chr2_CO$`CO Start`),]
chr2_bin <- binning(chr2_CO$midpoint, nbins = round(max(chr2_snp$`SNP End`)/1000000), type = "kmeans")
chr2_bin <- as.data.frame(summary(chr2_bin))
#chr2_bin$freq <- chr2_bin$freq*2/4713
chr2_bin <- within(chr2_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr2_bin$levels), ',', fixed=TRUE))))
chr2_bin <- do.call(data.frame, chr2_bin)
chr2_bin <- chr2_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr2_bin <- chr2_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr2_bin[1,4] <- 440104
chr2_bin$foo.X1 <- chr2_bin$foo.X1 - 440104
chr2_bin$foo.X2 <- chr2_bin$foo.X2 - 440104
chr2_bin[round(max(chr2_snp$`SNP End`)/1000000),5] <- max(chr2_snp$`SNP End`)
chr2_bin$length <- (chr2_bin$foo.X2-chr2_bin$foo.X1)/1000000
chr2_bin$rate <- ((chr2_bin$freq/4713)*100)/chr2_bin$length

#have not converted rest of chromosomes to what I did in 1 & 2
chr3_CO <- chr3_CO[order(chr3_CO$`CO Start`),]
chr3_bin <- binning(chr3_CO$midpoint, nbins = max(chr3_snp$`SNP End`)/1000000, type = "kmeans")
chr3_bin <- as.data.frame(summary(chr3_bin))
#chr3_bin$freq <- chr3_bin$freq*2/4713
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
#chr4_bin$freq <- chr4_bin$freq*2/4713
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
#chr5_bin$freq <- chr5_bin$freq*2/4713
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
#chr6_bin$freq <- chr6_bin$freq*2/4713
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
#chr7_bin$freq <- chr7_bin$freq*2/4713
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
#chr8_bin$freq <- chr8_bin$freq*2/4713
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
#chr9_bin$freq <- chr9_bin$freq*2/4713
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
#chr10_bin$freq <- chr10_bin$freq*2/4713
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

#using function, converted SNP start to Mb to get cM/Mb for final genetic position
chr1_snp2 <- snp_rate(chr1_bin, chr1_snp)
chr1_snp2$`SNP Start`<- chr1_snp2$`SNP Start`/1000000
chr1_snp2 <- chr1_snp2[order(chr1_snp2$`SNP Start`),]
#chr1_snp2 <- chr1_snp2[-c(278:300),]
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


chr2_snp2 <- snp_rate(chr2_bin, chr2_snp)
chr2_snp2$`SNP Start` <- chr2_snp2$`SNP Start`/1000000
chr2_snp2 <- chr2_snp2[-(228:237),]
chr2_spl <- smooth.spline(chr2_snp2$rate, spar = 1.2)
chr2_snp2$pos <- (chr2_snp2$`SNP Start`*chr2_spl$y)
plot(chr2_snp2$`SNP Start`, chr2_snp2$pos)
plot(chr2_snp2$`SNP Start`, chr2_snp2$pos/chr2_snp2$`SNP Start`, type = "l")
chr2_finalpos <- chr2_snp2[order(chr2_snp2$pos),]
is.unsorted(chr2_finalpos$pos)
plot(chr2_snp2$`SNP Start`, chr2_finalpos$pos/chr2_snp2$`SNP Start`, type = "l")

chr3_snp2 <- snp_rate(chr3_bin, chr3_snp)
chr3_snp2$`SNP Start` <- chr3_snp2$`SNP Start`/1000000
chr3_spl <- smooth.spline(chr3_snp2$rate, spar = 1.2)
chr3_snp2$pos <- (chr3_snp2$`SNP Start`*chr3_spl$y)
plot(chr3_snp2$`SNP Start`, chr3_snp2$pos)
plot(chr3_snp2$`SNP Start`, chr3_snp2$pos/chr3_snp2$`SNP Start`, type = "l")
chr3_finalpos <- chr3_snp2[order(chr3_snp2$pos),]
is.unsorted(chr3_finalpos$pos)
plot(chr3_snp2$`SNP Start`, chr3_finalpos$pos/chr3_snp2$`SNP Start`, type = "l")

chr4_snp2 <- snp_rate(chr4_bin, chr4_snp)
chr4_snp2$`SNP Start` <- chr4_snp2$`SNP Start`/1000000
chr4_spl <- smooth.spline(chr4_snp2$rate, spar = 1.15)
chr4_snp2$pos <- (chr4_snp2$`SNP Start`*chr4_spl$y)
plot(chr4_snp2$`SNP Start`, chr4_snp2$pos)
plot(chr4_snp2$`SNP Start`, chr4_snp2$pos/chr4_snp2$`SNP Start`, type = "l")
chr4_finalpos <- chr4_snp2[order(chr4_snp2$pos),]
is.unsorted(chr4_finalpos$pos)
plot(chr4_snp2$`SNP Start`, chr4_finalpos$pos/chr4_snp2$`SNP Start`, type = "l")

chr5_snp2 <- snp_rate(chr5_bin, chr5_snp)
chr5_snp2$`SNP Start` <- chr5_snp2$`SNP Start`/1000000
chr5_spl <- smooth.spline(chr5_snp2$rate, spar = 1.2)
chr5_snp2$pos <- (chr5_snp2$`SNP Start`*chr5_spl$y)
plot(chr5_snp2$`SNP Start`, chr5_snp2$pos)
plot(chr5_snp2$`SNP Start`, chr5_snp2$pos/chr5_snp2$`SNP Start`, type = "l")
chr5_finalpos <- chr5_snp2[order(chr5_snp2$pos),]
is.unsorted(chr5_finalpos$pos)
plot(chr5_snp2$`SNP Start`, chr5_finalpos$pos/chr5_snp2$`SNP Start`, type = "l")

#chr 6 is lowkey fuked up
chr6_snp2 <- snp_rate(chr6_bin, chr6_snp)
chr6_snp2$`SNP Start` <- chr6_snp2$`SNP Start`/1000000
chr6_spl <- smooth.spline(chr6_snp2$rate, spar = 1)
chr6_snp2$pos <- (chr6_snp2$`SNP Start`*chr6_spl$y)
plot(chr6_snp2$`SNP Start`, chr6_snp2$pos)
plot(chr6_snp2$`SNP Start`, chr6_snp2$pos/chr6_snp2$`SNP Start`, type = "l")
chr6_finalpos <- chr6_snp2[order(chr6_snp2$pos),]
is.unsorted(chr6_finalpos$pos)
plot(chr6_snp2$`SNP Start`, chr6_finalpos$pos/chr6_snp2$`SNP Start`, type = "l")

chr7_snp2 <- snp_rate(chr7_bin, chr7_snp)
chr7_snp2$`SNP Start` <- chr7_snp2$`SNP Start`/1000000
chr7_spl <- smooth.spline(chr7_snp2$rate, spar = 1.15)
chr7_snp2$pos <- (chr7_snp2$`SNP Start`*chr7_spl$y)
plot(chr7_snp2$`SNP Start`, chr7_snp2$pos)
plot(chr7_snp2$`SNP Start`, chr7_snp2$pos/chr7_snp2$`SNP Start`, type = "l")
chr7_finalpos <- chr7_snp2[order(chr7_snp2$pos),]
is.unsorted(chr7_finalpos$pos)
plot(chr7_snp2$`SNP Start`, chr7_finalpos$pos/chr7_snp2$`SNP Start`, type = "l")

chr8_snp2 <- snp_rate(chr8_bin, chr8_snp)
chr8_snp2$`SNP Start` <- chr8_snp2$`SNP Start`/1000000
chr8_spl <- smooth.spline(chr8_snp2$rate, spar = 1.2)
chr8_snp2$pos <- (chr8_snp2$`SNP Start`*chr8_spl$y)
plot(chr8_snp2$`SNP Start`, chr8_snp2$pos)
plot(chr8_snp2$`SNP Start`, chr8_snp2$pos/chr8_snp2$`SNP Start`, type = "l")
chr8_finalpos <- chr8_snp2[order(chr8_snp2$pos),]
is.unsorted(chr8_finalpos$pos)
plot(chr8_snp2$`SNP Start`, chr8_finalpos$pos/chr8_snp2$`SNP Start`, type = "l")

chr9_snp2 <- snp_rate(chr9_bin, chr9_snp)
chr9_snp2$`SNP Start` <- chr9_snp2$`SNP Start`/1000000
chr9_spl <- smooth.spline(chr9_snp2$rate, spar = 1.1)
chr9_snp2$pos <- (chr9_snp2$`SNP Start`*chr9_spl$y)
plot(chr9_snp2$`SNP Start`, chr9_snp2$pos)
plot(chr9_snp2$`SNP Start`, chr9_snp2$pos/chr9_snp2$`SNP Start`, type = "l")
chr9_finalpos <- chr9_snp2[order(chr9_snp2$pos),]
is.unsorted(chr9_finalpos$pos)
plot(chr9_snp2$`SNP Start`, chr9_finalpos$pos/chr9_snp2$`SNP Start`, type = "l")

chr10_snp2 <- snp_rate(chr10_bin, chr10_snp)
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

final_map <- list(chr1[[1]], chr2[[1]], 
                  chr3[[1]], chr4[[1]], chr5[[1]], 
                  chr6[[1]], chr7[[1]], chr8[[1]], 
                  chr9[[1]], chr10[[1]])



#Creating vector of centromere positions
real_centromere <- c(69.19425, 29.31842, 43.34131, 42.1754, 47.11507,
                     43.40838, 30.01601, 43.69602, 47.12911, 21.17685)
real_centromere <- real_centromere/100

###Simulating a realistic breeding program in maize

##Burning in populations before use
#One population is "good", other is "bad" GVs
#Use 10 generations of random crossing & selecting best/worst to diverge them

#Creating the "good"/elite pop
pop_good_sel10 <- vector(mode = "list", length = 100)
addEff_coupling <- runif(212, min = 0, max = 0.1)
addEff_repulsion <- runif(209, min = -0.1, max = 0.1)
for(i in 1:100){
  founderPop <- quickHaplo(nInd = 200, nChr = 10, inbred = TRUE, ploidy = 2L, segSites = c(nrow(chr1_finalpos), nrow(chr2_finalpos), 
                                                                                           nrow(chr3_finalpos), nrow(chr4_finalpos), nrow(chr5_finalpos), 
                                                                                           nrow(chr6_finalpos), nrow(chr7_finalpos), nrow(chr8_finalpos),
                                                                                           nrow(chr9_finalpos), nrow(chr10_finalpos)))
  founderPop@genMap <- final_map
  founderPop@centromere <- real_centromere
  SP = SimParam$new(founderPop)
  SP$setTrackRec(TRUE)
  SP$p = 0.15
  #trait_yield <- new("TraitA", nLoci = 209L, lociPerChr= c(32L, 24L, 24L, 24L, 22L, 16L, 18L, 18L, 17L, 14L), 
   #                  lociLoc = c(as.integer(runif(32L, min = 4L, max = 300L)), as.integer(runif(24L, min = 1L, max = 227L)), as.integer(runif(24L, min = 1L, max = 219L)),
    #                             as.integer(runif(24L, min = 1L, max = 256L)), as.integer(runif(22L, min = 1L, max = 146L)), as.integer(runif(16L, min = 1L, max = 92L)),
     #                            as.integer(runif(18L, min = 1L, max = 233L)), as.integer(runif(18L, min = 1L, max = 173L)), as.integer(runif(17L, min = 1L, max = 149L)),
      #                           as.integer(runif(14L, min = 1L, max = 184L))), addEff = addEff, intercept = 0.1)
  trait_yield <- new("TraitA", nLoci = 209L, lociPerChr= c(35L, 24L, 24L, 21L, 22L, 16L, 18L, 18L, 17L, 14L),
                                 lociLoc = c(1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,13L,14L,15L,16L,17L,18L,19L,20L,45L,46L,47L,
                                            48L,49L,50L,51L,52L,53L,54L,55L,56L,58L,59L,60L,1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,
                                           35L,36L,37L,38L,39L,40L,41L,42L,43L,44L,45L,46L, 1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,
                                          34L,35L,36L,37L,38L,39L,40L,41L,42L,43L,44L,45L, 1L,2L,3L,4L,5L,9L,10L,11L,12L,
                                         36L,37L,38L,39L,40L,41L,42L,43L,44L,45L,46L,47L, 1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,
                                        32L,33L,34L,35L,36L,37L,38L,39L,40L,41L,42L, 1L,2L,3L,4L,5L,6L,7L,8L,
                                       25L,26L,27L,28L,29L,30L,31L,32L, 1L,2L,3L,4L,5L,6L,7L,8L,9L, 26L,27L,28L,29L,30L,31L,32L,33L,
                                      34L, 1L,2L,3L,4L,5L,6L,7L,8L,9L, 26L,27L,28L,29L,30L,31L,32L,33L,
                                     34L, 1L,2L,3L,4L,5L,6L,7L,8L,9L, 23L,24L,25L,26L,27L,28L,29L,30L, 
                                    1L,2L,3L,4L,5L,6L,7L, 22L,23L,24L,25L,26L,27L,28L), addEff = addEff_repulsion, intercept = 0.1)
  #SP$addTraitA(nQtlPerChr = 21)
  SP$manAddTrait(trait_yield)
  pop_good <- newPop(founderPop, simParam = SP)
  pop_good <- setPheno(pop_good, h2 = 0.8, simParam = SP)
  
  pop_good1 <- randCross(pop_good, nCrosses = 10, nProgeny=10, simParam = SP)
  pop_good1 <- setPheno(pop_good1, h2 = 0.8, simParam = SP)
  
  pop_good_sel <- selectInd(pop_good1, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good2 <- randCross(pop_good_sel, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good2 <- setPheno(pop_good2, h2 = 0.8, simParam = SP)
  
  pop_good_sel2 <- selectInd(pop_good2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good3 <- randCross(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good3 <- setPheno(pop_good3, h2 = 0.8, simParam = SP)
  
  pop_good_sel3 <- selectInd(pop_good3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good4 <- randCross(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good4 <- setPheno(pop_good4, h2 = 0.8, simParam = SP)
  
  pop_good_sel4 <- selectInd(pop_good4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good5 <- randCross(pop_good_sel4, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good5 <- setPheno(pop_good5, h2 = 0.8, simParam = SP)
  
  pop_good_sel5 <- selectInd(pop_good5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good6 <- randCross(pop_good_sel5, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good6 <- setPheno(pop_good6, h2 = 0.8, simParam = SP)
  
  pop_good_sel6 <- selectInd(pop_good6, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good7 <- randCross(pop_good_sel6, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good7 <- setPheno(pop_good7, h2 = 0.8, simParam = SP)
  
  pop_good_sel7 <- selectInd(pop_good7, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good8 <- randCross(pop_good_sel7, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good8 <- setPheno(pop_good8, h2 = 0.8, simParam = SP)
  
  pop_good_sel8 <- selectInd(pop_good8, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good9 <- randCross(pop_good_sel8, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good9 <- setPheno(pop_good9, h2 = 0.8, simParam = SP)
  
  pop_good_sel9 <- selectInd(pop_good9, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good10 <- randCross(pop_good_sel9, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good10 <- setPheno(pop_good10, h2 = 0.8, simParam = SP)
  
  pop_good_sel10[[i]] <- selectInd(pop_good10, nInd = 2, use = "gv", trait = 1, selectop = TRUE, returnPop = TRUE, simParam = SP)
}

#trait_yield <- new("TraitAD", nLoci = 207L, lociPerChr= c(30L, 24L, 24L, 24L, 22L, 16L, 18L, 18L, 17L, 14L),
#            lociLoc = c(1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,13L,14L,15L,45L,46L,47L,
#                       48L,49L,50L,51L,52L,53L,54L,55L,56L,58L,59L,60L,1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,
#                      35L,36L,37L,38L,39L,40L,41L,42L,43L,44L,45L,46L, 1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,
#                     34L,35L,36L,37L,38L,39L,40L,41L,42L,43L,44L,45L, 1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,
#                    36L,37L,38L,39L,40L,41L,42L,43L,44L,45L,46L,47L, 1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,
#                   32L,33L,34L,35L,36L,37L,38L,39L,40L,41L,42L, 1L,2L,3L,4L,5L,6L,7L,8L,
#                  25L,26L,27L,28L,29L,30L,31L,32L, 1L,2L,3L,4L,5L,6L,7L,8L,9L, 26L,27L,28L,29L,30L,31L,32L,33L,
#                 34L, 1L,2L,3L,4L,5L,6L,7L,8L,9L, 26L,27L,28L,29L,30L,31L,32L,33L,
#                34L, 1L,2L,3L,4L,5L,6L,7L,8L,9L, 23L,24L,25L,26L,27L,28L,29L,30L, 
#               1L,2L,3L,4L,5L,6L,7L, 22L,23L,24L,25L,26L,27L,28L), addEff = rep(0.5, 207), domEff = rep(0.5, 207), intercept = 0)

#Bad population
pop_bad_sel10 <- vector(mode = "list", length = 100)
for(i in 1:100){
  founderPop <- quickHaplo(nInd = 200, nChr = 10, inbred = TRUE, 
                           ploidy = 2L, segSites = c(nrow(chr1_finalpos), 
                                                     nrow(chr2_finalpos), nrow(chr3_finalpos), 
                                                     nrow(chr4_finalpos), nrow(chr5_finalpos), 
                                                     nrow(chr6_finalpos), nrow(chr7_finalpos), 
                                                     nrow(chr8_finalpos),nrow(chr9_finalpos), 
                                                     nrow(chr10_finalpos)))
  founderPop@genMap <- final_map
  founderPop@centromere <- real_centromere
  SP = SimParam$new(founderPop)
  SP$setTrackRec(TRUE)
  SP$p = 0.15
  #SP$addTraitA(nQtlPerChr = 21)
  #trait_yield <- new("TraitA", nLoci = 209L, lociPerChr= c(32L, 24L, 24L, 24L, 22L, 16L, 18L, 18L, 17L, 14L), 
   #                  lociLoc = c(as.integer(runif(32L, min = 4L, max = 300L)), as.integer(runif(24L, min = 1L, max = 227L)), as.integer(runif(24L, min = 1L, max = 219L)),
    #                             as.integer(runif(24L, min = 1L, max = 256L)), as.integer(runif(22L, min = 1L, max = 146L)), as.integer(runif(16L, min = 1L, max = 92L)),
     #                            as.integer(runif(18L, min = 1L, max = 233L)), as.integer(runif(18L, min = 1L, max = 173L)), as.integer(runif(17L, min = 1L, max = 149L)),
      #                           as.integer(runif(14L, min = 1L, max = 184L))), addEff = addEff, intercept = 0.1)
  trait_yield <- new("TraitA", nLoci = 209L, lociPerChr= c(35L, 24L, 24L, 21L, 22L, 16L, 18L, 18L, 17L, 14L),
                     lociLoc = c(1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,13L,14L,15L,16L,17L,18L,19L,20L,45L,46L,47L,
                                  48L,49L,50L,51L,52L,53L,54L,55L,56L,58L,59L,60L,1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,
                                 35L,36L,37L,38L,39L,40L,41L,42L,43L,44L,45L,46L, 1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,
                                 34L,35L,36L,37L,38L,39L,40L,41L,42L,43L,44L,45L, 1L,2L,3L,4L,5L,9L,10L,11L,12L,
                                 36L,37L,38L,39L,40L,41L,42L,43L,44L,45L,46L,47L, 1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,
                                 32L,33L,34L,35L,36L,37L,38L,39L,40L,41L,42L, 1L,2L,3L,4L,5L,6L,7L,8L,
                                 25L,26L,27L,28L,29L,30L,31L,32L, 1L,2L,3L,4L,5L,6L,7L,8L,9L, 26L,27L,28L,29L,30L,31L,32L,33L,
                                 34L, 1L,2L,3L,4L,5L,6L,7L,8L,9L, 26L,27L,28L,29L,30L,31L,32L,33L,
                                 34L, 1L,2L,3L,4L,5L,6L,7L,8L,9L, 23L,24L,25L,26L,27L,28L,29L,30L, 
                                 1L,2L,3L,4L,5L,6L,7L, 22L,23L,24L,25L,26L,27L,28L), addEff = addEff_repulsion, intercept = 0.1)
  SP$manAddTrait(trait_yield)
  #additive effects of resistance loci add to 1 to get maintain variance
  trait <- new("TraitA", nLoci = 3L, lociPerChr = c(0L, 0L, 0L, 3L, 0L, 0L, 0L, 0L, 0L, 0L),
               lociLoc = c(6L,7L,8L), addEff = c(0.33,0.33,0.33), intercept = 0)
  SP$manAddTrait(trait)
  
  pop_bad <- newPop(founderPop, simParam = SP)
  pop_bad <- setPheno(pop_bad, h2 = 0.8, simParam = SP)
  
  pop_bad1 <- randCross(pop_bad, nCrosses = 10, nProgeny=10, simParam = SP)
  pop_bad1 <- setPheno(pop_bad1, h2 = 0.8, simParam = SP)
  
  #pop_bad_sel <- selectInd(pop_bad1, nInd = 10, use = "gv", trait = 1, selectTop = FALSE, returnPop = TRUE, simParam = SP)
  pop_bad2 <- selectCross(pop_bad1, nInd = 5, use = "gv", trait = 2, selectTop = TRUE, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_bad2 <- setPheno(pop_bad2, h2 = 0.8, simParam = SP)
  
  #pop_bad_sel2 <- selectInd(pop_bad2, nInd = 10, use = "gv", trait = 1, selectTop = FALSE, returnPop = TRUE, simParam = SP)
  pop_bad3 <- selectCross(pop_bad2, nInd = 5, use = "gv", trait = 2, selectTop = TRUE, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_bad3 <- setPheno(pop_bad3, h2 = 0.8, simParam = SP)
  
  #pop_bad_sel3 <- selectInd(pop_bad3, nInd = 10, use = "gv", trait = 1, selectTop = FALSE, returnPop = TRUE, simParam = SP)
  pop_bad4 <- selectCross(pop_bad3, nInd = 5, use = "gv", trait = 2, selectTop = TRUE, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_bad4 <- setPheno(pop_bad4, h2 = 0.8, simParam = SP)
  
  #pop_bad_sel4 <- selectInd(pop_bad4, nInd = 10, use = "gv", trait = 1, selectTop = FALSE, returnPop = TRUE, simParam = SP)
  pop_bad5 <- selectCross(pop_bad4, nInd = 5, use = "gv", trait = 2, selectTop = TRUE, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_bad5 <- setPheno(pop_bad5, h2 = 0.8, simParam = SP)
  
  #pop_bad_sel5 <- selectInd(pop_bad5, nInd = 10, use = "gv", trait = 1, selectTop = FALSE, returnPop = TRUE, simParam = SP)
  pop_bad6 <- selectCross(pop_bad5, nInd = 5, use = "gv", trait = 2, selectTop = TRUE, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_bad6 <- setPheno(pop_bad6, h2 = 0.8, simParam = SP)
  
  #pop_bad_sel6 <- selectInd(pop_bad6, nInd = 10, use = "gv", trait = 1, selectTop = FALSE, returnPop = TRUE, simParam = SP)
  pop_bad7 <- selectCross(pop_bad6, nInd = 5, use = "gv", trait = 2, selectTop = TRUE, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_bad7 <- setPheno(pop_bad7, h2 = 0.8, simParam = SP)
  
  #pop_bad_sel7 <- selectInd(pop_bad7, nInd = 10, use = "gv", trait = 1, selectTop = FALSE, returnPop = TRUE, simParam = SP)
  pop_bad8 <- selectCross(pop_bad7, nInd = 5, use = "gv", trait = 2, selectTop = TRUE, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_bad8 <- setPheno(pop_bad8, h2 = 0.8, simParam = SP)
  
  #pop_bad_sel8 <- selectInd(pop_bad8, nInd = 10, use = "gv", trait = 1, selectTop = FALSE, returnPop = TRUE, simParam = SP)
  pop_bad9 <- selectCross(pop_bad8, nInd = 5, use = "gv", trait = 2, selectTop = TRUE, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_bad9 <- setPheno(pop_bad9, h2 = 0.8, simParam = SP)
  
  #pop_bad_sel9 <- selectInd(pop_bad9, nInd = 10, use = "gv", trait = 1, selectTop = FALSE, returnPop = TRUE, simParam = SP)
  pop_bad10 <- selectCross(pop_bad9, nInd = 5, use = "gv", trait = 2, selectTop = TRUE, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_bad10 <- setPheno(pop_bad10, h2 = 0.8, simParam = SP)
  
  #pop_bad_sel10_1 <- selectInd(pop_bad10, nInd = 5, use = "gv", trait = 1, selecTop =FALSE, returnPop = TRUE, simParam = SP)
  pop_bad_sel10[[i]] <- selectInd(pop_bad10, nInd = 2, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
}

goodpop = mergePops(pop_good_sel10)
badpop = mergePops(pop_bad_sel10)

#use genParam to track variance of all types

pop1_sel2_gv <- matrix(data = NA, nrow = 100, ncol = 20)
pop1_sel3_gv <- matrix(data = NA, nrow = 100, ncol = 20)
pop1_sel4_gv <- matrix(data = NA, nrow = 100, ncol = 20)
pop1_sel5_gv <- matrix(data = NA, nrow = 100, ncol = 20)
pop1_sel6_gv <- matrix(data = NA, nrow = 100, ncol = 20)
pop1_sel7_gv <- matrix(data = NA, nrow = 100, ncol = 20)
wtvar1 <- matrix(data = NA, ncol = 2, nrow = 200)
wtvar2 <- matrix(data = NA, ncol = 2, nrow = 200)
wtvar3 <- matrix(data = NA, ncol = 2, nrow = 200)
wtvar4 <- matrix(data = NA, ncol = 2, nrow = 200)
wtvar5 <- matrix(data = NA, ncol = 2, nrow = 200)
wtvar6 <- matrix(data = NA, ncol = 2, nrow = 200)
wtvar7 <- matrix(data = NA, ncol = 2, nrow = 200)
for(i in 1:100){
  SP$switchGenMap(final_map, centromere = real_centromere)
  pop <- randCross2(goodpop[49,], badpop[155,], nCrosses = 20, nProgeny = 10, simParam = SP)
  pop <- setPheno(pop = pop, h2 = 0.8, simParam = SP)
  wtvar1[i] <- varA(pop)
  
  pop1_sel <- selectInd(pop, nInd = 80, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_2 <- selectInd(pop1_sel, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_cross <- randCross2(pop1_sel1_2, goodpop[49,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel1_cross <- setPheno(pop1_sel1_cross, h2 = 0.8, simParam = SP)
  wtvar2[i] <- varA(pop1_sel1_2)
  
  pop1_sel2_gv[i,] <- gv(pop1_sel1_2)
  pop1_sel2 <- selectInd(pop1_sel1_cross, nInd = 20, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_2 <- selectInd(pop1_sel2, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_cross <- randCross2(pop1_sel2_2, goodpop[49,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel2_cross <- setPheno(pop1_sel2_cross, h2 = 0.8, simParam = SP)
  wtvar3[i] <- varA(pop1_sel2_2)
  
  pop1_sel3_gv[i,] <- gv(pop1_sel2_2)
  pop1_sel3 <- selectInd(pop1_sel2_cross, nInd = 20, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_2 <- selectInd(pop1_sel3, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_cross <- randCross2(pop1_sel3_2, goodpop[49,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel3_cross <- setPheno(pop1_sel3_cross, h2 = 0.8, simParam = SP)
  wtvar4[i] <- varA(pop1_sel3_2)
  
  pop1_sel4_gv[i,] <- gv(pop1_sel3_2)
  pop1_sel4 <- selectInd(pop1_sel3_cross, nInd = 20, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_2 <- selectInd(pop1_sel4, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_cross <- randCross2(pop1_sel4_2, goodpop[49,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel4_cross <- setPheno(pop1_sel4_cross, h2 = 0.8, simParam = SP)
  wtvar5[i] <- varA(pop1_sel4_2)
  
  pop1_sel5_gv[i,] <- gv(pop1_sel4_2)
  pop1_sel5 <- selectInd(pop1_sel4_cross, nInd = 20, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_2 <- selectInd(pop1_sel5, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_cross <- randCross2(pop1_sel5_2, goodpop[49,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel5_cross <- setPheno(pop1_sel5_cross, h2 = 0.8, simParam = SP)
  wtvar6[i] <- varA(pop1_sel5_2)
  
  pop1_sel6_gv[i,] <- gv(pop1_sel5_2)
}

##Starting intermating schemes for different genetic maps

pop_good_selwt <- vector(mode = "list", length = 100)
wtvar1 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar2 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar3 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar4 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar5 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar6 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar7 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar8 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar9 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar10 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar11 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar12 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar13 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar14 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar15 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar16 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar17 <- matrix(data = NA, ncol = 1, nrow = 100)
wtgv1 <- matrix(data = NA, ncol = 5, nrow = 100)
wtgv2 <- matrix(data = NA, ncol = 5, nrow = 100)
wtgv3 <- matrix(data = NA, ncol = 5, nrow = 100)
wtgv4 <- matrix(data = NA, ncol = 5, nrow = 100)
wtgv5 <- matrix(data = NA, ncol = 5, nrow = 100)
wtgv6 <- matrix(data = NA, ncol = 5, nrow = 100)
wtgv7 <- matrix(data = NA, ncol = 5, nrow = 100)
wtgv8 <- matrix(data = NA, ncol = 5, nrow = 100)
wtgv9 <- matrix(data = NA, ncol = 5, nrow = 100)
wtgv10 <- matrix(data = NA, ncol = 5, nrow = 100)
wtgv11 <- matrix(data = NA, ncol = 5, nrow = 100)
wtgv12 <- matrix(data = NA, ncol = 5, nrow = 100)
wtgv13 <- matrix(data = NA, ncol = 5, nrow = 100)
wtgv14 <- matrix(data = NA, ncol = 5, nrow = 100)
wtgv15 <- matrix(data = NA, ncol = 5, nrow = 100)
#so all simulations have same marker effects
addEff <- runif(207, min= -0.1, max = 0.1)
for(i in 1:100){
  founderPop <- quickHaplo(nInd = 200, nChr = 10, inbred = TRUE, ploidy = 2L, segSites = c(nrow(chr1_finalpos), nrow(chr2_finalpos), 
                                                                                           nrow(chr3_finalpos), nrow(chr4_finalpos), nrow(chr5_finalpos), 
                                                                                           nrow(chr6_finalpos), nrow(chr7_finalpos), nrow(chr8_finalpos),
                                                                                           nrow(chr9_finalpos), nrow(chr10_finalpos)))
  founderPop@genMap <- final_map
  founderPop@centromere <- real_centromere
  SP = SimParam$new(founderPop)
  SP$setTrackRec(TRUE)
  SP$p = 0.15
  #SP$addTraitADE(nQtlPerChr = 92)
  trait_yield <- new("TraitAD", nLoci = 207L, lociPerChr= c(30L, 24L, 24L, 24L, 22L, 16L, 18L, 18L, 17L, 14L),
              lociLoc = c(1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,13L,14L,15L,45L,46L,47L,
                         48L,49L,50L,51L,52L,53L,54L,55L,56L,58L,59L,60L,1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,
                        35L,36L,37L,38L,39L,40L,41L,42L,43L,44L,45L,46L, 1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,
                       34L,35L,36L,37L,38L,39L,40L,41L,42L,43L,44L,45L, 1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,
                      36L,37L,38L,39L,40L,41L,42L,43L,44L,45L,46L,47L, 1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,
                     32L,33L,34L,35L,36L,37L,38L,39L,40L,41L,42L, 1L,2L,3L,4L,5L,6L,7L,8L,
                    25L,26L,27L,28L,29L,30L,31L,32L, 1L,2L,3L,4L,5L,6L,7L,8L,9L, 26L,27L,28L,29L,30L,31L,32L,33L,
                   34L, 1L,2L,3L,4L,5L,6L,7L,8L,9L, 26L,27L,28L,29L,30L,31L,32L,33L,
                  34L, 1L,2L,3L,4L,5L,6L,7L,8L,9L, 23L,24L,25L,26L,27L,28L,29L,30L, 
                 1L,2L,3L,4L,5L,6L,7L, 22L,23L,24L,25L,26L,27L,28L), addEff = addEff, domEff = rep(0, 207), intercept = 0.1)
  SP$manAddTrait(trait_yield)
  pop_good <- newPop(founderPop, simParam = SP)
  pop_good <- setPheno(pop_good, h2 = 0.8, simParam = SP)
  wtvar1[i,] <- genicVarA(pop_good)
  
  pop_good1 <- randCross(pop_good, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good1 <- setPheno(pop_good1, h2 = 0.8, simParam = SP)
  wtvar2[i,] <- genicVarA(pop_good1)
  
  pop_good_sel <- selectInd(pop_good1, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good2 <- randCross(pop_good_sel, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good2 <- setPheno(pop_good2, h2 = 0.8, simParam = SP)
  wtvar3[i,] <- genicVarA(pop_good2)
  wtgv1[i,] <- gv(pop_good_sel)
  
  pop_good_sel2 <- selectInd(pop_good2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good3 <- randCross(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good3 <- setPheno(pop_good3, h2 = 0.8, simParam = SP)
  wtvar4[i,] <- genicVarA(pop_good3)
  wtgv2[i,] <- gv(pop_good_sel2)
  
  pop_good_sel3 <- selectInd(pop_good3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good4 <- randCross(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good4 <- setPheno(pop_good4, h2 = 0.8, simParam = SP)
  wtvar5[i,] <- genicVarA(pop_good4)
  wtgv3[i,] <- gv(pop_good_sel3)
  
  pop_good_sel4 <- selectInd(pop_good4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good5 <- randCross(pop_good_sel4, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good5 <- setPheno(pop_good5, h2 = 0.8, simParam = SP)
  wtvar6[i,] <- genicVarA(pop_good5)
  wtgv4[i,] <- gv(pop_good_sel4)
  
  pop_good_sel5 <- selectInd(pop_good5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good6 <- randCross(pop_good_sel5, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good6 <- setPheno(pop_good6, h2 = 0.8, simParam = SP)
  wtvar7[i,] <- genicVarA(pop_good6)
  wtgv5[i,] <- gv(pop_good_sel5)
  
  pop_good_sel6 <- selectInd(pop_good6, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good7 <- randCross(pop_good_sel6, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good7 <- setPheno(pop_good7, h2 = 0.8, simParam = SP)
  wtvar8[i,] <- genicVarA(pop_good7)
  wtgv6[i,] <- gv(pop_good_sel6)
  
  pop_good_sel7 <- selectInd(pop_good7, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good8 <- randCross(pop_good_sel7, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good8 <- setPheno(pop_good8, h2 = 0.8, simParam = SP)
  wtvar9[i,] <- genicVarA(pop_good8)
  wtgv7[i,] <- gv(pop_good_sel7)
  
  pop_good_sel8 <- selectInd(pop_good8, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good9 <- randCross(pop_good_sel8, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good9 <- setPheno(pop_good9, h2 = 0.8, simParam = SP)
  wtvar10[i,] <- genicVarA(pop_good9)
  wtgv8[i,] <- gv(pop_good_sel8)
  
  pop_good_sel9 <- selectInd(pop_good9, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good10 <- randCross(pop_good_sel9, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good10 <- setPheno(pop_good10, h2 = 0.8, simParam = SP)
  wtvar11[i,] <- genicVarA(pop_good10)
  wtgv9[i,] <- gv(pop_good_sel9)
  
  pop_good_sel10 <- selectInd(pop_good10, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good11 <- randCross(pop_good_sel10, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good11 <- setPheno(pop_good11, h2 = 0.8, simParam = SP)
  wtvar12[i,] <- genicVarA(pop_good11)
  wtgv10[i,] <- gv(pop_good_sel10)
  
  pop_good_sel11 <- selectInd(pop_good11, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good12 <- randCross(pop_good_sel11, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good12 <- setPheno(pop_good12, h2 = 0.8, simParam = SP)
  wtvar13[i,] <- genicVarA(pop_good12)
  wtgv11[i,] <- gv(pop_good_sel11)
  
  pop_good_sel12 <- selectInd(pop_good12, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good13 <- randCross(pop_good_sel12, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good13 <- setPheno(pop_good13, h2 = 0.8, simParam = SP)
  wtvar14[i,] <- genicVarA(pop_good13)
  wtgv12[i,] <- gv(pop_good_sel12)
  
  pop_good_sel13 <- selectInd(pop_good13, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good14 <- randCross(pop_good_sel13, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good14 <- setPheno(pop_good14, h2 = 0.8, simParam = SP)
  wtvar15[i,] <- genicVarA(pop_good14)
  wtgv13[i,] <- gv(pop_good_sel13)
  
  pop_good_sel14 <- selectInd(pop_good14, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good15 <- randCross(pop_good_sel14, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good15 <- setPheno(pop_good15, h2 = 0.8, simParam = SP)
  wtvar16[i,] <- genicVarA(pop_good15)
  wtgv14[i,] <- gv(pop_good_sel14)
  
  pop_good_sel15 <- selectInd(pop_good15, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good16 <- randCross(pop_good_sel14, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good16 <- setPheno(pop_good15, h2 = 0.8, simParam = SP)
  wtvar17[i,] <- genicVarA(pop_good16)
  wtgv15[i,] <- gv(pop_good_sel15)
  
  pop_good_selwt[[i]] <- selectInd(pop_good10, nInd = 2, use = "gv", trait = 1, selectop = TRUE, returnPop = TRUE, simParam = SP)
}

ddm1var1 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var2 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var3 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var4 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var5 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var6 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var7 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var8 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var9 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var10 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var11 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var12 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var13 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var14 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var15 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var16 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var17 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1gv1 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm1gv2 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm1gv3 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm1gv4 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm1gv5 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm1gv6 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm1gv7 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm1gv8 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm1gv9 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm1gv10 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm1gv11 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm1gv12 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm1gv13 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm1gv14 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm1gv15 <- matrix(data = NA, ncol = 5, nrow = 100)
pop_good_selddm1 <- vector(mode = "list", length = 100)
for(i in 1:100){
  founderPop <- quickHaplo(nInd = 200, nChr = 10, inbred = TRUE, ploidy = 2L, segSites = c(nrow(chr1_finalpos), nrow(chr2_finalpos), 
                                                                                           nrow(chr3_finalpos), nrow(chr4_finalpos), nrow(chr5_finalpos), 
                                                                                           nrow(chr6_finalpos), nrow(chr7_finalpos), nrow(chr8_finalpos),
                                                                                           nrow(chr9_finalpos), nrow(chr10_finalpos)))
  founderPop@genMap <- ddm1_map
  founderPop@centromere <- ddm1_centromere
  SP = SimParam$new(founderPop)
  SP$setTrackRec(TRUE)
  SP$p = 0.15
  #SP$addTraitADE(nQtlPerChr = 92)
  trait_yield <- new("TraitAD", nLoci = 207L, lociPerChr= c(30L, 24L, 24L, 24L, 22L, 16L, 18L, 18L, 17L, 14L),
                     lociLoc = c(1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,13L,14L,15L,45L,46L,47L,
                                 48L,49L,50L,51L,52L,53L,54L,55L,56L,58L,59L,60L,1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,
                                 35L,36L,37L,38L,39L,40L,41L,42L,43L,44L,45L,46L, 1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,
                                 34L,35L,36L,37L,38L,39L,40L,41L,42L,43L,44L,45L, 1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,
                                 36L,37L,38L,39L,40L,41L,42L,43L,44L,45L,46L,47L, 1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,
                                 32L,33L,34L,35L,36L,37L,38L,39L,40L,41L,42L, 1L,2L,3L,4L,5L,6L,7L,8L,
                                 25L,26L,27L,28L,29L,30L,31L,32L, 1L,2L,3L,4L,5L,6L,7L,8L,9L, 26L,27L,28L,29L,30L,31L,32L,33L,
                                 34L, 1L,2L,3L,4L,5L,6L,7L,8L,9L, 26L,27L,28L,29L,30L,31L,32L,33L,
                                 34L, 1L,2L,3L,4L,5L,6L,7L,8L,9L, 23L,24L,25L,26L,27L,28L,29L,30L, 
                                 1L,2L,3L,4L,5L,6L,7L, 22L,23L,24L,25L,26L,27L,28L), addEff = addEff, domEff = rep(0, 207), intercept = 0.1)
  SP$manAddTrait(trait_yield)
  pop_good <- newPop(founderPop, simParam = SP)
  pop_good <- setPheno(pop_good, h2 = 0.8, simParam = SP)
  ddm1var1[i,] <- genicVarA(pop_good)

  pop_good1 <- randCross(pop_good, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good1 <- setPheno(pop_good1, h2 = 0.8, simParam = SP)
  ddm1var2[i,] <- genicVarA(pop_good1)

  pop_good_sel <- selectInd(pop_good1, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good2 <- randCross(pop_good_sel, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good2 <- setPheno(pop_good2, h2 = 0.8, simParam = SP)
  ddm1var3[i,] <- genicVarA(pop_good2)
  ddm1gv1[i,] <- gv(pop_good_sel)
  
  pop_good_sel2 <- selectInd(pop_good2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good3 <- randCross(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good3 <- setPheno(pop_good3, h2 = 0.8, simParam = SP)
  ddm1var4[i,] <- genicVarA(pop_good3)
  ddm1gv2[i,] <- gv(pop_good_sel2)
  
  pop_good_sel3 <- selectInd(pop_good3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good4 <- randCross(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good4 <- setPheno(pop_good4, h2 = 0.8, simParam = SP)
  ddm1var5[i,] <- genicVarA(pop_good4)
  ddm1gv3[i,] <- gv(pop_good_sel3)
  
  pop_good_sel4 <- selectInd(pop_good4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good5 <- randCross(pop_good_sel4, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good5 <- setPheno(pop_good5, h2 = 0.8, simParam = SP)
  ddm1var6[i,] <- genicVarA(pop_good5)
  ddm1gv4[i,] <- gv(pop_good_sel4)
  
  pop_good_sel5 <- selectInd(pop_good5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good6 <- randCross(pop_good_sel5, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good6 <- setPheno(pop_good6, h2 = 0.8, simParam = SP)
  ddm1var7[i,] <- genicVarA(pop_good6)
  ddm1gv5[i,] <- gv(pop_good_sel5)
  
  pop_good_sel6 <- selectInd(pop_good6, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good7 <- randCross(pop_good_sel6, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good7 <- setPheno(pop_good7, h2 = 0.8, simParam = SP)
  ddm1var8[i,] <- genicVarA(pop_good7)
  ddm1gv6[i,] <- gv(pop_good_sel6)
  
  pop_good_sel7 <- selectInd(pop_good7, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good8 <- randCross(pop_good_sel7, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good8 <- setPheno(pop_good8, h2 = 0.8, simParam = SP)
  ddm1var9[i,] <-genicVarA(pop_good8)
  ddm1gv7[i,] <- gv(pop_good_sel7)
  
  pop_good_sel8 <- selectInd(pop_good8, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good9 <- randCross(pop_good_sel8, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good9 <- setPheno(pop_good9, h2 = 0.8, simParam = SP)
  ddm1var10[i,] <- genicVarA(pop_good9)
  ddm1gv8[i,] <- gv(pop_good_sel8)
  
  pop_good_sel9 <- selectInd(pop_good9, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good10 <- randCross(pop_good_sel9, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good10 <- setPheno(pop_good10, h2 = 0.8, simParam = SP)
  ddm1var11[i,] <- genicVarA(pop_good10)
  ddm1gv9[i,] <- gv(pop_good_sel9)
  
  pop_good_sel10 <- selectInd(pop_good10, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good11 <- randCross(pop_good_sel10, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good11 <- setPheno(pop_good11, h2 = 0.8, simParam = SP)
  ddm1var12[i,] <- genicVarA(pop_good11)
  ddm1gv10[i,] <- gv(pop_good_sel10)
  
  pop_good_sel11 <- selectInd(pop_good11, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good12 <- randCross(pop_good_sel11, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good12 <- setPheno(pop_good12, h2 = 0.8, simParam = SP)
  ddm1var13[i,] <- genicVarA(pop_good12)
  ddm1gv11[i,] <- gv(pop_good_sel11)
  
  pop_good_sel12 <- selectInd(pop_good12, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good13 <- randCross(pop_good_sel12, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good13 <- setPheno(pop_good13, h2 = 0.8, simParam = SP)
  ddm1var14[i,] <- genicVarA(pop_good13)
  ddm1gv12[i,] <- gv(pop_good_sel12)
  
  pop_good_sel13 <- selectInd(pop_good13, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good14 <- randCross(pop_good_sel13, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good14 <- setPheno(pop_good14, h2 = 0.8, simParam = SP)
  ddm1var15[i,] <- genicVarA(pop_good14)
  ddm1gv13[i,] <- gv(pop_good_sel13)
  
  pop_good_sel14 <- selectInd(pop_good14, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good15 <- randCross(pop_good_sel14, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good15 <- setPheno(pop_good15, h2 = 0.8, simParam = SP)
  ddm1var16[i,] <- genicVarA(pop_good15)
  ddm1gv14[i,] <- gv(pop_good_sel14)
  
  pop_good_sel15 <- selectInd(pop_good15, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good16 <- randCross(pop_good_sel15, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good16 <- setPheno(pop_good16, h2 = 0.8, simParam = SP)
  ddm1var17[i,] <- genicVarA(pop_good16)
  ddm1gv15[i,] <- gv(pop_good_sel15)
  
  pop_good_selddm1[[i]] <- selectInd(pop_good13, nInd = 2, use = "gv", trait = 1, selectop = TRUE, returnPop = TRUE, simParam = SP)
}

zmet2var1 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var2 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var3 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var4 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var5 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var6 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var7 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var8 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var9 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var10 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var11 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var12 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var13 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var14 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var15 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var16 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var17 <- matrix(data = NA, ncol = 1, nrow = 100)
pop_good_selzmet2 <- vector(mode = "list", length = 100)
zmet2gv1 <- matrix(data = NA, ncol = 5, nrow = 100)
zmet2gv2 <- matrix(data = NA, ncol = 5, nrow = 100)
zmet2gv3 <- matrix(data = NA, ncol = 5, nrow = 100)
zmet2gv4 <- matrix(data = NA, ncol = 5, nrow = 100)
zmet2gv5 <- matrix(data = NA, ncol = 5, nrow = 100)
zmet2gv6 <- matrix(data = NA, ncol = 5, nrow = 100)
zmet2gv7 <- matrix(data = NA, ncol = 5, nrow = 100)
zmet2gv8 <- matrix(data = NA, ncol = 5, nrow = 100)
zmet2gv9 <- matrix(data = NA, ncol = 5, nrow = 100)
zmet2gv10 <- matrix(data = NA, ncol = 5, nrow = 100)
zmet2gv11 <- matrix(data = NA, ncol = 5, nrow = 100)
zmet2gv12 <- matrix(data = NA, ncol = 5, nrow = 100)
zmet2gv13 <- matrix(data = NA, ncol = 5, nrow = 100)
zmet2gv14 <- matrix(data = NA, ncol = 5, nrow = 100)
zmet2gv15 <- matrix(data = NA, ncol = 5, nrow = 100)
for(i in 1:100){
  founderPop <- quickHaplo(nInd = 200, nChr = 10, inbred = TRUE, ploidy = 2L, segSites = c(nrow(chr1_finalpos), nrow(chr2_finalpos), 
                                                                                           nrow(chr3_finalpos), nrow(chr4_finalpos), nrow(chr5_finalpos), 
                                                                                           nrow(chr6_finalpos), nrow(chr7_finalpos), nrow(chr8_finalpos),
                                                                                           nrow(chr9_finalpos), nrow(chr10_finalpos)))
  founderPop@genMap <- zmet2_map
  founderPop@centromere <- zmet2_centromere
  SP = SimParam$new(founderPop)
  SP$setTrackRec(TRUE)
  SP$p = 0.15
  #SP$addTraitADE(nQtlPerChr = 92)
  trait_yield <- new("TraitAD", nLoci = 207L, lociPerChr= c(30L, 24L, 24L, 24L, 22L, 16L, 18L, 18L, 17L, 14L),
                     lociLoc = c(1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,13L,14L,15L,45L,46L,47L,
                                 48L,49L,50L,51L,52L,53L,54L,55L,56L,58L,59L,60L,1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,
                                 35L,36L,37L,38L,39L,40L,41L,42L,43L,44L,45L,46L, 1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,
                                 34L,35L,36L,37L,38L,39L,40L,41L,42L,43L,44L,45L, 1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,
                                 36L,37L,38L,39L,40L,41L,42L,43L,44L,45L,46L,47L, 1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,
                                 32L,33L,34L,35L,36L,37L,38L,39L,40L,41L,42L, 1L,2L,3L,4L,5L,6L,7L,8L,
                                 25L,26L,27L,28L,29L,30L,31L,32L, 1L,2L,3L,4L,5L,6L,7L,8L,9L, 26L,27L,28L,29L,30L,31L,32L,33L,
                                 34L, 1L,2L,3L,4L,5L,6L,7L,8L,9L, 26L,27L,28L,29L,30L,31L,32L,33L,
                                 34L, 1L,2L,3L,4L,5L,6L,7L,8L,9L, 23L,24L,25L,26L,27L,28L,29L,30L, 
                                 1L,2L,3L,4L,5L,6L,7L, 22L,23L,24L,25L,26L,27L,28L), addEff = addEff, domEff = rep(0, 207), intercept = 0.1)
  SP$manAddTrait(trait_yield)
  pop_good <- newPop(founderPop, simParam = SP)
  pop_good <- setPheno(pop_good, h2 = 0.8, simParam = SP)
  zmet2var1[i,] <- genicVarA(pop_good)
  
  pop_good1 <- randCross(pop_good, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good1 <- setPheno(pop_good1, h2 = 0.8, simParam = SP)
  zmet2var2[i,] <- genicVarA(pop_good1)
  
  pop_good_sel <- selectInd(pop_good1, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good2 <- randCross(pop_good_sel, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good2 <- setPheno(pop_good2, h2 = 0.8, simParam = SP)
  zmet2var3[i,] <- genicVarA(pop_good2)
  zmet2gv1[i,] <- gv(pop_good_sel)
  
  pop_good_sel2 <- selectInd(pop_good2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good3 <- randCross(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good3 <- setPheno(pop_good3, h2 = 0.8, simParam = SP)
  zmet2var4[i,] <- genicVarA(pop_good3)
  zmet2gv2[i,] <- gv(pop_good_sel2)
  
  pop_good_sel3 <- selectInd(pop_good3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good4 <- randCross(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good4 <- setPheno(pop_good4, h2 = 0.8, simParam = SP)
  zmet2var5[i,] <- genicVarA(pop_good4)
  zmet2gv3[i,] <- gv(pop_good_sel3)
  
  pop_good_sel4 <- selectInd(pop_good4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good5 <- randCross(pop_good_sel4, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good5 <- setPheno(pop_good5, h2 = 0.8, simParam = SP)
  zmet2var6[i,] <- genicVarA(pop_good5)
  zmet2gv4[i,] <- gv(pop_good_sel4)
  
  pop_good_sel5 <- selectInd(pop_good5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good6 <- randCross(pop_good_sel5, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good6 <- setPheno(pop_good6, h2 = 0.8, simParam = SP)
  zmet2var7[i,] <- genicVarA(pop_good6)
  zmet2gv5[i,] <- gv(pop_good_sel5)
  
  pop_good_sel6 <- selectInd(pop_good6, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good7 <- randCross(pop_good_sel6, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good7 <- setPheno(pop_good7, h2 = 0.8, simParam = SP)
  zmet2var8[i,] <- genicVarA(pop_good7)
  zmet2gv6[i,] <- gv(pop_good_sel6)
  
  pop_good_sel7 <- selectInd(pop_good7, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good8 <- randCross(pop_good_sel7, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good8 <- setPheno(pop_good8, h2 = 0.8, simParam = SP)
  zmet2var9[i,] <- genicVarA(pop_good8)
  zmet2gv7[i,] <- gv(pop_good_sel7)
  
  pop_good_sel8 <- selectInd(pop_good8, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good9 <- randCross(pop_good_sel8, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good9 <- setPheno(pop_good9, h2 = 0.8, simParam = SP)
  zmet2var10[i,] <- genicVarA(pop_good9)
  zmet2gv8[i,] <- gv(pop_good_sel8)
  
  pop_good_sel9 <- selectInd(pop_good9, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good10 <- randCross(pop_good_sel9, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good10 <- setPheno(pop_good10, h2 = 0.8, simParam = SP)
  zmet2var11[i,] <- genicVarA(pop_good10)
  zmet2gv9[i,] <- gv(pop_good_sel9)
  
  pop_good_sel10 <- selectInd(pop_good10, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good11 <- randCross(pop_good_sel10, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good11 <- setPheno(pop_good11, h2 = 0.8, simParam = SP)
  zmet2var12[i,] <- genicVarA(pop_good11)
  zmet2gv10[i,] <- gv(pop_good_sel10)
  
  pop_good_sel11 <- selectInd(pop_good11, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good12 <- randCross(pop_good_sel11, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good12 <- setPheno(pop_good12, h2 = 0.8, simParam = SP)
  zmet2var13[i,] <- genicVarA(pop_good12)
  zmet2gv11[i,] <- gv(pop_good_sel11)
  
  pop_good_sel12 <- selectInd(pop_good12, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good13 <- randCross(pop_good_sel12, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good13 <- setPheno(pop_good13, h2 = 0.8, simParam = SP)
  zmet2var14[i,] <- genicVarA(pop_good13)
  zmet2gv12[i,] <- gv(pop_good_sel12)
  
  pop_good_sel13 <- selectInd(pop_good13, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good14 <- randCross(pop_good_sel13, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good14 <- setPheno(pop_good14, h2 = 0.8, simParam = SP)
  zmet2var15[i,] <- genicVarA(pop_good14)
  zmet2gv13[i,] <- gv(pop_good_sel13)
  
  pop_good_sel14 <- selectInd(pop_good14, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good15 <- randCross(pop_good_sel14, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good15 <- setPheno(pop_good15, h2 = 0.8, simParam = SP)
  zmet2var16[i,] <- genicVarA(pop_good15)
  zmet2gv14[i,] <- gv(pop_good_sel14)
  
  pop_good_sel15 <- selectInd(pop_good15, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good16 <- randCross(pop_good_sel15, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good16 <- setPheno(pop_good16, h2 = 0.8, simParam = SP)
  zmet2var17[i,] <- genicVarA(pop_good16)
  zmet2gv15[i,] <- gv(pop_good_sel15)
  
  pop_good_selzmet2[[i]] <- selectInd(pop_good10, nInd = 2, use = "gv", trait = 1, selectop = TRUE, returnPop = TRUE, simParam = SP)
}

recq4var1 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var2 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var3 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var4 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var5 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var6 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var7 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var8 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var9 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var10 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var11 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var12 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var13 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var14 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var15 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var16 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var17 <- matrix(data = NA, ncol = 1, nrow = 100)
pop_good_selrecq4 <- vector(mode = "list", length = 100)
recq4gv1 <- matrix(data = NA, ncol = 5, nrow = 100)
recq4gv2 <- matrix(data = NA, ncol = 5, nrow = 100)
recq4gv3 <- matrix(data = NA, ncol = 5, nrow = 100)
recq4gv4 <- matrix(data = NA, ncol = 5, nrow = 100)
recq4gv5 <- matrix(data = NA, ncol = 5, nrow = 100)
recq4gv6 <- matrix(data = NA, ncol = 5, nrow = 100)
recq4gv7 <- matrix(data = NA, ncol = 5, nrow = 100)
recq4gv8 <- matrix(data = NA, ncol = 5, nrow = 100)
recq4gv9 <- matrix(data = NA, ncol = 5, nrow = 100)
recq4gv10 <- matrix(data = NA, ncol = 5, nrow = 100)
recq4gv11 <- matrix(data = NA, ncol = 5, nrow = 100)
recq4gv12 <- matrix(data = NA, ncol = 5, nrow = 100)
recq4gv13 <- matrix(data = NA, ncol = 5, nrow = 100)
recq4gv14 <- matrix(data = NA, ncol = 5, nrow = 100)
recq4gv15 <- matrix(data = NA, ncol = 5, nrow = 100)
for(i in 1:100){
  founderPop <- quickHaplo(nInd = 200, nChr = 10, inbred = TRUE, ploidy = 2L, segSites = c(nrow(chr1_finalpos), nrow(chr2_finalpos), 
                                                                                           nrow(chr3_finalpos), nrow(chr4_finalpos), nrow(chr5_finalpos), 
                                                                                           nrow(chr6_finalpos), nrow(chr7_finalpos), nrow(chr8_finalpos),
                                                                                           nrow(chr9_finalpos), nrow(chr10_finalpos)))
  founderPop@genMap <- recq4_map
  founderPop@centromere <- recq4_centromere
  SP = SimParam$new(founderPop)
  SP$setTrackRec(TRUE)
  SP$p = 0.15
  #SP$addTraitADE(nQtlPerChr = 92)
  trait_yield <- new("TraitA", nLoci = 207L, lociPerChr= c(30L, 24L, 24L, 24L, 22L, 16L, 18L, 18L, 17L, 14L),
                     lociLoc = c(1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,13L,14L,15L,45L,46L,47L,
                                 48L,49L,50L,51L,52L,53L,54L,55L,56L,58L,59L,60L,1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,
                                 35L,36L,37L,38L,39L,40L,41L,42L,43L,44L,45L,46L, 1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,
                                 34L,35L,36L,37L,38L,39L,40L,41L,42L,43L,44L,45L, 1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,
                                 36L,37L,38L,39L,40L,41L,42L,43L,44L,45L,46L,47L, 1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,
                                 32L,33L,34L,35L,36L,37L,38L,39L,40L,41L,42L, 1L,2L,3L,4L,5L,6L,7L,8L,
                                 25L,26L,27L,28L,29L,30L,31L,32L, 1L,2L,3L,4L,5L,6L,7L,8L,9L, 26L,27L,28L,29L,30L,31L,32L,33L,
                                 34L, 1L,2L,3L,4L,5L,6L,7L,8L,9L, 26L,27L,28L,29L,30L,31L,32L,33L,
                                 34L, 1L,2L,3L,4L,5L,6L,7L,8L,9L, 23L,24L,25L,26L,27L,28L,29L,30L, 
                                 1L,2L,3L,4L,5L,6L,7L, 22L,23L,24L,25L,26L,27L,28L), addEff = addEff, intercept = 0.1)
  SP$manAddTrait(trait_yield)
  pop_good <- newPop(founderPop, simParam = SP)
  pop_good <- setPheno(pop_good, h2 = 0.8, simParam = SP)
  recq4var1[i,] <- genicVarA(pop_good)
  
  pop_good1 <- randCross(pop_good, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good1 <- setPheno(pop_good1, h2 = 0.8, simParam = SP)
  recq4var2[i,] <- genicVarA(pop_good1)
  
  pop_good_sel <- selectInd(pop_good1, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good2 <- randCross(pop_good_sel, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good2 <- setPheno(pop_good2, h2 = 0.8, simParam = SP)
  recq4var3[i,] <- genicVarA(pop_good2)
  recq4gv1[i,] <- gv(pop_good_sel)
  
  pop_good_sel2 <- selectInd(pop_good2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good3 <- randCross(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good3 <- setPheno(pop_good3, h2 = 0.8, simParam = SP)
  recq4var4[i,] <- genicVarA(pop_good3)
  recq4gv2[i,] <- gv(pop_good_sel2)
  
  pop_good_sel3 <- selectInd(pop_good3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good4 <- randCross(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good4 <- setPheno(pop_good4, h2 = 0.8, simParam = SP)
  recq4var5[i,] <- genicVarA(pop_good4)
  recq4gv3[i,] <- gv(pop_good_sel3)
  
  pop_good_sel4 <- selectInd(pop_good4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good5 <- randCross(pop_good_sel4, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good5 <- setPheno(pop_good5, h2 = 0.8, simParam = SP)
  recq4var6[i,] <- genicVarA(pop_good5)
  recq4gv4[i,] <- gv(pop_good_sel4)
  
  pop_good_sel5 <- selectInd(pop_good5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good6 <- randCross(pop_good_sel5, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good6 <- setPheno(pop_good6, h2 = 0.8, simParam = SP)
  recq4var7[i,] <- genicVarA(pop_good6)
  recq4gv5[i,] <- gv(pop_good_sel5)
  
  pop_good_sel6 <- selectInd(pop_good6, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good7 <- randCross(pop_good_sel6, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good7 <- setPheno(pop_good7, h2 = 0.8, simParam = SP)
  recq4var8[i,] <- genicVarA(pop_good7)
  recq4gv6[i,] <- gv(pop_good_sel6)
  
  pop_good_sel7 <- selectInd(pop_good7, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good8 <- randCross(pop_good_sel7, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good8 <- setPheno(pop_good8, h2 = 0.8, simParam = SP)
  recq4var9[i,] <- genicVarA(pop_good8)
  recq4gv7[i,] <- gv(pop_good_sel7)
  
  pop_good_sel8 <- selectInd(pop_good8, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good9 <- randCross(pop_good_sel8, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good9 <- setPheno(pop_good9, h2 = 0.8, simParam = SP)
  recq4var10[i,] <- genicVarA(pop_good9)
  recq4gv8[i,] <- gv(pop_good_sel8)
  
  pop_good_sel9 <- selectInd(pop_good9, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good10 <- randCross(pop_good_sel9, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good10 <- setPheno(pop_good10, h2 = 0.8, simParam = SP)
  recq4var11[i,] <- genicVarA(pop_good10)
  recq4gv9[i,] <- gv(pop_good_sel9)
  
  pop_good_sel10 <- selectInd(pop_good10, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good11 <- randCross(pop_good_sel10, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good11 <- setPheno(pop_good11, h2 = 0.8, simParam = SP)
  recq4var12[i,] <- genicVarA(pop_good11)
  recq4gv10[i,] <- gv(pop_good_sel10)
  
  pop_good_sel11 <- selectInd(pop_good11, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good12 <- randCross(pop_good_sel11, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good12 <- setPheno(pop_good12, h2 = 0.8, simParam = SP)
  recq4var13[i,] <- genicVarA(pop_good12)
  recq4gv11[i,] <- gv(pop_good_sel11)
  
  pop_good_sel12 <- selectInd(pop_good12, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good13 <- randCross(pop_good_sel12, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good13 <- setPheno(pop_good13, h2 = 0.8, simParam = SP)
  recq4var14[i,] <- genicVarA(pop_good13)
  recq4gv12[i,] <- gv(pop_good_sel12)
  
  pop_good_sel13 <- selectInd(pop_good13, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good14 <- randCross(pop_good_sel13, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good14 <- setPheno(pop_good14, h2 = 0.8, simParam = SP)
  recq4var15[i,] <- genicVarA(pop_good14)
  recq4gv13[i,] <- gv(pop_good_sel13)
  
  pop_good_sel14 <- selectInd(pop_good14, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good15 <- randCross(pop_good_sel14, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good15 <- setPheno(pop_good15, h2 = 0.8, simParam = SP)
  recq4var16[i,] <- genicVarA(pop_good15)
  recq4gv14[i,] <- gv(pop_good_sel14)
  
  pop_good_sel15 <- selectInd(pop_good15, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good16 <- randCross(pop_good_sel15, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good16 <- setPheno(pop_good16, h2 = 0.8, simParam = SP)
  recq4var17[i,] <- genicVarA(pop_good16)
  recq4gv15[i,] <- gv(pop_good_sel15)
  
  pop_good_selrecq4[[i]] <- selectInd(pop_good10, nInd = 2, use = "gv", trait = 1, selectop = TRUE, returnPop = TRUE, simParam = SP)
}

wtvarall <- c(mean(wtvar1), mean(wtvar2), mean(wtvar3), mean(wtvar4), mean(wtvar5),
        mean(wtvar6),  mean(wtvar7),  mean(wtvar8),  mean(wtvar9),  mean(wtvar10), mean(wtvar11),
        mean(wtvar12), mean(wtvar13),  mean(wtvar14),  mean(wtvar15),  mean(wtvar16),  mean(wtvar17))

ddm1varall <- c(mean(ddm1var1), mean(ddm1var2), mean(ddm1var3), mean(ddm1var4), mean(ddm1var5),
              mean(ddm1var6),  mean(ddm1var7),  mean(ddm1var8),  mean(ddm1var9),  mean(ddm1var10), 
              mean(ddm1var11), mean(ddm1var12), mean(ddm1var13), mean(ddm1var14), mean(ddm1var15), mean(ddm1var16), mean(ddm1var17))

zmet2varall <- c(mean(zmet2var1), mean(zmet2var2), mean(zmet2var3), mean(zmet2var4), mean(zmet2var5),
                 mean(zmet2var6),  mean(zmet2var7),  mean(zmet2var8),  mean(zmet2var9),  mean(zmet2var10),
                 mean(zmet2var11), mean(zmet2var12), mean(zmet2var13), mean(zmet2var14), mean(zmet2var15), mean(zmet2var16), mean(zmet2var17))

fanmcvarall <- c(mean(recq4var1), mean(recq4var2), mean(recq4var3), mean(recq4var4), mean(recq4var5),
                 mean(recq4var6),  mean(recq4var7),  mean(recq4var8),  mean(recq4var9),  mean(recq4var10),
                 mean(recq4var11), mean(recq4var12), mean(recq4var13), mean(recq4var14), mean(recq4var15), mean(recq4var16), mean(recq4var17))

wt <- c(wtgv1[,1:5], wtgv2[,1:5], wtgv3[,1:5], wtgv4[,1:5], wtgv5[,1:5],
        wtgv6[,1:5], wtgv7[,1:5], wtgv8[,1:5], wtgv9[,1:5], wtgv10[,1:5], wtgv11[,1:5],
        wtgv12[,1:5], wtgv13[,1:5], wtgv14[,1:5], wtgv15[,1:5])

ddm1 <- c(ddm1gv1[,1:5], ddm1gv2[,1:5], ddm1gv3[,1:5], ddm1gv4[,1:5], ddm1gv5[,1:5],
          ddm1gv6[,1:5], ddm1gv7[,1:5], ddm1gv8[,1:5], ddm1gv9[,1:5], ddm1gv10[,1:5], ddm1gv11[,1:5],
          ddm1gv12[,1:5], ddm1gv13[,1:5], ddm1gv14[,1:5], ddm1gv15[,1:5])

zmet2 <- c(zmet2gv1[,1:5], zmet2gv2[,1:5], zmet2gv3[,1:5], zmet2gv4[,1:5], zmet2gv5[,1:5],
           zmet2gv6[,1:5], zmet2gv7[,1:5], zmet2gv8[,1:5], zmet2gv9[,1:5], zmet2gv10[,1:5], 
           zmet2gv11[,1:5], zmet2gv12[,1:5], zmet2gv13[,1:5], zmet2gv14[,1:5], zmet2gv15[,1:5])

recq4 <- c(recq4gv1[,1:5], recq4gv2[,1:5], recq4gv3[,1:5], recq4gv4[,1:5], recq4gv5[,1:5],
           recq4gv6[,1:5], recq4gv7[,1:5], recq4gv8[,1:5], recq4gv9[,1:5], recq4gv10[,1:5], 
           recq4gv11[,1:5], recq4gv12[,1:5], recq4gv13[,1:5], recq4gv14[,1:5], recq4gv15[,1:5])

all <- cbind(wt, zmet2, ddm1, recq4)
all <- as.data.frame(all)
all2 <- as.data.frame(matrix(data = NA, nrow = 30000))
all2$gv <- c(all$wt, all$ddm1, all$zmet2, all$recq4)
all2$generation <- rep(1:15, size = 4, each = 500)
all2$gen_map <- rep(c("wt","ddm1", "zmet", "recq4"), size = 4, each = 7500)
lastgen <- all2[which(all2$generation == 9),]

intermating <- aov(gv ~ as.factor(gen_map), data = lastgen)
intermating_gen <- aov(gv ~ as.factor(gen_map)*as.factor(generation), data = all2)
#full <- anova(intermating, intermating_gen)
summary(intermating)
summary(intermating_gen)
TukeyHSD(intermating, conf.level = .95)
TukeyHSD(intermating_gen, conf.level=.95)
plot(TukeyHSD(intermating_gen, conf.level=.95), las = 2)

ggqqplot(residuals(intermating_gen))

ggplot(all2, aes(x=as.factor(generation), y=gv, fill=gen_map)) + 
  geom_boxplot() + theme_bw() + xlab("Generation") + ylab("GVs") + 
  scale_fill_manual(values=c("#aae4c2", "#fff87a", "#35c7ff", '#ff8b77')) + ggtitle("Intermating GVs with accurate gene space") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=14), legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14), plot.title = element_text(size=20))

wtvarall <- c(wtvar1, wtvar2, wtvar3, wtvar4, wtvar5,
              wtvar6,  wtvar7,  wtvar8,  wtvar9,  wtvar10, wtvar11,
              wtvar12, wtvar13,  wtvar14,  wtvar15,  wtvar16,  wtvar17)

ddm1varall <- c(ddm1var1, ddm1var2, ddm1var3, ddm1var4, ddm1var5,
                ddm1var6,  ddm1var7,  ddm1var8,  ddm1var9,  ddm1var10, 
                ddm1var11, ddm1var12, ddm1var13, ddm1var14, ddm1var15, ddm1var16, ddm1var17)

zmet2varall <- c(zmet2var1, zmet2var2, zmet2var3, zmet2var4, zmet2var5,
                 zmet2var6,  zmet2var7,  zmet2var8,  zmet2var9,  zmet2var10,
                 zmet2var11, zmet2var12, zmet2var13, zmet2var14, zmet2var15, zmet2var16, zmet2var17)

recq4varall <- c(recq4var1, recq4var2, recq4var3, recq4var4, recq4var5,
                 recq4var6,  recq4var7,  recq4var8,  recq4var9,  recq4var10,
                 recq4var11, recq4var12, recq4var13, recq4var14, recq4var15, recq4var16, recq4var17)

allvar <- cbind(wtvarall, ddm1varall, zmet2varall, recq4varall)
allvar <- as.data.frame(allvar)
allvar2 <- as.data.frame(matrix(data = NA, nrow = 6800))
allvar2$var <- c(allvar$wtvarall, allvar$ddm1varall, allvar$zmet2varall, allvar$recq4varall)
allvar2$gen <- rep(1:17, size = 4, each = 100)
allvar2$gen_map <- rep(c("wt","ddm1", "zmet", "recq4"), size = 4, each = 1700)
lastvar <- allvar2[which(allvar2$gen == 9),]

intermatingvar <- aov(var ~ as.factor(gen_map), data = lastvar)
summary(intermatingvar)
TukeyHSD(intermatingvar, conf.level = .95)

ggplot(allvar2, aes(x=as.factor(gen), y=var, fill=gen_map)) + 
  geom_boxplot() + theme_bw() + xlab("Generation") + ylab("Genetic Variances") + 
  scale_fill_manual(values=c("#aae4c2", "#fff87a", "#35c7ff",'#ff8b77')) + ggtitle("Intermating Genetic Variance through 15 Generations") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=14), legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14), plot.title = element_text(size=20))


save(wt, file = "wt_gene_gv.RData")
save(ddm1, file = "ddm1_gene_gv.RData")
save(zmet2, file = "zmet2_gene_gv.RData")
save(recq4, file = "recq4_gene_gv.RData")
save(zmet2varall, file = "zmet2_gene_var.RData")
save(ddm1varall, file = "ddm1_gene_var.RData")
save(wtvarall, file = "wt_gene_var.RData")
