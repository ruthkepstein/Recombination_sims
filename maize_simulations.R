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
hist(chr5_snp$`SNP Start`, breaks = 100)
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
NAM <- read.table("NAM_US_COs_v4.txt", header = TRUE)
#NAM <- NAM[,-c(2:4)]
#colnames(NAM) <- c("Chr", "CO Start", "CO End")
NAM <- NAM[order(NAM$Chr,NAM$Start),]
#pop_size <- read.table("pop_size_nam.txt", header = TRUE)
#US_pop_size <- pop_size[1:4713,]

chr1_CO <- NAM[ which(NAM$Chr == 'chr1'),]
chr1_CO$midpoint <- (chr1_CO$Start + chr1_CO$End)/2
#hist(chr1_CO$`CO Start`, breaks = 300, xlab = "Physical Pos (bp)", main = "WT Chr.1 CO Density", col = "black")

chr2_CO <- NAM[ which(NAM$Chr == 'chr2'),]
#hist(chr2_CO$`CO Start`, breaks = 300)
chr2_CO$midpoint <- (chr2_CO$Start + chr2_CO$End)/2

chr3_CO <- NAM[ which(NAM$Chr == 'chr3'),]
#hist(chr3_CO$`CO Start`, breaks = 300)
chr3_CO$midpoint <- (chr3_CO$Start + chr3_CO$End)/2

chr4_CO <- NAM[ which(NAM$Chr == 'chr4'),]
#hist(chr4_CO$`CO Start`, breaks = 300)
chr4_CO$midpoint <- (chr4_CO$Start + chr4_CO$End)/2

chr5_CO <- NAM[ which(NAM$Chr == 'chr5'),]
#hist(chr5_CO$`CO Start`, breaks = 300)
chr5_CO$midpoint <- (chr5_CO$Start + chr5_CO$End)/2

chr6_CO <- NAM[ which(NAM$Chr == 'chr6'),]
#hist(chr6_CO$`CO Start`, breaks = 300)
chr6_CO$midpoint <- (chr6_CO$Start + chr6_CO$End)/2

chr7_CO <- NAM[ which(NAM$Chr == 'chr7'),]
#hist(chr7_CO$`CO Start`, breaks = 300)
chr7_CO$midpoint <- (chr7_CO$Start + chr7_CO$End)/2

chr8_CO <- NAM[ which(NAM$Chr == 'chr8'),]
#hist(chr8_CO$`CO Start`, breaks = 300)
chr8_CO$midpoint <- (chr8_CO$Start + chr8_CO$End)/2

chr9_CO <- NAM[ which(NAM$Chr == 'chr9'),]
#hist(chr9_CO$`CO Start`, breaks = 300)
chr9_CO$midpoint <- (chr9_CO$Start + chr9_CO$End)/2

chr10_CO <- NAM[ which(NAM$Chr == 'chr10'),]
#hist(chr10_CO$`CO Start`, breaks = 300)
chr10_CO$midpoint <- (chr10_CO$Start + chr10_CO$End)/2

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
chr1_CO <- chr1_CO[order(chr1_CO$Start),]
chr1_bin <- binning(chr1_CO$midpoint, nbins = max(chr1_snp$`SNP End`)/1000000, type = "kmeans")
chr1_bin <- as.data.frame(summary(chr1_bin))
#transforming data; making bin interval into 2 columns
chr1_bin <- within(chr1_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr1_bin$levels), ',', fixed=TRUE))))
chr1_bin <- do.call(data.frame, chr1_bin)
chr1_bin <- chr1_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr1_bin <- chr1_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr1_bin[1,4] <- 547913.5
#making intervals start at 0
chr1_bin$foo.X1 <- chr1_bin$foo.X1 - 547913.5
chr1_bin$foo.X2 <- chr1_bin$foo.X2 - 547913.5
#expanding last bin to include last SNP site to avoid NAs in future
chr1_bin[max(chr1_snp$`SNP End`)/1000000,5] <- max(chr1_snp$`SNP End`)
#adding length of bin as column and making in Mb
chr1_bin$length <- (chr1_bin$foo.X2-chr1_bin$foo.X1)/1000000
chr1_bin$rate <- ((chr1_bin$freq/4713)*100)/chr1_bin$length

chr2_CO <- chr2_CO[order(chr2_CO$Start),]
chr2_bin <- binning(chr2_CO$midpoint, nbins = round(max(chr2_snp$`SNP End`)/1000000), type = "kmeans")
chr2_bin <- as.data.frame(summary(chr2_bin))
#chr2_bin$freq <- chr2_bin$freq*2/4713
chr2_bin <- within(chr2_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr2_bin$levels), ',', fixed=TRUE))))
chr2_bin <- do.call(data.frame, chr2_bin)
chr2_bin <- chr2_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr2_bin <- chr2_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr2_bin[1,4] <- 448683.5
chr2_bin$foo.X1 <- chr2_bin$foo.X1 - 448683.5
chr2_bin$foo.X2 <- chr2_bin$foo.X2 - 448683.5
chr2_bin[round(max(chr2_snp$`SNP End`)/1000000),5] <- max(chr2_snp$`SNP End`)
chr2_bin$length <- (chr2_bin$foo.X2-chr2_bin$foo.X1)/1000000
chr2_bin$rate <- ((chr2_bin$freq/4713)*100)/chr2_bin$length

#have not converted rest of chromosomes to what I did in 1 & 2
chr3_CO <- chr3_CO[order(chr3_CO$Start),]
chr3_bin <- binning(chr3_CO$midpoint, nbins = max(chr3_snp$`SNP End`)/1000000, type = "kmeans")
chr3_bin <- as.data.frame(summary(chr3_bin))
#chr3_bin$freq <- chr3_bin$freq*2/4713
chr3_bin <- within(chr3_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr3_bin$levels), ',', fixed=TRUE))))
chr3_bin <- do.call(data.frame, chr3_bin)
chr3_bin <- chr3_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr3_bin <- chr3_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr3_bin[1,4] <- 144065.5
chr3_bin$foo.X1 <- chr3_bin$foo.X1 - 144065.5
chr3_bin$foo.X2 <- chr3_bin$foo.X2 - 144065.5
chr3_bin[max(chr3_snp$`SNP End`)/1000000,5] <- max(chr3_snp$`SNP End`)
chr3_bin$length <- (chr3_bin$foo.X2-chr3_bin$foo.X1)/1000000
chr3_bin$rate <- ((chr3_bin$freq/4713)*100)/chr3_bin$length

chr4_CO <- chr4_CO[order(chr4_CO$Start),]
chr4_bin <- binning(chr4_CO$midpoint, nbins = max(chr4_snp$`SNP End`)/1000000, type = "kmeans")
chr4_bin <- as.data.frame(summary(chr4_bin))
#chr4_bin$freq <- chr4_bin$freq*2/4713
chr4_bin <- within(chr4_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr4_bin$levels), ',', fixed=TRUE))))
chr4_bin <- do.call(data.frame, chr4_bin)
chr4_bin <- chr4_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr4_bin <- chr4_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr4_bin[1,4] <- 545527.2
chr4_bin$foo.X1 <- chr4_bin$foo.X1 - 545527.2
chr4_bin$foo.X2 <- chr4_bin$foo.X2 - 545527.2
chr4_bin[max(chr2_snp$`SNP End`)/1000000,5] <- max(chr4_snp$`SNP End`)
chr4_bin$length <- (chr4_bin$foo.X2-chr4_bin$foo.X1)/1000000
chr4_bin$rate <- ((chr4_bin$freq/4713)*100)/chr4_bin$length
#chr4_bin <- chr4_bin[-242,]

chr5_CO <- chr5_CO[order(chr5_CO$Start),]
chr5_bin <- binning(chr5_CO$midpoint, nbins = max(chr5_snp$`SNP End`)/1000000, type = "kmeans")
chr5_bin <- as.data.frame(summary(chr5_bin))
#chr5_bin$freq <- chr5_bin$freq*2/4713
chr5_bin <- within(chr5_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr5_bin$levels), ',', fixed=TRUE))))
chr5_bin <- do.call(data.frame, chr5_bin)
chr5_bin <- chr5_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr5_bin <- chr5_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr5_bin[1,4] <- 268345.8
chr5_bin$foo.X1 <- chr5_bin$foo.X1 - 268345.8
chr5_bin$foo.X2 <- chr5_bin$foo.X2 - 268345.8
chr5_bin[max(chr2_snp$`SNP End`)/1000000,5] <- max(chr5_snp$`SNP End`)
chr5_bin$length <- (chr5_bin$foo.X2-chr5_bin$foo.X1)/1000000
chr5_bin$rate <- ((chr5_bin$freq/4713)*100)/chr5_bin$length

chr6_CO <- chr6_CO[order(chr6_CO$Start),]
chr6_bin <- binning(chr6_CO$midpoint, nbins = max(chr6_snp$`SNP End`)/1000000, type = "kmeans")
chr6_bin <- as.data.frame(summary(chr6_bin))
#chr6_bin$freq <- chr6_bin$freq*2/4713
chr6_bin <- within(chr6_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr6_bin$levels), ',', fixed=TRUE))))
chr6_bin <- do.call(data.frame, chr6_bin)
chr6_bin <- chr6_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr6_bin <- chr6_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr6_bin[1,4] <- 215965.5
chr6_bin$foo.X1 <- chr6_bin$foo.X1 - 215965.5
chr6_bin$foo.X2 <- chr6_bin$foo.X2 - 215965.5
chr6_bin[max(chr6_snp$`SNP End`)/1000000,5] <- max(chr6_snp$`SNP End`)
chr6_bin$length <- (chr6_bin$foo.X2-chr6_bin$foo.X1)/1000000
chr6_bin$rate <- ((chr6_bin$freq/4713)*100)/chr6_bin$length

chr7_CO <- chr7_CO[order(chr7_CO$Start),]
chr7_bin <- binning(chr7_CO$midpoint, nbins = max(chr7_snp$`SNP End`)/1000000, type = "kmeans")
chr7_bin <- as.data.frame(summary(chr7_bin))
#chr7_bin$freq <- chr7_bin$freq*2/4713
chr7_bin <- within(chr7_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr7_bin$levels), ',', fixed=TRUE))))
chr7_bin <- do.call(data.frame, chr7_bin)
chr7_bin <- chr7_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr7_bin <- chr7_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr7_bin[1,4] <- 60293.5
chr7_bin$foo.X1 <- chr7_bin$foo.X1 - 60293.5
chr7_bin$foo.X2 <- chr7_bin$foo.X2 - 60293.5
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
chr8_bin[1,4] <- 292875.5
chr8_bin$foo.X1 <- chr8_bin$foo.X1 - 292875.5
chr8_bin$foo.X2 <- chr8_bin$foo.X2 - 292875.5
chr8_bin[max(chr8_snp$`SNP End`)/1000000,5] <- max(chr8_snp$`SNP End`)
chr8_bin$length <- (chr8_bin$foo.X2-chr8_bin$foo.X1)/1000000
chr8_bin$rate <- ((chr8_bin$freq/4713)*100)/chr8_bin$length

chr9_CO <- chr9_CO[order(chr9_CO$Start),]
chr9_bin <- binning(chr9_CO$midpoint, nbins = max(chr9_snp$`SNP End`)/1000000, type = "kmeans")
chr9_bin <- as.data.frame(summary(chr9_bin))
#chr9_bin$freq <- chr9_bin$freq*2/4713
chr9_bin <- within(chr9_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr9_bin$levels), ',', fixed=TRUE))))
chr9_bin <- do.call(data.frame, chr9_bin)
chr9_bin <- chr9_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr9_bin <- chr9_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr9_bin[1,4] <- 65379.5
chr9_bin$foo.X1 <- chr9_bin$foo.X1 - 65379.5
chr9_bin$foo.X2 <- chr9_bin$foo.X2 - 65379.5
chr9_bin[max(chr9_snp$`SNP End`)/1000000,5] <- max(chr9_snp$`SNP End`)
chr9_bin$length <- (chr9_bin$foo.X2-chr9_bin$foo.X1)/1000000
chr9_bin$rate <- ((chr9_bin$freq/4713)*100)/chr9_bin$length

chr10_CO <- chr10_CO[order(chr10_CO$Start),]
chr10_bin <- binning(chr10_CO$midpoint, nbins = max(chr10_snp$`SNP End`)/1000000, type = "kmeans")
chr10_bin <- as.data.frame(summary(chr10_bin))
#chr10_bin$freq <- chr10_bin$freq*2/4713
chr10_bin <- within(chr10_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr10_bin$levels), ',', fixed=TRUE))))
chr10_bin <- do.call(data.frame, chr10_bin)
chr10_bin <- chr10_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr10_bin <- chr10_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr10_bin[1,4] <- 199814.5
chr10_bin$foo.X1 <- chr10_bin$foo.X1 - 199814.5
chr10_bin$foo.X2 <- chr10_bin$foo.X2 - 199814.5
chr10_bin[max(chr10_snp$`SNP End`)/1000000,5] <- max(chr10_snp$`SNP End`)
chr10_bin$length <- (chr10_bin$foo.X2-chr10_bin$foo.X1)/1000000
chr10_bin$rate <- ((chr10_bin$freq/4713)*100)/chr10_bin$length

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

#function to assign genetic positions to each SNP given SNP before
gen_pos <- function(SNP){
  SNP$pos <- NA
  SNP$pos[1]<-SNP$`SNP Start`[1]*SNP$rate[1]
  for(i in 1:nrow(SNP)){
    if(i>1){
      SNP$pos[i]<- SNP$pos[i-1] + (SNP$`SNP Start`[i] - SNP$`SNP Start`[i-1])*SNP$rate[i]
    }
  }
  print(SNP$pos)
}
#function to graph recombination rate after SNP position has been found
graph_recomb <- function(SNP){
  SNP$pos2[1]<-SNP$`SNP Start`[1]*SNP$rate[1]
  for(i in 1:nrow(SNP)){
    if(i>1){
      SNP$pos2[i]<- SNP$pos[i] - (SNP$pos[i-1] - (SNP$`SNP Start`[i] + SNP$`SNP Start`[i-1])*SNP$rate[i])
    }
  }
  print(SNP$pos2)
}

#Calculating genetic maps for all chromosomes
chr1_snp2 <- snp_rate(chr1_bin, chr1_snp)
chr1_snp2$`SNP Start`<- chr1_snp2$`SNP Start`/1000000
chr1_snp2 <- chr1_snp2[order(chr1_snp2$`SNP Start`),]
#creation of genetic positions from smoothed recombination rate
chr1_snp2$pos <- gen_pos(chr1_snp2)
#graph to look at Mb vs. cM along chromosome
plot(chr1_snp2$`SNP Start`, chr1_snp2$pos, type = "l")
#graph to look at Mb vs. cM/Mb to see recombination rate along chromosome
chr1_snp2$pos2 <- graph_recomb(chr1_snp2)
plot(chr1_snp2$`SNP Start`, chr1_snp2$pos2/chr1_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 1 Recombination Distribution")
chr1_finalpos <- chr1_snp2[order(chr1_snp2$pos),]
#want False to input into AlphaSimR
is.unsorted(chr1_finalpos$pos)
#plot again to make sure it looks the same
plot(chr1_finalpos$`SNP Start`, chr1_finalpos$pos)


chr2_snp2 <- snp_rate(chr2_bin, chr2_snp)
chr2_snp2$`SNP Start` <- chr2_snp2$`SNP Start`/1000000
chr2_snp2$pos <- gen_pos(chr2_snp2)
plot(chr2_snp2$`SNP Start`, chr2_snp2$pos)
graph_wt2 <- graph_recomb(chr2_snp2)
plot(chr2_snp2$`SNP Start`, graph_wt2/chr2_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 2 Recombination Distribution")
chr2_finalpos <- chr2_snp2[order(chr2_snp2$pos),]
is.unsorted(chr2_finalpos$pos)

chr3_snp2 <- snp_rate(chr3_bin, chr3_snp)
chr3_snp2$`SNP Start` <- chr3_snp2$`SNP Start`/1000000
chr3_snp2$pos <- gen_pos(chr3_snp2)
plot(chr3_snp2$`SNP Start`, chr3_snp2$pos)
graph_wt3 <- graph_recomb(chr3_snp2)
plot(chr3_snp2$`SNP Start`, graph_wt3/chr3_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 3 Recombination Distribution")
chr3_finalpos <- chr3_snp2[order(chr3_snp2$pos),]
is.unsorted(chr3_finalpos$pos)

chr4_snp2 <- snp_rate(chr4_bin, chr4_snp)
chr4_snp2$`SNP Start` <- chr4_snp2$`SNP Start`/1000000
chr4_snp2$pos <- gen_pos(chr4_snp2)
plot(chr4_snp2$`SNP Start`, chr4_snp2$pos)
graph_wt4 <- graph_recomb(chr4_snp2)
plot(chr4_snp2$`SNP Start`, graph_wt4/chr4_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 4 Recombination Distribution")
chr4_finalpos <- chr4_snp2[order(chr4_snp2$pos),]
is.unsorted(chr4_finalpos$pos)

chr5_snp2 <- snp_rate(chr5_bin, chr5_snp)
chr5_snp2$`SNP Start` <- chr5_snp2$`SNP Start`/1000000
chr5_snp2$pos <- gen_pos(chr5_snp2)
plot(chr5_snp2$`SNP Start`, chr5_snp2$pos)
graph_wt5 <- graph_recomb(chr5_snp2)
plot(chr5_snp2$`SNP Start`, graph_wt5/chr5_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 5 Recombination Distribution")
chr5_finalpos <- chr5_snp2[order(chr5_snp2$pos),]
is.unsorted(chr5_finalpos$pos)

chr6_snp2 <- snp_rate(chr6_bin, chr6_snp)
chr6_snp2$`SNP Start` <- chr6_snp2$`SNP Start`/1000000
chr6_snp2$pos <- gen_pos(chr6_snp2)
plot(chr6_snp2$`SNP Start`, chr6_snp2$pos)
graph_wt6 <- graph_recomb(chr6_snp2)
plot(chr6_snp2$`SNP Start`, graph_wt6/chr6_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 6 Recombination Distribution")
chr6_finalpos <- chr6_snp2[order(chr6_snp2$pos),]
is.unsorted(chr6_finalpos$pos)

chr7_snp2 <- snp_rate(chr7_bin, chr7_snp)
chr7_snp2$`SNP Start` <- chr7_snp2$`SNP Start`/1000000
chr7_snp2$pos <- gen_pos(chr7_snp2)
plot(chr7_snp2$`SNP Start`, chr7_snp2$pos)
graph_wt7 <- graph_recomb(chr7_snp2)
plot(chr7_snp2$`SNP Start`, graph_wt7/chr7_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 7 Recombination Distribution")
chr7_finalpos <- chr7_snp2[order(chr7_snp2$pos),]
is.unsorted(chr7_finalpos$pos)

chr8_snp2 <- snp_rate(chr8_bin, chr8_snp)
chr8_snp2$`SNP Start` <- chr8_snp2$`SNP Start`/1000000
chr8_snp2$pos <- gen_pos(chr8_snp2)
plot(chr8_snp2$`SNP Start`, chr8_snp2$pos)
graph_wt8 <- graph_recomb(chr8_snp2)
plot(chr8_snp2$`SNP Start`, graph_wt8/chr8_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 8 Recombination Distribution")
chr8_finalpos <- chr8_snp2[order(chr8_snp2$pos),]
is.unsorted(chr8_finalpos$pos)

chr9_snp2 <- snp_rate(chr9_bin, chr9_snp)
chr9_snp2$`SNP Start` <- chr9_snp2$`SNP Start`/1000000
chr9_snp2$pos <- gen_pos(chr9_snp2)
plot(chr9_snp2$`SNP Start`, chr9_snp2$pos)
graph_wt9 <- graph_recomb(chr9_snp2)
plot(chr9_snp2$`SNP Start`, graph_wt9/chr9_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 9 Recombination Distribution")
chr9_finalpos <- chr9_snp2[order(chr9_snp2$pos),]
is.unsorted(chr9_finalpos$pos)

chr10_snp2 <- snp_rate(chr10_bin, chr10_snp)
chr10_snp2$`SNP Start` <- chr10_snp2$`SNP Start`/1000000
chr10_snp2$pos <- gen_pos(chr10_snp2)
plot(chr10_snp2$`SNP Start`, chr10_snp2$pos)
graph_wt10 <- graph_recomb(chr10_snp2)
plot(chr10_snp2$`SNP Start`, graph_wt10/chr10_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 10 Recombination Distribution")
chr10_finalpos <- chr10_snp2[order(chr10_snp2$pos),]
is.unsorted(chr10_finalpos$pos)

saveRDS(chr1_finalpos, "chr1_finalpos.RData")
saveRDS(chr2_finalpos, "chr2_finalpos.RData")
saveRDS(chr3_finalpos, "chr3_finalpos.RData")
saveRDS(chr4_finalpos, "chr4_finalpos.RData")
saveRDS(chr5_finalpos, "chr5_finalpos.RData")
saveRDS(chr6_finalpos, "chr6_finalpos.RData")
saveRDS(chr7_finalpos, "chr7_finalpos.RData")
saveRDS(chr8_finalpos, "chr8_finalpos.RData")
saveRDS(chr9_finalpos, "chr9_finalpos.RData")
saveRDS(chr10_finalpos, "chr10_finalpos.RData")

#Putting the final genetic map together
chr1 <- chr1_finalpos$pos/100

chr2 <- chr2_finalpos$pos/100

chr3 <- chr3_finalpos$pos/100

chr4 <- chr4_finalpos$pos/100

chr10 <- chr10_finalpos$pos/100

chr5 <- chr5_finalpos$pos/100

chr6 <- chr6_finalpos$pos/100

chr7 <- chr7_finalpos$pos/100

chr8 <- chr8_finalpos$pos/100

chr9 <- chr9_finalpos$pos/100

final_map = vector("list",10)
final_map[[1]] = chr1
final_map[[2]] = chr2
final_map[[3]] = chr3
final_map[[4]] = chr4
final_map[[5]] = chr5
final_map[[6]] = chr6
final_map[[7]] = chr7
final_map[[8]] = chr8
final_map[[9]] = chr9
final_map[[10]] = chr10
for(i in 1:10){
  names(final_map[[i]]) = paste(i, 1:segSites[i], sep="_")
}

#adding names to each individual SNP based on loci and chromosome #

saveRDS(final_map, file = "final_map.RData")
