library(AlphaSimR)
library(Rcpp)

setwd("C:/Users/16192/Documents/PNAS_Simulations")

##reading in SNPs from B73xMo17 based on v4 B73 ref
final_snps <- read.table("SNP_V4.bed", header = FALSE)
colnames(final_snps) <- c("Chr#", "SNP Start", "SNP End")
chr1_snp <- final_snps[ which(final_snps$`Chr#` == "chr1"),]
hist(chr1_snp$`SNP Start`, breaks = 100)
chr1_snp$rate <- NA
#making SNPs start at 0
chr1_snp$`SNP Start` <- chr1_snp$`SNP Start`- 46746
chr1_snp$`SNP End` <- chr1_snp$`SNP End` - 46746
chr1_snp <- chr1_snp[order(chr1_snp$`SNP Start`),]

chr2_snp <- final_snps[ which(final_snps$`Chr#` == "chr2"),]
hist(chr2_snp$`SNP Start`, breaks = 100)
chr2_snp$rate <- NA
chr2_snp$`SNP Start` <- chr2_snp$`SNP Start`- 366079
chr2_snp$`SNP End` <- chr2_snp$`SNP End` - 366079
chr2_snp <- chr2_snp[order(chr2_snp$`SNP Start`),]

chr3_snp <- final_snps[ which(final_snps$`Chr#` == "chr3"),]
hist(chr3_snp$`SNP Start`, breaks = 100)
chr3_snp$rate <- NA
chr3_snp$`SNP Start` <- chr3_snp$`SNP Start`- 10863
chr3_snp$`SNP End` <- chr3_snp$`SNP End` - 10863
min3 <-min(chr3_snp[,2])
chr3_snp <- chr3_snp[order(chr3_snp$`SNP Start`),]

chr4_snp <- final_snps[ which(final_snps$`Chr#` == "chr4"),]
hist(chr4_snp$`SNP Start`, breaks = 100)
chr4_snp$rate <- NA
min4 <-min(chr4_snp[,2])
chr4_snp$`SNP Start` <- chr4_snp$`SNP Start`- 157899
chr4_snp$`SNP End` <- chr4_snp$`SNP End` - 157899
chr4_snp <- chr4_snp[order(chr4_snp$`SNP Start`),]

chr5_snp <- final_snps[ which(final_snps$`Chr#` == "chr5"),]
hist(chr5_snp$`SNP Start`, breaks = 100)
chr5_snp$rate <- NA
min5 <-min(chr5_snp[,2])
chr5_snp$`SNP Start` <- chr5_snp$`SNP Start`- 3291566
chr5_snp$`SNP End` <- chr5_snp$`SNP End` - 3291566
chr5_snp <- chr5_snp[order(chr5_snp$`SNP Start`),]

chr6_snp <- final_snps[ which(final_snps$`Chr#` == "chr6"),]
hist(chr6_snp$`SNP Start`, breaks = 100)
chr6_snp$rate <- NA
min6 <-min(chr6_snp[,2])
chr6_snp$`SNP Start` <- chr6_snp$`SNP Start`- 151963
chr6_snp$`SNP End` <- chr6_snp$`SNP End` - 151963
chr6_snp <- chr6_snp[order(chr6_snp$`SNP Start`),]

chr7_snp <- final_snps[ which(final_snps$`Chr#` == "chr7"),]
hist(chr7_snp$`SNP Start`, breaks = 100)
chr7_snp$rate <- NA
min7 <-min(chr7_snp[,2])
chr7_snp$`SNP Start` <- chr7_snp$`SNP Start`- 43560
chr7_snp$`SNP End` <- chr7_snp$`SNP End` - 43560
chr7_snp <- chr7_snp[order(chr7_snp$`SNP Start`),]

chr8_snp <- final_snps[ which(final_snps$`Chr#` == "chr8"),]
hist(chr8_snp$`SNP Start`, breaks = 100)
chr8_snp$rate <- NA
min8 <-min(chr8_snp[,2])
chr8_snp$`SNP Start` <- chr8_snp$`SNP Start`- 174210
chr8_snp$`SNP End` <- chr8_snp$`SNP End` - 174210
chr8_snp <- chr8_snp[order(chr8_snp$`SNP Start`),]

chr9_snp <- final_snps[ which(final_snps$`Chr#` == "chr9"),]
hist(chr9_snp$`SNP Start`, breaks = 100)
chr9_snp$rate <- NA
min9<-min(chr9_snp[,2])
chr9_snp$`SNP Start` <- chr9_snp$`SNP Start`- 2812
chr9_snp$`SNP End` <- chr9_snp$`SNP End` - 2812
chr9_snp <- chr9_snp[order(chr9_snp$`SNP Start`),]

chr10_snp <- final_snps[ which(final_snps$`Chr#` == "chr10"),]
hist(chr10_snp$`SNP Start`, breaks = 100)
chr10_snp$rate <- NA
min10 <-min(chr10_snp[,2])
chr10_snp$`SNP Start` <- chr10_snp$`SNP Start`- 17232
chr10_snp$`SNP End` <- chr10_snp$`SNP End` - 17232
chr10_snp <- chr10_snp[order(chr10_snp$`SNP Start`),]

##reading in COs in B73xMo17 uplifted to v4 & finding midpoint of CO interval
final_COs <- read.table("Final_CO.csv", header = TRUE, sep = ",")
colnames(final_COs) <- c("Chr#", "CO length", "CO Start", "CO End")
final_COs$`CO length` <- final_COs$`CO length`/1000

chr1_CO <- final_COs[ which(final_COs$`Chr#` == "chr1"),]
chr1_CO$midpoint <- (chr1_CO$`CO Start`+ chr1_CO$`CO End`)/2
hist(chr1_CO$`CO Start`, breaks = 100)

chr2_CO <- final_COs[ which(final_COs$`Chr#` == "chr2"),]
hist(chr2_CO$`CO Start`, breaks = 100)
chr2_CO$midpoint <- (chr2_CO$`CO Start`+ chr2_CO$`CO End`)/2

chr3_CO <- final_COs[ which(final_COs$`Chr#` == "chr3"),]
hist(chr3_CO$`CO Start`, breaks = 100)
chr3_CO$midpoint <- (chr3_CO$`CO Start`+ chr3_CO$`CO End`)/2

chr4_CO <- final_COs[ which(final_COs$`Chr#` == "chr4"),]
hist(chr4_CO$`CO Start`, breaks = 100)
chr4_CO$midpoint <- (chr4_CO$`CO Start`+ chr4_CO$`CO End`)/2

chr5_CO <- final_COs[ which(final_COs$`Chr#` == "chr5"),]
hist(chr5_CO$`CO Start`, breaks = 100)
chr5_CO$midpoint <- (chr5_CO$`CO Start`+ chr5_CO$`CO End`)/2

chr6_CO <- final_COs[ which(final_COs$`Chr#` == "chr6"),]
hist(chr6_CO$`CO Start`, breaks = 100)
chr6_CO$midpoint <- (chr6_CO$`CO Start`+ chr6_CO$`CO End`)/2

chr7_CO <- final_COs[ which(final_COs$`Chr#` == "chr7"),]
hist(chr7_CO$`CO Start`, breaks = 100)
chr7_CO$midpoint <- (chr7_CO$`CO Start`+ chr7_CO$`CO End`)/2

chr8_CO <- final_COs[ which(final_COs$`Chr#` == "chr8"),]
hist(chr8_CO$`CO Start`, breaks = 100)
chr8_CO$midpoint <- (chr8_CO$`CO Start`+ chr8_CO$`CO End`)/2

chr9_CO <- final_COs[ which(final_COs$`Chr#` == "chr9"),]
hist(chr9_CO$`CO Start`, breaks = 100)
chr9_CO$midpoint <- (chr9_CO$`CO Start`+ chr9_CO$`CO End`)/2

chr10_CO <- final_COs[ which(final_COs$`Chr#` == "chr10"),]
hist(chr10_CO$`CO Start`, breaks = 100)
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
chr1_bin <- binning(chr1_CO$midpoint, nbins = 100, type = "kmeans")
chr1_bin <- as.data.frame(summary(chr1_bin))
#transforming data; making bin interval into 2 columns
chr1_bin <- within(chr1_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr1_bin$levels), ',', fixed=TRUE))))
chr1_bin <- do.call(data.frame, chr1_bin)
chr1_bin <- chr1_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr1_bin <- chr1_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr1_bin[1,4] <- 33878.5
#making intervals start at 0
chr1_bin$foo.X1 <- chr1_bin$foo.X1 - 33878.5
chr1_bin$foo.X2 <- chr1_bin$foo.X2 - 33878.5
#adding length of bin as column and making in Mb
chr1_bin$length <- (chr1_bin$foo.X2-chr1_bin$foo.X1)/1000000
chr1_bin$rate <- ((chr1_bin$freq/257)*100)/chr1_bin$length

chr2_CO <- chr2_CO[order(chr2_CO$`CO Start`),]
chr2_bin <- binning(chr2_CO$midpoint, nbins = 120, type = "kmeans")
chr2_bin <- as.data.frame(summary(chr2_bin))
chr2_bin <- within(chr2_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr2_bin$levels), ',', fixed=TRUE))))
chr2_bin <- do.call(data.frame, chr2_bin)
chr2_bin <- chr2_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr2_bin <- chr2_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr2_bin[1,4] <- 1968777
chr2_bin$foo.X1 <- chr2_bin$foo.X1 - 1968777
chr2_bin$foo.X2 <- chr2_bin$foo.X2 - 1968777
chr2_bin$length <- (chr2_bin$foo.X2-chr2_bin$foo.X1)/1000000
chr2_bin$rate <- ((chr2_bin$freq/257)*100)/chr2_bin$length

#have not converted rest of chromosomes to what I did in 1 & 2
chr3_CO <- chr3_CO[order(chr3_CO$`CO Start`),]
chr3_bin <- binning(chr3_CO$midpoint, nbins = 120, type = "kmeans")
chr3_bin <- as.data.frame(summary(chr3_bin))
chr3_bin <- within(chr3_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr3_bin$levels), ',', fixed=TRUE))))
chr3_bin <- do.call(data.frame, chr3_bin)
chr3_bin <- chr3_bin %>% mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr3_bin <- chr3_bin %>% mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr3_bin[1,4] <- 142443
chr3_bin$foo.X1 <- chr3_bin$foo.X1 - 142443
chr3_bin$foo.X2 <- chr3_bin$foo.X2 - 142443
chr3_bin$length <- (chr3_bin$foo.X2-chr3_bin$foo.X1)/1000000
chr3_bin$rate <- ((chr3_bin$freq/257)*100)/chr3_bin$length

chr4_CO <- chr4_CO[order(chr4_CO$`CO Start`),]
chr4_bin <- binning(chr4_CO$midpoint, nbins = 100, type = "kmeans")
chr4_bin <- as.data.frame(summary(chr4_bin))
chr4_bin <- within(chr4_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr4_bin$levels), ',', fixed=TRUE))))
chr4_bin <- do.call(data.frame, chr4_bin)
chr4_bin <- chr4_bin %>% mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr4_bin <- chr4_bin %>% mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr4_bin[1,4] <- 2742511
chr4_bin$foo.X1 <- chr4_bin$foo.X1 - 2742511
chr4_bin$foo.X2 <- chr4_bin$foo.X2 - 2742511
chr4_bin$length <- (chr4_bin$foo.X2-chr4_bin$foo.X1)/1000000
chr4_bin$rate <- ((chr4_bin$freq/257)*100)/chr4_bin$length

chr5_CO <- chr5_CO[order(chr5_CO$`CO Start`),]
chr5_bin <- binning(chr5_CO$midpoint, nbins = 80, type = "kmeans")
chr5_bin <- as.data.frame(summary(chr5_bin))
chr5_bin <- within(chr5_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr5_bin$levels), ',', fixed=TRUE))))
chr5_bin <- do.call(data.frame, chr5_bin)
chr5_bin <- chr5_bin %>% mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr5_bin <- chr5_bin %>% mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr5_bin[1,4] <- 3493080
chr5_bin$foo.X1 <- chr5_bin$foo.X1 - 3493080
chr5_bin$foo.X2 <- chr5_bin$foo.X2 - 3493080
chr5_bin$length <- (chr5_bin$foo.X2-chr5_bin$foo.X1)/1000000
chr5_bin$rate <- ((chr5_bin$freq/257)*100)/chr5_bin$length

chr6_CO <- chr6_CO[order(chr6_CO$`CO Start`),]
chr6_bin <- binning(chr6_CO$midpoint, nbins = 75, type = "kmeans")
chr6_bin <- as.data.frame(summary(chr6_bin))
chr6_bin <- within(chr6_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr6_bin$levels), ',', fixed=TRUE))))
chr6_bin <- do.call(data.frame, chr6_bin)
chr6_bin <- chr6_bin %>% mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr6_bin <- chr6_bin %>% mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr6_bin[1,4] <- 1.551805e+07
chr6_bin$foo.X1 <- chr6_bin$foo.X1 - 1.551805e+07
chr6_bin$foo.X2 <- chr6_bin$foo.X2 - 1.551805e+07
chr6_bin$length <- (chr6_bin$foo.X2-chr6_bin$foo.X1)/1000000
chr6_bin$rate <- ((chr6_bin$freq/257)*100)/chr6_bin$length

chr7_CO <- chr7_CO[order(chr7_CO$`CO Start`),]
chr7_bin <- binning(chr7_CO$midpoint, nbins = 75, type = "kmeans")
chr7_bin <- as.data.frame(summary(chr7_bin))
chr7_bin <- within(chr7_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr7_bin$levels), ',', fixed=TRUE))))
chr7_bin <- do.call(data.frame, chr7_bin)
chr7_bin <- chr7_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr7_bin <- chr7_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr7_bin[1,4] <- 1418491
chr7_bin$foo.X1 <- chr7_bin$foo.X1 - 1418491
chr7_bin$foo.X2 <- chr7_bin$foo.X2 - 1418491
chr7_bin$length <- (chr7_bin$foo.X2-chr7_bin$foo.X1)/1000000
chr7_bin$rate <- ((chr7_bin$freq/257)*100)/chr7_bin$length

chr8_CO <- chr8_CO[order(chr8_CO$`CO Start`),]
chr8_bin <- binning(chr8_CO$midpoint, nbins = 75, type = "kmeans")
chr8_bin <- as.data.frame(summary(chr8_bin))
chr8_bin <- within(chr8_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr8_bin$levels), ',', fixed=TRUE))))
chr8_bin <- do.call(data.frame, chr8_bin)
chr8_bin <- chr8_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr8_bin <- chr8_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr8_bin[1,4] <- 1655452
chr8_bin$foo.X1 <- chr8_bin$foo.X1 - 1655452
chr8_bin$foo.X2 <- chr8_bin$foo.X2 - 1655452
chr8_bin$length <- (chr8_bin$foo.X2-chr8_bin$foo.X1)/1000000
chr8_bin$rate <- ((chr8_bin$freq/257)*100)/chr8_bin$length

chr9_CO <- chr9_CO[order(chr9_CO$`CO Start`),]
chr9_bin <- binning(chr9_CO$midpoint, nbins = 70, type = "kmeans")
chr9_bin <- as.data.frame(summary(chr9_bin))
chr9_bin <- within(chr9_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr9_bin$levels), ',', fixed=TRUE))))
chr9_bin <- do.call(data.frame, chr9_bin)
chr9_bin <- chr9_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr9_bin <- chr9_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr9_bin[1,4] <- 22966.5
chr9_bin$foo.X1 <- chr9_bin$foo.X1 - 22966.5
chr9_bin$foo.X2 <- chr9_bin$foo.X2 - 22966.5
chr9_bin$length <- (chr9_bin$foo.X2-chr9_bin$foo.X1)/1000000
chr9_bin$rate <- ((chr9_bin$freq/257)*100)/chr9_bin$length

chr10_CO <- chr10_CO[order(chr10_CO$`CO Start`),]
chr10_bin <- binning(chr10_CO$midpoint, nbins = 70, type = "kmeans")
chr10_bin <- as.data.frame(summary(chr10_bin))
chr10_bin <- within(chr10_bin, foo<-data.frame(do.call('rbind', strsplit(as.character(chr10_bin$levels), ',', fixed=TRUE))))
chr10_bin <- do.call(data.frame, chr10_bin)
chr10_bin <- chr10_bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
chr10_bin <- chr10_bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
chr10_bin[1,4] <- 33839.5
chr10_bin$foo.X1 <- chr10_bin$foo.X1 - 33839.5
chr10_bin$foo.X2 <- chr10_bin$foo.X2 - 33839.5
chr10_bin$length <- (chr10_bin$foo.X2-chr10_bin$foo.X1)/1000000
chr10_bin$rate <- ((chr10_bin$freq/257)*100)/chr10_bin$length
#use plot to look at distribution of k-means

##assigning frequency to SNPs based on recombination frequency in each bin
snp_rate <- function(chr_bin, chr_snp){
  for(i in 1:nrow(chr_snp)){
    for(k in 1:nrow(chr_bin)){
      if((chr_snp$`SNP Start`[i] >= chr_bin$foo.X1[k]) && (chr_snp$`SNP Start`[i] <= chr_bin$foo.X2[k])){
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
#some SNPs were outside of the interval; have to manually apply a recombination freq
chr1_snp2[is.na(chr1_snp2)] <- 0.64581752
#smoothing the recombination rate so transitions between bins are not so abrupt
chr1_spl <- smooth.spline(chr1_snp2$rate, spar = 0.6)
#creation of genetic positions from smoothed recombination rate
chr1_snp2$pos <- (chr1_snp2$`SNP Start`*chr1_spl$y)
#graph to look at Mb vs. cM along chromosome
plot(chr1_snp2$`SNP Start`, chr1_snp2$pos)
ggplot(chr1_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
#graph to look at Mb vs. cM/Mb to see recombination rate along chromosome
plot(chr1_snp2$`SNP Start`, chr1_snp2$pos/chr1_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 1 Recombination Distribution")

chr2_snp2 <- snp_rate(chr2_bin, chr2_snp)
chr2_snp2$`SNP Start` <- chr2_snp2$`SNP Start`/1000000
chr2_snp2[is.na(chr2_snp2)] <- 0.58768322
chr2_spl <- smooth.spline(chr2_snp2$rate, spar = 1)
chr2_snp2$pos <- (chr2_snp2$`SNP Start`*chr2_spl$y)
plot(chr2_snp2$`SNP Start`, chr2_snp2$pos)
plot(chr2_snp2$`SNP Start`, chr2_snp2$pos/chr2_snp2$`SNP Start`, type = "l")

chr3_snp2 <- snp_rate(chr3_bin, chr3_snp)
chr3_snp2$`SNP Start` <- chr3_snp2$`SNP Start`/1000000
chr3_snp2[is.na(chr3_snp2)] <- 0.14400631
chr3_spl <- smooth.spline(chr3_snp2$rate, spar = 0.6)
chr3_snp2$pos <- (chr3_snp2$`SNP Start`*chr3_spl$y)
plot(chr3_snp2$`SNP Start`, chr3_snp2$pos)
plot(chr3_snp2$`SNP Start`, chr3_snp2$pos/chr3_snp2$`SNP Start`, type = "l")

chr4_snp2 <- snp_rate(chr4_bin, chr4_snp)
chr4_snp2$`SNP Start` <- chr4_snp2$`SNP Start`/1000000
chr4_snp2[is.na(chr4_snp2)] <- 1.59884286
chr4_spl <- smooth.spline(chr4_snp2$rate, spar = 0.6)
chr4_snp2$pos <- (chr4_snp2$`SNP Start`*chr4_spl$y)
plot(chr4_snp2$`SNP Start`, chr4_snp2$pos)
plot(chr4_snp2$`SNP Start`, chr4_snp2$pos/chr4_snp2$`SNP Start`, type = "l")

chr5_snp2 <- snp_rate(chr5_bin, chr5_snp)
chr5_snp2$`SNP Start` <- chr5_snp2$`SNP Start`/1000000
chr5_snp2[is.na(chr5_snp2)] <- 1.09145879
chr5_spl <- smooth.spline(chr5_snp2$rate, spar = 0.6)
chr5_snp2$pos <- (chr5_snp2$`SNP Start`*chr5_spl$y)
plot(chr5_snp2$`SNP Start`, chr5_snp2$pos)
plot(chr5_snp2$`SNP Start`, chr5_snp2$pos/chr5_snp2$`SNP Start`, type = "l")

#chr 6 is lowkey fuked up
chr6_snp2 <- snp_rate(chr6_bin, chr6_snp)
chr6_snp2$`SNP Start` <- chr6_snp2$`SNP Start`/1000000
chr6_snp2[is.na(chr6_snp2)] <- 0.011111111
chr6_snp2$pos <- (chr6_snp2$`SNP Start`*chr6_snp2$rate)
plot(chr6_snp2$`SNP Start`, chr6_snp2$pos)
plot(chr6_snp2$`SNP Start`, chr6_snp2$pos/chr6_snp2$`SNP Start`, type = "l")

chr7_snp2 <- snp_rate(chr7_bin, chr7_snp)
chr7_snp2$`SNP Start` <- chr7_snp2$`SNP Start`/1000000
chr7_snp2[is.na(chr7_snp2)] <- 0.041025641
chr7_snp2$pos <- (chr7_snp2$`SNP Start`*chr7_snp2$rate)
plot(chr7_snp2$`SNP Start`, chr7_snp2$pos)
plot(chr7_snp2$`SNP Start`, chr7_snp2$pos/chr7_snp2$`SNP Start`, type = "l")

chr8_snp2 <- snp_rate(chr8_bin, chr8_snp)
chr8_snp2$`SNP Start` <- chr8_snp2$`SNP Start`/1000000
chr8_snp2[is.na(chr8_snp2)] <- 0.005208333
chr8_snp2$pos <- (chr8_snp2$`SNP Start`*chr8_snp2$rate)
plot(chr8_snp2$`SNP Start`, chr8_snp2$pos)
plot(chr8_snp2$`SNP Start`, chr8_snp2$pos/chr8_snp2$`SNP Start`, type = "l")

chr9_snp2 <- snp_rate(chr9_bin, chr9_snp)
chr9_snp2$`SNP Start` <- chr9_snp2$`SNP Start`/1000000
#normalized between rates at 2 ends to replace NAs
chr9_snp2[is.na(chr9_snp2)] <- 0.01569506
chr9_snp2$pos <- (chr9_snp2$`SNP Start`*chr9_snp2$rate)
plot(chr9_snp2$`SNP Start`, chr9_snp2$pos)
plot(chr9_snp2$`SNP Start`, chr9_snp2$pos/chr9_snp2$`SNP Start`, type = "l")

chr10_snp2 <- snp_rate(chr10_bin, chr10_snp)
chr10_snp2$`SNP Start` <- chr10_snp2$`SNP Start`/1000000
chr10_snp2[is.na(chr10_snp2)] <- 0.01569506
chr10_snp2$pos <- (chr10_snp2$`SNP Start`*chr10_snp2$rate)
plot(chr10_snp2$`SNP Start`, chr10_snp2$pos)
plot(chr10_snp2$`SNP Start`, chr10_snp2$pos/chr9_snp2$`SNP Start`, type = "l")

##Reading in haplotype data
haplotypes <- read.table("B73_mol17_haplotypes.vcf", header = FALSE)


###Simulating a realistic breeding program in maize

##Setting up program with altered map 
#need genetic map, haplotypes
#columns of haplotypes corresponds to number of segregating sites & must be equal to length of genetic map
#if genetic map is 100 then there must be 100 columns of haplotype matrix
founderpop <- newMapPop(genMap, haplotypes, inbred = FALSE, ploidy = 2L)

#creating simulation parameters with the founder population
#tracking recombination as well
SP = SimParam$new(founderPop)$setTrackRec(TRUE)
#assuming crossover interference
SP$v = 2.6
#adding an additive trait
SimParam$addTraitA(
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
pop <- newPop(founderpop, simParam = SP)
#find phenotypes for the first population in order to do phenotypic selection
pop <- setPheno(
  pop,
  h2 = NULL,
  H2 = NULL,
  onlyPheno = FALSE,
  simParam = SP
)
#selecting top individuals from first population
#select top 10 using phenotype for one trait
pop_sel <- selectInd(pop, nInd = 10, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)

#put it together to iterate one program 40 times with 20 generations of selection
#after 20 gen of selection, inbred or DH
for(i in 1:40){
  founderpop <- newMapPop(genMap, haplotypes, inbred = FALSE, ploidy = 2L)
  SP = SimParam$new(founderPop)$setTrackRec(TRUE)
  SP$v = 2.6
  SimParam$addTraitA(
    nQtlPerChr,
    mean = 0,
    var = 1,
    corA = NULL,
    gamma = FALSE,
    shape = 1,
    force = FALSE
  )
  SP$setVarE(h2=0.5)
  
  pop <- newPop(founderpop, simParam = SP)
  
  pop_F1 <- randCross(pop, nCrosses = 100, nProgeny = 100, simParam = SP)
  pop_DH <- makeDH(pop_F1, nDH = 1, simParam = NULL)
  pop1 <- setPheno(pop_DH, varE = 1, simParam = SP)
  
  pop1_sel <- selectInd(pop1, nInd = 10, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel_cross <- randCross(pop1_sel, 10, nProgeny = 100, simParam = SP)
  pop1_sel_cross <- setPheno(pop1_sel_cross, simParam = SP)
  
  pop1_sel2 <- selectInd(pop1_sel_cross, nInd = 5, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_cross <- randCross(pop1_sel2, 10, nProgeny = 100, simParam = SP)
  pop1_sel2_cross <- setPheno(pop1_sel2_cross, simParam = SP)
  
  pop1_sel3 <- selectInd(pop1_sel2_cross, nInd = 5, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_cross <- randCross(pop1_sel3, 10, nProgeny = 100, simParam = SP)
  pop1_sel3_cross <- setPheno(pop1_sel3_cross, simParam = SP)
  
  pop1_sel4 <- selectInd(pop1_sel3_cross, nInd = 5, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_cross <- randCross(pop1_sel4, 10, nProgeny = 100, simParam = SP)
  pop1_sel4_cross <- setPheno(pop1_sel4_cross, simParam = SP)
  
  pop1_sel5 <- selectInd(pop1_sel4_cross, nInd = 5, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_cross <- randCross(pop1_sel5, 10, nProgeny = 100, simParam = SP)
  pop1_sel5_cross <- setPheno(pop1_sel5_cross, simParam = SP)
  
  pop1_sel6 <- selectInd(pop1_sel5_cross, nInd = 5, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel6_cross <- randCross(pop1_sel6, 10, nProgeny = 100, simParam = SP)
  pop1_sel6_cross <- setPheno(pop1_sel6_cross, simParam = SP)
  
  pop1_sel7 <- selectInd(pop1_sel6_cross, nInd = 5, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel7_cross <- randCross(pop1_sel7, 10, nProgeny = 100, simParam = SP)
  pop1_sel7_cross <- setPheno(pop1_sel7_cross, simParam = SP)
  
  pop1_sel8 <- selectInd(pop1_sel7_cross, nInd = 5, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel8_cross <- randCross(pop1_sel8, 10, nProgeny = 100, simParam = SP)
  pop1_sel8_cross <- setPheno(pop1_sel8_cross, simParam = SP)
  
  pop1_sel9 <- selectInd(pop1_sel8_cross, nInd = 5, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel9_cross <- randCross(pop1_sel9, 10, nProgeny = 100, simParam = SP)
  pop1_sel9_cross <- setPheno(pop1_sel9_cross, simParam = SP)
  
  pop1_sel10 <- selectInd(pop1_sel9_cross, nInd = 5, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel10_cross <- randCross(pop1_sel10, 10, nProgeny = 100, simParam = SP)
  pop1_sel10_cross <- setPheno(pop1_sel10_cross, simParam = SP)
  
  pop1_sel11 <- selectInd(pop1_sel10_cross, nInd = 5, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel11_cross <- randCross(pop1_sel11, 10, nProgeny = 100, simParam = SP)
  pop1_sel11_cross <- setPheno(pop1_sel11_cross, simParam = SP)
  
  pop1_sel12 <- selectInd(pop1_sel11_cross, nInd = 5, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel12_cross <- randCross(pop1_sel12, 10, nProgeny = 100, simParam = SP)
  pop1_sel12_cross <- setPheno(pop1_sel12_cross, simParam = SP)
  
  pop1_sel13 <- selectInd(pop1_sel12_cross, nInd = 5, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel13_cross <- randCross(pop1_sel13, 10, nProgeny = 100, simParam = SP)
  pop1_sel13_cross <- setPheno(pop1_sel13_cross, simParam = SP)
  
  pop1_sel14 <- selectInd(pop1_sel13_cross, nInd = 5, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel14_cross <- randCross(pop1_sel14, 10, nProgeny = 100, simParam = SP)
  pop1_sel14_cross <- setPheno(pop1_sel14_cross, simParam = SP)
  
  pop1_sel15 <- selectInd(pop1_sel14_cross, nInd = 5, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel15_cross <- randCross(pop1_sel15, 10, nProgeny = 100, simParam = SP)
  pop1_sel15_cross <- setPheno(pop1_sel15_cross, simParam = SP)
  
  pop1_sel16 <- selectInd(pop1_sel15_cross, nInd = 5, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel16_cross <- randCross(pop1_sel16, 10, nProgeny = 100, simParam = SP)
  pop1_sel16_cross <- setPheno(pop1_sel16_cross, simParam = SP)
  
  pop1_sel17 <- selectInd(pop1_sel16_cross, nInd = 5, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel17_cross <- randCross(pop1_sel17, 10, nProgeny = 100, simParam = SP)
  pop1_sel17_cross <- setPheno(pop1_sel17_cross, simParam = SP)
  
  pop1_sel18 <- selectInd(pop1_sel17_cross, nInd = 5, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel18_cross <- randCross(pop1_sel18, 10, nProgeny = 100, simParam = SP)
  pop1_sel18_cross <- setPheno(pop1_sel18_cross, simParam = SP)
  
  pop1_sel19 <- selectInd(pop1_sel18_cross, nInd = 5, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel19_cross <- randCross(pop1_sel19, 10, nProgeny = 100, simParam = SP)
  pop1_sel19_cross <- setPheno(pop1_sel19_cross, simParam = SP)
  
  pop1_sel20 <- selectInd(pop1_sel19_cross, nInd = 5, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel20_cross <- randCross(pop1_sel20, 10, nProgeny = 100, simParam = SP)
  pop1_sel20_cross <- setPheno(pop1_sel20_cross, simParam = SP)
  
  final_sel <- selectInd(pop1_sel20_cross, nInd = 5, use = "pheno", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  finalpop_gv[i] <- gv(final_sel)
}
