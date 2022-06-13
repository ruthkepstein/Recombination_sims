library(AlphaSimR)
library(Rcpp)
library(ggplot2)
library(dplyr)

#setwd("C:/Users/16192/Documents/PNAS_Simulations")
set.seed(420)

#Reading in the SNP dataset from japonica subspecies genome and randomly sampling 2000 of these SNPs
japonica_snps <- read.table("japonica_SNPs.bed", header =FALSE)
ind_snps <- read.table("indica_SNPs.bed", header =FALSE)
colnames(japonica_snps) <- c("Chr#", "SNP Start", "SNP End")
colnames(ind_snps) <- c("Chr#", "SNP Start", "SNP End")
japonica_snps <- sample_n(japonica_snps, 2000)
ind_snps <- sample_n(ind_snps, 2000)
japonica_snps <- japonica_snps[order(japonica_snps$`Chr#`,japonica_snps$`SNP Start`),]
ind_snps <- ind_snps[order(ind_snps$`Chr#`,ind_snps$`SNP Start`),]

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

#recombination rates
japonica_CO <- read.table("japonica_rec_rate.bed", header = FALSE)
colnames(japonica_CO) <- c("Chr", "CO Start", "CO End", "rate")
japonica_CO <- japonica_CO[order(japonica_CO$Chr,japonica_CO$`CO Start`),]

japonica_chr1_CO <- japonica_CO[ which(japonica_CO$Chr == "chr01"),]
japonica_chr1_CO$midpoint <- (japonica_chr1_CO$`CO Start`+ japonica_chr1_CO$`CO End`)/2
japonica_chr1_CO <- japonica_chr1_CO[order(japonica_chr1_CO$`CO Start`),]

japonica_chr2_CO <- japonica_CO[ which(japonica_CO$Chr == "chr02"),]
japonica_chr2_CO$midpoint <- (japonica_chr2_CO$`CO Start`+ japonica_chr2_CO$`CO End`)/2
japonica_chr2_CO <- japonica_chr2_CO[order(japonica_chr2_CO$`CO Start`),]

japonica_chr3_CO <- japonica_CO[ which(japonica_CO$Chr == "chr03"),]
japonica_chr3_CO$midpoint <- (japonica_chr3_CO$`CO Start`+ japonica_chr3_CO$`CO End`)/2
japonica_chr3_CO <- japonica_chr3_CO[order(japonica_chr3_CO$`CO Start`),]

japonica_chr4_CO <- japonica_CO[ which(japonica_CO$Chr == "chr04"),]
japonica_chr4_CO$midpoint <- (japonica_chr4_CO$`CO Start`+ japonica_chr4_CO$`CO End`)/2
japonica_chr4_CO <- japonica_chr4_CO[order(japonica_chr4_CO$`CO Start`),]

japonica_chr5_CO <- japonica_CO[ which(japonica_CO$Chr == "chr05"),]
japonica_chr5_CO$midpoint <- (japonica_chr5_CO$`CO Start`+ japonica_chr5_CO$`CO End`)/2
japonica_chr5_CO <- japonica_chr5_CO[order(japonica_chr5_CO$`CO Start`),]

japonica_chr6_CO <- japonica_CO[ which(japonica_CO$Chr == "chr06"),]
japonica_chr6_CO$midpoint <- (japonica_chr6_CO$`CO Start`+ japonica_chr6_CO$`CO End`)/2
japonica_chr6_CO <- japonica_chr6_CO[order(japonica_chr6_CO$`CO Start`),]

japonica_chr7_CO <- japonica_CO[ which(japonica_CO$Chr == "chr07"),]
japonica_chr7_CO$midpoint <- (japonica_chr7_CO$`CO Start`+ japonica_chr7_CO$`CO End`)/2
japonica_chr7_CO <- japonica_chr7_CO[order(japonica_chr7_CO$`CO Start`),]

japonica_chr8_CO <- japonica_CO[ which(japonica_CO$Chr == "chr08"),]
japonica_chr8_CO$midpoint <- (japonica_chr8_CO$`CO Start`+ japonica_chr8_CO$`CO End`)/2
japonica_chr8_CO <- japonica_chr8_CO[order(japonica_chr8_CO$`CO Start`),]

japonica_chr9_CO <- japonica_CO[ which(japonica_CO$Chr == "chr09"),]
japonica_chr9_CO$midpoint <- (japonica_chr9_CO$`CO Start`+ japonica_chr9_CO$`CO End`)/2
japonica_chr9_CO <- japonica_chr9_CO[order(japonica_chr9_CO$`CO Start`),]

japonica_chr10_CO <- japonica_CO[ which(japonica_CO$Chr == "chr10"),]
japonica_chr10_CO$midpoint <- (japonica_chr10_CO$`CO Start`+ japonica_chr10_CO$`CO End`)/2
japonica_chr10_CO <- japonica_chr10_CO[order(japonica_chr10_CO$`CO Start`),]

japonica_chr11_CO <- japonica_CO[ which(japonica_CO$Chr == "chr11"),]
japonica_chr11_CO$midpoint <- (japonica_chr11_CO$`CO Start`+ japonica_chr11_CO$`CO End`)/2
japonica_chr11_CO <- japonica_chr11_CO[order(japonica_chr11_CO$`CO Start`),]

japonica_chr12_CO <- japonica_CO[ which(japonica_CO$Chr == "chr12"),]
japonica_chr12_CO$midpoint <- (japonica_chr12_CO$`CO Start`+ japonica_chr12_CO$`CO End`)/2
japonica_chr12_CO <- japonica_chr12_CO[order(japonica_chr12_CO$`CO Start`),]


ind_CO <- read.table("indica_rec_rates.bed", header = FALSE)
colnames(ind_CO) <- c("Chr", "CO Start", "CO End", "rate")
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

#making intervals start at 0
japonica_chr1_CO$`CO End` <- japonica_chr1_CO$`CO End` - min(japonica_chr1_CO$`CO Start`)
japonica_chr1_CO$`CO Start` <- japonica_chr1_CO$`CO Start` - min(japonica_chr1_CO$`CO Start`)

japonica_chr2_CO$`CO End` <- japonica_chr2_CO$`CO End` - min(japonica_chr2_CO$`CO Start`)
japonica_chr2_CO$`CO Start` <- japonica_chr2_CO$`CO Start` - min(japonica_chr2_CO$`CO Start`)

japonica_chr3_CO$`CO End` <- japonica_chr3_CO$`CO End` - min(japonica_chr3_CO$`CO Start`)
japonica_chr3_CO$`CO Start` <- japonica_chr3_CO$`CO Start` - min(japonica_chr3_CO$`CO Start`)

japonica_chr4_CO$`CO End` <- japonica_chr4_CO$`CO End` - min(japonica_chr4_CO$`CO Start`)
japonica_chr4_CO$`CO Start` <- japonica_chr4_CO$`CO Start` - min(japonica_chr4_CO$`CO Start`)

japonica_chr5_CO$`CO End` <- japonica_chr5_CO$`CO End` - min(japonica_chr5_CO$`CO Start`)
japonica_chr5_CO$`CO Start` <- japonica_chr5_CO$`CO Start` - min(japonica_chr5_CO$`CO Start`)

japonica_chr6_CO$`CO End` <- japonica_chr6_CO$`CO End` - min(japonica_chr6_CO$`CO Start`)
japonica_chr6_CO$`CO Start` <- japonica_chr6_CO$`CO Start` - min(japonica_chr6_CO$`CO Start`)

japonica_chr7_CO$`CO End` <- japonica_chr7_CO$`CO End` - min(japonica_chr7_CO$`CO Start`)
japonica_chr7_CO$`CO Start` <- japonica_chr7_CO$`CO Start` - min(japonica_chr7_CO$`CO Start`)

japonica_chr8_CO$`CO End` <- japonica_chr8_CO$`CO End` - min(japonica_chr8_CO$`CO Start`)
japonica_chr8_CO$`CO Start` <- japonica_chr8_CO$`CO Start` - min(japonica_chr8_CO$`CO Start`)

japonica_chr9_CO$`CO End` <- japonica_chr9_CO$`CO End` - min(japonica_chr9_CO$`CO Start`)
japonica_chr9_CO$`CO Start` <- japonica_chr9_CO$`CO Start` - min(japonica_chr9_CO$`CO Start`)

japonica_chr10_CO$`CO End` <- japonica_chr10_CO$`CO End` - min(japonica_chr10_CO$`CO Start`)
japonica_chr10_CO$`CO Start` <- japonica_chr10_CO$`CO Start` - min(japonica_chr10_CO$`CO Start`)

japonica_chr11_CO$`CO End` <- japonica_chr11_CO$`CO End` - min(japonica_chr11_CO$`CO Start`)
japonica_chr11_CO$`CO Start` <- japonica_chr11_CO$`CO Start` - min(japonica_chr11_CO$`CO Start`)

japonica_chr12_CO$`CO End` <- japonica_chr12_CO$`CO End` - min(japonica_chr12_CO$`CO Start`)
japonica_chr12_CO$`CO Start` <- japonica_chr12_CO$`CO Start` - min(japonica_chr12_CO$`CO Start`)

ind_chr1_CO$`CO End` <- ind_chr1_CO$`CO End` - min(ind_chr1_CO$`CO Start`)
ind_chr1_CO$`CO Start` <- ind_chr1_CO$`CO Start` - min(ind_chr1_CO$`CO Start`)

ind_chr2_CO$`CO End` <- ind_chr2_CO$`CO End` - min(ind_chr2_CO$`CO Start`)
ind_chr2_CO$`CO Start` <- ind_chr2_CO$`CO Start` - min(ind_chr2_CO$`CO Start`)

ind_chr3_CO$`CO End` <- ind_chr3_CO$`CO End` - min(ind_chr3_CO$`CO Start`)
ind_chr3_CO$`CO Start` <- ind_chr3_CO$`CO Start` - min(ind_chr3_CO$`CO Start`)

ind_chr4_CO$`CO End` <- ind_chr4_CO$`CO End` - min(ind_chr4_CO$`CO Start`)
ind_chr4_CO$`CO Start` <- ind_chr4_CO$`CO Start` - min(ind_chr4_CO$`CO Start`)

ind_chr5_CO$`CO End` <- ind_chr5_CO$`CO End` - min(ind_chr5_CO$`CO Start`)
ind_chr5_CO$`CO Start` <- ind_chr5_CO$`CO Start` - min(ind_chr5_CO$`CO Start`)

ind_chr6_CO$`CO End` <- ind_chr6_CO$`CO End` - min(ind_chr6_CO$`CO Start`)
ind_chr6_CO$`CO Start` <- ind_chr6_CO$`CO Start` - min(ind_chr6_CO$`CO Start`)

ind_chr7_CO$`CO End` <- ind_chr7_CO$`CO End` - min(ind_chr7_CO$`CO Start`)
ind_chr7_CO$`CO Start` <- ind_chr7_CO$`CO Start` - min(ind_chr7_CO$`CO Start`)

ind_chr8_CO$`CO End` <- ind_chr8_CO$`CO End` - min(ind_chr8_CO$`CO Start`)
ind_chr8_CO$`CO Start` <- ind_chr8_CO$`CO Start` - min(ind_chr8_CO$`CO Start`)

ind_chr9_CO$`CO End` <- ind_chr9_CO$`CO End` - min(ind_chr9_CO$`CO Start`)
ind_chr9_CO$`CO Start` <- ind_chr9_CO$`CO Start` - min(ind_chr9_CO$`CO Start`)

ind_chr10_CO$`CO End` <- ind_chr10_CO$`CO End` - min(ind_chr10_CO$`CO Start`)
ind_chr10_CO$`CO Start` <- ind_chr10_CO$`CO Start` - min(ind_chr10_CO$`CO Start`)

ind_chr11_CO$`CO End` <- ind_chr11_CO$`CO End` - min(ind_chr11_CO$`CO Start`)
ind_chr11_CO$`CO Start` <- ind_chr11_CO$`CO Start` - min(ind_chr11_CO$`CO Start`)

ind_chr12_CO$`CO End` <- ind_chr12_CO$`CO End` - min(ind_chr12_CO$`CO Start`)
ind_chr12_CO$`CO Start` <- ind_chr12_CO$`CO Start` - min(ind_chr12_CO$`CO Start`)

isTRUE(ind_chr1_CO[1,2] == 0)
isTRUE(ind_chr2_CO[1,2] == 0)
isTRUE(ind_chr3_CO[1,2] == 0)
isTRUE(ind_chr4_CO[1,2] == 0)
isTRUE(ind_chr5_CO[1,2] == 0)
isTRUE(ind_chr6_CO[1,2] == 0)
isTRUE(ind_chr7_CO[1,2] == 0)
isTRUE(ind_chr8_CO[1,2] == 0)
isTRUE(ind_chr9_CO[1,2] == 0)
isTRUE(ind_chr10_CO[1,2] == 0)
isTRUE(ind_chr11_CO[1,2] == 0)
isTRUE(ind_chr12_CO[1,2] == 0)

isTRUE(japonica_chr1_CO[1,2] == 0)
isTRUE(japonica_chr2_CO[1,2] == 0)
isTRUE(japonica_chr3_CO[1,2] == 0)
isTRUE(japonica_chr4_CO[1,2] == 0)
isTRUE(japonica_chr5_CO[1,2] == 0)
isTRUE(japonica_chr6_CO[1,2] == 0)
isTRUE(japonica_chr7_CO[1,2] == 0)
isTRUE(japonica_chr8_CO[1,2] == 0)
isTRUE(japonica_chr9_CO[1,2] == 0)
isTRUE(japonica_chr10_CO[1,2] == 0)
isTRUE(japonica_chr11_CO[1,2] == 0)
isTRUE(japonica_chr12_CO[1,2] == 0)

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
japonica_chr1_CO_2 <- japonica_chr1_CO
bins<-as.integer(nrow(japonica_chr1_CO)/40)
japonica_chr1_CO_2$rates<- rollapply(japonica_chr1_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr1_CO_2<-fill_start(japonica_chr1_CO_2)
japonica_chr1_CO_2<- japonica_chr1_CO_2 %>% drop_na(rates)

japonica_chr2_CO_2 <- japonica_chr2_CO
bins<-as.integer(nrow(japonica_chr2_CO)/40)
japonica_chr2_CO_2$rates<- rollapply(japonica_chr2_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr2_CO_2<-fill_start(japonica_chr2_CO_2)
japonica_chr2_CO_2<- japonica_chr2_CO_2 %>% drop_na(rates)

japonica_chr3_CO_2 <- japonica_chr3_CO
bins<-as.integer(nrow(japonica_chr3_CO)/40)
japonica_chr3_CO_2$rates<- rollapply(japonica_chr3_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr3_CO_2<-fill_start(japonica_chr3_CO_2)
japonica_chr3_CO_2<- japonica_chr3_CO_2 %>% drop_na(rates)

japonica_chr4_CO_2 <- japonica_chr4_CO
bins<-as.integer(nrow(japonica_chr4_CO)/40)
japonica_chr4_CO_2$rates<- rollapply(japonica_chr4_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr4_CO_2<-fill_start(japonica_chr4_CO_2)
japonica_chr4_CO_2<- japonica_chr4_CO_2 %>% drop_na(rates)

japonica_chr5_CO_2 <- japonica_chr5_CO
bins<-as.integer(nrow(japonica_chr5_CO)/40)
japonica_chr5_CO_2$rates<- rollapply(japonica_chr5_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr5_CO_2<-fill_start(japonica_chr5_CO_2)
japonica_chr5_CO_2<- japonica_chr5_CO_2 %>% drop_na(rates)

japonica_chr6_CO_2 <- japonica_chr6_CO
bins<-as.integer(nrow(japonica_chr6_CO)/40)
japonica_chr6_CO_2$rates<- rollapply(japonica_chr6_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr6_CO_2<-fill_start(japonica_chr6_CO_2)
japonica_chr6_CO_2<- japonica_chr6_CO_2 %>% drop_na(rates)

japonica_chr7_CO_2 <- japonica_chr7_CO
bins<-as.integer(nrow(japonica_chr7_CO)/40)
japonica_chr7_CO_2$rates<- rollapply(japonica_chr7_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr7_CO_2<-fill_start(japonica_chr7_CO_2)
japonica_chr7_CO_2<- japonica_chr7_CO_2 %>% drop_na(rates)

japonica_chr8_CO_2 <- japonica_chr8_CO
bins<-as.integer(nrow(japonica_chr8_CO)/40)
japonica_chr8_CO_2$rates<- rollapply(japonica_chr8_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr8_CO_2<-fill_start(japonica_chr8_CO_2)
japonica_chr8_CO_2<- japonica_chr8_CO_2 %>% drop_na(rates)

japonica_chr9_CO_2 <- japonica_chr9_CO
bins<-as.integer(nrow(japonica_chr9_CO)/40)
japonica_chr9_CO_2$rates<- rollapply(japonica_chr9_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr9_CO_2<-fill_start(japonica_chr9_CO_2)
japonica_chr9_CO_2<- japonica_chr9_CO_2 %>% drop_na(rates)

japonica_chr10_CO_2 <- japonica_chr10_CO
bins<-as.integer(nrow(japonica_chr10_CO)/40)
japonica_chr10_CO_2$rates<- rollapply(japonica_chr10_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr10_CO_2<-fill_start(japonica_chr10_CO_2)
japonica_chr10_CO_2<- japonica_chr10_CO_2 %>% drop_na(rates)

japonica_chr11_CO_2 <- japonica_chr11_CO
bins<-as.integer(nrow(japonica_chr11_CO)/40)
japonica_chr11_CO_2$rates<- rollapply(japonica_chr11_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr11_CO_2<-fill_start(japonica_chr11_CO_2)
japonica_chr11_CO_2<- japonica_chr11_CO_2 %>% drop_na(rates)

japonica_chr12_CO_2 <- japonica_chr12_CO
bins<-as.integer(nrow(japonica_chr12_CO)/40)
japonica_chr12_CO_2$rates<- rollapply(japonica_chr12_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr12_CO_2<-fill_start(japonica_chr12_CO_2)
japonica_chr12_CO_2<- japonica_chr12_CO_2 %>% drop_na(rates)


ind_chr1_CO_2 <- ind_chr1_CO
bins<-as.integer(nrow(ind_chr1_CO)/40)
ind_chr1_CO_2$rates<- rollapply(ind_chr1_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ind_chr1_CO_2<-fill_start(ind_chr1_CO_2)
ind_chr1_CO_2<- ind_chr1_CO_2 %>% drop_na(rates)

ind_chr2_CO_2 <- ind_chr2_CO
bins<-as.integer(nrow(ind_chr2_CO)/40)
ind_chr2_CO_2$rates<- rollapply(ind_chr2_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ind_chr2_CO_2<-fill_start(ind_chr2_CO_2)
ind_chr2_CO_2<- ind_chr2_CO_2 %>% drop_na(rates)

ind_chr3_CO_2 <- ind_chr3_CO
bins<-as.integer(nrow(ind_chr3_CO)/40)
ind_chr3_CO_2$rates<- rollapply(ind_chr3_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ind_chr3_CO_2<-fill_start(ind_chr3_CO_2)
ind_chr3_CO_2<- ind_chr3_CO_2 %>% drop_na(rates)

ind_chr4_CO_2 <- ind_chr4_CO
bins<-as.integer(nrow(ind_chr4_CO)/40)
ind_chr4_CO_2$rates<- rollapply(ind_chr4_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ind_chr4_CO_2<-fill_start(ind_chr4_CO_2)
ind_chr4_CO_2<- ind_chr4_CO_2 %>% drop_na(rates)

ind_chr5_CO_2 <- ind_chr5_CO
bins<-as.integer(nrow(ind_chr5_CO)/40)
ind_chr5_CO_2$rates<- rollapply(ind_chr5_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ind_chr5_CO_2<-fill_start(ind_chr5_CO_2)
ind_chr5_CO_2<- ind_chr5_CO_2 %>% drop_na(rates)

ind_chr6_CO_2 <- ind_chr6_CO
bins<-as.integer(nrow(ind_chr6_CO)/40)
ind_chr6_CO_2$rates<- rollapply(ind_chr6_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ind_chr6_CO_2<-fill_start(ind_chr6_CO_2)
ind_chr6_CO_2<- ind_chr6_CO_2 %>% drop_na(rates)

ind_chr7_CO_2 <- ind_chr7_CO
bins<-as.integer(nrow(ind_chr7_CO)/40)
ind_chr7_CO_2$rates<- rollapply(ind_chr7_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ind_chr7_CO_2<-fill_start(ind_chr7_CO_2)
ind_chr7_CO_2<- ind_chr7_CO_2 %>% drop_na(rates)

ind_chr8_CO_2 <- ind_chr8_CO
bins<-as.integer(nrow(ind_chr8_CO)/40)
ind_chr8_CO_2$rates<- rollapply(ind_chr8_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ind_chr8_CO_2<-fill_start(ind_chr8_CO_2)
ind_chr8_CO_2<- ind_chr8_CO_2 %>% drop_na(rates)

ind_chr9_CO_2 <- ind_chr9_CO
bins<-as.integer(nrow(ind_chr9_CO)/40)
ind_chr9_CO_2$rates<- rollapply(ind_chr9_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ind_chr9_CO_2<-fill_start(ind_chr9_CO_2)
ind_chr9_CO_2<- ind_chr9_CO_2 %>% drop_na(rates)

ind_chr10_CO_2 <- ind_chr10_CO
bins<-as.integer(nrow(ind_chr10_CO)/40)
ind_chr10_CO_2$rates<- rollapply(ind_chr10_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ind_chr10_CO_2<-fill_start(ind_chr10_CO_2)
ind_chr10_CO_2<- ind_chr10_CO_2 %>% drop_na(rates)

ind_chr11_CO_2 <- ind_chr11_CO
bins<-as.integer(nrow(ind_chr11_CO)/40)
ind_chr11_CO_2$rates<- rollapply(ind_chr11_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ind_chr11_CO_2<-fill_start(ind_chr11_CO_2)
ind_chr11_CO_2<- ind_chr11_CO_2 %>% drop_na(rates)

ind_chr12_CO_2 <- ind_chr12_CO
bins<-as.integer(nrow(ind_chr12_CO)/40)
ind_chr12_CO_2$rates<- rollapply(ind_chr12_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ind_chr12_CO_2<-fill_start(ind_chr12_CO_2)
ind_chr12_CO_2<- ind_chr12_CO_2 %>% drop_na(rates)


##assigning frequency to SNPs based on avg recombination rate in each bin
snp_rate <- function(chr_rate, chr_snp){
  for(i in 1:nrow(chr_snp)){
    for(k in 1:nrow(chr_rate)){
      if(isTRUE((chr_snp$`SNP Start`[i] >= chr_rate$`CO Start`[k]) && (chr_snp$`SNP End`[i] <= chr_rate$`CO End`[k]))){
        chr_snp$rate[i] <- chr_rate$rates[k]
      }
    }
  }
  print(chr_snp)
}

##JAPONICA
#using function, converted SNP start to Mb to get cM/Mb for final genetic position - assign rates
japonica_chr1_snp2 <- snp_rate(japonica_chr1_CO_2, japonica_chr1_snp)
japonica_chr2_snp2 <- snp_rate(japonica_chr2_CO_2, japonica_chr2_snp)
japonica_chr3_snp2 <- snp_rate(japonica_chr3_CO_2, japonica_chr3_snp)
japonica_chr4_snp2 <- snp_rate(japonica_chr4_CO_2, japonica_chr4_snp)
japonica_chr5_snp2 <- snp_rate(japonica_chr5_CO_2, japonica_chr5_snp)
japonica_chr6_snp2 <- snp_rate(japonica_chr6_CO_2, japonica_chr6_snp)
japonica_chr7_snp2 <- snp_rate(japonica_chr7_CO_2, japonica_chr7_snp)
japonica_chr8_snp2 <- snp_rate(japonica_chr8_CO_2, japonica_chr8_snp)
japonica_chr9_snp2 <- snp_rate(japonica_chr9_CO_2, japonica_chr9_snp)
japonica_chr10_snp2 <- snp_rate(japonica_chr10_CO_2, japonica_chr10_snp)
japonica_chr11_snp2 <- snp_rate(japonica_chr11_CO_2, japonica_chr11_snp)
japonica_chr12_snp2 <- snp_rate(japonica_chr12_CO_2, japonica_chr12_snp)

#make mutable copies of data
japonica_chr1_snp3 <- japonica_chr1_snp2
japonica_chr2_snp3 <- japonica_chr2_snp2
japonica_chr3_snp3 <- japonica_chr3_snp2
japonica_chr4_snp3 <- japonica_chr4_snp2
japonica_chr5_snp3 <- japonica_chr5_snp2
japonica_chr6_snp3 <- japonica_chr6_snp2
japonica_chr7_snp3 <- japonica_chr7_snp2
japonica_chr8_snp3 <- japonica_chr8_snp2
japonica_chr9_snp3 <- japonica_chr9_snp2
japonica_chr10_snp3 <- japonica_chr10_snp2
japonica_chr11_snp3 <- japonica_chr11_snp2
japonica_chr12_snp3 <- japonica_chr12_snp2

#adjust start
japonica_chr1_snp3$`SNP Start`<- japonica_chr1_snp3$`SNP Start`/1000000
japonica_chr2_snp3$`SNP Start` <- japonica_chr2_snp3$`SNP Start`/1000000
japonica_chr3_snp3$`SNP Start` <- japonica_chr3_snp3$`SNP Start`/1000000
japonica_chr4_snp3$`SNP Start` <- japonica_chr4_snp3$`SNP Start`/1000000
japonica_chr5_snp3$`SNP Start` <- japonica_chr5_snp3$`SNP Start`/1000000
japonica_chr6_snp3$`SNP Start` <- japonica_chr6_snp3$`SNP Start`/1000000
japonica_chr7_snp3$`SNP Start` <- japonica_chr7_snp3$`SNP Start`/1000000
japonica_chr8_snp3$`SNP Start` <- japonica_chr8_snp3$`SNP Start`/1000000
japonica_chr9_snp3$`SNP Start` <- japonica_chr9_snp3$`SNP Start`/1000000
japonica_chr10_snp3$`SNP Start` <- japonica_chr10_snp3$`SNP Start`/1000000
japonica_chr11_snp3$`SNP Start` <- japonica_chr11_snp3$`SNP Start`/1000000
japonica_chr12_snp3$`SNP Start` <- japonica_chr12_snp3$`SNP Start`/1000000

fill_NA<-function(SNP){
  for(i in 1:nrow(SNP)){
    if(is.na(SNP$rate[i])){
      SNP$rate[i]<-SNP$rate[i-1]
    }
  }
  print(SNP)
}
japonica_chr1_snp3<-fill_NA(japonica_chr1_snp3)
japonica_chr2_snp3<-fill_NA(japonica_chr2_snp3)
japonica_chr3_snp3<-fill_NA(japonica_chr3_snp3)
japonica_chr4_snp3<-fill_NA(japonica_chr4_snp3)
japonica_chr5_snp3<-fill_NA(japonica_chr5_snp3)
japonica_chr6_snp3<-fill_NA(japonica_chr6_snp3)
japonica_chr7_snp3<-fill_NA(japonica_chr7_snp3)
japonica_chr8_snp3<-fill_NA(japonica_chr8_snp3)
japonica_chr9_snp3<-fill_NA(japonica_chr9_snp3)
japonica_chr10_snp3<-fill_NA(japonica_chr10_snp3)
japonica_chr11_snp3<-fill_NA(japonica_chr11_snp3)
japonica_chr12_snp3<-fill_NA(japonica_chr12_snp3)

gen_pos <- function(SNP){
  SNP$pos <- NA
  SNP$pos[1]<-SNP$`SNP Start`[1]*SNP$rate[1]
  for(i in 1:nrow(SNP)){
    if(i>1){
      #SNP$pos[i]<-SNP$`SNP Start`[i]*spl$y[i] + SNP$`SNP Start`[i-1]
      SNP$pos[i]<- SNP$pos[i-1] + (SNP$`SNP Start`[i] - SNP$`SNP Start`[i-1])*SNP$rate[i]
    }
  }
  print(SNP$pos)
}

graph_recomb <- function(SNP){
  SNP$pos2[1]<-SNP$`SNP Start`[1]*SNP$rate[1]
  for(i in 1:nrow(SNP)){
    if(i>1){
      SNP$pos2[i]<- (SNP$pos[i] - SNP$pos[i-1])/ (SNP$`SNP Start`[i] - SNP$`SNP Start`[i-1])
    }
  }
  print(SNP$pos2)
}

japonica_chr1_spl <- smooth.spline(japonica_chr1_snp3$rate, spar = 0)
japonica_chr1_snp3$pos <- gen_pos(japonica_chr1_snp3)
plot(japonica_chr1_snp3$`SNP Start`, japonica_chr1_snp3$rate, type = "l")
ggplot(japonica_chr1_snp3, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
japonica_chr1_finalpos <- japonica_chr1_snp3[order(japonica_chr1_snp3$pos),]
is.unsorted(japonica_chr1_finalpos$pos)
rates <-graph_recomb(japonica_chr1_snp3)
plot(japonica_chr1_snp3$`SNP Start`, rates, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 1 Recombination Distribution")
plot(japonica_chr1_snp3$`SNP Start`, japonica_chr1_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Chromosome 1 Genetic Map")
plot(japonica_chr1_finalpos$`SNP Start`, japonica_chr1_finalpos$pos)


japonica_chr2_spl <- smooth.spline(japonica_chr2_snp3$rate, spar = 0.1)
japonica_chr2_snp3$pos <- gen_pos(japonica_chr2_snp3)
plot(japonica_chr2_snp3$`SNP Start`, japonica_chr2_snp3$pos)
plot(japonica_chr2_snp3$`SNP Start`, japonica_chr2_snp3$pos/japonica_chr2_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Chromosome 2 Recombination Distribution")
japonica_chr2_finalpos <- japonica_chr2_snp3[order(japonica_chr2_snp3$pos),]
is.unsorted(japonica_chr2_finalpos$pos)
plot(japonica_chr2_snp3$`SNP Start`, japonica_chr2_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Chromosome 2 Genetic Map")

japonica_chr3_spl <- smooth.spline(japonica_chr3_snp3$rate, spar = 0)
japonica_chr3_snp3$pos <- gen_pos(japonica_chr3_snp3)
plot(japonica_chr3_snp3$`SNP Start`, japonica_chr3_snp3$pos)
plot(japonica_chr3_snp3$`SNP Start`, japonica_chr3_snp3$pos/japonica_chr3_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Chromosome 3 Recombination Distribution")
japonica_chr3_finalpos <- japonica_chr3_snp3[order(japonica_chr3_snp3$pos),]
is.unsorted(japonica_chr3_finalpos$pos)
plot(japonica_chr3_snp3$`SNP Start`, japonica_chr3_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Chromosome 3 Genetic Map")

japonica_chr4_spl <- smooth.spline(japonica_chr4_snp3$rate, spar = 0)
japonica_chr4_snp3$pos <- gen_pos(japonica_chr4_snp3)
plot(japonica_chr4_snp3$`SNP Start`, japonica_chr4_snp3$pos)
plot(japonica_chr4_snp3$`SNP Start`, japonica_chr4_snp3$pos/japonica_chr4_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Chromosome 4 Recombination Distribution")
japonica_chr4_finalpos <- japonica_chr4_snp3[order(japonica_chr4_snp3$pos),]
is.unsorted(japonica_chr4_finalpos$pos)
plot(japonica_chr4_snp3$`SNP Start`, japonica_chr4_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Chromosome 4 Genetic Map")

japonica_chr5_spl <- smooth.spline(japonica_chr5_snp3$rate, spar =0)
japonica_chr5_snp3$pos <- gen_pos(japonica_chr5_snp3)
plot(japonica_chr5_snp3$`SNP Start`, japonica_chr5_snp3$pos)
plot(japonica_chr5_snp3$`SNP Start`, japonica_chr5_snp3$pos/japonica_chr5_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Chromosome 5 Recombination Distribution")
japonica_chr5_finalpos <- japonica_chr5_snp3[order(japonica_chr5_snp3$pos),]
is.unsorted(japonica_chr5_finalpos$pos)
plot(japonica_chr5_snp3$`SNP Start`, japonica_chr5_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Chromosome 5 Genetic Map")

japonica_chr6_spl <- smooth.spline(japonica_chr6_snp3$rate, spar =0)
japonica_chr6_snp3$pos <- gen_pos(japonica_chr6_snp3)
plot(japonica_chr6_snp3$`SNP Start`, japonica_chr6_snp3$pos)
plot(japonica_chr6_snp3$`SNP Start`, japonica_chr6_snp3$pos/japonica_chr6_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Chromosome 6 Recombination Distribution")
japonica_chr6_finalpos <- japonica_chr6_snp3[order(japonica_chr6_snp3$pos),]
is.unsorted(japonica_chr6_finalpos$pos)
plot(japonica_chr6_snp3$`SNP Start`, japonica_chr6_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Chromosome 6 Genetic Map")

japonica_chr7_spl <- smooth.spline(japonica_chr7_snp3$rate, spar = 0)
japonica_chr7_snp3$pos <- gen_pos(japonica_chr7_snp3)
plot(japonica_chr7_snp3$`SNP Start`, japonica_chr7_snp3$pos)
plot(japonica_chr7_snp3$`SNP Start`, japonica_chr7_snp3$pos/japonica_chr7_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Chromosome 7 Recombination Distribution")
japonica_chr7_finalpos <- japonica_chr7_snp3[order(japonica_chr7_snp3$pos),]
is.unsorted(japonica_chr7_finalpos$pos)
plot(japonica_chr7_snp3$`SNP Start`, japonica_chr7_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Chromosome 7 Genetic Map")

japonica_chr8_spl <- smooth.spline(japonica_chr8_snp3$rate, spar = 0)
japonica_chr8_snp3$pos <- gen_pos(japonica_chr8_snp3)
plot(japonica_chr8_snp3$`SNP Start`, japonica_chr8_snp3$pos)
plot(japonica_chr8_snp3$`SNP Start`, japonica_chr8_snp3$pos/japonica_chr8_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Chromosome 8 Recombination Distribution")
japonica_chr8_finalpos <- japonica_chr8_snp3[order(japonica_chr8_snp3$pos),]
is.unsorted(japonica_chr8_finalpos$pos)
plot(japonica_chr8_snp3$`SNP Start`, japonica_chr8_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Chromosome 8 Genetic Map")

japonica_chr9_spl <- smooth.spline(japonica_chr9_snp3$rate, spar = 0)
japonica_chr9_snp3$pos <- gen_pos(japonica_chr9_snp3)
plot(japonica_chr9_snp3$`SNP Start`, japonica_chr9_snp3$pos)
plot(japonica_chr9_snp3$`SNP Start`, japonica_chr9_snp3$pos/japonica_chr9_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Chromosome 9 Recombination Distribution")
japonica_chr9_finalpos <- japonica_chr9_snp3[order(japonica_chr9_snp3$pos),]
is.unsorted(japonica_chr9_finalpos$pos)
plot(japonica_chr9_snp3$`SNP Start`, japonica_chr9_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Chromosome 9 Genetic Map")

japonica_chr10_spl <- smooth.spline(japonica_chr10_snp3$rate, spar =0)
japonica_chr10_snp3$pos <- gen_pos(japonica_chr10_snp3)
plot(japonica_chr10_snp3$`SNP Start`, japonica_chr10_snp3$pos)
plot(japonica_chr10_snp3$`SNP Start`, japonica_chr10_snp3$pos/japonica_chr10_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Chromosome 10 Recombination Distribution")
japonica_chr10_finalpos <- japonica_chr10_snp3[order(japonica_chr10_snp3$pos),]
is.unsorted(japonica_chr10_finalpos$pos)
plot(japonica_chr10_snp3$`SNP Start`, japonica_chr10_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Chromosome 10 Genetic Map")

japonica_chr11_spl <- smooth.spline(japonica_chr11_snp3$rate, spar = 0)
japonica_chr11_snp3$pos <- gen_pos(japonica_chr11_snp3)
plot(japonica_chr11_snp3$`SNP Start`, japonica_chr11_snp3$pos)
plot(japonica_chr11_snp3$`SNP Start`, japonica_chr11_snp3$pos/japonica_chr11_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Chromosome 11 Recombination Distribution")
japonica_chr11_finalpos <- japonica_chr11_snp3[order(japonica_chr11_snp3$pos),]
is.unsorted(japonica_chr11_finalpos$pos)
plot(japonica_chr11_snp3$`SNP Start`, japonica_chr11_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Chromosome 11 Genetic Map")

japonica_chr12_spl <- smooth.spline(japonica_chr12_snp3$rate, spar = 0)
japonica_chr12_snp3$pos <- gen_pos(japonica_chr12_snp3)
plot(japonica_chr12_snp3$`SNP Start`, japonica_chr12_snp3$pos)
plot(japonica_chr12_snp3$`SNP Start`, japonica_chr12_snp3$pos/japonica_chr12_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Chromosome 12 Recombination Distribution")
japonica_chr12_finalpos <- japonica_chr12_snp3[order(japonica_chr12_snp3$pos),]
is.unsorted(japonica_chr12_finalpos$pos)
plot(japonica_chr12_snp3$`SNP Start`, japonica_chr12_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Chromosome 12 Genetic Map")

##INDICA 
#assign recomb rates
ind_chr1_snp2 <- snp_rate(ind_chr1_CO_2, ind_chr1_snp)
ind_chr2_snp2 <- snp_rate(ind_chr2_CO_2, ind_chr2_snp)
ind_chr3_snp2 <- snp_rate(ind_chr3_CO_2, ind_chr3_snp)
ind_chr4_snp2 <- snp_rate(ind_chr4_CO_2, ind_chr4_snp)
ind_chr5_snp2 <- snp_rate(ind_chr5_CO_2, ind_chr5_snp)
ind_chr6_snp2 <- snp_rate(ind_chr6_CO_2, ind_chr6_snp)
ind_chr7_snp2 <- snp_rate(ind_chr7_CO_2, ind_chr7_snp)
ind_chr8_snp2 <- snp_rate(ind_chr8_CO_2, ind_chr8_snp)
ind_chr9_snp2 <- snp_rate(ind_chr9_CO_2, ind_chr9_snp)
ind_chr10_snp2 <- snp_rate(ind_chr10_CO_2, ind_chr10_snp)
ind_chr11_snp2 <- snp_rate(ind_chr11_CO_2, ind_chr11_snp)
ind_chr12_snp2 <- snp_rate(ind_chr12_CO_2, ind_chr12_snp)

#make mutable copies of data
ind_chr1_snp3 <- ind_chr1_snp2
ind_chr2_snp3 <- ind_chr2_snp2
ind_chr3_snp3 <- ind_chr3_snp2
ind_chr4_snp3 <- ind_chr4_snp2
ind_chr5_snp3 <- ind_chr5_snp2
ind_chr6_snp3 <- ind_chr6_snp2
ind_chr7_snp3 <- ind_chr7_snp2
ind_chr8_snp3 <- ind_chr8_snp2
ind_chr9_snp3 <- ind_chr9_snp2
ind_chr10_snp3 <- ind_chr10_snp2
ind_chr11_snp3 <- ind_chr11_snp2
ind_chr12_snp3 <- ind_chr12_snp2

#adjust start
ind_chr1_snp3$`SNP Start`<- ind_chr1_snp3$`SNP Start`/1000000
ind_chr2_snp3$`SNP Start` <- ind_chr2_snp3$`SNP Start`/1000000
ind_chr3_snp3$`SNP Start` <- ind_chr3_snp3$`SNP Start`/1000000
ind_chr4_snp3$`SNP Start` <- ind_chr4_snp3$`SNP Start`/1000000
ind_chr5_snp3$`SNP Start` <- ind_chr5_snp3$`SNP Start`/1000000
ind_chr6_snp3$`SNP Start` <- ind_chr6_snp3$`SNP Start`/1000000
ind_chr7_snp3$`SNP Start` <- ind_chr7_snp3$`SNP Start`/1000000
ind_chr8_snp3$`SNP Start` <- ind_chr8_snp3$`SNP Start`/1000000
ind_chr9_snp3$`SNP Start` <- ind_chr9_snp3$`SNP Start`/1000000
ind_chr10_snp3$`SNP Start` <- ind_chr10_snp3$`SNP Start`/1000000
ind_chr11_snp3$`SNP Start` <- ind_chr11_snp3$`SNP Start`/1000000
ind_chr12_snp3$`SNP Start` <- ind_chr12_snp3$`SNP Start`/1000000

ind_chr1_snp3<-fill_NA(ind_chr1_snp3)
ind_chr2_snp3<-fill_NA(ind_chr2_snp3)
ind_chr3_snp3<-fill_NA(ind_chr3_snp3)
ind_chr4_snp3<-fill_NA(ind_chr4_snp3)
ind_chr5_snp3<-fill_NA(ind_chr5_snp3)
ind_chr6_snp3<-fill_NA(ind_chr6_snp3)
ind_chr7_snp3<-fill_NA(ind_chr7_snp3)
ind_chr8_snp3<-fill_NA(ind_chr8_snp3)
ind_chr9_snp3<-fill_NA(ind_chr9_snp3)
ind_chr10_snp3<-fill_NA(ind_chr10_snp3)
ind_chr11_snp3<-fill_NA(ind_chr11_snp3)
ind_chr12_snp3<-fill_NA(ind_chr12_snp3)

ind_chr1_snp3 <- ind_chr1_snp3[order(ind_chr1_snp3$`SNP Start`),]
ind_chr1_spl <- smooth.spline(ind_chr1_snp3$rate, spar =0)
ind_chr1_snp3$pos <- gen_pos(ind_chr1_snp3)
plot(ind_chr1_snp3$`SNP Start`, ind_chr1_snp3$pos)
ggplot(ind_chr1_snp3, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
plot(ind_chr1_snp3$`SNP Start`, ind_chr1_snp3$pos/ind_chr1_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Indica Chromosome 1 Recombination Distribution")
ind_chr1_finalpos <- ind_chr1_snp3[order(ind_chr1_snp3$pos),]
is.unsorted(ind_chr1_finalpos$pos)
plot(ind_chr1_snp3$`SNP Start`, ind_chr1_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Indica Chromosome 1 Genetic Map")

ind_chr2_spl <- smooth.spline(ind_chr2_snp3$rate, spar = 0)
ind_chr2_snp3$pos <- gen_pos(ind_chr2_snp3)
plot(ind_chr2_snp3$`SNP Start`, ind_chr2_snp3$pos)
plot(ind_chr2_snp3$`SNP Start`, ind_chr2_snp3$pos/ind_chr2_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Indica Chromosome 2 Recombination Distribution")
ind_chr2_finalpos <- ind_chr2_snp3[order(ind_chr2_snp3$pos),]
is.unsorted(ind_chr2_finalpos$pos)
plot(ind_chr2_snp3$`SNP Start`, ind_chr2_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Indica Chromosome 2 Genetic Map")

ind_chr3_spl <- smooth.spline(ind_chr3_snp3$rate, spar = 0)
ind_chr3_snp3$pos <- gen_pos(ind_chr3_snp3)
plot(ind_chr3_snp3$`SNP Start`, ind_chr3_snp3$pos)
plot(ind_chr3_snp3$`SNP Start`, ind_chr3_snp3$pos/ind_chr3_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Indica Chromosome 3 Recombination Distribution")
ind_chr3_finalpos <- ind_chr3_snp3[order(ind_chr3_snp3$pos),]
is.unsorted(ind_chr3_finalpos$pos)
plot(ind_chr3_snp3$`SNP Start`, ind_chr3_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Indica Chromosome 3 Genetic Map")

ind_chr4_spl <- smooth.spline(ind_chr4_snp3$rate, spar = 0)
ind_chr4_snp3$pos <- gen_pos(ind_chr4_snp3)
plot(ind_chr4_snp3$`SNP Start`, ind_chr4_snp3$pos)
plot(ind_chr4_snp3$`SNP Start`, ind_chr4_snp3$pos/ind_chr4_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Indica Chromosome 4 Recombination Distribution")
ind_chr4_finalpos <- ind_chr4_snp3[order(ind_chr4_snp3$pos),]
is.unsorted(ind_chr4_finalpos$pos)
plot(ind_chr4_snp3$`SNP Start`, ind_chr4_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Indica Chromosome 4 Genetic Map")

ind_chr5_spl <- smooth.spline(ind_chr5_snp3$rate, spar = 0)
ind_chr5_snp3$pos <- gen_pos(ind_chr5_snp3)
plot(ind_chr5_snp3$`SNP Start`, ind_chr5_snp3$pos)
plot(ind_chr5_snp3$`SNP Start`, ind_chr5_snp3$pos/ind_chr5_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Indica Chromosome 5 Recombination Distribution")
ind_chr5_finalpos <- ind_chr5_snp3[order(ind_chr5_snp3$pos),]
is.unsorted(ind_chr5_finalpos$pos)
plot(ind_chr5_snp3$`SNP Start`, ind_chr5_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Indica Chromosome 5 Genetic Map")

ind_chr6_spl <- smooth.spline(ind_chr6_snp3$rate, spar = 0)
ind_chr6_snp3$pos <- gen_pos(ind_chr6_snp3)
plot(ind_chr6_snp3$`SNP Start`, ind_chr6_snp3$pos)
plot(ind_chr6_snp3$`SNP Start`, ind_chr6_snp3$pos/ind_chr6_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Indica Chromosome 6 Recombination Distribution")
ind_chr6_finalpos <- ind_chr6_snp3[order(ind_chr6_snp3$pos),]
is.unsorted(ind_chr6_finalpos$pos)
plot(ind_chr6_snp3$`SNP Start`, ind_chr6_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Indica Chromosome 6 Genetic Map")

ind_chr7_spl <- smooth.spline(ind_chr7_snp3$rate, spar =0)
ind_chr7_snp3$pos <- gen_pos(ind_chr7_snp3)
plot(ind_chr7_snp3$`SNP Start`, ind_chr7_snp3$pos)
plot(ind_chr7_snp3$`SNP Start`, ind_chr7_snp3$pos/ind_chr7_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Indica Chromosome 7 Recombination Distribution")
ind_chr7_finalpos <- ind_chr7_snp3[order(ind_chr7_snp3$pos),]
is.unsorted(ind_chr7_finalpos$pos)
plot(ind_chr7_snp3$`SNP Start`, ind_chr7_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Indica Chromosome 7 Genetic Map")

ind_chr8_spl <- smooth.spline(ind_chr8_snp3$rate, spar = 0)
ind_chr8_snp3$pos <- gen_pos(ind_chr8_snp3)
plot(ind_chr8_snp3$`SNP Start`, ind_chr8_snp3$pos)
plot(ind_chr8_snp3$`SNP Start`, ind_chr8_snp3$pos/ind_chr8_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Indica Chromosome 8 Recombination Distribution")
ind_chr8_finalpos <- ind_chr8_snp3[order(ind_chr8_snp3$pos),]
is.unsorted(ind_chr8_finalpos$pos)
plot(ind_chr8_snp3$`SNP Start`, ind_chr8_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Indica Chromosome 8 Genetic Map")

ind_chr9_spl <- smooth.spline(ind_chr9_snp3$rate, spar =0)
ind_chr9_snp3$pos <- gen_pos(ind_chr9_snp3)
plot(ind_chr9_snp3$`SNP Start`, ind_chr9_snp3$pos)
plot(ind_chr9_snp3$`SNP Start`, ind_chr9_snp3$pos/ind_chr9_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Indica Chromosome 9 Recombination Distribution")
ind_chr9_finalpos <- ind_chr9_snp3[order(ind_chr9_snp3$pos),]
is.unsorted(ind_chr9_finalpos$pos)
plot(ind_chr9_snp3$`SNP Start`, ind_chr9_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Indica Chromosome 9 Genetic Map")

ind_chr10_spl <- smooth.spline(ind_chr10_snp3$rate, spar = 0)
ind_chr10_snp3$pos <- gen_pos(ind_chr10_snp3)
plot(ind_chr10_snp3$`SNP Start`, ind_chr10_snp3$pos)
plot(ind_chr10_snp3$`SNP Start`, ind_chr10_snp3$pos/ind_chr10_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Indica Chromosome 10 Recombination Distribution")
ind_chr10_finalpos <- ind_chr10_snp3[order(ind_chr10_snp3$pos),]
is.unsorted(ind_chr10_finalpos$pos)
plot(ind_chr10_snp3$`SNP Start`, ind_chr10_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Indica Chromosome 10 Genetic Map")

ind_chr11_spl <- smooth.spline(ind_chr11_snp3$rate, spar =0)
ind_chr11_snp3$pos <- gen_pos(ind_chr11_snp3)
plot(ind_chr11_snp3$`SNP Start`, ind_chr11_snp3$pos)
plot(ind_chr11_snp3$`SNP Start`, ind_chr11_snp3$pos/ind_chr11_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Indica Chromosome 11 Recombination Distribution")
ind_chr11_finalpos <- ind_chr11_snp3[order(ind_chr11_snp3$pos),]
is.unsorted(ind_chr11_finalpos$pos)
plot(ind_chr11_snp3$`SNP Start`, ind_chr11_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Indica Chromosome 11 Genetic Map")

ind_chr12_spl <- smooth.spline(ind_chr12_snp3$rate, spar = 0)
ind_chr12_snp3$pos <- gen_pos(ind_chr12_snp3)
plot(ind_chr12_snp3$`SNP Start`, ind_chr12_snp3$pos)
plot(ind_chr12_snp3$`SNP Start`, ind_chr12_snp3$pos/ind_chr12_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Indica Chromosome 12 Genetic Map")
ind_chr12_finalpos <- ind_chr12_snp3[order(ind_chr12_snp3$pos),]
is.unsorted(ind_chr12_finalpos$pos)
plot(ind_chr12_snp3$`SNP Start`, ind_chr12_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Indica Chromosome 12 Genetic Map")

#Japonica final genetic map
chr1 <- japonica_chr1_finalpos$pos/100
chr2 <- japonica_chr2_finalpos$pos/100
chr3 <- japonica_chr3_finalpos$pos/100
chr4 <- japonica_chr4_finalpos$pos/100
chr5 <- japonica_chr5_finalpos$pos/100
chr6 <- japonica_chr6_finalpos$pos/100
chr7 <- japonica_chr7_finalpos$pos/100
chr8 <- japonica_chr8_finalpos$pos/100
chr9 <- japonica_chr9_finalpos$pos/100
chr10<- japonica_chr10_finalpos$pos/100
chr11 <- japonica_chr11_finalpos$pos/100
chr12<- japonica_chr12_finalpos$pos/100

segSites<-readRDS("japonica_num_SNP.RData")
japonica_map = vector("list",10)
japonica_map[[1]] = chr1
japonica_map[[2]] = chr2
japonica_map[[3]] = chr3
japonica_map[[4]] = chr4
japonica_map[[5]] = chr5
japonica_map[[6]] = chr6
japonica_map[[7]] = chr7
japonica_map[[8]] = chr8
japonica_map[[9]] = chr9
japonica_map[[10]] = chr10
japonica_map[[11]] = chr11
japonica_map[[12]] = chr12
for(i in 1:12){
  names(japonica_map[[i]]) = paste(i, 1:segSites[i], sep="_")
}

saveRDS(japonica_map, file="japonica_final_map.RData")

#actual positions:http://rice.uga.edu/annotation_pseudo_centromeres.shtml
# 1- 16.7
# 2- 13.6 
# 3- 19.4
# 4- 9.7
# 5- 12.4
# 6- 15.3
# 7- 12.1
# 8- 12.9
# 9- 2.8
# 10- 8.2
# 11- 12
# 12- 11.9
find_centromere<-function(centromere,finalpos){
  row.names(finalpos) <- NULL
  index<- finalpos[which.min(abs(centromere-finalpos$`SNP Start`)),]
  row<-as.integer(rownames(index))
  print(finalpos$pos[row])
}
c1 <-find_centromere(16.7,japonica_chr1_finalpos)
c2 <-find_centromere(13.6,japonica_chr2_finalpos)
c3 <-find_centromere(19.4,japonica_chr3_finalpos)
c4 <-find_centromere(9.7,japonica_chr4_finalpos)
c5 <-find_centromere(12.4,japonica_chr5_finalpos)
c6 <-find_centromere(15.3,japonica_chr6_finalpos)
c7 <-find_centromere(12.1,japonica_chr7_finalpos)
c8 <-find_centromere(12.9,japonica_chr8_finalpos)
c9 <-find_centromere(2.8,japonica_chr9_finalpos)
c10 <-find_centromere(8.2,japonica_chr10_finalpos)
c11 <-find_centromere(12,japonica_chr11_finalpos)
c12 <-find_centromere(11.9,japonica_chr12_finalpos)

japonica_centromere <- c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12)
japonica_centromere <- japonica_centromere/100

c1 <-find_centromere(16.7,ind_chr1_finalpos)
c2 <-find_centromere(13.6,ind_chr2_finalpos)
c3 <-find_centromere(19.4,ind_chr3_finalpos)
c4 <-find_centromere(9.7,ind_chr4_finalpos)
c5 <-find_centromere(12.4,ind_chr5_finalpos)
c6 <-find_centromere(15.3,ind_chr6_finalpos)
c7 <-find_centromere(12.1,ind_chr7_finalpos)
c8 <-find_centromere(12.9,ind_chr8_finalpos)
c9 <-find_centromere(2.8,ind_chr9_finalpos)
c10 <-find_centromere(8.2,ind_chr10_finalpos)
c11 <-find_centromere(12,ind_chr11_finalpos)
c12 <-find_centromere(11.9,ind_chr12_finalpos)

ind_centromere <- c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12)
ind_centromere <- ind_centromere/100

saveRDS(ind_centromere, file="indica_centromeres.RData")
saveRDS(japonica_centromere, file="japonica_centromeres.RData")

saveRDS(japonica_chr1_finalpos, file="japonica_chr1_finalpos.RData")
saveRDS(japonica_chr2_finalpos, file="japonica_chr2_finalpos.RData")
saveRDS(japonica_chr3_finalpos, file="japonica_chr3_finalpos.RData")
saveRDS(japonica_chr4_finalpos, file="japonica_chr4_finalpos.RData")
saveRDS(japonica_chr5_finalpos, file="japonica_chr5_finalpos.RData")
saveRDS(japonica_chr6_finalpos, file="japonica_chr6_finalpos.RData")
saveRDS(japonica_chr7_finalpos, file="japonica_chr7_finalpos.RData")
saveRDS(japonica_chr8_finalpos, file="japonica_chr8_finalpos.RData")
saveRDS(japonica_chr9_finalpos, file="japonica_chr9_finalpos.RData")
saveRDS(japonica_chr10_finalpos, file="japonica_chr10_finalpos.RData")
saveRDS(japonica_chr11_finalpos, file="japonica_chr11_finalpos.RData")
saveRDS(japonica_chr12_finalpos, file="japonica_chr12_finalpos.RData")

japonica_nSNP<-c(nrow(japonica_chr1_finalpos),nrow(japonica_chr2_finalpos),nrow(japonica_chr3_finalpos),nrow(japonica_chr4_finalpos),nrow(japonica_chr5_finalpos),nrow(japonica_chr6_finalpos),nrow(japonica_chr7_finalpos),nrow(japonica_chr8_finalpos),nrow(japonica_chr9_finalpos),nrow(japonica_chr10_finalpos),nrow(japonica_chr11_finalpos),nrow(japonica_chr12_finalpos))
saveRDS(japonica_nSNP, file="japonica_num_SNP.RData")
