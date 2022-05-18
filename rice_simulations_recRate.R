library(AlphaSimR)
library(Rcpp)
library(ggplot2)
library(dplyr)

#setwd("C:/Users/16192/Documents/PNAS_Simulations")
set.seed(420)

jap_snps <- read.table("japonica_SNPs.bed", header =FALSE)
ind_snps <- read.table("indica_SNPs.bed", header =FALSE)
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

#recombination rates
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

isTRUE(jap_chr1_CO[1,2] == 0)
isTRUE(jap_chr2_CO[1,2] == 0)
isTRUE(jap_chr3_CO[1,2] == 0)
isTRUE(jap_chr4_CO[1,2] == 0)
isTRUE(jap_chr5_CO[1,2] == 0)
isTRUE(jap_chr6_CO[1,2] == 0)
isTRUE(jap_chr7_CO[1,2] == 0)
isTRUE(jap_chr8_CO[1,2] == 0)
isTRUE(jap_chr9_CO[1,2] == 0)
isTRUE(jap_chr10_CO[1,2] == 0)
isTRUE(jap_chr11_CO[1,2] == 0)
isTRUE(jap_chr12_CO[1,2] == 0)

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
jap_chr1_snp2 <- snp_rate(jap_chr1_CO_2, jap_chr1_snp)
jap_chr2_snp2 <- snp_rate(jap_chr2_CO_2, jap_chr2_snp)
jap_chr3_snp2 <- snp_rate(jap_chr3_CO_2, jap_chr3_snp)
jap_chr4_snp2 <- snp_rate(jap_chr4_CO_2, jap_chr4_snp)
jap_chr5_snp2 <- snp_rate(jap_chr5_CO_2, jap_chr5_snp)
jap_chr6_snp2 <- snp_rate(jap_chr6_CO_2, jap_chr6_snp)
jap_chr7_snp2 <- snp_rate(jap_chr7_CO_2, jap_chr7_snp)
jap_chr8_snp2 <- snp_rate(jap_chr8_CO_2, jap_chr8_snp)
jap_chr9_snp2 <- snp_rate(jap_chr9_CO_2, jap_chr9_snp)
jap_chr10_snp2 <- snp_rate(jap_chr10_CO_2, jap_chr10_snp)
jap_chr11_snp2 <- snp_rate(jap_chr11_CO_2, jap_chr11_snp)
jap_chr12_snp2 <- snp_rate(jap_chr12_CO_2, jap_chr12_snp)

#make mutable copies of data
jap_chr1_snp3 <- jap_chr1_snp2
jap_chr2_snp3 <- jap_chr2_snp2
jap_chr3_snp3 <- jap_chr3_snp2
jap_chr4_snp3 <- jap_chr4_snp2
jap_chr5_snp3 <- jap_chr5_snp2
jap_chr6_snp3 <- jap_chr6_snp2
jap_chr7_snp3 <- jap_chr7_snp2
jap_chr8_snp3 <- jap_chr8_snp2
jap_chr9_snp3 <- jap_chr9_snp2
jap_chr10_snp3 <- jap_chr10_snp2
jap_chr11_snp3 <- jap_chr11_snp2
jap_chr12_snp3 <- jap_chr12_snp2

#adjust start
jap_chr1_snp3$`SNP Start`<- jap_chr1_snp3$`SNP Start`/1000000
jap_chr2_snp3$`SNP Start` <- jap_chr2_snp3$`SNP Start`/1000000
jap_chr3_snp3$`SNP Start` <- jap_chr3_snp3$`SNP Start`/1000000
jap_chr4_snp3$`SNP Start` <- jap_chr4_snp3$`SNP Start`/1000000
jap_chr5_snp3$`SNP Start` <- jap_chr5_snp3$`SNP Start`/1000000
jap_chr6_snp3$`SNP Start` <- jap_chr6_snp3$`SNP Start`/1000000
jap_chr7_snp3$`SNP Start` <- jap_chr7_snp3$`SNP Start`/1000000
jap_chr8_snp3$`SNP Start` <- jap_chr8_snp3$`SNP Start`/1000000
jap_chr9_snp3$`SNP Start` <- jap_chr9_snp3$`SNP Start`/1000000
jap_chr10_snp3$`SNP Start` <- jap_chr10_snp3$`SNP Start`/1000000
jap_chr11_snp3$`SNP Start` <- jap_chr11_snp3$`SNP Start`/1000000
jap_chr12_snp3$`SNP Start` <- jap_chr12_snp3$`SNP Start`/1000000

fill_NA<-function(SNP){
  for(i in 1:nrow(SNP)){
    if(is.na(SNP$rate[i])){
      SNP$rate[i]<-SNP$rate[i-1]
    }
  }
  print(SNP)
}
jap_chr1_snp3<-fill_NA(jap_chr1_snp3)
jap_chr2_snp3<-fill_NA(jap_chr2_snp3)
jap_chr3_snp3<-fill_NA(jap_chr3_snp3)
jap_chr4_snp3<-fill_NA(jap_chr4_snp3)
jap_chr5_snp3<-fill_NA(jap_chr5_snp3)
jap_chr6_snp3<-fill_NA(jap_chr6_snp3)
jap_chr7_snp3<-fill_NA(jap_chr7_snp3)
jap_chr8_snp3<-fill_NA(jap_chr8_snp3)
jap_chr9_snp3<-fill_NA(jap_chr9_snp3)
jap_chr10_snp3<-fill_NA(jap_chr10_snp3)
jap_chr11_snp3<-fill_NA(jap_chr11_snp3)
jap_chr12_snp3<-fill_NA(jap_chr12_snp3)

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

jap_chr1_spl <- smooth.spline(jap_chr1_snp3$rate, spar=0.1)
jap_chr1_snp3$pos <- gen_pos(jap_chr1_snp3)
plot(jap_chr1_snp3$`SNP Start`, jap_chr1_snp3$pos, type = "l")
ggplot(jap_chr1_snp3, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
plot(jap_chr1_snp3$`SNP Start`, jap_chr1_snp3$pos/jap_chr1_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica jap Chromosome 1 Recombination Distribution")
jap_chr1_finalpos <- jap_chr1_snp3[order(jap_chr1_snp3$pos),]
is.unsorted(jap_chr1_finalpos$pos)
plot(jap_chr1_snp3$`SNP Start`, jap_chr1_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica jap Chromosome 1 Genetic Map")
plot(jap_chr1_finalpos$`SNP Start`, jap_chr1_finalpos$pos)
saveRDS(jap_chr1_finalpos,"chr1_snp2jap.RData")

jap_chr2_spl <- smooth.spline(jap_chr2_snp3$rate, spar = 0.1)
jap_chr2_snp3$pos <- gen_pos(jap_chr2_snp3)
plot(jap_chr2_snp3$`SNP Start`, jap_chr2_snp3$pos)
plot(jap_chr2_snp3$`SNP Start`, jap_chr2_snp3$pos/jap_chr2_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica jap Chromosome 2 Recombination Distribution")
jap_chr2_finalpos <- jap_chr2_snp3[order(jap_chr2_snp3$pos),]
is.unsorted(jap_chr2_finalpos$pos)
plot(jap_chr2_snp3$`SNP Start`, jap_chr2_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica jap Chromosome 2 Genetic Map")
saveRDS(jap_chr2_finalpos,"chr2_snp2jap.RData")

jap_chr3_spl <- smooth.spline(jap_chr3_snp3$rate, spar =.1)
jap_chr3_snp3$pos <- gen_pos(jap_chr3_snp3)
plot(jap_chr3_snp3$`SNP Start`, jap_chr3_snp3$pos, type = "l")
plot(jap_chr3_snp3$`SNP Start`, jap_chr3_snp3$pos/jap_chr3_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica jap Chromosome 3 Recombination Distribution")
jap_chr3_finalpos <- jap_chr3_snp3[order(jap_chr3_snp3$pos),]
is.unsorted(jap_chr3_finalpos$pos)
#jap_chr3_finalpos$pos <- jap_chr3_finalpos$pos + abs(min(jap_chr3_finalpos$pos))
plot(jap_chr3_snp3$`SNP Start`, jap_chr3_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica jap Chromosome 3 Genetic Map")
saveRDS(jap_chr3_finalpos,"chr3_snp2jap.RData")

jap_chr4_spl <- smooth.spline(jap_chr4_snp3$rate, spar =.1)
jap_chr4_snp3$pos <- gen_pos(jap_chr4_snp3)
plot(jap_chr4_snp3$`SNP Start`, jap_chr4_snp3$pos)
plot(jap_chr4_snp3$`SNP Start`, jap_chr4_snp3$pos/jap_chr4_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica jap Chromosome 4 Recombination Distribution")
jap_chr4_finalpos <- jap_chr4_snp3[order(jap_chr4_snp3$pos),]
is.unsorted(jap_chr4_finalpos$pos)
plot(jap_chr4_snp3$`SNP Start`, jap_chr4_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica jap Chromosome 4 Genetic Map")
saveRDS(jap_chr4_finalpos,"chr4_snp2jap.RData")

jap_chr5_spl <- smooth.spline(jap_chr5_snp3$rate, spar =.1)
jap_chr5_snp3$pos <- gen_pos(jap_chr5_snp3)
plot(jap_chr5_snp3$`SNP Start`, jap_chr5_snp3$pos)
plot(jap_chr5_snp3$`SNP Start`, jap_chr5_snp3$pos/jap_chr5_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica jap Chromosome 5 Recombination Distribution")
jap_chr5_finalpos <- jap_chr5_snp3[order(jap_chr5_snp3$pos),]
is.unsorted(jap_chr5_finalpos$pos)
plot(jap_chr5_snp3$`SNP Start`, jap_chr5_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica jap Chromosome 5 Genetic Map")
saveRDS(jap_chr5_finalpos,"chr5_snp2jap.RData")

jap_chr6_spl <- smooth.spline(jap_chr6_snp3$rate, spar = .1)
jap_chr6_snp3$pos <- gen_pos(jap_chr6_snp3)
plot(jap_chr6_snp3$`SNP Start`, jap_chr6_snp3$pos)
plot(jap_chr6_snp3$`SNP Start`, jap_chr6_snp3$pos/jap_chr6_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica jap Chromosome 6 Recombination Distribution")
jap_chr6_finalpos <- jap_chr6_snp3[order(jap_chr6_snp3$pos),]
is.unsorted(jap_chr6_finalpos$pos)
plot(jap_chr6_snp3$`SNP Start`, jap_chr6_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica jap Chromosome 6 Genetic Map")
saveRDS(jap_chr6_finalpos,"chr6_snp2jap.RData")

jap_chr7_spl <- smooth.spline(jap_chr7_snp3$rate, spar = .1)
jap_chr7_snp3$pos <- gen_pos(jap_chr7_snp3)
plot(jap_chr7_snp3$`SNP Start`, jap_chr7_snp3$pos)
plot(jap_chr7_snp3$`SNP Start`, jap_chr7_snp3$pos/jap_chr7_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica jap Chromosome 7 Recombination Distribution")
jap_chr7_finalpos <- jap_chr7_snp3[order(jap_chr7_snp3$pos),]
is.unsorted(jap_chr7_finalpos$pos)
plot(jap_chr7_snp3$`SNP Start`, jap_chr7_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica jap Chromosome 7 Genetic Map")
saveRDS(jap_chr7_finalpos,"chr7_snp2jap.RData")

jap_chr8_spl <- smooth.spline(jap_chr8_snp3$rate, spar = .1)
jap_chr8_snp3$pos <- gen_pos(jap_chr8_snp3)
plot(jap_chr8_snp3$`SNP Start`, jap_chr8_snp3$pos)
plot(jap_chr8_snp3$`SNP Start`, jap_chr8_snp3$pos/jap_chr8_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica jap Chromosome 8 Recombination Distribution")
jap_chr8_finalpos <- jap_chr8_snp3[order(jap_chr8_snp3$pos),]
is.unsorted(jap_chr8_finalpos$pos)
plot(jap_chr8_snp3$`SNP Start`, jap_chr8_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica jap Chromosome 8 Genetic Map")
saveRDS(jap_chr8_finalpos,"chr8_snp2jap.RData")

jap_chr9_spl <- smooth.spline(jap_chr9_snp3$rate, spar = .1)
jap_chr9_snp3$pos <- gen_pos(jap_chr9_snp3)
plot(jap_chr9_snp3$`SNP Start`, jap_chr9_snp3$pos)
plot(jap_chr9_snp3$`SNP Start`, jap_chr9_snp3$pos/jap_chr9_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica jap Chromosome 9 Recombination Distribution")
jap_chr9_finalpos <- jap_chr9_snp3[order(jap_chr9_snp3$pos),]
is.unsorted(jap_chr9_finalpos$pos)
plot(jap_chr9_snp3$`SNP Start`, jap_chr9_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica jap Chromosome 9 Genetic Map")
saveRDS(jap_chr9_finalpos,"chr9_snp2jap.RData")

jap_chr10_spl <- smooth.spline(jap_chr10_snp3$rate, spar =.1)
jap_chr10_snp3$pos <- gen_pos(jap_chr10_snp3)
plot(jap_chr10_snp3$`SNP Start`, jap_chr10_snp3$pos)
plot(jap_chr10_snp3$`SNP Start`, jap_chr10_snp3$pos/jap_chr10_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica jap Chromosome 10 Recombination Distribution")
jap_chr10_finalpos <- jap_chr10_snp3[order(jap_chr10_snp3$pos),]
is.unsorted(jap_chr10_finalpos$pos)
plot(jap_chr10_snp3$`SNP Start`, jap_chr10_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica jap Chromosome 10 Genetic Map")
saveRDS(jap_chr10_finalpos,"chr10_snp2jap.RData")

jap_chr11_spl <- smooth.spline(jap_chr11_snp3$rate, spar = .1)
jap_chr11_snp3$pos <- gen_pos(jap_chr11_snp3)
plot(jap_chr11_snp3$`SNP Start`, jap_chr11_snp3$pos)
plot(jap_chr11_snp3$`SNP Start`, jap_chr11_snp3$pos/jap_chr11_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica jap Chromosome 11 Recombination Distribution")
jap_chr11_finalpos <- jap_chr11_snp3[order(jap_chr11_snp3$pos),]
is.unsorted(jap_chr11_finalpos$pos)
plot(jap_chr11_snp3$`SNP Start`, jap_chr11_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica jap Chromosome 11 Genetic Map")
saveRDS(jap_chr11_finalpos,"chr11_snp2jap.RData")

jap_chr12_spl <- smooth.spline(jap_chr12_snp3$rate, spar = .1)
jap_chr12_snp3$pos <- gen_pos(jap_chr12_snp3)
plot(jap_chr12_snp3$`SNP Start`, jap_chr12_snp3$pos)
plot(jap_chr12_snp3$`SNP Start`, jap_chr12_snp3$pos/jap_chr12_snp3$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica jap Chromosome 12 Recombination Distribution")
jap_chr12_finalpos <- jap_chr12_snp3[order(jap_chr12_snp3$pos),]
is.unsorted(jap_chr12_finalpos$pos)
plot(jap_chr12_snp3$`SNP Start`, jap_chr12_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica jap Chromosome 12 Genetic Map")
saveRDS(jap_chr12_finalpos,"chr12_snp2jap.RData")

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

jap_chr6 <- jap_chr6_finalpos$pos/100
jap_chr6len <- length(jap_chr6)
dim(jap_chr6) <- c(jap_chr6len,1)
jap_chr6 <- list(jap_chr6)

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

japonica_final_map <- list(jap_chr1[[1]], jap_chr2[[1]], 
                  jap_chr3[[1]], jap_chr4[[1]], jap_chr5[[1]], 
                  jap_chr6[[1]], jap_chr7[[1]], jap_chr8[[1]], 
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

indica_final_map <- list(ind_chr1[[1]], ind_chr2[[1]], 
                  ind_chr3[[1]], ind_chr4[[1]], ind_chr5[[1]], 
                  ind_chr6[[1]], ind_chr7[[1]], ind_chr8[[1]], 
                  ind_chr9[[1]], ind_chr10[[1]], ind_chr11[[1]], ind_chr12[[1]])
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
c1 <-find_centromere(16.7,jap_chr1_finalpos)
c2 <-find_centromere(13.6,jap_chr2_finalpos)
c3 <-find_centromere(19.4,jap_chr3_finalpos)
c4 <-find_centromere(9.7,jap_chr4_finalpos)
c5 <-find_centromere(12.4,jap_chr5_finalpos)
c6 <-find_centromere(15.3,jap_chr6_finalpos)
c7 <-find_centromere(12.1,jap_chr7_finalpos)
c8 <-find_centromere(12.9,jap_chr8_finalpos)
c9 <-find_centromere(2.8,jap_chr9_finalpos)
c10 <-find_centromere(8.2,jap_chr10_finalpos)
c11 <-find_centromere(12,jap_chr11_finalpos)
c12 <-find_centromere(11.9,jap_chr12_finalpos)

jap_centromere <- c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12)
jap_centromere <- jap_centromere/100

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

saveRDS(japonica_final_map, file="japonica_final_map.RData")
saveRDS(indica_final_map, file="indica_final_map.RData")

saveRDS(ind_centromere, file="indica_centromeres.RData")
saveRDS(jap_centromere, file="japonica_centromeres.RData")

saveRDS(jap_chr1_finalpos, file="jap_chr1_finalpos.RData")
saveRDS(jap_chr2_finalpos, file="jap_chr2_finalpos.RData")
saveRDS(jap_chr3_finalpos, file="jap_chr3_finalpos.RData")
saveRDS(jap_chr4_finalpos, file="jap_chr4_finalpos.RData")
saveRDS(jap_chr5_finalpos, file="jap_chr5_finalpos.RData")
saveRDS(jap_chr6_finalpos, file="jap_chr6_finalpos.RData")
saveRDS(jap_chr7_finalpos, file="jap_chr7_finalpos.RData")
saveRDS(jap_chr8_finalpos, file="jap_chr8_finalpos.RData")
saveRDS(jap_chr9_finalpos, file="jap_chr9_finalpos.RData")
saveRDS(jap_chr10_finalpos, file="jap_chr10_finalpos.RData")
saveRDS(jap_chr11_finalpos, file="jap_chr11_finalpos.RData")
saveRDS(jap_chr12_finalpos, file="jap_chr12_finalpos.RData")

jap_nSNP<-c(nrow(jap_chr1_finalpos),nrow(jap_chr2_finalpos),nrow(jap_chr3_finalpos),nrow(jap_chr4_finalpos),nrow(jap_chr5_finalpos),nrow(jap_chr6_finalpos),nrow(jap_chr7_finalpos),nrow(jap_chr8_finalpos),nrow(jap_chr9_finalpos),nrow(jap_chr10_finalpos),nrow(jap_chr11_finalpos),nrow(jap_chr12_finalpos))
saveRDS(jap_nSNP, file="japonica_num_SNP.RData")
