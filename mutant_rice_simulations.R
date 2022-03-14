library(AlphaSimR)
library(Rcpp)
library(ggplot2)
library(dplyr)

set.seed(420)

jap_snps <- read.table("japonica_SNPs.bed", header =FALSE)
colnames(jap_snps) <- c("Chr#", "SNP Start", "SNP End")
colnames(ind_snps) <- c("Chr#", "SNP Start", "SNP End")
#sample SNPs?
jap_snps <- sample_n(jap_snps, 4000)
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


#jap mutant recomb rates

jap_CO <- read.csv("jap_WT_rate.csv", header = TRUE)
colnames(jap_CO) <- c("Chr", "CO Start", "CO End", "rate")
jap_CO <- jap_CO[order(jap_CO$Chr,jap_CO$`CO Start`),]

jap_chr1_CO <- jap_CO[ which(jap_CO$Chr == "1"),]
jap_chr1_CO$midpoint <- (jap_chr1_CO$`CO Start`+ jap_chr1_CO$`CO End`)/2
jap_chr1_CO <- jap_chr1_CO[order(jap_chr1_CO$`CO Start`),]

jap_chr2_CO <- jap_CO[ which(jap_CO$Chr == "2"),]
jap_chr2_CO$midpoint <- (jap_chr2_CO$`CO Start`+ jap_chr2_CO$`CO End`)/2
jap_chr2_CO <- jap_chr2_CO[order(jap_chr2_CO$`CO Start`),]

jap_chr3_CO <- jap_CO[ which(jap_CO$Chr == "3"),]
jap_chr3_CO$midpoint <- (jap_chr3_CO$`CO Start`+ jap_chr3_CO$`CO End`)/2
jap_chr3_CO <- jap_chr3_CO[order(jap_chr3_CO$`CO Start`),]

jap_chr4_CO <- jap_CO[ which(jap_CO$Chr == "4"),]
jap_chr4_CO$midpoint <- (jap_chr4_CO$`CO Start`+ jap_chr4_CO$`CO End`)/2
jap_chr4_CO <- jap_chr4_CO[order(jap_chr4_CO$`CO Start`),]

jap_chr5_CO <- jap_CO[ which(jap_CO$Chr == "5"),]
jap_chr5_CO$midpoint <- (jap_chr5_CO$`CO Start`+ jap_chr5_CO$`CO End`)/2
jap_chr5_CO <- jap_chr5_CO[order(jap_chr5_CO$`CO Start`),]

jap_chr6_CO <- jap_CO[ which(jap_CO$Chr == "6"),]
jap_chr6_CO$midpoint <- (jap_chr6_CO$`CO Start`+ jap_chr6_CO$`CO End`)/2
jap_chr6_CO <- jap_chr6_CO[order(jap_chr6_CO$`CO Start`),]

jap_chr7_CO <- jap_CO[ which(jap_CO$Chr == "7"),]
jap_chr7_CO$midpoint <- (jap_chr7_CO$`CO Start`+ jap_chr7_CO$`CO End`)/2
jap_chr7_CO <- jap_chr7_CO[order(jap_chr7_CO$`CO Start`),]

jap_chr8_CO <- jap_CO[ which(jap_CO$Chr == "8"),]
jap_chr8_CO$midpoint <- (jap_chr8_CO$`CO Start`+ jap_chr8_CO$`CO End`)/2
jap_chr8_CO <- jap_chr8_CO[order(jap_chr8_CO$`CO Start`),]

jap_chr9_CO <- jap_CO[ which(jap_CO$Chr == "9"),]
jap_chr9_CO$midpoint <- (jap_chr9_CO$`CO Start`+ jap_chr9_CO$`CO End`)/2
jap_chr9_CO <- jap_chr9_CO[order(jap_chr9_CO$`CO Start`),]

jap_chr10_CO <- jap_CO[ which(jap_CO$Chr == "10"),]
jap_chr10_CO$midpoint <- (jap_chr10_CO$`CO Start`+ jap_chr10_CO$`CO End`)/2
jap_chr10_CO <- jap_chr10_CO[order(jap_chr10_CO$`CO Start`),]

jap_chr11_CO <- jap_CO[ which(jap_CO$Chr == "11"),]
jap_chr11_CO$midpoint <- (jap_chr11_CO$`CO Start`+ jap_chr11_CO$`CO End`)/2
jap_chr11_CO <- jap_chr11_CO[order(jap_chr11_CO$`CO Start`),]

jap_chr12_CO <- jap_CO[ which(jap_CO$Chr == "12"),]
jap_chr12_CO$midpoint <- (jap_chr12_CO$`CO Start`+ jap_chr12_CO$`CO End`)/2
jap_chr12_CO <- jap_chr12_CO[order(jap_chr12_CO$`CO Start`),]


#calculating recombination rate per bin of CO data
library(dlookr)
library(tidyverse)
library(OneR)

#making intervals start at 0
jap_chr1_CO$`CO Start` <- jap_chr1_CO$`CO Start` - min(jap_chr1_CO$`CO Start`)
jap_chr1_CO$`CO End` <- jap_chr1_CO$`CO End` - min(jap_chr1_CO$`CO Start`)

jap_chr2_CO$`CO Start` <- jap_chr2_CO$`CO Start` - min(jap_chr2_CO$`CO Start`)
jap_chr2_CO$`CO End` <- jap_chr2_CO$`CO End` - min(jap_chr2_CO$`CO Start`)

jap_chr3_CO$`CO Start` <- jap_chr3_CO$`CO Start` - min(jap_chr3_CO$`CO Start`)
jap_chr3_CO$`CO End` <- jap_chr3_CO$`CO End` - min(jap_chr3_CO$`CO Start`)

jap_chr4_CO$`CO Start` <- jap_chr4_CO$`CO Start` - min(jap_chr4_CO$`CO Start`)
jap_chr4_CO$`CO End` <- jap_chr4_CO$`CO End` - min(jap_chr4_CO$`CO Start`)

jap_chr5_CO$`CO Start` <- jap_chr5_CO$`CO Start` - min(jap_chr5_CO$`CO Start`)
jap_chr5_CO$`CO End` <- jap_chr5_CO$`CO End` - min(jap_chr5_CO$`CO Start`)

jap_chr6_CO$`CO Start` <- jap_chr6_CO$`CO Start` - min(jap_chr6_CO$`CO Start`)
jap_chr6_CO$`CO End` <- jap_chr6_CO$`CO End` - min(jap_chr6_CO$`CO Start`)

jap_chr7_CO$`CO Start` <- jap_chr7_CO$`CO Start` - min(jap_chr7_CO$`CO Start`)
jap_chr7_CO$`CO End` <- jap_chr7_CO$`CO End` - min(jap_chr7_CO$`CO Start`)

jap_chr8_CO$`CO Start` <- jap_chr8_CO$`CO Start` - min(jap_chr8_CO$`CO Start`)
jap_chr8_CO$`CO End` <- jap_chr8_CO$`CO End` - min(jap_chr8_CO$`CO Start`)

jap_chr9_CO$`CO Start` <- jap_chr9_CO$`CO Start` - min(jap_chr9_CO$`CO Start`)
jap_chr9_CO$`CO End` <- jap_chr9_CO$`CO End` - min(jap_chr9_CO$`CO Start`)

jap_chr10_CO$`CO Start` <- jap_chr10_CO$`CO Start` - min(jap_chr10_CO$`CO Start`)
jap_chr10_CO$`CO End` <- jap_chr10_CO$`CO End` - min(jap_chr10_CO$`CO Start`)

jap_chr11_CO$`CO Start` <- jap_chr11_CO$`CO Start` - min(jap_chr11_CO$`CO Start`)
jap_chr11_CO$`CO End` <- jap_chr11_CO$`CO End` - min(jap_chr11_CO$`CO Start`)

jap_chr12_CO$`CO Start` <- jap_chr12_CO$`CO Start` - min(jap_chr12_CO$`CO Start`)
jap_chr12_CO$`CO End` <- jap_chr12_CO$`CO End` - min(jap_chr12_CO$`CO Start`)

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


##assigning frequency to SNPs based on recombination frequency in each bin
snp_rate <- function(chr_rate, chr_snp){
  for(i in 1:nrow(chr_snp)){
    for(k in 1:nrow(chr_rate)){
      if(isTRUE((chr_snp$`SNP Start`[i] >= chr_rate$`CO Start`[k]) && (chr_snp$`SNP Start`[i] <= chr_rate$`CO End`[k]))){
        chr_snp$rate[i] <- chr_rate$rate[k]
      }
    }
  }
  print(chr_snp)
}

#converted SNP start to Mb
jap_chr1_snp$`SNP Start`<- jap_chr1_snp$`SNP Start`/1000000
jap_chr2_snp$`SNP Start` <- jap_chr2_snp$`SNP Start`/1000000
jap_chr3_snp$`SNP Start` <- jap_chr3_snp$`SNP Start`/1000000
jap_chr4_snp$`SNP Start` <- jap_chr4_snp$`SNP Start`/1000000
jap_chr5_snp$`SNP Start` <- jap_chr5_snp$`SNP Start`/1000000
jap_chr6_snp$`SNP Start` <- jap_chr6_snp$`SNP Start`/1000000
jap_chr7_snp$`SNP Start` <- jap_chr7_snp$`SNP Start`/1000000
jap_chr8_snp$`SNP Start` <- jap_chr8_snp$`SNP Start`/1000000
jap_chr9_snp$`SNP Start` <- jap_chr9_snp$`SNP Start`/1000000
jap_chr10_snp$`SNP Start` <- jap_chr10_snp$`SNP Start`/1000000
jap_chr11_snp$`SNP Start` <- jap_chr11_snp$`SNP Start`/1000000
jap_chr12_snp$`SNP Start` <- jap_chr12_snp$`SNP Start`/1000000

#using function,  get cM/Mb for final genetic position - assign rates
jap_chr1_snp2 <- snp_rate(jap_chr1_CO, jap_chr1_snp)
jap_chr2_snp2 <- snp_rate(jap_chr2_CO, jap_chr2_snp)
jap_chr3_snp2 <- snp_rate(jap_chr3_CO, jap_chr3_snp)
jap_chr4_snp2 <- snp_rate(jap_chr4_CO, jap_chr4_snp)
jap_chr5_snp2 <- snp_rate(jap_chr5_CO, jap_chr5_snp)
jap_chr6_snp2 <- snp_rate(jap_chr6_CO, jap_chr6_snp)
jap_chr7_snp2 <- snp_rate(jap_chr7_CO, jap_chr7_snp)
jap_chr8_snp2 <- snp_rate(jap_chr8_CO, jap_chr8_snp)
jap_chr9_snp2 <- snp_rate(jap_chr9_CO, jap_chr9_snp)
jap_chr10_snp2 <- snp_rate(jap_chr10_CO, jap_chr10_snp)
jap_chr11_snp2 <- snp_rate(jap_chr11_CO, jap_chr11_snp)
jap_chr12_snp2 <- snp_rate(jap_chr12_CO, jap_chr12_snp)

jap_chr1_snp2 <- jap_chr1_snp2[order(jap_chr1_snp2$`SNP Start`),]
jap_chr1_spl <- smooth.spline(jap_chr1_snp2$rate, spar = .4)
jap_chr1_snp2$pos <- (jap_chr1_snp2$`SNP Start`*jap_chr1_spl$y)
plot(jap_chr1_snp2$`SNP Start`, jap_chr1_snp2$pos)
ggplot(jap_chr1_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
plot(jap_chr1_snp2$`SNP Start`, jap_chr1_snp2$pos/jap_chr1_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 1 Recombination Distribution")
jap_chr1_finalpos <- jap_chr1_snp2[order(jap_chr1_snp2$pos),]
is.unsorted(jap_chr1_finalpos$pos)
plot(jap_chr1_snp2$`SNP Start`, jap_chr1_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Chromosome 1 Genetic Map")
plot(jap_chr1_finalpos$`SNP Start`, jap_chr1_finalpos$pos)
