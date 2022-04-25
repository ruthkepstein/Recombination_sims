library(AlphaSimR)
library(Rcpp)
library(ggplot2)
library(dplyr)

set.seed(420)

japonica_snps <- read.table("japonica_SNPs.bed", header =FALSE)
colnames(japonica_snps) <- c("Chr#", "SNP Start", "SNP End")
recq4l_snps <- sample_n(japonica_snps, 2000)
recq4l_snps <- recq4l_snps[order(recq4l_snps$`Chr#`,recq4l_snps$`SNP Start`),]

#splitting recq4l snps
recq4l_chr1_snp <- recq4l_snps[ which(recq4l_snps$`Chr#` == "Chr1"),]
recq4l_chr1_snp$rate <- NA
recq4l_chr1_snp$`SNP End` <- recq4l_chr1_snp$`SNP End` - min(recq4l_chr1_snp$`SNP Start`)
recq4l_chr1_snp$`SNP Start` <- recq4l_chr1_snp$`SNP Start`- min(recq4l_chr1_snp$`SNP Start`)

recq4l_chr2_snp <- recq4l_snps[ which(recq4l_snps$`Chr#` == "Chr2"),]
recq4l_chr2_snp$rate <- NA
recq4l_chr2_snp$`SNP End` <- recq4l_chr2_snp$`SNP End` - min(recq4l_chr2_snp$`SNP Start`)
recq4l_chr2_snp$`SNP Start` <- recq4l_chr2_snp$`SNP Start`- min(recq4l_chr2_snp$`SNP Start`)

recq4l_chr3_snp <- recq4l_snps[ which(recq4l_snps$`Chr#` == "Chr3"),]
recq4l_chr3_snp$rate <- NA
recq4l_chr3_snp$`SNP End` <- recq4l_chr3_snp$`SNP End` - min(recq4l_chr3_snp$`SNP Start`)
recq4l_chr3_snp$`SNP Start` <- recq4l_chr3_snp$`SNP Start`- min(recq4l_chr3_snp$`SNP Start`)

recq4l_chr4_snp <- recq4l_snps[ which(recq4l_snps$`Chr#` == "Chr4"),]
recq4l_chr4_snp$rate <- NA
recq4l_chr4_snp$`SNP End` <- recq4l_chr4_snp$`SNP End` - min(recq4l_chr4_snp$`SNP Start`)
recq4l_chr4_snp$`SNP Start` <- recq4l_chr4_snp$`SNP Start`- min(recq4l_chr4_snp$`SNP Start`)

recq4l_chr5_snp <- recq4l_snps[ which(recq4l_snps$`Chr#` == "Chr5"),]
recq4l_chr5_snp$rate <- NA
recq4l_chr5_snp$`SNP End` <- recq4l_chr5_snp$`SNP End` - min(recq4l_chr5_snp$`SNP Start`)
recq4l_chr5_snp$`SNP Start` <- recq4l_chr5_snp$`SNP Start`- min(recq4l_chr5_snp$`SNP Start`)

recq4l_chr6_snp <- recq4l_snps[ which(recq4l_snps$`Chr#` == "Chr6"),]
recq4l_chr6_snp$rate <- NA
recq4l_chr6_snp$`SNP End` <- recq4l_chr6_snp$`SNP End` - min(recq4l_chr6_snp$`SNP Start`)
recq4l_chr6_snp$`SNP Start` <- recq4l_chr6_snp$`SNP Start`- min(recq4l_chr6_snp$`SNP Start`)

recq4l_chr7_snp <- recq4l_snps[ which(recq4l_snps$`Chr#` == "Chr7"),]
recq4l_chr7_snp$rate <- NA
recq4l_chr7_snp$`SNP End` <- recq4l_chr7_snp$`SNP End` - min(recq4l_chr7_snp$`SNP Start`)
recq4l_chr7_snp$`SNP Start` <- recq4l_chr7_snp$`SNP Start`- min(recq4l_chr7_snp$`SNP Start`)

recq4l_chr8_snp <- recq4l_snps[ which(recq4l_snps$`Chr#` == "Chr8"),]
recq4l_chr8_snp$rate <- NA
recq4l_chr8_snp$`SNP End` <- recq4l_chr8_snp$`SNP End` - min(recq4l_chr8_snp$`SNP Start`)
recq4l_chr8_snp$`SNP Start` <- recq4l_chr8_snp$`SNP Start`- min(recq4l_chr8_snp$`SNP Start`)

recq4l_chr9_snp <- recq4l_snps[ which(recq4l_snps$`Chr#` == "Chr9"),]
recq4l_chr9_snp$rate <- NA
recq4l_chr9_snp$`SNP End` <- recq4l_chr9_snp$`SNP End` - min(recq4l_chr9_snp$`SNP Start`)
recq4l_chr9_snp$`SNP Start` <- recq4l_chr9_snp$`SNP Start`- min(recq4l_chr9_snp$`SNP Start`)

recq4l_chr10_snp <- recq4l_snps[ which(recq4l_snps$`Chr#` == "Chr10"),]
recq4l_chr10_snp$rate <- NA
recq4l_chr10_snp$`SNP End` <- recq4l_chr10_snp$`SNP End` - min(recq4l_chr10_snp$`SNP Start`)
recq4l_chr10_snp$`SNP Start` <- recq4l_chr10_snp$`SNP Start`- min(recq4l_chr10_snp$`SNP Start`)

recq4l_chr11_snp <- recq4l_snps[ which(recq4l_snps$`Chr#` == "Chr11"),]
recq4l_chr11_snp$rate <- NA
recq4l_chr11_snp$`SNP End` <- recq4l_chr11_snp$`SNP End` - min(recq4l_chr11_snp$`SNP Start`)
recq4l_chr11_snp$`SNP Start` <- recq4l_chr11_snp$`SNP Start`- min(recq4l_chr11_snp$`SNP Start`)

recq4l_chr12_snp <- recq4l_snps[ which(recq4l_snps$`Chr#` == "Chr12"),]
recq4l_chr12_snp$rate <- NA
recq4l_chr12_snp$`SNP End` <- recq4l_chr12_snp$`SNP End` - min(recq4l_chr12_snp$`SNP Start`)
recq4l_chr12_snp$`SNP Start` <- recq4l_chr12_snp$`SNP Start`- min(recq4l_chr12_snp$`SNP Start`)


#recq4l mutant recomb rates from "unleashing meoitic ... paper"
recq4l_CO <- read.csv("jap_mut_recq4l.csv", header = TRUE)
colnames(recq4l_CO) <- c("Chr", "CO Start", "CO End", "recq4l_rate")
recq4l_CO <- recq4l_CO[order(recq4l_CO$Chr,recq4l_CO$`CO Start`),]

recq4l_chr1_CO <- recq4l_CO[ which(recq4l_CO$Chr == "1"),]
recq4l_chr1_CO$midpoint <- (recq4l_chr1_CO$`CO Start`+ recq4l_chr1_CO$`CO End`)/2
recq4l_chr1_CO <- recq4l_chr1_CO[order(recq4l_chr1_CO$`CO Start`),]

recq4l_chr2_CO <- recq4l_CO[ which(recq4l_CO$Chr == "2"),]
recq4l_chr2_CO$midpoint <- (recq4l_chr2_CO$`CO Start`+ recq4l_chr2_CO$`CO End`)/2
recq4l_chr2_CO <- recq4l_chr2_CO[order(recq4l_chr2_CO$`CO Start`),]

recq4l_chr3_CO <- recq4l_CO[ which(recq4l_CO$Chr == "3"),]
recq4l_chr3_CO$midpoint <- (recq4l_chr3_CO$`CO Start`+ recq4l_chr3_CO$`CO End`)/2
recq4l_chr3_CO <- recq4l_chr3_CO[order(recq4l_chr3_CO$`CO Start`),]

recq4l_chr4_CO <- recq4l_CO[ which(recq4l_CO$Chr == "4"),]
recq4l_chr4_CO$midpoint <- (recq4l_chr4_CO$`CO Start`+ recq4l_chr4_CO$`CO End`)/2
recq4l_chr4_CO <- recq4l_chr4_CO[order(recq4l_chr4_CO$`CO Start`),]

recq4l_chr5_CO <- recq4l_CO[ which(recq4l_CO$Chr == "5"),]
recq4l_chr5_CO$midpoint <- (recq4l_chr5_CO$`CO Start`+ recq4l_chr5_CO$`CO End`)/2
recq4l_chr5_CO <- recq4l_chr5_CO[order(recq4l_chr5_CO$`CO Start`),]

recq4l_chr6_CO <- recq4l_CO[ which(recq4l_CO$Chr == "6"),]
recq4l_chr6_CO$midpoint <- (recq4l_chr6_CO$`CO Start`+ recq4l_chr6_CO$`CO End`)/2
recq4l_chr6_CO <- recq4l_chr6_CO[order(recq4l_chr6_CO$`CO Start`),]

recq4l_chr7_CO <- recq4l_CO[ which(recq4l_CO$Chr == "7"),]
recq4l_chr7_CO$midpoint <- (recq4l_chr7_CO$`CO Start`+ recq4l_chr7_CO$`CO End`)/2
recq4l_chr7_CO <- recq4l_chr7_CO[order(recq4l_chr7_CO$`CO Start`),]

recq4l_chr8_CO <- recq4l_CO[ which(recq4l_CO$Chr == "8"),]
recq4l_chr8_CO$midpoint <- (recq4l_chr8_CO$`CO Start`+ recq4l_chr8_CO$`CO End`)/2
recq4l_chr8_CO <- recq4l_chr8_CO[order(recq4l_chr8_CO$`CO Start`),]

recq4l_chr9_CO <- recq4l_CO[ which(recq4l_CO$Chr == "9"),]
recq4l_chr9_CO$midpoint <- (recq4l_chr9_CO$`CO Start`+ recq4l_chr9_CO$`CO End`)/2
recq4l_chr9_CO <- recq4l_chr9_CO[order(recq4l_chr9_CO$`CO Start`),]

recq4l_chr10_CO <- recq4l_CO[ which(recq4l_CO$Chr == "10"),]
recq4l_chr10_CO$midpoint <- (recq4l_chr10_CO$`CO Start`+ recq4l_chr10_CO$`CO End`)/2
recq4l_chr10_CO <- recq4l_chr10_CO[order(recq4l_chr10_CO$`CO Start`),]

recq4l_chr11_CO <- recq4l_CO[ which(recq4l_CO$Chr == "11"),]
recq4l_chr11_CO$midpoint <- (recq4l_chr11_CO$`CO Start`+ recq4l_chr11_CO$`CO End`)/2
recq4l_chr11_CO <- recq4l_chr11_CO[order(recq4l_chr11_CO$`CO Start`),]

recq4l_chr12_CO <- recq4l_CO[ which(recq4l_CO$Chr == "12"),]
recq4l_chr12_CO$midpoint <- (recq4l_chr12_CO$`CO Start`+ recq4l_chr12_CO$`CO End`)/2
recq4l_chr12_CO <- recq4l_chr12_CO[order(recq4l_chr12_CO$`CO Start`),]


#calculating recombination rate per bin of CO data
library(dlookr)
library(tidyverse)
library(OneR)

#making intervals start at 0
recq4l_chr1_CO$`CO End` <- recq4l_chr1_CO$`CO End` - min(recq4l_chr1_CO$`CO Start`)
recq4l_chr1_CO$`CO Start` <- recq4l_chr1_CO$`CO Start` - min(recq4l_chr1_CO$`CO Start`)

recq4l_chr2_CO$`CO End` <- recq4l_chr2_CO$`CO End` - min(recq4l_chr2_CO$`CO Start`)
recq4l_chr2_CO$`CO Start` <- recq4l_chr2_CO$`CO Start` - min(recq4l_chr2_CO$`CO Start`)

recq4l_chr3_CO$`CO End` <- recq4l_chr3_CO$`CO End` - min(recq4l_chr3_CO$`CO Start`)
recq4l_chr3_CO$`CO Start` <- recq4l_chr3_CO$`CO Start` - min(recq4l_chr3_CO$`CO Start`)

recq4l_chr4_CO$`CO End` <- recq4l_chr4_CO$`CO End` - min(recq4l_chr4_CO$`CO Start`)
recq4l_chr4_CO$`CO Start` <- recq4l_chr4_CO$`CO Start` - min(recq4l_chr4_CO$`CO Start`)

recq4l_chr5_CO$`CO End` <- recq4l_chr5_CO$`CO End` - min(recq4l_chr5_CO$`CO Start`)
recq4l_chr5_CO$`CO Start` <- recq4l_chr5_CO$`CO Start` - min(recq4l_chr5_CO$`CO Start`)

recq4l_chr6_CO$`CO End` <- recq4l_chr6_CO$`CO End` - min(recq4l_chr6_CO$`CO Start`)
recq4l_chr6_CO$`CO Start` <- recq4l_chr6_CO$`CO Start` - min(recq4l_chr6_CO$`CO Start`)

recq4l_chr7_CO$`CO End` <- recq4l_chr7_CO$`CO End` - min(recq4l_chr7_CO$`CO Start`)
recq4l_chr7_CO$`CO Start` <- recq4l_chr7_CO$`CO Start` - min(recq4l_chr7_CO$`CO Start`)

recq4l_chr8_CO$`CO End` <- recq4l_chr8_CO$`CO End` - min(recq4l_chr8_CO$`CO Start`)
recq4l_chr8_CO$`CO Start` <- recq4l_chr8_CO$`CO Start` - min(recq4l_chr8_CO$`CO Start`)

recq4l_chr9_CO$`CO End` <- recq4l_chr9_CO$`CO End` - min(recq4l_chr9_CO$`CO Start`)
recq4l_chr9_CO$`CO Start` <- recq4l_chr9_CO$`CO Start` - min(recq4l_chr9_CO$`CO Start`)

recq4l_chr10_CO$`CO End` <- recq4l_chr10_CO$`CO End` - min(recq4l_chr10_CO$`CO Start`)
recq4l_chr10_CO$`CO Start` <- recq4l_chr10_CO$`CO Start` - min(recq4l_chr10_CO$`CO Start`)

recq4l_chr11_CO$`CO End` <- recq4l_chr11_CO$`CO End` - min(recq4l_chr11_CO$`CO Start`)
recq4l_chr11_CO$`CO Start` <- recq4l_chr11_CO$`CO Start` - min(recq4l_chr11_CO$`CO Start`)

recq4l_chr12_CO$`CO End` <- recq4l_chr12_CO$`CO End` - min(recq4l_chr12_CO$`CO Start`)
recq4l_chr12_CO$`CO Start` <- recq4l_chr12_CO$`CO Start` - min(recq4l_chr12_CO$`CO Start`)


## multiplying WT recombination rates by the avg difference
# exclude pericentromeric regions (suppresion regions)
# 1. create avg diff column (supression region = 0)
# 2. loop through, multiply wt rate from other paper (fine scale recombination rate) by avg rate
#japonica wildtype recomb rates
jap_CO <- read.csv("jap_WT_rate.csv", header = TRUE)
colnames(jap_CO) <- c("Chr", "CO Start", "CO End", "WT_rate")
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

#apply avg difference to telomeric regions (divide chromosome into fifths and apply avg diff to first & last fifth)
avg_diff <- 1.581055
pericentromeric <- function(CO){
  rownames(CO)<-c(1:nrow(CO))
  CO$avg_rate <- avg_diff
  fifth<- max(CO$`CO End`)/5
  start<-fifth*2
  end<-fifth*4
  for(i in 1:nrow(CO)){
    if(CO$`CO Start`[i]>= start && CO$`CO End`[i]<= end ){
      CO$avg_rate[i] <- 0
    }
  }
  print(CO)
}
jap_chr1_CO<- pericentromeric(jap_chr1_CO)
jap_chr2_CO <- pericentromeric(jap_chr2_CO)
jap_chr3_CO <- pericentromeric(jap_chr3_CO)
jap_chr4_CO <- pericentromeric(jap_chr4_CO)
jap_chr5_CO <- pericentromeric(jap_chr5_CO)
jap_chr6_CO <- pericentromeric(jap_chr6_CO)
jap_chr7_CO <- pericentromeric(jap_chr7_CO)
jap_chr8_CO <- pericentromeric(jap_chr8_CO)
jap_chr9_CO <- pericentromeric(jap_chr9_CO)
jap_chr10_CO <- pericentromeric(jap_chr10_CO)
jap_chr11_CO <- pericentromeric(jap_chr11_CO)
jap_chr12_CO <- pericentromeric(jap_chr12_CO)
#import fine scale recombination rates
WTJap_CO <- read.table("japonica_rec_rate.bed", header = FALSE)
colnames(WTJap_CO) <- c("Chr", "CO Start", "CO End", "rate")
WTJap_CO <- WTJap_CO[order(WTJap_CO$Chr,WTJap_CO$`CO Start`),]

WTJap_chr1_CO <- WTJap_CO[ which(WTJap_CO$Chr == "chr01"),]
WTJap_chr1_CO$midpoint <- (WTJap_chr1_CO$`CO Start`+ WTJap_chr1_CO$`CO End`)/2
WTJap_chr1_CO <- WTJap_chr1_CO[order(WTJap_chr1_CO$`CO Start`),]

WTJap_chr2_CO <- WTJap_CO[ which(WTJap_CO$Chr == "chr02"),]
WTJap_chr2_CO$midpoint <- (WTJap_chr2_CO$`CO Start`+ WTJap_chr2_CO$`CO End`)/2
WTJap_chr2_CO <- WTJap_chr2_CO[order(WTJap_chr2_CO$`CO Start`),]

WTJap_chr3_CO <- WTJap_CO[ which(WTJap_CO$Chr == "chr03"),]
WTJap_chr3_CO$midpoint <- (WTJap_chr3_CO$`CO Start`+ WTJap_chr3_CO$`CO End`)/2
WTJap_chr3_CO <- WTJap_chr3_CO[order(WTJap_chr3_CO$`CO Start`),]

WTJap_chr4_CO <- WTJap_CO[ which(WTJap_CO$Chr == "chr04"),]
WTJap_chr4_CO$midpoint <- (WTJap_chr4_CO$`CO Start`+ WTJap_chr4_CO$`CO End`)/2
WTJap_chr4_CO <- WTJap_chr4_CO[order(WTJap_chr4_CO$`CO Start`),]

WTJap_chr5_CO <- WTJap_CO[ which(WTJap_CO$Chr == "chr05"),]
WTJap_chr5_CO$midpoint <- (WTJap_chr5_CO$`CO Start`+ WTJap_chr5_CO$`CO End`)/2
WTJap_chr5_CO <- WTJap_chr5_CO[order(WTJap_chr5_CO$`CO Start`),]

WTJap_chr6_CO <- WTJap_CO[ which(WTJap_CO$Chr == "chr06"),]
WTJap_chr6_CO$midpoint <- (WTJap_chr6_CO$`CO Start`+ WTJap_chr6_CO$`CO End`)/2
WTJap_chr6_CO <- WTJap_chr6_CO[order(WTJap_chr6_CO$`CO Start`),]

WTJap_chr7_CO <- WTJap_CO[ which(WTJap_CO$Chr == "chr07"),]
WTJap_chr7_CO$midpoint <- (WTJap_chr7_CO$`CO Start`+ WTJap_chr7_CO$`CO End`)/2
WTJap_chr7_CO <- WTJap_chr7_CO[order(WTJap_chr7_CO$`CO Start`),]

WTJap_chr8_CO <- WTJap_CO[ which(WTJap_CO$Chr == "chr08"),]
WTJap_chr8_CO$midpoint <- (WTJap_chr8_CO$`CO Start`+ WTJap_chr8_CO$`CO End`)/2
WTJap_chr8_CO <- WTJap_chr8_CO[order(WTJap_chr8_CO$`CO Start`),]

WTJap_chr9_CO <- WTJap_CO[ which(WTJap_CO$Chr == "chr09"),]
WTJap_chr9_CO$midpoint <- (WTJap_chr9_CO$`CO Start`+ WTJap_chr9_CO$`CO End`)/2
WTJap_chr9_CO <- WTJap_chr9_CO[order(WTJap_chr9_CO$`CO Start`),]

WTJap_chr10_CO <- WTJap_CO[ which(WTJap_CO$Chr == "chr10"),]
WTJap_chr10_CO$midpoint <- (WTJap_chr10_CO$`CO Start`+ WTJap_chr10_CO$`CO End`)/2
WTJap_chr10_CO <- WTJap_chr10_CO[order(WTJap_chr10_CO$`CO Start`),]

WTJap_chr11_CO <- WTJap_CO[ which(WTJap_CO$Chr == "chr11"),]
WTJap_chr11_CO$midpoint <- (WTJap_chr11_CO$`CO Start`+ WTJap_chr11_CO$`CO End`)/2
WTJap_chr11_CO <- WTJap_chr11_CO[order(WTJap_chr11_CO$`CO Start`),]

WTJap_chr12_CO <- WTJap_CO[ which(WTJap_CO$Chr == "chr12"),]
WTJap_chr12_CO$midpoint <- (WTJap_chr12_CO$`CO Start`+ WTJap_chr12_CO$`CO End`)/2
WTJap_chr12_CO <- WTJap_chr12_CO[order(WTJap_chr12_CO$`CO Start`),]

WTJap_chr1_CO$`CO End` <- WTJap_chr1_CO$`CO End` - min(WTJap_chr1_CO$`CO Start`)
WTJap_chr1_CO$`CO Start` <- WTJap_chr1_CO$`CO Start` - min(WTJap_chr1_CO$`CO Start`)

WTJap_chr2_CO$`CO End` <- WTJap_chr2_CO$`CO End` - min(WTJap_chr2_CO$`CO Start`)
WTJap_chr2_CO$`CO Start` <- WTJap_chr2_CO$`CO Start` - min(WTJap_chr2_CO$`CO Start`)

WTJap_chr3_CO$`CO End` <- WTJap_chr3_CO$`CO End` - min(WTJap_chr3_CO$`CO Start`)
WTJap_chr3_CO$`CO Start` <- WTJap_chr3_CO$`CO Start` - min(WTJap_chr3_CO$`CO Start`)

WTJap_chr4_CO$`CO End` <- WTJap_chr4_CO$`CO End` - min(WTJap_chr4_CO$`CO Start`)
WTJap_chr4_CO$`CO Start` <- WTJap_chr4_CO$`CO Start` - min(WTJap_chr4_CO$`CO Start`)

WTJap_chr5_CO$`CO End` <- WTJap_chr5_CO$`CO End` - min(WTJap_chr5_CO$`CO Start`)
WTJap_chr5_CO$`CO Start` <- WTJap_chr5_CO$`CO Start` - min(WTJap_chr5_CO$`CO Start`)

WTJap_chr6_CO$`CO End` <- WTJap_chr6_CO$`CO End` - min(WTJap_chr6_CO$`CO Start`)
WTJap_chr6_CO$`CO Start` <- WTJap_chr6_CO$`CO Start` - min(WTJap_chr6_CO$`CO Start`)

WTJap_chr7_CO$`CO End` <- WTJap_chr7_CO$`CO End` - min(WTJap_chr7_CO$`CO Start`)
WTJap_chr7_CO$`CO Start` <- WTJap_chr7_CO$`CO Start` - min(WTJap_chr7_CO$`CO Start`)

WTJap_chr8_CO$`CO End` <- WTJap_chr8_CO$`CO End` - min(WTJap_chr8_CO$`CO Start`)
WTJap_chr8_CO$`CO Start` <- WTJap_chr8_CO$`CO Start` - min(WTJap_chr8_CO$`CO Start`)

WTJap_chr9_CO$`CO End` <- WTJap_chr9_CO$`CO End` - min(WTJap_chr9_CO$`CO Start`)
WTJap_chr9_CO$`CO Start` <- WTJap_chr9_CO$`CO Start` - min(WTJap_chr9_CO$`CO Start`)

WTJap_chr10_CO$`CO End` <- WTJap_chr10_CO$`CO End` - min(WTJap_chr10_CO$`CO Start`)
WTJap_chr10_CO$`CO Start` <- WTJap_chr10_CO$`CO Start` - min(WTJap_chr10_CO$`CO Start`)

WTJap_chr11_CO$`CO End` <- WTJap_chr11_CO$`CO End` - min(WTJap_chr11_CO$`CO Start`)
WTJap_chr11_CO$`CO Start` <- WTJap_chr11_CO$`CO Start` - min(WTJap_chr11_CO$`CO Start`)

WTJap_chr12_CO$`CO End` <- WTJap_chr12_CO$`CO End` - min(WTJap_chr12_CO$`CO Start`)
WTJap_chr12_CO$`CO Start` <- WTJap_chr12_CO$`CO Start` - min(WTJap_chr12_CO$`CO Start`)

## Multiply recombination fine scale data by the avg rate
new_rates <- function(old_rate){
  for(i in 1:nrow(old_rate)){
    old_rate$rate[i] <- old_rate$rate[i] + (avg_diff*old_rate$rate[i])
  }
  print(old_rate)
}
recq4l_chr1_CO_2<-new_rates(WTJap_chr1_CO)
recq4l_chr2_CO_2<-new_rates(WTJap_chr2_CO)
recq4l_chr3_CO_2<-new_rates(WTJap_chr3_CO)
recq4l_chr4_CO_2<-new_rates(WTJap_chr4_CO)
recq4l_chr5_CO_2<-new_rates(WTJap_chr5_CO)
recq4l_chr6_CO_2<-new_rates(WTJap_chr6_CO)
recq4l_chr7_CO_2<-new_rates(WTJap_chr7_CO)
recq4l_chr8_CO_2<-new_rates(WTJap_chr8_CO)
recq4l_chr9_CO_2<-new_rates(WTJap_chr9_CO)
recq4l_chr10_CO_2<-new_rates(WTJap_chr10_CO)
recq4l_chr11_CO_2<-new_rates(WTJap_chr11_CO)
recq4l_chr12_CO_2<-new_rates(WTJap_chr12_CO)

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
recq4l_chr1_CO_3 <- recq4l_chr1_CO_2
bins<-as.integer(nrow(recq4l_chr1_CO_2)/40)
recq4l_chr1_CO_3$rates<- rollapply(recq4l_chr1_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
recq4l_chr1_CO_3<-fill_start(recq4l_chr1_CO_3)
recq4l_chr1_CO_3<- recq4l_chr1_CO_3 %>% drop_na(rates)

recq4l_chr2_CO_3 <- recq4l_chr2_CO_2
bins<-as.integer(nrow(recq4l_chr2_CO_2)/40)
recq4l_chr2_CO_3$rates<- rollapply(recq4l_chr2_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
recq4l_chr2_CO_3<-fill_start(recq4l_chr2_CO_3)
recq4l_chr2_CO_3<- recq4l_chr2_CO_3 %>% drop_na(rates)

recq4l_chr3_CO_3 <- recq4l_chr3_CO_2
bins<-as.integer(nrow(recq4l_chr3_CO_2)/40)
recq4l_chr3_CO_3$rates<- rollapply(recq4l_chr3_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
recq4l_chr3_CO_3<-fill_start(recq4l_chr3_CO_3)
recq4l_chr3_CO_3<- recq4l_chr3_CO_3 %>% drop_na(rates)

recq4l_chr4_CO_3 <- recq4l_chr4_CO_2
bins<-as.integer(nrow(recq4l_chr4_CO_2)/40)
recq4l_chr4_CO_3$rates<- rollapply(recq4l_chr4_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
recq4l_chr4_CO_3<-fill_start(recq4l_chr4_CO_3)
recq4l_chr4_CO_3<- recq4l_chr4_CO_3 %>% drop_na(rates)

recq4l_chr5_CO_3 <- recq4l_chr5_CO_2
bins<-as.integer(nrow(recq4l_chr5_CO_2)/40)
recq4l_chr5_CO_3$rates<- rollapply(recq4l_chr5_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
recq4l_chr5_CO_3<-fill_start(recq4l_chr5_CO_3)
recq4l_chr5_CO_3<- recq4l_chr5_CO_3 %>% drop_na(rates)

recq4l_chr6_CO_3 <- recq4l_chr6_CO_2
bins<-as.integer(nrow(recq4l_chr6_CO_2)/40)
recq4l_chr6_CO_3$rates<- rollapply(recq4l_chr6_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
recq4l_chr6_CO_3<-fill_start(recq4l_chr6_CO_3)
recq4l_chr6_CO_3<- recq4l_chr6_CO_3 %>% drop_na(rates)

recq4l_chr7_CO_3 <- recq4l_chr7_CO_2
bins<-as.integer(nrow(recq4l_chr7_CO_2)/40)
recq4l_chr7_CO_3$rates<- rollapply(recq4l_chr7_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
recq4l_chr7_CO_3<-fill_start(recq4l_chr7_CO_3)
recq4l_chr7_CO_3<- recq4l_chr7_CO_3 %>% drop_na(rates)

recq4l_chr8_CO_3 <- recq4l_chr8_CO_2
bins<-as.integer(nrow(recq4l_chr8_CO_2)/40)
recq4l_chr8_CO_3$rates<- rollapply(recq4l_chr8_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
recq4l_chr8_CO_3<-fill_start(recq4l_chr8_CO_3)
recq4l_chr8_CO_3<- recq4l_chr8_CO_3 %>% drop_na(rates)

recq4l_chr9_CO_3 <- recq4l_chr9_CO_2
bins<-as.integer(nrow(recq4l_chr9_CO_2)/40)
recq4l_chr9_CO_3$rates<- rollapply(recq4l_chr9_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
recq4l_chr9_CO_3<-fill_start(recq4l_chr9_CO_3)
recq4l_chr9_CO_3<- recq4l_chr9_CO_3 %>% drop_na(rates)

recq4l_chr10_CO_3 <- recq4l_chr10_CO_2
bins<-as.integer(nrow(recq4l_chr10_CO_2)/40)
recq4l_chr10_CO_3$rates<- rollapply(recq4l_chr10_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
recq4l_chr10_CO_3<-fill_start(recq4l_chr10_CO_3)
recq4l_chr10_CO_3<- recq4l_chr10_CO_3 %>% drop_na(rates)

recq4l_chr11_CO_3 <- recq4l_chr11_CO_2
bins<-as.integer(nrow(recq4l_chr11_CO_2)/40)
recq4l_chr11_CO_3$rates<- rollapply(recq4l_chr11_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
recq4l_chr11_CO_3<-fill_start(recq4l_chr11_CO_3)
recq4l_chr11_CO_3<- recq4l_chr11_CO_3 %>% drop_na(rates)

recq4l_chr12_CO_3 <- recq4l_chr12_CO_2
bins<-as.integer(nrow(recq4l_chr12_CO_2)/40)
recq4l_chr12_CO_3$rates<- rollapply(recq4l_chr12_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
recq4l_chr12_CO_3<-fill_start(recq4l_chr12_CO_3)
recq4l_chr12_CO_3<- recq4l_chr12_CO_3 %>% drop_na(rates)


## assigning frequency to SNPs based on recombination frequency in each bin
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


#using function,  get cM/Mb for final genetic position - assign rates
recq4l_chr1_snp2 <- snp_rate(recq4l_chr1_CO_3, recq4l_chr1_snp)
recq4l_chr2_snp2 <- snp_rate(recq4l_chr2_CO_3, recq4l_chr2_snp)
recq4l_chr3_snp2 <- snp_rate(recq4l_chr3_CO_3, recq4l_chr3_snp)
recq4l_chr4_snp2 <- snp_rate(recq4l_chr4_CO_3, recq4l_chr4_snp)
recq4l_chr5_snp2 <- snp_rate(recq4l_chr5_CO_3, recq4l_chr5_snp)
recq4l_chr6_snp2 <- snp_rate(recq4l_chr6_CO_3, recq4l_chr6_snp)
recq4l_chr7_snp2 <- snp_rate(recq4l_chr7_CO_3, recq4l_chr7_snp)
recq4l_chr8_snp2 <- snp_rate(recq4l_chr8_CO_3, recq4l_chr8_snp)
recq4l_chr9_snp2 <- snp_rate(recq4l_chr9_CO_3, recq4l_chr9_snp)
recq4l_chr10_snp2 <- snp_rate(recq4l_chr10_CO_3, recq4l_chr10_snp)
recq4l_chr11_snp2 <- snp_rate(recq4l_chr11_CO_3, recq4l_chr11_snp)
recq4l_chr12_snp2 <- snp_rate(recq4l_chr12_CO_3, recq4l_chr12_snp)

#converted SNP start to Mb
recq4l_chr1_snp2$`SNP Start` <- recq4l_chr1_snp2$`SNP Start`/1000000
recq4l_chr2_snp2$`SNP Start` <- recq4l_chr2_snp2$`SNP Start`/1000000
recq4l_chr3_snp2$`SNP Start` <- recq4l_chr3_snp2$`SNP Start`/1000000
recq4l_chr4_snp2$`SNP Start` <- recq4l_chr4_snp2$`SNP Start`/1000000
recq4l_chr5_snp2$`SNP Start` <- recq4l_chr5_snp2$`SNP Start`/1000000
recq4l_chr6_snp2$`SNP Start` <- recq4l_chr6_snp2$`SNP Start`/1000000
recq4l_chr7_snp2$`SNP Start` <- recq4l_chr7_snp2$`SNP Start`/1000000
recq4l_chr8_snp2$`SNP Start` <- recq4l_chr8_snp2$`SNP Start`/1000000
recq4l_chr9_snp2$`SNP Start` <- recq4l_chr9_snp2$`SNP Start`/1000000
recq4l_chr10_snp2$`SNP Start` <- recq4l_chr10_snp2$`SNP Start`/1000000
recq4l_chr11_snp2$`SNP Start` <- recq4l_chr11_snp2$`SNP Start`/1000000
recq4l_chr12_snp2$`SNP Start` <- recq4l_chr12_snp2$`SNP Start`/1000000

recq4l_chr1_snp2$`SNP End` <- recq4l_chr1_snp2$`SNP End`/1000000
recq4l_chr2_snp2$`SNP End` <- recq4l_chr2_snp2$`SNP End`/1000000
recq4l_chr3_snp2$`SNP End` <- recq4l_chr3_snp2$`SNP End`/1000000
recq4l_chr4_snp2$`SNP End` <- recq4l_chr4_snp2$`SNP End`/1000000
recq4l_chr5_snp2$`SNP End` <- recq4l_chr5_snp2$`SNP End`/1000000
recq4l_chr6_snp2$`SNP End` <- recq4l_chr6_snp2$`SNP End`/1000000
recq4l_chr7_snp2$`SNP End` <- recq4l_chr7_snp2$`SNP End`/1000000
recq4l_chr8_snp2$`SNP End` <- recq4l_chr8_snp2$`SNP End`/1000000
recq4l_chr9_snp2$`SNP End` <- recq4l_chr9_snp2$`SNP End`/1000000
recq4l_chr10_snp2$`SNP End` <- recq4l_chr10_snp2$`SNP End`/1000000
recq4l_chr11_snp2$`SNP End` <- recq4l_chr11_snp2$`SNP End`/1000000
recq4l_chr12_snp2$`SNP End` <- recq4l_chr12_snp2$`SNP End`/1000000

#omit empty col
recq4l_chr1_snp2<-na.omit(recq4l_chr1_snp2)
recq4l_chr2_snp2<-na.omit(recq4l_chr2_snp2)
recq4l_chr3_snp2<-na.omit(recq4l_chr3_snp2)
recq4l_chr4_snp2<-na.omit(recq4l_chr4_snp2)
recq4l_chr5_snp2<-na.omit(recq4l_chr5_snp2)
recq4l_chr6_snp2<-na.omit(recq4l_chr6_snp2)
recq4l_chr7_snp2<-na.omit(recq4l_chr7_snp2)
recq4l_chr8_snp2<-na.omit(recq4l_chr8_snp2)
recq4l_chr9_snp2<-na.omit(recq4l_chr9_snp2)
recq4l_chr10_snp2<-na.omit(recq4l_chr10_snp2)
recq4l_chr11_snp2<-na.omit(recq4l_chr11_snp2)
recq4l_chr12_snp2<-na.omit(recq4l_chr12_snp2)


#gen maps
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


recq4l_chr1_spl <- smooth.spline(recq4l_chr1_snp2$rate, spar = 0.1)
recq4l_chr1_snp2$pos <- gen_pos(recq4l_chr1_snp2)
plot(recq4l_chr1_snp2$`SNP Start`, recq4l_chr1_snp2$pos)
ggplot(recq4l_chr1_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
plot(recq4l_chr1_snp2$`SNP Start`, recq4l_chr1_snp2$pos/recq4l_chr1_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Recq4l Chromosome 1 Recombination Distribution")
recq4l_chr1_finalpos <- recq4l_chr1_snp2[order(recq4l_chr1_snp2$pos),]
is.unsorted(recq4l_chr1_finalpos$pos)
plot(recq4l_chr1_snp2$`SNP Start`, recq4l_chr1_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Recq4l Chromosome 1 Genetic Map")
plot(recq4l_chr1_finalpos$`SNP Start`, recq4l_chr1_finalpos$pos)

recq4l_chr2_spl <- smooth.spline(recq4l_chr2_snp2$rate, spar = 0.1)
recq4l_chr2_snp2$pos <- gen_pos(recq4l_chr2_snp2)
plot(recq4l_chr2_snp2$`SNP Start`, recq4l_chr2_snp2$pos)
plot(recq4l_chr2_snp2$`SNP Start`, recq4l_chr2_snp2$pos/recq4l_chr2_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Recq4l Chromosome 2 Recombination Distribution")
recq4l_chr2_finalpos <- recq4l_chr2_snp2[order(recq4l_chr2_snp2$pos),]
is.unsorted(recq4l_chr2_finalpos$pos)
plot(recq4l_chr2_snp2$`SNP Start`, recq4l_chr2_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Recq4l Chromosome 2 Genetic Map")

recq4l_chr3_spl <- smooth.spline(recq4l_chr3_snp2$rate, spar = 0.1)
recq4l_chr3_snp2$pos <- gen_pos(recq4l_chr3_snp2)
plot(recq4l_chr3_snp2$`SNP Start`, recq4l_chr3_snp2$pos)
plot(recq4l_chr3_snp2$`SNP Start`, recq4l_chr3_snp2$pos/recq4l_chr3_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Recq4l Chromosome 3 Recombination Distribution")
recq4l_chr3_finalpos <- recq4l_chr3_snp2[order(recq4l_chr3_snp2$pos),]
is.unsorted(recq4l_chr3_finalpos$pos)
plot(recq4l_chr3_snp2$`SNP Start`, recq4l_chr3_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Recq4l Chromosome 3 Genetic Map")

recq4l_chr4_spl <- smooth.spline(recq4l_chr4_snp2$rate, spar =0.1)
recq4l_chr4_snp2$pos <- gen_pos(recq4l_chr4_snp2)
plot(recq4l_chr4_snp2$`SNP Start`, recq4l_chr4_snp2$pos)
plot(recq4l_chr4_snp2$`SNP Start`, recq4l_chr4_snp2$pos/recq4l_chr4_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Recq4l Chromosome 4 Recombination Distribution")
recq4l_chr4_finalpos <- recq4l_chr4_snp2[order(recq4l_chr4_snp2$pos),]
is.unsorted(recq4l_chr4_finalpos$pos)
plot(recq4l_chr4_snp2$`SNP Start`, recq4l_chr4_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Recq4l Chromosome 4 Genetic Map")

recq4l_chr5_spl <- smooth.spline(recq4l_chr5_snp2$rate, spar =0.1)
recq4l_chr5_snp2$pos <- gen_pos(recq4l_chr5_snp2)
plot(recq4l_chr5_snp2$`SNP Start`, recq4l_chr5_snp2$pos)
plot(recq4l_chr5_snp2$`SNP Start`, recq4l_chr5_snp2$pos/recq4l_chr5_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Recq4l Chromosome 5 Recombination Distribution")
recq4l_chr5_finalpos <- recq4l_chr5_snp2[order(recq4l_chr5_snp2$pos),]
is.unsorted(recq4l_chr5_finalpos$pos)
plot(recq4l_chr5_snp2$`SNP Start`, recq4l_chr5_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Recq4l Chromosome 5 Genetic Map")

recq4l_chr6_spl <- smooth.spline(recq4l_chr6_snp2$rate, spar = 0.1)
recq4l_chr6_snp2$pos <- gen_pos(recq4l_chr6_snp2)
plot(recq4l_chr6_snp2$`SNP Start`, recq4l_chr6_snp2$pos)
plot(recq4l_chr6_snp2$`SNP Start`, recq4l_chr6_snp2$pos/recq4l_chr6_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Recq4l Chromosome 6 Recombination Distribution")
recq4l_chr6_finalpos <- recq4l_chr6_snp2[order(recq4l_chr6_snp2$pos),]
is.unsorted(recq4l_chr6_finalpos$pos)
plot(recq4l_chr6_snp2$`SNP Start`, recq4l_chr6_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Recq4l Chromosome 6 Genetic Map")

recq4l_chr7_spl <- smooth.spline(recq4l_chr7_snp2$rate, spar =0.1)
recq4l_chr7_snp2$pos <- gen_pos(recq4l_chr7_snp2)
plot(recq4l_chr7_snp2$`SNP Start`, recq4l_chr7_snp2$pos)
plot(recq4l_chr7_snp2$`SNP Start`, recq4l_chr7_snp2$pos/recq4l_chr7_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Recq4l Chromosome 7 Recombination Distribution")
recq4l_chr7_finalpos <- recq4l_chr7_snp2[order(recq4l_chr7_snp2$pos),]
is.unsorted(recq4l_chr7_finalpos$pos)
plot(recq4l_chr7_snp2$`SNP Start`, recq4l_chr7_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Recq4l Chromosome 7 Genetic Map")

recq4l_chr8_spl <- smooth.spline(recq4l_chr8_snp2$rate, spar = 0.1)
recq4l_chr8_snp2$pos <- gen_pos(recq4l_chr8_snp2)
plot(recq4l_chr8_snp2$`SNP Start`, recq4l_chr8_snp2$pos)
plot(recq4l_chr8_snp2$`SNP Start`, recq4l_chr8_snp2$pos/recq4l_chr8_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Recq4l Chromosome 8 Recombination Distribution")
recq4l_chr8_finalpos <- recq4l_chr8_snp2[order(recq4l_chr8_snp2$pos),]
is.unsorted(recq4l_chr8_finalpos$pos)
plot(recq4l_chr8_snp2$`SNP Start`, recq4l_chr8_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Recq4l Chromosome 8 Genetic Map")

recq4l_chr9_spl <- smooth.spline(recq4l_chr9_snp2$rate, spar = 0.1)
recq4l_chr9_snp2$pos <- gen_pos(recq4l_chr9_snp2)
plot(recq4l_chr9_snp2$`SNP Start`, recq4l_chr9_snp2$pos)
plot(recq4l_chr9_snp2$`SNP Start`, recq4l_chr9_snp2$pos/recq4l_chr9_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Recq4l Chromosome 9 Recombination Distribution")
recq4l_chr9_finalpos <- recq4l_chr9_snp2[order(recq4l_chr9_snp2$pos),]
is.unsorted(recq4l_chr9_finalpos$pos)
plot(recq4l_chr9_snp2$`SNP Start`, recq4l_chr9_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Recq4l Chromosome 9 Genetic Map")

recq4l_chr10_spl <- smooth.spline(recq4l_chr10_snp2$rate, spar =0.1)
recq4l_chr10_snp2$pos <- gen_pos(recq4l_chr10_snp2)
plot(recq4l_chr10_snp2$`SNP Start`, recq4l_chr10_snp2$pos)
plot(recq4l_chr10_snp2$`SNP Start`, recq4l_chr10_snp2$pos/recq4l_chr10_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Recq4l Chromosome 10 Recombination Distribution")
recq4l_chr10_finalpos <- recq4l_chr10_snp2[order(recq4l_chr10_snp2$pos),]
is.unsorted(recq4l_chr10_finalpos$pos)
plot(recq4l_chr10_snp2$`SNP Start`, recq4l_chr10_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Recq4l Chromosome 10 Genetic Map")

recq4l_chr11_spl <- smooth.spline(recq4l_chr11_snp2$rate, spar = 0.1)
recq4l_chr11_snp2$pos <- gen_pos(recq4l_chr11_snp2)
plot(recq4l_chr11_snp2$`SNP Start`, recq4l_chr11_snp2$pos)
plot(recq4l_chr11_snp2$`SNP Start`, recq4l_chr11_snp2$pos/recq4l_chr11_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Recq4l Chromosome 11 Recombination Distribution")
recq4l_chr11_finalpos <- recq4l_chr11_snp2[order(recq4l_chr11_snp2$pos),]
is.unsorted(recq4l_chr11_finalpos$pos)
plot(recq4l_chr11_snp2$`SNP Start`, recq4l_chr11_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Recq4l Chromosome 11 Genetic Map")

recq4l_chr12_spl <- smooth.spline(recq4l_chr12_snp2$rate, spar = 0.1)
recq4l_chr12_snp2$pos <- gen_pos(recq4l_chr12_snp2)
plot(recq4l_chr12_snp2$`SNP Start`, recq4l_chr12_snp2$pos)
plot(recq4l_chr12_snp2$`SNP Start`, recq4l_chr12_snp2$pos/recq4l_chr12_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Recq4l Chromosome 12 Recombination Distribution")
recq4l_chr12_finalpos <- recq4l_chr12_snp2[order(recq4l_chr12_snp2$pos),]
is.unsorted(recq4l_chr12_finalpos$pos)
plot(recq4l_chr12_snp2$`SNP Start`, recq4l_chr12_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Recq4l Chromosome 12 Genetic Map")

#Final genetic map
recq4l_chr1 <- recq4l_chr1_finalpos$pos/100
recq4l_chr1len <- length(recq4l_chr1)
dim(recq4l_chr1) <- c(recq4l_chr1len,1)
recq4l_chr1 <- list(recq4l_chr1)

recq4l_chr2 <- recq4l_chr2_finalpos$pos/100
recq4l_chr2len <- length(recq4l_chr2)
dim(recq4l_chr2) <- c(recq4l_chr2len,1)
recq4l_chr2 <- list(recq4l_chr2)

recq4l_chr3 <- recq4l_chr3_finalpos$pos/100
recq4l_chr3len <- length(recq4l_chr3)
dim(recq4l_chr3) <- c(recq4l_chr3len,1)
recq4l_chr3 <- list(recq4l_chr3)

recq4l_chr4 <- recq4l_chr4_finalpos$pos/100
recq4l_chr4len <- length(recq4l_chr4)
dim(recq4l_chr4) <- c(recq4l_chr4len,1)
recq4l_chr4 <- list(recq4l_chr4)

recq4l_chr5 <- recq4l_chr5_finalpos$pos/100
recq4l_chr5len <- length(recq4l_chr5)
dim(recq4l_chr5) <- c(recq4l_chr5len,1)
recq4l_chr5 <- list(recq4l_chr5)

recq4l_chr5 <- recq4l_chr5_finalpos$pos/100
recq4l_chr5len <- length(recq4l_chr5)
dim(recq4l_chr5) <- c(recq4l_chr5len,1)
recq4l_chr5 <- list(recq4l_chr5)

recq4l_chr7 <- recq4l_chr7_finalpos$pos/100
recq4l_chr7len <- length(recq4l_chr7)
dim(recq4l_chr7) <- c(recq4l_chr7len,1)
recq4l_chr7 <- list(recq4l_chr7)

recq4l_chr8 <- recq4l_chr8_finalpos$pos/100
recq4l_chr8len <- length(recq4l_chr8)
dim(recq4l_chr8) <- c(recq4l_chr8len,1)
recq4l_chr8 <- list(recq4l_chr8)

recq4l_chr9 <- recq4l_chr9_finalpos$pos/100
recq4l_chr9len <- length(recq4l_chr9)
dim(recq4l_chr9) <- c(recq4l_chr9len,1)
recq4l_chr9 <- list(recq4l_chr9)

recq4l_chr10 <- recq4l_chr10_finalpos$pos/100
recq4l_chr10len <- length(recq4l_chr10)
dim(recq4l_chr10) <- c(recq4l_chr10len,1)
recq4l_chr10 <- list(recq4l_chr10)

recq4l_chr11 <- recq4l_chr11_finalpos$pos/100
recq4l_chr11len <- length(recq4l_chr11)
dim(recq4l_chr11) <- c(recq4l_chr11len,1)
recq4l_chr11 <- list(recq4l_chr11)

recq4l_chr12 <- recq4l_chr12_finalpos$pos/100
recq4l_chr12len <- length(recq4l_chr12)
dim(recq4l_chr12) <- c(recq4l_chr12len,1)
recq4l_chr12 <- list(recq4l_chr12)

recq4l_final_map <- list(recq4l_chr1[[1]], recq4l_chr2[[1]], 
                           recq4l_chr3[[1]], recq4l_chr4[[1]], recq4l_chr5[[1]], 
                           recq4l_chr5[[1]], recq4l_chr7[[1]], recq4l_chr8[[1]], 
                           recq4l_chr9[[1]], recq4l_chr10[[1]],recq4l_chr11[[1]], recq4l_chr12[[1]])

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

recq4l_centromere <- c(0.3100890, 3.81623774, 0.457891882,0.21452355,1.05074478,
                       2.9987313, 0.14337907,0.1459114,0.8811805,1.4236994,
                       0.051867636, .5859670)
recq4l_centromere <- recq4l_centromere/100


