library(AlphaSimR)
library(Rcpp)
library(ggplot2)
library(dplyr)

set.seed(420)

japonica_snps <- read.table("japonica_SNPs.bed", header =FALSE)
colnames(japonica_snps) <- c("Chr#", "SNP Start", "SNP End")
ddm1_snps <- sample_n(japonica_snps, 4000)
ddm1_snps <- ddm1_snps[order(ddm1_snps$`Chr#`,ddm1_snps$`SNP Start`),]

#splitting ddm1 snps
ddm1_chr1_snp <- ddm1_snps[ which(ddm1_snps$`Chr#` == "Chr1"),]
ddm1_chr1_snp$rate <- NA
ddm1_chr1_snp$`SNP End` <- ddm1_chr1_snp$`SNP End` - min(ddm1_chr1_snp$`SNP Start`)
ddm1_chr1_snp$`SNP Start` <- ddm1_chr1_snp$`SNP Start`- min(ddm1_chr1_snp$`SNP Start`)

ddm1_chr2_snp <- ddm1_snps[ which(ddm1_snps$`Chr#` == "Chr2"),]
ddm1_chr2_snp$rate <- NA
ddm1_chr2_snp$`SNP End` <- ddm1_chr2_snp$`SNP End` - min(ddm1_chr2_snp$`SNP Start`)
ddm1_chr2_snp$`SNP Start` <- ddm1_chr2_snp$`SNP Start`- min(ddm1_chr2_snp$`SNP Start`)

ddm1_chr3_snp <- ddm1_snps[ which(ddm1_snps$`Chr#` == "Chr3"),]
ddm1_chr3_snp$rate <- NA
ddm1_chr3_snp$`SNP End` <- ddm1_chr3_snp$`SNP End` - min(ddm1_chr3_snp$`SNP Start`)
ddm1_chr3_snp$`SNP Start` <- ddm1_chr3_snp$`SNP Start`- min(ddm1_chr3_snp$`SNP Start`)

ddm1_chr4_snp <- ddm1_snps[ which(ddm1_snps$`Chr#` == "Chr4"),]
ddm1_chr4_snp$rate <- NA
ddm1_chr4_snp$`SNP End` <- ddm1_chr4_snp$`SNP End` - min(ddm1_chr4_snp$`SNP Start`)
ddm1_chr4_snp$`SNP Start` <- ddm1_chr4_snp$`SNP Start`- min(ddm1_chr4_snp$`SNP Start`)

ddm1_chr5_snp <- ddm1_snps[ which(ddm1_snps$`Chr#` == "Chr5"),]
ddm1_chr5_snp$rate <- NA
ddm1_chr5_snp$`SNP End` <- ddm1_chr5_snp$`SNP End` - min(ddm1_chr5_snp$`SNP Start`)
ddm1_chr5_snp$`SNP Start` <- ddm1_chr5_snp$`SNP Start`- min(ddm1_chr5_snp$`SNP Start`)

ddm1_chr6_snp <- ddm1_snps[ which(ddm1_snps$`Chr#` == "Chr6"),]
ddm1_chr6_snp$rate <- NA
ddm1_chr6_snp$`SNP End` <- ddm1_chr6_snp$`SNP End` - min(ddm1_chr6_snp$`SNP Start`)
ddm1_chr6_snp$`SNP Start` <- ddm1_chr6_snp$`SNP Start`- min(ddm1_chr6_snp$`SNP Start`)

ddm1_chr7_snp <- ddm1_snps[ which(ddm1_snps$`Chr#` == "Chr7"),]
ddm1_chr7_snp$rate <- NA
ddm1_chr7_snp$`SNP End` <- ddm1_chr7_snp$`SNP End` - min(ddm1_chr7_snp$`SNP Start`)
ddm1_chr7_snp$`SNP Start` <- ddm1_chr7_snp$`SNP Start`- min(ddm1_chr7_snp$`SNP Start`)

ddm1_chr8_snp <- ddm1_snps[ which(ddm1_snps$`Chr#` == "Chr8"),]
ddm1_chr8_snp$rate <- NA
ddm1_chr8_snp$`SNP End` <- ddm1_chr8_snp$`SNP End` - min(ddm1_chr8_snp$`SNP Start`)
ddm1_chr8_snp$`SNP Start` <- ddm1_chr8_snp$`SNP Start`- min(ddm1_chr8_snp$`SNP Start`)

ddm1_chr9_snp <- ddm1_snps[ which(ddm1_snps$`Chr#` == "Chr9"),]
ddm1_chr9_snp$rate <- NA
ddm1_chr9_snp$`SNP End` <- ddm1_chr9_snp$`SNP End` - min(ddm1_chr9_snp$`SNP Start`)
ddm1_chr9_snp$`SNP Start` <- ddm1_chr9_snp$`SNP Start`- min(ddm1_chr9_snp$`SNP Start`)

ddm1_chr10_snp <- ddm1_snps[ which(ddm1_snps$`Chr#` == "Chr10"),]
ddm1_chr10_snp$rate <- NA
ddm1_chr10_snp$`SNP End` <- ddm1_chr10_snp$`SNP End` - min(ddm1_chr10_snp$`SNP Start`)
ddm1_chr10_snp$`SNP Start` <- ddm1_chr10_snp$`SNP Start`- min(ddm1_chr10_snp$`SNP Start`)

ddm1_chr11_snp <- ddm1_snps[ which(ddm1_snps$`Chr#` == "Chr11"),]
ddm1_chr11_snp$rate <- NA
ddm1_chr11_snp$`SNP End` <- ddm1_chr11_snp$`SNP End` - min(ddm1_chr11_snp$`SNP Start`)
ddm1_chr11_snp$`SNP Start` <- ddm1_chr11_snp$`SNP Start`- min(ddm1_chr11_snp$`SNP Start`)

ddm1_chr12_snp <- ddm1_snps[ which(ddm1_snps$`Chr#` == "Chr12"),]
ddm1_chr12_snp$rate <- NA
ddm1_chr12_snp$`SNP End` <- ddm1_chr12_snp$`SNP End` - min(ddm1_chr12_snp$`SNP Start`)
ddm1_chr12_snp$`SNP Start` <- ddm1_chr12_snp$`SNP Start`- min(ddm1_chr12_snp$`SNP Start`)

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

#apply avg difference to telomeric regions (divide chromosome into fifths and apply avg diff to first and last fifth)
avg_diff <-2.346557
jap_chr1_CO$avg_rate <- avg_diff
jap_chr1_CO[8:23,6] <- 0

rownames(jap_chr2_CO)<-c(1:20)
jap_chr2_CO$avg_rate <- avg_diff
jap_chr2_CO[7:15,6] <- 0

rownames(jap_chr3_CO)<-c(1:13)
jap_chr3_CO$avg_rate <- avg_diff
jap_chr3_CO[6:12,6] <- 0

rownames(jap_chr4_CO)<-c(1:17)
jap_chr4_CO$avg_rate <- avg_diff
jap_chr4_CO[7:16,6] <- 0

rownames(jap_chr5_CO)<-c(1:15)
jap_chr5_CO$avg_rate <- avg_diff
jap_chr5_CO[5:13,6] <- 0

rownames(jap_chr6_CO)<-c(1:18)
jap_chr6_CO$avg_rate <- avg_diff
jap_chr6_CO[3:15,6] <- 0

rownames(jap_chr7_CO)<-c(1:20)
jap_chr7_CO$avg_rate <- avg_diff
jap_chr7_CO[5:19,6] <- 0

rownames(jap_chr8_CO)<-c(1:19)
jap_chr8_CO$avg_rate <- avg_diff
jap_chr8_CO[7:17,6] <- 0

rownames(jap_chr9_CO)<-c(1:17)
jap_chr9_CO$avg_rate <- avg_diff
jap_chr9_CO[5:14,6] <- 0

rownames(jap_chr10_CO)<-c(1:21)
jap_chr10_CO$avg_rate <- avg_diff
jap_chr10_CO[5:20,6] <- 0

rownames(jap_chr11_CO)<-c(1:20)
jap_chr11_CO$avg_rate <- avg_diff
jap_chr11_CO[5:16,6] <- 0

rownames(jap_chr12_CO)<-c(1:21)
jap_chr12_CO$avg_rate <- avg_diff
jap_chr12_CO[6:19,6] <- 0

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

#make intervals start at 0
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
new_rates <- function(avg_rate, old_rate){
  for(i in 1:nrow(old_rate)){
    for(k in 1:nrow(avg_rate)){
      if(isTRUE((old_rate$`CO Start`[i] >= (avg_rate$`CO Start`[k]*1000000)) && (old_rate$`CO End`[i] <= (avg_rate$`CO End`[k] *1000000)))){
        old_rate$rate[i] <- old_rate$rate[i] + (avg_rate$avg_rate[k]*old_rate$rate[i])
      }
    }
  }
  print(old_rate)
}
ddm1_chr1_CO_2<-new_rates(jap_chr1_CO, WTJap_chr1_CO)
ddm1_chr2_CO_2<-new_rates(jap_chr2_CO, WTJap_chr2_CO)
ddm1_chr3_CO_2<-new_rates(jap_chr3_CO, WTJap_chr3_CO)
ddm1_chr4_CO_2<-new_rates(jap_chr4_CO, WTJap_chr4_CO)
ddm1_chr5_CO_2<-new_rates(jap_chr5_CO, WTJap_chr5_CO)
ddm1_chr6_CO_2<-new_rates(jap_chr6_CO, WTJap_chr6_CO)
ddm1_chr7_CO_2<-new_rates(jap_chr7_CO, WTJap_chr7_CO)
ddm1_chr8_CO_2<-new_rates(jap_chr8_CO, WTJap_chr8_CO)
ddm1_chr9_CO_2<-new_rates(jap_chr9_CO, WTJap_chr9_CO)
ddm1_chr10_CO_2<-new_rates(jap_chr10_CO, WTJap_chr10_CO)
ddm1_chr11_CO_2<-new_rates(jap_chr11_CO, WTJap_chr11_CO)
ddm1_chr12_CO_2<-new_rates(jap_chr12_CO, WTJap_chr12_CO)

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
ddm1_chr1_CO_3 <- ddm1_chr1_CO_2
bins<-as.integer(nrow(ddm1_chr1_CO_2)/440)
ddm1_chr1_CO_3$rates<- rollapply(ddm1_chr1_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ddm1_chr1_CO_3<-fill_start(ddm1_chr1_CO_3)
ddm1_chr1_CO_3<- ddm1_chr1_CO_3 %>% drop_na(rates)

ddm1_chr2_CO_3 <- ddm1_chr2_CO_2
bins<-as.integer(nrow(ddm1_chr2_CO_2)/400)
ddm1_chr2_CO_3$rates<- rollapply(ddm1_chr2_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ddm1_chr2_CO_3<-fill_start(ddm1_chr2_CO_3)
ddm1_chr2_CO_3<- ddm1_chr2_CO_3 %>% drop_na(rates)

ddm1_chr3_CO_3 <- ddm1_chr3_CO_2
bins<-as.integer(nrow(ddm1_chr3_CO_2)/410)
ddm1_chr3_CO_3$rates<- rollapply(ddm1_chr3_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ddm1_chr3_CO_3<-fill_start(ddm1_chr3_CO_3)
ddm1_chr3_CO_3<- ddm1_chr3_CO_3 %>% drop_na(rates)

ddm1_chr4_CO_3 <- ddm1_chr4_CO_2
bins<-as.integer(nrow(ddm1_chr4_CO_2)/390)
ddm1_chr4_CO_3$rates<- rollapply(ddm1_chr4_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ddm1_chr4_CO_3<-fill_start(ddm1_chr4_CO_3)
ddm1_chr4_CO_3<- ddm1_chr4_CO_3 %>% drop_na(rates)

ddm1_chr5_CO_3 <- ddm1_chr5_CO_2
bins<-as.integer(nrow(ddm1_chr5_CO_2)/330)
ddm1_chr5_CO_3$rates<- rollapply(ddm1_chr5_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ddm1_chr5_CO_3<-fill_start(ddm1_chr5_CO_3)
ddm1_chr5_CO_3<- ddm1_chr5_CO_3 %>% drop_na(rates)

ddm1_chr6_CO_3 <- ddm1_chr6_CO_2
bins<-as.integer(nrow(ddm1_chr6_CO_2)/320)
ddm1_chr6_CO_3$rates<- rollapply(ddm1_chr6_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ddm1_chr6_CO_3<-fill_start(ddm1_chr6_CO_3)
ddm1_chr6_CO_3<- ddm1_chr6_CO_3 %>% drop_na(rates)

ddm1_chr7_CO_3 <- ddm1_chr7_CO_2
bins<-as.integer(nrow(ddm1_chr7_CO_2)/350)
ddm1_chr7_CO_3$rates<- rollapply(ddm1_chr7_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ddm1_chr7_CO_3<-fill_start(ddm1_chr7_CO_3)
ddm1_chr7_CO_3<- ddm1_chr7_CO_3 %>% drop_na(rates)

ddm1_chr8_CO_3 <- ddm1_chr8_CO_2
bins<-as.integer(nrow(ddm1_chr8_CO_2)/280)
ddm1_chr8_CO_3$rates<- rollapply(ddm1_chr8_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ddm1_chr8_CO_3<-fill_start(ddm1_chr8_CO_3)
ddm1_chr8_CO_3<- ddm1_chr8_CO_3 %>% drop_na(rates)

ddm1_chr9_CO_3 <- ddm1_chr9_CO_2
bins<-as.integer(nrow(ddm1_chr9_CO_2)/220)
ddm1_chr9_CO_3$rates<- rollapply(ddm1_chr9_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ddm1_chr9_CO_3<-fill_start(ddm1_chr9_CO_3)
ddm1_chr9_CO_3<- ddm1_chr9_CO_3 %>% drop_na(rates)

ddm1_chr10_CO_3 <- ddm1_chr10_CO_2
bins<-as.integer(nrow(ddm1_chr10_CO_2)/270)
ddm1_chr10_CO_3$rates<- rollapply(ddm1_chr10_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ddm1_chr10_CO_3<-fill_start(ddm1_chr10_CO_3)
ddm1_chr10_CO_3<- ddm1_chr10_CO_3 %>% drop_na(rates)

ddm1_chr11_CO_3 <- ddm1_chr11_CO_2
bins<-as.integer(nrow(ddm1_chr11_CO_2)/300)
ddm1_chr11_CO_3$rates<- rollapply(ddm1_chr11_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ddm1_chr11_CO_3<-fill_start(ddm1_chr11_CO_3)
ddm1_chr11_CO_3<- ddm1_chr11_CO_3 %>% drop_na(rates)

ddm1_chr12_CO_3 <- ddm1_chr12_CO_2
bins<-as.integer(nrow(ddm1_chr12_CO_2)/310)
ddm1_chr12_CO_3$rates<- rollapply(ddm1_chr12_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ddm1_chr12_CO_3<-fill_start(ddm1_chr12_CO_3)
ddm1_chr12_CO_3<- ddm1_chr12_CO_3 %>% drop_na(rates)

## assigning frequency to SNPs based on recombination frequency in each bin
snp_rate <- function(chr_rate, chr_snp){
  for(i in 1:nrow(chr_snp)){
    for(k in 1:nrow(chr_rate)){
      if(isTRUE((chr_snp$`SNP Start`[i] >= chr_rate$`CO Start`[k]) && (chr_snp$`SNP End`[i] <= chr_rate$`CO End`[k]))){
        chr_snp$rate[i] <- chr_rate$rate[k]
      }
    }
  }
  print(chr_snp)
}


#using function,  get cM/Mb for final genetic position - assign rates
ddm1_chr1_snp2 <- snp_rate(ddm1_chr1_CO_3, ddm1_chr1_snp)
ddm1_chr2_snp2 <- snp_rate(ddm1_chr2_CO_3, ddm1_chr2_snp)
ddm1_chr3_snp2 <- snp_rate(ddm1_chr3_CO_3, ddm1_chr3_snp)
ddm1_chr4_snp2 <- snp_rate(ddm1_chr4_CO_3, ddm1_chr4_snp)
ddm1_chr5_snp2 <- snp_rate(ddm1_chr5_CO_3, ddm1_chr5_snp)
ddm1_chr6_snp2 <- snp_rate(ddm1_chr6_CO_3, ddm1_chr6_snp)
ddm1_chr7_snp2 <- snp_rate(ddm1_chr7_CO_3, ddm1_chr7_snp)
ddm1_chr8_snp2 <- snp_rate(ddm1_chr8_CO_3, ddm1_chr8_snp)
ddm1_chr9_snp2 <- snp_rate(ddm1_chr9_CO_3, ddm1_chr9_snp)
ddm1_chr10_snp2 <- snp_rate(ddm1_chr10_CO_3, ddm1_chr10_snp)
ddm1_chr11_snp2 <- snp_rate(ddm1_chr11_CO_3, ddm1_chr11_snp)
ddm1_chr12_snp2 <- snp_rate(ddm1_chr12_CO_3, ddm1_chr12_snp)

#converted SNP start to Mb
ddm1_chr1_snp2$`SNP Start` <- ddm1_chr1_snp2$`SNP Start`/1000000
ddm1_chr2_snp2$`SNP Start` <- ddm1_chr2_snp2$`SNP Start`/1000000
ddm1_chr3_snp2$`SNP Start` <- ddm1_chr3_snp2$`SNP Start`/1000000
ddm1_chr4_snp2$`SNP Start` <- ddm1_chr4_snp2$`SNP Start`/1000000
ddm1_chr5_snp2$`SNP Start` <- ddm1_chr5_snp2$`SNP Start`/1000000
ddm1_chr6_snp2$`SNP Start` <- ddm1_chr6_snp2$`SNP Start`/1000000
ddm1_chr7_snp2$`SNP Start` <- ddm1_chr7_snp2$`SNP Start`/1000000
ddm1_chr8_snp2$`SNP Start` <- ddm1_chr8_snp2$`SNP Start`/1000000
ddm1_chr9_snp2$`SNP Start` <- ddm1_chr9_snp2$`SNP Start`/1000000
ddm1_chr10_snp2$`SNP Start` <- ddm1_chr10_snp2$`SNP Start`/1000000
ddm1_chr11_snp2$`SNP Start` <- ddm1_chr11_snp2$`SNP Start`/1000000
ddm1_chr12_snp2$`SNP Start` <- ddm1_chr12_snp2$`SNP Start`/1000000

ddm1_chr1_snp2$`SNP End` <- ddm1_chr1_snp2$`SNP End`/1000000
ddm1_chr2_snp2$`SNP End` <- ddm1_chr2_snp2$`SNP End`/1000000
ddm1_chr3_snp2$`SNP End` <- ddm1_chr3_snp2$`SNP End`/1000000
ddm1_chr4_snp2$`SNP End` <- ddm1_chr4_snp2$`SNP End`/1000000
ddm1_chr5_snp2$`SNP End` <- ddm1_chr5_snp2$`SNP End`/1000000
ddm1_chr6_snp2$`SNP End` <- ddm1_chr6_snp2$`SNP End`/1000000
ddm1_chr7_snp2$`SNP End` <- ddm1_chr7_snp2$`SNP End`/1000000
ddm1_chr8_snp2$`SNP End` <- ddm1_chr8_snp2$`SNP End`/1000000
ddm1_chr9_snp2$`SNP End` <- ddm1_chr9_snp2$`SNP End`/1000000
ddm1_chr10_snp2$`SNP End` <- ddm1_chr10_snp2$`SNP End`/1000000
ddm1_chr11_snp2$`SNP End` <- ddm1_chr11_snp2$`SNP End`/1000000
ddm1_chr12_snp2$`SNP End` <- ddm1_chr12_snp2$`SNP End`/1000000

#omit empty col
ddm1_chr1_snp2<-na.omit(ddm1_chr1_snp2)
ddm1_chr2_snp2<-na.omit(ddm1_chr2_snp2)
ddm1_chr3_snp2<-na.omit(ddm1_chr3_snp2)
ddm1_chr4_snp2<-na.omit(ddm1_chr4_snp2)
ddm1_chr5_snp2<-na.omit(ddm1_chr5_snp2)
ddm1_chr6_snp2<-na.omit(ddm1_chr6_snp2)
ddm1_chr7_snp2<-na.omit(ddm1_chr7_snp2)
ddm1_chr8_snp2<-na.omit(ddm1_chr8_snp2)
ddm1_chr9_snp2<-na.omit(ddm1_chr9_snp2)
ddm1_chr10_snp2<-na.omit(ddm1_chr10_snp2)
ddm1_chr11_snp2<-na.omit(ddm1_chr11_snp2)
ddm1_chr12_snp2<-na.omit(ddm1_chr12_snp2)

#gen maps
ddm1_chr1_spl <- smooth.spline(ddm1_chr1_snp2$rate, spar = 1.2)
ddm1_chr1_snp2$pos <- (ddm1_chr1_snp2$`SNP Start`*ddm1_chr1_spl$y)
#ddm1_chr1_snp2$pos <- (ddm1_chr1_snp2$`SNP Start`*ddm1_chr1_snp2$rate)
plot(ddm1_chr1_snp2$`SNP Start`, ddm1_chr1_snp2$pos)
ggplot(ddm1_chr1_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
plot(ddm1_chr1_snp2$`SNP Start`, ddm1_chr1_snp2$pos/ddm1_chr1_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica ddm1 Chromosome 1 Recombination Distribution")
ddm1_chr1_finalpos <- ddm1_chr1_snp2[order(ddm1_chr1_snp2$pos),]
is.unsorted(ddm1_chr1_finalpos$pos)
#ddm1_chr1_spl <- smooth.spline(ddm1_chr1_finalpos$pos, spar = .5)
plot(ddm1_chr1_snp2$`SNP Start`, ddm1_chr1_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 1 Genetic Map")
plot(ddm1_chr1_finalpos$`SNP Start`, ddm1_chr1_finalpos$pos)

ddm1_chr2_spl <- smooth.spline(ddm1_chr2_snp2$rate, spar = 1.2)
ddm1_chr2_snp2$pos <- (ddm1_chr2_snp2$`SNP Start`*ddm1_chr2_spl$y)
plot(ddm1_chr2_snp2$`SNP Start`, ddm1_chr2_snp2$pos)
plot(ddm1_chr2_snp2$`SNP Start`, ddm1_chr2_snp2$pos/ddm1_chr2_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica ddm1 Chromosome 2 Recombination Distribution")
ddm1_chr2_finalpos <- ddm1_chr2_snp2[order(ddm1_chr2_snp2$pos),]
is.unsorted(ddm1_chr2_finalpos$pos)
plot(ddm1_chr2_snp2$`SNP Start`, ddm1_chr2_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 2 Genetic Map")

ddm1_chr3_spl <- smooth.spline(ddm1_chr3_snp2$rate, spar = .7)
ddm1_chr3_snp2$pos <- (ddm1_chr3_snp2$`SNP Start`*ddm1_chr3_spl$y)
plot(ddm1_chr3_snp2$`SNP Start`, ddm1_chr3_snp2$pos)
plot(ddm1_chr3_snp2$`SNP Start`, ddm1_chr3_snp2$pos/ddm1_chr3_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica ddm1 Chromosome 3 Recombination Distribution")
ddm1_chr3_finalpos <- ddm1_chr3_snp2[order(ddm1_chr3_snp2$pos),]
is.unsorted(ddm1_chr3_finalpos$pos)
plot(ddm1_chr3_snp2$`SNP Start`, ddm1_chr3_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 3 Genetic Map")

ddm1_chr4_spl <- smooth.spline(ddm1_chr4_snp2$rate, spar =.8)
ddm1_chr4_snp2$pos <- (ddm1_chr4_snp2$`SNP Start`*ddm1_chr4_spl$y)
plot(ddm1_chr4_snp2$`SNP Start`, ddm1_chr4_snp2$pos)
plot(ddm1_chr4_snp2$`SNP Start`, ddm1_chr4_snp2$pos/ddm1_chr4_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica ddm1 Chromosome 4 Recombination Distribution")
ddm1_chr4_finalpos <- ddm1_chr4_snp2[order(ddm1_chr4_snp2$pos),]
is.unsorted(ddm1_chr4_finalpos$pos)
plot(ddm1_chr4_snp2$`SNP Start`, ddm1_chr4_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 4 Genetic Map")

ddm1_chr5_spl <- smooth.spline(ddm1_chr5_snp2$rate, spar =.9)
ddm1_chr5_snp2$pos <- (ddm1_chr5_snp2$`SNP Start`*ddm1_chr5_spl$y)
plot(ddm1_chr5_snp2$`SNP Start`, ddm1_chr5_snp2$pos)
plot(ddm1_chr5_snp2$`SNP Start`, ddm1_chr5_snp2$pos/ddm1_chr5_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica ddm1 Chromosome 5 Recombination Distribution")
ddm1_chr5_finalpos <- ddm1_chr5_snp2[order(ddm1_chr5_snp2$pos),]
is.unsorted(ddm1_chr5_finalpos$pos)
ddm1_chr5_finalpos$pos <- ddm1_chr5_finalpos$pos + abs(min(ddm1_chr5_finalpos$pos))
plot(ddm1_chr5_snp2$`SNP Start`, ddm1_chr5_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 5 Genetic Map")

ddm1_chr6_spl <- smooth.spline(ddm1_chr6_snp2$rate, spar = .8)
ddm1_chr6_snp2$pos <- (ddm1_chr6_snp2$`SNP Start`*ddm1_chr6_spl$y)
plot(ddm1_chr6_snp2$`SNP Start`, ddm1_chr6_snp2$pos)
plot(ddm1_chr6_snp2$`SNP Start`, ddm1_chr6_snp2$pos/ddm1_chr6_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica ddm1 Chromosome 6 Recombination Distribution")
ddm1_chr6_finalpos <- ddm1_chr6_snp2[order(ddm1_chr6_snp2$pos),]
is.unsorted(ddm1_chr6_finalpos$pos)
plot(ddm1_chr6_snp2$`SNP Start`, ddm1_chr6_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 6 Genetic Map")

ddm1_chr7_spl <- smooth.spline(ddm1_chr7_snp2$rate, spar = 1.2)
ddm1_chr7_snp2$pos <- (ddm1_chr7_snp2$`SNP Start`*ddm1_chr7_spl$y)
plot(ddm1_chr7_snp2$`SNP Start`, ddm1_chr7_snp2$pos)
plot(ddm1_chr7_snp2$`SNP Start`, ddm1_chr7_snp2$pos/ddm1_chr7_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica ddm1 Chromosome 7 Recombination Distribution")
ddm1_chr7_finalpos <- ddm1_chr7_snp2[order(ddm1_chr7_snp2$pos),]
is.unsorted(ddm1_chr7_finalpos$pos)
plot(ddm1_chr7_snp2$`SNP Start`, ddm1_chr7_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 7 Genetic Map")

ddm1_chr8_spl <- smooth.spline(ddm1_chr8_snp2$rate, spar = .8)
ddm1_chr8_snp2$pos <- (ddm1_chr8_snp2$`SNP Start`*ddm1_chr8_spl$y)
plot(ddm1_chr8_snp2$`SNP Start`, ddm1_chr8_snp2$pos)
plot(ddm1_chr8_snp2$`SNP Start`, ddm1_chr8_snp2$pos/ddm1_chr8_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica ddm1 Chromosome 8 Recombination Distribution")
ddm1_chr8_finalpos <- ddm1_chr8_snp2[order(ddm1_chr8_snp2$pos),]
is.unsorted(ddm1_chr8_finalpos$pos)
plot(ddm1_chr8_snp2$`SNP Start`, ddm1_chr8_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 8 Genetic Map")

ddm1_chr9_spl <- smooth.spline(ddm1_chr9_snp2$rate, spar = .9)
ddm1_chr9_snp2$pos <- (ddm1_chr9_snp2$`SNP Start`*ddm1_chr9_spl$y)
plot(ddm1_chr9_snp2$`SNP Start`, ddm1_chr9_snp2$pos)
plot(ddm1_chr9_snp2$`SNP Start`, ddm1_chr9_snp2$pos/ddm1_chr9_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica ddm1 Chromosome 9 Recombination Distribution")
ddm1_chr9_finalpos <- ddm1_chr9_snp2[order(ddm1_chr9_snp2$pos),]
is.unsorted(ddm1_chr9_finalpos$pos)
plot(ddm1_chr9_snp2$`SNP Start`, ddm1_chr9_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 9 Genetic Map")

ddm1_chr10_spl <- smooth.spline(ddm1_chr10_snp2$rate, spar =.6)
ddm1_chr10_snp2$pos <- (ddm1_chr10_snp2$`SNP Start`*ddm1_chr10_spl$y)
plot(ddm1_chr10_snp2$`SNP Start`, ddm1_chr10_snp2$pos)
plot(ddm1_chr10_snp2$`SNP Start`, ddm1_chr10_snp2$pos/ddm1_chr10_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica ddm1 Chromosome 10 Recombination Distribution")
ddm1_chr10_finalpos <- ddm1_chr10_snp2[order(ddm1_chr10_snp2$pos),]
is.unsorted(ddm1_chr10_finalpos$pos)
plot(ddm1_chr10_snp2$`SNP Start`, ddm1_chr10_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 10 Genetic Map")

ddm1_chr11_spl <- smooth.spline(ddm1_chr11_snp2$rate, spar = .45)
ddm1_chr11_snp2$pos <- (ddm1_chr11_snp2$`SNP Start`*ddm1_chr11_spl$y)
plot(ddm1_chr11_snp2$`SNP Start`, ddm1_chr11_snp2$pos)
plot(ddm1_chr11_snp2$`SNP Start`, ddm1_chr11_snp2$pos/ddm1_chr11_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica ddm1 Chromosome 11 Recombination Distribution")
ddm1_chr11_finalpos <- ddm1_chr11_snp2[order(ddm1_chr11_snp2$pos),]
is.unsorted(ddm1_chr11_finalpos$pos)
plot(ddm1_chr11_snp2$`SNP Start`, ddm1_chr11_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 11 Genetic Map")

ddm1_chr12_spl <- smooth.spline(ddm1_chr12_snp2$rate, spar = .415)
ddm1_chr12_snp2$pos <- (ddm1_chr12_snp2$`SNP Start`*ddm1_chr12_spl$y)
plot(ddm1_chr12_snp2$`SNP Start`, ddm1_chr12_snp2$pos)
plot(ddm1_chr12_snp2$`SNP Start`, ddm1_chr12_snp2$pos/ddm1_chr12_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica ddm1 Chromosome 12 Recombination Distribution")
ddm1_chr12_finalpos <- ddm1_chr12_snp2[order(ddm1_chr12_snp2$pos),]
is.unsorted(ddm1_chr12_finalpos$pos)
plot(ddm1_chr12_snp2$`SNP Start`, ddm1_chr12_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 12 Genetic Map")

# ddm1_chr2_snp2$pos <- (ddm1_chr2_snp2$`SNP Start`*ddm1_chr2_snp2$rate)
# plot(ddm1_chr2_snp2$`SNP Start`, ddm1_chr2_snp2$pos)
# plot(ddm1_chr2_snp2$`SNP Start`, ddm1_chr2_snp2$pos/ddm1_chr2_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
#      ylab = "Recombination rate (cM/Mb)", main = "Japonica ddm1 Chromosome 2 Recombination Distribution")
# ddm1_chr2_finalpos <- ddm1_chr2_snp2[order(ddm1_chr2_snp2$pos),]
# is.unsorted(ddm1_chr2_finalpos$pos)
# ddm1_chr2_spl <- smooth.spline(ddm1_chr2_finalpos$pos, spar = .5)
# plot(ddm1_chr2_snp2$`SNP Start`, ddm1_chr2_spl$y, type = "l", xlab = "Physical Positions (Mb)",
#      ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 2 Genetic Map")
# 
# ddm1_chr3_spl <- smooth.spline(ddm1_chr3_snp2$rate, spar = .7)
# ddm1_chr3_snp2$pos <- (ddm1_chr3_snp2$`SNP Start`*ddm1_chr3_spl$y)
# #ddm1_chr3_snp2$pos <- (ddm1_chr3_snp2$`SNP Start`*ddm1_chr3_snp2$rate)
# plot(ddm1_chr3_snp2$`SNP Start`, ddm1_chr3_snp2$pos)
# plot(ddm1_chr3_snp2$`SNP Start`, ddm1_chr3_snp2$pos/ddm1_chr3_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
#      ylab = "Recombination rate (cM/Mb)", main = "Japonica ddm1 Chromosome 3 Recombination Distribution")
# ddm1_chr3_finalpos <- ddm1_chr3_snp2[order(ddm1_chr3_snp2$pos),]
# is.unsorted(ddm1_chr3_finalpos$pos)
# #ddm1_chr3_spl <- smooth.spline(ddm1_chr3_finalpos$pos, spar = .5)
# plot(ddm1_chr3_snp2$`SNP Start`, ddm1_chr3_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
#      ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 3 Genetic Map")
# 
# #ddm1_chr4_spl <- smooth.spline(ddm1_chr4_snp2$rate, spar = .8)
# #ddm1_chr4_snp2$pos <- (ddm1_chr4_snp2$`SNP Start`*ddm1_chr4_spl$y)
# ddm1_chr4_snp2$pos <- (ddm1_chr4_snp2$`SNP Start`*ddm1_chr4_snp2$rate)
# plot(ddm1_chr4_snp2$`SNP Start`, ddm1_chr4_snp2$pos)
# plot(ddm1_chr4_snp2$`SNP Start`, ddm1_chr4_snp2$pos/ddm1_chr4_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
#      ylab = "Recombination rate (cM/Mb)", main = "Japonica ddm1 Chromosome 4 Recombination Distribution")
# ddm1_chr4_finalpos <- ddm1_chr4_snp2[order(ddm1_chr4_snp2$pos),]
# is.unsorted(ddm1_chr4_finalpos$pos)
# ddm1_chr4_spl <- smooth.spline(ddm1_chr4_finalpos$pos, spar = .5)
# plot(ddm1_chr4_snp2$`SNP Start`, ddm1_chr4_spl$y, type = "l", xlab = "Physical Positions (Mb)",
#      ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 4 Genetic Map")
# 
# #ddm1_chr5_spl <- smooth.spline(ddm1_chr5_snp2$rate, spar =.66)
# #ddm1_chr5_snp2$pos <- (ddm1_chr5_snp2$`SNP Start`*ddm1_chr5_spl$y)
# ddm1_chr5_snp2$pos <- (ddm1_chr5_snp2$`SNP Start`*ddm1_chr5_snp2$rate)
# plot(ddm1_chr5_snp2$`SNP Start`, ddm1_chr5_snp2$pos)
# plot(ddm1_chr5_snp2$`SNP Start`, ddm1_chr5_snp2$pos/ddm1_chr5_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
#      ylab = "Recombination rate (cM/Mb)", main = "Japonica ddm1 Chromosome 5 Recombination Distribution")
# ddm1_chr5_finalpos <- ddm1_chr5_snp2[order(ddm1_chr5_snp2$pos),]
# is.unsorted(ddm1_chr5_finalpos$pos)
# ddm1_chr5_spl <- smooth.spline(ddm1_chr5_finalpos$pos, spar = .5)
# plot(ddm1_chr5_snp2$`SNP Start`, ddm1_chr5_spl$y, type = "l", xlab = "Physical Positions (Mb)",
#      ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 5 Genetic Map")
# 
# ddm1_chr6_spl <- smooth.spline(ddm1_chr6_snp2$rate, spar = 1.2)
# ddm1_chr6_snp2$pos <- (ddm1_chr6_snp2$`SNP Start`*ddm1_chr6_spl$y)
# #ddm1_chr6_snp2$pos <- (ddm1_chr6_snp2$`SNP Start`*ddm1_chr6_snp2$rate)
# plot(ddm1_chr6_snp2$`SNP Start`, ddm1_chr6_snp2$pos)
# plot(ddm1_chr6_snp2$`SNP Start`, ddm1_chr6_snp2$pos/ddm1_chr6_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
#      ylab = "Recombination rate (cM/Mb)", main = "Japonica ddm1 Chromosome 6 Recombination Distribution")
# ddm1_chr6_finalpos <- ddm1_chr6_snp2[order(ddm1_chr6_snp2$pos),]
# is.unsorted(ddm1_chr6_finalpos$pos)
# #ddm1_chr6_spl <- smooth.spline(ddm1_chr6_finalpos$pos, spar = 1.3)
# plot(ddm1_chr6_snp2$`SNP Start`, ddm1_chr6_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
#      ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 6 Genetic Map")
# 
# #ddm1_chr7_spl <- smooth.spline(ddm1_chr7_snp2$rate, spar = 0.7)
# #ddm1_chr7_snp2$pos <- (ddm1_chr7_snp2$`SNP Start`*ddm1_chr7_spl$y)
# ddm1_chr7_snp2$pos <- (ddm1_chr7_snp2$`SNP Start`*ddm1_chr7_snp2$rate)
# plot(ddm1_chr7_snp2$`SNP Start`, ddm1_chr7_snp2$pos)
# plot(ddm1_chr7_snp2$`SNP Start`, ddm1_chr7_snp2$pos/ddm1_chr7_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
#      ylab = "Recombination rate (cM/Mb)", main = "Japonica ddm1 Chromosome 7 Recombination Distribution")
# ddm1_chr7_finalpos <- ddm1_chr7_snp2[order(ddm1_chr7_snp2$pos),]
# is.unsorted(ddm1_chr7_finalpos$pos)
# ddm1_chr7_spl <- smooth.spline(ddm1_chr7_finalpos$pos, spar = .5)
# plot(ddm1_chr7_snp2$`SNP Start`, ddm1_chr7_spl$y, type = "l", xlab = "Physical Positions (Mb)",
#      ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 7 Genetic Map")
# 
# #ddm1_chr8_spl <- smooth.spline(ddm1_chr8_snp2$rate, spar = 1)
# #ddm1_chr8_snp2$pos <- (ddm1_chr8_snp2$`SNP Start`*ddm1_chr8_spl$y)
# ddm1_chr8_snp2$pos <- (ddm1_chr8_snp2$`SNP Start`*ddm1_chr8_snp2$rate)
# plot(ddm1_chr8_snp2$`SNP Start`, ddm1_chr8_snp2$pos)
# plot(ddm1_chr8_snp2$`SNP Start`, ddm1_chr8_snp2$pos/ddm1_chr8_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
#      ylab = "Recombination rate (cM/Mb)", main = "Japonica ddm1 Chromosome 8 Recombination Distribution")
# ddm1_chr8_finalpos <- ddm1_chr8_snp2[order(ddm1_chr8_snp2$pos),]
# is.unsorted(ddm1_chr8_finalpos$pos)
# ddm1_chr8_spl <- smooth.spline(ddm1_chr8_finalpos$pos, spar = .5)
# plot(ddm1_chr8_snp2$`SNP Start`, ddm1_chr8_spl$y, type = "l", xlab = "Physical Positions (Mb)",
#      ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 8 Genetic Map")
# 
# #ddm1_chr9_spl <- smooth.spline(ddm1_chr9_snp2$rate, spar = .9)
# #ddm1_chr9_snp2$pos <- (ddm1_chr9_snp2$`SNP Start`*ddm1_chr9_spl$y)
# ddm1_chr9_snp2$pos <- (ddm1_chr9_snp2$`SNP Start`*ddm1_chr9_snp2$rate)
# plot(ddm1_chr9_snp2$`SNP Start`, ddm1_chr9_snp2$pos)
# plot(ddm1_chr9_snp2$`SNP Start`, ddm1_chr9_snp2$pos/ddm1_chr9_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
#      ylab = "Recombination rate (cM/Mb)", main = "Japonica ddm1 Chromosome 9 Recombination Distribution")
# ddm1_chr9_finalpos <- ddm1_chr9_snp2[order(ddm1_chr9_snp2$pos),]
# is.unsorted(ddm1_chr9_finalpos$pos)
# ddm1_chr9_spl <- smooth.spline(ddm1_chr9_finalpos$pos, spar = .5)
# plot(ddm1_chr9_snp2$`SNP Start`, ddm1_chr9_spl$y, type = "l", xlab = "Physical Positions (Mb)",
#      ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 9 Genetic Map")
# 
# #ddm1_chr10_spl <- smooth.spline(ddm1_chr10_snp2$rate, spar =1)
# #ddm1_chr10_snp2$pos <- (ddm1_chr10_snp2$`SNP Start`*ddm1_chr10_spl$y)
# ddm1_chr10_snp2$pos <- (ddm1_chr10_snp2$`SNP Start`*ddm1_chr10_snp2$rate)
# plot(ddm1_chr10_snp2$`SNP Start`, ddm1_chr10_snp2$pos)
# plot(ddm1_chr10_snp2$`SNP Start`, ddm1_chr10_snp2$pos/ddm1_chr10_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
#      ylab = "Recombination rate (cM/Mb)", main = "Japonica ddm1 Chromosome 10 Recombination Distribution")
# ddm1_chr10_finalpos <- ddm1_chr10_snp2[order(ddm1_chr10_snp2$pos),]
# is.unsorted(ddm1_chr10_finalpos$pos)
# ddm1_chr10_spl <- smooth.spline(ddm1_chr10_finalpos$pos, spar = .5)
# plot(ddm1_chr10_snp2$`SNP Start`, ddm1_chr10_spl$y, type = "l", xlab = "Physical Positions (Mb)",
#      ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 10 Genetic Map")
# 
# #ddm1_chr11_spl <- smooth.spline(ddm1_chr11_snp2$rate, spar = 1)
# #ddm1_chr11_snp2$pos <- (ddm1_chr11_snp2$`SNP Start`*ddm1_chr11_spl$y)
# ddm1_chr11_snp2$pos <- (ddm1_chr11_snp2$`SNP Start`*ddm1_chr11_snp2$rate)
# plot(ddm1_chr11_snp2$`SNP Start`, ddm1_chr11_snp2$pos)
# plot(ddm1_chr11_snp2$`SNP Start`, ddm1_chr11_snp2$pos/ddm1_chr11_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
#      ylab = "Recombination rate (cM/Mb)", main = "Japonica ddm1 Chromosome 11 Recombination Distribution")
# ddm1_chr11_finalpos <- ddm1_chr11_snp2[order(ddm1_chr11_snp2$pos),]
# is.unsorted(ddm1_chr11_finalpos$pos)
# ddm1_chr11_spl <- smooth.spline(ddm1_chr11_finalpos$pos, spar = .5)
# plot(ddm1_chr11_snp2$`SNP Start`, ddm1_chr11_spl$y, type = "l", xlab = "Physical Positions (Mb)",
#      ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 11 Genetic Map")
# 
# #ddm1_chr12_spl <- smooth.spline(ddm1_chr12_snp2$rate, spar = .99)
# #ddm1_chr12_snp2$pos <- (ddm1_chr12_snp2$`SNP Start`*ddm1_chr12_spl$y)
# ddm1_chr12_snp2$pos <- (ddm1_chr12_snp2$`SNP Start`*ddm1_chr12_snp2$rate)
# plot(ddm1_chr12_snp2$`SNP Start`, ddm1_chr12_snp2$pos)
# plot(ddm1_chr12_snp2$`SNP Start`, ddm1_chr12_snp2$pos/ddm1_chr12_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
#      ylab = "Recombination rate (cM/Mb)", main = "Japonica ddm1 Chromosome 12 Recombination Distribution")
# ddm1_chr12_finalpos <- ddm1_chr12_snp2[order(ddm1_chr12_snp2$pos),]
# is.unsorted(ddm1_chr12_finalpos$pos)
# ddm1_chr12_spl <- smooth.spline(ddm1_chr12_finalpos$pos, spar = .5)
# plot(ddm1_chr12_snp2$`SNP Start`, ddm1_chr12_spl$y, type = "l", xlab = "Physical Positions (Mb)",
#      ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 12 Genetic Map")

#Final genetic map
ddm1_chr1 <- ddm1_chr1_finalpos$pos/100
ddm1_chr1len <- length(ddm1_chr1)
dim(ddm1_chr1) <- c(ddm1_chr1len,1)
ddm1_chr1 <- list(ddm1_chr1)

ddm1_chr2 <- ddm1_chr2_finalpos$pos/100
ddm1_chr2len <- length(ddm1_chr2)
dim(ddm1_chr2) <- c(ddm1_chr2len,1)
ddm1_chr2 <- list(ddm1_chr2)

ddm1_chr3 <- ddm1_chr3_finalpos$pos/100
ddm1_chr3len <- length(ddm1_chr3)
dim(ddm1_chr3) <- c(ddm1_chr3len,1)
ddm1_chr3 <- list(ddm1_chr3)

ddm1_chr4 <- ddm1_chr4_finalpos$pos/100
ddm1_chr4len <- length(ddm1_chr4)
dim(ddm1_chr4) <- c(ddm1_chr4len,1)
ddm1_chr4 <- list(ddm1_chr4)

ddm1_chr5 <- ddm1_chr5_finalpos$pos/100
ddm1_chr5len <- length(ddm1_chr5)
dim(ddm1_chr5) <- c(ddm1_chr5len,1)
ddm1_chr5 <- list(ddm1_chr5)

ddm1_chr5 <- ddm1_chr5_finalpos$pos/100
ddm1_chr5len <- length(ddm1_chr5)
dim(ddm1_chr5) <- c(ddm1_chr5len,1)
ddm1_chr5 <- list(ddm1_chr5)

ddm1_chr7 <- ddm1_chr7_finalpos$pos/100
ddm1_chr7len <- length(ddm1_chr7)
dim(ddm1_chr7) <- c(ddm1_chr7len,1)
ddm1_chr7 <- list(ddm1_chr7)

ddm1_chr8 <- ddm1_chr8_finalpos$pos/100
ddm1_chr8len <- length(ddm1_chr8)
dim(ddm1_chr8) <- c(ddm1_chr8len,1)
ddm1_chr8 <- list(ddm1_chr8)

ddm1_chr9 <- ddm1_chr9_finalpos$pos/100
ddm1_chr9len <- length(ddm1_chr9)
dim(ddm1_chr9) <- c(ddm1_chr9len,1)
ddm1_chr9 <- list(ddm1_chr9)

ddm1_chr10 <- ddm1_chr10_finalpos$pos/100
ddm1_chr10len <- length(ddm1_chr10)
dim(ddm1_chr10) <- c(ddm1_chr10len,1)
ddm1_chr10 <- list(ddm1_chr10)

ddm1_chr11 <- ddm1_chr11_finalpos$pos/100
ddm1_chr11len <- length(ddm1_chr11)
dim(ddm1_chr11) <- c(ddm1_chr11len,1)
ddm1_chr11 <- list(ddm1_chr11)

ddm1_chr12 <- ddm1_chr12_finalpos$pos/100
ddm1_chr12len <- length(ddm1_chr12)
dim(ddm1_chr12) <- c(ddm1_chr12len,1)
ddm1_chr12 <- list(ddm1_chr12)

ddm1_final_map <- list(ddm1_chr1[[1]], ddm1_chr2[[1]], 
                        ddm1_chr3[[1]], ddm1_chr4[[1]], ddm1_chr5[[1]], 
                        ddm1_chr5[[1]], ddm1_chr7[[1]], ddm1_chr8[[1]], 
                        ddm1_chr9[[1]], ddm1_chr10[[1]],ddm1_chr11[[1]], ddm1_chr12[[1]])

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
c1 <-find_centromere(16.7,ddm1_chr1_finalpos)
c2 <-find_centromere(13.6,ddm1_chr2_finalpos)
c3 <-find_centromere(19.4,ddm1_chr3_finalpos)
c4 <-find_centromere(9.7,ddm1_chr4_finalpos)
c5 <-find_centromere(12.4,ddm1_chr5_finalpos)
c6 <-find_centromere(15.3,ddm1_chr6_finalpos)
c7 <-find_centromere(12.1,ddm1_chr7_finalpos)
c8 <-find_centromere(12.9,ddm1_chr8_finalpos)
c9 <-find_centromere(2.8,ddm1_chr9_finalpos)
c10 <-find_centromere(8.2,ddm1_chr10_finalpos)
c11 <-find_centromere(12,ddm1_chr11_finalpos)
c12 <-find_centromere(11.9,ddm1_chr12_finalpos)

ddm1_centromere <- c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12)
ddm1_centromere <- ddm1_centromere/100



