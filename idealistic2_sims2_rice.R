library(AlphaSimR)
library(Rcpp)
library(ggplot2)
library(dplyr)

set.seed(420)

japonica_snps <- read.table("japonica_SNPs.bed", header =FALSE)
colnames(japonica_snps) <- c("Chr#", "SNP Start", "SNP End")
ideal2_snps <- sample_n(japonica_snps, 2000)
ideal2_snps <- ideal2_snps[order(ideal2_snps$`Chr#`,ideal2_snps$`SNP Start`),]

#splitting ideal2 snps
ideal2_chr1_snp <- ideal2_snps[ which(ideal2_snps$`Chr#` == "Chr1"),]
ideal2_chr1_snp$rate <- NA
ideal2_chr1_snp$`SNP End` <- ideal2_chr1_snp$`SNP End` - min(ideal2_chr1_snp$`SNP Start`)
ideal2_chr1_snp$`SNP Start` <- ideal2_chr1_snp$`SNP Start`- min(ideal2_chr1_snp$`SNP Start`)

ideal2_chr2_snp <- ideal2_snps[ which(ideal2_snps$`Chr#` == "Chr2"),]
ideal2_chr2_snp$rate <- NA
ideal2_chr2_snp$`SNP End` <- ideal2_chr2_snp$`SNP End` - min(ideal2_chr2_snp$`SNP Start`)
ideal2_chr2_snp$`SNP Start` <- ideal2_chr2_snp$`SNP Start`- min(ideal2_chr2_snp$`SNP Start`)

ideal2_chr3_snp <- ideal2_snps[ which(ideal2_snps$`Chr#` == "Chr3"),]
ideal2_chr3_snp$rate <- NA
ideal2_chr3_snp$`SNP End` <- ideal2_chr3_snp$`SNP End` - min(ideal2_chr3_snp$`SNP Start`)
ideal2_chr3_snp$`SNP Start` <- ideal2_chr3_snp$`SNP Start`- min(ideal2_chr3_snp$`SNP Start`)

ideal2_chr4_snp <- ideal2_snps[ which(ideal2_snps$`Chr#` == "Chr4"),]
ideal2_chr4_snp$rate <- NA
ideal2_chr4_snp$`SNP End` <- ideal2_chr4_snp$`SNP End` - min(ideal2_chr4_snp$`SNP Start`)
ideal2_chr4_snp$`SNP Start` <- ideal2_chr4_snp$`SNP Start`- min(ideal2_chr4_snp$`SNP Start`)

ideal2_chr5_snp <- ideal2_snps[ which(ideal2_snps$`Chr#` == "Chr5"),]
ideal2_chr5_snp$rate <- NA
ideal2_chr5_snp$`SNP End` <- ideal2_chr5_snp$`SNP End` - min(ideal2_chr5_snp$`SNP Start`)
ideal2_chr5_snp$`SNP Start` <- ideal2_chr5_snp$`SNP Start`- min(ideal2_chr5_snp$`SNP Start`)

ideal2_chr6_snp <- ideal2_snps[ which(ideal2_snps$`Chr#` == "Chr6"),]
ideal2_chr6_snp$rate <- NA
ideal2_chr6_snp$`SNP End` <- ideal2_chr6_snp$`SNP End` - min(ideal2_chr6_snp$`SNP Start`)
ideal2_chr6_snp$`SNP Start` <- ideal2_chr6_snp$`SNP Start`- min(ideal2_chr6_snp$`SNP Start`)

ideal2_chr7_snp <- ideal2_snps[ which(ideal2_snps$`Chr#` == "Chr7"),]
ideal2_chr7_snp$rate <- NA
ideal2_chr7_snp$`SNP End` <- ideal2_chr7_snp$`SNP End` - min(ideal2_chr7_snp$`SNP Start`)
ideal2_chr7_snp$`SNP Start` <- ideal2_chr7_snp$`SNP Start`- min(ideal2_chr7_snp$`SNP Start`)

ideal2_chr8_snp <- ideal2_snps[ which(ideal2_snps$`Chr#` == "Chr8"),]
ideal2_chr8_snp$rate <- NA
ideal2_chr8_snp$`SNP End` <- ideal2_chr8_snp$`SNP End` - min(ideal2_chr8_snp$`SNP Start`)
ideal2_chr8_snp$`SNP Start` <- ideal2_chr8_snp$`SNP Start`- min(ideal2_chr8_snp$`SNP Start`)

ideal2_chr9_snp <- ideal2_snps[ which(ideal2_snps$`Chr#` == "Chr9"),]
ideal2_chr9_snp$rate <- NA
ideal2_chr9_snp$`SNP End` <- ideal2_chr9_snp$`SNP End` - min(ideal2_chr9_snp$`SNP Start`)
ideal2_chr9_snp$`SNP Start` <- ideal2_chr9_snp$`SNP Start`- min(ideal2_chr9_snp$`SNP Start`)

ideal2_chr10_snp <- ideal2_snps[ which(ideal2_snps$`Chr#` == "Chr10"),]
ideal2_chr10_snp$rate <- NA
ideal2_chr10_snp$`SNP End` <- ideal2_chr10_snp$`SNP End` - min(ideal2_chr10_snp$`SNP Start`)
ideal2_chr10_snp$`SNP Start` <- ideal2_chr10_snp$`SNP Start`- min(ideal2_chr10_snp$`SNP Start`)

ideal2_chr11_snp <- ideal2_snps[ which(ideal2_snps$`Chr#` == "Chr11"),]
ideal2_chr11_snp$rate <- NA
ideal2_chr11_snp$`SNP End` <- ideal2_chr11_snp$`SNP End` - min(ideal2_chr11_snp$`SNP Start`)
ideal2_chr11_snp$`SNP Start` <- ideal2_chr11_snp$`SNP Start`- min(ideal2_chr11_snp$`SNP Start`)

ideal2_chr12_snp <- ideal2_snps[ which(ideal2_snps$`Chr#` == "Chr12"),]
ideal2_chr12_snp$rate <- NA
ideal2_chr12_snp$`SNP End` <- ideal2_chr12_snp$`SNP End` - min(ideal2_chr12_snp$`SNP Start`)
ideal2_chr12_snp$`SNP Start` <- ideal2_chr12_snp$`SNP Start`- min(ideal2_chr12_snp$`SNP Start`)

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

#apply avg difference to pericentromeric regions (divide chromosome into fifths and apply avg diff to middle fifth)
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
zmet2_avg_diff <-1
pericentromeric <- function(CO){
  rownames(CO)<-c(1:nrow(CO))
  CO$avg_rate <- 0
  fifth<- max(CO$`CO End`)/5
  start<-fifth*2
  end<-fifth*4
  for(i in 1:nrow(CO)){
    if(CO$`CO Start`[i]>= start && CO$`CO End`[i]<= end ){
      CO$avg_rate[i] <- zmet2_avg_diff
    }
  }
  print(CO)
}
ddm1_avg_diff <-2.346557
telomeric <- function(CO){
  rownames(CO)<-c(1:nrow(CO))
  fifth<- max(CO$`CO End`)/5
  start<-fifth*2
  end<-fifth*4
  for(i in 1:nrow(CO)){
    if(CO$`CO Start`[i]<= start){
      CO$avg_rate[i] <- ddm1_avg_diff
    }
    else if(CO$`CO Start`[i]>=end ){
      CO$avg_rate[i] <- ddm1_avg_diff
    }
  }
  print(CO)
}

#creating ideal zmet2/ddm1 double mutant
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
jap_chr1_CO<- telomeric(jap_chr1_CO)
jap_chr2_CO <- telomeric(jap_chr2_CO)
jap_chr3_CO <- telomeric(jap_chr3_CO)
jap_chr4_CO <- telomeric(jap_chr4_CO)
jap_chr5_CO <- telomeric(jap_chr5_CO)
jap_chr6_CO <- telomeric(jap_chr6_CO)
jap_chr7_CO <- telomeric(jap_chr7_CO)
jap_chr8_CO <- telomeric(jap_chr8_CO)
jap_chr9_CO <- telomeric(jap_chr9_CO)
jap_chr10_CO <- telomeric(jap_chr10_CO)
jap_chr11_CO <- telomeric(jap_chr11_CO)
jap_chr12_CO <- telomeric(jap_chr12_CO)

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
ideal2_chr1_CO_2<-new_rates(jap_chr1_CO, WTJap_chr1_CO)
ideal2_chr2_CO_2<-new_rates(jap_chr2_CO, WTJap_chr2_CO)
ideal2_chr3_CO_2<-new_rates(jap_chr3_CO, WTJap_chr3_CO)
ideal2_chr4_CO_2<-new_rates(jap_chr4_CO, WTJap_chr4_CO)
ideal2_chr5_CO_2<-new_rates(jap_chr5_CO, WTJap_chr5_CO)
ideal2_chr6_CO_2<-new_rates(jap_chr6_CO, WTJap_chr6_CO)
ideal2_chr7_CO_2<-new_rates(jap_chr7_CO, WTJap_chr7_CO)
ideal2_chr8_CO_2<-new_rates(jap_chr8_CO, WTJap_chr8_CO)
ideal2_chr9_CO_2<-new_rates(jap_chr9_CO, WTJap_chr9_CO)
ideal2_chr10_CO_2<-new_rates(jap_chr10_CO, WTJap_chr10_CO)
ideal2_chr11_CO_2<-new_rates(jap_chr11_CO, WTJap_chr11_CO)
ideal2_chr12_CO_2<-new_rates(jap_chr12_CO, WTJap_chr12_CO)

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
ideal2_chr1_CO_3 <- ideal2_chr1_CO_2
bins<-as.integer(nrow(ideal2_chr1_CO_2)/40)
ideal2_chr1_CO_3$rates<- rollapply(ideal2_chr1_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ideal2_chr1_CO_3<-fill_start(ideal2_chr1_CO_3)
ideal2_chr1_CO_3<- ideal2_chr1_CO_3 %>% drop_na(rates)

ideal2_chr2_CO_3 <- ideal2_chr2_CO_2
bins<-as.integer(nrow(ideal2_chr2_CO_2)/40)
ideal2_chr2_CO_3$rates<- rollapply(ideal2_chr2_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ideal2_chr2_CO_3<-fill_start(ideal2_chr2_CO_3)
ideal2_chr2_CO_3<- ideal2_chr2_CO_3 %>% drop_na(rates)

ideal2_chr3_CO_3 <- ideal2_chr3_CO_2
bins<-as.integer(nrow(ideal2_chr3_CO_2)/40)
ideal2_chr3_CO_3$rates<- rollapply(ideal2_chr3_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ideal2_chr3_CO_3<-fill_start(ideal2_chr3_CO_3)
ideal2_chr3_CO_3<- ideal2_chr3_CO_3 %>% drop_na(rates)

ideal2_chr4_CO_3 <- ideal2_chr4_CO_2
bins<-as.integer(nrow(ideal2_chr4_CO_2)/40)
ideal2_chr4_CO_3$rates<- rollapply(ideal2_chr4_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ideal2_chr4_CO_3<-fill_start(ideal2_chr4_CO_3)
ideal2_chr4_CO_3<- ideal2_chr4_CO_3 %>% drop_na(rates)

ideal2_chr5_CO_3 <- ideal2_chr5_CO_2
bins<-as.integer(nrow(ideal2_chr5_CO_2)/40)
ideal2_chr5_CO_3$rates<- rollapply(ideal2_chr5_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ideal2_chr5_CO_3<-fill_start(ideal2_chr5_CO_3)
ideal2_chr5_CO_3<- ideal2_chr5_CO_3 %>% drop_na(rates)

ideal2_chr6_CO_3 <- ideal2_chr6_CO_2
bins<-as.integer(nrow(ideal2_chr6_CO_2)/40)
ideal2_chr6_CO_3$rates<- rollapply(ideal2_chr6_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ideal2_chr6_CO_3<-fill_start(ideal2_chr6_CO_3)
ideal2_chr6_CO_3<- ideal2_chr6_CO_3 %>% drop_na(rates)

ideal2_chr7_CO_3 <- ideal2_chr7_CO_2
bins<-as.integer(nrow(ideal2_chr7_CO_2)/40)
ideal2_chr7_CO_3$rates<- rollapply(ideal2_chr7_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ideal2_chr7_CO_3<-fill_start(ideal2_chr7_CO_3)
ideal2_chr7_CO_3<- ideal2_chr7_CO_3 %>% drop_na(rates)

ideal2_chr8_CO_3 <- ideal2_chr8_CO_2
bins<-as.integer(nrow(ideal2_chr8_CO_2)/40)
ideal2_chr8_CO_3$rates<- rollapply(ideal2_chr8_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ideal2_chr8_CO_3<-fill_start(ideal2_chr8_CO_3)
ideal2_chr8_CO_3<- ideal2_chr8_CO_3 %>% drop_na(rates)

ideal2_chr9_CO_3 <- ideal2_chr9_CO_2
bins<-as.integer(nrow(ideal2_chr9_CO_2)/40)
ideal2_chr9_CO_3$rates<- rollapply(ideal2_chr9_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ideal2_chr9_CO_3<-fill_start(ideal2_chr9_CO_3)
ideal2_chr9_CO_3<- ideal2_chr9_CO_3 %>% drop_na(rates)

ideal2_chr10_CO_3 <- ideal2_chr10_CO_2
bins<-as.integer(nrow(ideal2_chr10_CO_2)/40)
ideal2_chr10_CO_3$rates<- rollapply(ideal2_chr10_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ideal2_chr10_CO_3<-fill_start(ideal2_chr10_CO_3)
ideal2_chr10_CO_3<- ideal2_chr10_CO_3 %>% drop_na(rates)

ideal2_chr11_CO_3 <- ideal2_chr11_CO_2
bins<-as.integer(nrow(ideal2_chr11_CO_2)/40)
ideal2_chr11_CO_3$rates<- rollapply(ideal2_chr11_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ideal2_chr11_CO_3<-fill_start(ideal2_chr11_CO_3)
ideal2_chr11_CO_3<- ideal2_chr11_CO_3 %>% drop_na(rates)

ideal2_chr12_CO_3 <- ideal2_chr12_CO_2
bins<-as.integer(nrow(ideal2_chr12_CO_2)/40)
ideal2_chr12_CO_3$rates<- rollapply(ideal2_chr12_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ideal2_chr12_CO_3<-fill_start(ideal2_chr12_CO_3)
ideal2_chr12_CO_3<- ideal2_chr12_CO_3 %>% drop_na(rates)

## assigning frequency to SNPs based on recombination frequency in each bin
snp_rate <- function(chr_rate, chr_snp){
  for(i in 1:nrow(chr_snp)){
    for(k in 1:nrow(chr_rate)){
      if(isTRUE((chr_snp$`SNP Start`[i] >= chr_rate$`CO Start`[k]) && (chr_snp$`SNP Start`[i] <= chr_rate$`CO End`[k]))){
        chr_snp$rate[i] <- chr_rate$rates[k]
      }
    }
  }
  print(chr_snp)
}


#using function,  get cM/Mb for final genetic position - assign rates
ideal2_chr1_snp2 <- snp_rate(ideal2_chr1_CO_3, ideal2_chr1_snp)
ideal2_chr2_snp2 <- snp_rate(ideal2_chr2_CO_3, ideal2_chr2_snp)
ideal2_chr3_snp2 <- snp_rate(ideal2_chr3_CO_3, ideal2_chr3_snp)
ideal2_chr4_snp2 <- snp_rate(ideal2_chr4_CO_3, ideal2_chr4_snp)
ideal2_chr5_snp2 <- snp_rate(ideal2_chr5_CO_3, ideal2_chr5_snp)
ideal2_chr6_snp2 <- snp_rate(ideal2_chr6_CO_3, ideal2_chr6_snp)
ideal2_chr7_snp2 <- snp_rate(ideal2_chr7_CO_3, ideal2_chr7_snp)
ideal2_chr8_snp2 <- snp_rate(ideal2_chr8_CO_3, ideal2_chr8_snp)
ideal2_chr9_snp2 <- snp_rate(ideal2_chr9_CO_3, ideal2_chr9_snp)
ideal2_chr10_snp2 <- snp_rate(ideal2_chr10_CO_3, ideal2_chr10_snp)
ideal2_chr11_snp2 <- snp_rate(ideal2_chr11_CO_3, ideal2_chr11_snp)
ideal2_chr12_snp2 <- snp_rate(ideal2_chr12_CO_3, ideal2_chr12_snp)

#converted SNP start to Mb
ideal2_chr1_snp2$`SNP Start` <- ideal2_chr1_snp2$`SNP Start`/1000000
ideal2_chr2_snp2$`SNP Start` <- ideal2_chr2_snp2$`SNP Start`/1000000
ideal2_chr3_snp2$`SNP Start` <- ideal2_chr3_snp2$`SNP Start`/1000000
ideal2_chr4_snp2$`SNP Start` <- ideal2_chr4_snp2$`SNP Start`/1000000
ideal2_chr5_snp2$`SNP Start` <- ideal2_chr5_snp2$`SNP Start`/1000000
ideal2_chr6_snp2$`SNP Start` <- ideal2_chr6_snp2$`SNP Start`/1000000
ideal2_chr7_snp2$`SNP Start` <- ideal2_chr7_snp2$`SNP Start`/1000000
ideal2_chr8_snp2$`SNP Start` <- ideal2_chr8_snp2$`SNP Start`/1000000
ideal2_chr9_snp2$`SNP Start` <- ideal2_chr9_snp2$`SNP Start`/1000000
ideal2_chr10_snp2$`SNP Start` <- ideal2_chr10_snp2$`SNP Start`/1000000
ideal2_chr11_snp2$`SNP Start` <- ideal2_chr11_snp2$`SNP Start`/1000000
ideal2_chr12_snp2$`SNP Start` <- ideal2_chr12_snp2$`SNP Start`/1000000

ideal2_chr1_snp2$`SNP End` <- ideal2_chr1_snp2$`SNP End`/1000000
ideal2_chr2_snp2$`SNP End` <- ideal2_chr2_snp2$`SNP End`/1000000
ideal2_chr3_snp2$`SNP End` <- ideal2_chr3_snp2$`SNP End`/1000000
ideal2_chr4_snp2$`SNP End` <- ideal2_chr4_snp2$`SNP End`/1000000
ideal2_chr5_snp2$`SNP End` <- ideal2_chr5_snp2$`SNP End`/1000000
ideal2_chr6_snp2$`SNP End` <- ideal2_chr6_snp2$`SNP End`/1000000
ideal2_chr7_snp2$`SNP End` <- ideal2_chr7_snp2$`SNP End`/1000000
ideal2_chr8_snp2$`SNP End` <- ideal2_chr8_snp2$`SNP End`/1000000
ideal2_chr9_snp2$`SNP End` <- ideal2_chr9_snp2$`SNP End`/1000000
ideal2_chr10_snp2$`SNP End` <- ideal2_chr10_snp2$`SNP End`/1000000
ideal2_chr11_snp2$`SNP End` <- ideal2_chr11_snp2$`SNP End`/1000000
ideal2_chr12_snp2$`SNP End` <- ideal2_chr12_snp2$`SNP End`/1000000

#omit empty col
ideal2_chr1_snp2<-na.omit(ideal2_chr1_snp2)
ideal2_chr2_snp2<-na.omit(ideal2_chr2_snp2)
ideal2_chr3_snp2<-na.omit(ideal2_chr3_snp2)
ideal2_chr4_snp2<-na.omit(ideal2_chr4_snp2)
ideal2_chr5_snp2<-na.omit(ideal2_chr5_snp2)
ideal2_chr6_snp2<-na.omit(ideal2_chr6_snp2)
ideal2_chr7_snp2<-na.omit(ideal2_chr7_snp2)
ideal2_chr8_snp2<-na.omit(ideal2_chr8_snp2)
ideal2_chr9_snp2<-na.omit(ideal2_chr9_snp2)
ideal2_chr10_snp2<-na.omit(ideal2_chr10_snp2)
ideal2_chr11_snp2<-na.omit(ideal2_chr11_snp2)
ideal2_chr12_snp2<-na.omit(ideal2_chr12_snp2)

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


ideal2_chr1_spl <- smooth.spline(ideal2_chr1_snp2$rate, spar = .1)
ideal2_chr1_snp2$pos <-gen_pos(ideal2_chr1_snp2)
plot(ideal2_chr1_snp2$`SNP Start`, ideal2_chr1_snp2$pos)
ggplot(ideal2_chr1_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
plot(ideal2_chr1_snp2$`SNP Start`, ideal2_chr1_snp2$pos/ideal2_chr1_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica ideal2 Chromosome 1 Recombination Distribution")
ideal2_chr1_finalpos <- ideal2_chr1_snp2[order(ideal2_chr1_snp2$pos),]
is.unsorted(ideal2_chr1_finalpos$pos)
plot(ideal2_chr1_snp2$`SNP Start`, ideal2_chr1_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ideal2 Chromosome 1 Genetic Map")
plot(ideal2_chr1_finalpos$`SNP Start`, ideal2_chr1_finalpos$pos)


ideal2_chr2_spl <- smooth.spline(ideal2_chr2_snp2$rate, spar = .1)
ideal2_chr2_snp2$pos <-gen_pos(ideal2_chr2_snp2)
plot(ideal2_chr2_snp2$`SNP Start`, ideal2_chr2_snp2$pos)
plot(ideal2_chr2_snp2$`SNP Start`, ideal2_chr2_snp2$pos/ideal2_chr2_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica ideal2 Chromosome 2 Recombination Distribution")
ideal2_chr2_finalpos <- ideal2_chr2_snp2[order(ideal2_chr2_snp2$pos),]
is.unsorted(ideal2_chr2_finalpos$pos)
plot(ideal2_chr2_snp2$`SNP Start`, ideal2_chr2_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ideal2 Chromosome 2 Genetic Map")

ideal2_chr3_spl <- smooth.spline(ideal2_chr3_snp2$rate, spar = .1)
ideal2_chr3_snp2$pos <-gen_pos(ideal2_chr3_snp2)
plot(ideal2_chr3_snp2$`SNP Start`, ideal2_chr3_snp2$pos)
plot(ideal2_chr3_snp2$`SNP Start`, ideal2_chr3_snp2$pos/ideal2_chr3_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica ideal2 Chromosome 3 Recombination Distribution")
ideal2_chr3_finalpos <- ideal2_chr3_snp2[order(ideal2_chr3_snp2$pos),]
is.unsorted(ideal2_chr3_finalpos$pos)
plot(ideal2_chr3_snp2$`SNP Start`, ideal2_chr3_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ideal2 Chromosome 3 Genetic Map")

ideal2_chr4_spl <- smooth.spline(ideal2_chr4_snp2$rate, spar = .1)
ideal2_chr4_snp2$pos <-gen_pos(ideal2_chr4_snp2)
plot(ideal2_chr4_snp2$`SNP Start`, ideal2_chr4_snp2$pos)
plot(ideal2_chr4_snp2$`SNP Start`, ideal2_chr4_snp2$pos/ideal2_chr4_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica ideal2 Chromosome 4 Recombination Distribution")
ideal2_chr4_finalpos <- ideal2_chr4_snp2[order(ideal2_chr4_snp2$pos),]
is.unsorted(ideal2_chr4_finalpos$pos)
plot(ideal2_chr4_snp2$`SNP Start`, ideal2_chr4_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ideal2 Chromosome 4 Genetic Map")

ideal2_chr5_spl <- smooth.spline(ideal2_chr5_snp2$rate, spar =.1)
ideal2_chr5_snp2$pos <-gen_pos(ideal2_chr5_snp2)
plot(ideal2_chr5_snp2$`SNP Start`, ideal2_chr5_snp2$pos)
plot(ideal2_chr5_snp2$`SNP Start`, ideal2_chr5_snp2$pos/ideal2_chr5_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica ideal2 Chromosome 5 Recombination Distribution")
ideal2_chr5_finalpos <- ideal2_chr5_snp2[order(ideal2_chr5_snp2$pos),]
is.unsorted(ideal2_chr5_finalpos$pos)
plot(ideal2_chr5_snp2$`SNP Start`, ideal2_chr5_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ideal2 Chromosome 5 Genetic Map")

ideal2_chr6_spl <- smooth.spline(ideal2_chr6_snp2$rate, spar = .1)
ideal2_chr6_snp2$pos <-gen_pos(ideal2_chr6_snp2)
plot(ideal2_chr6_snp2$`SNP Start`, ideal2_chr6_snp2$pos)
plot(ideal2_chr6_snp2$`SNP Start`, ideal2_chr6_snp2$pos/ideal2_chr6_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica ideal2 Chromosome 6 Recombination Distribution")
ideal2_chr6_finalpos <- ideal2_chr6_snp2[order(ideal2_chr6_snp2$pos),]
is.unsorted(ideal2_chr6_finalpos$pos)
plot(ideal2_chr6_snp2$`SNP Start`, ideal2_chr6_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ideal2 Chromosome 6 Genetic Map")

ideal2_chr7_spl <- smooth.spline(ideal2_chr7_snp2$rate, spar = .1)
ideal2_chr7_snp2$pos <-gen_pos(ideal2_chr7_snp2)
plot(ideal2_chr7_snp2$`SNP Start`, ideal2_chr7_snp2$pos)
plot(ideal2_chr7_snp2$`SNP Start`, ideal2_chr7_snp2$pos/ideal2_chr7_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica ideal2 Chromosome 7 Recombination Distribution")
ideal2_chr7_finalpos <- ideal2_chr7_snp2[order(ideal2_chr7_snp2$pos),]
is.unsorted(ideal2_chr7_finalpos$pos)
plot(ideal2_chr7_snp2$`SNP Start`, ideal2_chr7_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ideal2 Chromosome 7 Genetic Map")

ideal2_chr8_spl <- smooth.spline(ideal2_chr8_snp2$rate, spar = .1)
ideal2_chr8_snp2$pos <-gen_pos(ideal2_chr8_snp2)
plot(ideal2_chr8_snp2$`SNP Start`, ideal2_chr8_snp2$pos)
plot(ideal2_chr8_snp2$`SNP Start`, ideal2_chr8_snp2$pos/ideal2_chr8_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica ideal2 Chromosome 8 Recombination Distribution")
ideal2_chr8_finalpos <- ideal2_chr8_snp2[order(ideal2_chr8_snp2$pos),]
is.unsorted(ideal2_chr8_finalpos$pos)
plot(ideal2_chr8_snp2$`SNP Start`, ideal2_chr8_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ideal2 Chromosome 8 Genetic Map")

ideal2_chr9_spl <- smooth.spline(ideal2_chr9_snp2$rate, spar = .1)
ideal2_chr9_snp2$pos <-gen_pos(ideal2_chr9_snp2)
plot(ideal2_chr9_snp2$`SNP Start`, ideal2_chr9_snp2$pos)
plot(ideal2_chr9_snp2$`SNP Start`, ideal2_chr9_snp2$pos/ideal2_chr9_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica ideal2 Chromosome 9 Recombination Distribution")
ideal2_chr9_finalpos <- ideal2_chr9_snp2[order(ideal2_chr9_snp2$pos),]
is.unsorted(ideal2_chr9_finalpos$pos)
plot(ideal2_chr9_snp2$`SNP Start`, ideal2_chr9_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ideal2 Chromosome 9 Genetic Map")

ideal2_chr10_spl <- smooth.spline(ideal2_chr10_snp2$rate, spar =.1)
ideal2_chr10_snp2$pos <-gen_pos(ideal2_chr10_snp2)
plot(ideal2_chr10_snp2$`SNP Start`, ideal2_chr10_snp2$pos)
plot(ideal2_chr10_snp2$`SNP Start`, ideal2_chr10_snp2$pos/ideal2_chr10_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica ideal2 Chromosome 10 Recombination Distribution")
ideal2_chr10_finalpos <- ideal2_chr10_snp2[order(ideal2_chr10_snp2$pos),]
is.unsorted(ideal2_chr10_finalpos$pos)
plot(ideal2_chr10_snp2$`SNP Start`, ideal2_chr10_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ideal2 Chromosome 10 Genetic Map")

ideal2_chr11_spl <- smooth.spline(ideal2_chr11_snp2$rate, spar = .1)
ideal2_chr11_snp2$pos <-gen_pos(ideal2_chr11_snp2)
plot(ideal2_chr11_snp2$`SNP Start`, ideal2_chr11_snp2$pos)
plot(ideal2_chr11_snp2$`SNP Start`, ideal2_chr11_snp2$pos/ideal2_chr11_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica ideal2 Chromosome 11 Recombination Distribution")
ideal2_chr11_finalpos <- ideal2_chr11_snp2[order(ideal2_chr11_snp2$pos),]
is.unsorted(ideal2_chr11_finalpos$pos)
plot(ideal2_chr11_snp2$`SNP Start`, ideal2_chr11_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ideal2 Chromosome 11 Genetic Map")

ideal2_chr12_spl <- smooth.spline(ideal2_chr12_snp2$rate, spar = 0.1)
ideal2_chr12_snp2$pos <-gen_pos(ideal2_chr12_snp2)
plot(ideal2_chr12_snp2$`SNP Start`, ideal2_chr12_snp2$pos)
plot(ideal2_chr12_snp2$`SNP Start`, ideal2_chr12_snp2$pos/ideal2_chr12_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica ideal2 Chromosome 12 Recombination Distribution")
ideal2_chr12_finalpos <- ideal2_chr12_snp2[order(ideal2_chr12_snp2$pos),]
is.unsorted(ideal2_chr12_finalpos$pos)
plot(ideal2_chr12_snp2$`SNP Start`, ideal2_chr12_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ideal2 Chromosome 12 Genetic Map")

#Final genetic map
ideal2_chr1 <- ideal2_chr1_finalpos$pos/100
ideal2_chr1len <- length(ideal2_chr1)
dim(ideal2_chr1) <- c(ideal2_chr1len,1)
ideal2_chr1 <- list(ideal2_chr1)

ideal2_chr2 <- ideal2_chr2_finalpos$pos/100
ideal2_chr2len <- length(ideal2_chr2)
dim(ideal2_chr2) <- c(ideal2_chr2len,1)
ideal2_chr2 <- list(ideal2_chr2)

ideal2_chr3 <- ideal2_chr3_finalpos$pos/100
ideal2_chr3len <- length(ideal2_chr3)
dim(ideal2_chr3) <- c(ideal2_chr3len,1)
ideal2_chr3 <- list(ideal2_chr3)

ideal2_chr4 <- ideal2_chr4_finalpos$pos/100
ideal2_chr4len <- length(ideal2_chr4)
dim(ideal2_chr4) <- c(ideal2_chr4len,1)
ideal2_chr4 <- list(ideal2_chr4)

ideal2_chr5 <- ideal2_chr5_finalpos$pos/100
ideal2_chr5len <- length(ideal2_chr5)
dim(ideal2_chr5) <- c(ideal2_chr5len,1)
ideal2_chr5 <- list(ideal2_chr5)

ideal2_chr5 <- ideal2_chr5_finalpos$pos/100
ideal2_chr5len <- length(ideal2_chr5)
dim(ideal2_chr5) <- c(ideal2_chr5len,1)
ideal2_chr5 <- list(ideal2_chr5)

ideal2_chr7 <- ideal2_chr7_finalpos$pos/100
ideal2_chr7len <- length(ideal2_chr7)
dim(ideal2_chr7) <- c(ideal2_chr7len,1)
ideal2_chr7 <- list(ideal2_chr7)

ideal2_chr8 <- ideal2_chr8_finalpos$pos/100
ideal2_chr8len <- length(ideal2_chr8)
dim(ideal2_chr8) <- c(ideal2_chr8len,1)
ideal2_chr8 <- list(ideal2_chr8)

ideal2_chr9 <- ideal2_chr9_finalpos$pos/100
ideal2_chr9len <- length(ideal2_chr9)
dim(ideal2_chr9) <- c(ideal2_chr9len,1)
ideal2_chr9 <- list(ideal2_chr9)

ideal2_chr10 <- ideal2_chr10_finalpos$pos/100
ideal2_chr10len <- length(ideal2_chr10)
dim(ideal2_chr10) <- c(ideal2_chr10len,1)
ideal2_chr10 <- list(ideal2_chr10)

ideal2_chr11 <- ideal2_chr11_finalpos$pos/100
ideal2_chr11len <- length(ideal2_chr11)
dim(ideal2_chr11) <- c(ideal2_chr11len,1)
ideal2_chr11 <- list(ideal2_chr11)

ideal2_chr12 <- ideal2_chr12_finalpos$pos/100
ideal2_chr12len <- length(ideal2_chr12)
dim(ideal2_chr12) <- c(ideal2_chr12len,1)
ideal2_chr12 <- list(ideal2_chr12)

ideal2_final_map <- list(ideal2_chr1[[1]], ideal2_chr2[[1]], 
                        ideal2_chr3[[1]], ideal2_chr4[[1]], ideal2_chr5[[1]], 
                        ideal2_chr5[[1]], ideal2_chr7[[1]], ideal2_chr8[[1]], 
                        ideal2_chr9[[1]], ideal2_chr10[[1]],ideal2_chr11[[1]], ideal2_chr12[[1]])

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
c1 <-find_centromere(16.7,ideal2_chr1_finalpos)
c2 <-find_centromere(13.6,ideal2_chr2_finalpos)
c3 <-find_centromere(19.4,ideal2_chr3_finalpos)
c4 <-find_centromere(9.7,ideal2_chr4_finalpos)
c5 <-find_centromere(12.4,ideal2_chr5_finalpos)
c6 <-find_centromere(15.3,ideal2_chr6_finalpos)
c7 <-find_centromere(12.1,ideal2_chr7_finalpos)
c8 <-find_centromere(12.9,ideal2_chr8_finalpos)
c9 <-find_centromere(2.8,ideal2_chr9_finalpos)
c10 <-find_centromere(8.2,ideal2_chr10_finalpos)
c11 <-find_centromere(12,ideal2_chr11_finalpos)
c12 <-find_centromere(11.9,ideal2_chr12_finalpos)

ideal2_centromere <- c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12)
ideal2_centromere <- ideal2_centromere/100


