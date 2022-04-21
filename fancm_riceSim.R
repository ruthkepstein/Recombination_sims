library(AlphaSimR)
library(Rcpp)
library(ggplot2)
library(dplyr)

set.seed(420)

japonica_snps <- read.table("japonica_SNPs.bed", header =FALSE)
colnames(japonica_snps) <- c("Chr#", "SNP Start", "SNP End")
fancm_snps <- sample_n(japonica_snps, 2000)
fancm_snps <- fancm_snps[order(fancm_snps$`Chr#`,fancm_snps$`SNP Start`),]

#splitting fancm snps
fancm_chr1_snp <- fancm_snps[ which(fancm_snps$`Chr#` == "Chr1"),]
fancm_chr1_snp$rate <- NA
fancm_chr1_snp$`SNP End` <- fancm_chr1_snp$`SNP End` - min(fancm_chr1_snp$`SNP Start`)
fancm_chr1_snp$`SNP Start` <- fancm_chr1_snp$`SNP Start`- min(fancm_chr1_snp$`SNP Start`)

fancm_chr2_snp <- fancm_snps[ which(fancm_snps$`Chr#` == "Chr2"),]
fancm_chr2_snp$rate <- NA
fancm_chr2_snp$`SNP End` <- fancm_chr2_snp$`SNP End` - min(fancm_chr2_snp$`SNP Start`)
fancm_chr2_snp$`SNP Start` <- fancm_chr2_snp$`SNP Start`- min(fancm_chr2_snp$`SNP Start`)

fancm_chr3_snp <- fancm_snps[ which(fancm_snps$`Chr#` == "Chr3"),]
fancm_chr3_snp$rate <- NA
fancm_chr3_snp$`SNP End` <- fancm_chr3_snp$`SNP End` - min(fancm_chr3_snp$`SNP Start`)
fancm_chr3_snp$`SNP Start` <- fancm_chr3_snp$`SNP Start`- min(fancm_chr3_snp$`SNP Start`)

fancm_chr4_snp <- fancm_snps[ which(fancm_snps$`Chr#` == "Chr4"),]
fancm_chr4_snp$rate <- NA
fancm_chr4_snp$`SNP End` <- fancm_chr4_snp$`SNP End` - min(fancm_chr4_snp$`SNP Start`)
fancm_chr4_snp$`SNP Start` <- fancm_chr4_snp$`SNP Start`- min(fancm_chr4_snp$`SNP Start`)

fancm_chr5_snp <- fancm_snps[ which(fancm_snps$`Chr#` == "Chr5"),]
fancm_chr5_snp$rate <- NA
fancm_chr5_snp$`SNP End` <- fancm_chr5_snp$`SNP End` - min(fancm_chr5_snp$`SNP Start`)
fancm_chr5_snp$`SNP Start` <- fancm_chr5_snp$`SNP Start`- min(fancm_chr5_snp$`SNP Start`)

fancm_chr6_snp <- fancm_snps[ which(fancm_snps$`Chr#` == "Chr6"),]
fancm_chr6_snp$rate <- NA
fancm_chr6_snp$`SNP End` <- fancm_chr6_snp$`SNP End` - min(fancm_chr6_snp$`SNP Start`)
fancm_chr6_snp$`SNP Start` <- fancm_chr6_snp$`SNP Start`- min(fancm_chr6_snp$`SNP Start`)

fancm_chr7_snp <- fancm_snps[ which(fancm_snps$`Chr#` == "Chr7"),]
fancm_chr7_snp$rate <- NA
fancm_chr7_snp$`SNP End` <- fancm_chr7_snp$`SNP End` - min(fancm_chr7_snp$`SNP Start`)
fancm_chr7_snp$`SNP Start` <- fancm_chr7_snp$`SNP Start`- min(fancm_chr7_snp$`SNP Start`)

fancm_chr8_snp <- fancm_snps[ which(fancm_snps$`Chr#` == "Chr8"),]
fancm_chr8_snp$rate <- NA
fancm_chr8_snp$`SNP End` <- fancm_chr8_snp$`SNP End` - min(fancm_chr8_snp$`SNP Start`)
fancm_chr8_snp$`SNP Start` <- fancm_chr8_snp$`SNP Start`- min(fancm_chr8_snp$`SNP Start`)

fancm_chr9_snp <- fancm_snps[ which(fancm_snps$`Chr#` == "Chr9"),]
fancm_chr9_snp$rate <- NA
fancm_chr9_snp$`SNP End` <- fancm_chr9_snp$`SNP End` - min(fancm_chr9_snp$`SNP Start`)
fancm_chr9_snp$`SNP Start` <- fancm_chr9_snp$`SNP Start`- min(fancm_chr9_snp$`SNP Start`)

fancm_chr10_snp <- fancm_snps[ which(fancm_snps$`Chr#` == "Chr10"),]
fancm_chr10_snp$rate <- NA
fancm_chr10_snp$`SNP End` <- fancm_chr10_snp$`SNP End` - min(fancm_chr10_snp$`SNP Start`)
fancm_chr10_snp$`SNP Start` <- fancm_chr10_snp$`SNP Start`- min(fancm_chr10_snp$`SNP Start`)

fancm_chr11_snp <- fancm_snps[ which(fancm_snps$`Chr#` == "Chr11"),]
fancm_chr11_snp$rate <- NA
fancm_chr11_snp$`SNP End` <- fancm_chr11_snp$`SNP End` - min(fancm_chr11_snp$`SNP Start`)
fancm_chr11_snp$`SNP Start` <- fancm_chr11_snp$`SNP Start`- min(fancm_chr11_snp$`SNP Start`)

fancm_chr12_snp <- fancm_snps[ which(fancm_snps$`Chr#` == "Chr12"),]
fancm_chr12_snp$rate <- NA
fancm_chr12_snp$`SNP End` <- fancm_chr12_snp$`SNP End` - min(fancm_chr12_snp$`SNP Start`)
fancm_chr12_snp$`SNP Start` <- fancm_chr12_snp$`SNP Start`- min(fancm_chr12_snp$`SNP Start`)


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

#apply avg difference to telomeric regions (divide chromosome into fifths and apply avg diff to middle fifth)
avg_diff <-1.252264
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
fancm_chr1_CO_2<-new_rates(WTJap_chr1_CO)
fancm_chr2_CO_2<-new_rates(WTJap_chr2_CO)
fancm_chr3_CO_2<-new_rates(WTJap_chr3_CO)
fancm_chr4_CO_2<-new_rates(WTJap_chr4_CO)
fancm_chr5_CO_2<-new_rates(WTJap_chr5_CO)
fancm_chr6_CO_2<-new_rates(WTJap_chr6_CO)
fancm_chr7_CO_2<-new_rates(WTJap_chr7_CO)
fancm_chr8_CO_2<-new_rates(WTJap_chr8_CO)
fancm_chr9_CO_2<-new_rates(WTJap_chr9_CO)
fancm_chr10_CO_2<-new_rates(WTJap_chr10_CO)
fancm_chr11_CO_2<-new_rates(WTJap_chr11_CO)
fancm_chr12_CO_2<-new_rates(WTJap_chr12_CO)

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
fancm_chr1_CO_3 <- fancm_chr1_CO_2
bins<-as.integer(nrow(fancm_chr1_CO_2)/10)
fancm_chr1_CO_3$rates<- rollapply(fancm_chr1_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
fancm_chr1_CO_3<-fill_start(fancm_chr1_CO_3)
fancm_chr1_CO_3<- fancm_chr1_CO_3 %>% drop_na(rates)

fancm_chr2_CO_3 <- fancm_chr2_CO_2
bins<-as.integer(nrow(fancm_chr2_CO_2)/10)
fancm_chr2_CO_3$rates<- rollapply(fancm_chr2_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
fancm_chr2_CO_3<-fill_start(fancm_chr2_CO_3)
fancm_chr2_CO_3<- fancm_chr2_CO_3 %>% drop_na(rates)

fancm_chr3_CO_3 <- fancm_chr3_CO_2
bins<-as.integer(nrow(fancm_chr3_CO_2)/10)
fancm_chr3_CO_3$rates<- rollapply(fancm_chr3_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
fancm_chr3_CO_3<-fill_start(fancm_chr3_CO_3)
fancm_chr3_CO_3<- fancm_chr3_CO_3 %>% drop_na(rates)

fancm_chr4_CO_3 <- fancm_chr4_CO_2
bins<-as.integer(nrow(fancm_chr4_CO_2)/10)
fancm_chr4_CO_3$rates<- rollapply(fancm_chr4_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
fancm_chr4_CO_3<-fill_start(fancm_chr4_CO_3)
fancm_chr4_CO_3<- fancm_chr4_CO_3 %>% drop_na(rates)

fancm_chr5_CO_3 <- fancm_chr5_CO_2
bins<-as.integer(nrow(fancm_chr5_CO_2)/10)
fancm_chr5_CO_3$rates<- rollapply(fancm_chr5_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
fancm_chr5_CO_3<-fill_start(fancm_chr5_CO_3)
fancm_chr5_CO_3<- fancm_chr5_CO_3 %>% drop_na(rates)

fancm_chr6_CO_3 <- fancm_chr6_CO_2
bins<-as.integer(nrow(fancm_chr6_CO_2)/10)
fancm_chr6_CO_3$rates<- rollapply(fancm_chr6_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
fancm_chr6_CO_3<-fill_start(fancm_chr6_CO_3)
fancm_chr6_CO_3<- fancm_chr6_CO_3 %>% drop_na(rates)

fancm_chr7_CO_3 <- fancm_chr7_CO_2
bins<-as.integer(nrow(fancm_chr7_CO_2)/10)
fancm_chr7_CO_3$rates<- rollapply(fancm_chr7_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
fancm_chr7_CO_3<-fill_start(fancm_chr7_CO_3)
fancm_chr7_CO_3<- fancm_chr7_CO_3 %>% drop_na(rates)

fancm_chr8_CO_3 <- fancm_chr8_CO_2
bins<-as.integer(nrow(fancm_chr8_CO_2)/10)
fancm_chr8_CO_3$rates<- rollapply(fancm_chr8_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
fancm_chr8_CO_3<-fill_start(fancm_chr8_CO_3)
fancm_chr8_CO_3<- fancm_chr8_CO_3 %>% drop_na(rates)

fancm_chr9_CO_3 <- fancm_chr9_CO_2
bins<-as.integer(nrow(fancm_chr9_CO_2)/10)
fancm_chr9_CO_3$rates<- rollapply(fancm_chr9_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
fancm_chr9_CO_3<-fill_start(fancm_chr9_CO_3)
fancm_chr9_CO_3<- fancm_chr9_CO_3 %>% drop_na(rates)

fancm_chr10_CO_3 <- fancm_chr10_CO_2
bins<-as.integer(nrow(fancm_chr10_CO_2)/10)
fancm_chr10_CO_3$rates<- rollapply(fancm_chr10_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
fancm_chr10_CO_3<-fill_start(fancm_chr10_CO_3)
fancm_chr10_CO_3<- fancm_chr10_CO_3 %>% drop_na(rates)

fancm_chr11_CO_3 <- fancm_chr11_CO_2
bins<-as.integer(nrow(fancm_chr11_CO_2)/10)
fancm_chr11_CO_3$rates<- rollapply(fancm_chr11_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
fancm_chr11_CO_3<-fill_start(fancm_chr11_CO_3)
fancm_chr11_CO_3<- fancm_chr11_CO_3 %>% drop_na(rates)

fancm_chr12_CO_3 <- fancm_chr12_CO_2
bins<-as.integer(nrow(fancm_chr12_CO_2)/10)
fancm_chr12_CO_3$rates<- rollapply(fancm_chr12_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
fancm_chr12_CO_3<-fill_start(fancm_chr12_CO_3)
fancm_chr12_CO_3<- fancm_chr12_CO_3 %>% drop_na(rates)

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
fancm_chr1_snp2 <- snp_rate(fancm_chr1_CO_3, fancm_chr1_snp)
fancm_chr2_snp2 <- snp_rate(fancm_chr2_CO_3, fancm_chr2_snp)
fancm_chr3_snp2 <- snp_rate(fancm_chr3_CO_3, fancm_chr3_snp)
fancm_chr4_snp2 <- snp_rate(fancm_chr4_CO_3, fancm_chr4_snp)
fancm_chr5_snp2 <- snp_rate(fancm_chr5_CO_3, fancm_chr5_snp)
fancm_chr6_snp2 <- snp_rate(fancm_chr6_CO_3, fancm_chr6_snp)
fancm_chr7_snp2 <- snp_rate(fancm_chr7_CO_3, fancm_chr7_snp)
fancm_chr8_snp2 <- snp_rate(fancm_chr8_CO_3, fancm_chr8_snp)
fancm_chr9_snp2 <- snp_rate(fancm_chr9_CO_3, fancm_chr9_snp)
fancm_chr10_snp2 <- snp_rate(fancm_chr10_CO_3, fancm_chr10_snp)
fancm_chr11_snp2 <- snp_rate(fancm_chr11_CO_3, fancm_chr11_snp)
fancm_chr12_snp2 <- snp_rate(fancm_chr12_CO_3, fancm_chr12_snp)

#converted SNP start to Mb
fancm_chr1_snp2$`SNP Start` <- fancm_chr1_snp2$`SNP Start`/1000000
fancm_chr2_snp2$`SNP Start` <- fancm_chr2_snp2$`SNP Start`/1000000
fancm_chr3_snp2$`SNP Start` <- fancm_chr3_snp2$`SNP Start`/1000000
fancm_chr4_snp2$`SNP Start` <- fancm_chr4_snp2$`SNP Start`/1000000
fancm_chr5_snp2$`SNP Start` <- fancm_chr5_snp2$`SNP Start`/1000000
fancm_chr6_snp2$`SNP Start` <- fancm_chr6_snp2$`SNP Start`/1000000
fancm_chr7_snp2$`SNP Start` <- fancm_chr7_snp2$`SNP Start`/1000000
fancm_chr8_snp2$`SNP Start` <- fancm_chr8_snp2$`SNP Start`/1000000
fancm_chr9_snp2$`SNP Start` <- fancm_chr9_snp2$`SNP Start`/1000000
fancm_chr10_snp2$`SNP Start` <- fancm_chr10_snp2$`SNP Start`/1000000
fancm_chr11_snp2$`SNP Start` <- fancm_chr11_snp2$`SNP Start`/1000000
fancm_chr12_snp2$`SNP Start` <- fancm_chr12_snp2$`SNP Start`/1000000

fancm_chr1_snp2$`SNP End` <- fancm_chr1_snp2$`SNP End`/1000000
fancm_chr2_snp2$`SNP End` <- fancm_chr2_snp2$`SNP End`/1000000
fancm_chr3_snp2$`SNP End` <- fancm_chr3_snp2$`SNP End`/1000000
fancm_chr4_snp2$`SNP End` <- fancm_chr4_snp2$`SNP End`/1000000
fancm_chr5_snp2$`SNP End` <- fancm_chr5_snp2$`SNP End`/1000000
fancm_chr6_snp2$`SNP End` <- fancm_chr6_snp2$`SNP End`/1000000
fancm_chr7_snp2$`SNP End` <- fancm_chr7_snp2$`SNP End`/1000000
fancm_chr8_snp2$`SNP End` <- fancm_chr8_snp2$`SNP End`/1000000
fancm_chr9_snp2$`SNP End` <- fancm_chr9_snp2$`SNP End`/1000000
fancm_chr10_snp2$`SNP End` <- fancm_chr10_snp2$`SNP End`/1000000
fancm_chr11_snp2$`SNP End` <- fancm_chr11_snp2$`SNP End`/1000000
fancm_chr12_snp2$`SNP End` <- fancm_chr12_snp2$`SNP End`/1000000

#omit empty col
fancm_chr1_snp2<-na.omit(fancm_chr1_snp2)
fancm_chr2_snp2<-na.omit(fancm_chr2_snp2)
fancm_chr3_snp2<-na.omit(fancm_chr3_snp2)
fancm_chr4_snp2<-na.omit(fancm_chr4_snp2)
fancm_chr5_snp2<-na.omit(fancm_chr5_snp2)
fancm_chr6_snp2<-na.omit(fancm_chr6_snp2)
fancm_chr7_snp2<-na.omit(fancm_chr7_snp2)
fancm_chr8_snp2<-na.omit(fancm_chr8_snp2)
fancm_chr9_snp2<-na.omit(fancm_chr9_snp2)
fancm_chr10_snp2<-na.omit(fancm_chr10_snp2)
fancm_chr11_snp2<-na.omit(fancm_chr11_snp2)
fancm_chr12_snp2<-na.omit(fancm_chr12_snp2)

#gen maps
fancm_chr1_spl <- smooth.spline(fancm_chr1_snp2$rate, spar = .4)
fancm_chr1_snp2$pos <- (fancm_chr1_snp2$`SNP Start`*fancm_chr1_spl$y)
plot(fancm_chr1_snp2$`SNP Start`, fancm_chr1_snp2$pos)
ggplot(fancm_chr1_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
plot(fancm_chr1_snp2$`SNP Start`, fancm_chr1_snp2$pos/fancm_chr1_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica fancm Chromosome 1 Recombination Distribution")
fancm_chr1_finalpos <- fancm_chr1_snp2[order(fancm_chr1_snp2$pos),]
is.unsorted(fancm_chr1_finalpos$pos)
plot(fancm_chr1_snp2$`SNP Start`, fancm_chr1_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica fancm Chromosome 1 Genetic Map")
plot(fancm_chr1_finalpos$`SNP Start`, fancm_chr1_finalpos$pos)

fancm_chr2_spl <- smooth.spline(fancm_chr2_snp2$rate, spar = .4)
fancm_chr2_snp2$pos <- (fancm_chr2_snp2$`SNP Start`*fancm_chr2_spl$y)
plot(fancm_chr2_snp2$`SNP Start`, fancm_chr2_snp2$pos)
plot(fancm_chr2_snp2$`SNP Start`, fancm_chr2_snp2$pos/fancm_chr2_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica fancm Chromosome 2 Recombination Distribution")
fancm_chr2_finalpos <- fancm_chr2_snp2[order(fancm_chr2_snp2$pos),]
is.unsorted(fancm_chr2_finalpos$pos)
plot(fancm_chr2_snp2$`SNP Start`, fancm_chr2_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica fancm Chromosome 2 Genetic Map")

fancm_chr3_spl <- smooth.spline(fancm_chr3_snp2$rate, spar = .4)
fancm_chr3_snp2$pos <- (fancm_chr3_snp2$`SNP Start`*fancm_chr3_spl$y)
plot(fancm_chr3_snp2$`SNP Start`, fancm_chr3_snp2$pos)
plot(fancm_chr3_snp2$`SNP Start`, fancm_chr3_snp2$pos/fancm_chr3_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica fancm Chromosome 3 Recombination Distribution")
fancm_chr3_finalpos <- fancm_chr3_snp2[order(fancm_chr3_snp2$pos),]
is.unsorted(fancm_chr3_finalpos$pos)
plot(fancm_chr3_snp2$`SNP Start`, fancm_chr3_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica fancm Chromosome 3 Genetic Map")

fancm_chr4_spl <- smooth.spline(fancm_chr4_snp2$rate, spar = .4)
fancm_chr4_snp2$pos <- (fancm_chr4_snp2$`SNP Start`*fancm_chr4_spl$y)
plot(fancm_chr4_snp2$`SNP Start`, fancm_chr4_snp2$pos)
plot(fancm_chr4_snp2$`SNP Start`, fancm_chr4_snp2$pos/fancm_chr4_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica fancm Chromosome 4 Recombination Distribution")
fancm_chr4_finalpos <- fancm_chr4_snp2[order(fancm_chr4_snp2$pos),]
is.unsorted(fancm_chr4_finalpos$pos)
fancm_chr4_finalpos$pos <- fancm_chr4_finalpos$pos + abs(min(fancm_chr4_finalpos$pos))
plot(fancm_chr4_snp2$`SNP Start`, fancm_chr4_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica fancm Chromosome 4 Genetic Map")

fancm_chr5_spl <- smooth.spline(fancm_chr5_snp2$rate, spar =.4)
fancm_chr5_snp2$pos <- (fancm_chr5_snp2$`SNP Start`*fancm_chr5_spl$y)
plot(fancm_chr5_snp2$`SNP Start`, fancm_chr5_snp2$pos)
plot(fancm_chr5_snp2$`SNP Start`, fancm_chr5_snp2$pos/fancm_chr5_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica fancm Chromosome 5 Recombination Distribution")
fancm_chr5_finalpos <- fancm_chr5_snp2[order(fancm_chr5_snp2$pos),]
is.unsorted(fancm_chr5_finalpos$pos)
plot(fancm_chr5_snp2$`SNP Start`, fancm_chr5_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica fancm Chromosome 5 Genetic Map")

fancm_chr6_spl <- smooth.spline(fancm_chr6_snp2$rate, spar = .4)
fancm_chr6_snp2$pos <- (fancm_chr6_snp2$`SNP Start`*fancm_chr6_spl$y)
plot(fancm_chr6_snp2$`SNP Start`, fancm_chr6_snp2$pos)
plot(fancm_chr6_snp2$`SNP Start`, fancm_chr6_snp2$pos/fancm_chr6_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica fancm Chromosome 6 Recombination Distribution")
fancm_chr6_finalpos <- fancm_chr6_snp2[order(fancm_chr6_snp2$pos),]
is.unsorted(fancm_chr6_finalpos$pos)
plot(fancm_chr6_snp2$`SNP Start`, fancm_chr6_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica fancm Chromosome 6 Genetic Map")

fancm_chr7_spl <- smooth.spline(fancm_chr7_snp2$rate, spar = 0.4)
fancm_chr7_snp2$pos <- (fancm_chr7_snp2$`SNP Start`*fancm_chr7_spl$y)
plot(fancm_chr7_snp2$`SNP Start`, fancm_chr7_snp2$pos)
plot(fancm_chr7_snp2$`SNP Start`, fancm_chr7_snp2$pos/fancm_chr7_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica fancm Chromosome 7 Recombination Distribution")
fancm_chr7_finalpos <- fancm_chr7_snp2[order(fancm_chr7_snp2$pos),]
is.unsorted(fancm_chr7_finalpos$pos)
plot(fancm_chr7_snp2$`SNP Start`, fancm_chr7_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica fancm Chromosome 7 Genetic Map")

fancm_chr8_spl <- smooth.spline(fancm_chr8_snp2$rate, spar = .4)
fancm_chr8_snp2$pos <- (fancm_chr8_snp2$`SNP Start`*fancm_chr8_spl$y)
plot(fancm_chr8_snp2$`SNP Start`, fancm_chr8_snp2$pos)
plot(fancm_chr8_snp2$`SNP Start`, fancm_chr8_snp2$pos/fancm_chr8_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica fancm Chromosome 8 Recombination Distribution")
fancm_chr8_finalpos <- fancm_chr8_snp2[order(fancm_chr8_snp2$pos),]
is.unsorted(fancm_chr8_finalpos$pos)
plot(fancm_chr8_snp2$`SNP Start`, fancm_chr8_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica fancm Chromosome 8 Genetic Map")

fancm_chr9_spl <- smooth.spline(fancm_chr9_snp2$rate, spar = .4)
fancm_chr9_snp2$pos <- (fancm_chr9_snp2$`SNP Start`*fancm_chr9_spl$y)
plot(fancm_chr9_snp2$`SNP Start`, fancm_chr9_snp2$pos)
plot(fancm_chr9_snp2$`SNP Start`, fancm_chr9_snp2$pos/fancm_chr9_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica fancm Chromosome 9 Recombination Distribution")
fancm_chr9_finalpos <- fancm_chr9_snp2[order(fancm_chr9_snp2$pos),]
is.unsorted(fancm_chr9_finalpos$pos)
fancm_chr9_finalpos$pos <- fancm_chr9_finalpos$pos + abs(min(fancm_chr9_finalpos$pos))
plot(fancm_chr9_snp2$`SNP Start`, fancm_chr9_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica fancm Chromosome 9 Genetic Map")

fancm_chr10_spl <- smooth.spline(fancm_chr10_snp2$rate, spar =.4)
fancm_chr10_snp2$pos <- (fancm_chr10_snp2$`SNP Start`*fancm_chr10_spl$y)
plot(fancm_chr10_snp2$`SNP Start`, fancm_chr10_snp2$pos)
plot(fancm_chr10_snp2$`SNP Start`, fancm_chr10_snp2$pos/fancm_chr10_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica fancm Chromosome 10 Recombination Distribution")
fancm_chr10_finalpos <- fancm_chr10_snp2[order(fancm_chr10_snp2$pos),]
is.unsorted(fancm_chr10_finalpos$pos)
fancm_chr10_finalpos$pos <- fancm_chr10_finalpos$pos + abs(min(fancm_chr10_finalpos$pos))
plot(fancm_chr10_snp2$`SNP Start`, fancm_chr10_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica fancm Chromosome 10 Genetic Map")

fancm_chr11_spl <- smooth.spline(fancm_chr11_snp2$rate, spar = .4)
fancm_chr11_snp2$pos <- (fancm_chr11_snp2$`SNP Start`*fancm_chr11_spl$y)
plot(fancm_chr11_snp2$`SNP Start`, fancm_chr11_snp2$pos)
plot(fancm_chr11_snp2$`SNP Start`, fancm_chr11_snp2$pos/fancm_chr11_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica fancm Chromosome 11 Recombination Distribution")
fancm_chr11_finalpos <- fancm_chr11_snp2[order(fancm_chr11_snp2$pos),]
is.unsorted(fancm_chr11_finalpos$pos)
plot(fancm_chr11_snp2$`SNP Start`, fancm_chr11_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica fancm Chromosome 11 Genetic Map")

fancm_chr12_spl <- smooth.spline(fancm_chr12_snp2$rate, spar = .4)
fancm_chr12_snp2$pos <- (fancm_chr12_snp2$`SNP Start`*fancm_chr12_spl$y)
plot(fancm_chr12_snp2$`SNP Start`, fancm_chr12_snp2$pos)
plot(fancm_chr12_snp2$`SNP Start`, fancm_chr12_snp2$pos/fancm_chr12_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica fancm Chromosome 12 Recombination Distribution")
fancm_chr12_finalpos <- fancm_chr12_snp2[order(fancm_chr12_snp2$pos),]
is.unsorted(fancm_chr12_finalpos$pos)
plot(fancm_chr12_snp2$`SNP Start`, fancm_chr12_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica fancm Chromosome 12 Genetic Map")

#Final genetic map
fancm_chr1 <- fancm_chr1_finalpos$pos/100
fancm_chr1len <- length(fancm_chr1)
dim(fancm_chr1) <- c(fancm_chr1len,1)
fancm_chr1 <- list(fancm_chr1)

fancm_chr2 <- fancm_chr2_finalpos$pos/100
fancm_chr2len <- length(fancm_chr2)
dim(fancm_chr2) <- c(fancm_chr2len,1)
fancm_chr2 <- list(fancm_chr2)

fancm_chr3 <- fancm_chr3_finalpos$pos/100
fancm_chr3len <- length(fancm_chr3)
dim(fancm_chr3) <- c(fancm_chr3len,1)
fancm_chr3 <- list(fancm_chr3)

fancm_chr4 <- fancm_chr4_finalpos$pos/100
fancm_chr4len <- length(fancm_chr4)
dim(fancm_chr4) <- c(fancm_chr4len,1)
fancm_chr4 <- list(fancm_chr4)

fancm_chr5 <- fancm_chr5_finalpos$pos/100
fancm_chr5len <- length(fancm_chr5)
dim(fancm_chr5) <- c(fancm_chr5len,1)
fancm_chr5 <- list(fancm_chr5)

fancm_chr5 <- fancm_chr5_finalpos$pos/100
fancm_chr5len <- length(fancm_chr5)
dim(fancm_chr5) <- c(fancm_chr5len,1)
fancm_chr5 <- list(fancm_chr5)

fancm_chr7 <- fancm_chr7_finalpos$pos/100
fancm_chr7len <- length(fancm_chr7)
dim(fancm_chr7) <- c(fancm_chr7len,1)
fancm_chr7 <- list(fancm_chr7)

fancm_chr8 <- fancm_chr8_finalpos$pos/100
fancm_chr8len <- length(fancm_chr8)
dim(fancm_chr8) <- c(fancm_chr8len,1)
fancm_chr8 <- list(fancm_chr8)

fancm_chr9 <- fancm_chr9_finalpos$pos/100
fancm_chr9len <- length(fancm_chr9)
dim(fancm_chr9) <- c(fancm_chr9len,1)
fancm_chr9 <- list(fancm_chr9)

fancm_chr10 <- fancm_chr10_finalpos$pos/100
fancm_chr10len <- length(fancm_chr10)
dim(fancm_chr10) <- c(fancm_chr10len,1)
fancm_chr10 <- list(fancm_chr10)

fancm_chr11 <- fancm_chr11_finalpos$pos/100
fancm_chr11len <- length(fancm_chr11)
dim(fancm_chr11) <- c(fancm_chr11len,1)
fancm_chr11 <- list(fancm_chr11)

fancm_chr12 <- fancm_chr12_finalpos$pos/100
fancm_chr12len <- length(fancm_chr12)
dim(fancm_chr12) <- c(fancm_chr12len,1)
fancm_chr12 <- list(fancm_chr12)

fancm_final_map <- list(fancm_chr1[[1]], fancm_chr2[[1]], 
                         fancm_chr3[[1]], fancm_chr4[[1]], fancm_chr5[[1]], 
                         fancm_chr5[[1]], fancm_chr7[[1]], fancm_chr8[[1]], 
                         fancm_chr9[[1]], fancm_chr10[[1]],fancm_chr11[[1]], fancm_chr12[[1]])

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
c1 <-find_centromere(16.7,fancm_chr1_finalpos)
c2 <-find_centromere(13.6,fancm_chr2_finalpos)
c3 <-find_centromere(19.4,fancm_chr3_finalpos)
c4 <-find_centromere(9.7,fancm_chr4_finalpos)
c5 <-find_centromere(12.4,fancm_chr5_finalpos)
c6 <-find_centromere(15.3,fancm_chr6_finalpos)
c7 <-find_centromere(12.1,fancm_chr7_finalpos)
c8 <-find_centromere(12.9,fancm_chr8_finalpos)
c9 <-find_centromere(2.8,fancm_chr9_finalpos)
c10 <-find_centromere(8.2,fancm_chr10_finalpos)
c11 <-find_centromere(12,fancm_chr11_finalpos)
c12 <-find_centromere(11.9,fancm_chr12_finalpos)

fancm_centromere <- c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12)
fancm_centromere <- fancm_centromere/100


