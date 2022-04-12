library(AlphaSimR)
library(Rcpp)
library(ggplot2)
library(dplyr)

set.seed(420)

japonica_snps <- read.table("japonica_SNPs.bed", header =FALSE)
colnames(japonica_snps) <- c("Chr#", "SNP Start", "SNP End")
zmet2_snps <- sample_n(japonica_snps, 4000)
zmet2_snps <- zmet2_snps[order(zmet2_snps$`Chr#`,zmet2_snps$`SNP Start`),]

#splitting zmet2 snps
zmet2_chr1_snp <- zmet2_snps[ which(zmet2_snps$`Chr#` == "Chr1"),]
zmet2_chr1_snp$rate <- NA
zmet2_chr1_snp$`SNP End` <- zmet2_chr1_snp$`SNP End` - min(zmet2_chr1_snp$`SNP Start`)
zmet2_chr1_snp$`SNP Start` <- zmet2_chr1_snp$`SNP Start`- min(zmet2_chr1_snp$`SNP Start`)

zmet2_chr2_snp <- zmet2_snps[ which(zmet2_snps$`Chr#` == "Chr2"),]
zmet2_chr2_snp$rate <- NA
zmet2_chr2_snp$`SNP End` <- zmet2_chr2_snp$`SNP End` - min(zmet2_chr2_snp$`SNP Start`)
zmet2_chr2_snp$`SNP Start` <- zmet2_chr2_snp$`SNP Start`- min(zmet2_chr2_snp$`SNP Start`)

zmet2_chr3_snp <- zmet2_snps[ which(zmet2_snps$`Chr#` == "Chr3"),]
zmet2_chr3_snp$rate <- NA
zmet2_chr3_snp$`SNP End` <- zmet2_chr3_snp$`SNP End` - min(zmet2_chr3_snp$`SNP Start`)
zmet2_chr3_snp$`SNP Start` <- zmet2_chr3_snp$`SNP Start`- min(zmet2_chr3_snp$`SNP Start`)

zmet2_chr4_snp <- zmet2_snps[ which(zmet2_snps$`Chr#` == "Chr4"),]
zmet2_chr4_snp$rate <- NA
zmet2_chr4_snp$`SNP End` <- zmet2_chr4_snp$`SNP End` - min(zmet2_chr4_snp$`SNP Start`)
zmet2_chr4_snp$`SNP Start` <- zmet2_chr4_snp$`SNP Start`- min(zmet2_chr4_snp$`SNP Start`)

zmet2_chr5_snp <- zmet2_snps[ which(zmet2_snps$`Chr#` == "Chr5"),]
zmet2_chr5_snp$rate <- NA
zmet2_chr5_snp$`SNP End` <- zmet2_chr5_snp$`SNP End` - min(zmet2_chr5_snp$`SNP Start`)
zmet2_chr5_snp$`SNP Start` <- zmet2_chr5_snp$`SNP Start`- min(zmet2_chr5_snp$`SNP Start`)

zmet2_chr6_snp <- zmet2_snps[ which(zmet2_snps$`Chr#` == "Chr6"),]
zmet2_chr6_snp$rate <- NA
zmet2_chr6_snp$`SNP End` <- zmet2_chr6_snp$`SNP End` - min(zmet2_chr6_snp$`SNP Start`)
zmet2_chr6_snp$`SNP Start` <- zmet2_chr6_snp$`SNP Start`- min(zmet2_chr6_snp$`SNP Start`)

zmet2_chr7_snp <- zmet2_snps[ which(zmet2_snps$`Chr#` == "Chr7"),]
zmet2_chr7_snp$rate <- NA
zmet2_chr7_snp$`SNP End` <- zmet2_chr7_snp$`SNP End` - min(zmet2_chr7_snp$`SNP Start`)
zmet2_chr7_snp$`SNP Start` <- zmet2_chr7_snp$`SNP Start`- min(zmet2_chr7_snp$`SNP Start`)

zmet2_chr8_snp <- zmet2_snps[ which(zmet2_snps$`Chr#` == "Chr8"),]
zmet2_chr8_snp$rate <- NA
zmet2_chr8_snp$`SNP End` <- zmet2_chr8_snp$`SNP End` - min(zmet2_chr8_snp$`SNP Start`)
zmet2_chr8_snp$`SNP Start` <- zmet2_chr8_snp$`SNP Start`- min(zmet2_chr8_snp$`SNP Start`)

zmet2_chr9_snp <- zmet2_snps[ which(zmet2_snps$`Chr#` == "Chr9"),]
zmet2_chr9_snp$rate <- NA
zmet2_chr9_snp$`SNP End` <- zmet2_chr9_snp$`SNP End` - min(zmet2_chr9_snp$`SNP Start`)
zmet2_chr9_snp$`SNP Start` <- zmet2_chr9_snp$`SNP Start`- min(zmet2_chr9_snp$`SNP Start`)

zmet2_chr10_snp <- zmet2_snps[ which(zmet2_snps$`Chr#` == "Chr10"),]
zmet2_chr10_snp$rate <- NA
zmet2_chr10_snp$`SNP End` <- zmet2_chr10_snp$`SNP End` - min(zmet2_chr10_snp$`SNP Start`)
zmet2_chr10_snp$`SNP Start` <- zmet2_chr10_snp$`SNP Start`- min(zmet2_chr10_snp$`SNP Start`)

zmet2_chr11_snp <- zmet2_snps[ which(zmet2_snps$`Chr#` == "Chr11"),]
zmet2_chr11_snp$rate <- NA
zmet2_chr11_snp$`SNP End` <- zmet2_chr11_snp$`SNP End` - min(zmet2_chr11_snp$`SNP Start`)
zmet2_chr11_snp$`SNP Start` <- zmet2_chr11_snp$`SNP Start`- min(zmet2_chr11_snp$`SNP Start`)

zmet2_chr12_snp <- zmet2_snps[ which(zmet2_snps$`Chr#` == "Chr12"),]
zmet2_chr12_snp$rate <- NA
zmet2_chr12_snp$`SNP End` <- zmet2_chr12_snp$`SNP End` - min(zmet2_chr12_snp$`SNP Start`)
zmet2_chr12_snp$`SNP Start` <- zmet2_chr12_snp$`SNP Start`- min(zmet2_chr12_snp$`SNP Start`)

## multiplying WT recombination rates by the avg difference
# exclude pericentromeric regions (suppresion regions)
# 1. create avg diff column (supression region = 0)
# 2. loop through, multiply wt rate from other paper (fine scale recombination rate) by avg rate

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


##Using zmet2 recombination landscape--> 20% increase in COs
avg_diff <-1.2 #in telomeres? should I apply this generally or ... specify telomeric regions?

## Multiply recombination fine scale data by the avg rate
new_rates <- function(old_rate){
  for(i in 1:nrow(old_rate)){
    old_rate$rate[i] <-(avg_diff*old_rate$rate[i]) + old_rate$rate[i]
  }
  print(old_rate)
}
zmet2_chr1_CO_2<-new_rates(WTJap_chr1_CO)
zmet2_chr2_CO_2<-new_rates(WTJap_chr2_CO)
zmet2_chr3_CO_2<-new_rates(WTJap_chr3_CO)
zmet2_chr4_CO_2<-new_rates(WTJap_chr4_CO)
zmet2_chr5_CO_2<-new_rates(WTJap_chr5_CO)
zmet2_chr6_CO_2<-new_rates(WTJap_chr6_CO)
zmet2_chr7_CO_2<-new_rates(WTJap_chr7_CO)
zmet2_chr8_CO_2<-new_rates(WTJap_chr8_CO)
zmet2_chr9_CO_2<-new_rates(WTJap_chr9_CO)
zmet2_chr10_CO_2<-new_rates(WTJap_chr10_CO)
zmet2_chr11_CO_2<-new_rates(WTJap_chr11_CO)
zmet2_chr12_CO_2<-new_rates(WTJap_chr12_CO)

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
zmet2_chr1_CO_3 <- zmet2_chr1_CO_2
bins<-as.integer(nrow(zmet2_chr1_CO_2)/44)
zmet2_chr1_CO_3$rates<- rollapply(zmet2_chr1_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
zmet2_chr1_CO_3<-fill_start(zmet2_chr1_CO_3)
zmet2_chr1_CO_3<- zmet2_chr1_CO_3 %>% drop_na(rates)

zmet2_chr2_CO_3 <- zmet2_chr2_CO_2
bins<-as.integer(nrow(zmet2_chr2_CO_2)/40)
zmet2_chr2_CO_3$rates<- rollapply(zmet2_chr2_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
zmet2_chr2_CO_3<-fill_start(zmet2_chr2_CO_3)
zmet2_chr2_CO_3<- zmet2_chr2_CO_3 %>% drop_na(rates)

zmet2_chr3_CO_3 <- zmet2_chr3_CO_2
bins<-as.integer(nrow(zmet2_chr3_CO_2)/41)
zmet2_chr3_CO_3$rates<- rollapply(zmet2_chr3_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
zmet2_chr3_CO_3<-fill_start(zmet2_chr3_CO_3)
zmet2_chr3_CO_3<- zmet2_chr3_CO_3 %>% drop_na(rates)

zmet2_chr4_CO_3 <- zmet2_chr4_CO_2
bins<-as.integer(nrow(zmet2_chr4_CO_2)/39)
zmet2_chr4_CO_3$rates<- rollapply(zmet2_chr4_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
zmet2_chr4_CO_3<-fill_start(zmet2_chr4_CO_3)
zmet2_chr4_CO_3<- zmet2_chr4_CO_3 %>% drop_na(rates)

zmet2_chr5_CO_3 <- zmet2_chr5_CO_2
bins<-as.integer(nrow(zmet2_chr5_CO_2)/33)
zmet2_chr5_CO_3$rates<- rollapply(zmet2_chr5_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
zmet2_chr5_CO_3<-fill_start(zmet2_chr5_CO_3)
zmet2_chr5_CO_3<- zmet2_chr5_CO_3 %>% drop_na(rates)

zmet2_chr6_CO_3 <- zmet2_chr6_CO_2
bins<-as.integer(nrow(zmet2_chr6_CO_2)/32)
zmet2_chr6_CO_3$rates<- rollapply(zmet2_chr6_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
zmet2_chr6_CO_3<-fill_start(zmet2_chr6_CO_3)
zmet2_chr6_CO_3<- zmet2_chr6_CO_3 %>% drop_na(rates)

zmet2_chr7_CO_3 <- zmet2_chr7_CO_2
bins<-as.integer(nrow(zmet2_chr7_CO_2)/35)
zmet2_chr7_CO_3$rates<- rollapply(zmet2_chr7_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
zmet2_chr7_CO_3<-fill_start(zmet2_chr7_CO_3)
zmet2_chr7_CO_3<- zmet2_chr7_CO_3 %>% drop_na(rates)

zmet2_chr8_CO_3 <- zmet2_chr8_CO_2
bins<-as.integer(nrow(zmet2_chr8_CO_2)/28)
zmet2_chr8_CO_3$rates<- rollapply(zmet2_chr8_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
zmet2_chr8_CO_3<-fill_start(zmet2_chr8_CO_3)
zmet2_chr8_CO_3<- zmet2_chr8_CO_3 %>% drop_na(rates)

zmet2_chr9_CO_3 <- zmet2_chr9_CO_2
bins<-as.integer(nrow(zmet2_chr9_CO_2)/22)
zmet2_chr9_CO_3$rates<- rollapply(zmet2_chr9_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
zmet2_chr9_CO_3<-fill_start(zmet2_chr9_CO_3)
zmet2_chr9_CO_3<- zmet2_chr9_CO_3 %>% drop_na(rates)

zmet2_chr10_CO_3 <- zmet2_chr10_CO_2
bins<-as.integer(nrow(zmet2_chr10_CO_2)/27)
zmet2_chr10_CO_3$rates<- rollapply(zmet2_chr10_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
zmet2_chr10_CO_3<-fill_start(zmet2_chr10_CO_3)
zmet2_chr10_CO_3<- zmet2_chr10_CO_3 %>% drop_na(rates)

zmet2_chr11_CO_3 <- zmet2_chr11_CO_2
bins<-as.integer(nrow(zmet2_chr11_CO_2)/30)
zmet2_chr11_CO_3$rates<- rollapply(zmet2_chr11_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
zmet2_chr11_CO_3<-fill_start(zmet2_chr11_CO_3)
zmet2_chr11_CO_3<- zmet2_chr11_CO_3 %>% drop_na(rates)

zmet2_chr12_CO_3 <- zmet2_chr12_CO_2
bins<-as.integer(nrow(zmet2_chr12_CO_2)/31)
zmet2_chr12_CO_3$rates<- rollapply(zmet2_chr12_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
zmet2_chr12_CO_3<-fill_start(zmet2_chr12_CO_3)
zmet2_chr12_CO_3<- zmet2_chr12_CO_3 %>% drop_na(rates)

## assigning frequency to SNPs based on recombination frequency in each bin
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


#using function,  get cM/Mb for final genetic position - assign rates
zmet2_chr1_snp2 <- snp_rate(zmet2_chr1_CO_3, zmet2_chr1_snp)
zmet2_chr2_snp2 <- snp_rate(zmet2_chr2_CO_3, zmet2_chr2_snp)
zmet2_chr3_snp2 <- snp_rate(zmet2_chr3_CO_3, zmet2_chr3_snp)
zmet2_chr4_snp2 <- snp_rate(zmet2_chr4_CO_3, zmet2_chr4_snp)
zmet2_chr5_snp2 <- snp_rate(zmet2_chr5_CO_3, zmet2_chr5_snp)
zmet2_chr6_snp2 <- snp_rate(zmet2_chr6_CO_3, zmet2_chr6_snp)
zmet2_chr7_snp2 <- snp_rate(zmet2_chr7_CO_3, zmet2_chr7_snp)
zmet2_chr8_snp2 <- snp_rate(zmet2_chr8_CO_3, zmet2_chr8_snp)
zmet2_chr9_snp2 <- snp_rate(zmet2_chr9_CO_3, zmet2_chr9_snp)
zmet2_chr10_snp2 <- snp_rate(zmet2_chr10_CO_3, zmet2_chr10_snp)
zmet2_chr11_snp2 <- snp_rate(zmet2_chr11_CO_3, zmet2_chr11_snp)
zmet2_chr12_snp2 <- snp_rate(zmet2_chr12_CO_3, zmet2_chr12_snp)

#converted SNP start to Mb
zmet2_chr1_snp2$`SNP Start` <- zmet2_chr1_snp2$`SNP Start`/1000000
zmet2_chr2_snp2$`SNP Start` <- zmet2_chr2_snp2$`SNP Start`/1000000
zmet2_chr3_snp2$`SNP Start` <- zmet2_chr3_snp2$`SNP Start`/1000000
zmet2_chr4_snp2$`SNP Start` <- zmet2_chr4_snp2$`SNP Start`/1000000
zmet2_chr5_snp2$`SNP Start` <- zmet2_chr5_snp2$`SNP Start`/1000000
zmet2_chr6_snp2$`SNP Start` <- zmet2_chr6_snp2$`SNP Start`/1000000
zmet2_chr7_snp2$`SNP Start` <- zmet2_chr7_snp2$`SNP Start`/1000000
zmet2_chr8_snp2$`SNP Start` <- zmet2_chr8_snp2$`SNP Start`/1000000
zmet2_chr9_snp2$`SNP Start` <- zmet2_chr9_snp2$`SNP Start`/1000000
zmet2_chr10_snp2$`SNP Start` <- zmet2_chr10_snp2$`SNP Start`/1000000
zmet2_chr11_snp2$`SNP Start` <- zmet2_chr11_snp2$`SNP Start`/1000000
zmet2_chr12_snp2$`SNP Start` <- zmet2_chr12_snp2$`SNP Start`/1000000

zmet2_chr1_snp2$`SNP End` <- zmet2_chr1_snp2$`SNP End`/1000000
zmet2_chr2_snp2$`SNP End` <- zmet2_chr2_snp2$`SNP End`/1000000
zmet2_chr3_snp2$`SNP End` <- zmet2_chr3_snp2$`SNP End`/1000000
zmet2_chr4_snp2$`SNP End` <- zmet2_chr4_snp2$`SNP End`/1000000
zmet2_chr5_snp2$`SNP End` <- zmet2_chr5_snp2$`SNP End`/1000000
zmet2_chr6_snp2$`SNP End` <- zmet2_chr6_snp2$`SNP End`/1000000
zmet2_chr7_snp2$`SNP End` <- zmet2_chr7_snp2$`SNP End`/1000000
zmet2_chr8_snp2$`SNP End` <- zmet2_chr8_snp2$`SNP End`/1000000
zmet2_chr9_snp2$`SNP End` <- zmet2_chr9_snp2$`SNP End`/1000000
zmet2_chr10_snp2$`SNP End` <- zmet2_chr10_snp2$`SNP End`/1000000
zmet2_chr11_snp2$`SNP End` <- zmet2_chr11_snp2$`SNP End`/1000000
zmet2_chr12_snp2$`SNP End` <- zmet2_chr12_snp2$`SNP End`/1000000

#omit empty col
zmet2_chr1_snp2<-na.omit(zmet2_chr1_snp2)
zmet2_chr2_snp2<-na.omit(zmet2_chr2_snp2)
zmet2_chr3_snp2<-na.omit(zmet2_chr3_snp2)
zmet2_chr4_snp2<-na.omit(zmet2_chr4_snp2)
zmet2_chr5_snp2<-na.omit(zmet2_chr5_snp2)
zmet2_chr6_snp2<-na.omit(zmet2_chr6_snp2)
zmet2_chr7_snp2<-na.omit(zmet2_chr7_snp2)
zmet2_chr8_snp2<-na.omit(zmet2_chr8_snp2)
zmet2_chr9_snp2<-na.omit(zmet2_chr9_snp2)
zmet2_chr10_snp2<-na.omit(zmet2_chr10_snp2)
zmet2_chr11_snp2<-na.omit(zmet2_chr11_snp2)
zmet2_chr12_snp2<-na.omit(zmet2_chr12_snp2)

#gen maps
zmet2_chr1_snp2 <- zmet2_chr1_snp2[order(zmet2_chr1_snp2$`SNP Start`),]
zmet2_chr1_spl <- smooth.spline(zmet2_chr1_snp2$rate, spar = .7)
zmet2_chr1_snp2$pos <- (zmet2_chr1_snp2$`SNP Start`*zmet2_chr1_spl$y)
plot(zmet2_chr1_snp2$`SNP Start`, zmet2_chr1_snp2$pos)
ggplot(zmet2_chr1_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
plot(zmet2_chr1_snp2$`SNP Start`, zmet2_chr1_snp2$pos/zmet2_chr1_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica zmet2 Chromosome 1 Recombination Distribution")
zmet2_chr1_finalpos <- zmet2_chr1_snp2[order(zmet2_chr1_snp2$pos),]
is.unsorted(zmet2_chr1_finalpos$pos)
plot(zmet2_chr1_snp2$`SNP Start`, zmet2_chr1_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica zmet2 Chromosome 1 Genetic Map")
plot(zmet2_chr1_finalpos$`SNP Start`, zmet2_chr1_finalpos$pos)


zmet2_chr2_spl <- smooth.spline(zmet2_chr2_snp2$rate, spar = .67)
zmet2_chr2_snp2$pos <- (zmet2_chr2_snp2$`SNP Start`*zmet2_chr2_spl$y)
plot(zmet2_chr2_snp2$`SNP Start`, zmet2_chr2_snp2$pos)
plot(zmet2_chr2_snp2$`SNP Start`, zmet2_chr2_snp2$pos/zmet2_chr2_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica zmet2 Chromosome 2 Recombination Distribution")
zmet2_chr2_finalpos <- zmet2_chr2_snp2[order(zmet2_chr2_snp2$pos),]
is.unsorted(zmet2_chr2_finalpos$pos)
plot(zmet2_chr2_snp2$`SNP Start`, zmet2_chr2_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica zmet2 Chromosome 2 Genetic Map")

zmet2_chr3_spl <- smooth.spline(zmet2_chr3_snp2$rate, spar = .6)
zmet2_chr3_snp2$pos <- (zmet2_chr3_snp2$`SNP Start`*zmet2_chr3_spl$y)
plot(zmet2_chr3_snp2$`SNP Start`, zmet2_chr3_snp2$pos)
plot(zmet2_chr3_snp2$`SNP Start`, zmet2_chr3_snp2$pos/zmet2_chr3_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica zmet2 Chromosome 3 Recombination Distribution")
zmet2_chr3_finalpos <- zmet2_chr3_snp2[order(zmet2_chr3_snp2$pos),]
is.unsorted(zmet2_chr3_finalpos$pos)
plot(zmet2_chr3_snp2$`SNP Start`, zmet2_chr3_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica zmet2 Chromosome 3 Genetic Map")

zmet2_chr4_spl <- smooth.spline(zmet2_chr4_snp2$rate, spar = .75)
zmet2_chr4_snp2$pos <- (zmet2_chr4_snp2$`SNP Start`*zmet2_chr4_spl$y)
plot(zmet2_chr4_snp2$`SNP Start`, zmet2_chr4_snp2$pos)
plot(zmet2_chr4_snp2$`SNP Start`, zmet2_chr4_snp2$pos/zmet2_chr4_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica zmet2 Chromosome 4 Recombination Distribution")
zmet2_chr4_finalpos <- zmet2_chr4_snp2[order(zmet2_chr4_snp2$pos),]
is.unsorted(zmet2_chr4_finalpos$pos)
plot(zmet2_chr4_snp2$`SNP Start`, zmet2_chr4_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica zmet2 Chromosome 4 Genetic Map")

zmet2_chr5_spl <- smooth.spline(zmet2_chr5_snp2$rate, spar =.66)
zmet2_chr5_snp2$pos <- (zmet2_chr5_snp2$`SNP Start`*zmet2_chr5_spl$y)
plot(zmet2_chr5_snp2$`SNP Start`, zmet2_chr5_snp2$pos)
plot(zmet2_chr5_snp2$`SNP Start`, zmet2_chr5_snp2$pos/zmet2_chr5_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica zmet2 Chromosome 5 Recombination Distribution")
zmet2_chr5_finalpos <- zmet2_chr5_snp2[order(zmet2_chr5_snp2$pos),]
is.unsorted(zmet2_chr5_finalpos$pos)
plot(zmet2_chr5_snp2$`SNP Start`, zmet2_chr5_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica zmet2 Chromosome 5 Genetic Map")

zmet2_chr6_spl <- smooth.spline(zmet2_chr6_snp2$rate, spar = .745)
zmet2_chr6_snp2$pos <- (zmet2_chr6_snp2$`SNP Start`*zmet2_chr6_spl$y)
plot(zmet2_chr6_snp2$`SNP Start`, zmet2_chr6_snp2$pos)
plot(zmet2_chr6_snp2$`SNP Start`, zmet2_chr6_snp2$pos/zmet2_chr6_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica zmet2 Chromosome 6 Recombination Distribution")
zmet2_chr6_finalpos <- zmet2_chr6_snp2[order(zmet2_chr6_snp2$pos),]
is.unsorted(zmet2_chr6_finalpos$pos)
plot(zmet2_chr6_snp2$`SNP Start`, zmet2_chr6_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica zmet2 Chromosome 6 Genetic Map")

zmet2_chr7_spl <- smooth.spline(zmet2_chr7_snp2$rate, spar = 0.7)
zmet2_chr7_snp2$pos <- (zmet2_chr7_snp2$`SNP Start`*zmet2_chr7_spl$y)
plot(zmet2_chr7_snp2$`SNP Start`, zmet2_chr7_snp2$pos)
plot(zmet2_chr7_snp2$`SNP Start`, zmet2_chr7_snp2$pos/zmet2_chr7_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica zmet2 Chromosome 7 Recombination Distribution")
zmet2_chr7_finalpos <- zmet2_chr7_snp2[order(zmet2_chr7_snp2$pos),]
is.unsorted(zmet2_chr7_finalpos$pos)
plot(zmet2_chr7_snp2$`SNP Start`, zmet2_chr7_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica zmet2 Chromosome 7 Genetic Map")

zmet2_chr8_spl <- smooth.spline(zmet2_chr8_snp2$rate, spar = 1)
zmet2_chr8_snp2$pos <- (zmet2_chr8_snp2$`SNP Start`*zmet2_chr8_spl$y)
plot(zmet2_chr8_snp2$`SNP Start`, zmet2_chr8_snp2$pos)
plot(zmet2_chr8_snp2$`SNP Start`, zmet2_chr8_snp2$pos/zmet2_chr8_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica zmet2 Chromosome 8 Recombination Distribution")
zmet2_chr8_finalpos <- zmet2_chr8_snp2[order(zmet2_chr8_snp2$pos),]
is.unsorted(zmet2_chr8_finalpos$pos)
plot(zmet2_chr8_snp2$`SNP Start`, zmet2_chr8_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica zmet2 Chromosome 8 Genetic Map")

zmet2_chr9_spl <- smooth.spline(zmet2_chr9_snp2$rate, spar = .9)
zmet2_chr9_snp2$pos <- (zmet2_chr9_snp2$`SNP Start`*zmet2_chr9_spl$y)
plot(zmet2_chr9_snp2$`SNP Start`, zmet2_chr9_snp2$pos)
plot(zmet2_chr9_snp2$`SNP Start`, zmet2_chr9_snp2$pos/zmet2_chr9_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica zmet2 Chromosome 9 Recombination Distribution")
zmet2_chr9_finalpos <- zmet2_chr9_snp2[order(zmet2_chr9_snp2$pos),]
is.unsorted(zmet2_chr9_finalpos$pos)
plot(zmet2_chr9_snp2$`SNP Start`, zmet2_chr9_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica zmet2 Chromosome 9 Genetic Map")

zmet2_chr10_spl <- smooth.spline(zmet2_chr10_snp2$rate, spar =.7)
zmet2_chr10_snp2$pos <- (zmet2_chr10_snp2$`SNP Start`*zmet2_chr10_spl$y)
plot(zmet2_chr10_snp2$`SNP Start`, zmet2_chr10_snp2$pos)
plot(zmet2_chr10_snp2$`SNP Start`, zmet2_chr10_snp2$pos/zmet2_chr10_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica zmet2 Chromosome 10 Recombination Distribution")
zmet2_chr10_finalpos <- zmet2_chr10_snp2[order(zmet2_chr10_snp2$pos),]
is.unsorted(zmet2_chr10_finalpos$pos)
plot(zmet2_chr10_snp2$`SNP Start`, zmet2_chr10_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica zmet2 Chromosome 10 Genetic Map")

zmet2_chr11_spl <- smooth.spline(zmet2_chr11_snp2$rate, spar = .7)
zmet2_chr11_snp2$pos <- (zmet2_chr11_snp2$`SNP Start`*zmet2_chr11_spl$y)
plot(zmet2_chr11_snp2$`SNP Start`, zmet2_chr11_snp2$pos)
plot(zmet2_chr11_snp2$`SNP Start`, zmet2_chr11_snp2$pos/zmet2_chr11_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica zmet2 Chromosome 11 Recombination Distribution")
zmet2_chr11_finalpos <- zmet2_chr11_snp2[order(zmet2_chr11_snp2$pos),]
is.unsorted(zmet2_chr11_finalpos$pos)
plot(zmet2_chr11_snp2$`SNP Start`, zmet2_chr11_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica zmet2 Chromosome 11 Genetic Map")

zmet2_chr12_spl <- smooth.spline(zmet2_chr12_snp2$rate, spar = .9)
zmet2_chr12_snp2$pos <- (zmet2_chr12_snp2$`SNP Start`*zmet2_chr12_spl$y)
plot(zmet2_chr12_snp2$`SNP Start`, zmet2_chr12_snp2$pos)
plot(zmet2_chr12_snp2$`SNP Start`, zmet2_chr12_snp2$pos/zmet2_chr12_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica zmet2 Chromosome 12 Recombination Distribution")
zmet2_chr12_finalpos <- zmet2_chr12_snp2[order(zmet2_chr12_snp2$pos),]
is.unsorted(zmet2_chr12_finalpos$pos)
plot(zmet2_chr12_snp2$`SNP Start`, zmet2_chr12_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica zmet2 Chromosome 12 Genetic Map")

#Final genetic map
zmet2_chr1 <- zmet2_chr1_finalpos$pos/100
zmet2_chr1len <- length(zmet2_chr1)
dim(zmet2_chr1) <- c(zmet2_chr1len,1)
zmet2_chr1 <- list(zmet2_chr1)

zmet2_chr2 <- zmet2_chr2_finalpos$pos/100
zmet2_chr2len <- length(zmet2_chr2)
dim(zmet2_chr2) <- c(zmet2_chr2len,1)
zmet2_chr2 <- list(zmet2_chr2)

zmet2_chr3 <- zmet2_chr3_finalpos$pos/100
zmet2_chr3len <- length(zmet2_chr3)
dim(zmet2_chr3) <- c(zmet2_chr3len,1)
zmet2_chr3 <- list(zmet2_chr3)

zmet2_chr4 <- zmet2_chr4_finalpos$pos/100
zmet2_chr4len <- length(zmet2_chr4)
dim(zmet2_chr4) <- c(zmet2_chr4len,1)
zmet2_chr4 <- list(zmet2_chr4)

zmet2_chr5 <- zmet2_chr5_finalpos$pos/100
zmet2_chr5len <- length(zmet2_chr5)
dim(zmet2_chr5) <- c(zmet2_chr5len,1)
zmet2_chr5 <- list(zmet2_chr5)

zmet2_chr5 <- zmet2_chr5_finalpos$pos/100
zmet2_chr5len <- length(zmet2_chr5)
dim(zmet2_chr5) <- c(zmet2_chr5len,1)
zmet2_chr5 <- list(zmet2_chr5)

zmet2_chr7 <- zmet2_chr7_finalpos$pos/100
zmet2_chr7len <- length(zmet2_chr7)
dim(zmet2_chr7) <- c(zmet2_chr7len,1)
zmet2_chr7 <- list(zmet2_chr7)

zmet2_chr8 <- zmet2_chr8_finalpos$pos/100
zmet2_chr8len <- length(zmet2_chr8)
dim(zmet2_chr8) <- c(zmet2_chr8len,1)
zmet2_chr8 <- list(zmet2_chr8)

zmet2_chr9 <- zmet2_chr9_finalpos$pos/100
zmet2_chr9len <- length(zmet2_chr9)
dim(zmet2_chr9) <- c(zmet2_chr9len,1)
zmet2_chr9 <- list(zmet2_chr9)

zmet2_chr10 <- zmet2_chr10_finalpos$pos/100
zmet2_chr10len <- length(zmet2_chr10)
dim(zmet2_chr10) <- c(zmet2_chr10len,1)
zmet2_chr10 <- list(zmet2_chr10)

zmet2_chr11 <- zmet2_chr11_finalpos$pos/100
zmet2_chr11len <- length(zmet2_chr11)
dim(zmet2_chr11) <- c(zmet2_chr11len,1)
zmet2_chr11 <- list(zmet2_chr11)

zmet2_chr12 <- zmet2_chr12_finalpos$pos/100
zmet2_chr12len <- length(zmet2_chr12)
dim(zmet2_chr12) <- c(zmet2_chr12len,1)
zmet2_chr12 <- list(zmet2_chr12)

zmet2_final_map <- list(zmet2_chr1[[1]], zmet2_chr2[[1]], 
                       zmet2_chr3[[1]], zmet2_chr4[[1]], zmet2_chr5[[1]], 
                       zmet2_chr5[[1]], zmet2_chr7[[1]], zmet2_chr8[[1]], 
                       zmet2_chr9[[1]], zmet2_chr10[[1]],zmet2_chr11[[1]], zmet2_chr12[[1]])

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
c1 <-find_centromere(16.7,zmet2_chr1_finalpos)
c2 <-find_centromere(13.6,zmet2_chr2_finalpos)
c3 <-find_centromere(19.4,zmet2_chr3_finalpos)
c4 <-find_centromere(9.7,zmet2_chr4_finalpos)
c5 <-find_centromere(12.4,zmet2_chr5_finalpos)
c6 <-find_centromere(15.3,zmet2_chr6_finalpos)
c7 <-find_centromere(12.1,zmet2_chr7_finalpos)
c8 <-find_centromere(12.9,zmet2_chr8_finalpos)
c9 <-find_centromere(2.8,zmet2_chr9_finalpos)
c10 <-find_centromere(8.2,zmet2_chr10_finalpos)
c11 <-find_centromere(12,zmet2_chr11_finalpos)
c12 <-find_centromere(11.9,zmet2_chr12_finalpos)

zmet2_centromere <- c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12)
zmet2_centromere <- zmet2_centromere/100


