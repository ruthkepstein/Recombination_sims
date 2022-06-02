library(AlphaSimR)
library(Rcpp)
library(ggplot2)
library(dplyr)

#Set seed to ensure sampling is the same each time
set.seed(420)

#Read in the SNP dataset from japonica subspecies genome, randomly sample 2000 SNPs, and order them by SNP start positions
japonica_snps <- read.table("japonica_SNPs.bed", header =FALSE)
colnames(japonica_snps) <- c("Chr#", "SNP Start", "SNP End")
recq4_snps <- sample_n(japonica_snps, 2000)
recq4_snps <- recq4_snps[order(recq4_snps$`Chr#`,recq4_snps$`SNP Start`),]

#Create separate dataframes for SNPs on each chromosome, starting the SNP start positions at zero
recq4_chr1_snp <- recq4_snps[ which(recq4_snps$`Chr#` == "Chr1"),]
recq4_chr1_snp$rate <- NA
recq4_chr1_snp$`SNP End` <- recq4_chr1_snp$`SNP End` - min(recq4_chr1_snp$`SNP Start`)
recq4_chr1_snp$`SNP Start` <- recq4_chr1_snp$`SNP Start`- min(recq4_chr1_snp$`SNP Start`)

recq4_chr2_snp <- recq4_snps[ which(recq4_snps$`Chr#` == "Chr2"),]
recq4_chr2_snp$rate <- NA
recq4_chr2_snp$`SNP End` <- recq4_chr2_snp$`SNP End` - min(recq4_chr2_snp$`SNP Start`)
recq4_chr2_snp$`SNP Start` <- recq4_chr2_snp$`SNP Start`- min(recq4_chr2_snp$`SNP Start`)

recq4_chr3_snp <- recq4_snps[ which(recq4_snps$`Chr#` == "Chr3"),]
recq4_chr3_snp$rate <- NA
recq4_chr3_snp$`SNP End` <- recq4_chr3_snp$`SNP End` - min(recq4_chr3_snp$`SNP Start`)
recq4_chr3_snp$`SNP Start` <- recq4_chr3_snp$`SNP Start`- min(recq4_chr3_snp$`SNP Start`)

recq4_chr4_snp <- recq4_snps[ which(recq4_snps$`Chr#` == "Chr4"),]
recq4_chr4_snp$rate <- NA
recq4_chr4_snp$`SNP End` <- recq4_chr4_snp$`SNP End` - min(recq4_chr4_snp$`SNP Start`)
recq4_chr4_snp$`SNP Start` <- recq4_chr4_snp$`SNP Start`- min(recq4_chr4_snp$`SNP Start`)

recq4_chr5_snp <- recq4_snps[ which(recq4_snps$`Chr#` == "Chr5"),]
recq4_chr5_snp$rate <- NA
recq4_chr5_snp$`SNP End` <- recq4_chr5_snp$`SNP End` - min(recq4_chr5_snp$`SNP Start`)
recq4_chr5_snp$`SNP Start` <- recq4_chr5_snp$`SNP Start`- min(recq4_chr5_snp$`SNP Start`)

recq4_chr6_snp <- recq4_snps[ which(recq4_snps$`Chr#` == "Chr6"),]
recq4_chr6_snp$rate <- NA
recq4_chr6_snp$`SNP End` <- recq4_chr6_snp$`SNP End` - min(recq4_chr6_snp$`SNP Start`)
recq4_chr6_snp$`SNP Start` <- recq4_chr6_snp$`SNP Start`- min(recq4_chr6_snp$`SNP Start`)

recq4_chr7_snp <- recq4_snps[ which(recq4_snps$`Chr#` == "Chr7"),]
recq4_chr7_snp$rate <- NA
recq4_chr7_snp$`SNP End` <- recq4_chr7_snp$`SNP End` - min(recq4_chr7_snp$`SNP Start`)
recq4_chr7_snp$`SNP Start` <- recq4_chr7_snp$`SNP Start`- min(recq4_chr7_snp$`SNP Start`)

recq4_chr8_snp <- recq4_snps[ which(recq4_snps$`Chr#` == "Chr8"),]
recq4_chr8_snp$rate <- NA
recq4_chr8_snp$`SNP End` <- recq4_chr8_snp$`SNP End` - min(recq4_chr8_snp$`SNP Start`)
recq4_chr8_snp$`SNP Start` <- recq4_chr8_snp$`SNP Start`- min(recq4_chr8_snp$`SNP Start`)

recq4_chr9_snp <- recq4_snps[ which(recq4_snps$`Chr#` == "Chr9"),]
recq4_chr9_snp$rate <- NA
recq4_chr9_snp$`SNP End` <- recq4_chr9_snp$`SNP End` - min(recq4_chr9_snp$`SNP Start`)
recq4_chr9_snp$`SNP Start` <- recq4_chr9_snp$`SNP Start`- min(recq4_chr9_snp$`SNP Start`)

recq4_chr10_snp <- recq4_snps[ which(recq4_snps$`Chr#` == "Chr10"),]
recq4_chr10_snp$rate <- NA
recq4_chr10_snp$`SNP End` <- recq4_chr10_snp$`SNP End` - min(recq4_chr10_snp$`SNP Start`)
recq4_chr10_snp$`SNP Start` <- recq4_chr10_snp$`SNP Start`- min(recq4_chr10_snp$`SNP Start`)

recq4_chr11_snp <- recq4_snps[ which(recq4_snps$`Chr#` == "Chr11"),]
recq4_chr11_snp$rate <- NA
recq4_chr11_snp$`SNP End` <- recq4_chr11_snp$`SNP End` - min(recq4_chr11_snp$`SNP Start`)
recq4_chr11_snp$`SNP Start` <- recq4_chr11_snp$`SNP Start`- min(recq4_chr11_snp$`SNP Start`)

recq4_chr12_snp <- recq4_snps[ which(recq4_snps$`Chr#` == "Chr12"),]
recq4_chr12_snp$rate <- NA
recq4_chr12_snp$`SNP End` <- recq4_chr12_snp$`SNP End` - min(recq4_chr12_snp$`SNP Start`)
recq4_chr12_snp$`SNP Start` <- recq4_chr12_snp$`SNP Start`- min(recq4_chr12_snp$`SNP Start`)

library(dlookr)
library(tidyverse)
library(OneR)

##Read in wildtype recombination rate data and create separate dataframe for each chromosome
japonica_CO <- read.csv("japonica_wt_rate.csv", header = TRUE)
colnames(japonica_CO) <- c("Chr", "CO Start", "CO End", "WT_rate")
japonica_CO <- japonica_CO[order(japonica_CO$Chr,japonica_CO$`CO Start`),]

japonica_chr1_CO <- japonica_CO[ which(japonica_CO$Chr == "1"),]
japonica_chr1_CO <- japonica_chr1_CO[order(japonica_chr1_CO$`CO Start`),]

japonica_chr2_CO <- japonica_CO[ which(japonica_CO$Chr == "2"),]
japonica_chr2_CO <- japonica_chr2_CO[order(japonica_chr2_CO$`CO Start`),]

japonica_chr3_CO <- japonica_CO[ which(japonica_CO$Chr == "3"),]
japonica_chr3_CO <- japonica_chr3_CO[order(japonica_chr3_CO$`CO Start`),]

japonica_chr4_CO <- japonica_CO[ which(japonica_CO$Chr == "4"),]
japonica_chr4_CO <- japonica_chr4_CO[order(japonica_chr4_CO$`CO Start`),]

japonica_chr5_CO <- japonica_CO[ which(japonica_CO$Chr == "5"),]
japonica_chr5_CO <- japonica_chr5_CO[order(japonica_chr5_CO$`CO Start`),]

japonica_chr6_CO <- japonica_CO[ which(japonica_CO$Chr == "6"),]
japonica_chr6_CO <- japonica_chr6_CO[order(japonica_chr6_CO$`CO Start`),]

japonica_chr7_CO <- japonica_CO[ which(japonica_CO$Chr == "7"),]
japonica_chr7_CO <- japonica_chr7_CO[order(japonica_chr7_CO$`CO Start`),]

japonica_chr8_CO <- japonica_CO[ which(japonica_CO$Chr == "8"),]
japonica_chr8_CO <- japonica_chr8_CO[order(japonica_chr8_CO$`CO Start`),]

japonica_chr9_CO <- japonica_CO[ which(japonica_CO$Chr == "9"),]
japonica_chr9_CO <- japonica_chr9_CO[order(japonica_chr9_CO$`CO Start`),]

japonica_chr10_CO <- japonica_CO[ which(japonica_CO$Chr == "10"),]
japonica_chr10_CO <- japonica_chr10_CO[order(japonica_chr10_CO$`CO Start`),]

japonica_chr11_CO <- japonica_CO[ which(japonica_CO$Chr == "11"),]
japonica_chr11_CO <- japonica_chr11_CO[order(japonica_chr11_CO$`CO Start`),]

japonica_chr12_CO <- japonica_CO[ which(japonica_CO$Chr == "12"),]
japonica_chr12_CO <- japonica_chr12_CO[order(japonica_chr12_CO$`CO Start`),]

#Making each recombination rate intervals start at zero
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

#Average difference (variable: avg_diff) is avg recombination rate difference between recq4 mutant and wildtype rates
#Function to apply avg difference to chromosome arms: divide chromosome into five segments and apply avg_diff to arms (first, second, fourth, fifth segment)
avg_diff <- 1.581055
telomeric <- function(CO){
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
japonica_chr1_CO<- telomeric(japonica_chr1_CO)
japonica_chr2_CO <- telomeric(japonica_chr2_CO)
japonica_chr3_CO <- telomeric(japonica_chr3_CO)
japonica_chr4_CO <- telomeric(japonica_chr4_CO)
japonica_chr5_CO <- telomeric(japonica_chr5_CO)
japonica_chr6_CO <- telomeric(japonica_chr6_CO)
japonica_chr7_CO <- telomeric(japonica_chr7_CO)
japonica_chr8_CO <- telomeric(japonica_chr8_CO)
japonica_chr9_CO <- telomeric(japonica_chr9_CO)
japonica_chr10_CO <- telomeric(japonica_chr10_CO)
japonica_chr11_CO <- telomeric(japonica_chr11_CO)
japonica_chr12_CO <- telomeric(japonica_chr12_CO)

#Import fine scale recombination rates and create separate dataframe for each chromosome
WTjaponica_CO <- read.table("japonica_rec_rate.bed", header = FALSE)
colnames(WTjaponica_CO) <- c("Chr", "CO Start", "CO End", "rate")
WTjaponica_CO <- WTjaponica_CO[order(WTjaponica_CO$Chr,WTjaponica_CO$`CO Start`),]

WTjaponica_chr1_CO <- WTjaponica_CO[ which(WTjaponica_CO$Chr == "chr01"),]
WTjaponica_chr1_CO <- WTjaponica_chr1_CO[order(WTjaponica_chr1_CO$`CO Start`),]

WTjaponica_chr2_CO <- WTjaponica_CO[ which(WTjaponica_CO$Chr == "chr02"),]
WTjaponica_chr2_CO <- WTjaponica_chr2_CO[order(WTjaponica_chr2_CO$`CO Start`),]

WTjaponica_chr3_CO <- WTjaponica_CO[ which(WTjaponica_CO$Chr == "chr03"),]
WTjaponica_chr3_CO <- WTjaponica_chr3_CO[order(WTjaponica_chr3_CO$`CO Start`),]

WTjaponica_chr4_CO <- WTjaponica_CO[ which(WTjaponica_CO$Chr == "chr04"),]
WTjaponica_chr4_CO <- WTjaponica_chr4_CO[order(WTjaponica_chr4_CO$`CO Start`),]

WTjaponica_chr5_CO <- WTjaponica_CO[ which(WTjaponica_CO$Chr == "chr05"),]
WTjaponica_chr5_CO <- WTjaponica_chr5_CO[order(WTjaponica_chr5_CO$`CO Start`),]

WTjaponica_chr6_CO <- WTjaponica_CO[ which(WTjaponica_CO$Chr == "chr06"),]
WTjaponica_chr6_CO <- WTjaponica_chr6_CO[order(WTjaponica_chr6_CO$`CO Start`),]

WTjaponica_chr7_CO <- WTjaponica_CO[ which(WTjaponica_CO$Chr == "chr07"),]
WTjaponica_chr7_CO <- WTjaponica_chr7_CO[order(WTjaponica_chr7_CO$`CO Start`),]

WTjaponica_chr8_CO <- WTjaponica_CO[ which(WTjaponica_CO$Chr == "chr08"),]
WTjaponica_chr8_CO <- WTjaponica_chr8_CO[order(WTjaponica_chr8_CO$`CO Start`),]

WTjaponica_chr9_CO <- WTjaponica_CO[ which(WTjaponica_CO$Chr == "chr09"),]
WTjaponica_chr9_CO <- WTjaponica_chr9_CO[order(WTjaponica_chr9_CO$`CO Start`),]

WTjaponica_chr10_CO <- WTjaponica_CO[ which(WTjaponica_CO$Chr == "chr10"),]
WTjaponica_chr10_CO <- WTjaponica_chr10_CO[order(WTjaponica_chr10_CO$`CO Start`),]

WTjaponica_chr11_CO <- WTjaponica_CO[ which(WTjaponica_CO$Chr == "chr11"),]
WTjaponica_chr11_CO <- WTjaponica_chr11_CO[order(WTjaponica_chr11_CO$`CO Start`),]

WTjaponica_chr12_CO <- WTjaponica_CO[ which(WTjaponica_CO$Chr == "chr12"),]
WTjaponica_chr12_CO <- WTjaponica_chr12_CO[order(WTjaponica_chr12_CO$`CO Start`),]

#Make intervals start at 0
WTjaponica_chr1_CO$`CO End` <- WTjaponica_chr1_CO$`CO End` - min(WTjaponica_chr1_CO$`CO Start`)
WTjaponica_chr1_CO$`CO Start` <- WTjaponica_chr1_CO$`CO Start` - min(WTjaponica_chr1_CO$`CO Start`)

WTjaponica_chr2_CO$`CO End` <- WTjaponica_chr2_CO$`CO End` - min(WTjaponica_chr2_CO$`CO Start`)
WTjaponica_chr2_CO$`CO Start` <- WTjaponica_chr2_CO$`CO Start` - min(WTjaponica_chr2_CO$`CO Start`)

WTjaponica_chr3_CO$`CO End` <- WTjaponica_chr3_CO$`CO End` - min(WTjaponica_chr3_CO$`CO Start`)
WTjaponica_chr3_CO$`CO Start` <- WTjaponica_chr3_CO$`CO Start` - min(WTjaponica_chr3_CO$`CO Start`)

WTjaponica_chr4_CO$`CO End` <- WTjaponica_chr4_CO$`CO End` - min(WTjaponica_chr4_CO$`CO Start`)
WTjaponica_chr4_CO$`CO Start` <- WTjaponica_chr4_CO$`CO Start` - min(WTjaponica_chr4_CO$`CO Start`)

WTjaponica_chr5_CO$`CO End` <- WTjaponica_chr5_CO$`CO End` - min(WTjaponica_chr5_CO$`CO Start`)
WTjaponica_chr5_CO$`CO Start` <- WTjaponica_chr5_CO$`CO Start` - min(WTjaponica_chr5_CO$`CO Start`)

WTjaponica_chr6_CO$`CO End` <- WTjaponica_chr6_CO$`CO End` - min(WTjaponica_chr6_CO$`CO Start`)
WTjaponica_chr6_CO$`CO Start` <- WTjaponica_chr6_CO$`CO Start` - min(WTjaponica_chr6_CO$`CO Start`)

WTjaponica_chr7_CO$`CO End` <- WTjaponica_chr7_CO$`CO End` - min(WTjaponica_chr7_CO$`CO Start`)
WTjaponica_chr7_CO$`CO Start` <- WTjaponica_chr7_CO$`CO Start` - min(WTjaponica_chr7_CO$`CO Start`)

WTjaponica_chr8_CO$`CO End` <- WTjaponica_chr8_CO$`CO End` - min(WTjaponica_chr8_CO$`CO Start`)
WTjaponica_chr8_CO$`CO Start` <- WTjaponica_chr8_CO$`CO Start` - min(WTjaponica_chr8_CO$`CO Start`)

WTjaponica_chr9_CO$`CO End` <- WTjaponica_chr9_CO$`CO End` - min(WTjaponica_chr9_CO$`CO Start`)
WTjaponica_chr9_CO$`CO Start` <- WTjaponica_chr9_CO$`CO Start` - min(WTjaponica_chr9_CO$`CO Start`)

WTjaponica_chr10_CO$`CO End` <- WTjaponica_chr10_CO$`CO End` - min(WTjaponica_chr10_CO$`CO Start`)
WTjaponica_chr10_CO$`CO Start` <- WTjaponica_chr10_CO$`CO Start` - min(WTjaponica_chr10_CO$`CO Start`)

WTjaponica_chr11_CO$`CO End` <- WTjaponica_chr11_CO$`CO End` - min(WTjaponica_chr11_CO$`CO Start`)
WTjaponica_chr11_CO$`CO Start` <- WTjaponica_chr11_CO$`CO Start` - min(WTjaponica_chr11_CO$`CO Start`)

WTjaponica_chr12_CO$`CO End` <- WTjaponica_chr12_CO$`CO End` - min(WTjaponica_chr12_CO$`CO Start`)
WTjaponica_chr12_CO$`CO Start` <- WTjaponica_chr12_CO$`CO Start` - min(WTjaponica_chr12_CO$`CO Start`)

## Function to calculate new rates: loop through fine scale recombination rate data and multiply by avg rate data, then add to old rate
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
recq4_chr1_CO_2<-new_rates(japonica_chr1_CO, WTjaponica_chr1_CO)
recq4_chr2_CO_2<-new_rates(japonica_chr2_CO, WTjaponica_chr2_CO)
recq4_chr3_CO_2<-new_rates(japonica_chr3_CO, WTjaponica_chr3_CO)
recq4_chr4_CO_2<-new_rates(japonica_chr4_CO, WTjaponica_chr4_CO)
recq4_chr5_CO_2<-new_rates(japonica_chr5_CO, WTjaponica_chr5_CO)
recq4_chr6_CO_2<-new_rates(japonica_chr6_CO, WTjaponica_chr6_CO)
recq4_chr7_CO_2<-new_rates(japonica_chr7_CO, WTjaponica_chr7_CO)
recq4_chr8_CO_2<-new_rates(japonica_chr8_CO, WTjaponica_chr8_CO)
recq4_chr9_CO_2<-new_rates(japonica_chr9_CO, WTjaponica_chr9_CO)
recq4_chr10_CO_2<-new_rates(japonica_chr10_CO, WTjaponica_chr10_CO)
recq4_chr11_CO_2<-new_rates(japonica_chr11_CO, WTjaponica_chr11_CO)
recq4_chr12_CO_2<-new_rates(japonica_chr12_CO, WTjaponica_chr12_CO)

##Bin the recombination rate data into 40 equally sized bins on each chromosome
#Function to label the start SNP positions of each bin
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

#Bin size (variable: bins) calculated by dividing the total size of dataframe by 40
#Rollapply function "rolls" by interval of bin size and calculates a moving average recombination rate for each bin
#Use drop_NA function to condense table so it only contains the avg recombination rate

library(zoo)
library(tidyr)
recq4_chr1_CO_3 <- recq4_chr1_CO_2
bins<-as.integer(nrow(recq4_chr1_CO_2)/40)
recq4_chr1_CO_3$rates<- rollapply(recq4_chr1_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
recq4_chr1_CO_3<-fill_start(recq4_chr1_CO_3)
recq4_chr1_CO_3<- recq4_chr1_CO_3 %>% drop_na(rates)

recq4_chr2_CO_3 <- recq4_chr2_CO_2
bins<-as.integer(nrow(recq4_chr2_CO_2)/40)
recq4_chr2_CO_3$rates<- rollapply(recq4_chr2_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
recq4_chr2_CO_3<-fill_start(recq4_chr2_CO_3)
recq4_chr2_CO_3<- recq4_chr2_CO_3 %>% drop_na(rates)

recq4_chr3_CO_3 <- recq4_chr3_CO_2
bins<-as.integer(nrow(recq4_chr3_CO_2)/40)
recq4_chr3_CO_3$rates<- rollapply(recq4_chr3_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
recq4_chr3_CO_3<-fill_start(recq4_chr3_CO_3)
recq4_chr3_CO_3<- recq4_chr3_CO_3 %>% drop_na(rates)

recq4_chr4_CO_3 <- recq4_chr4_CO_2
bins<-as.integer(nrow(recq4_chr4_CO_2)/40)
recq4_chr4_CO_3$rates<- rollapply(recq4_chr4_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
recq4_chr4_CO_3<-fill_start(recq4_chr4_CO_3)
recq4_chr4_CO_3<- recq4_chr4_CO_3 %>% drop_na(rates)

recq4_chr5_CO_3 <- recq4_chr5_CO_2
bins<-as.integer(nrow(recq4_chr5_CO_2)/40)
recq4_chr5_CO_3$rates<- rollapply(recq4_chr5_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
recq4_chr5_CO_3<-fill_start(recq4_chr5_CO_3)
recq4_chr5_CO_3<- recq4_chr5_CO_3 %>% drop_na(rates)

recq4_chr6_CO_3 <- recq4_chr6_CO_2
bins<-as.integer(nrow(recq4_chr6_CO_2)/40)
recq4_chr6_CO_3$rates<- rollapply(recq4_chr6_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
recq4_chr6_CO_3<-fill_start(recq4_chr6_CO_3)
recq4_chr6_CO_3<- recq4_chr6_CO_3 %>% drop_na(rates)

recq4_chr7_CO_3 <- recq4_chr7_CO_2
bins<-as.integer(nrow(recq4_chr7_CO_2)/40)
recq4_chr7_CO_3$rates<- rollapply(recq4_chr7_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
recq4_chr7_CO_3<-fill_start(recq4_chr7_CO_3)
recq4_chr7_CO_3<- recq4_chr7_CO_3 %>% drop_na(rates)

recq4_chr8_CO_3 <- recq4_chr8_CO_2
bins<-as.integer(nrow(recq4_chr8_CO_2)/40)
recq4_chr8_CO_3$rates<- rollapply(recq4_chr8_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
recq4_chr8_CO_3<-fill_start(recq4_chr8_CO_3)
recq4_chr8_CO_3<- recq4_chr8_CO_3 %>% drop_na(rates)

recq4_chr9_CO_3 <- recq4_chr9_CO_2
bins<-as.integer(nrow(recq4_chr9_CO_2)/40)
recq4_chr9_CO_3$rates<- rollapply(recq4_chr9_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
recq4_chr9_CO_3<-fill_start(recq4_chr9_CO_3)
recq4_chr9_CO_3<- recq4_chr9_CO_3 %>% drop_na(rates)

recq4_chr10_CO_3 <- recq4_chr10_CO_2
bins<-as.integer(nrow(recq4_chr10_CO_2)/40)
recq4_chr10_CO_3$rates<- rollapply(recq4_chr10_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
recq4_chr10_CO_3<-fill_start(recq4_chr10_CO_3)
recq4_chr10_CO_3<- recq4_chr10_CO_3 %>% drop_na(rates)

recq4_chr11_CO_3 <- recq4_chr11_CO_2
bins<-as.integer(nrow(recq4_chr11_CO_2)/40)
recq4_chr11_CO_3$rates<- rollapply(recq4_chr11_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
recq4_chr11_CO_3<-fill_start(recq4_chr11_CO_3)
recq4_chr11_CO_3<- recq4_chr11_CO_3 %>% drop_na(rates)

recq4_chr12_CO_3 <- recq4_chr12_CO_2
bins<-as.integer(nrow(recq4_chr12_CO_2)/40)
recq4_chr12_CO_3$rates<- rollapply(recq4_chr12_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
recq4_chr12_CO_3<-fill_start(recq4_chr12_CO_3)
recq4_chr12_CO_3<- recq4_chr12_CO_3 %>% drop_na(rates)


##Assign recombination rate to each SNP
#Function simultaneously loops through recombination rate data and SNPs and assigns recombination rate to SNP based on closest relative location (i.e. if SNP start and end position falls within a recombination rate interval, it is assigned that rate)
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

recq4_chr1_snp2 <- snp_rate(recq4_chr1_CO_3, recq4_chr1_snp)
recq4_chr2_snp2 <- snp_rate(recq4_chr2_CO_3, recq4_chr2_snp)
recq4_chr3_snp2 <- snp_rate(recq4_chr3_CO_3, recq4_chr3_snp)
recq4_chr4_snp2 <- snp_rate(recq4_chr4_CO_3, recq4_chr4_snp)
recq4_chr5_snp2 <- snp_rate(recq4_chr5_CO_3, recq4_chr5_snp)
recq4_chr6_snp2 <- snp_rate(recq4_chr6_CO_3, recq4_chr6_snp)
recq4_chr7_snp2 <- snp_rate(recq4_chr7_CO_3, recq4_chr7_snp)
recq4_chr8_snp2 <- snp_rate(recq4_chr8_CO_3, recq4_chr8_snp)
recq4_chr9_snp2 <- snp_rate(recq4_chr9_CO_3, recq4_chr9_snp)
recq4_chr10_snp2 <- snp_rate(recq4_chr10_CO_3, recq4_chr10_snp)
recq4_chr11_snp2 <- snp_rate(recq4_chr11_CO_3, recq4_chr11_snp)
recq4_chr12_snp2 <- snp_rate(recq4_chr12_CO_3, recq4_chr12_snp)

#Convert SNP start/end positions from bp to Mb
recq4_chr1_snp2$`SNP Start` <- recq4_chr1_snp2$`SNP Start`/1000000
recq4_chr2_snp2$`SNP Start` <- recq4_chr2_snp2$`SNP Start`/1000000
recq4_chr3_snp2$`SNP Start` <- recq4_chr3_snp2$`SNP Start`/1000000
recq4_chr4_snp2$`SNP Start` <- recq4_chr4_snp2$`SNP Start`/1000000
recq4_chr5_snp2$`SNP Start` <- recq4_chr5_snp2$`SNP Start`/1000000
recq4_chr6_snp2$`SNP Start` <- recq4_chr6_snp2$`SNP Start`/1000000
recq4_chr7_snp2$`SNP Start` <- recq4_chr7_snp2$`SNP Start`/1000000
recq4_chr8_snp2$`SNP Start` <- recq4_chr8_snp2$`SNP Start`/1000000
recq4_chr9_snp2$`SNP Start` <- recq4_chr9_snp2$`SNP Start`/1000000
recq4_chr10_snp2$`SNP Start` <- recq4_chr10_snp2$`SNP Start`/1000000
recq4_chr11_snp2$`SNP Start` <- recq4_chr11_snp2$`SNP Start`/1000000
recq4_chr12_snp2$`SNP Start` <- recq4_chr12_snp2$`SNP Start`/1000000

recq4_chr1_snp2$`SNP End` <- recq4_chr1_snp2$`SNP End`/1000000
recq4_chr2_snp2$`SNP End` <- recq4_chr2_snp2$`SNP End`/1000000
recq4_chr3_snp2$`SNP End` <- recq4_chr3_snp2$`SNP End`/1000000
recq4_chr4_snp2$`SNP End` <- recq4_chr4_snp2$`SNP End`/1000000
recq4_chr5_snp2$`SNP End` <- recq4_chr5_snp2$`SNP End`/1000000
recq4_chr6_snp2$`SNP End` <- recq4_chr6_snp2$`SNP End`/1000000
recq4_chr7_snp2$`SNP End` <- recq4_chr7_snp2$`SNP End`/1000000
recq4_chr8_snp2$`SNP End` <- recq4_chr8_snp2$`SNP End`/1000000
recq4_chr9_snp2$`SNP End` <- recq4_chr9_snp2$`SNP End`/1000000
recq4_chr10_snp2$`SNP End` <- recq4_chr10_snp2$`SNP End`/1000000
recq4_chr11_snp2$`SNP End` <- recq4_chr11_snp2$`SNP End`/1000000
recq4_chr12_snp2$`SNP End` <- recq4_chr12_snp2$`SNP End`/1000000

#Function to fill in SNP positions without assigned rates: assign closest rate 
fill_NA<-function(SNP){
  for(i in 1:nrow(SNP)){
    if(is.na(SNP$rate[i])){
      SNP$rate[i]<-SNP$rate[i-1]
    }
  }
  print(SNP)
}
recq4_chr1_snp2<-fill_NA(recq4_chr1_snp2)
recq4_chr2_snp2<-fill_NA(recq4_chr2_snp2)
recq4_chr3_snp2<-fill_NA(recq4_chr3_snp2)
recq4_chr4_snp2<-fill_NA(recq4_chr4_snp2)
recq4_chr5_snp2<-fill_NA(recq4_chr5_snp2)
recq4_chr6_snp2<-fill_NA(recq4_chr6_snp2)
recq4_chr7_snp2<-fill_NA(recq4_chr7_snp2)
recq4_chr8_snp2<-fill_NA(recq4_chr8_snp2)
recq4_chr9_snp2<-fill_NA(recq4_chr9_snp2)
recq4_chr10_snp2<-fill_NA(recq4_chr10_snp2)
recq4_chr11_snp2<-fill_NA(recq4_chr11_snp2)
recq4_chr12_snp2<-fill_NA(recq4_chr12_snp2)


#Function to calculate genetic position: 
#(previous genetic position)+((current physical position - previous physical position)*(recombination rate at current position))
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

#Create genetic maps: use fxn to calculate genetic position, map the genetic position by physical position
recq4_chr1_snp2$pos <- gen_pos(recq4_chr1_snp2)
plot(recq4_chr1_snp2$`SNP Start`, recq4_chr1_snp2$pos)
recq4_chr1_finalpos <- recq4_chr1_snp2[order(recq4_chr1_snp2$pos),]
is.unsorted(recq4_chr1_finalpos$pos)
plot(recq4_chr1_snp2$`SNP Start`, recq4_chr1_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica recq4 Chromosome 1 Genetic Map")
plot(recq4_chr1_finalpos$`SNP Start`, recq4_chr1_finalpos$pos)

recq4_chr2_snp2$pos <- gen_pos(recq4_chr2_snp2)
plot(recq4_chr2_snp2$`SNP Start`, recq4_chr2_snp2$pos)
recq4_chr2_finalpos <- recq4_chr2_snp2[order(recq4_chr2_snp2$pos),]
is.unsorted(recq4_chr2_finalpos$pos)
plot(recq4_chr2_snp2$`SNP Start`, recq4_chr2_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica recq4 Chromosome 2 Genetic Map")

recq4_chr3_snp2$pos <- gen_pos(recq4_chr3_snp2)
plot(recq4_chr3_snp2$`SNP Start`, recq4_chr3_snp2$pos)
recq4_chr3_finalpos <- recq4_chr3_snp2[order(recq4_chr3_snp2$pos),]
is.unsorted(recq4_chr3_finalpos$pos)
plot(recq4_chr3_snp2$`SNP Start`, recq4_chr3_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica recq4 Chromosome 3 Genetic Map")

recq4_chr4_snp2$pos <- gen_pos(recq4_chr4_snp2)
plot(recq4_chr4_snp2$`SNP Start`, recq4_chr4_snp2$pos)
recq4_chr4_finalpos <- recq4_chr4_snp2[order(recq4_chr4_snp2$pos),]
is.unsorted(recq4_chr4_finalpos$pos)
plot(recq4_chr4_snp2$`SNP Start`, recq4_chr4_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica recq4 Chromosome 4 Genetic Map")

recq4_chr5_snp2$pos <- gen_pos(recq4_chr5_snp2)
plot(recq4_chr5_snp2$`SNP Start`, recq4_chr5_snp2$pos)
recq4_chr5_finalpos <- recq4_chr5_snp2[order(recq4_chr5_snp2$pos),]
is.unsorted(recq4_chr5_finalpos$pos)
plot(recq4_chr5_snp2$`SNP Start`, recq4_chr5_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica recq4 Chromosome 5 Genetic Map")

recq4_chr6_snp2$pos <- gen_pos(recq4_chr6_snp2)
plot(recq4_chr6_snp2$`SNP Start`, recq4_chr6_snp2$pos)
recq4_chr6_finalpos <- recq4_chr6_snp2[order(recq4_chr6_snp2$pos),]
is.unsorted(recq4_chr6_finalpos$pos)
plot(recq4_chr6_snp2$`SNP Start`, recq4_chr6_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica recq4 Chromosome 6 Genetic Map")

recq4_chr7_snp2$pos <- gen_pos(recq4_chr7_snp2)
plot(recq4_chr7_snp2$`SNP Start`, recq4_chr7_snp2$pos)
recq4_chr7_finalpos <- recq4_chr7_snp2[order(recq4_chr7_snp2$pos),]
is.unsorted(recq4_chr7_finalpos$pos)
plot(recq4_chr7_snp2$`SNP Start`, recq4_chr7_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica recq4 Chromosome 7 Genetic Map")

recq4_chr8_snp2$pos <- gen_pos(recq4_chr8_snp2)
plot(recq4_chr8_snp2$`SNP Start`, recq4_chr8_snp2$pos)
recq4_chr8_finalpos <- recq4_chr8_snp2[order(recq4_chr8_snp2$pos),]
is.unsorted(recq4_chr8_finalpos$pos)
plot(recq4_chr8_snp2$`SNP Start`, recq4_chr8_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica recq4 Chromosome 8 Genetic Map")

recq4_chr9_snp2$pos <- gen_pos(recq4_chr9_snp2)
plot(recq4_chr9_snp2$`SNP Start`, recq4_chr9_snp2$pos)
recq4_chr9_finalpos <- recq4_chr9_snp2[order(recq4_chr9_snp2$pos),]
is.unsorted(recq4_chr9_finalpos$pos)
plot(recq4_chr9_snp2$`SNP Start`, recq4_chr9_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica recq4 Chromosome 9 Genetic Map")

recq4_chr10_snp2$pos <- gen_pos(recq4_chr10_snp2)
plot(recq4_chr10_snp2$`SNP Start`, recq4_chr10_snp2$pos)
recq4_chr10_finalpos <- recq4_chr10_snp2[order(recq4_chr10_snp2$pos),]
is.unsorted(recq4_chr10_finalpos$pos)
plot(recq4_chr10_snp2$`SNP Start`, recq4_chr10_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica recq4 Chromosome 10 Genetic Map")

recq4_chr11_snp2$pos <- gen_pos(recq4_chr11_snp2)
plot(recq4_chr11_snp2$`SNP Start`, recq4_chr11_snp2$pos)
recq4_chr11_finalpos <- recq4_chr11_snp2[order(recq4_chr11_snp2$pos),]
is.unsorted(recq4_chr11_finalpos$pos)
plot(recq4_chr11_snp2$`SNP Start`, recq4_chr11_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica recq4 Chromosome 11 Genetic Map")

recq4_chr12_snp2$pos <- gen_pos(recq4_chr12_snp2)
plot(recq4_chr12_snp2$`SNP Start`, recq4_chr12_snp2$pos)
recq4_chr12_finalpos <- recq4_chr12_snp2[order(recq4_chr12_snp2$pos),]
is.unsorted(recq4_chr12_finalpos$pos)
plot(recq4_chr12_snp2$`SNP Start`, recq4_chr12_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica recq4 Chromosome 12 Genetic Map")

##Create final recq4 genetic map
chr1 <- recq4_chr1_finalpos$pos/100
chr2 <- recq4_chr2_finalpos$pos/100
chr3 <- recq4_chr3_finalpos$pos/100
chr4 <- recq4_chr4_finalpos$pos/100
chr5 <- recq4_chr5_finalpos$pos/100
chr6 <- recq4_chr6_finalpos$pos/100
chr7 <- recq4_chr7_finalpos$pos/100
chr8 <- recq4_chr8_finalpos$pos/100
chr9 <- recq4_chr9_finalpos$pos/100
chr10<- recq4_chr10_finalpos$pos/100
chr11 <- recq4_chr11_finalpos$pos/100
chr12<- recq4_chr12_finalpos$pos/100

#Add genetic positions to a list (variable:japonica_map) and label each position by  chromosome and location
segSites<-readRDS("japonica_num_SNP.RData")
recq4_map = vector("list",10)
recq4_map[[1]] = chr1
recq4_map[[2]] = chr2
recq4_map[[3]] = chr3
recq4_map[[4]] = chr4
recq4_map[[5]] = chr5
recq4_map[[6]] = chr6
recq4_map[[7]] = chr7
recq4_map[[8]] = chr8
recq4_map[[9]] = chr9
recq4_map[[10]] = chr10
recq4_map[[11]] = chr11
recq4_map[[12]] = chr12
for(i in 1:12){
  names(recq4_map[[i]]) = paste(i, 1:segSites[i], sep="_")
}

saveRDS(recq4_map, file="recq4_final_map.RData")
#Actual physical centromere positions:http://rice.uga.edu/annotation_pseudo_centromeres.shtml
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

#Function to find centromere genetic position: find closest physical position in table to centromere physical position and save the genetic position in that row
find_centromere<-function(centromere,finalpos){
  row.names(finalpos) <- NULL
  index<- finalpos[which.min(abs(centromere-finalpos$`SNP Start`)),]
  row<-as.integer(rownames(index))
  print(finalpos$pos[row])
}
c1 <-find_centromere(16.7,recq4_chr1_finalpos)
c2 <-find_centromere(13.6,recq4_chr2_finalpos)
c3 <-find_centromere(19.4,recq4_chr3_finalpos)
c4 <-find_centromere(9.7,recq4_chr4_finalpos)
c5 <-find_centromere(12.4,recq4_chr5_finalpos)
c6 <-find_centromere(15.3,recq4_chr6_finalpos)
c7 <-find_centromere(12.1,recq4_chr7_finalpos)
c8 <-find_centromere(12.9,recq4_chr8_finalpos)
c9 <-find_centromere(2.8,recq4_chr9_finalpos)
c10 <-find_centromere(8.2,recq4_chr10_finalpos)
c11 <-find_centromere(12,recq4_chr11_finalpos)
c12 <-find_centromere(11.9,recq4_chr12_finalpos)

recq4_centromere <- c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12)
recq4_centromere <- recq4_centromere/100

saveRDS(recq4_centromere, file="recq4_centromeres.RData")
saveRDS(recq4_chr1_finalpos, file="recq4_chr1_finalpos.RData")
saveRDS(recq4_chr2_finalpos, file="recq4_chr2_finalpos.RData")
saveRDS(recq4_chr3_finalpos, file="recq4_chr3_finalpos.RData")
saveRDS(recq4_chr4_finalpos, file="recq4_chr4_finalpos.RData")
saveRDS(recq4_chr5_finalpos, file="recq4_chr5_finalpos.RData")
saveRDS(recq4_chr6_finalpos, file="recq4_chr6_finalpos.RData")
saveRDS(recq4_chr7_finalpos, file="recq4_chr7_finalpos.RData")
saveRDS(recq4_chr8_finalpos, file="recq4_chr8_finalpos.RData")
saveRDS(recq4_chr9_finalpos, file="recq4_chr9_finalpos.RData")
saveRDS(recq4_chr10_finalpos, file="recq4_chr10_finalpos.RData")
saveRDS(recq4_chr11_finalpos, file="recq4_chr11_finalpos.RData")
saveRDS(recq4_chr12_finalpos, file="recq4_chr12_finalpos.RData")

recq4_nSNP<-c(nrow(recq4_chr1_finalpos),nrow(recq4_chr2_finalpos),nrow(recq4_chr3_finalpos),nrow(recq4_chr4_finalpos),nrow(recq4_chr5_finalpos),nrow(recq4_chr6_finalpos),nrow(recq4_chr7_finalpos),nrow(recq4_chr8_finalpos),nrow(recq4_chr9_finalpos),nrow(recq4_chr10_finalpos),nrow(recq4_chr11_finalpos),nrow(recq4_chr12_finalpos))
saveRDS(recq4_nSNP, file="recq4_num_SNP.RData")


