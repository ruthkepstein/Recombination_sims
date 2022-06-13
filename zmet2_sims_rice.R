library(AlphaSimR)
library(Rcpp)
library(ggplot2)
library(dplyr)

#Set seed to ensure sampling is the same each time
set.seed(420)

#Read in the SNP dataset from japonica subspecies genome, randomly sample 2000 SNPs, and order them by SNP start positions
japonica_snps <- read.table("japonica_SNPs.bed", header =FALSE)
colnames(japonica_snps) <- c("Chr#", "SNP Start", "SNP End")
zmet2_snps <- sample_n(japonica_snps, 2000)
zmet2_snps <- zmet2_snps[order(zmet2_snps$`Chr#`,zmet2_snps$`SNP Start`),]

#Create separate dataframes for SNPs on each chromosome, starting the SNP start positions at zero
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

##Read in wildtype recombination rate data and create separate dataframe for each chromosome
jap_CO <- read.csv("jap_wt_rate.csv", header = TRUE)
colnames(jap_CO) <- c("Chr", "CO Start", "CO End", "WT_rate")
jap_CO <- jap_CO[order(jap_CO$Chr,jap_CO$`CO Start`),]

jap_chr1_CO <- jap_CO[ which(jap_CO$Chr == "1"),]
jap_chr1_CO <- jap_chr1_CO[order(jap_chr1_CO$`CO Start`),]

jap_chr2_CO <- jap_CO[ which(jap_CO$Chr == "2"),]
jap_chr2_CO <- jap_chr2_CO[order(jap_chr2_CO$`CO Start`),]

jap_chr3_CO <- jap_CO[ which(jap_CO$Chr == "3"),]
jap_chr3_CO <- jap_chr3_CO[order(jap_chr3_CO$`CO Start`),]

jap_chr4_CO <- jap_CO[ which(jap_CO$Chr == "4"),]
jap_chr4_CO <- jap_chr4_CO[order(jap_chr4_CO$`CO Start`),]

jap_chr5_CO <- jap_CO[ which(jap_CO$Chr == "5"),]
jap_chr5_CO <- jap_chr5_CO[order(jap_chr5_CO$`CO Start`),]

jap_chr6_CO <- jap_CO[ which(jap_CO$Chr == "6"),]
jap_chr6_CO <- jap_chr6_CO[order(jap_chr6_CO$`CO Start`),]

jap_chr7_CO <- jap_CO[ which(jap_CO$Chr == "7"),]
jap_chr7_CO <- jap_chr7_CO[order(jap_chr7_CO$`CO Start`),]

jap_chr8_CO <- jap_CO[ which(jap_CO$Chr == "8"),]
jap_chr8_CO <- jap_chr8_CO[order(jap_chr8_CO$`CO Start`),]

jap_chr9_CO <- jap_CO[ which(jap_CO$Chr == "9"),]
jap_chr9_CO <- jap_chr9_CO[order(jap_chr9_CO$`CO Start`),]

jap_chr10_CO <- jap_CO[ which(jap_CO$Chr == "10"),]
jap_chr10_CO <- jap_chr10_CO[order(jap_chr10_CO$`CO Start`),]

jap_chr11_CO <- jap_CO[ which(jap_CO$Chr == "11"),]
jap_chr11_CO <- jap_chr11_CO[order(jap_chr11_CO$`CO Start`),]

jap_chr12_CO <- jap_CO[ which(jap_CO$Chr == "12"),]
jap_chr12_CO <- jap_chr12_CO[order(jap_chr12_CO$`CO Start`),]

#Making each recombination rate intervals start at zero
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

#Average difference (variable: avg_diff) is avg recombination rate difference between zmet2 mutant and wildtype rates
#Function to apply avg difference to pericentromeric region: divide chromosome into five segments and apply avg_diff to middle fifth segment
avg_diff <-1
pericentromeric <- function(CO){
  rownames(CO)<-c(1:nrow(CO))
  CO$avg_rate <- 0
  fifth<- max(CO$`CO End`)/5
  start<-fifth
  end<-fifth*3
  for(i in 1:nrow(CO)){
    if(CO$`CO Start`[i]>= start && CO$`CO End`[i]<= end ){
      CO$avg_rate[i] <- avg_diff
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

#Import fine scale recombination rates and create separate dataframe for each chromosome
WTJap_CO <- read.table("japonica_rec_rate.bed", header = FALSE)
colnames(WTJap_CO) <- c("Chr", "CO Start", "CO End", "rate")
WTJap_CO <- WTJap_CO[order(WTJap_CO$Chr,WTJap_CO$`CO Start`),]

WTJap_chr1_CO <- WTJap_CO[ which(WTJap_CO$Chr == "chr01"),]
WTJap_chr1_CO <- WTJap_chr1_CO[order(WTJap_chr1_CO$`CO Start`),]

WTJap_chr2_CO <- WTJap_CO[ which(WTJap_CO$Chr == "chr02"),]
WTJap_chr2_CO <- WTJap_chr2_CO[order(WTJap_chr2_CO$`CO Start`),]

WTJap_chr3_CO <- WTJap_CO[ which(WTJap_CO$Chr == "chr03"),]
WTJap_chr3_CO <- WTJap_chr3_CO[order(WTJap_chr3_CO$`CO Start`),]

WTJap_chr4_CO <- WTJap_CO[ which(WTJap_CO$Chr == "chr04"),]
WTJap_chr4_CO <- WTJap_chr4_CO[order(WTJap_chr4_CO$`CO Start`),]

WTJap_chr5_CO <- WTJap_CO[ which(WTJap_CO$Chr == "chr05"),]
WTJap_chr5_CO <- WTJap_chr5_CO[order(WTJap_chr5_CO$`CO Start`),]

WTJap_chr6_CO <- WTJap_CO[ which(WTJap_CO$Chr == "chr06"),]
WTJap_chr6_CO <- WTJap_chr6_CO[order(WTJap_chr6_CO$`CO Start`),]

WTJap_chr7_CO <- WTJap_CO[ which(WTJap_CO$Chr == "chr07"),]
WTJap_chr7_CO <- WTJap_chr7_CO[order(WTJap_chr7_CO$`CO Start`),]

WTJap_chr8_CO <- WTJap_CO[ which(WTJap_CO$Chr == "chr08"),]
WTJap_chr8_CO <- WTJap_chr8_CO[order(WTJap_chr8_CO$`CO Start`),]

WTJap_chr9_CO <- WTJap_CO[ which(WTJap_CO$Chr == "chr09"),]
WTJap_chr9_CO <- WTJap_chr9_CO[order(WTJap_chr9_CO$`CO Start`),]

WTJap_chr10_CO <- WTJap_CO[ which(WTJap_CO$Chr == "chr10"),]
WTJap_chr10_CO <- WTJap_chr10_CO[order(WTJap_chr10_CO$`CO Start`),]

WTJap_chr11_CO <- WTJap_CO[ which(WTJap_CO$Chr == "chr11"),]
WTJap_chr11_CO <- WTJap_chr11_CO[order(WTJap_chr11_CO$`CO Start`),]

WTJap_chr12_CO <- WTJap_CO[ which(WTJap_CO$Chr == "chr12"),]
WTJap_chr12_CO <- WTJap_chr12_CO[order(WTJap_chr12_CO$`CO Start`),]

#Make intervals start at 0
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
zmet2_chr1_CO_2<-new_rates(jap_chr1_CO, WTJap_chr1_CO)
zmet2_chr2_CO_2<-new_rates(jap_chr2_CO, WTJap_chr2_CO)
zmet2_chr3_CO_2<-new_rates(jap_chr3_CO, WTJap_chr3_CO)
zmet2_chr4_CO_2<-new_rates(jap_chr4_CO, WTJap_chr4_CO)
zmet2_chr5_CO_2<-new_rates(jap_chr5_CO, WTJap_chr5_CO)
zmet2_chr6_CO_2<-new_rates(jap_chr6_CO, WTJap_chr6_CO)
zmet2_chr7_CO_2<-new_rates(jap_chr7_CO, WTJap_chr7_CO)
zmet2_chr8_CO_2<-new_rates(jap_chr8_CO, WTJap_chr8_CO)
zmet2_chr9_CO_2<-new_rates(jap_chr9_CO, WTJap_chr9_CO)
zmet2_chr10_CO_2<-new_rates(jap_chr10_CO, WTJap_chr10_CO)
zmet2_chr11_CO_2<-new_rates(jap_chr11_CO, WTJap_chr11_CO)
zmet2_chr12_CO_2<-new_rates(jap_chr12_CO, WTJap_chr12_CO)

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
zmet2_chr1_CO_3 <- zmet2_chr1_CO_2
bins<-as.integer(nrow(zmet2_chr1_CO_2)/40)
zmet2_chr1_CO_3$rates<- rollapply(zmet2_chr1_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
zmet2_chr1_CO_3<-fill_start(zmet2_chr1_CO_3)
zmet2_chr1_CO_3<- zmet2_chr1_CO_3 %>% drop_na(rates)

zmet2_chr2_CO_3 <- zmet2_chr2_CO_2
bins<-as.integer(nrow(zmet2_chr2_CO_2)/40)
zmet2_chr2_CO_3$rates<- rollapply(zmet2_chr2_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
zmet2_chr2_CO_3<-fill_start(zmet2_chr2_CO_3)
zmet2_chr2_CO_3<- zmet2_chr2_CO_3 %>% drop_na(rates)

zmet2_chr3_CO_3 <- zmet2_chr3_CO_2
bins<-as.integer(nrow(zmet2_chr3_CO_2)/40)
zmet2_chr3_CO_3$rates<- rollapply(zmet2_chr3_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
zmet2_chr3_CO_3<-fill_start(zmet2_chr3_CO_3)
zmet2_chr3_CO_3<- zmet2_chr3_CO_3 %>% drop_na(rates)

zmet2_chr4_CO_3 <- zmet2_chr4_CO_2
bins<-as.integer(nrow(zmet2_chr4_CO_2)/40)
zmet2_chr4_CO_3$rates<- rollapply(zmet2_chr4_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
zmet2_chr4_CO_3<-fill_start(zmet2_chr4_CO_3)
zmet2_chr4_CO_3<- zmet2_chr4_CO_3 %>% drop_na(rates)

zmet2_chr5_CO_3 <- zmet2_chr5_CO_2
bins<-as.integer(nrow(zmet2_chr5_CO_2)/40)
zmet2_chr5_CO_3$rates<- rollapply(zmet2_chr5_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
zmet2_chr5_CO_3<-fill_start(zmet2_chr5_CO_3)
zmet2_chr5_CO_3<- zmet2_chr5_CO_3 %>% drop_na(rates)

zmet2_chr6_CO_3 <- zmet2_chr6_CO_2
bins<-as.integer(nrow(zmet2_chr6_CO_2)/40)
zmet2_chr6_CO_3$rates<- rollapply(zmet2_chr6_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
zmet2_chr6_CO_3<-fill_start(zmet2_chr6_CO_3)
zmet2_chr6_CO_3<- zmet2_chr6_CO_3 %>% drop_na(rates)

zmet2_chr7_CO_3 <- zmet2_chr7_CO_2
bins<-as.integer(nrow(zmet2_chr7_CO_2)/40)
zmet2_chr7_CO_3$rates<- rollapply(zmet2_chr7_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
zmet2_chr7_CO_3<-fill_start(zmet2_chr7_CO_3)
zmet2_chr7_CO_3<- zmet2_chr7_CO_3 %>% drop_na(rates)

zmet2_chr8_CO_3 <- zmet2_chr8_CO_2
bins<-as.integer(nrow(zmet2_chr8_CO_2)/40)
zmet2_chr8_CO_3$rates<- rollapply(zmet2_chr8_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
zmet2_chr8_CO_3<-fill_start(zmet2_chr8_CO_3)
zmet2_chr8_CO_3<- zmet2_chr8_CO_3 %>% drop_na(rates)

zmet2_chr9_CO_3 <- zmet2_chr9_CO_2
bins<-as.integer(nrow(zmet2_chr9_CO_2)/40)
zmet2_chr9_CO_3$rates<- rollapply(zmet2_chr9_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
zmet2_chr9_CO_3<-fill_start(zmet2_chr9_CO_3)
zmet2_chr9_CO_3<- zmet2_chr9_CO_3 %>% drop_na(rates)

zmet2_chr10_CO_3 <- zmet2_chr10_CO_2
bins<-as.integer(nrow(zmet2_chr10_CO_2)/40)
zmet2_chr10_CO_3$rates<- rollapply(zmet2_chr10_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
zmet2_chr10_CO_3<-fill_start(zmet2_chr10_CO_3)
zmet2_chr10_CO_3<- zmet2_chr10_CO_3 %>% drop_na(rates)

zmet2_chr11_CO_3 <- zmet2_chr11_CO_2
bins<-as.integer(nrow(zmet2_chr11_CO_2)/40)
zmet2_chr11_CO_3$rates<- rollapply(zmet2_chr11_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
zmet2_chr11_CO_3<-fill_start(zmet2_chr11_CO_3)
zmet2_chr11_CO_3<- zmet2_chr11_CO_3 %>% drop_na(rates)

zmet2_chr12_CO_3 <- zmet2_chr12_CO_2
bins<-as.integer(nrow(zmet2_chr12_CO_2)/40)
zmet2_chr12_CO_3$rates<- rollapply(zmet2_chr12_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
zmet2_chr12_CO_3<-fill_start(zmet2_chr12_CO_3)
zmet2_chr12_CO_3<- zmet2_chr12_CO_3 %>% drop_na(rates)

##Assign recombination rate to each SNP
#Function simulataneously loops through recombination rate data and SNPs and assigns recombination rate to SNP based on closest relative location (i.e. if SNP start and end position falls within a recombination rate interval, it is assigned that rate)
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

#Convert SNP start/end positions from bp to Mb
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

#Function to fill in SNP positions without assigned rates: assign closest rate 
fill_NA<-function(SNP){
  for(i in 1:nrow(SNP)){
    if(is.na(SNP$rate[i])){
      SNP$rate[i]<-SNP$rate[i-1]
    }
  }
  print(SNP)
}
zmet2_chr1_snp2<-fill_NA(zmet2_chr1_snp2)
zmet2_chr2_snp2<-fill_NA(zmet2_chr2_snp2)
zmet2_chr3_snp2<-fill_NA(zmet2_chr3_snp2)
zmet2_chr4_snp2<-fill_NA(zmet2_chr4_snp2)
zmet2_chr5_snp2<-fill_NA(zmet2_chr5_snp2)
zmet2_chr6_snp2<-fill_NA(zmet2_chr6_snp2)
zmet2_chr7_snp2<-fill_NA(zmet2_chr7_snp2)
zmet2_chr8_snp2<-fill_NA(zmet2_chr8_snp2)
zmet2_chr9_snp2<-fill_NA(zmet2_chr9_snp2)
zmet2_chr10_snp2<-fill_NA(zmet2_chr10_snp2)
zmet2_chr11_snp2<-fill_NA(zmet2_chr11_snp2)
zmet2_chr12_snp2<-fill_NA(zmet2_chr12_snp2)

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
zmet2_chr1_snp2$pos <-gen_pos(zmet2_chr1_snp2)
plot(zmet2_chr1_snp2$`SNP Start`, zmet2_chr1_snp2$pos)
ggplot(zmet2_chr1_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
zmet2_chr1_finalpos <- zmet2_chr1_snp2[order(zmet2_chr1_snp2$pos),]
is.unsorted(zmet2_chr1_finalpos$pos)
plot(zmet2_chr1_snp2$`SNP Start`, zmet2_chr1_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica zmet2 Chromosome 1 Genetic Map")
plot(zmet2_chr1_finalpos$`SNP Start`, zmet2_chr1_finalpos$pos)

zmet2_chr2_snp2$pos <-gen_pos(zmet2_chr2_snp2)
plot(zmet2_chr2_snp2$`SNP Start`, zmet2_chr2_snp2$pos)
zmet2_chr2_finalpos <- zmet2_chr2_snp2[order(zmet2_chr2_snp2$pos),]
is.unsorted(zmet2_chr2_finalpos$pos)
plot(zmet2_chr2_snp2$`SNP Start`, zmet2_chr2_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica zmet2 Chromosome 2 Genetic Map")

zmet2_chr3_snp2$pos <-gen_pos(zmet2_chr3_snp2)
plot(zmet2_chr3_snp2$`SNP Start`, zmet2_chr3_snp2$pos)
zmet2_chr3_finalpos <- zmet2_chr3_snp2[order(zmet2_chr3_snp2$pos),]
is.unsorted(zmet2_chr3_finalpos$pos)
plot(zmet2_chr3_snp2$`SNP Start`, zmet2_chr3_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica zmet2 Chromosome 3 Genetic Map")

zmet2_chr4_snp2$pos <-gen_pos(zmet2_chr4_snp2)
plot(zmet2_chr4_snp2$`SNP Start`, zmet2_chr4_snp2$pos)
zmet2_chr4_finalpos <- zmet2_chr4_snp2[order(zmet2_chr4_snp2$pos),]
is.unsorted(zmet2_chr4_finalpos$pos)
plot(zmet2_chr4_snp2$`SNP Start`, zmet2_chr4_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica zmet2 Chromosome 4 Genetic Map")

zmet2_chr5_snp2$pos <-gen_pos(zmet2_chr5_snp2)
plot(zmet2_chr5_snp2$`SNP Start`, zmet2_chr5_snp2$pos)
zmet2_chr5_finalpos <- zmet2_chr5_snp2[order(zmet2_chr5_snp2$pos),]
is.unsorted(zmet2_chr5_finalpos$pos)
plot(zmet2_chr5_snp2$`SNP Start`, zmet2_chr5_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica zmet2 Chromosome 5 Genetic Map")

zmet2_chr6_snp2$pos <-gen_pos(zmet2_chr6_snp2)
plot(zmet2_chr6_snp2$`SNP Start`, zmet2_chr6_snp2$pos)
zmet2_chr6_finalpos <- zmet2_chr6_snp2[order(zmet2_chr6_snp2$pos),]
is.unsorted(zmet2_chr6_finalpos$pos)
plot(zmet2_chr6_snp2$`SNP Start`, zmet2_chr6_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica zmet2 Chromosome 6 Genetic Map")

zmet2_chr7_snp2$pos <-gen_pos(zmet2_chr7_snp2)
plot(zmet2_chr7_snp2$`SNP Start`, zmet2_chr7_snp2$pos)
zmet2_chr7_finalpos <- zmet2_chr7_snp2[order(zmet2_chr7_snp2$pos),]
is.unsorted(zmet2_chr7_finalpos$pos)
plot(zmet2_chr7_snp2$`SNP Start`, zmet2_chr7_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica zmet2 Chromosome 7 Genetic Map")

zmet2_chr8_snp2$pos <-gen_pos(zmet2_chr8_snp2)
plot(zmet2_chr8_snp2$`SNP Start`, zmet2_chr8_snp2$pos)
zmet2_chr8_finalpos <- zmet2_chr8_snp2[order(zmet2_chr8_snp2$pos),]
is.unsorted(zmet2_chr8_finalpos$pos)
plot(zmet2_chr8_snp2$`SNP Start`, zmet2_chr8_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica zmet2 Chromosome 8 Genetic Map")

zmet2_chr9_snp2$pos <-gen_pos(zmet2_chr9_snp2)
plot(zmet2_chr9_snp2$`SNP Start`, zmet2_chr9_snp2$pos)
zmet2_chr9_finalpos <- zmet2_chr9_snp2[order(zmet2_chr9_snp2$pos),]
is.unsorted(zmet2_chr9_finalpos$pos)
plot(zmet2_chr9_snp2$`SNP Start`, zmet2_chr9_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica zmet2 Chromosome 9 Genetic Map")

zmet2_chr10_snp2$pos <-gen_pos(zmet2_chr10_snp2)
plot(zmet2_chr10_snp2$`SNP Start`, zmet2_chr10_snp2$pos)
zmet2_chr10_finalpos <- zmet2_chr10_snp2[order(zmet2_chr10_snp2$pos),]
is.unsorted(zmet2_chr10_finalpos$pos)
plot(zmet2_chr10_snp2$`SNP Start`, zmet2_chr10_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica zmet2 Chromosome 10 Genetic Map")

zmet2_chr11_snp2$pos <-gen_pos(zmet2_chr11_snp2)
plot(zmet2_chr11_snp2$`SNP Start`, zmet2_chr11_snp2$pos)
zmet2_chr11_finalpos <- zmet2_chr11_snp2[order(zmet2_chr11_snp2$pos),]
is.unsorted(zmet2_chr11_finalpos$pos)
plot(zmet2_chr11_snp2$`SNP Start`, zmet2_chr11_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica zmet2 Chromosome 11 Genetic Map")

zmet2_chr12_snp2$pos <-gen_pos(zmet2_chr12_snp2)
plot(zmet2_chr12_snp2$`SNP Start`, zmet2_chr12_snp2$pos)
zmet2_chr12_finalpos <- zmet2_chr12_snp2[order(zmet2_chr12_snp2$pos),]
is.unsorted(zmet2_chr12_finalpos$pos)
plot(zmet2_chr12_snp2$`SNP Start`, zmet2_chr12_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica zmet2 Chromosome 12 Genetic Map")

##Create final zmet2 genetic map
chr1 <- zmet2_chr1_finalpos$pos/100
chr2 <- zmet2_chr2_finalpos$pos/100
chr3 <- zmet2_chr3_finalpos$pos/100
chr4 <- zmet2_chr4_finalpos$pos/100
chr5 <- zmet2_chr5_finalpos$pos/100
chr6 <- zmet2_chr6_finalpos$pos/100
chr7 <- zmet2_chr7_finalpos$pos/100
chr8 <- zmet2_chr8_finalpos$pos/100
chr9 <- zmet2_chr9_finalpos$pos/100
chr10<- zmet2_chr10_finalpos$pos/100
chr11 <- zmet2_chr11_finalpos$pos/100
chr12<- zmet2_chr12_finalpos$pos/100

#Add genetic positions to a list (variable:japonica_map) and label each position by  chromosome and location
segSites<-readRDS("japonica_num_SNP.RData")
zmet2_map = vector("list",10)
zmet2_map[[1]] = chr1
zmet2_map[[2]] = chr2
zmet2_map[[3]] = chr3
zmet2_map[[4]] = chr4
zmet2_map[[5]] = chr5
zmet2_map[[6]] = chr6
zmet2_map[[7]] = chr7
zmet2_map[[8]] = chr8
zmet2_map[[9]] = chr9
zmet2_map[[10]] = chr10
zmet2_map[[11]] = chr11
zmet2_map[[12]] = chr12
for(i in 1:12){
  names(zmet2_map[[i]]) = paste(i, 1:segSites[i], sep="_")
}

saveRDS(zmet2_map, file="zmet2_final_map.RData")

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

saveRDS(zmet2_centromere, file="zmet2_centromeres.RData")

saveRDS(zmet2_chr1_finalpos, file="zmet2_chr1_finalpos.RData")
saveRDS(zmet2_chr2_finalpos, file="zmet2_chr2_finalpos.RData")
saveRDS(zmet2_chr3_finalpos, file="zmet2_chr3_finalpos.RData")
saveRDS(zmet2_chr4_finalpos, file="zmet2_chr4_finalpos.RData")
saveRDS(zmet2_chr5_finalpos, file="zmet2_chr5_finalpos.RData")
saveRDS(zmet2_chr6_finalpos, file="zmet2_chr6_finalpos.RData")
saveRDS(zmet2_chr7_finalpos, file="zmet2_chr7_finalpos.RData")
saveRDS(zmet2_chr8_finalpos, file="zmet2_chr8_finalpos.RData")
saveRDS(zmet2_chr9_finalpos, file="zmet2_chr9_finalpos.RData")
saveRDS(zmet2_chr10_finalpos, file="zmet2_chr10_finalpos.RData")
saveRDS(zmet2_chr11_finalpos, file="zmet2_chr11_finalpos.RData")
saveRDS(zmet2_chr12_finalpos, file="zmet2_chr12_finalpos.RData")

zmet2_nSNP<-c(nrow(zmet2_chr1_finalpos),nrow(zmet2_chr2_finalpos),nrow(zmet2_chr3_finalpos),nrow(zmet2_chr4_finalpos),nrow(zmet2_chr5_finalpos),nrow(zmet2_chr6_finalpos),nrow(zmet2_chr7_finalpos),nrow(zmet2_chr8_finalpos),nrow(zmet2_chr9_finalpos),nrow(zmet2_chr10_finalpos),nrow(zmet2_chr11_finalpos),nrow(zmet2_chr12_finalpos))
saveRDS(zmet2_nSNP, file="zmet2_num_SNP.RData")


