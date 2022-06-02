library(AlphaSimR)
library(Rcpp)
library(ggplot2)
library(dplyr)

#Set seed to ensure sampling is the same each time
set.seed(420)

#Read in the SNP dataset from japonica subspecies genome, randomly sample 2000 SNPs, and order them by SNP start positions
japonica_snps <- read.table("japonica_SNPs.bed", header =FALSE)
colnames(japonica_snps) <- c("Chr#", "SNP Start", "SNP End")
ddm1_snps <- sample_n(japonica_snps, 2000)
ddm1_snps <- ddm1_snps[order(ddm1_snps$`Chr#`,ddm1_snps$`SNP Start`),]

#Create separate dataframes for SNPs on each chromosome, starting the SNP start positions at zero
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

#Average difference (variable: avg_diff) is avg recombination rate difference between recq4 mutant and wildtype rates
#Function to apply avg difference to chromosome arms: divide chromosome into five segments and apply avg_diff to arms (first, second, fourth, fifth segment)avg_diff <-2.346557
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
ddm1_chr1_CO_3 <- ddm1_chr1_CO_2
bins<-as.integer(nrow(ddm1_chr1_CO_2)/40)
ddm1_chr1_CO_3$rates<- rollapply(ddm1_chr1_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ddm1_chr1_CO_3<-fill_start(ddm1_chr1_CO_3)
ddm1_chr1_CO_3<- ddm1_chr1_CO_3 %>% drop_na(rates)

ddm1_chr2_CO_3 <- ddm1_chr2_CO_2
bins<-as.integer(nrow(ddm1_chr2_CO_2)/40)
ddm1_chr2_CO_3$rates<- rollapply(ddm1_chr2_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ddm1_chr2_CO_3<-fill_start(ddm1_chr2_CO_3)
ddm1_chr2_CO_3<- ddm1_chr2_CO_3 %>% drop_na(rates)

ddm1_chr3_CO_3 <- ddm1_chr3_CO_2
bins<-as.integer(nrow(ddm1_chr3_CO_2)/40)
ddm1_chr3_CO_3$rates<- rollapply(ddm1_chr3_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ddm1_chr3_CO_3<-fill_start(ddm1_chr3_CO_3)
ddm1_chr3_CO_3<- ddm1_chr3_CO_3 %>% drop_na(rates)

ddm1_chr4_CO_3 <- ddm1_chr4_CO_2
bins<-as.integer(nrow(ddm1_chr4_CO_2)/40)
ddm1_chr4_CO_3$rates<- rollapply(ddm1_chr4_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ddm1_chr4_CO_3<-fill_start(ddm1_chr4_CO_3)
ddm1_chr4_CO_3<- ddm1_chr4_CO_3 %>% drop_na(rates)

ddm1_chr5_CO_3 <- ddm1_chr5_CO_2
bins<-as.integer(nrow(ddm1_chr5_CO_2)/40)
ddm1_chr5_CO_3$rates<- rollapply(ddm1_chr5_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ddm1_chr5_CO_3<-fill_start(ddm1_chr5_CO_3)
ddm1_chr5_CO_3<- ddm1_chr5_CO_3 %>% drop_na(rates)

ddm1_chr6_CO_3 <- ddm1_chr6_CO_2
bins<-as.integer(nrow(ddm1_chr6_CO_2)/40)
ddm1_chr6_CO_3$rates<- rollapply(ddm1_chr6_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ddm1_chr6_CO_3<-fill_start(ddm1_chr6_CO_3)
ddm1_chr6_CO_3<- ddm1_chr6_CO_3 %>% drop_na(rates)

ddm1_chr7_CO_3 <- ddm1_chr7_CO_2
bins<-as.integer(nrow(ddm1_chr7_CO_2)/40)
ddm1_chr7_CO_3$rates<- rollapply(ddm1_chr7_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ddm1_chr7_CO_3<-fill_start(ddm1_chr7_CO_3)
ddm1_chr7_CO_3<- ddm1_chr7_CO_3 %>% drop_na(rates)

ddm1_chr8_CO_3 <- ddm1_chr8_CO_2
bins<-as.integer(nrow(ddm1_chr8_CO_2)/40)
ddm1_chr8_CO_3$rates<- rollapply(ddm1_chr8_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ddm1_chr8_CO_3<-fill_start(ddm1_chr8_CO_3)
ddm1_chr8_CO_3<- ddm1_chr8_CO_3 %>% drop_na(rates)

ddm1_chr9_CO_3 <- ddm1_chr9_CO_2
bins<-as.integer(nrow(ddm1_chr9_CO_2)/40)
ddm1_chr9_CO_3$rates<- rollapply(ddm1_chr9_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ddm1_chr9_CO_3<-fill_start(ddm1_chr9_CO_3)
ddm1_chr9_CO_3<- ddm1_chr9_CO_3 %>% drop_na(rates)

ddm1_chr10_CO_3 <- ddm1_chr10_CO_2
bins<-as.integer(nrow(ddm1_chr10_CO_2)/40)
ddm1_chr10_CO_3$rates<- rollapply(ddm1_chr10_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ddm1_chr10_CO_3<-fill_start(ddm1_chr10_CO_3)
ddm1_chr10_CO_3<- ddm1_chr10_CO_3 %>% drop_na(rates)

ddm1_chr11_CO_3 <- ddm1_chr11_CO_2
bins<-as.integer(nrow(ddm1_chr11_CO_2)/40)
ddm1_chr11_CO_3$rates<- rollapply(ddm1_chr11_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ddm1_chr11_CO_3<-fill_start(ddm1_chr11_CO_3)
ddm1_chr11_CO_3<- ddm1_chr11_CO_3 %>% drop_na(rates)

ddm1_chr12_CO_3 <- ddm1_chr12_CO_2
bins<-as.integer(nrow(ddm1_chr12_CO_2)/40)
ddm1_chr12_CO_3$rates<- rollapply(ddm1_chr12_CO_2$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
ddm1_chr12_CO_3<-fill_start(ddm1_chr12_CO_3)
ddm1_chr12_CO_3<- ddm1_chr12_CO_3 %>% drop_na(rates)


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

#Convert SNP start/end positions from bp to Mb
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

#Function to fill in SNP positions without assigned rates: assign closest rate 
fill_NA<-function(SNP){
  for(i in 1:nrow(SNP)){
    if(is.na(SNP$rate[i])){
     SNP$rate[i]<-SNP$rate[i-1]
    }
  }
  print(SNP)
}

ddm1_chr1_snp2<-fill_NA(ddm1_chr1_snp2)
ddm1_chr2_snp2<-fill_NA(ddm1_chr2_snp2)
ddm1_chr3_snp2<-fill_NA(ddm1_chr3_snp2)
ddm1_chr4_snp2<-fill_NA(ddm1_chr4_snp2)
ddm1_chr5_snp2<-fill_NA(ddm1_chr5_snp2)
ddm1_chr6_snp2<-fill_NA(ddm1_chr6_snp2)
ddm1_chr7_snp2<-fill_NA(ddm1_chr7_snp2)
ddm1_chr8_snp2<-fill_NA(ddm1_chr8_snp2)
ddm1_chr9_snp2<-fill_NA(ddm1_chr9_snp2)
ddm1_chr10_snp2<-fill_NA(ddm1_chr10_snp2)
ddm1_chr11_snp2<-fill_NA(ddm1_chr11_snp2)
ddm1_chr12_snp2<-fill_NA(ddm1_chr12_snp2)

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
ddm1_chr1_snp2$pos <- gen_pos(ddm1_chr1_snp2)
plot(ddm1_chr1_snp2$`SNP Start`, ddm1_chr1_snp2$pos, type = "l")
ggplot(ddm1_chr1_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
ddm1_chr1_finalpos <- ddm1_chr1_snp2[order(ddm1_chr1_snp2$pos),]
is.unsorted(ddm1_chr1_finalpos$pos)
plot(ddm1_chr1_snp2$`SNP Start`, ddm1_chr1_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 1 Genetic Map")
plot(ddm1_chr1_finalpos$`SNP Start`, ddm1_chr1_finalpos$pos)

ddm1_chr2_snp2$pos <- gen_pos(ddm1_chr2_snp2)
plot(ddm1_chr2_snp2$`SNP Start`, ddm1_chr2_snp2$pos)
ddm1_chr2_finalpos <- ddm1_chr2_snp2[order(ddm1_chr2_snp2$pos),]
is.unsorted(ddm1_chr2_finalpos$pos)
plot(ddm1_chr2_snp2$`SNP Start`, ddm1_chr2_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 2 Genetic Map")

ddm1_chr3_snp2$pos <- gen_pos(ddm1_chr3_snp2)
plot(ddm1_chr3_snp2$`SNP Start`, ddm1_chr3_snp2$pos, type = "l")
ddm1_chr3_finalpos <- ddm1_chr3_snp2[order(ddm1_chr3_snp2$pos),]
is.unsorted(ddm1_chr3_finalpos$pos)
#ddm1_chr3_finalpos$pos <- ddm1_chr3_finalpos$pos + abs(min(ddm1_chr3_finalpos$pos))
plot(ddm1_chr3_snp2$`SNP Start`, ddm1_chr3_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 3 Genetic Map")

ddm1_chr4_snp2$pos <- gen_pos(ddm1_chr4_snp2)
plot(ddm1_chr4_snp2$`SNP Start`, ddm1_chr4_snp2$pos)
ddm1_chr4_finalpos <- ddm1_chr4_snp2[order(ddm1_chr4_snp2$pos),]
is.unsorted(ddm1_chr4_finalpos$pos)
plot(ddm1_chr4_snp2$`SNP Start`, ddm1_chr4_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 4 Genetic Map")

ddm1_chr5_snp2$pos <- gen_pos(ddm1_chr5_snp2)
plot(ddm1_chr5_snp2$`SNP Start`, ddm1_chr5_snp2$pos)
ddm1_chr5_finalpos <- ddm1_chr5_snp2[order(ddm1_chr5_snp2$pos),]
is.unsorted(ddm1_chr5_finalpos$pos)
plot(ddm1_chr5_snp2$`SNP Start`, ddm1_chr5_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 5 Genetic Map")

ddm1_chr6_snp2$pos <- gen_pos(ddm1_chr6_snp2)
plot(ddm1_chr6_snp2$`SNP Start`, ddm1_chr6_snp2$pos)
ddm1_chr6_finalpos <- ddm1_chr6_snp2[order(ddm1_chr6_snp2$pos),]
is.unsorted(ddm1_chr6_finalpos$pos)
plot(ddm1_chr6_snp2$`SNP Start`, ddm1_chr6_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 6 Genetic Map")

ddm1_chr7_snp2$pos <- gen_pos(ddm1_chr7_snp2)
plot(ddm1_chr7_snp2$`SNP Start`, ddm1_chr7_snp2$pos)
ddm1_chr7_finalpos <- ddm1_chr7_snp2[order(ddm1_chr7_snp2$pos),]
is.unsorted(ddm1_chr7_finalpos$pos)
plot(ddm1_chr7_snp2$`SNP Start`, ddm1_chr7_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 7 Genetic Map")

ddm1_chr8_snp2$pos <- gen_pos(ddm1_chr8_snp2)
plot(ddm1_chr8_snp2$`SNP Start`, ddm1_chr8_snp2$pos)
ddm1_chr8_finalpos <- ddm1_chr8_snp2[order(ddm1_chr8_snp2$pos),]
is.unsorted(ddm1_chr8_finalpos$pos)
plot(ddm1_chr8_snp2$`SNP Start`, ddm1_chr8_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 8 Genetic Map")

ddm1_chr9_snp2$pos <- gen_pos(ddm1_chr9_snp2)
plot(ddm1_chr9_snp2$`SNP Start`, ddm1_chr9_snp2$pos)
ddm1_chr9_finalpos <- ddm1_chr9_snp2[order(ddm1_chr9_snp2$pos),]
is.unsorted(ddm1_chr9_finalpos$pos)
plot(ddm1_chr9_snp2$`SNP Start`, ddm1_chr9_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 9 Genetic Map")

ddm1_chr10_snp2$pos <- gen_pos(ddm1_chr10_snp2)
plot(ddm1_chr10_snp2$`SNP Start`, ddm1_chr10_snp2$pos)
ddm1_chr10_finalpos <- ddm1_chr10_snp2[order(ddm1_chr10_snp2$pos),]
is.unsorted(ddm1_chr10_finalpos$pos)
plot(ddm1_chr10_snp2$`SNP Start`, ddm1_chr10_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 10 Genetic Map")

ddm1_chr11_snp2$pos <- gen_pos(ddm1_chr11_snp2)
plot(ddm1_chr11_snp2$`SNP Start`, ddm1_chr11_snp2$pos)
ddm1_chr11_finalpos <- ddm1_chr11_snp2[order(ddm1_chr11_snp2$pos),]
is.unsorted(ddm1_chr11_finalpos$pos)
plot(ddm1_chr11_snp2$`SNP Start`, ddm1_chr11_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 11 Genetic Map")

ddm1_chr12_snp2$pos <- gen_pos(ddm1_chr12_snp2)
plot(ddm1_chr12_snp2$`SNP Start`, ddm1_chr12_snp2$pos)
ddm1_chr12_finalpos <- ddm1_chr12_snp2[order(ddm1_chr12_snp2$pos),]
is.unsorted(ddm1_chr12_finalpos$pos)
plot(ddm1_chr12_snp2$`SNP Start`, ddm1_chr12_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica ddm1 Chromosome 12 Genetic Map")

##Create final recq4 genetic map
chr1 <- ddm1_chr1_finalpos$pos/100
chr2 <- ddm1_chr2_finalpos$pos/100
chr3 <- ddm1_chr3_finalpos$pos/100
chr4 <- ddm1_chr4_finalpos$pos/100
chr5 <- ddm1_chr5_finalpos$pos/100
chr6 <- ddm1_chr6_finalpos$pos/100
chr7 <- ddm1_chr7_finalpos$pos/100
chr8 <- ddm1_chr8_finalpos$pos/100
chr9 <- ddm1_chr9_finalpos$pos/100
chr10<- ddm1_chr10_finalpos$pos/100
chr11 <- ddm1_chr11_finalpos$pos/100
chr12<- ddm1_chr12_finalpos$pos/100

#Add genetic positions to a list (variable:japonica_map) and label each position by chromosome and location
segSites<-readRDS("japonica_num_SNP.RData")
ddm1_map = vector("list",10)
ddm1_map[[1]] = chr1
ddm1_map[[2]] = chr2
ddm1_map[[3]] = chr3
ddm1_map[[4]] = chr4
ddm1_map[[5]] = chr5
ddm1_map[[6]] = chr6
ddm1_map[[7]] = chr7
ddm1_map[[8]] = chr8
ddm1_map[[9]] = chr9
ddm1_map[[10]] = chr10
ddm1_map[[11]] = chr11
ddm1_map[[12]] = chr12
for(i in 1:12){
  names(ddm1_map[[i]]) = paste(i, 1:segSites[i], sep="_")
}

saveRDS(ddm1_map, file="ddm1_final_map.RData")

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

saveRDS(ddm1_centromere, file="ddm1_centromeres.RData")
saveRDS(ddm1_chr1_finalpos, file="ddm1_chr1_finalpos.RData")
saveRDS(ddm1_chr2_finalpos, file="ddm1_chr2_finalpos.RData")
saveRDS(ddm1_chr3_finalpos, file="ddm1_chr3_finalpos.RData")
saveRDS(ddm1_chr4_finalpos, file="ddm1_chr4_finalpos.RData")
saveRDS(ddm1_chr5_finalpos, file="ddm1_chr5_finalpos.RData")
saveRDS(ddm1_chr6_finalpos, file="ddm1_chr6_finalpos.RData")
saveRDS(ddm1_chr7_finalpos, file="ddm1_chr7_finalpos.RData")
saveRDS(ddm1_chr8_finalpos, file="ddm1_chr8_finalpos.RData")
saveRDS(ddm1_chr9_finalpos, file="ddm1_chr9_finalpos.RData")
saveRDS(ddm1_chr10_finalpos, file="ddm1_chr10_finalpos.RData")
saveRDS(ddm1_chr11_finalpos, file="ddm1_chr11_finalpos.RData")
saveRDS(ddm1_chr12_finalpos, file="ddm1_chr12_finalpos.RData")

ddm1_nSNP<-c(nrow(ddm1_chr1_finalpos),nrow(ddm1_chr2_finalpos),nrow(ddm1_chr3_finalpos),nrow(ddm1_chr4_finalpos),nrow(ddm1_chr5_finalpos),nrow(ddm1_chr6_finalpos),nrow(ddm1_chr7_finalpos),nrow(ddm1_chr8_finalpos),nrow(ddm1_chr9_finalpos),nrow(ddm1_chr10_finalpos),nrow(ddm1_chr11_finalpos),nrow(ddm1_chr12_finalpos))
saveRDS(ddm1_nSNP, file="ddm1_num_SNP.RData")


