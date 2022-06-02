library(AlphaSimR)
library(Rcpp)
library(ggplot2)
library(dplyr)

#Set seed to ensure sampling is the same each time
set.seed(420)

#Read in the SNP dataset from japonica subspecies genome, randomly sample 2000 SNPs, and order them by SNP start positions
japonica_snps <- read.table("japonica_SNPs.bed", header =FALSE)
colnames(japonica_snps) <- c("Chr#", "SNP Start", "SNP End")
japonica_snps <- sample_n(japonica_snps, 2000)
japonica_snps <- japonica_snps[order(japonica_snps$`Chr#`,japonica_snps$`SNP Start`),]

#Create separate dataframes for SNPs on each chromosome, starting the SNP start positions at zero
japonica_chr1_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr1"),]
japonica_chr1_snp$rate <- NA
japonica_chr1_snp$`SNP End` <- japonica_chr1_snp$`SNP End` - min(japonica_chr1_snp$`SNP Start`)
japonica_chr1_snp$`SNP Start` <- japonica_chr1_snp$`SNP Start`- min(japonica_chr1_snp$`SNP Start`)

japonica_chr2_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr2"),]
japonica_chr2_snp$rate <- NA
japonica_chr2_snp$`SNP End` <- japonica_chr2_snp$`SNP End` - min(japonica_chr2_snp$`SNP Start`)
japonica_chr2_snp$`SNP Start` <- japonica_chr2_snp$`SNP Start`- min(japonica_chr2_snp$`SNP Start`)

japonica_chr3_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr3"),]
japonica_chr3_snp$rate <- NA
japonica_chr3_snp$`SNP End` <- japonica_chr3_snp$`SNP End` - min(japonica_chr3_snp$`SNP Start`)
japonica_chr3_snp$`SNP Start` <- japonica_chr3_snp$`SNP Start`- min(japonica_chr3_snp$`SNP Start`)

japonica_chr4_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr4"),]
japonica_chr4_snp$rate <- NA
japonica_chr4_snp$`SNP End` <- japonica_chr4_snp$`SNP End` - min(japonica_chr4_snp$`SNP Start`)
japonica_chr4_snp$`SNP Start` <- japonica_chr4_snp$`SNP Start`- min(japonica_chr4_snp$`SNP Start`)

japonica_chr5_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr5"),]
japonica_chr5_snp$rate <- NA
japonica_chr5_snp$`SNP End` <- japonica_chr5_snp$`SNP End` - min(japonica_chr5_snp$`SNP Start`)
japonica_chr5_snp$`SNP Start` <- japonica_chr5_snp$`SNP Start`- min(japonica_chr5_snp$`SNP Start`)

japonica_chr6_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr6"),]
japonica_chr6_snp$rate <- NA
japonica_chr6_snp$`SNP End` <- japonica_chr6_snp$`SNP End` - min(japonica_chr6_snp$`SNP Start`)
japonica_chr6_snp$`SNP Start` <- japonica_chr6_snp$`SNP Start`- min(japonica_chr6_snp$`SNP Start`)

japonica_chr7_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr7"),]
japonica_chr7_snp$rate <- NA
japonica_chr7_snp$`SNP End` <- japonica_chr7_snp$`SNP End` - min(japonica_chr7_snp$`SNP Start`)
japonica_chr7_snp$`SNP Start` <- japonica_chr7_snp$`SNP Start`- min(japonica_chr7_snp$`SNP Start`)

japonica_chr8_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr8"),]
japonica_chr8_snp$rate <- NA
japonica_chr8_snp$`SNP End` <- japonica_chr8_snp$`SNP End` - min(japonica_chr8_snp$`SNP Start`)
japonica_chr8_snp$`SNP Start` <- japonica_chr8_snp$`SNP Start`- min(japonica_chr8_snp$`SNP Start`)

japonica_chr9_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr9"),]
japonica_chr9_snp$rate <- NA
japonica_chr9_snp$`SNP End` <- japonica_chr9_snp$`SNP End` - min(japonica_chr9_snp$`SNP Start`)
japonica_chr9_snp$`SNP Start` <- japonica_chr9_snp$`SNP Start`- min(japonica_chr9_snp$`SNP Start`)

japonica_chr10_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr10"),]
japonica_chr10_snp$rate <- NA
japonica_chr10_snp$`SNP End` <- japonica_chr10_snp$`SNP End` - min(japonica_chr10_snp$`SNP Start`)
japonica_chr10_snp$`SNP Start` <- japonica_chr10_snp$`SNP Start`- min(japonica_chr10_snp$`SNP Start`)

japonica_chr11_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr11"),]
japonica_chr11_snp$rate <- NA
japonica_chr11_snp$`SNP End` <- japonica_chr11_snp$`SNP End` - min(japonica_chr11_snp$`SNP Start`)
japonica_chr11_snp$`SNP Start` <- japonica_chr11_snp$`SNP Start`- min(japonica_chr11_snp$`SNP Start`)

japonica_chr12_snp <- japonica_snps[ which(japonica_snps$`Chr#` == "Chr12"),]
japonica_chr12_snp$rate <- NA
japonica_chr12_snp$`SNP End` <- japonica_chr12_snp$`SNP End` - min(japonica_chr12_snp$`SNP Start`)
japonica_chr12_snp$`SNP Start` <- japonica_chr12_snp$`SNP Start`- min(japonica_chr12_snp$`SNP Start`)

#Read in genome-wide recombination rate data and creating separate dataframes for each chromosomes
japonica_CO <- read.table("japonica_rec_rate.bed", header = FALSE)
colnames(japonica_CO) <- c("Chr", "CO Start", "CO End", "rate")
japonica_CO <- japonica_CO[order(japonica_CO$Chr,japonica_CO$`CO Start`),]

japonica_chr1_CO <- japonica_CO[ which(japonica_CO$Chr == "chr01"),]
japonica_chr1_CO <- japonica_chr1_CO[order(japonica_chr1_CO$`CO Start`),]

japonica_chr2_CO <- japonica_CO[ which(japonica_CO$Chr == "chr02"),]
japonica_chr2_CO <- japonica_chr2_CO[order(japonica_chr2_CO$`CO Start`),]

japonica_chr3_CO <- japonica_CO[ which(japonica_CO$Chr == "chr03"),]
japonica_chr3_CO <- japonica_chr3_CO[order(japonica_chr3_CO$`CO Start`),]

japonica_chr4_CO <- japonica_CO[ which(japonica_CO$Chr == "chr04"),]
japonica_chr4_CO <- japonica_chr4_CO[order(japonica_chr4_CO$`CO Start`),]

japonica_chr5_CO <- japonica_CO[ which(japonica_CO$Chr == "chr05"),]
japonica_chr5_CO <- japonica_chr5_CO[order(japonica_chr5_CO$`CO Start`),]

japonica_chr6_CO <- japonica_CO[ which(japonica_CO$Chr == "chr06"),]
japonica_chr6_CO <- japonica_chr6_CO[order(japonica_chr6_CO$`CO Start`),]

japonica_chr7_CO <- japonica_CO[ which(japonica_CO$Chr == "chr07"),]
japonica_chr7_CO <- japonica_chr7_CO[order(japonica_chr7_CO$`CO Start`),]

japonica_chr8_CO <- japonica_CO[ which(japonica_CO$Chr == "chr08"),]
japonica_chr8_CO <- japonica_chr8_CO[order(japonica_chr8_CO$`CO Start`),]

japonica_chr9_CO <- japonica_CO[ which(japonica_CO$Chr == "chr09"),]
japonica_chr9_CO <- japonica_chr9_CO[order(japonica_chr9_CO$`CO Start`),]

japonica_chr10_CO <- japonica_CO[ which(japonica_CO$Chr == "chr10"),]
japonica_chr10_CO <- japonica_chr10_CO[order(japonica_chr10_CO$`CO Start`),]

japonica_chr11_CO <- japonica_CO[ which(japonica_CO$Chr == "chr11"),]
japonica_chr11_CO <- japonica_chr11_CO[order(japonica_chr11_CO$`CO Start`),]

japonica_chr12_CO <- japonica_CO[ which(japonica_CO$Chr == "chr12"),]
japonica_chr12_CO <- japonica_chr12_CO[order(japonica_chr12_CO$`CO Start`),]

library(dlookr)
library(tidyverse)
library(OneR)

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

#Verifying that recombination rate intervals start at zero
isTRUE(japonica_chr1_CO[1,2] == 0)
isTRUE(japonica_chr2_CO[1,2] == 0)
isTRUE(japonica_chr3_CO[1,2] == 0)
isTRUE(japonica_chr4_CO[1,2] == 0)
isTRUE(japonica_chr5_CO[1,2] == 0)
isTRUE(japonica_chr6_CO[1,2] == 0)
isTRUE(japonica_chr7_CO[1,2] == 0)
isTRUE(japonica_chr8_CO[1,2] == 0)
isTRUE(japonica_chr9_CO[1,2] == 0)
isTRUE(japonica_chr10_CO[1,2] == 0)
isTRUE(japonica_chr11_CO[1,2] == 0)
isTRUE(japonica_chr12_CO[1,2] == 0)

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
japonica_chr1_CO_2 <- japonica_chr1_CO
bins<-as.integer(nrow(japonica_chr1_CO)/40)
japonica_chr1_CO_2$rates<- rollapply(japonica_chr1_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr1_CO_2<-fill_start(japonica_chr1_CO_2)
japonica_chr1_CO_2<- japonica_chr1_CO_2 %>% drop_na(rates)

japonica_chr2_CO_2 <- japonica_chr2_CO
bins<-as.integer(nrow(japonica_chr2_CO)/40)
japonica_chr2_CO_2$rates<- rollapply(japonica_chr2_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr2_CO_2<-fill_start(japonica_chr2_CO_2)
japonica_chr2_CO_2<- japonica_chr2_CO_2 %>% drop_na(rates)

japonica_chr3_CO_2 <- japonica_chr3_CO
bins<-as.integer(nrow(japonica_chr3_CO)/40)
japonica_chr3_CO_2$rates<- rollapply(japonica_chr3_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr3_CO_2<-fill_start(japonica_chr3_CO_2)
japonica_chr3_CO_2<- japonica_chr3_CO_2 %>% drop_na(rates)

japonica_chr4_CO_2 <- japonica_chr4_CO
bins<-as.integer(nrow(japonica_chr4_CO)/40)
japonica_chr4_CO_2$rates<- rollapply(japonica_chr4_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr4_CO_2<-fill_start(japonica_chr4_CO_2)
japonica_chr4_CO_2<- japonica_chr4_CO_2 %>% drop_na(rates)

japonica_chr5_CO_2 <- japonica_chr5_CO
bins<-as.integer(nrow(japonica_chr5_CO)/40)
japonica_chr5_CO_2$rates<- rollapply(japonica_chr5_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr5_CO_2<-fill_start(japonica_chr5_CO_2)
japonica_chr5_CO_2<- japonica_chr5_CO_2 %>% drop_na(rates)

japonica_chr6_CO_2 <- japonica_chr6_CO
bins<-as.integer(nrow(japonica_chr6_CO)/40)
japonica_chr6_CO_2$rates<- rollapply(japonica_chr6_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr6_CO_2<-fill_start(japonica_chr6_CO_2)
japonica_chr6_CO_2<- japonica_chr6_CO_2 %>% drop_na(rates)

japonica_chr7_CO_2 <- japonica_chr7_CO
bins<-as.integer(nrow(japonica_chr7_CO)/40)
japonica_chr7_CO_2$rates<- rollapply(japonica_chr7_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr7_CO_2<-fill_start(japonica_chr7_CO_2)
japonica_chr7_CO_2<- japonica_chr7_CO_2 %>% drop_na(rates)

japonica_chr8_CO_2 <- japonica_chr8_CO
bins<-as.integer(nrow(japonica_chr8_CO)/40)
japonica_chr8_CO_2$rates<- rollapply(japonica_chr8_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr8_CO_2<-fill_start(japonica_chr8_CO_2)
japonica_chr8_CO_2<- japonica_chr8_CO_2 %>% drop_na(rates)

japonica_chr9_CO_2 <- japonica_chr9_CO
bins<-as.integer(nrow(japonica_chr9_CO)/40)
japonica_chr9_CO_2$rates<- rollapply(japonica_chr9_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr9_CO_2<-fill_start(japonica_chr9_CO_2)
japonica_chr9_CO_2<- japonica_chr9_CO_2 %>% drop_na(rates)

japonica_chr10_CO_2 <- japonica_chr10_CO
bins<-as.integer(nrow(japonica_chr10_CO)/40)
japonica_chr10_CO_2$rates<- rollapply(japonica_chr10_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr10_CO_2<-fill_start(japonica_chr10_CO_2)
japonica_chr10_CO_2<- japonica_chr10_CO_2 %>% drop_na(rates)

japonica_chr11_CO_2 <- japonica_chr11_CO
bins<-as.integer(nrow(japonica_chr11_CO)/40)
japonica_chr11_CO_2$rates<- rollapply(japonica_chr11_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr11_CO_2<-fill_start(japonica_chr11_CO_2)
japonica_chr11_CO_2<- japonica_chr11_CO_2 %>% drop_na(rates)

japonica_chr12_CO_2 <- japonica_chr12_CO
bins<-as.integer(nrow(japonica_chr12_CO)/40)
japonica_chr12_CO_2$rates<- rollapply(japonica_chr12_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr12_CO_2<-fill_start(japonica_chr12_CO_2)
japonica_chr12_CO_2<- japonica_chr12_CO_2 %>% drop_na(rates)

##Assign recombination rate to each SNP
#Function simulataneously loops through recombination rate data and SNPs and assigns recombination rate to SNP based on closest relative location (i.e. if SNP start and end position falls within a recombination rate interval, it is assigned that rate)
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

japonica_chr1_snp2 <- snp_rate(japonica_chr1_CO_2, japonica_chr1_snp)
japonica_chr2_snp2 <- snp_rate(japonica_chr2_CO_2, japonica_chr2_snp)
japonica_chr3_snp2 <- snp_rate(japonica_chr3_CO_2, japonica_chr3_snp)
japonica_chr4_snp2 <- snp_rate(japonica_chr4_CO_2, japonica_chr4_snp)
japonica_chr5_snp2 <- snp_rate(japonica_chr5_CO_2, japonica_chr5_snp)
japonica_chr6_snp2 <- snp_rate(japonica_chr6_CO_2, japonica_chr6_snp)
japonica_chr7_snp2 <- snp_rate(japonica_chr7_CO_2, japonica_chr7_snp)
japonica_chr8_snp2 <- snp_rate(japonica_chr8_CO_2, japonica_chr8_snp)
japonica_chr9_snp2 <- snp_rate(japonica_chr9_CO_2, japonica_chr9_snp)
japonica_chr10_snp2 <- snp_rate(japonica_chr10_CO_2, japonica_chr10_snp)
japonica_chr11_snp2 <- snp_rate(japonica_chr11_CO_2, japonica_chr11_snp)
japonica_chr12_snp2 <- snp_rate(japonica_chr12_CO_2, japonica_chr12_snp)

#Make mutable copies of data
japonica_chr1_snp3 <- japonica_chr1_snp2
japonica_chr2_snp3 <- japonica_chr2_snp2
japonica_chr3_snp3 <- japonica_chr3_snp2
japonica_chr4_snp3 <- japonica_chr4_snp2
japonica_chr5_snp3 <- japonica_chr5_snp2
japonica_chr6_snp3 <- japonica_chr6_snp2
japonica_chr7_snp3 <- japonica_chr7_snp2
japonica_chr8_snp3 <- japonica_chr8_snp2
japonica_chr9_snp3 <- japonica_chr9_snp2
japonica_chr10_snp3 <- japonica_chr10_snp2
japonica_chr11_snp3 <- japonica_chr11_snp2
japonica_chr12_snp3 <- japonica_chr12_snp2

#Convert start position from bp to Mb
japonica_chr1_snp3$`SNP Start`<- japonica_chr1_snp3$`SNP Start`/1000000
japonica_chr2_snp3$`SNP Start` <- japonica_chr2_snp3$`SNP Start`/1000000
japonica_chr3_snp3$`SNP Start` <- japonica_chr3_snp3$`SNP Start`/1000000
japonica_chr4_snp3$`SNP Start` <- japonica_chr4_snp3$`SNP Start`/1000000
japonica_chr5_snp3$`SNP Start` <- japonica_chr5_snp3$`SNP Start`/1000000
japonica_chr6_snp3$`SNP Start` <- japonica_chr6_snp3$`SNP Start`/1000000
japonica_chr7_snp3$`SNP Start` <- japonica_chr7_snp3$`SNP Start`/1000000
japonica_chr8_snp3$`SNP Start` <- japonica_chr8_snp3$`SNP Start`/1000000
japonica_chr9_snp3$`SNP Start` <- japonica_chr9_snp3$`SNP Start`/1000000
japonica_chr10_snp3$`SNP Start` <- japonica_chr10_snp3$`SNP Start`/1000000
japonica_chr11_snp3$`SNP Start` <- japonica_chr11_snp3$`SNP Start`/1000000
japonica_chr12_snp3$`SNP Start` <- japonica_chr12_snp3$`SNP Start`/1000000

#Function to fill in SNP positions without assigned rates: assign closest rate 
fill_NA<-function(SNP){
  for(i in 1:nrow(SNP)){
    if(is.na(SNP$rate[i])){
      SNP$rate[i]<-SNP$rate[i-1]
    }
  }
  print(SNP)
}
japonica_chr1_snp3<-fill_NA(japonica_chr1_snp3)
japonica_chr2_snp3<-fill_NA(japonica_chr2_snp3)
japonica_chr3_snp3<-fill_NA(japonica_chr3_snp3)
japonica_chr4_snp3<-fill_NA(japonica_chr4_snp3)
japonica_chr5_snp3<-fill_NA(japonica_chr5_snp3)
japonica_chr6_snp3<-fill_NA(japonica_chr6_snp3)
japonica_chr7_snp3<-fill_NA(japonica_chr7_snp3)
japonica_chr8_snp3<-fill_NA(japonica_chr8_snp3)
japonica_chr9_snp3<-fill_NA(japonica_chr9_snp3)
japonica_chr10_snp3<-fill_NA(japonica_chr10_snp3)
japonica_chr11_snp3<-fill_NA(japonica_chr11_snp3)
japonica_chr12_snp3<-fill_NA(japonica_chr12_snp3)

#Function to calculate genetic position: 
#(previous genetic position)+((current physical position - previous physical position)*(recombination rate at current position))

gen_pos <- function(SNP){
  SNP$pos <- NA
  SNP$pos[1]<-SNP$`SNP Start`[1]*SNP$rate[1]
  for(i in 1:nrow(SNP)){
    if(i>1){
      SNP$pos[i]<- SNP$pos[i-1] + (SNP$`SNP Start`[i] - SNP$`SNP Start`[i-1])*SNP$rate[i]
    }
  }
  print(SNP$pos)
}

#Create genetic maps: use fxn to calculate genetic position, map the genetic position by physical position
japonica_chr1_snp3$pos <- gen_pos(japonica_chr1_snp3)
plot(japonica_chr1_snp3$`SNP Start`, japonica_chr1_snp3$rate, type = "l")
ggplot(japonica_chr1_snp3, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
japonica_chr1_finalpos <- japonica_chr1_snp3[order(japonica_chr1_snp3$pos),]
is.unsorted(japonica_chr1_finalpos$pos)
plot(japonica_chr1_snp3$`SNP Start`, japonica_chr1_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Chromosome 1 Genetic Map")

japonica_chr2_snp3$pos <- gen_pos(japonica_chr2_snp3)
plot(japonica_chr2_snp3$`SNP Start`, japonica_chr2_snp3$pos)
japonica_chr2_finalpos <- japonica_chr2_snp3[order(japonica_chr2_snp3$pos),]
is.unsorted(japonica_chr2_finalpos$pos)
plot(japonica_chr2_snp3$`SNP Start`, japonica_chr2_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Chromosome 2 Genetic Map")

japonica_chr3_snp3$pos <- gen_pos(japonica_chr3_snp3)
plot(japonica_chr3_snp3$`SNP Start`, japonica_chr3_snp3$pos)
japonica_chr3_finalpos <- japonica_chr3_snp3[order(japonica_chr3_snp3$pos),]
is.unsorted(japonica_chr3_finalpos$pos)
plot(japonica_chr3_snp3$`SNP Start`, japonica_chr3_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Chromosome 3 Genetic Map")

japonica_chr4_spl <- smooth.spline(japonica_chr4_snp3$rate, spar = 0)
japonica_chr4_snp3$pos <- gen_pos(japonica_chr4_snp3)
plot(japonica_chr4_snp3$`SNP Start`, japonica_chr4_snp3$pos)
japonica_chr4_finalpos <- japonica_chr4_snp3[order(japonica_chr4_snp3$pos),]
is.unsorted(japonica_chr4_finalpos$pos)
plot(japonica_chr4_snp3$`SNP Start`, japonica_chr4_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Chromosome 4 Genetic Map")

japonica_chr5_snp3$pos <- gen_pos(japonica_chr5_snp3)
plot(japonica_chr5_snp3$`SNP Start`, japonica_chr5_snp3$pos)
japonica_chr5_finalpos <- japonica_chr5_snp3[order(japonica_chr5_snp3$pos),]
is.unsorted(japonica_chr5_finalpos$pos)
plot(japonica_chr5_snp3$`SNP Start`, japonica_chr5_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Chromosome 5 Genetic Map")

japonica_chr6_snp3$pos <- gen_pos(japonica_chr6_snp3)
plot(japonica_chr6_snp3$`SNP Start`, japonica_chr6_snp3$pos)
japonica_chr6_finalpos <- japonica_chr6_snp3[order(japonica_chr6_snp3$pos),]
is.unsorted(japonica_chr6_finalpos$pos)
plot(japonica_chr6_snp3$`SNP Start`, japonica_chr6_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Chromosome 6 Genetic Map")

japonica_chr7_snp3$pos <- gen_pos(japonica_chr7_snp3)
plot(japonica_chr7_snp3$`SNP Start`, japonica_chr7_snp3$pos)
japonica_chr7_finalpos <- japonica_chr7_snp3[order(japonica_chr7_snp3$pos),]
is.unsorted(japonica_chr7_finalpos$pos)
plot(japonica_chr7_snp3$`SNP Start`, japonica_chr7_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Chromosome 7 Genetic Map")

japonica_chr8_snp3$pos <- gen_pos(japonica_chr8_snp3)
plot(japonica_chr8_snp3$`SNP Start`, japonica_chr8_snp3$pos)
japonica_chr8_finalpos <- japonica_chr8_snp3[order(japonica_chr8_snp3$pos),]
is.unsorted(japonica_chr8_finalpos$pos)
plot(japonica_chr8_snp3$`SNP Start`, japonica_chr8_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Chromosome 8 Genetic Map")

japonica_chr9_snp3$pos <- gen_pos(japonica_chr9_snp3)
plot(japonica_chr9_snp3$`SNP Start`, japonica_chr9_snp3$pos)
japonica_chr9_finalpos <- japonica_chr9_snp3[order(japonica_chr9_snp3$pos),]
is.unsorted(japonica_chr9_finalpos$pos)
plot(japonica_chr9_snp3$`SNP Start`, japonica_chr9_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Chromosome 9 Genetic Map")

japonica_chr10_snp3$pos <- gen_pos(japonica_chr10_snp3)
plot(japonica_chr10_snp3$`SNP Start`, japonica_chr10_snp3$pos)
japonica_chr10_finalpos <- japonica_chr10_snp3[order(japonica_chr10_snp3$pos),]
is.unsorted(japonica_chr10_finalpos$pos)
plot(japonica_chr10_snp3$`SNP Start`, japonica_chr10_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Chromosome 10 Genetic Map")

japonica_chr11_snp3$pos <- gen_pos(japonica_chr11_snp3)
plot(japonica_chr11_snp3$`SNP Start`, japonica_chr11_snp3$pos)
japonica_chr11_finalpos <- japonica_chr11_snp3[order(japonica_chr11_snp3$pos),]
is.unsorted(japonica_chr11_finalpos$pos)
plot(japonica_chr11_snp3$`SNP Start`, japonica_chr11_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Chromosome 11 Genetic Map")

japonica_chr12_snp3$pos <- gen_pos(japonica_chr12_snp3)
plot(japonica_chr12_snp3$`SNP Start`, japonica_chr12_snp3$pos)
japonica_chr12_finalpos <- japonica_chr12_snp3[order(japonica_chr12_snp3$pos),]
is.unsorted(japonica_chr12_finalpos$pos)
plot(japonica_chr12_snp3$`SNP Start`, japonica_chr12_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Chromosome 12 Genetic Map")


#Create final WT Japonica genetic map
chr1 <- japonica_chr1_finalpos$pos/100
chr2 <- japonica_chr2_finalpos$pos/100
chr3 <- japonica_chr3_finalpos$pos/100
chr4 <- japonica_chr4_finalpos$pos/100
chr5 <- japonica_chr5_finalpos$pos/100
chr6 <- japonica_chr6_finalpos$pos/100
chr7 <- japonica_chr7_finalpos$pos/100
chr8 <- japonica_chr8_finalpos$pos/100
chr9 <- japonica_chr9_finalpos$pos/100
chr10<- japonica_chr10_finalpos$pos/100
chr11 <- japonica_chr11_finalpos$pos/100
chr12<- japonica_chr12_finalpos$pos/100

#Add genetic positions to a list (variable:japonica_map) and label each position by  chromosome and location
segSites<-readRDS("japonica_num_SNP.RData")
japonica_map = vector("list",10)
japonica_map[[1]] = chr1
japonica_map[[2]] = chr2
japonica_map[[3]] = chr3
japonica_map[[4]] = chr4
japonica_map[[5]] = chr5
japonica_map[[6]] = chr6
japonica_map[[7]] = chr7
japonica_map[[8]] = chr8
japonica_map[[9]] = chr9
japonica_map[[10]] = chr10
japonica_map[[11]] = chr11
japonica_map[[12]] = chr12
for(i in 1:12){
  names(japonica_map[[i]]) = paste(i, 1:segSites[i], sep="_")
}

saveRDS(japonica_map, file="japonica_final_map.RData")

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
c1 <-find_centromere(16.7,japonica_chr1_finalpos)
c2 <-find_centromere(13.6,japonica_chr2_finalpos)
c3 <-find_centromere(19.4,japonica_chr3_finalpos)
c4 <-find_centromere(9.7,japonica_chr4_finalpos)
c5 <-find_centromere(12.4,japonica_chr5_finalpos)
c6 <-find_centromere(15.3,japonica_chr6_finalpos)
c7 <-find_centromere(12.1,japonica_chr7_finalpos)
c8 <-find_centromere(12.9,japonica_chr8_finalpos)
c9 <-find_centromere(2.8,japonica_chr9_finalpos)
c10 <-find_centromere(8.2,japonica_chr10_finalpos)
c11 <-find_centromere(12,japonica_chr11_finalpos)
c12 <-find_centromere(11.9,japonica_chr12_finalpos)

japonica_centromere <- c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12)
japonica_centromere <- japonica_centromere/100

saveRDS(japonica_centromere, file="japonica_centromeres.RData")

saveRDS(japonica_chr1_finalpos, file="japonica_chr1_finalpos.RData")
saveRDS(japonica_chr2_finalpos, file="japonica_chr2_finalpos.RData")
saveRDS(japonica_chr3_finalpos, file="japonica_chr3_finalpos.RData")
saveRDS(japonica_chr4_finalpos, file="japonica_chr4_finalpos.RData")
saveRDS(japonica_chr5_finalpos, file="japonica_chr5_finalpos.RData")
saveRDS(japonica_chr6_finalpos, file="japonica_chr6_finalpos.RData")
saveRDS(japonica_chr7_finalpos, file="japonica_chr7_finalpos.RData")
saveRDS(japonica_chr8_finalpos, file="japonica_chr8_finalpos.RData")
saveRDS(japonica_chr9_finalpos, file="japonica_chr9_finalpos.RData")
saveRDS(japonica_chr10_finalpos, file="japonica_chr10_finalpos.RData")
saveRDS(japonica_chr11_finalpos, file="japonica_chr11_finalpos.RData")
saveRDS(japonica_chr12_finalpos, file="japonica_chr12_finalpos.RData")

#Create vector of SNP count on each chromosome
japonica_nSNP<-c(nrow(japonica_chr1_finalpos),nrow(japonica_chr2_finalpos),nrow(japonica_chr3_finalpos),nrow(japonica_chr4_finalpos),nrow(japonica_chr5_finalpos),nrow(japonica_chr6_finalpos),nrow(japonica_chr7_finalpos),nrow(japonica_chr8_finalpos),nrow(japonica_chr9_finalpos),nrow(japonica_chr10_finalpos),nrow(japonica_chr11_finalpos),nrow(japonica_chr12_finalpos))
saveRDS(japonica_nSNP, file="japonica_num_SNP.RData")
