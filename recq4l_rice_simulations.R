library(AlphaSimR)
library(Rcpp)
library(ggplot2)
library(dplyr)

set.seed(420)

japonica_snps <- read.table("japonica_SNPs.bed", header =FALSE)
colnames(japonica_snps) <- c("Chr#", "SNP Start", "SNP End")
#sample SNPs?
recq4l_snps <- sample_n(japonica_snps, 4000)
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


#recq4l mutant recomb rates

recq4l_CO <- read.csv("jap_mut_recq4l.csv", header = TRUE)
colnames(recq4l_CO) <- c("Chr", "CO Start", "CO End", "rate")
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
recq4l_chr1_CO$`CO Start` <- recq4l_chr1_CO$`CO Start` - min(recq4l_chr1_CO$`CO Start`)
recq4l_chr1_CO$`CO End` <- recq4l_chr1_CO$`CO End` - min(recq4l_chr1_CO$`CO Start`)

recq4l_chr2_CO$`CO Start` <- recq4l_chr2_CO$`CO Start` - min(recq4l_chr2_CO$`CO Start`)
recq4l_chr2_CO$`CO End` <- recq4l_chr2_CO$`CO End` - min(recq4l_chr2_CO$`CO Start`)

recq4l_chr3_CO$`CO Start` <- recq4l_chr3_CO$`CO Start` - min(recq4l_chr3_CO$`CO Start`)
recq4l_chr3_CO$`CO End` <- recq4l_chr3_CO$`CO End` - min(recq4l_chr3_CO$`CO Start`)

recq4l_chr4_CO$`CO Start` <- recq4l_chr4_CO$`CO Start` - min(recq4l_chr4_CO$`CO Start`)
recq4l_chr4_CO$`CO End` <- recq4l_chr4_CO$`CO End` - min(recq4l_chr4_CO$`CO Start`)

recq4l_chr5_CO$`CO Start` <- recq4l_chr5_CO$`CO Start` - min(recq4l_chr5_CO$`CO Start`)
recq4l_chr5_CO$`CO End` <- recq4l_chr5_CO$`CO End` - min(recq4l_chr5_CO$`CO Start`)

recq4l_chr6_CO$`CO Start` <- recq4l_chr6_CO$`CO Start` - min(recq4l_chr6_CO$`CO Start`)
recq4l_chr6_CO$`CO End` <- recq4l_chr6_CO$`CO End` - min(recq4l_chr6_CO$`CO Start`)

recq4l_chr7_CO$`CO Start` <- recq4l_chr7_CO$`CO Start` - min(recq4l_chr7_CO$`CO Start`)
recq4l_chr7_CO$`CO End` <- recq4l_chr7_CO$`CO End` - min(recq4l_chr7_CO$`CO Start`)

recq4l_chr8_CO$`CO Start` <- recq4l_chr8_CO$`CO Start` - min(recq4l_chr8_CO$`CO Start`)
recq4l_chr8_CO$`CO End` <- recq4l_chr8_CO$`CO End` - min(recq4l_chr8_CO$`CO Start`)

recq4l_chr9_CO$`CO Start` <- recq4l_chr9_CO$`CO Start` - min(recq4l_chr9_CO$`CO Start`)
recq4l_chr9_CO$`CO End` <- recq4l_chr9_CO$`CO End` - min(recq4l_chr9_CO$`CO Start`)

recq4l_chr10_CO$`CO Start` <- recq4l_chr10_CO$`CO Start` - min(recq4l_chr10_CO$`CO Start`)
recq4l_chr10_CO$`CO End` <- recq4l_chr10_CO$`CO End` - min(recq4l_chr10_CO$`CO Start`)

recq4l_chr11_CO$`CO Start` <- recq4l_chr11_CO$`CO Start` - min(recq4l_chr11_CO$`CO Start`)
recq4l_chr11_CO$`CO End` <- recq4l_chr11_CO$`CO End` - min(recq4l_chr11_CO$`CO Start`)

recq4l_chr12_CO$`CO Start` <- recq4l_chr12_CO$`CO Start` - min(recq4l_chr12_CO$`CO Start`)
recq4l_chr12_CO$`CO End` <- recq4l_chr12_CO$`CO End` - min(recq4l_chr12_CO$`CO Start`)

isTRUE(recq4l_chr1_CO[1,2] == 0)
isTRUE(recq4l_chr2_CO[1,2] == 0)
isTRUE(recq4l_chr3_CO[1,2] == 0)
isTRUE(recq4l_chr4_CO[1,2] == 0)
isTRUE(recq4l_chr5_CO[1,2] == 0)
isTRUE(recq4l_chr6_CO[1,2] == 0)
isTRUE(recq4l_chr7_CO[1,2] == 0)
isTRUE(recq4l_chr8_CO[1,2] == 0)
isTRUE(recq4l_chr9_CO[1,2] == 0)
isTRUE(recq4l_chr10_CO[1,2] == 0)
isTRUE(recq4l_chr11_CO[1,2] == 0)
isTRUE(recq4l_chr12_CO[1,2] == 0)


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
recq4l_chr1_snp$`SNP Start`<- recq4l_chr1_snp$`SNP Start`/1000000
recq4l_chr2_snp$`SNP Start` <- recq4l_chr2_snp$`SNP Start`/1000000
recq4l_chr3_snp$`SNP Start` <- recq4l_chr3_snp$`SNP Start`/1000000
recq4l_chr4_snp$`SNP Start` <- recq4l_chr4_snp$`SNP Start`/1000000
recq4l_chr5_snp$`SNP Start` <- recq4l_chr5_snp$`SNP Start`/1000000
recq4l_chr6_snp$`SNP Start` <- recq4l_chr6_snp$`SNP Start`/1000000
recq4l_chr7_snp$`SNP Start` <- recq4l_chr7_snp$`SNP Start`/1000000
recq4l_chr8_snp$`SNP Start` <- recq4l_chr8_snp$`SNP Start`/1000000
recq4l_chr9_snp$`SNP Start` <- recq4l_chr9_snp$`SNP Start`/1000000
recq4l_chr10_snp$`SNP Start` <- recq4l_chr10_snp$`SNP Start`/1000000
recq4l_chr11_snp$`SNP Start` <- recq4l_chr11_snp$`SNP Start`/1000000
recq4l_chr12_snp$`SNP Start` <- recq4l_chr12_snp$`SNP Start`/1000000

recq4l_chr1_snp$`SNP End`<- recq4l_chr1_snp$`SNP End`/1000000
recq4l_chr2_snp$`SNP End` <- recq4l_chr2_snp$`SNP End`/1000000
recq4l_chr3_snp$`SNP End` <- recq4l_chr3_snp$`SNP End`/1000000
recq4l_chr4_snp$`SNP End` <- recq4l_chr4_snp$`SNP End`/1000000
recq4l_chr5_snp$`SNP End` <- recq4l_chr5_snp$`SNP End`/1000000
recq4l_chr6_snp$`SNP End` <- recq4l_chr6_snp$`SNP End`/1000000
recq4l_chr7_snp$`SNP End` <- recq4l_chr7_snp$`SNP End`/1000000
recq4l_chr8_snp$`SNP End` <- recq4l_chr8_snp$`SNP End`/1000000
recq4l_chr9_snp$`SNP End` <- recq4l_chr9_snp$`SNP End`/1000000
recq4l_chr10_snp$`SNP End` <- recq4l_chr10_snp$`SNP End`/1000000
recq4l_chr11_snp$`SNP End` <- recq4l_chr11_snp$`SNP End`/1000000
recq4l_chr12_snp$`SNP End` <- recq4l_chr12_snp$`SNP End`/1000000

#using function,  get cM/Mb for final genetic position - assign rates
recq4l_chr1_snp2 <- snp_rate(recq4l_chr1_CO, recq4l_chr1_snp)
recq4l_chr2_snp2 <- snp_rate(recq4l_chr2_CO, recq4l_chr2_snp)
recq4l_chr3_snp2 <- snp_rate(recq4l_chr3_CO, recq4l_chr3_snp)
recq4l_chr4_snp2 <- snp_rate(recq4l_chr4_CO, recq4l_chr4_snp)
recq4l_chr5_snp2 <- snp_rate(recq4l_chr5_CO, recq4l_chr5_snp)
recq4l_chr6_snp2 <- snp_rate(recq4l_chr6_CO, recq4l_chr6_snp)
recq4l_chr7_snp2 <- snp_rate(recq4l_chr7_CO, recq4l_chr7_snp)
recq4l_chr8_snp2 <- snp_rate(recq4l_chr8_CO, recq4l_chr8_snp)
recq4l_chr9_snp2 <- snp_rate(recq4l_chr9_CO, recq4l_chr9_snp)
recq4l_chr10_snp2 <- snp_rate(recq4l_chr10_CO, recq4l_chr10_snp)
recq4l_chr11_snp2 <- snp_rate(recq4l_chr11_CO, recq4l_chr11_snp)
recq4l_chr12_snp2 <- snp_rate(recq4l_chr12_CO, recq4l_chr12_snp)

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
#dataset too small to smooth...
recq4l_chr1_snp2 <- recq4l_chr1_snp2[order(recq4l_chr1_snp2$`SNP Start`),]
recq4l_chr1_spl <- smooth.spline(recq4l_chr1_snp2$rate, spar = .7)
recq4l_chr1_snp2$pos <- (recq4l_chr1_snp2$`SNP Start`*recq4l_chr1_spl$y)
plot(recq4l_chr1_snp2$`SNP Start`, recq4l_chr1_snp2$pos)
ggplot(recq4l_chr1_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
plot(recq4l_chr1_snp2$`SNP Start`, recq4l_chr1_snp2$pos/recq4l_chr1_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Recq4l Chromosome 1 Recombination Distribution")
recq4l_chr1_finalpos <- recq4l_chr1_snp2[order(recq4l_chr1_snp2$pos),]
is.unsorted(recq4l_chr1_finalpos$pos)
plot(recq4l_chr1_snp2$`SNP Start`, recq4l_chr1_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Recq4l Chromosome 1 Genetic Map")
plot(recq4l_chr1_finalpos$`SNP Start`, recq4l_chr1_finalpos$pos)


recq4l_chr2_spl <- smooth.spline(recq4l_chr2_snp2$rate, spar = .4)
recq4l_chr2_snp2$pos <- (recq4l_chr2_snp2$`SNP Start`*recq4l_chr2_spl$y)
plot(recq4l_chr2_snp2$`SNP Start`, recq4l_chr2_snp2$pos)
plot(recq4l_chr2_snp2$`SNP Start`, recq4l_chr2_snp2$pos/recq4l_chr2_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Recq4l Chromosome 2 Recombination Distribution")
recq4l_chr2_finalpos <- recq4l_chr2_snp2[order(recq4l_chr2_snp2$pos),]
is.unsorted(recq4l_chr2_finalpos$pos)
plot(recq4l_chr2_snp2$`SNP Start`, recq4l_chr2_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Recq4l Chromosome 2 Genetic Map")

recq4l_chr3_spl <- smooth.spline(recq4l_chr3_snp2$rate, spar = .7)
recq4l_chr3_snp2$pos <- (recq4l_chr3_snp2$`SNP Start`*recq4l_chr3_spl$y)
plot(recq4l_chr3_snp2$`SNP Start`, recq4l_chr3_snp2$pos)
plot(recq4l_chr3_snp2$`SNP Start`, recq4l_chr3_snp2$pos/recq4l_chr3_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Recq4l Chromosome 3 Recombination Distribution")
recq4l_chr3_finalpos <- recq4l_chr3_snp2[order(recq4l_chr3_snp2$pos),]
is.unsorted(recq4l_chr3_finalpos$pos)
plot(recq4l_chr3_snp2$`SNP Start`, recq4l_chr3_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Recq4l Chromosome 3 Genetic Map")

recq4l_chr4_spl <- smooth.spline(recq4l_chr4_snp2$rate, spar = .7)
recq4l_chr4_snp2$pos <- (recq4l_chr4_snp2$`SNP Start`*recq4l_chr4_spl$y)
plot(recq4l_chr4_snp2$`SNP Start`, recq4l_chr4_snp2$pos)
plot(recq4l_chr4_snp2$`SNP Start`, recq4l_chr4_snp2$pos/recq4l_chr4_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Recq4l Chromosome 4 Recombination Distribution")
recq4l_chr4_finalpos <- recq4l_chr4_snp2[order(recq4l_chr4_snp2$pos),]
is.unsorted(recq4l_chr4_finalpos$pos)
plot(recq4l_chr4_snp2$`SNP Start`, recq4l_chr4_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Recq4l Chromosome 4 Genetic Map")

recq4l_chr5_spl <- smooth.spline(recq4l_chr5_snp2$rate, spar =.7)
recq4l_chr5_snp2$pos <- (recq4l_chr5_snp2$`SNP Start`*recq4l_chr5_spl$y)
plot(recq4l_chr5_snp2$`SNP Start`, recq4l_chr5_snp2$pos)
plot(recq4l_chr5_snp2$`SNP Start`, recq4l_chr5_snp2$pos/recq4l_chr5_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Recq4l Chromosome 5 Recombination Distribution")
recq4l_chr5_finalpos <- recq4l_chr5_snp2[order(recq4l_chr5_snp2$pos),]
is.unsorted(recq4l_chr5_finalpos$pos)
plot(recq4l_chr5_snp2$`SNP Start`, recq4l_chr5_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Recq4l Chromosome 5 Genetic Map")

#slight increasing
recq4l_chr6_spl <- smooth.spline(recq4l_chr6_snp2$rate, spar = 1)
recq4l_chr6_snp2$pos <- (recq4l_chr6_snp2$`SNP Start`*recq4l_chr6_spl$y)
plot(recq4l_chr6_snp2$`SNP Start`, recq4l_chr6_snp2$pos)
plot(recq4l_chr6_snp2$`SNP Start`, recq4l_chr6_snp2$pos/recq4l_chr6_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Recq4l Chromosome 6 Recombination Distribution")
recq4l_chr6_finalpos <- recq4l_chr6_snp2[order(recq4l_chr6_snp2$pos),]
is.unsorted(recq4l_chr6_finalpos$pos)
plot(recq4l_chr6_snp2$`SNP Start`, recq4l_chr6_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Recq4l Chromosome 6 Genetic Map")

#slight increase
recq4l_chr7_spl <- smooth.spline(recq4l_chr7_snp2$rate, spar = 0.9)
recq4l_chr7_snp2$pos <- (recq4l_chr7_snp2$`SNP Start`*recq4l_chr7_spl$y)
plot(recq4l_chr7_snp2$`SNP Start`, recq4l_chr7_snp2$pos)
plot(recq4l_chr7_snp2$`SNP Start`, recq4l_chr7_snp2$pos/recq4l_chr7_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Recq4l Chromosome 7 Recombination Distribution")
recq4l_chr7_finalpos <- recq4l_chr7_snp2[order(recq4l_chr7_snp2$pos),]
is.unsorted(recq4l_chr7_finalpos$pos)
plot(recq4l_chr7_snp2$`SNP Start`, recq4l_chr7_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Recq4l Chromosome 7 Genetic Map")

recq4l_chr8_spl <- smooth.spline(recq4l_chr8_snp2$rate, spar = .7)
recq4l_chr8_snp2$pos <- (recq4l_chr8_snp2$`SNP Start`*recq4l_chr8_spl$y)
plot(recq4l_chr8_snp2$`SNP Start`, recq4l_chr8_snp2$pos)
plot(recq4l_chr8_snp2$`SNP Start`, recq4l_chr8_snp2$pos/recq4l_chr8_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Recq4l Chromosome 8 Recombination Distribution")
recq4l_chr8_finalpos <- recq4l_chr8_snp2[order(recq4l_chr8_snp2$pos),]
is.unsorted(recq4l_chr8_finalpos$pos)
plot(recq4l_chr8_snp2$`SNP Start`, recq4l_chr8_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Recq4l Chromosome 8 Genetic Map")

#slight increase
recq4l_chr9_spl <- smooth.spline(recq4l_chr9_snp2$rate, spar = .8)
recq4l_chr9_snp2$pos <- (recq4l_chr9_snp2$`SNP Start`*recq4l_chr9_spl$y)
plot(recq4l_chr9_snp2$`SNP Start`, recq4l_chr9_snp2$pos)
plot(recq4l_chr9_snp2$`SNP Start`, recq4l_chr9_snp2$pos/recq4l_chr9_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Recq4l Chromosome 9 Recombination Distribution")
recq4l_chr9_finalpos <- recq4l_chr9_snp2[order(recq4l_chr9_snp2$pos),]
is.unsorted(recq4l_chr9_finalpos$pos)
plot(recq4l_chr9_snp2$`SNP Start`, recq4l_chr9_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Recq4l Chromosome 9 Genetic Map")

#slight increase
recq4l_chr10_spl <- smooth.spline(recq4l_chr10_snp2$rate, spar =.7)
recq4l_chr10_snp2$pos <- (recq4l_chr10_snp2$`SNP Start`*recq4l_chr10_spl$y)
plot(recq4l_chr10_snp2$`SNP Start`, recq4l_chr10_snp2$pos)
plot(recq4l_chr10_snp2$`SNP Start`, recq4l_chr10_snp2$pos/recq4l_chr10_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Recq4l Chromosome 10 Recombination Distribution")
recq4l_chr10_finalpos <- recq4l_chr10_snp2[order(recq4l_chr10_snp2$pos),]
is.unsorted(recq4l_chr10_finalpos$pos)
plot(recq4l_chr10_snp2$`SNP Start`, recq4l_chr10_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Recq4l Chromosome 10 Genetic Map")

recq4l_chr11_spl <- smooth.spline(recq4l_chr11_snp2$rate, spar = .7)
recq4l_chr11_snp2$pos <- (recq4l_chr11_snp2$`SNP Start`*recq4l_chr11_spl$y)
plot(recq4l_chr11_snp2$`SNP Start`, recq4l_chr11_snp2$pos)
plot(recq4l_chr11_snp2$`SNP Start`, recq4l_chr11_snp2$pos/recq4l_chr11_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Recq4l Chromosome 11 Recombination Distribution")
recq4l_chr11_finalpos <- recq4l_chr11_snp2[order(recq4l_chr11_snp2$pos),]
is.unsorted(recq4l_chr11_finalpos$pos)
plot(recq4l_chr11_snp2$`SNP Start`, recq4l_chr11_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Recq4l Chromosome 11 Genetic Map")

recq4l_chr12_spl <- smooth.spline(recq4l_chr12_snp2$rate, spar = .9)
recq4l_chr12_snp2$pos <- (recq4l_chr12_snp2$`SNP Start`*recq4l_chr12_spl$y)
plot(recq4l_chr12_snp2$`SNP Start`, recq4l_chr12_snp2$pos)
plot(recq4l_chr12_snp2$`SNP Start`, recq4l_chr12_snp2$pos/recq4l_chr12_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Recq4l Chromosome 12 Recombination Distribution")
recq4l_chr12_finalpos <- recq4l_chr12_snp2[order(recq4l_chr12_snp2$pos),]
is.unsorted(recq4l_chr12_finalpos$pos)
plot(recq4l_chr12_snp2$`SNP Start`, recq4l_chr12_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Recq4l Chromosome 12 Genetic Map")
