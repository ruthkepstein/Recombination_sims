library(AlphaSimR)
library(Rcpp)
library(ggplot2)
library(dplyr)

set.seed(420)

japonica_snps <- read.table("japonica_SNPs.bed", header =FALSE)
colnames(japonica_snps) <- c("Chr#", "SNP Start", "SNP End")
#sample SNPs?
fancm_snps <- sample_n(japonica_snps, 4000)
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


#fancm mutant recomb rates

fancm_CO <- read.csv("jap_mut_fancm.csv", header = TRUE)
colnames(fancm_CO) <- c("Chr", "CO Start", "CO End", "rate")
fancm_CO <- fancm_CO[order(fancm_CO$Chr,fancm_CO$`CO Start`),]

fancm_chr1_CO <- fancm_CO[ which(fancm_CO$Chr == "1"),]
fancm_chr1_CO$midpoint <- (fancm_chr1_CO$`CO Start`+ fancm_chr1_CO$`CO End`)/2
fancm_chr1_CO <- fancm_chr1_CO[order(fancm_chr1_CO$`CO Start`),]

fancm_chr2_CO <- fancm_CO[ which(fancm_CO$Chr == "2"),]
fancm_chr2_CO$midpoint <- (fancm_chr2_CO$`CO Start`+ fancm_chr2_CO$`CO End`)/2
fancm_chr2_CO <- fancm_chr2_CO[order(fancm_chr2_CO$`CO Start`),]

fancm_chr3_CO <- fancm_CO[ which(fancm_CO$Chr == "3"),]
fancm_chr3_CO$midpoint <- (fancm_chr3_CO$`CO Start`+ fancm_chr3_CO$`CO End`)/2
fancm_chr3_CO <- fancm_chr3_CO[order(fancm_chr3_CO$`CO Start`),]

fancm_chr4_CO <- fancm_CO[ which(fancm_CO$Chr == "4"),]
fancm_chr4_CO$midpoint <- (fancm_chr4_CO$`CO Start`+ fancm_chr4_CO$`CO End`)/2
fancm_chr4_CO <- fancm_chr4_CO[order(fancm_chr4_CO$`CO Start`),]

fancm_chr5_CO <- fancm_CO[ which(fancm_CO$Chr == "5"),]
fancm_chr5_CO$midpoint <- (fancm_chr5_CO$`CO Start`+ fancm_chr5_CO$`CO End`)/2
fancm_chr5_CO <- fancm_chr5_CO[order(fancm_chr5_CO$`CO Start`),]

fancm_chr6_CO <- fancm_CO[ which(fancm_CO$Chr == "6"),]
fancm_chr6_CO$midpoint <- (fancm_chr6_CO$`CO Start`+ fancm_chr6_CO$`CO End`)/2
fancm_chr6_CO <- fancm_chr6_CO[order(fancm_chr6_CO$`CO Start`),]

fancm_chr7_CO <- fancm_CO[ which(fancm_CO$Chr == "7"),]
fancm_chr7_CO$midpoint <- (fancm_chr7_CO$`CO Start`+ fancm_chr7_CO$`CO End`)/2
fancm_chr7_CO <- fancm_chr7_CO[order(fancm_chr7_CO$`CO Start`),]

fancm_chr8_CO <- fancm_CO[ which(fancm_CO$Chr == "8"),]
fancm_chr8_CO$midpoint <- (fancm_chr8_CO$`CO Start`+ fancm_chr8_CO$`CO End`)/2
fancm_chr8_CO <- fancm_chr8_CO[order(fancm_chr8_CO$`CO Start`),]

fancm_chr9_CO <- fancm_CO[ which(fancm_CO$Chr == "9"),]
fancm_chr9_CO$midpoint <- (fancm_chr9_CO$`CO Start`+ fancm_chr9_CO$`CO End`)/2
fancm_chr9_CO <- fancm_chr9_CO[order(fancm_chr9_CO$`CO Start`),]

fancm_chr10_CO <- fancm_CO[ which(fancm_CO$Chr == "10"),]
fancm_chr10_CO$midpoint <- (fancm_chr10_CO$`CO Start`+ fancm_chr10_CO$`CO End`)/2
fancm_chr10_CO <- fancm_chr10_CO[order(fancm_chr10_CO$`CO Start`),]

fancm_chr11_CO <- fancm_CO[ which(fancm_CO$Chr == "11"),]
fancm_chr11_CO$midpoint <- (fancm_chr11_CO$`CO Start`+ fancm_chr11_CO$`CO End`)/2
fancm_chr11_CO <- fancm_chr11_CO[order(fancm_chr11_CO$`CO Start`),]

fancm_chr12_CO <- fancm_CO[ which(fancm_CO$Chr == "12"),]
fancm_chr12_CO$midpoint <- (fancm_chr12_CO$`CO Start`+ fancm_chr12_CO$`CO End`)/2
fancm_chr12_CO <- fancm_chr12_CO[order(fancm_chr12_CO$`CO Start`),]


#calculating recombination rate per bin of CO data
library(dlookr)
library(tidyverse)
library(OneR)

#making intervals start at 0
fancm_chr1_CO$`CO Start` <- fancm_chr1_CO$`CO Start` - min(fancm_chr1_CO$`CO Start`)
fancm_chr1_CO$`CO End` <- fancm_chr1_CO$`CO End` - min(fancm_chr1_CO$`CO Start`)

fancm_chr2_CO$`CO Start` <- fancm_chr2_CO$`CO Start` - min(fancm_chr2_CO$`CO Start`)
fancm_chr2_CO$`CO End` <- fancm_chr2_CO$`CO End` - min(fancm_chr2_CO$`CO Start`)

fancm_chr3_CO$`CO Start` <- fancm_chr3_CO$`CO Start` - min(fancm_chr3_CO$`CO Start`)
fancm_chr3_CO$`CO End` <- fancm_chr3_CO$`CO End` - min(fancm_chr3_CO$`CO Start`)

fancm_chr4_CO$`CO Start` <- fancm_chr4_CO$`CO Start` - min(fancm_chr4_CO$`CO Start`)
fancm_chr4_CO$`CO End` <- fancm_chr4_CO$`CO End` - min(fancm_chr4_CO$`CO Start`)

fancm_chr5_CO$`CO Start` <- fancm_chr5_CO$`CO Start` - min(fancm_chr5_CO$`CO Start`)
fancm_chr5_CO$`CO End` <- fancm_chr5_CO$`CO End` - min(fancm_chr5_CO$`CO Start`)

fancm_chr6_CO$`CO Start` <- fancm_chr6_CO$`CO Start` - min(fancm_chr6_CO$`CO Start`)
fancm_chr6_CO$`CO End` <- fancm_chr6_CO$`CO End` - min(fancm_chr6_CO$`CO Start`)

fancm_chr7_CO$`CO Start` <- fancm_chr7_CO$`CO Start` - min(fancm_chr7_CO$`CO Start`)
fancm_chr7_CO$`CO End` <- fancm_chr7_CO$`CO End` - min(fancm_chr7_CO$`CO Start`)

fancm_chr8_CO$`CO Start` <- fancm_chr8_CO$`CO Start` - min(fancm_chr8_CO$`CO Start`)
fancm_chr8_CO$`CO End` <- fancm_chr8_CO$`CO End` - min(fancm_chr8_CO$`CO Start`)

fancm_chr9_CO$`CO Start` <- fancm_chr9_CO$`CO Start` - min(fancm_chr9_CO$`CO Start`)
fancm_chr9_CO$`CO End` <- fancm_chr9_CO$`CO End` - min(fancm_chr9_CO$`CO Start`)

fancm_chr10_CO$`CO Start` <- fancm_chr10_CO$`CO Start` - min(fancm_chr10_CO$`CO Start`)
fancm_chr10_CO$`CO End` <- fancm_chr10_CO$`CO End` - min(fancm_chr10_CO$`CO Start`)

fancm_chr11_CO$`CO Start` <- fancm_chr11_CO$`CO Start` - min(fancm_chr11_CO$`CO Start`)
fancm_chr11_CO$`CO End` <- fancm_chr11_CO$`CO End` - min(fancm_chr11_CO$`CO Start`)

fancm_chr12_CO$`CO Start` <- fancm_chr12_CO$`CO Start` - min(fancm_chr12_CO$`CO Start`)
fancm_chr12_CO$`CO End` <- fancm_chr12_CO$`CO End` - min(fancm_chr12_CO$`CO Start`)

isTRUE(fancm_chr1_CO[1,2] == 0)
isTRUE(fancm_chr2_CO[1,2] == 0)
isTRUE(fancm_chr3_CO[1,2] == 0)
isTRUE(fancm_chr4_CO[1,2] == 0)
isTRUE(fancm_chr5_CO[1,2] == 0)
isTRUE(fancm_chr6_CO[1,2] == 0)
isTRUE(fancm_chr7_CO[1,2] == 0)
isTRUE(fancm_chr8_CO[1,2] == 0)
isTRUE(fancm_chr9_CO[1,2] == 0)
isTRUE(fancm_chr10_CO[1,2] == 0)
isTRUE(fancm_chr11_CO[1,2] == 0)
isTRUE(fancm_chr12_CO[1,2] == 0)


##assigning frequency to SNPs based on recombination frequency in each bin
snp_rate <- function(chr_rate, chr_snp){
  for(i in 1:nrow(chr_snp)){
    for(k in 1:nrow(chr_rate)){
      if(isTRUE((chr_snp$`SNP Start`[i] >= chr_rate$`CO Start`[k]) && (chr_snp$`SNP Start`[i] <= chr_rate$`CO End`[k]))){
        chr_snp$rate[i] <- chr_rate$rate[k]
      }
    }
  }
}

#converted SNP start to Mb
fancm_chr1_snp$`SNP Start`<- fancm_chr1_snp$`SNP Start`/1000000
fancm_chr2_snp$`SNP Start` <- fancm_chr2_snp$`SNP Start`/1000000
fancm_chr3_snp$`SNP Start` <- fancm_chr3_snp$`SNP Start`/1000000
fancm_chr4_snp$`SNP Start` <- fancm_chr4_snp$`SNP Start`/1000000
fancm_chr5_snp$`SNP Start` <- fancm_chr5_snp$`SNP Start`/1000000
fancm_chr6_snp$`SNP Start` <- fancm_chr6_snp$`SNP Start`/1000000
fancm_chr7_snp$`SNP Start` <- fancm_chr7_snp$`SNP Start`/1000000
fancm_chr8_snp$`SNP Start` <- fancm_chr8_snp$`SNP Start`/1000000
fancm_chr9_snp$`SNP Start` <- fancm_chr9_snp$`SNP Start`/1000000
fancm_chr10_snp$`SNP Start` <- fancm_chr10_snp$`SNP Start`/1000000
fancm_chr11_snp$`SNP Start` <- fancm_chr11_snp$`SNP Start`/1000000
fancm_chr12_snp$`SNP Start` <- fancm_chr12_snp$`SNP Start`/1000000

#using function,  get cM/Mb for final genetic position - assign rates
fancm_chr1_snp2 <- snp_rate(fancm_chr1_CO, fancm_chr1_snp)
fancm_chr2_snp2 <- snp_rate(fancm_chr2_CO, fancm_chr2_snp)
fancm_chr3_snp2 <- snp_rate(fancm_chr3_CO, fancm_chr3_snp)
fancm_chr4_snp2 <- snp_rate(fancm_chr4_CO, fancm_chr4_snp)
fancm_chr5_snp2 <- snp_rate(fancm_chr5_CO, fancm_chr5_snp)
fancm_chr6_snp2 <- snp_rate(fancm_chr6_CO, fancm_chr6_snp)
fancm_chr7_snp2 <- snp_rate(fancm_chr7_CO, fancm_chr7_snp)
fancm_chr8_snp2 <- snp_rate(fancm_chr8_CO, fancm_chr8_snp)
fancm_chr9_snp2 <- snp_rate(fancm_chr9_CO, fancm_chr9_snp)
fancm_chr10_snp2 <- snp_rate(fancm_chr10_CO, fancm_chr10_snp)
fancm_chr11_snp2 <- snp_rate(fancm_chr11_CO, fancm_chr11_snp)
fancm_chr12_snp2 <- snp_rate(fancm_chr12_CO, fancm_chr12_snp)

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
#dataset too small to smooth...
fancm_chr1_snp2 <- fancm_chr1_snp2[order(fancm_chr1_snp2$`SNP Start`),]
fancm_chr1_spl <- smooth.spline(fancm_chr1_snp2$rate, spar = .7)
fancm_chr1_snp2$pos <- (fancm_chr1_snp2$`SNP Start`*fancm_chr1_spl$y)
plot(fancm_chr1_snp2$`SNP Start`, fancm_chr1_snp2$pos)
ggplot(fancm_chr1_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
plot(fancm_chr1_snp2$`SNP Start`, fancm_chr1_snp2$pos/fancm_chr1_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Fancm Chromosome 1 Recombination Distribution")
fancm_chr1_finalpos <- fancm_chr1_snp2[order(fancm_chr1_snp2$pos),]
is.unsorted(fancm_chr1_finalpos$pos)
plot(fancm_chr1_snp2$`SNP Start`, fancm_chr1_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Fancm Chromosome 1 Genetic Map")
plot(fancm_chr1_finalpos$`SNP Start`, fancm_chr1_finalpos$pos)


fancm_chr2_spl <- smooth.spline(fancm_chr2_snp2$rate, spar = .4)
fancm_chr2_snp2$pos <- (fancm_chr2_snp2$`SNP Start`*fancm_chr2_spl$y)
plot(fancm_chr2_snp2$`SNP Start`, fancm_chr2_snp2$pos)
plot(fancm_chr2_snp2$`SNP Start`, fancm_chr2_snp2$pos/fancm_chr2_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Fancm Chromosome 2 Recombination Distribution")
fancm_chr2_finalpos <- fancm_chr2_snp2[order(fancm_chr2_snp2$pos),]
is.unsorted(fancm_chr2_finalpos$pos)
plot(fancm_chr2_snp2$`SNP Start`, fancm_chr2_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Fancm Chromosome 2 Genetic Map")

fancm_chr3_spl <- smooth.spline(fancm_chr3_snp2$rate, spar = .3)
fancm_chr3_snp2$pos <- (fancm_chr3_snp2$`SNP Start`*fancm_chr3_spl$y)
plot(fancm_chr3_snp2$`SNP Start`, fancm_chr3_snp2$pos)
plot(fancm_chr3_snp2$`SNP Start`, fancm_chr3_snp2$pos/fancm_chr3_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Fancm Chromosome 3 Recombination Distribution")
fancm_chr3_finalpos <- fancm_chr3_snp2[order(fancm_chr3_snp2$pos),]
is.unsorted(fancm_chr3_finalpos$pos)
plot(fancm_chr3_snp2$`SNP Start`, fancm_chr3_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Fancm Chromosome 3 Genetic Map")

fancm_chr4_spl <- smooth.spline(fancm_chr4_snp2$rate, spar = .3)
fancm_chr4_snp2$pos <- (fancm_chr4_snp2$`SNP Start`*fancm_chr4_spl$y)
plot(fancm_chr4_snp2$`SNP Start`, fancm_chr4_snp2$pos)
plot(fancm_chr4_snp2$`SNP Start`, fancm_chr4_snp2$pos/fancm_chr4_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Fancm Chromosome 4 Recombination Distribution")
fancm_chr4_finalpos <- fancm_chr4_snp2[order(fancm_chr4_snp2$pos),]
is.unsorted(fancm_chr4_finalpos$pos)
plot(fancm_chr4_snp2$`SNP Start`, fancm_chr4_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Fancm Chromosome 4 Genetic Map")

fancm_chr5_spl <- smooth.spline(fancm_chr5_snp2$rate, spar =.4)
fancm_chr5_snp2$pos <- (fancm_chr5_snp2$`SNP Start`*fancm_chr5_spl$y)
plot(fancm_chr5_snp2$`SNP Start`, fancm_chr5_snp2$pos)
plot(fancm_chr5_snp2$`SNP Start`, fancm_chr5_snp2$pos/fancm_chr5_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Fancm Chromosome 5 Recombination Distribution")
fancm_chr5_finalpos <- fancm_chr5_snp2[order(fancm_chr5_snp2$pos),]
is.unsorted(fancm_chr5_finalpos$pos)
plot(fancm_chr5_snp2$`SNP Start`, fancm_chr5_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Fancm Chromosome 5 Genetic Map")

#slight increasing
fancm_chr6_spl <- smooth.spline(fancm_chr6_snp2$rate, spar = .289)
fancm_chr6_snp2$pos <- (fancm_chr6_snp2$`SNP Start`*fancm_chr6_spl$y)
plot(fancm_chr6_snp2$`SNP Start`, fancm_chr6_snp2$pos)
plot(fancm_chr6_snp2$`SNP Start`, fancm_chr6_snp2$pos/fancm_chr6_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Fancm Chromosome 6 Recombination Distribution")
fancm_chr6_finalpos <- fancm_chr6_snp2[order(fancm_chr6_snp2$pos),]
is.unsorted(fancm_chr6_finalpos$pos)
plot(fancm_chr6_snp2$`SNP Start`, fancm_chr6_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Fancm Chromosome 6 Genetic Map")

#slight increase
fancm_chr7_spl <- smooth.spline(fancm_chr7_snp2$rate, spar = 0.27)
fancm_chr7_snp2$pos <- (fancm_chr7_snp2$`SNP Start`*fancm_chr7_spl$y)
plot(fancm_chr7_snp2$`SNP Start`, fancm_chr7_snp2$pos)
plot(fancm_chr7_snp2$`SNP Start`, fancm_chr7_snp2$pos/fancm_chr7_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Fancm Chromosome 7 Recombination Distribution")
fancm_chr7_finalpos <- fancm_chr7_snp2[order(fancm_chr7_snp2$pos),]
is.unsorted(fancm_chr7_finalpos$pos)
plot(fancm_chr7_snp2$`SNP Start`, fancm_chr7_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Fancm Chromosome 7 Genetic Map")

fancm_chr8_spl <- smooth.spline(fancm_chr8_snp2$rate, spar = .4)
fancm_chr8_snp2$pos <- (fancm_chr8_snp2$`SNP Start`*fancm_chr8_spl$y)
plot(fancm_chr8_snp2$`SNP Start`, fancm_chr8_snp2$pos)
plot(fancm_chr8_snp2$`SNP Start`, fancm_chr8_snp2$pos/fancm_chr8_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Fancm Chromosome 8 Recombination Distribution")
fancm_chr8_finalpos <- fancm_chr8_snp2[order(fancm_chr8_snp2$pos),]
is.unsorted(fancm_chr8_finalpos$pos)
plot(fancm_chr8_snp2$`SNP Start`, fancm_chr8_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Fancm Chromosome 8 Genetic Map")

#slight increase
fancm_chr9_spl <- smooth.spline(fancm_chr9_snp2$rate, spar = .35)
fancm_chr9_snp2$pos <- (fancm_chr9_snp2$`SNP Start`*fancm_chr9_spl$y)
plot(fancm_chr9_snp2$`SNP Start`, fancm_chr9_snp2$pos)
plot(fancm_chr9_snp2$`SNP Start`, fancm_chr9_snp2$pos/fancm_chr9_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Fancm Chromosome 9 Recombination Distribution")
fancm_chr9_finalpos <- fancm_chr9_snp2[order(fancm_chr9_snp2$pos),]
is.unsorted(fancm_chr9_finalpos$pos)
plot(fancm_chr9_snp2$`SNP Start`, fancm_chr9_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Fancm Chromosome 9 Genetic Map")

#slight increase
fancm_chr10_spl <- smooth.spline(fancm_chr10_snp2$rate, spar =.29)
fancm_chr10_snp2$pos <- (fancm_chr10_snp2$`SNP Start`*fancm_chr10_spl$y)
plot(fancm_chr10_snp2$`SNP Start`, fancm_chr10_snp2$pos)
plot(fancm_chr10_snp2$`SNP Start`, fancm_chr10_snp2$pos/fancm_chr10_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Fancm Chromosome 10 Recombination Distribution")
fancm_chr10_finalpos <- fancm_chr10_snp2[order(fancm_chr10_snp2$pos),]
is.unsorted(fancm_chr10_finalpos$pos)
plot(fancm_chr10_snp2$`SNP Start`, fancm_chr10_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Fancm Chromosome 10 Genetic Map")

fancm_chr11_spl <- smooth.spline(fancm_chr11_snp2$rate, spar = .29)
fancm_chr11_snp2$pos <- (fancm_chr11_snp2$`SNP Start`*fancm_chr11_spl$y)
plot(fancm_chr11_snp2$`SNP Start`, fancm_chr11_snp2$pos)
plot(fancm_chr11_snp2$`SNP Start`, fancm_chr11_snp2$pos/fancm_chr11_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Fancm Chromosome 11 Recombination Distribution")
fancm_chr11_finalpos <- fancm_chr11_snp2[order(fancm_chr11_snp2$pos),]
is.unsorted(fancm_chr11_finalpos$pos)
plot(fancm_chr11_snp2$`SNP Start`, fancm_chr11_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Fancm Chromosome 11 Genetic Map")

fancm_chr12_spl <- smooth.spline(fancm_chr12_snp2$rate, spar = .29)
fancm_chr12_snp2$pos <- (fancm_chr12_snp2$`SNP Start`*fancm_chr12_spl$y)
plot(fancm_chr12_snp2$`SNP Start`, fancm_chr12_snp2$pos)
plot(fancm_chr12_snp2$`SNP Start`, fancm_chr12_snp2$pos/fancm_chr12_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Japonica Fancm Chromosome 12 Recombination Distribution")
fancm_chr12_finalpos <- fancm_chr12_snp2[order(fancm_chr12_snp2$pos),]
is.unsorted(fancm_chr12_finalpos$pos)
plot(fancm_chr12_snp2$`SNP Start`, fancm_chr12_finalpos$pos, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Genetic Position (cM)", main = "Japonica Fancm Chromosome 12 Genetic Map")
