library(AlphaSimR)
library(Rcpp)
library(ggplot2)
library(dplyr)

setwd("C:/Users/16192/Documents/PNAS_Simulations")
set.seed(420)

##assigning frequency to SNPs based on recombination frequency in each bin
snp_rate_zmet2 <- function(chr_w_zmet2, chr_snp){
  for(i in 1:nrow(chr_snp)){
    for(k in 1:nrow(chr_w_zmet2)){
      if(isTRUE((chr_snp$`SNP Start`[i] >= chr_w_zmet2$foo.X1[k]) && (chr_snp$`SNP Start`[i] <= chr_w_zmet2$foo.X2[k]))){
        chr_snp$rate[i] <- chr_w_zmet2$final[k]
      }
    }
  }
  print(chr_snp)
}

##Using zmet2 recombination landscape--> 20% increase in COs

zmet2_dist <- read.table("maize_genome_ddm1_zmet2.txt", header = FALSE)
colnames(zmet2_dist) <- c("Chr", "Start", "End", "Female WT", "Male WT", "ddm1_1", "ddm1_2", "zmet2")
zmet2_dist$zmet2 <- zmet2_dist$zmet2*2/95
zmet2_dist$`Male WT`<- zmet2_dist$`Male WT`*2/135

#zmet2_dist$diffwt <- (zmet2_dist$zmet2-zmet2_dist$`Male WT`)/((zmet2_dist$zmet2+zmet2_dist$`Male WT`)/2)
#zmet2_dist$diffwt2 <- zmet2_dist$diffwt*1
#zmet2_dist$diffwt2[zmet2_dist$diffwt2 <= 0] <- 0
#zmet2_dist$diffwt2 <- abs(zmet2_dist$diffwt2)
zmet2_dist$diffwt <- 0

chr1_distzmet2 <- zmet2_dist[ which(zmet2_dist$Chr == 1),]
chr2_distzmet2 <- zmet2_dist[ which(zmet2_dist$Chr == 2),]
chr3_distzmet2 <- zmet2_dist[ which(zmet2_dist$Chr == 3),]
chr4_distzmet2 <- zmet2_dist[ which(zmet2_dist$Chr == 4),]
chr5_distzmet2 <- zmet2_dist[ which(zmet2_dist$Chr == 5),]
chr6_distzmet2 <- zmet2_dist[ which(zmet2_dist$Chr == 6),]
chr7_distzmet2 <- zmet2_dist[ which(zmet2_dist$Chr == 7),]
chr8_distzmet2 <- zmet2_dist[ which(zmet2_dist$Chr == 8),]
chr9_distzmet2 <- zmet2_dist[ which(zmet2_dist$Chr == 9),]
chr10_distzmet2 <- zmet2_dist[ which(zmet2_dist$Chr == 10),]

chr1_distzmet2[16:44,9] <- 1

chr2_distzmet2[13:34,9] <- 1

chr3_distzmet2[13:33,9] <- 1

chr4_distzmet2[13:35,9] <- 1

chr5_distzmet2[12:31,9] <- 1

chr6_distzmet2[9:24,9] <- 1

chr7_distzmet2[10:25,9] <- 1

chr8_distzmet2[10:25,9] <- 1

chr9_distzmet2[10:22,9] <- 1

chr10_distzmet2[8:21,9] <- 1

zmet2_wt <- function(chr_bin, zmet2_dist){
  for(i in 1:nrow(chr_bin)){
    for(k in 1:nrow(zmet2_dist)){
      if(isTRUE(chr_bin$foo.X1[i] >= zmet2_dist$Start[k] && chr_bin$foo.X2 <= zmet2_dist$End[k])){
        chr_bin$final[i] <- chr_bin$rate[i] + (chr_bin$rate[i]*zmet2_dist$diffwt[k])
      }
    }
  } 
  return(chr_bin)
}
chr1_w_zmet2 <- zmet2_wt(chr1_bin, chr1_distzmet2)
chr2_w_zmet2 <- zmet2_wt(chr2_bin, chr2_distzmet2)
chr3_w_zmet2 <- zmet2_wt(chr3_bin, chr3_distzmet2)
chr4_w_zmet2 <- zmet2_wt(chr4_bin, chr4_distzmet2)
chr5_w_zmet2 <- zmet2_wt(chr5_bin, chr5_distzmet2)
chr6_w_zmet2 <- zmet2_wt(chr6_bin, chr6_distzmet2)
chr7_w_zmet2 <- zmet2_wt(chr7_bin, chr7_distzmet2)
chr8_w_zmet2 <- zmet2_wt(chr8_bin, chr8_distzmet2)
chr9_w_zmet2 <- zmet2_wt(chr9_bin, chr9_distzmet2)
chr10_w_zmet2 <- zmet2_wt(chr10_bin, chr10_distzmet2)

#using function, converted SNP start to Mb to get cM/Mb for final genetic position
chr1_snp2zmet2 <- snp_rate_zmet2(chr1_w_zmet2, chr1_snp)
chr1_snp2zmet2$`SNP Start`<- chr1_snp2zmet2$`SNP Start`/1000000
chr1_snp2zmet2 <- chr1_snp2zmet2[order(chr1_snp2zmet2$`SNP Start`),]
#creation of genetic positions from smoothed recombination rate
chr1_snp2zmet2$pos <- gen_pos(chr1_snp2zmet2)
chr1_snp2zmet2$pos2 <- graph_recomb(chr1_snp2zmet2)
#graph to look at Mb vs. cM along chromosome
plot(chr1_snp2$`SNP Start`, chr1_snp2$pos)
#graph to look at Mb vs. cM/Mb to see recombination rate along chromosome
chr1_finalpos <- chr1_snp2[order(chr1_snp2$pos),]
#want False to input into AlphaSimR
is.unsorted(chr1_finalpos$pos)

chr2_snp2 <- snp_rate_zmet2(chr2_w_zmet2, chr2_snp)
chr2_snp2$`SNP Start` <- chr2_snp2$`SNP Start`/1000000
chr2_snp2$pos <- gen_pos(chr2_snp2)
plot(chr2_snp2$`SNP Start`, chr2_snp2$pos)
chr2_finalpos <- chr2_snp2[order(chr2_snp2$pos),]
is.unsorted(chr2_finalpos$pos)

chr3_snp2 <- snp_rate_zmet2(chr3_w_zmet2, chr3_snp)
chr3_snp2$`SNP Start` <- chr3_snp2$`SNP Start`/1000000
chr3_snp2$pos <- gen_pos(chr3_snp2)
plot(chr3_snp2$`SNP Start`, chr3_snp2$pos)
chr3_finalpos <- chr3_snp2[order(chr3_snp2$pos),]
is.unsorted(chr3_finalpos$pos)

chr4_snp2 <- snp_rate_zmet2(chr4_w_zmet2, chr4_snp)
chr4_snp2$`SNP Start` <- chr4_snp2$`SNP Start`/1000000
chr4_snp2$pos <- gen_pos(chr4_snp2)
plot(chr4_snp2$`SNP Start`, chr4_snp2$pos)
chr4_finalpos <- chr4_snp2[order(chr4_snp2$pos),]
is.unsorted(chr4_finalpos$pos)

chr5_snp2 <- snp_rate_zmet2(chr5_w_zmet2, chr5_snp)
chr5_snp2$`SNP Start` <- chr5_snp2$`SNP Start`/1000000
chr5_snp2$pos <- gen_pos(chr5_snp2)
plot(chr5_snp2$`SNP Start`, chr5_snp2$pos)
chr5_finalpos <- chr5_snp2[order(chr5_snp2$pos),]
is.unsorted(chr5_finalpos$pos)

chr6_snp2 <- snp_rate_zmet2(chr6_w_zmet2, chr6_snp)
chr6_snp2$`SNP Start` <- chr6_snp2$`SNP Start`/1000000
chr6_snp2$pos <- gen_pos(chr6_snp2)
plot(chr6_snp2$`SNP Start`, chr6_snp2$pos)
chr6_finalpos <- chr6_snp2[order(chr6_snp2$pos),]
is.unsorted(chr6_finalpos$pos)

chr7_snp2 <- snp_rate_zmet2(chr7_w_zmet2, chr7_snp)
chr7_snp2$`SNP Start` <- chr7_snp2$`SNP Start`/1000000
chr7_snp2$pos <- gen_pos(chr7_snp2)
plot(chr7_snp2$`SNP Start`, chr7_snp2$pos)
chr7_finalpos <- chr7_snp2[order(chr7_snp2$pos),]
is.unsorted(chr7_finalpos$pos)

chr8_snp2 <- snp_rate_zmet2(chr8_w_zmet2, chr8_snp)
chr8_snp2$`SNP Start` <- chr8_snp2$`SNP Start`/1000000
chr8_snp2$pos <- gen_pos(chr8_snp2)
plot(chr8_snp2$`SNP Start`, chr8_snp2$pos)
chr8_finalpos <- chr8_snp2[order(chr8_snp2$pos),]
is.unsorted(chr8_finalpos$pos)

chr9_snp2 <- snp_rate_zmet2(chr9_w_zmet2, chr9_snp)
chr9_snp2$`SNP Start` <- chr9_snp2$`SNP Start`/1000000
chr9_snp2$pos <- gen_pos(chr9_snp2)
plot(chr9_snp2$`SNP Start`, chr9_snp2$pos)
chr9_finalpos <- chr9_snp2[order(chr9_snp2$pos),]
is.unsorted(chr9_finalpos$pos)

chr10_snp2 <- snp_rate_zmet2(chr10_w_zmet2, chr10_snp)
chr10_snp2$`SNP Start` <- chr10_snp2$`SNP Start`/1000000
chr10_snp2$pos <- gen_pos(chr10_snp2)
plot(chr10_snp2$`SNP Start`, chr10_snp2$pos)
chr10_finalpos <- chr10_snp2[order(chr10_snp2$pos),]
is.unsorted(chr10_finalpos$pos)

#Putting the final genetic map together
chr1 <- chr1_finalpos$pos/100
chr1len <- length(chr1)
dim(chr1) <- c(chr1len,1)
chr1 <- list(chr1)

chr2 <- chr2_finalpos$pos/100
chr2len <- length(chr2)
dim(chr2) <- c(chr2len,1)
chr2 <- list(chr2)

chr3 <- chr3_finalpos$pos/100
chr3len <- length(chr3)
dim(chr3) <- c(chr3len,1)
chr3 <- list(chr3)

chr4 <- chr4_finalpos$pos/100
chr4len <- length(chr4)
dim(chr4) <- c(chr4len,1)
chr4 <- list(chr4)

chr10 <- chr10_finalpos$pos/100
chr10len <- length(chr10)
dim(chr10) <- c(chr10len,1)
chr10 <- list(chr10)

chr5 <- chr5_finalpos$pos/100
chr5len <- length(chr5)
dim(chr5) <- c(chr5len,1)
chr5 <- list(chr5)

chr6 <- chr6_finalpos$pos/100
chr6len <- length(chr6)
dim(chr6) <- c(chr6len,1)
chr6 <- list(chr6)

chr7 <- chr7_finalpos$pos/100
chr7len <- length(chr7)
dim(chr7) <- c(chr7len,1)
chr7 <- list(chr7)

chr8 <- chr8_finalpos$pos/100
chr8len <- length(chr8)
dim(chr8) <- c(chr8len,1)
chr8 <- list(chr8)

chr9 <- chr9_finalpos$pos/100
chr9len <- length(chr9)
dim(chr9) <- c(chr9len,1)
chr9 <- list(chr9)

zmet2_map <- list(chr1[[1]], chr2[[1]], 
                  chr3[[1]], chr4[[1]], chr5[[1]], 
                  chr6[[1]], chr7[[1]], chr8[[1]], 
                  chr9[[1]], chr10[[1]])

#saveRDS(zmet2_map, file="zmet2_map.RData")
zmet2_map <- readRDS("zmet2_map.RData")

#Creating vector of centromere positions for zmet2
zmet2_centromere <- c(163.7717, 113.67616, 69.41110, 67.51603, 112, 
                      20.881080, 78.00210, 73.49867, 68.812445, 67.53644)
zmet2_centromere <- zmet2_centromere/100
