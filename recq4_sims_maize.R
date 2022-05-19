library(AlphaSimR)
library(Rcpp)
library(ggplot2)
library(dplyr)

#setwd("C:/Users/16192/Documents/PNAS_Simulations")
set.seed(420)

#rice_snps <- read.table("japonica_SNPs.bed", header =FALSE)
#colnames(rice_snps) <- c("Chr#", "SNP Start", "SNP End")

rice_CO <- read.csv("jap_wt_rate.csv", header = TRUE)
rice_CO_recq4 <- read.csv("jap_mut_recq4.csv", header = TRUE)

colnames(rice_CO) <- c("Chr", "CO Start", "CO End", "rate_wt")
rice_CO <- rice_CO[order(rice_CO$Chr,rice_CO$`CO Start`),]
colnames(rice_CO_recq4) <- c("Chr", "CO Start", "CO End", "rate_recq4")
rice_CO_recq4 <- rice_CO_recq4[order(rice_CO_recq4$Chr,rice_CO_recq4$`CO Start`),]

all_rice <- cbind(rice_CO, rice_CO_recq4)
all_rice$diff <- abs(all_rice$rate_wt-all_rice$rate_recq4)/((all_rice$rate_wt+all_rice$rate_recq4)/2)
all_rice <- all_rice[-c(144,178,217),]
mean(all_rice$diff)

recq4_dist <- read.table("maize_genome_ddm1_zmet2.txt", header = FALSE)
colnames(recq4_dist) <- c("Chr", "Start", "End", "Female WT", "Male WT", "ddm1_1", "ddm1_2", "zmet2")

recq4_dist$diffwt <- 0

chr1_recq4_dist <- recq4_dist[ which(recq4_dist$Chr == 1),]
chr2_recq4_dist <- recq4_dist[ which(recq4_dist$Chr == 2),]
chr3_recq4_dist <- recq4_dist[ which(recq4_dist$Chr == 3),]
chr4_recq4_dist <- recq4_dist[ which(recq4_dist$Chr == 4),]
chr5_recq4_dist <- recq4_dist[ which(recq4_dist$Chr == 5),]
chr6_recq4_dist <- recq4_dist[ which(recq4_dist$Chr == 6),]
chr7_recq4_dist <- recq4_dist[ which(recq4_dist$Chr == 7),]
chr8_recq4_dist <- recq4_dist[ which(recq4_dist$Chr == 8),]
chr9_recq4_dist <- recq4_dist[ which(recq4_dist$Chr == 9),]
chr10_recq4_dist <- recq4_dist[ which(recq4_dist$Chr == 10),]

#Assigning change in recombination to chr arms
chr1_recq4_dist[1:15,9] <- 1.581055
chr1_recq4_dist[45:60,9] <- 1.581055

chr2_recq4_dist[1:12,9] <- 1.581055
chr2_recq4_dist[35:47,9] <- 1.581055

chr3_recq4_dist[1:12,9] <- 1.581055
chr3_recq4_dist[34:46,9] <- 1.581055

chr4_recq4_dist[1:12,9] <- 1.581055
chr4_recq4_dist[36:48,9] <- 1.581055

chr5_recq4_dist[1:11,9] <- 1.581055
chr5_recq4_dist[32:43,9] <- 1.581055

chr6_recq4_dist[1:8,9] <- 1.581055
chr6_recq4_dist[25:33,9] <- 1.581055

chr7_recq4_dist[1:9,9] <- 1.581055
chr7_recq4_dist[26:35,9] <- 1.581055

chr8_recq4_dist[1:9,9] <- 1.581055
chr8_recq4_dist[26:35,9] <- 1.581055

chr9_recq4_dist[1:9,9] <- 1.581055
chr9_recq4_dist[23:31,9] <- 1.581055

chr10_recq4_dist[1:7,9] <- 1.581055
chr10_recq4_dist[22:29,9] <- 1.581055

recq4_wt <- function(chr_bin, recq4_dist){
  for(i in 1:nrow(chr_bin)){
    for(k in 1:nrow(recq4_dist)){
      if(isTRUE(chr_bin$foo.X1[i] >= recq4_dist$Start[k] && chr_bin$foo.X2 <= recq4_dist$End[k])){
        chr_bin$final[i] <- chr_bin$rate[i] + (chr_bin$rate[i]*recq4_dist$diffwt[k])
      }
    }
  }
  return(chr_bin)
}
chr1_w_recq4 <- recq4_wt(chr1_bin, chr1_recq4_dist)
chr2_w_recq4 <- recq4_wt(chr2_bin, chr2_recq4_dist)
chr3_w_recq4 <- recq4_wt(chr3_bin, chr3_recq4_dist)
chr4_w_recq4 <- recq4_wt(chr4_bin, chr4_recq4_dist)
chr5_w_recq4 <- recq4_wt(chr5_bin, chr5_recq4_dist)
chr6_w_recq4 <- recq4_wt(chr6_bin, chr6_recq4_dist)
chr7_w_recq4 <- recq4_wt(chr7_bin, chr7_recq4_dist)
chr8_w_recq4 <- recq4_wt(chr8_bin, chr8_recq4_dist)
chr9_w_recq4 <- recq4_wt(chr9_bin, chr9_recq4_dist)
chr10_w_recq4 <- recq4_wt(chr10_bin, chr10_recq4_dist)

snp_rate_recq4 <- function(chr_w_recq4, chr_snp){
  for(i in 1:nrow(chr_snp)){
    for(k in 1:nrow(chr_w_recq4)){
      if(isTRUE((chr_snp$`SNP Start`[i] >= chr_w_recq4$foo.X1[k]) && (chr_snp$`SNP Start`[i] <= chr_w_recq4$foo.X2[k]))){
        chr_snp$rate[i] <- chr_w_recq4$final[k]
      }
    }
  }
  print(chr_snp)
}


chr1_snp2recq4 <- snp_rate_recq4(chr1_w_recq4, chr1_snp)
chr1_snp2recq4$`SNP Start`<- chr1_snp2recq4$`SNP Start`/1000000
chr1_snp2recq4 <- chr1_snp2recq4[order(chr1_snp2recq4$`SNP Start`),]
#creation of genetic positions from smoothed recombination rate
chr1_snp2recq4$pos <- gen_pos(chr1_snp2recq4)
chr1_snp2recq4$pos2 <- graph_recomb(chr1_snp2recq4)
#graph to look at Mb vs. cM along chromosome
plot(chr1_snp2$`SNP Start`, chr1_snp2$pos, type = "l", col = "blue", 
     main = "Chr1. Genetic Maps", xlab = "Physical Positions (Mb)", ylab = "Recombination rate (cM/Mb)")
chr1_finalpos <- chr1_snp2recq4[order(chr1_snp2recq4$pos),]
#want False to input into AlphaSimR
is.unsorted(chr1_finalpos$pos)

chr2_snp2 <- snp_rate_recq4(chr2_w_recq4, chr2_snp)
chr2_snp2$`SNP Start` <- chr2_snp2$`SNP Start`/1000000
chr2_snp2$pos <- gen_pos(chr2_snp2)
plot(chr2_snp2$`SNP Start`, chr2_snp2$pos)
chr2_finalpos <- chr2_snp2[order(chr2_snp2$pos),]
is.unsorted(chr2_finalpos$pos)

chr3_snp2 <- snp_rate_recq4(chr3_w_recq4, chr3_snp)
chr3_snp2$`SNP Start` <- chr3_snp2$`SNP Start`/1000000
chr3_snp2$pos <- gen_pos(chr3_snp2)
plot(chr3_snp2$`SNP Start`, chr3_snp2$pos)
chr3_finalpos <- chr3_snp2[order(chr3_snp2$pos),]
is.unsorted(chr3_finalpos$pos)

chr4_snp2 <- snp_rate_recq4(chr4_w_recq4, chr4_snp)
chr4_snp2$`SNP Start` <- chr4_snp2$`SNP Start`/1000000
chr4_snp2$pos <- gen_pos(chr4_snp2)
plot(chr4_snp2$`SNP Start`, chr4_snp2$pos)
chr4_finalpos <- chr4_snp2[order(chr4_snp2$pos),]
is.unsorted(chr4_finalpos$pos)

chr5_snp2 <- snp_rate_recq4(chr5_w_recq4, chr5_snp)
chr5_snp2$`SNP Start` <- chr5_snp2$`SNP Start`/1000000
chr5_snp2$pos <- gen_pos(chr5_snp2)
plot(chr5_snp2$`SNP Start`, chr5_snp2$pos)
chr5_finalpos <- chr5_snp2[order(chr5_snp2$pos),]
is.unsorted(chr5_finalpos$pos)

chr6_snp2 <- snp_rate_recq4(chr6_w_recq4, chr6_snp)
chr6_snp2$`SNP Start` <- chr6_snp2$`SNP Start`/1000000
chr6_snp2$pos <- gen_pos(chr6_snp2)
plot(chr6_snp2$`SNP Start`, chr6_snp2$pos)
chr6_finalpos <- chr6_snp2[order(chr6_snp2$pos),]
is.unsorted(chr6_finalpos$pos)

chr7_snp2 <- snp_rate_recq4(chr7_w_recq4, chr7_snp)
chr7_snp2$`SNP Start` <- chr7_snp2$`SNP Start`/1000000
chr7_snp2$pos <- gen_pos(chr7_snp2)
plot(chr7_snp2$`SNP Start`, chr7_snp2$pos)
chr7_finalpos <- chr7_snp2[order(chr7_snp2$pos),]
is.unsorted(chr7_finalpos$pos)

chr8_snp2 <- snp_rate_recq4(chr8_w_recq4, chr8_snp)
chr8_snp2$`SNP Start` <- chr8_snp2$`SNP Start`/1000000
chr8_snp2$pos <- gen_pos(chr8_snp2)
plot(chr8_snp2$`SNP Start`, chr8_snp2$pos)
chr8_finalpos <- chr8_snp2[order(chr8_snp2$pos),]
is.unsorted(chr8_finalpos$pos)

chr9_snp2 <- snp_rate_recq4(chr9_w_recq4, chr9_snp)
chr9_snp2$`SNP Start` <- chr9_snp2$`SNP Start`/1000000
chr9_snp2$pos <- gen_pos(chr9_snp2)
plot(chr9_snp2$`SNP Start`, chr9_snp2$pos)
chr9_finalpos <- chr9_snp2[order(chr9_snp2$pos),]
is.unsorted(chr9_finalpos$pos)

chr10_snp2 <- snp_rate_recq4(chr10_w_recq4, chr10_snp)
chr10_snp2$`SNP Start` <- chr10_snp2$`SNP Start`/1000000
chr10_snp2$pos <- gen_pos(chr10_snp2)
plot(chr10_snp2$`SNP Start`, chr10_snp2$pos)
chr10_finalpos <- chr10_snp2[order(chr10_snp2$pos),]
is.unsorted(chr10_finalpos$pos)

chr1 <- chr1_finalpos$pos/100

chr2 <- chr2_finalpos$pos/100

chr3 <- chr3_finalpos$pos/100

chr4 <- chr4_finalpos$pos/100

chr10 <- chr10_finalpos$pos/100

chr5 <- chr5_finalpos$pos/100

chr6 <- chr6_finalpos$pos/100

chr7 <- chr7_finalpos$pos/100

chr8 <- chr8_finalpos$pos/100

chr9 <- chr9_finalpos$pos/100

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
for(i in 1:10){
  names(recq4_map[[i]]) = paste(i, 1:segSites[i], sep="_")
}

saveRDS(recq4_map, file="recq4_map.RData")

#Creating vector of centromere positions
#change this
recq4_centromere <- c(355.4204, 285.0390, 168.4086, 162.7559, 239,
                      51.99450, 199.68464, 190.21911, 172.31881, 170.4290)
recq4_centromere <- (recq4_centromere/100)
