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
chr1_recq4_dist[1:15,9] <- 0.95
chr1_recq4_dist[45:60,9] <- 0.95

chr2_recq4_dist[1:12,9] <- 0.95
chr2_recq4_dist[35:47,9] <- 0.95

chr3_recq4_dist[1:12,9] <- 0.95
chr3_recq4_dist[34:46,9] <- 0.95

chr4_recq4_dist[1:12,9] <- 0.95
chr4_recq4_dist[36:48,9] <- 0.95

chr5_recq4_dist[1:11,9] <- 0.95
chr5_recq4_dist[32:43,9] <- 0.95

chr6_recq4_dist[1:8,9] <- 0.95
chr6_recq4_dist[25:33,9] <- 0.95

chr7_recq4_dist[1:9,9] <- 0.95
chr7_recq4_dist[26:35,9] <- 0.95

chr8_recq4_dist[1:9,9] <- 0.95
chr8_recq4_dist[26:35,9] <- 0.95

chr9_recq4_dist[1:9,9] <- 0.95
chr9_recq4_dist[23:31,9] <- 0.95

chr10_recq4_dist[1:7,9] <- 0.95
chr10_recq4_dist[22:29,9] <- 0.95

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


chr1_snp2 <- snp_rate_recq4(chr1_w_recq4, chr1_snp)
chr1_snp2$`SNP Start`<- chr1_snp2$`SNP Start`/1000000
chr1_snp2 <- chr1_snp2[order(chr1_snp2$`SNP Start`),]
#smoothing the recombination rate so transitions between bins are not so abrupt
chr1_spl2 <- smooth.spline(chr1_snp2$rate, spar = 1.2)
#creation of genetic positions from smoothed recombination rate
chr1_snp2$pos <- (chr1_snp2$`SNP Start`*chr1_spl2$y)
#graph to look at Mb vs. cM along chromosome
plot(chr1_snp2$`SNP Start`, chr1_snp2$pos, type = "l", col = "blue", 
     main = "Chr1. Genetic Maps", xlab = "Physical Positions (Mb)", ylab = "Recombination rate (cM/Mb)")
ggplot(chr1_snp2, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
#graph to look at Mb vs. cM/Mb to see recombination rate along chromosome
plot(chr1_snp2$`SNP Start`, chr1_snp2$pos/chr1_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 1 Recombination Distribution")
chr1_finalpos <- chr1_snp2[order(chr1_snp2$pos),]
#want False to input into AlphaSimR
is.unsorted(chr1_finalpos$pos)
#plot again to make sure it looks the same
plot(chr1_snp2$`SNP Start`, chr1_finalpos$pos/chr1_snp2$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 1 Recombination Distribution")
plot(chr1_finalpos$`SNP Start`, chr1_finalpos$pos)


chr2_snp2 <- snp_rate_recq4(chr2_w_recq4, chr2_snp)
chr2_snp2$`SNP Start` <- chr2_snp2$`SNP Start`/1000000
chr2_snp2 <- chr2_snp2[-(228:237),]
chr2_spl <- smooth.spline(chr2_snp2$rate, spar = 1.25)
chr2_snp2$pos <- (chr2_snp2$`SNP Start`*chr2_spl$y)
plot(chr2_snp2$`SNP Start`, chr2_snp2$pos)
plot(chr2_snp2$`SNP Start`, chr2_snp2$pos/chr2_snp2$`SNP Start`, type = "l")
chr2_finalpos <- chr2_snp2[order(chr2_snp2$pos),]
is.unsorted(chr2_finalpos$pos)
plot(chr2_snp2$`SNP Start`, chr2_finalpos$pos/chr2_snp2$`SNP Start`, type = "l")

chr3_snp2 <- snp_rate_recq4(chr3_w_recq4, chr3_snp)
chr3_snp2$`SNP Start` <- chr3_snp2$`SNP Start`/1000000
chr3_spl <- smooth.spline(chr3_snp2$rate, spar = 1.2)
chr3_snp2$pos <- (chr3_snp2$`SNP Start`*chr3_spl$y)
plot(chr3_snp2$`SNP Start`, chr3_snp2$pos)
plot(chr3_snp2$`SNP Start`, chr3_snp2$pos/chr3_snp2$`SNP Start`, type = "l")
chr3_finalpos <- chr3_snp2[order(chr3_snp2$pos),]
is.unsorted(chr3_finalpos$pos)
plot(chr3_snp2$`SNP Start`, chr3_finalpos$pos/chr3_snp2$`SNP Start`, type = "l")

chr4_snp2 <- snp_rate_recq4(chr4_w_recq4, chr4_snp)
chr4_snp2$`SNP Start` <- chr4_snp2$`SNP Start`/1000000
chr4_spl <- smooth.spline(chr4_snp2$rate, spar = 1.25)
chr4_snp2$pos <- (chr4_snp2$`SNP Start`*chr4_spl$y)
plot(chr4_snp2$`SNP Start`, chr4_snp2$pos)
plot(chr4_snp2$`SNP Start`, chr4_snp2$pos/chr4_snp2$`SNP Start`, type = "l")
chr4_finalpos <- chr4_snp2[order(chr4_snp2$pos),]
is.unsorted(chr4_finalpos$pos)
plot(chr4_snp2$`SNP Start`, chr4_finalpos$pos/chr4_snp2$`SNP Start`, type = "l")

chr5_snp2 <- snp_rate_recq4(chr5_w_recq4, chr5_snp)
chr5_snp2$`SNP Start` <- chr5_snp2$`SNP Start`/1000000
chr5_spl <- smooth.spline(chr5_snp2$rate, spar = 1.2)
chr5_snp2$pos <- (chr5_snp2$`SNP Start`*chr5_spl$y)
plot(chr5_snp2$`SNP Start`, chr5_snp2$pos)
plot(chr5_snp2$`SNP Start`, chr5_snp2$pos/chr5_snp2$`SNP Start`, type = "l")
chr5_finalpos <- chr5_snp2[order(chr5_snp2$pos),]
is.unsorted(chr5_finalpos$pos)
plot(chr5_snp2$`SNP Start`, chr5_finalpos$pos/chr5_snp2$`SNP Start`, type = "l")

#chr 6 is lowkey fuked up
chr6_snp2 <- snp_rate_recq4(chr6_w_recq4, chr6_snp)
chr6_snp2$`SNP Start` <- chr6_snp2$`SNP Start`/1000000
chr6_spl <- smooth.spline(chr6_snp2$rate, spar = 1)
chr6_snp2$pos <- (chr6_snp2$`SNP Start`*chr6_spl$y)
plot(chr6_snp2$`SNP Start`, chr6_snp2$pos)
plot(chr6_snp2$`SNP Start`, chr6_snp2$pos/chr6_snp2$`SNP Start`, type = "l")
chr6_finalpos <- chr6_snp2[order(chr6_snp2$pos),]
is.unsorted(chr6_finalpos$pos)
plot(chr6_snp2$`SNP Start`, chr6_finalpos$pos/chr6_snp2$`SNP Start`, type = "l")

chr7_snp2 <- snp_rate_recq4(chr7_w_recq4, chr7_snp)
chr7_snp2$`SNP Start` <- chr7_snp2$`SNP Start`/1000000
chr7_spl <- smooth.spline(chr7_snp2$rate, spar = 1.2)
chr7_snp2$pos <- (chr7_snp2$`SNP Start`*chr7_spl$y)
plot(chr7_snp2$`SNP Start`, chr7_snp2$pos)
plot(chr7_snp2$`SNP Start`, chr7_snp2$pos/chr7_snp2$`SNP Start`, type = "l")
chr7_finalpos <- chr7_snp2[order(chr7_snp2$pos),]
is.unsorted(chr7_finalpos$pos)
plot(chr7_snp2$`SNP Start`, chr7_finalpos$pos/chr7_snp2$`SNP Start`, type = "l")

chr8_snp2 <- snp_rate_recq4(chr8_w_recq4, chr8_snp)
chr8_snp2$`SNP Start` <- chr8_snp2$`SNP Start`/1000000
chr8_spl <- smooth.spline(chr8_snp2$rate, spar = 1.2)
chr8_snp2$pos <- (chr8_snp2$`SNP Start`*chr8_spl$y)
plot(chr8_snp2$`SNP Start`, chr8_snp2$pos)
plot(chr8_snp2$`SNP Start`, chr8_snp2$pos/chr8_snp2$`SNP Start`, type = "l")
chr8_finalpos <- chr8_snp2[order(chr8_snp2$pos),]
is.unsorted(chr8_finalpos$pos)
plot(chr8_snp2$`SNP Start`, chr8_finalpos$pos/chr8_snp2$`SNP Start`, type = "l")

chr9_snp2 <- snp_rate_recq4(chr9_w_recq4, chr9_snp)
chr9_snp2$`SNP Start` <- chr9_snp2$`SNP Start`/1000000
chr9_spl <- smooth.spline(chr9_snp2$rate, spar = 1.15)
chr9_snp2$pos <- (chr9_snp2$`SNP Start`*chr9_spl$y)
plot(chr9_snp2$`SNP Start`, chr9_snp2$pos)
plot(chr9_snp2$`SNP Start`, chr9_snp2$pos/chr9_snp2$`SNP Start`, type = "l")
chr9_finalpos <- chr9_snp2[order(chr9_snp2$pos),]
is.unsorted(chr9_finalpos$pos)
plot(chr9_snp2$`SNP Start`, chr9_finalpos$pos/chr9_snp2$`SNP Start`, type = "l")

chr10_snp2 <- snp_rate_recq4(chr10_w_recq4, chr10_snp)
chr10_snp2$`SNP Start` <- chr10_snp2$`SNP Start`/1000000
chr10_spl <- smooth.spline(chr10_snp2$rate, spar = 1.2)
chr10_snp2$pos <- (chr10_snp2$`SNP Start`*chr10_spl$y)
plot(chr10_snp2$`SNP Start`, chr10_snp2$pos)
plot(chr10_snp2$`SNP Start`, chr10_snp2$pos/chr10_snp2$`SNP Start`, type = "l")
chr10_finalpos <- chr10_snp2[order(chr10_snp2$pos),]
is.unsorted(chr10_finalpos$pos)
plot(chr10_snp2$`SNP Start`, chr10_finalpos$pos/chr10_snp2$`SNP Start`, type = "l")

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

recq4_map <- list(chr1[[1]], chr2[[1]], 
                 chr3[[1]], chr4[[1]], chr5[[1]], 
                 chr6[[1]], chr7[[1]], chr8[[1]], 
                 chr9[[1]], chr10[[1]])

#Creating vector of centromere positions
#change this
recq4_centromere <- c(176.093, 184.6176, 77, 85.56475, 189.9949,
                      59.8399, 77.5, 92.6381, 115.3071, 66.1423)
recq4_centromere <- (recq4_centromere/100)

recq4_sel2_gv <- matrix(data = NA, nrow = 100, ncol = 20)
recq4_sel3_gv <- matrix(data = NA, nrow = 100, ncol = 20)
recq4_sel4_gv <- matrix(data = NA, nrow = 100, ncol = 20)
recq4_sel5_gv <- matrix(data = NA, nrow = 100, ncol = 20)
recq4_sel6_gv <- matrix(data = NA, nrow = 100, ncol = 20)
recq4_sel7_gv <- matrix(data = NA, nrow = 100, ncol = 20)
recq4var1 <- matrix(data = NA, ncol = 2, nrow = 200)
recq4var2 <- matrix(data = NA, ncol = 2, nrow = 200)
recq4var3 <- matrix(data = NA, ncol = 2, nrow = 200)
recq4var4 <- matrix(data = NA, ncol = 2, nrow = 200)
recq4var5 <- matrix(data = NA, ncol = 2, nrow = 200)
recq4var6 <- matrix(data = NA, ncol = 2, nrow = 200)
recq4var7 <- matrix(data = NA, ncol = 2, nrow = 200)
#stronger selection coefficient
for(i in 1:100){
  SP$switchGenMap(recq4_map, centromere = recq4_centromere)
  pop <- randCross2(goodpop[49,], badpop[155,], nCrosses = 200, nProgeny = 1, simParam = SP)
  pop <- setPheno(pop = pop, h2 = 0.8, simParam = SP)
  recq4var1[i] <- varA(pop)
  
  pop1_sel <- selectInd(pop, nInd = 80, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_2 <- selectInd(pop1_sel, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_cross <- randCross2(pop1_sel1_2, goodpop[49,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel1_cross <- setPheno(pop1_sel1_cross, h2 = 0.8, simParam = SP)
  recq4var2[i] <- varA(pop1_sel1_2)
  
  recq4_sel2_gv[i,] <- gv(pop1_sel1_2)
  pop1_sel2 <- selectInd(pop1_sel1_cross, nInd = 20, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_2 <- selectInd(pop1_sel2, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_cross <- randCross2(pop1_sel2_2, goodpop[49,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel2_cross <- setPheno(pop1_sel2_cross, h2 = 0.8, simParam = SP)
  recq4var3[i] <- varA(pop1_sel2_2)
  
  recq4_sel3_gv[i,] <- gv(pop1_sel2_2)
  pop1_sel3 <- selectInd(pop1_sel2_cross, nInd = 20, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_2 <- selectInd(pop1_sel3, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_cross <- randCross2(pop1_sel3_2, goodpop[49,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel3_cross <- setPheno(pop1_sel3_cross, h2 = 0.8, simParam = SP)
  recq4var4[i] <- varA(pop1_sel3_2)
  
  recq4_sel4_gv[i,] <- gv(pop1_sel3_2)
  pop1_sel4 <- selectInd(pop1_sel3_cross, nInd = 20, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_2 <- selectInd(pop1_sel4, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_cross <- randCross2(pop1_sel4_2, goodpop[49,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel4_cross <- setPheno(pop1_sel4_cross, h2 = 0.8, simParam = SP)
  recq4var5[i] <- varA(pop1_sel4_2)
  
  recq4_sel5_gv[i,] <- gv(pop1_sel4_2)
  pop1_sel5 <- selectInd(pop1_sel4_cross, nInd = 20, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_2 <- selectInd(pop1_sel5, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_cross <- randCross2(pop1_sel5_2, goodpop[49,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel5_cross <- setPheno(pop1_sel5_cross, h2 = 0.8, simParam = SP)
  recq4var6[i] <- varA(pop1_sel5_2)
  
  recq4_sel6_gv[i,] <- gv(pop1_sel5_2)
}

