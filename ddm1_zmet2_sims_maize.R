#reading in data frame to easily change recombination rate
ideal2_dist <- read.table("maize_genome_ddm1_zmet2.txt", header = FALSE)
colnames(ideal2_dist) <- c("Chr", "Start", "End", "Female WT", "Male WT", "ddm1_1", "ddm1_2", "zmet2")

#difference starts at 0
ideal2_dist$diffwt <- 0

chr1_ideal2 <- ideal2_dist[ which(ideal2_dist$Chr == 1),]
chr2_ideal2 <- ideal2_dist[ which(ideal2_dist$Chr == 2),]
chr3_ideal2 <- ideal2_dist[ which(ideal2_dist$Chr == 3),]
chr4_ideal2 <- ideal2_dist[ which(ideal2_dist$Chr == 4),]
chr5_ideal2 <- ideal2_dist[ which(ideal2_dist$Chr == 5),]
chr6_ideal2 <- ideal2_dist[ which(ideal2_dist$Chr == 6),]
chr7_ideal2 <- ideal2_dist[ which(ideal2_dist$Chr == 7),]
chr8_ideal2 <- ideal2_dist[ which(ideal2_dist$Chr == 8),]
chr9_ideal2 <- ideal2_dist[ which(ideal2_dist$Chr == 9),]
chr10_ideal2 <- ideal2_dist[ which(ideal2_dist$Chr == 10),]

#ideal2 & zmet2 double mutant
#changing chromosome positions based on hypothetical double mutant
chr1_ideal2[1:9,9] <- 2.346557
chr1_ideal2[16:44,9] <- 1
chr1_ideal2[50:60,9] <- 2.346557

chr2_ideal2[1:8,9] <- 2.346557
chr2_ideal2[13:34,9] <- 1
chr2_ideal2[39:47,9] <- 2.346557

chr3_ideal2[1:7,9] <- 2.346557
chr3_ideal2[13:33,9] <- 1
chr3_ideal2[39:46,9] <- 2.346557

chr4_ideal2[1:8,9] <- 2.346557
chr4_ideal2[13:35,9] <- 1
chr4_ideal2[40:48,9] <- 2.346557

chr5_ideal2[1:7,9] <- 2.346557
chr5_ideal2[12:31,9] <- 1
chr5_ideal2[36:43,9] <- 2.346557

chr6_ideal2[1:5,9] <- 2.346557
chr6_ideal2[9:24,9] <- 1
chr6_ideal2[28:33,9] <- 2.346557

chr7_ideal2[1:6,9] <- 2.346557
chr7_ideal2[10:25,9] <- 1
chr7_ideal2[29:35,9] <- 2.346557

chr8_ideal2[1:6,9] <- 2.346557
chr8_ideal2[10:25,9] <- 1
chr8_ideal2[29:35,9] <- 2.346557

chr9_ideal2[1:5,9] <- 2.346557
chr9_ideal2[10:22,9] <- 1
chr9_ideal2[26:31,9] <- 2.346557

chr10_ideal2[1:5,9] <- 2.346557
chr10_ideal2[8:21,9] <- 1
chr10_ideal2[24:29,9] <- 2.346557

#function to multiple recombination difference and add back to wildtype recombination rate in certain genomic regions
ideal2_wt <- function(chr_bin, ideal2_dist){
  for(i in 1:nrow(chr_bin)){
    for(k in 1:nrow(ideal2_dist)){
      if(isTRUE(chr_bin$foo.X1[i] >= ideal2_dist$Start[k] && chr_bin$foo.X2 <= ideal2_dist$End[k])){
        chr_bin$final[i] <- chr_bin$rate[i] + (chr_bin$rate[i]*ideal2_dist$diffwt[k])
      }
    }
  }
  return(chr_bin)
}

#applying that function
chr1_w_ideal2 <- ideal2_wt(chr1_bin, chr1_ideal2)
chr2_w_ideal2 <- ideal2_wt(chr2_bin, chr2_ideal2)
chr3_w_ideal2 <- ideal2_wt(chr3_bin, chr3_ideal2)
chr4_w_ideal2 <- ideal2_wt(chr4_bin, chr4_ideal2)
chr5_w_ideal2 <- ideal2_wt(chr5_bin, chr5_ideal2)
chr6_w_ideal2 <- ideal2_wt(chr6_bin, chr6_ideal2)
chr7_w_ideal2 <- ideal2_wt(chr7_bin, chr7_ideal2)
chr8_w_ideal2 <- ideal2_wt(chr8_bin, chr8_ideal2)
chr9_w_ideal2 <- ideal2_wt(chr9_bin, chr9_ideal2)
chr10_w_ideal2 <- ideal2_wt(chr10_bin, chr10_ideal2)

#function to assign recombination interval to SNPs
snp_rate_ideal2 <- function(chr_w_ideal2, chr_snp){
  for(i in 1:nrow(chr_snp)){
    for(k in 1:nrow(chr_w_ideal2)){
      if(isTRUE((chr_snp$`SNP Start`[i] >= chr_w_ideal2$foo.X1[k]) && (chr_snp$`SNP Start`[i] <= chr_w_ideal2$foo.X2[k]))){
        chr_snp$rate[i] <- chr_w_ideal2$final[k]
      }
    }
  }
  print(chr_snp)
}

#assigning rate, finding genetic positions of each SNP on each chromosome
chr1_snp2ideal2 <- snp_rate_ideal2(chr1_w_ideal2, chr1_snp)
chr1_snp2ideal2$`SNP Start`<- chr1_snp2ideal2$`SNP Start`/1000000
chr1_snp2ideal2 <- chr1_snp2ideal2[order(chr1_snp2ideal2$`SNP Start`),]
#creation of genetic positions from smoothed recombination rate
chr1_snp2ideal2$pos <- gen_pos(chr1_snp2ideal2)
chr1_snp2ideal2$pos2 <- graph_recomb(chr1_snp2ideal2)
#graph to look at Mb vs. cM along chromosome
plot(chr1_snp2$`SNP Start`, chr1_snp2$pos, type = "l", col = "blue", 
     main = "Chr1. Genetic Maps", xlab = "Physical Positions (Mb)", ylab = "Recombination rate (cM/Mb)")
chr1_finalpos <- chr1_snp2ideal2[order(chr1_snp2ideal2$pos),]
#want False to input into AlphaSimR
is.unsorted(chr1_finalpos$pos)

chr2_snp2ideal2 <- snp_rate_ideal2(chr2_w_ideal2, chr2_snp)
chr2_snp2ideal2$`SNP Start` <- chr2_snp2ideal2$`SNP Start`/1000000
chr2_snp2ideal2$pos <- gen_pos(chr2_snp2ideal2)
plot(chr2_snp2ideal2$`SNP Start`, chr2_snp2ideal2$pos)
chr2_finalpos <- chr2_snp2ideal2[order(chr2_snp2ideal2$pos),]
is.unsorted(chr2_finalpos$pos)

chr3_snp2ideal2 <- snp_rate_ideal2(chr3_w_ideal2, chr3_snp)
chr3_snp2ideal2$`SNP Start` <- chr3_snp2ideal2$`SNP Start`/1000000
chr3_snp2ideal2$pos <- gen_pos(chr3_snp2ideal2)
plot(chr3_snp2ideal2$`SNP Start`, chr3_snp2ideal2$pos)
chr3_finalpos <- chr3_snp2ideal2[order(chr3_snp2ideal2$pos),]
is.unsorted(chr3_finalpos$pos)

chr4_snp2ideal2 <- snp_rate_ideal2(chr4_w_ideal2, chr4_snp)
chr4_snp2ideal2$`SNP Start` <- chr4_snp2ideal2$`SNP Start`/1000000
chr4_snp2ideal2$pos <- gen_pos(chr4_snp2ideal2)
plot(chr4_snp2ideal2$`SNP Start`, chr4_snp2ideal2$pos)
chr4_finalpos <- chr4_snp2ideal2[order(chr4_snp2ideal2$pos),]
is.unsorted(chr4_finalpos$pos)

chr5_snp2ideal2 <- snp_rate_ideal2(chr5_w_ideal2, chr5_snp)
chr5_snp2ideal2$`SNP Start` <- chr5_snp2ideal2$`SNP Start`/1000000
chr5_snp2ideal2$pos <- gen_pos(chr5_snp2ideal2)
plot(chr5_snp2ideal2$`SNP Start`, chr5_snp2ideal2$pos)
chr5_finalpos <- chr5_snp2ideal2[order(chr5_snp2ideal2$pos),]
is.unsorted(chr5_finalpos$pos)

chr6_snp2ideal2 <- snp_rate_ideal2(chr6_w_ideal2, chr6_snp)
chr6_snp2ideal2$`SNP Start` <- chr6_snp2ideal2$`SNP Start`/1000000
chr6_snp2ideal2$pos <- gen_pos(chr6_snp2ideal2)
plot(chr6_snp2ideal2$`SNP Start`, chr6_snp2ideal2$pos, type = "l")
chr6_finalpos <- chr6_snp2ideal2[order(chr6_snp2ideal2$pos),]
is.unsorted(chr6_finalpos$pos)

chr7_snp2ideal2 <- snp_rate_ideal2(chr7_w_ideal2, chr7_snp)
chr7_snp2ideal2$`SNP Start` <- chr7_snp2ideal2$`SNP Start`/1000000
chr7_snp2ideal2$pos <- gen_pos(chr7_snp2ideal2)
plot(chr7_snp2ideal2$`SNP Start`, chr7_snp2ideal2$pos)
chr7_finalpos <- chr7_snp2ideal2[order(chr7_snp2ideal2$pos),]
is.unsorted(chr7_finalpos$pos)

chr8_snp2ideal2 <- snp_rate_ideal2(chr8_w_ideal2, chr8_snp)
chr8_snp2ideal2$`SNP Start` <- chr8_snp2ideal2$`SNP Start`/1000000
chr8_snp2ideal2$pos <- gen_pos(chr8_snp2ideal2)
plot(chr8_snp2ideal2$`SNP Start`, chr8_snp2ideal2$pos)
chr8_finalpos <- chr8_snp2ideal2[order(chr8_snp2ideal2$pos),]
is.unsorted(chr8_finalpos$pos)

chr9_snp2ideal2 <- snp_rate_ideal2(chr9_w_ideal2, chr9_snp)
chr9_snp2ideal2$`SNP Start` <- chr9_snp2ideal2$`SNP Start`/1000000
chr9_snp2ideal2$pos <- gen_pos(chr9_snp2ideal2)
plot(chr9_snp2ideal2$`SNP Start`, chr9_snp2ideal2$pos)
chr9_finalpos <- chr9_snp2ideal2[order(chr9_snp2ideal2$pos),]
is.unsorted(chr9_finalpos$pos)

chr10_snp2ideal2 <- snp_rate_ideal2(chr10_w_ideal2, chr10_snp)
chr10_snp2ideal2$`SNP Start` <- chr10_snp2ideal2$`SNP Start`/1000000
chr10_snp2ideal2$pos <- gen_pos(chr10_snp2ideal2)
plot(chr10_snp2ideal2$`SNP Start`, chr10_snp2ideal2$pos)
chr10_finalpos <- chr10_snp2ideal2[order(chr10_snp2ideal2$pos),]
is.unsorted(chr10_finalpos$pos)

#putting genetic map together with cM 
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

ideal2_map = vector("list",10)
ideal2_map[[1]] = chr1
ideal2_map[[2]] = chr2
ideal2_map[[3]] = chr3
ideal2_map[[4]] = chr4
ideal2_map[[5]] = chr5
ideal2_map[[6]] = chr6
ideal2_map[[7]] = chr7
ideal2_map[[8]] = chr8
ideal2_map[[9]] = chr9
ideal2_map[[10]] = chr10
for(i in 1:10){
  names(ideal2_map[[i]]) = paste(i, 1:segSites[i], sep="_")
}
saveRDS(ideal2_map, file = "ideal2_map.RData")
#Creating vector of centromere positions
ideal2_centromere <- c(414.7823, 342.8922, 211.6182, 206.9404, 112,
                       63.02667, 257.1550, 244.20528, 216.75260, 220.07058)
ideal2_centromere <- (ideal2_centromere/100)
