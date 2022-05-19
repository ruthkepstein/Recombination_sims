#reading in data frame to easily change recombination rate
ideal1_dist <- read.table("maize_genome_ddm1_zmet2.txt", header = FALSE)
colnames(ideal1_dist) <- c("Chr", "Start", "End", "Female WT", "Male WT", "ddm1_1", "ddm1_2", "zmet2")

#difference from wildtype starts at 0
ideal1_dist$diffwt <- 0

chr1_ideal1_dist <- ideal1_dist[ which(ideal1_dist$Chr == 1),]
chr2_ideal1_dist <- ideal1_dist[ which(ideal1_dist$Chr == 2),]
chr3_ideal1_dist <- ideal1_dist[ which(ideal1_dist$Chr == 3),]
chr4_ideal1_dist <- ideal1_dist[ which(ideal1_dist$Chr == 4),]
chr5_ideal1_dist <- ideal1_dist[ which(ideal1_dist$Chr == 5),]
chr6_ideal1_dist <- ideal1_dist[ which(ideal1_dist$Chr == 6),]
chr7_ideal1_dist <- ideal1_dist[ which(ideal1_dist$Chr == 7),]
chr8_ideal1_dist <- ideal1_dist[ which(ideal1_dist$Chr == 8),]
chr9_ideal1_dist <- ideal1_dist[ which(ideal1_dist$Chr == 9),]
chr10_ideal1_dist <- ideal1_dist[ which(ideal1_dist$Chr == 10),]

#changing chromosome positions based on hypothetical 10X increase
chr1_ideal1_dist[1:60,9] <- 10

chr2_ideal1_dist[1:47,9] <- 10

chr3_ideal1_dist[1:46,9] <- 10

chr4_ideal1_dist[1:48,9] <- 10

chr5_ideal1_dist[1:43,9] <- 10

chr6_ideal1_dist[1:33,9] <- 10

chr7_ideal1_dist[1:35,9] <- 10

chr8_ideal1_dist[1:35,9] <- 10

chr9_ideal1_dist[1:31,9] <- 10

chr10_ideal1_dist[1:29,9] <- 10

ideal1_wt <- function(chr_bin, ideal1_dist){
  for(i in 1:nrow(chr_bin)){
    for(k in 1:nrow(ideal1_dist)){
      if(isTRUE(chr_bin$foo.X1[i] >= ideal1_dist$Start[k] && chr_bin$foo.X2 <= ideal1_dist$End[k])){
        chr_bin$final[i] <- chr_bin$rate[i] + (chr_bin$rate[i]*ideal1_dist$diffwt[k])
      }
    }
  }
  return(chr_bin)
}
chr1_w_ideal1 <- ideal1_wt(chr1_bin, chr1_ideal1_dist)
chr2_w_ideal1 <- ideal1_wt(chr2_bin, chr2_ideal1_dist)
chr3_w_ideal1 <- ideal1_wt(chr3_bin, chr3_ideal1_dist)
chr4_w_ideal1 <- ideal1_wt(chr4_bin, chr4_ideal1_dist)
chr5_w_ideal1 <- ideal1_wt(chr5_bin, chr5_ideal1_dist)
chr6_w_ideal1 <- ideal1_wt(chr6_bin, chr6_ideal1_dist)
chr7_w_ideal1 <- ideal1_wt(chr7_bin, chr7_ideal1_dist)
chr8_w_ideal1 <- ideal1_wt(chr8_bin, chr8_ideal1_dist)
chr9_w_ideal1 <- ideal1_wt(chr9_bin, chr9_ideal1_dist)
chr10_w_ideal1 <- ideal1_wt(chr10_bin, chr10_ideal1_dist)

snp_rate_ideal1 <- function(chr_w_ideal1, chr_snp){
  for(i in 1:nrow(chr_snp)){
    for(k in 1:nrow(chr_w_ideal1)){
      if(isTRUE((chr_snp$`SNP Start`[i] >= chr_w_ideal1$foo.X1[k]) && (chr_snp$`SNP Start`[i] <= chr_w_ideal1$foo.X2[k]))){
        chr_snp$rate[i] <- chr_w_ideal1$final[k]
      }
    }
  }
  print(chr_snp)
}


chr1_snp2ideal1 <- snp_rate_ideal1(chr1_w_ideal1, chr1_snp)
chr1_snp2ideal1$`SNP Start`<- chr1_snp2ideal1$`SNP Start`/1000000
chr1_snp2ideal1 <- chr1_snp2ideal1[order(chr1_snp2ideal1$`SNP Start`),]
#creation of genetic positions from smoothed recombination rate
chr1_snp2ideal1$pos <- gen_pos(chr1_snp2ideal1)
chr1_snp2ideal1$pos2 <- graph_recomb(chr1_snp2ideal1)
plot(chr1_snp2$`SNP Start`, chr1_snp2$pos, type = "l")
plot(chr1_snp2ideal1$`SNP Start`, chr1_snp2ideal1$pos2/chr1_snp2ideal1$`SNP Start`, type = "l")
#graph to look at Mb vs. cM/Mb to see recombination rate along chromosome
chr1_finalpos <- chr1_snp2ideal1[order(chr1_snp2ideal1$pos),]
#want False to input into AlphaSimR
is.unsorted(chr1_finalpos$pos)

chr2_snp2 <- snp_rate_ideal1(chr2_w_ideal1, chr2_snp)
chr2_snp2$`SNP Start` <- chr2_snp2$`SNP Start`/1000000
chr2_snp2$pos <- gen_pos(chr2_snp2)
plot(chr2_snp2$`SNP Start`, chr2_snp2$pos)
chr2_finalpos <- chr2_snp2[order(chr2_snp2$pos),]
is.unsorted(chr2_finalpos$pos)

chr3_snp2 <- snp_rate_ideal1(chr3_w_ideal1, chr3_snp)
chr3_snp2$`SNP Start` <- chr3_snp2$`SNP Start`/1000000
chr3_snp2$pos <- gen_pos(chr3_snp2)
plot(chr3_snp2$`SNP Start`, chr3_snp2$pos)
chr3_finalpos <- chr3_snp2[order(chr3_snp2$pos),]
is.unsorted(chr3_finalpos$pos)

chr4_snp2 <- snp_rate_ideal1(chr4_w_ideal1, chr4_snp)
chr4_snp2$`SNP Start` <- chr4_snp2$`SNP Start`/1000000
chr4_snp2$pos <- gen_pos(chr4_snp2)
plot(chr4_snp2$`SNP Start`, chr4_snp2$pos)
chr4_finalpos <- chr4_snp2[order(chr4_snp2$pos),]
is.unsorted(chr4_finalpos$pos)

chr5_snp2 <- snp_rate_ideal1(chr5_w_ideal1, chr5_snp)
chr5_snp2$`SNP Start` <- chr5_snp2$`SNP Start`/1000000
chr5_snp2$pos <- gen_pos(chr5_snp2)
plot(chr5_snp2$`SNP Start`, chr5_snp2$pos)
chr5_finalpos <- chr5_snp2[order(chr5_snp2$pos),]
is.unsorted(chr5_finalpos$pos)

chr6_snp2 <- snp_rate_ideal1(chr6_w_ideal1, chr6_snp)
chr6_snp2$`SNP Start` <- chr6_snp2$`SNP Start`/1000000
chr6_snp2$pos <- gen_pos(chr6_snp2)
plot(chr6_snp2$`SNP Start`, chr6_snp2$pos)
chr6_finalpos <- chr6_snp2[order(chr6_snp2$pos),]
is.unsorted(chr6_finalpos$pos)

chr7_snp2 <- snp_rate_ideal1(chr7_w_ideal1, chr7_snp)
chr7_snp2$`SNP Start` <- chr7_snp2$`SNP Start`/1000000
chr7_snp2$pos <- gen_pos(chr7_snp2)
plot(chr7_snp2$`SNP Start`, chr7_snp2$pos)
chr7_finalpos <- chr7_snp2[order(chr7_snp2$pos),]
is.unsorted(chr7_finalpos$pos)

chr8_snp2 <- snp_rate_ideal1(chr8_w_ideal1, chr8_snp)
chr8_snp2$`SNP Start` <- chr8_snp2$`SNP Start`/1000000
chr8_snp2$pos <- gen_pos(chr8_snp2)
plot(chr8_snp2$`SNP Start`, chr8_snp2$pos)
chr8_finalpos <- chr8_snp2[order(chr8_snp2$pos),]
is.unsorted(chr8_finalpos$pos)

chr9_snp2 <- snp_rate_ideal1(chr9_w_ideal1, chr9_snp)
chr9_snp2$`SNP Start` <- chr9_snp2$`SNP Start`/1000000
chr9_snp2$pos <- gen_pos(chr9_snp2)
plot(chr9_snp2$`SNP Start`, chr9_snp2$pos)
chr9_finalpos <- chr9_snp2[order(chr9_snp2$pos),]
is.unsorted(chr9_finalpos$pos)

chr10_snp2 <- snp_rate_ideal1(chr10_w_ideal1, chr10_snp)
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

ideal1_map = vector("list",10)
ideal1_map[[1]] = chr1
ideal1_map[[2]] = chr2
ideal1_map[[3]] = chr3
ideal1_map[[4]] = chr4
ideal1_map[[5]] = chr5
ideal1_map[[6]] = chr6
ideal1_map[[7]] = chr7
ideal1_map[[8]] = chr8
ideal1_map[[9]] = chr9
ideal1_map[[10]] = chr10
for(i in 1:10){
  names(ideal1_map[[i]]) = paste(i, 1:segSites[i], sep="_")
}

saveRDS(ideal1_map, file="ideal1_map.RData")
#Creating vector of centromere positions

ideal1_centromere <- real_centromere*10
