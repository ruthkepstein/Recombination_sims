fancm <- read.table("fancm_rice_wt_mu.csv", header = TRUE, sep = ",")
fancm <- na.omit(fancm)
fancm$diff <- (fancm$r-fancm$wt_r)/fancm$wt_r
fancm <- na.omit(fancm)
fancm <- fancm[-which(fancm$diff == Inf),]
#fancm <- fancm[-which(fancm$diff2 == Inf),]
mean(fancm$diff)
sum(fancm$diff)
(sum(fancm$r) - sum(fancm$wt_r))/sum(fancm$wt_r)

fancm_dist <- read.table("maize_genome_ddm1_zmet2.txt", header = FALSE)
colnames(fancm_dist) <- c("Chr", "Start", "End", "Female WT", "Male WT", "fancm_1", "fancm_2", "zmet2")

fancm_dist$diffwt <- 0

chr1_fancm_dist <- fancm_dist[ which(fancm_dist$Chr == 1),]
chr2_fancm_dist <- fancm_dist[ which(fancm_dist$Chr == 2),]
chr3_fancm_dist <- fancm_dist[ which(fancm_dist$Chr == 3),]
chr4_fancm_dist <- fancm_dist[ which(fancm_dist$Chr == 4),]
chr5_fancm_dist <- fancm_dist[ which(fancm_dist$Chr == 5),]
chr6_fancm_dist <- fancm_dist[ which(fancm_dist$Chr == 6),]
chr7_fancm_dist <- fancm_dist[ which(fancm_dist$Chr == 7),]
chr8_fancm_dist <- fancm_dist[ which(fancm_dist$Chr == 8),]
chr9_fancm_dist <- fancm_dist[ which(fancm_dist$Chr == 9),]
chr10_fancm_dist <- fancm_dist[ which(fancm_dist$Chr == 10),]

#fancm dist
chr1_fancm_dist[1:15,9] <- 1.252264
chr1_fancm_dist[45:60,9] <- 1.252264

chr2_fancm_dist[1:12,9] <- 1.252264
chr2_fancm_dist[35:47,9] <- 1.252264

chr3_fancm_dist[1:12,9] <- 1.252264
chr3_fancm_dist[34:46,9] <- 1.252264

chr4_fancm_dist[1:12,9] <- 1.252264
chr4_fancm_dist[36:48,9] <- 1.252264

chr5_fancm_dist[1:11,9] <- 1.252264
chr5_fancm_dist[32:43,9] <- 1.252264

chr6_fancm_dist[1:8,9] <- 1.252264
chr6_fancm_dist[25:33,9] <- 1.252264

chr7_fancm_dist[1:9,9] <- 1.252264
chr7_fancm_dist[26:35,9] <- 1.252264

chr8_fancm_dist[1:9,9] <- 1.252264
chr8_fancm_dist[26:35,9] <- 1.252264

chr9_fancm_dist[1:9,9] <- 1.252264
chr9_fancm_dist[23:31,9] <- 1.252264

chr10_fancm_dist[1:7,9] <- 1.252264
chr10_fancm_dist[22:29,9] <- 1.252264

fancm_wt <- function(chr_bin, fancm_dist){
  for(i in 1:nrow(chr_bin)){
    for(k in 1:nrow(fancm_dist)){
      if(isTRUE(chr_bin$foo.X1[i] >= fancm_dist$Start[k] && chr_bin$foo.X2 <= fancm_dist$End[k])){
        chr_bin$final[i] <- chr_bin$rate[i] + (chr_bin$rate[i]*fancm_dist$diffwt[k])
      }
    }
  }
  return(chr_bin)
}
chr1_w_fancm <- fancm_wt(chr1_bin, chr1_fancm_dist)
chr2_w_fancm <- fancm_wt(chr2_bin, chr2_fancm_dist)
chr3_w_fancm <- fancm_wt(chr3_bin, chr3_fancm_dist)
chr4_w_fancm <- fancm_wt(chr4_bin, chr4_fancm_dist)
chr5_w_fancm <- fancm_wt(chr5_bin, chr5_fancm_dist)
chr6_w_fancm <- fancm_wt(chr6_bin, chr6_fancm_dist)
chr7_w_fancm <- fancm_wt(chr7_bin, chr7_fancm_dist)
chr8_w_fancm <- fancm_wt(chr8_bin, chr8_fancm_dist)
chr9_w_fancm <- fancm_wt(chr9_bin, chr9_fancm_dist)
chr10_w_fancm <- fancm_wt(chr10_bin, chr10_fancm_dist)

snp_rate_fancm <- function(chr_w_fancm, chr_snp){
  for(i in 1:nrow(chr_snp)){
    for(k in 1:nrow(chr_w_fancm)){
      if(isTRUE((chr_snp$`SNP Start`[i] >= chr_w_fancm$foo.X1[k]) && (chr_snp$`SNP Start`[i] <= chr_w_fancm$foo.X2[k]))){
        chr_snp$rate[i] <- chr_w_fancm$final[k]
      }
    }
  }
  print(chr_snp)
}


chr1_snp2fancm <- snp_rate_fancm(chr1_w_fancm, chr1_snp)
chr1_snp2fancm$`SNP Start`<- chr1_snp2fancm$`SNP Start`/1000000
chr1_snp2fancm <- chr1_snp2fancm[order(chr1_snp2fancm$`SNP Start`),]
#creation of genetic positions from smoothed recombination rate
chr1_snp2fancm$pos <- gen_pos(chr1_snp2fancm)
chr1_snp2fancm$pos2 <- graph_recomb(chr1_snp2fancm)
plot(chr1_snp2$`SNP Start`, chr1_snp2$pos, type = "l")
#graph to look at Mb vs. cM/Mb to see recombination rate along chromosome
chr1_finalpos <- chr1_snp2fancm[order(chr1_snp2fancm$pos),]
#want False to input into AlphaSimR
is.unsorted(chr1_finalpos$pos)

chr2_snp2 <- snp_rate_fancm(chr2_w_fancm, chr2_snp)
chr2_snp2$`SNP Start` <- chr2_snp2$`SNP Start`/1000000
chr2_snp2$pos <- gen_pos(chr2_snp2)
plot(chr2_snp2$`SNP Start`, chr2_snp2$pos, type = "l")
chr2_finalpos <- chr2_snp2[order(chr2_snp2$pos),]
is.unsorted(chr2_finalpos$pos)

chr3_snp2 <- snp_rate_fancm(chr3_w_fancm, chr3_snp)
chr3_snp2$`SNP Start` <- chr3_snp2$`SNP Start`/1000000
chr3_snp2$pos <- gen_pos(chr3_snp2)
plot(chr3_snp2$`SNP Start`, chr3_snp2$pos)
chr3_finalpos <- chr3_snp2[order(chr3_snp2$pos),]
is.unsorted(chr3_finalpos$pos)

chr4_snp2 <- snp_rate_fancm(chr4_w_fancm, chr4_snp)
chr4_snp2$`SNP Start` <- chr4_snp2$`SNP Start`/1000000
chr4_snp2$pos <- gen_pos(chr4_snp2)
plot(chr4_snp2$`SNP Start`, chr4_snp2$pos, type = "l")
chr4_finalpos <- chr4_snp2[order(chr4_snp2$pos),]
is.unsorted(chr4_finalpos$pos)

chr5_snp2 <- snp_rate_fancm(chr5_w_fancm, chr5_snp)
chr5_snp2$`SNP Start` <- chr5_snp2$`SNP Start`/1000000
chr5_snp2$pos <- gen_pos(chr5_snp2)
plot(chr5_snp2$`SNP Start`, chr5_snp2$pos)
chr5_finalpos <- chr5_snp2[order(chr5_snp2$pos),]
is.unsorted(chr5_finalpos$pos)

chr6_snp2 <- snp_rate_fancm(chr6_w_fancm, chr6_snp)
chr6_snp2$`SNP Start` <- chr6_snp2$`SNP Start`/1000000
chr6_snp2$pos <- gen_pos(chr6_snp2)
plot(chr6_snp2$`SNP Start`, chr6_snp2$pos)
chr6_finalpos <- chr6_snp2[order(chr6_snp2$pos),]
is.unsorted(chr6_finalpos$pos)

chr7_snp2 <- snp_rate_fancm(chr7_w_fancm, chr7_snp)
chr7_snp2$`SNP Start` <- chr7_snp2$`SNP Start`/1000000
chr7_snp2$pos <- gen_pos(chr7_snp2)
plot(chr7_snp2$`SNP Start`, chr7_snp2$pos, type = "l")
chr7_finalpos <- chr7_snp2[order(chr7_snp2$pos),]
is.unsorted(chr7_finalpos$pos)

chr8_snp2 <- snp_rate_fancm(chr8_w_fancm, chr8_snp)
chr8_snp2$`SNP Start` <- chr8_snp2$`SNP Start`/1000000
chr8_snp2$pos <- gen_pos(chr8_snp2)
plot(chr8_snp2$`SNP Start`, chr8_snp2$pos, type = "l")
chr8_finalpos <- chr8_snp2[order(chr8_snp2$pos),]
is.unsorted(chr8_finalpos$pos)

chr9_snp2 <- snp_rate_fancm(chr9_w_fancm, chr9_snp)
chr9_snp2$`SNP Start` <- chr9_snp2$`SNP Start`/1000000
chr9_snp2$pos <- gen_pos(chr9_snp2)
plot(chr9_snp2$`SNP Start`, chr9_snp2$pos)
chr9_finalpos <- chr9_snp2[order(chr9_snp2$pos),]
is.unsorted(chr9_finalpos$pos)

chr10_snp2 <- snp_rate_fancm(chr10_w_fancm, chr10_snp)
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

fancm_map = vector("list",10)
fancm_map[[1]] = chr1
fancm_map[[2]] = chr2
fancm_map[[3]] = chr3
fancm_map[[4]] = chr4
fancm_map[[5]] = chr5
fancm_map[[6]] = chr6
fancm_map[[7]] = chr7
fancm_map[[8]] = chr8
fancm_map[[9]] = chr9
fancm_map[[10]] = chr10
for(i in 1:10){
  names(fancm_map[[i]]) = paste(i, 1:segSites[i], sep="_")
}

saveRDS(fancm_map, file="fancm_map.RData")

#Creating vector of centromere positions
fancm_centromere <- c(312.2040, 248.9850, 147.2845, 142.3752, 210,
                      45.38585, 174.29785, 165.53841, 150.52966, 148.83764)
fancm_centromere <- (fancm_centromere/100)
