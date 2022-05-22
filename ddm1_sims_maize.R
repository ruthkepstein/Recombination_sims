##assigning frequency to SNPs based on recombination frequency in each bin
snp_rate_ddm1 <- function(chr_w_ddm1, chr_snp){
  for(i in 1:nrow(chr_snp)){
    for(k in 1:nrow(chr_w_ddm1)){
      if(isTRUE((chr_snp$`SNP Start`[i] >= chr_w_ddm1$foo.X1[k]) && (chr_snp$`SNP Start`[i] <= chr_w_ddm1$foo.X2[k]))){
        chr_snp$rate[i] <- chr_w_ddm1$final[k]
      }
    }
  }
  print(chr_snp)
}

##Using ddm1 recombination landscape--> 6% increase in COs
ddm1_dist <- read.table("maize_genome_ddm1_zmet2.txt", header = FALSE)
colnames(ddm1_dist) <- c("Chr", "Start", "End", "Female WT", "Male WT", "ddm1_1", "ddm1_2", "zmet2")
#normalizing the data
ddm1_dist$ddm1_1 <- ddm1_dist$ddm1_1*2/39
ddm1_dist$ddm1_2 <- ddm1_dist$ddm1_2*2/40
#ddm1_dist$`Female WT` <- ddm1_dist$`Female WT`*2/122
ddm1_dist$`Male WT`<- ddm1_dist$`Male WT`*2/135
#ddm1_dist$WT <- (ddm1_dist$`Male WT`)/2
ddm1_dist$ddm1 <- (ddm1_dist$ddm1_1+ddm1_dist$ddm1_2)/2

ddm1_dist$diffwt <- 0
#ddm1_dist$diffwt2 <- ddm1_dist$diffwt*1
#ddm1_dist$diffwt2[ddm1_dist$diffwt2 <= 0] <- 0
#ddm1_dist$diffwt2 <- abs(ddm1_dist$diffwt2)
#~160% increase in telomeres

#make generalization about ddm1

chr1_distddm1 <- ddm1_dist[ which(ddm1_dist$Chr == 1),]
chr2_distddm1 <- ddm1_dist[ which(ddm1_dist$Chr == 2),]
chr3_distddm1 <- ddm1_dist[ which(ddm1_dist$Chr == 3),]
chr4_distddm1 <- ddm1_dist[ which(ddm1_dist$Chr == 4),]
chr5_distddm1 <- ddm1_dist[ which(ddm1_dist$Chr == 5),]
chr6_distddm1 <- ddm1_dist[ which(ddm1_dist$Chr == 6),]
chr7_distddm1 <- ddm1_dist[ which(ddm1_dist$Chr == 7),]
chr8_distddm1 <- ddm1_dist[ which(ddm1_dist$Chr == 8),]
chr9_distddm1 <- ddm1_dist[ which(ddm1_dist$Chr == 9),]
chr10_distddm1 <- ddm1_dist[ which(ddm1_dist$Chr == 10),]

chr1_distddm1[1:9,10] <- 2.346557
chr1_distddm1[50:60,10] <- 2.346557

chr2_distddm1[1:8,10] <- 2.346557
chr2_distddm1[39:47,10] <- 2.346557

chr3_distddm1[1:7,10] <- 2.346557
chr3_distddm1[39:46,10] <- 2.346557

chr4_distddm1[1:8,10] <- 2.346557
chr4_distddm1[40:48,10] <- 2.346557

chr5_distddm1[1:7,10] <- 2.346557
chr5_distddm1[36:43,10] <- 2.346557

chr6_distddm1[1:5,10] <- 2.346557
chr6_distddm1[28:33,10] <- 2.346557

chr7_distddm1[1:6,10] <- 2.346557
chr7_distddm1[29:35,10] <- 2.346557

chr8_distddm1[1:6,10] <- 2.346557
chr8_distddm1[29:35,10] <- 2.346557

chr9_distddm1[1:5,10] <- 2.346557
chr9_distddm1[26:31,10] <- 2.346557

chr10_distddm1[1:5,10] <- 2.346557
chr10_distddm1[24:29,10] <- 2.346557

ddm1_wt <- function(chr_bin, ddm1_dist){
  for(i in 1:nrow(chr_bin)){
    for(k in 1:nrow(ddm1_dist)){
      if(isTRUE(chr_bin$foo.X1[i] >= ddm1_dist$Start[k] && chr_bin$foo.X2 <= ddm1_dist$End[k])){
        chr_bin$final[i] <- chr_bin$rate[i] + (chr_bin$rate[i]*ddm1_dist$diffwt[k])
      }
    }
  }
  return(chr_bin)
}
chr1_w_ddm1 <- ddm1_wt(chr1_bin, chr1_distddm1)
chr2_w_ddm1 <- ddm1_wt(chr2_bin, chr2_distddm1)
chr3_w_ddm1 <- ddm1_wt(chr3_bin, chr3_distddm1)
chr4_w_ddm1 <- ddm1_wt(chr4_bin, chr4_distddm1)
chr5_w_ddm1 <- ddm1_wt(chr5_bin, chr5_distddm1)
chr6_w_ddm1 <- ddm1_wt(chr6_bin, chr6_distddm1)
chr7_w_ddm1 <- ddm1_wt(chr7_bin, chr7_distddm1)
chr8_w_ddm1 <- ddm1_wt(chr8_bin, chr8_distddm1)
chr9_w_ddm1 <- ddm1_wt(chr9_bin, chr9_distddm1)
chr10_w_ddm1 <- ddm1_wt(chr10_bin, chr10_distddm1)

#using function, converted SNP start to Mb to get cM/Mb for final genetic position
chr1_snp2ddm1 <- snp_rate_ddm1(chr1_w_ddm1, chr1_snp)
chr1_snp2ddm1$`SNP Start`<- chr1_snp2ddm1$`SNP Start`/1000000
chr1_snp2ddm1 <- chr1_snp2ddm1[order(chr1_snp2ddm1$`SNP Start`),]
#creation of genetic positions from smoothed recombination rate
chr1_snp2ddm1$pos <- gen_pos(chr1_snp2ddm1)
chr1_snp2ddm1$pos2 <- graph_recomb(chr1_snp2ddm1)
#graph to look at Mb vs. cM along chromosome
plot(chr1_snp2$`SNP Start`, chr1_snp2$pos, type = "l", col = "blue", 
     main = "Chr1. Genetic Maps", xlab = "Physical Positions (Mb)", ylab = "Recombination rate (cM/Mb)")
chr1_finalpos <- chr1_snp2ddm1[order(chr1_snp2ddm1$pos),]
#want False to input into AlphaSimR
is.unsorted(chr1_finalpos$pos)

chr2_snp2ddm1 <- snp_rate_ddm1(chr2_w_ddm1, chr2_snp)
chr2_snp2ddm1$`SNP Start` <- chr2_snp2ddm1$`SNP Start`/1000000
chr2_snp2ddm1$pos <- gen_pos(chr2_snp2ddm1)
plot(chr2_snp2ddm1$`SNP Start`, chr2_snp2ddm1$pos)
chr2_finalpos <- chr2_snp2ddm1[order(chr2_snp2ddm1$pos),]
is.unsorted(chr2_finalpos$pos)

chr3_snp2ddm1 <- snp_rate_ddm1(chr3_w_ddm1, chr3_snp)
chr3_snp2ddm1$`SNP Start` <- chr3_snp2ddm1$`SNP Start`/1000000
chr3_snp2ddm1$pos <- gen_pos(chr3_snp2ddm1)
plot(chr3_snp2ddm1$`SNP Start`, chr3_snp2ddm1$pos)
chr3_finalpos <- chr3_snp2ddm1[order(chr3_snp2ddm1$pos),]
is.unsorted(chr3_finalpos$pos)

chr4_snp2ddm1 <- snp_rate_ddm1(chr4_w_ddm1, chr4_snp)
chr4_snp2ddm1$`SNP Start` <- chr4_snp2ddm1$`SNP Start`/1000000
chr4_snp2ddm1$pos <- gen_pos(chr4_snp2ddm1)
plot(chr4_snp2$`SNP Start`, chr4_snp2ddm1$pos)
chr4_finalpos <- chr4_snp2ddm1[order(chr4_snp2ddm1$pos),]
is.unsorted(chr4_finalpos$pos)

chr5_snp2ddm1 <- snp_rate_ddm1(chr5_w_ddm1, chr5_snp)
chr5_snp2ddm1$`SNP Start` <- chr5_snp2ddm1$`SNP Start`/1000000
chr5_snp2ddm1$pos <- gen_pos(chr5_snp2ddm1)
plot(chr5_snp2ddm1$`SNP Start`, chr5_snp2ddm1$pos)
chr5_finalpos <- chr5_snp2ddm1[order(chr5_snp2ddm1$pos),]
is.unsorted(chr5_finalpos$pos)

chr6_snp2ddm1 <- snp_rate_ddm1(chr6_w_ddm1, chr6_snp)
chr6_snp2ddm1$`SNP Start` <- chr6_snp2ddm1$`SNP Start`/1000000
chr6_snp2ddm1$pos <- gen_pos(chr6_snp2ddm1)
plot(chr6_snp2ddm1$`SNP Start`, chr6_snp2ddm1$pos)
chr6_finalpos <- chr6_snp2ddm1[order(chr6_snp2ddm1$pos),]
is.unsorted(chr6_finalpos$pos)

chr7_snp2ddm1 <- snp_rate_ddm1(chr7_w_ddm1, chr7_snp)
chr7_snp2ddm1$`SNP Start` <- chr7_snp2ddm1$`SNP Start`/1000000
chr7_snp2ddm1$pos <- gen_pos(chr7_snp2ddm1)
plot(chr7_snp2ddm1$`SNP Start`, chr7_snp2ddm1$pos)
chr7_finalpos <- chr7_snp2ddm1[order(chr7_snp2ddm1$pos),]
is.unsorted(chr7_finalpos$pos)

chr8_snp2ddm1 <- snp_rate_ddm1(chr8_w_ddm1, chr8_snp)
chr8_snp2ddm1$`SNP Start` <- chr8_snp2ddm1$`SNP Start`/1000000
chr8_snp2ddm1$pos <- gen_pos(chr8_snp2ddm1)
plot(chr8_snp2ddm1$`SNP Start`, chr8_snp2ddm1$pos)
chr8_finalpos <- chr8_snp2ddm1[order(chr8_snp2ddm1$pos),]
is.unsorted(chr8_finalpos$pos)

chr9_snp2ddm1 <- snp_rate_ddm1(chr9_w_ddm1, chr9_snp)
chr9_snp2ddm1$`SNP Start` <- chr9_snp2ddm1$`SNP Start`/1000000
chr9_snp2ddm1$pos <- gen_pos(chr9_snp2ddm1)
plot(chr9_snp2$`SNP Start`, chr9_snp2ddm1$pos)
chr9_finalpos <- chr9_snp2ddm1[order(chr9_snp2ddm1$pos),]
is.unsorted(chr9_finalpos$pos)

chr10_snp2ddm1 <- snp_rate_ddm1(chr10_w_ddm1, chr10_snp)
chr10_snp2ddm1$`SNP Start` <- chr10_snp2ddm1$`SNP Start`/1000000
chr10_snp2ddm1$pos <- gen_pos(chr10_snp2ddm1)
plot(chr10_snp2ddm1$`SNP Start`, chr10_snp2ddm1$pos)
chr10_finalpos <- chr10_snp2ddm1[order(chr10_snp2ddm1$pos),]
is.unsorted(chr10_finalpos$pos)

#Putting the final genetic map together
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
for(i in 1:10){
  names(ddm1_map[[i]]) = paste(i, 1:segSites[i], sep="_")
}

saveRDS(ddm1_map, file="ddm1_map.RData")

#Creating vector of centromere positions
#chnage this to reflect new pos
ddm1_centromere <- c(398.6166, 340.9768, 209.0365, 204.1758, 291,
                     63.18590, 256.7602, 244.40481, 215.48164, 219.13686)
ddm1_centromere <- (ddm1_centromere/100)
