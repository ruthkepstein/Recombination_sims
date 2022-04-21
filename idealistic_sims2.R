ideal1_dist <- read.table("maize_genome_ddm1_zmet2.txt", header = FALSE)
colnames(ideal1_dist) <- c("Chr", "Start", "End", "Female WT", "Male WT", "ddm1_1", "ddm1_2", "zmet2")

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

#10X increase through the genome
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
#smoothing the recombination rate so transitions between bins are not so abrupt
chr1_splideal12 <- smooth.spline(chr1_snp2ideal1$rate, spar = 1.2)
#creation of genetic positions from smoothed recombination rate
chr1_snp2ideal1$pos <- (chr1_snp2ideal1$`SNP Start`*chr1_splideal12$y)
#graph to look at Mb vs. cM along chromosome
plot(chr1_snp2ideal1$`SNP Start`, chr1_snp2ideal1$pos, type = "l", col = "blue", 
     main = "Chr1. Genetic Maps", xlab = "Physical Positions (Mb)", ylab = "Recombination rate (cM/Mb)")
ggplot(chr1_snp2ideal1, aes(`SNP Start`,pos)) + geom_point() + geom_smooth()
#graph to look at Mb vs. cM/Mb to see recombination rate along chromosome
plot(chr1_snp2ideal1$`SNP Start`, chr1_snp2ideal1$pos/chr1_snp2ideal1$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 1 Recombination Distribution")
chr1_finalpos <- chr1_snp2ideal1[order(chr1_snp2ideal1$pos),]
#want False to input into AlphaSimR
is.unsorted(chr1_finalpos$pos)
#plot again to make sure it looks the same
plot(chr1_snp2ideal1$`SNP Start`, chr1_finalpos$pos/chr1_snp2ideal1$`SNP Start`, type = "l", xlab = "Physical Positions (Mb)",
     ylab = "Recombination rate (cM/Mb)", main = "Chromosome 1 Recombination Distribution")
plot(chr1_finalpos$`SNP Start`, chr1_finalpos$pos)


chr2_snp2ideal1 <- snp_rate_ideal1(chr2_w_ideal1, chr2_snp)
chr2_snp2ideal1$`SNP Start` <- chr2_snp2ideal1$`SNP Start`/1000000
chr2_snp2ideal1 <- chr2_snp2ideal1[-(225:237),]
chr2_splideal1 <- smooth.spline(chr2_snp2ideal1$rate, spar = 1.3)
chr2_snp2ideal1$pos <- (chr2_snp2ideal1$`SNP Start`*chr2_splideal1$y)
plot(chr2_snp2ideal1$`SNP Start`, chr2_snp2ideal1$pos)
plot(chr2_snp2ideal1$`SNP Start`, chr2_snp2ideal1$pos/chr2_snp2ideal1$`SNP Start`, type = "l")
chr2_finalpos <- chr2_snp2ideal1[order(chr2_snp2ideal1$pos),]
is.unsorted(chr2_finalpos$pos)
plot(chr2_snp2ideal1$`SNP Start`, chr2_finalpos$pos/chr2_snp2ideal1$`SNP Start`, type = "l")

chr3_snp2ideal1 <- snp_rate_ideal1(chr3_w_ideal1, chr3_snp)
chr3_snp2ideal1$`SNP Start` <- chr3_snp2ideal1$`SNP Start`/1000000
chr3_splideal1 <- smooth.spline(chr3_snp2ideal1$rate, spar = 1.25)
chr3_snp2ideal1$pos <- (chr3_snp2ideal1$`SNP Start`*chr3_splideal1$y)
plot(chr3_snp2ideal1$`SNP Start`, chr3_snp2ideal1$pos)
plot(chr3_snp2ideal1$`SNP Start`, chr3_snp2ideal1$pos/chr3_snp2ideal1$`SNP Start`, type = "l")
chr3_finalpos <- chr3_snp2ideal1[order(chr3_snp2ideal1$pos),]
is.unsorted(chr3_finalpos$pos)
plot(chr3_snp2ideal1$`SNP Start`, chr3_finalpos$pos/chr3_snp2ideal1$`SNP Start`, type = "l")

chr4_snp2ideal1 <- snp_rate_ideal1(chr4_w_ideal1, chr4_snp)
chr4_snp2ideal1$`SNP Start` <- chr4_snp2ideal1$`SNP Start`/1000000
chr4_splideal1 <- smooth.spline(chr4_snp2ideal1$rate, spar = 1.25)
chr4_snp2ideal1$pos <- (chr4_snp2ideal1$`SNP Start`*chr4_splideal1$y)
plot(chr4_snp2ideal1$`SNP Start`, chr4_snp2ideal1$pos)
plot(chr4_snp2ideal1$`SNP Start`, chr4_snp2ideal1$pos/chr4_snp2ideal1$`SNP Start`, type = "l")
chr4_finalpos <- chr4_snp2ideal1[order(chr4_snp2ideal1$pos),]
is.unsorted(chr4_finalpos$pos)
plot(chr4_snp2ideal1$`SNP Start`, chr4_finalpos$pos/chr4_snp2ideal1$`SNP Start`, type = "l")

chr5_snp2ideal1 <- snp_rate_ideal1(chr5_w_ideal1, chr5_snp)
chr5_snp2ideal1$`SNP Start` <- chr5_snp2ideal1$`SNP Start`/1000000
chr5_splideal1 <- smooth.spline(chr5_snp2ideal1$rate, spar = 1.22)
chr5_snp2ideal1$pos <- (chr5_snp2ideal1$`SNP Start`*chr5_splideal1$y)
plot(chr5_snp2ideal1$`SNP Start`, chr5_snp2ideal1$pos)
plot(chr5_snp2ideal1$`SNP Start`, chr5_snp2ideal1$pos/chr5_snp2ideal1$`SNP Start`, type = "l")
chr5_finalpos <- chr5_snp2ideal1[order(chr5_snp2ideal1$pos),]
is.unsorted(chr5_finalpos$pos)
plot(chr5_snp2ideal1$`SNP Start`, chr5_finalpos$pos/chr5_snp2ideal1$`SNP Start`, type = "l")

#chr 6 is lowkey fuked up
chr6_snp2ideal1 <- snp_rate_ideal1(chr6_w_ideal1, chr6_snp)
chr6_snp2ideal1$`SNP Start` <- chr6_snp2ideal1$`SNP Start`/1000000
chr6_splideal1 <- smooth.spline(chr6_snp2ideal1$rate, spar = 1)
chr6_snp2ideal1$pos <- (chr6_snp2ideal1$`SNP Start`*chr6_splideal1$y)
plot(chr6_snp2ideal1$`SNP Start`, chr6_snp2ideal1$pos)
plot(chr6_snp2ideal1$`SNP Start`, chr6_snp2ideal1$pos/chr6_snp2ideal1$`SNP Start`, type = "l")
chr6_finalpos <- chr6_snp2ideal1[order(chr6_snp2ideal1$pos),]
is.unsorted(chr6_finalpos$pos)
plot(chr6_snp2ideal1$`SNP Start`, chr6_finalpos$pos/chr6_snp2ideal1$`SNP Start`, type = "l")

chr7_snp2ideal1 <- snp_rate_ideal1(chr7_w_ideal1, chr7_snp)
chr7_snp2ideal1$`SNP Start` <- chr7_snp2ideal1$`SNP Start`/1000000
chr7_splideal1 <- smooth.spline(chr7_snp2ideal1$rate, spar = 1.2)
chr7_snp2ideal1$pos <- (chr7_snp2ideal1$`SNP Start`*chr7_splideal1$y)
plot(chr7_snp2ideal1$`SNP Start`, chr7_snp2ideal1$pos)
plot(chr7_snp2ideal1$`SNP Start`, chr7_snp2ideal1$pos/chr7_snp2ideal1$`SNP Start`, type = "l")
chr7_finalpos <- chr7_snp2ideal1[order(chr7_snp2ideal1$pos),]
is.unsorted(chr7_finalpos$pos)
plot(chr7_snp2ideal1$`SNP Start`, chr7_finalpos$pos/chr7_snp2ideal1$`SNP Start`, type = "l")

chr8_snp2ideal1 <- snp_rate_ideal1(chr8_w_ideal1, chr8_snp)
chr8_snp2ideal1$`SNP Start` <- chr8_snp2ideal1$`SNP Start`/1000000
chr8_splideal1 <- smooth.spline(chr8_snp2ideal1$rate, spar = 1.2)
chr8_snp2ideal1$pos <- (chr8_snp2ideal1$`SNP Start`*chr8_splideal1$y)
plot(chr8_snp2ideal1$`SNP Start`, chr8_snp2ideal1$pos)
plot(chr8_snp2ideal1$`SNP Start`, chr8_snp2ideal1$pos/chr8_snp2ideal1$`SNP Start`, type = "l")
chr8_finalpos <- chr8_snp2ideal1[order(chr8_snp2ideal1$pos),]
is.unsorted(chr8_finalpos$pos)
plot(chr8_snp2ideal1$`SNP Start`, chr8_finalpos$pos/chr8_snp2ideal1$`SNP Start`, type = "l")

chr9_snp2ideal1 <- snp_rate_ideal1(chr9_w_ideal1, chr9_snp)
chr9_snp2ideal1$`SNP Start` <- chr9_snp2ideal1$`SNP Start`/1000000
chr9_splideal1 <- smooth.spline(chr9_snp2ideal1$rate, spar = 1.15)
chr9_snp2ideal1$pos <- (chr9_snp2ideal1$`SNP Start`*chr9_splideal1$y)
plot(chr9_snp2ideal1$`SNP Start`, chr9_snp2ideal1$pos)
plot(chr9_snp2ideal1$`SNP Start`, chr9_snp2ideal1$pos/chr9_snp2ideal1$`SNP Start`, type = "l")
chr9_finalpos <- chr9_snp2ideal1[order(chr9_snp2ideal1$pos),]
is.unsorted(chr9_finalpos$pos)
plot(chr9_snp2ideal1$`SNP Start`, chr9_finalpos$pos/chr9_snp2ideal1$`SNP Start`, type = "l")

chr10_snp2ideal1 <- snp_rate_ideal1(chr10_w_ideal1, chr10_snp)
chr10_snp2ideal1$`SNP Start` <- chr10_snp2ideal1$`SNP Start`/1000000
chr10_splideal1 <- smooth.spline(chr10_snp2ideal1$rate, spar = 1.2)
chr10_snp2ideal1$pos <- (chr10_snp2ideal1$`SNP Start`*chr10_splideal1$y)
plot(chr10_snp2ideal1$`SNP Start`, chr10_snp2ideal1$pos)
plot(chr10_snp2ideal1$`SNP Start`, chr10_snp2ideal1$pos/chr10_snp2ideal1$`SNP Start`, type = "l")
chr10_finalpos <- chr10_snp2ideal1[order(chr10_snp2ideal1$pos),]
is.unsorted(chr10_finalpos$pos)
plot(chr10_snp2ideal1$`SNP Start`, chr10_finalpos$pos/chr10_snp2ideal1$`SNP Start`, type = "l")

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

ideal1_map <- list(chr1[[1]], chr2[[1]], 
                  chr3[[1]], chr4[[1]], chr5[[1]], 
                  chr6[[1]], chr7[[1]], chr8[[1]], 
                  chr9[[1]], chr10[[1]])

#saveRDS(ideal1_map, file="ideal1_map.RData")
#Creating vector of centromere positions
#change this
ideal1_centromere <- c(69.19425, 29.31842, 43.34131, 42.1754, 47.11507,
                     43.40838, 30.01601, 43.69602, 47.12911, 21.17685)
ideal1_centromere <- (ideal1_centromere*10)/100
