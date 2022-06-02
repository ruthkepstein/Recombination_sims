snp_rate_targeted <- function(chr_w_targeted, chr_snp){
  for(i in 1:nrow(chr_snp)){
    for(k in 1:nrow(chr_w_targeted)){
      if(isTRUE((chr_snp$`SNP Start`[i] >= chr_w_targeted$foo.X1[k]) && (chr_snp$`SNP Start`[i] <= chr_w_targeted$foo.X2[k]))){
        chr_snp$rate[i] <- chr_w_targeted$final[k]
      }
    }
  }
  print(chr_snp)
}

##Using targeted recombination landscape--> 6% increase in COs
targeted_dist <- read.table("maize_genome_ddm1_zmet2.txt", header = FALSE)
colnames(targeted_dist) <- c("Chr", "Start", "End", "Female WT", "Male WT", "ddm1_1", "ddm1_2", "zmet2")
#targeted_dist$`Female WT` <- targeted_dist$`Female WT`*2/122
targeted_dist$`Male WT`<- targeted_dist$`Male WT`*2/135

targeted_dist$diffwt <- 0
#make generalization about targeted

chr1_disttargeted <- targeted_dist[ which(targeted_dist$Chr == 1),]
chr2_disttargeted <- targeted_dist[ which(targeted_dist$Chr == 2),]
chr3_disttargeted <- targeted_dist[ which(targeted_dist$Chr == 3),]
chr4_disttargeted <- targeted_dist[ which(targeted_dist$Chr == 4),]
chr5_disttargeted <- targeted_dist[ which(targeted_dist$Chr == 5),]
chr6_disttargeted <- targeted_dist[ which(targeted_dist$Chr == 6),]
chr7_disttargeted <- targeted_dist[ which(targeted_dist$Chr == 7),]
chr8_disttargeted <- targeted_dist[ which(targeted_dist$Chr == 8),]
chr9_disttargeted <- targeted_dist[ which(targeted_dist$Chr == 9),]
chr10_disttargeted <- targeted_dist[ which(targeted_dist$Chr == 10),]

chr1_disttargeted[55:60,9] <- 20
#chr1_disttargeted[50:60,10] <- 2.346557

targeted_wt <- function(chr_bin, targeted_dist){
  for(i in 1:nrow(chr_bin)){
    for(k in 1:nrow(targeted_dist)){
      if(isTRUE(chr_bin$foo.X1[i] >= targeted_dist$Start[k] && chr_bin$foo.X2 <= targeted_dist$End[k])){
        chr_bin$final[i] <- chr_bin$rate[i] + (chr_bin$rate[i]*targeted_dist$diffwt[k])
      }
    }
  }
  return(chr_bin)
}
chr1_w_targeted <- targeted_wt(chr1_bin, chr1_disttargeted)
chr2_w_targeted <- targeted_wt(chr2_bin, chr2_disttargeted)
chr3_w_targeted <- targeted_wt(chr3_bin, chr3_disttargeted)
chr4_w_targeted <- targeted_wt(chr4_bin, chr4_disttargeted)
chr5_w_targeted <- targeted_wt(chr5_bin, chr5_disttargeted)
chr6_w_targeted <- targeted_wt(chr6_bin, chr6_disttargeted)
chr7_w_targeted <- targeted_wt(chr7_bin, chr7_disttargeted)
chr8_w_targeted <- targeted_wt(chr8_bin, chr8_disttargeted)
chr9_w_targeted <- targeted_wt(chr9_bin, chr9_disttargeted)
chr10_w_targeted <- targeted_wt(chr10_bin, chr10_disttargeted)

#using function, converted SNP start to Mb to get cM/Mb for final genetic position
chr1_snp2 <- snp_rate_targeted(chr1_w_targeted, chr1_snp)
chr1_snp2$`SNP Start`<- chr1_snp2$`SNP Start`/1000000
chr1_snp2 <- chr1_snp2[order(chr1_snp2$`SNP Start`),]
#creation of genetic positions from smoothed recombination rate
chr1_snp2$pos <- gen_pos(chr1_snp2)
chr1_snp2$pos2 <- graph_recomb(chr1_snp2)
#graph to look at Mb vs. cM along chromosome
plot(chr1_snp2$`SNP Start`, chr1_snp2$pos2/chr1_snp2$`SNP Start`, type = "l", col = "blue", 
     main = "Chr1. Genetic Maps", xlab = "Physical Positions (Mb)", ylab = "Recombination rate (cM/Mb)")
chr1_finalpos <- chr1_snp2[order(chr1_snp2$pos),]
#want False to input into AlphaSimR
is.unsorted(chr1_finalpos$pos)

chr2_snp2 <- snp_rate_targeted(chr2_w_targeted, chr2_snp)
chr2_snp2$`SNP Start` <- chr2_snp2$`SNP Start`/1000000
chr2_snp2$pos <- gen_pos(chr2_snp2)
plot(chr2_snp2$`SNP Start`, chr2_snp2$pos)
chr2_finalpos <- chr2_snp2[order(chr2_snp2$pos),]
is.unsorted(chr2_finalpos$pos)

chr3_snp2 <- snp_rate_targeted(chr3_w_targeted, chr3_snp)
chr3_snp2$`SNP Start` <- chr3_snp2$`SNP Start`/1000000
chr3_snp2$pos <- gen_pos(chr3_snp2)
plot(chr3_snp2$`SNP Start`, chr3_snp2$pos)
chr3_finalpos <- chr3_snp2[order(chr3_snp2$pos),]
is.unsorted(chr3_finalpos$pos)

chr4_snp2 <- snp_rate_targeted(chr4_w_targeted, chr4_snp)
chr4_snp2$`SNP Start` <- chr4_snp2$`SNP Start`/1000000
chr4_snp2$pos <- gen_pos(chr4_snp2)
plot(chr4_snp2$`SNP Start`, chr4_snp2$pos)
chr4_finalpos <- chr4_snp2[order(chr4_snp2$pos),]
is.unsorted(chr4_finalpos$pos)

chr5_snp2 <- snp_rate_targeted(chr5_w_targeted, chr5_snp)
chr5_snp2$`SNP Start` <- chr5_snp2$`SNP Start`/1000000
chr5_snp2$pos <- gen_pos(chr5_snp2)
plot(chr5_snp2$`SNP Start`, chr5_snp2$pos)
chr5_finalpos <- chr5_snp2[order(chr5_snp2$pos),]
is.unsorted(chr5_finalpos$pos)

chr6_snp2 <- snp_rate_targeted(chr6_w_targeted, chr6_snp)
chr6_snp2$`SNP Start` <- chr6_snp2$`SNP Start`/1000000
chr6_snp2$pos <- gen_pos(chr6_snp2)
plot(chr6_snp2$`SNP Start`, chr6_snp2$pos)
chr6_finalpos <- chr6_snp2[order(chr6_snp2$pos),]
is.unsorted(chr6_finalpos$pos)

chr7_snp2 <- snp_rate_targeted(chr7_w_targeted, chr7_snp)
chr7_snp2$`SNP Start` <- chr7_snp2$`SNP Start`/1000000
chr7_snp2$pos <- gen_pos(chr7_snp2)
plot(chr7_snp2$`SNP Start`, chr7_snp2$pos)
chr7_finalpos <- chr7_snp2[order(chr7_snp2$pos),]
is.unsorted(chr7_finalpos$pos)

chr8_snp2 <- snp_rate_targeted(chr8_w_targeted, chr8_snp)
chr8_snp2$`SNP Start` <- chr8_snp2$`SNP Start`/1000000
chr8_snp2$pos <- gen_pos(chr8_snp2)
plot(chr8_snp2$`SNP Start`, chr8_snp2$pos)
chr8_finalpos <- chr8_snp2[order(chr8_snp2$pos),]
is.unsorted(chr8_finalpos$pos)

chr9_snp2 <- snp_rate_targeted(chr9_w_targeted, chr9_snp)
chr9_snp2$`SNP Start` <- chr9_snp2$`SNP Start`/1000000
chr9_snp2$pos <- gen_pos(chr9_snp2)
plot(chr9_snp2$`SNP Start`, chr9_snp2$pos)
chr9_finalpos <- chr9_snp2[order(chr9_snp2$pos),]
is.unsorted(chr9_finalpos$pos)

chr10_snp2 <- snp_rate_targeted(chr10_w_targeted, chr10_snp)
chr10_snp2$`SNP Start` <- chr10_snp2$`SNP Start`/1000000
chr10_snp2$pos <- gen_pos(chr10_snp2)
plot(chr10_snp2$`SNP Start`, chr10_snp2$pos)
chr10_finalpos <- chr10_snp2[order(chr10_snp2$pos),]
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

targeted2_map = vector("list",10)
targeted2_map[[1]] = chr1
targeted2_map[[2]] = chr2
targeted2_map[[3]] = chr3
targeted2_map[[4]] = chr4
targeted2_map[[5]] = chr5
targeted2_map[[6]] = chr6
targeted2_map[[7]] = chr7
targeted2_map[[8]] = chr8
targeted2_map[[9]] = chr9
targeted2_map[[10]] = chr10
for(i in 1:10){
  names(targeted2_map[[i]]) = paste(i, 1:segSites[i], sep="_")
}

#adding names to each individual SNP based on loci and chromosome #

saveRDS(targeted2_map, file = "targeted2_map.RData")

targeted_centromere <- real_centromere/100

library(AlphaSimR)

#generating mix of repulsion & coupling linkages between QTL
addEff_mix <- rnorm(300, mean = 0, sd = 0.1)
sum(addEff_mix)

#Getting QTLs from accurate gene space
chr1_QTL <- readRDS("chr1_QTLpos.Rdata")
chr2_QTL <- readRDS("chr2_QTLpos.Rdata")
chr3_QTL <- readRDS("chr3_QTLpos.Rdata")
chr4_QTL <- readRDS("chr4_QTLpos.Rdata")
chr5_QTL <- readRDS("chr5_QTLpos.Rdata")
chr6_QTL <- readRDS("chr6_QTLpos.Rdata")
chr7_QTL <- readRDS("chr7_QTLpos.Rdata")
chr8_QTL <- readRDS("chr8_QTLpos.Rdata")
chr9_QTL <- readRDS("chr9_QTLpos.Rdata")
chr10_QTL <- readRDS("chr10_QTLpos.Rdata")

#Setting up founder population
founderPop <- quickHaplo(nInd = 200, nChr = 10, inbred = TRUE, ploidy = 2L, segSites = segSites)
founderPop@genMap <- final_map
founderPop@centromere <- real_centromere
SP = SimParam$new(founderPop)
SP$setTrackRec(TRUE)
SP$p = 0.15
trait_yield <- new("TraitA", nLoci = 300L, lociPerChr= c(40L, 35L, 35L, 35L, 35L, 30L, 25L, 25L, 20L, 20L),
                   lociLoc = as.integer(c(chr1_QTL, chr2_QTL, chr3_QTL, chr4_QTL, chr5_QTL,
                                          chr6_QTL, chr7_QTL, chr8_QTL, chr9_QTL, chr10_QTL)), addEff = addEff_mix, intercept = 0.1)
SP$manAddTrait(trait_yield)

#Diverging founder population into 2 populations; "good pop" & "bad pop"
#Good pop = selecting FOR polygenic trait; no "resistance" QTLs
pop_good_sel10 <- vector(mode = "list", length = 100)
for(i in 1:100){
  pop_good <- newPop(founderPop, simParam = SP)
  pop_good <- setPheno(pop_good, h2 = 0.8, simParam = SP)
  
  pop_good1 <- randCross(pop_good, nCrosses = 10, nProgeny=10, simParam = SP)
  pop_good1 <- setPheno(pop_good1, h2 = 0.8, simParam = SP)
  
  pop_good_sel <- selectInd(pop_good1, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good2 <- randCross(pop_good_sel, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good2 <- setPheno(pop_good2, h2 = 0.8, simParam = SP)
  
  pop_good_sel2 <- selectInd(pop_good2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good3 <- randCross(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good3 <- setPheno(pop_good3, h2 = 0.8, simParam = SP)
  
  pop_good_sel3 <- selectInd(pop_good3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good4 <- randCross(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good4 <- setPheno(pop_good4, h2 = 0.8, simParam = SP)
  
  pop_good_sel4 <- selectInd(pop_good4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good5 <- randCross(pop_good_sel4, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good5 <- setPheno(pop_good5, h2 = 0.8, simParam = SP)
  
  pop_good_sel5 <- selectInd(pop_good5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good6 <- randCross(pop_good_sel5, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good6 <- setPheno(pop_good6, h2 = 0.8, simParam = SP)
  
  pop_good_sel6 <- selectInd(pop_good6, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good7 <- randCross(pop_good_sel6, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good7 <- setPheno(pop_good7, h2 = 0.8, simParam = SP)
  
  pop_good_sel7 <- selectInd(pop_good7, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good8 <- randCross(pop_good_sel7, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good8 <- setPheno(pop_good8, h2 = 0.8, simParam = SP)
  
  pop_good_sel8 <- selectInd(pop_good8, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good9 <- randCross(pop_good_sel8, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good9 <- setPheno(pop_good9, h2 = 0.8, simParam = SP)
  
  pop_good_sel9 <- selectInd(pop_good9, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good10 <- randCross(pop_good_sel9, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good10 <- setPheno(pop_good10, h2 = 0.8, simParam = SP)
  
  pop_good_sel10[[i]] <- selectInd(pop_good10, nInd = 5, use = "gv", trait = 1, selectop = TRUE, returnPop = TRUE, simParam = SP)
}

#Bad pop = not selecting for polygenic trait even though its present; "resistance" QTLs present & selected FOR
trait <- new("TraitAD", nLoci = 3L, lociPerChr = c(3L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
             lociLoc = c(260L, 272L, 288L), addEff = c(0.2, 0.2, 0.2), domEff = c(0.2, 0.2, 0.2), intercept = 0.1)
SP$resetPed()
SP$manAddTrait(trait)
pop_bad_sel10 <- vector(mode = "list", length = 100)
for(i in 1:100){
  pop_bad <- newPop(founderPop, simParam = SP)
  pop_bad <- setPheno(pop_bad, h2 = 0.8, simParam = SP)
  
  pop_bad1 <- randCross(pop_bad, nCrosses = 10, nProgeny=10, simParam = SP)
  pop_bad1 <- setPheno(pop_bad1, h2 = 0.8, simParam = SP)
  
  pop_bad2 <- selectInd(pop_bad1, nInd = 10, use = "gv", trait = 1, selectTop = FALSE, simParam = SP)
  pop_bad2 <- selectCross(pop_bad2, nInd = 5, use = "gv", trait = 2, selectTop = TRUE, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_bad2 <- setPheno(pop_bad2, h2 = 0.8, simParam = SP)
  
  pop_bad3 <- selectInd(pop_bad2, nInd = 10, use = "gv", trait = 1, selectTop = FALSE, simParam = SP)
  pop_bad3 <- selectCross(pop_bad3, nInd = 5, use = "gv", trait = 2, selectTop = TRUE, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_bad3 <- setPheno(pop_bad3, h2 = 0.8, simParam = SP)
  
  pop_bad4 <- selectInd(pop_bad3, nInd = 10, use = "gv", trait = 1, selectTop = FALSE, simParam = SP)
  pop_bad4 <- selectCross(pop_bad4, nInd = 5, use = "gv", trait = 2, selectTop = TRUE, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_bad4 <- setPheno(pop_bad4, h2 = 0.8, simParam = SP)
  
  pop_bad5 <- selectInd(pop_bad4, nInd = 10, use = "gv", trait = 1, selectTop = FALSE, simParam = SP)
  pop_bad5 <- selectCross(pop_bad5, nInd = 5, use = "gv", trait = 2, selectTop = TRUE, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_bad5 <- setPheno(pop_bad5, h2 = 0.8, simParam = SP)
  
  pop_bad6 <- selectInd(pop_bad5, nInd = 10, use = "gv", trait = 1, selectTop = FALSE, simParam = SP)
  pop_bad6 <- selectCross(pop_bad6, nInd = 5, use = "gv", trait = 2, selectTop = TRUE, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_bad6 <- setPheno(pop_bad6, h2 = 0.8, simParam = SP)
  
  pop_bad7 <- selectInd(pop_bad6, nInd = 10, use = "gv", trait = 1, selectTop = FALSE, simParam = SP)
  pop_bad7 <- selectCross(pop_bad7, nInd = 5, use = "gv", trait = 2, selectTop = TRUE, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_bad7 <- setPheno(pop_bad7, h2 = 0.8, simParam = SP)
  
  pop_bad8 <- selectInd(pop_bad7, nInd = 10, use = "gv", trait = 1, selectTop = FALSE, simParam = SP)
  pop_bad8 <- selectCross(pop_bad8, nInd = 5, use = "gv", trait = 2, selectTop = TRUE, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_bad8 <- setPheno(pop_bad8, h2 = 0.8, simParam = SP)
  
  pop_bad9 <- selectInd(pop_bad8, nInd = 10, use = "gv", trait = 1, selectTop = FALSE, simParam = SP)
  pop_bad9 <- selectCross(pop_bad9, nInd = 5, use = "gv", trait = 2, selectTop = TRUE, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_bad9 <- setPheno(pop_bad9, h2 = 0.8, simParam = SP)
  
  pop_bad10 <- selectInd(pop_bad9, nInd = 10, use = "gv", trait = 1, selectTop = FALSE, simParam = SP)
  pop_bad10 <- selectCross(pop_bad10, nInd = 5, use = "gv", trait = 2, selectTop = TRUE, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_bad10 <- setPheno(pop_bad10, h2 = 0.8, simParam = SP)
  
  pop_bad_sel10[[i]] <- selectInd(pop_bad10, nInd = 5, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
}

goodpop = mergePops(pop_good_sel10)
badpop = mergePops(pop_bad_sel10)

##Backcrossing scheme
#selecting best individual from "good pop" and individual with resistance QTLs in "bad pop"
best_good <- which(goodpop@gv == max(goodpop@gv), arr.ind = FALSE)
best_good <- 403
#this changes every time depending on which row has highest GVs
best_bad <- 463

good_trait1_geno <- pullQtlGeno(goodpop[best_good,], trait = 1)
#elite parent has 0, 0, 0 at all resistance loci bc absent
good_trait2_geno <- pullQtlGeno(goodpop[best_good,], trait = 2)
bad_trait1_geno <- pullQtlGeno(badpop[best_bad,], trait = 1)
#wild parent has 2, 2, 2 at all resistance loci 
bad_trait2_geno <- pullQtlGeno(badpop[best_bad,], trait = 2)

#function for elucidating % of elite parent of progeny from genotype data
percent_recurrent <- function(elite, progeny, recurrent_count){
  for(i in 1:(length(progeny)/2)){
    if(isTRUE((progeny[i,1] == elite[i,1]) && (progeny[i,2] == elite[i,2]))){
      recurrent_count = recurrent_count + 1
    } else if(isTRUE((progeny[i,1] == elite[i,1]) && (progeny[i,2] != elite[i,2]))){
      recurrent_count = recurrent_count + .5
    } else if(isTRUE((progeny[i,1] != elite[i,1]) && (progeny[i,2] == elite[i,2]))){
      recurrent_count = recurrent_count + .5
    } else if(isTRUE((progeny[i,1] == elite[i,2]) || (progeny[i,2] == elite[i,1]))){
      recurrent_count = recurrent_count + .5
    }
  }
  return(recurrent_count)
}

#iterating over all chromosomes & individuals
recurrent_all <- function(elite, progeny){
  recurrent_count = 0
  for(j in 1:10){
    for(k in 1:5){
      recurrent_count = percent_recurrent(elite[[j]][,,1], progeny[[j]][,,k], recurrent_count)
    }
  }
  return(recurrent_count)
}

#function to look at loci around resistance loci to find if elite or wild parents are present
#to assess linkage drag
drag <- c()
range_around <- function(progeny_loci, good, bad){
  bad_ratio = 0
  good_ratio = 0
  for(k in 1:5){
    for(i in 1:40){
      if(good[i] != bad[i]){
        if(progeny_loci[k,i] == bad[i]){
          bad_ratio = bad_ratio + 1
        }
        if(progeny_loci[k,i] == good[i]){
          good_ratio = good_ratio + 1
        }
      }
      if(progeny_loci[k,i] == 1){
        good_ratio = good_ratio + .5
        bad_ratio = bad_ratio + .5
      }
    }
    drag[k] <- bad_ratio/(good_ratio + bad_ratio)
  }
  return(mean(drag))
}

#wild-type first
wtvar1 <- matrix(data = NA, ncol = 2, nrow = 100)
wtvar2 <- matrix(data = NA, ncol = 2, nrow = 100)
wtvar3 <- matrix(data = NA, ncol = 2, nrow = 100)
wtvar4 <- matrix(data = NA, ncol = 2, nrow = 100)
wtvar5 <- matrix(data = NA, ncol = 2, nrow = 100)
wtvar6 <- matrix(data = NA, ncol = 2, nrow = 100)
F1_elite_perc <- c()
F1_wild_perc <- c()
S1_elite_perc <- c()
S1_wild_perc <- c()
S2_elite_perc <- c()
S2_wild_perc <- c()
S3_elite_perc <- c()
S3_wild_perc <- c()
S4_elite_perc <- c()
S4_wild_perc <- c()
S5_elite_perc <- c()
S5_wild_perc <- c()
F1drag <- c()
S1drag <- c()
S2drag <- c()
S3drag <- c()
S4drag <- c()
S5drag <- c()
for(i in 1:100){
  SP$switchGenMap(final_map, centromere = real_centromere)
  pop <- randCross2(goodpop[best_good,], badpop[best_bad,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop <- setPheno(pop = pop, h2 = 0.8, simParam = SP)
  wtvar1[i,] <- genicVarA(pop)
  
  pop1_sel <- selectInd(pop, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_2 <- selectInd(pop1_sel, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  F1_elite_perc[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel1_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel1_2@geno)) + recurrent_all(badpop[best_bad,]@geno, pop1_sel1_2@geno))
  F1_wild_perc[i] <- recurrent_all(badpop[best_bad,]@geno, pop1_sel1_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel1_2@geno)) + recurrent_all(badpop[best_bad,]@geno, pop1_sel1_2@geno))
  F1drag[i] <- range_around(pullQtlGeno(pop1_sel1_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel1_cross <- randCross2(pop1_sel1_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel1_cross <- setPheno(pop1_sel1_cross, h2 = 0.8, simParam = SP)
  wtvar2[i,] <- genicVarA(pop1_sel1_2)
  
  pop1_sel2 <- selectInd(pop1_sel1_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_2 <- selectInd(pop1_sel2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S1_elite_perc[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel2_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel2_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel2_2@geno))
  S1_wild_perc[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel2_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel2_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel2_2@geno))
  S1drag[i] <- range_around(pullQtlGeno(pop1_sel2_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel2_cross <- randCross2(pop1_sel2_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel2_cross <- setPheno(pop1_sel2_cross, h2 = 0.8, simParam = SP)
  wtvar3[i,] <- genicVarA(pop1_sel2_2)
  
  pop1_sel3 <- selectInd(pop1_sel2_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_2 <- selectInd(pop1_sel3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S2_elite_perc[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel3_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel3_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel3_2@geno))
  S2_wild_perc[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel3_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel3_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel3_2@geno))
  S2drag[i] <- range_around(pullQtlGeno(pop1_sel3_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel3_cross <- randCross2(pop1_sel3_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel3_cross <- setPheno(pop1_sel3_cross, h2 = 0.8, simParam = SP)
  wtvar4[i,] <- genicVarA(pop1_sel3_2)
  
  pop1_sel4 <- selectInd(pop1_sel3_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_2 <- selectInd(pop1_sel4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S3_elite_perc[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel4_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_2@geno))
  S3_wild_perc[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_2@geno))
  S3drag[i] <- range_around(pullQtlGeno(pop1_sel4_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel4_cross <- randCross2(pop1_sel4_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel4_cross <- setPheno(pop1_sel4_cross, h2 = 0.8, simParam = SP)
  wtvar5[i,] <- genicVarA(pop1_sel4_2)
  
  pop1_sel5 <- selectInd(pop1_sel4_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_2 <- selectInd(pop1_sel5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S4_elite_perc[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel5_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel5_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel5_2@geno))
  S4_wild_perc[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel5_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel5_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel5_2@geno))
  S4drag[i] <- range_around(pullQtlGeno(pop1_sel5_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel5_cross <- randCross2(pop1_sel5_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel5_cross <- setPheno(pop1_sel5_cross, h2 = 0.8, simParam = SP)
  wtvar6[i,] <- genicVarA(pop1_sel5_2)
  
  pop1_sel6 <- selectInd(pop1_sel5_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel6_2 <- selectInd(pop1_sel6, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S5_elite_perc[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel6_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel6_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel6_2@geno))
  S5_wild_perc[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel6_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel6_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel6_2@geno))
  S5drag[i] <- range_around(pullQtlGeno(pop1_sel6_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
}

#targeted with 6-times
targetedvar1 <- matrix(data = NA, ncol = 2, nrow = 100)
targetedvar2 <- matrix(data = NA, ncol = 2, nrow = 100)
targetedvar3 <- matrix(data = NA, ncol = 2, nrow = 100)
targetedvar4 <- matrix(data = NA, ncol = 2, nrow = 100)
targetedvar5 <- matrix(data = NA, ncol = 2, nrow = 100)
targetedvar6 <- matrix(data = NA, ncol = 2, nrow = 100)
targetedvar7 <- matrix(data = NA, ncol = 2, nrow = 100)
F1_elite_perctargeted <- c()
F1_wild_perctargeted <- c()
S1_elite_perctargeted <- c()
S1_wild_perctargeted <- c()
S2_elite_perctargeted <- c()
S2_wild_perctargeted <- c()
S3_elite_perctargeted <- c()
S3_wild_perctargeted <- c()
S4_elite_perctargeted <- c()
S4_wild_perctargeted <- c()
S5_elite_perctargeted <- c()
S5_wild_perctargeted <- c()
F1dragtargeted <- c()
S1dragtargeted <- c()
S2dragtargeted <- c()
S3dragtargeted <- c()
S4dragtargeted <- c()
S5dragtargeted <- c()
for(i in 1:100){
  SP$switchGenMap(targeted_map, centromere = targeted_centromere)
  pop <- randCross2(goodpop[best_good,], badpop[best_bad,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop <- setPheno(pop = pop, h2 = 0.8, simParam = SP)
  targetedvar1[i,] <- genicVarA(pop)
  
  pop1_sel <- selectInd(pop, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_2 <- selectInd(pop1_sel, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  F1_elite_perctargeted[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel1_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel1_2@geno)) + recurrent_all(badpop[best_bad,]@geno, pop1_sel1_2@geno))
  F1_wild_perctargeted[i] <- recurrent_all(badpop[best_bad,]@geno, pop1_sel1_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel1_2@geno)) + recurrent_all(badpop[best_bad,]@geno, pop1_sel1_2@geno))
  F1dragtargeted[i] <- range_around(pullQtlGeno(pop1_sel1_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel1_cross <- randCross2(pop1_sel1_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel1_cross <- setPheno(pop1_sel1_cross, h2 = 0.8, simParam = SP)
  targetedvar2[i,] <- genicVarA(pop1_sel1_2)
  
  pop1_sel2 <- selectInd(pop1_sel1_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_2 <- selectInd(pop1_sel2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S1_elite_perctargeted[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel2_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel2_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel2_2@geno))
  S1_wild_perctargeted[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel2_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel2_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel2_2@geno))
  S1dragtargeted[i] <- range_around(pullQtlGeno(pop1_sel2_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel2_cross <- randCross2(pop1_sel2_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel2_cross <- setPheno(pop1_sel2_cross, h2 = 0.8, simParam = SP)
  targetedvar3[i,] <- genicVarA(pop1_sel2_2)
  
  pop1_sel3 <- selectInd(pop1_sel2_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_2 <- selectInd(pop1_sel3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S2_elite_perctargeted[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel3_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel3_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel3_2@geno))
  S2_wild_perctargeted[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel3_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel3_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel3_2@geno))
  S2dragtargeted[i] <- range_around(pullQtlGeno(pop1_sel3_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel3_cross <- randCross2(pop1_sel3_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel3_cross <- setPheno(pop1_sel3_cross, h2 = 0.8, simParam = SP)
  targetedvar4[i,] <- genicVarA(pop1_sel3_2)
  
  pop1_sel4 <- selectInd(pop1_sel3_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_2 <- selectInd(pop1_sel4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S3_elite_perctargeted[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel4_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_2@geno))
  S3_wild_perctargeted[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_2@geno))
  S3dragtargeted[i] <- range_around(pullQtlGeno(pop1_sel4_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel4_cross <- randCross2(pop1_sel4_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel4_cross <- setPheno(pop1_sel4_cross, h2 = 0.8, simParam = SP)
  targetedvar5[i,] <- genicVarA(pop1_sel4_2)
  
  pop1_sel5 <- selectInd(pop1_sel4_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_2 <- selectInd(pop1_sel5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S4_elite_perctargeted[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel5_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel5_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel5_2@geno))
  S4_wild_perctargeted[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel5_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel5_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel5_2@geno))
  S4dragtargeted[i] <- range_around(pullQtlGeno(pop1_sel5_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel5_cross <- randCross2(pop1_sel5_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel5_cross <- setPheno(pop1_sel5_cross, h2 = 0.8, simParam = SP)
  targetedvar6[i,] <- genicVarA(pop1_sel5_2)
  
  pop1_sel6 <- selectInd(pop1_sel5_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel6_2 <- selectInd(pop1_sel6, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S5_elite_perctargeted[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel6_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel6_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel6_2@geno))
  S5_wild_perctargeted[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel6_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel6_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel6_2@geno))
  S5dragtargeted[i] <- range_around(pullQtlGeno(pop1_sel6_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
}

#targeted with 20-times
targeted2var1 <- matrix(data = NA, ncol = 2, nrow = 100)
targeted2var2 <- matrix(data = NA, ncol = 2, nrow = 100)
targeted2var3 <- matrix(data = NA, ncol = 2, nrow = 100)
targeted2var4 <- matrix(data = NA, ncol = 2, nrow = 100)
targeted2var5 <- matrix(data = NA, ncol = 2, nrow = 100)
targeted2var6 <- matrix(data = NA, ncol = 2, nrow = 100)
targeted2var7 <- matrix(data = NA, ncol = 2, nrow = 100)
F1_elite_perctargeted2 <- c()
F1_wild_perctargeted2 <- c()
S1_elite_perctargeted2 <- c()
S1_wild_perctargeted2 <- c()
S2_elite_perctargeted2 <- c()
S2_wild_perctargeted2 <- c()
S3_elite_perctargeted2 <- c()
S3_wild_perctargeted2 <- c()
S4_elite_perctargeted2 <- c()
S4_wild_perctargeted2 <- c()
S5_elite_perctargeted2 <- c()
S5_wild_perctargeted2 <- c()
F1dragtargeted2 <- c()
S1dragtargeted2 <- c()
S2dragtargeted2 <- c()
S3dragtargeted2 <- c()
S4dragtargeted2 <- c()
S5dragtargeted2 <- c()
for(i in 1:100){
  SP$switchGenMap(targeted2_map, centromere = targeted_centromere)
  pop <- randCross2(goodpop[best_good,], badpop[best_bad,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop <- setPheno(pop = pop, h2 = 0.8, simParam = SP)
  targeted2var1[i,] <- genicVarA(pop)
  
  pop1_sel <- selectInd(pop, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_2 <- selectInd(pop1_sel, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  F1_elite_perctargeted2[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel1_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel1_2@geno)) + recurrent_all(badpop[best_bad,]@geno, pop1_sel1_2@geno))
  F1_wild_perctargeted2[i] <- recurrent_all(badpop[best_bad,]@geno, pop1_sel1_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel1_2@geno)) + recurrent_all(badpop[best_bad,]@geno, pop1_sel1_2@geno))
  F1dragtargeted2[i] <- range_around(pullQtlGeno(pop1_sel1_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel1_cross <- randCross2(pop1_sel1_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel1_cross <- setPheno(pop1_sel1_cross, h2 = 0.8, simParam = SP)
  targeted2var2[i,] <- genicVarA(pop1_sel1_2)
  
  pop1_sel2 <- selectInd(pop1_sel1_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_2 <- selectInd(pop1_sel2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S1_elite_perctargeted2[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel2_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel2_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel2_2@geno))
  S1_wild_perctargeted2[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel2_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel2_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel2_2@geno))
  S1dragtargeted2[i] <- range_around(pullQtlGeno(pop1_sel2_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel2_cross <- randCross2(pop1_sel2_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel2_cross <- setPheno(pop1_sel2_cross, h2 = 0.8, simParam = SP)
  targeted2var3[i,] <- genicVarA(pop1_sel2_2)
  
  pop1_sel3 <- selectInd(pop1_sel2_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_2 <- selectInd(pop1_sel3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S2_elite_perctargeted2[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel3_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel3_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel3_2@geno))
  S2_wild_perctargeted2[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel3_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel3_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel3_2@geno))
  S2dragtargeted2[i] <- range_around(pullQtlGeno(pop1_sel3_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel3_cross <- randCross2(pop1_sel3_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel3_cross <- setPheno(pop1_sel3_cross, h2 = 0.8, simParam = SP)
  targeted2var4[i,] <- genicVarA(pop1_sel3_2)
  
  pop1_sel4 <- selectInd(pop1_sel3_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_2 <- selectInd(pop1_sel4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S3_elite_perctargeted2[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel4_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_2@geno))
  S3_wild_perctargeted2[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_2@geno))
  S3dragtargeted2[i] <- range_around(pullQtlGeno(pop1_sel4_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel4_cross <- randCross2(pop1_sel4_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel4_cross <- setPheno(pop1_sel4_cross, h2 = 0.8, simParam = SP)
  targeted2var5[i,] <- genicVarA(pop1_sel4_2)
  
  pop1_sel5 <- selectInd(pop1_sel4_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_2 <- selectInd(pop1_sel5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S4_elite_perctargeted2[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel5_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel5_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel5_2@geno))
  S4_wild_perctargeted2[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel5_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel5_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel5_2@geno))
  S4dragtargeted2[i] <- range_around(pullQtlGeno(pop1_sel5_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel5_cross <- randCross2(pop1_sel5_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel5_cross <- setPheno(pop1_sel5_cross, h2 = 0.8, simParam = SP)
  targeted2var6[i,] <- genicVarA(pop1_sel5_2)
  
  pop1_sel6 <- selectInd(pop1_sel5_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel6_2 <- selectInd(pop1_sel6, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S5_elite_perctargeted2[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel6_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel6_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel6_2@geno))
  S5_wild_perctargeted2[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel6_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel6_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel6_2@geno))
  S5dragtargeted2[i] <- range_around(pullQtlGeno(pop1_sel6_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
}

#introgression scheme analysis

#looking at recurrent parent % in progeny at each generation
F1 <- c(mean(F1_elite_perc), mean(F1_elite_perctargeted), mean(F1_elite_perctargeted2))
BC1 <- c(mean(S1_elite_perc), mean(S1_elite_perctargeted), mean(S1_elite_perctargeted2))
BC2 <- c(mean(S2_elite_perc), mean(S2_elite_perctargeted), mean(S2_elite_perctargeted2))
BC3 <- c(mean(S3_elite_perc), mean(S3_elite_perctargeted), mean(S3_elite_perctargeted2))
BC4 <- c(mean(S4_elite_perc), mean(S4_elite_perctargeted), mean(S4_elite_perctargeted2))

#looking at donor parent
F1_d <- c(mean(F1_wild_perc), mean(F1_wild_perctargeted), mean(F1_wild_perctargeted2))
BC1_d <- c(mean(S1_wild_perc), mean(S1_wild_perctargeted), mean(S1_wild_perctargeted2))
BC2_d <- c(mean(S2_wild_perc), mean(S2_wild_perctargeted), mean(S2_wild_perctargeted2))
BC3_d <- c(mean(S3_wild_perc), mean(S3_wild_perctargeted), mean(S3_wild_perctargeted2))
BC4_d <- c(mean(S4_wild_perc), mean(S4_wild_perctargeted), mean(S4_wild_perctargeted2))

#looking at linkage drag
F1_drag <- c(mean(F1drag), mean(F1dragtargeted), mean(F1dragtargeted2))
BC1_drag <- c(mean(S1drag), mean(S1dragtargeted), mean(S1dragtargeted2))
BC2_drag <- c(mean(S2drag), mean(S2dragtargeted), mean(S2dragtargeted2))
BC3_drag <- c(mean(S3drag), mean(S3dragtargeted), mean(S3dragtargeted2))
BC4_drag <- c(mean(S4drag), mean(S4dragtargeted), mean(S4dragtargeted2))

recur_percent <- cbind(F1, BC1, BC2, BC3, BC4)
all <- as.data.frame(recur_percent)
all2 <- as.data.frame(matrix(data = NA, nrow = 15))
all2$mean <- c(all$F1, all$BC1, all$BC2, all$BC3, all$BC4)
all2$generation <- rep(1:5, size = 5, each = 3)
all2$gen_map <- rep(c("wt","targeted 6-times", "targeted 20-times"), size = 3, each = 1)

donor_percent <- cbind(F1_d, BC1_d, BC2_d, BC3_d, BC4_d)
donor <- as.data.frame(donor_percent)
donor2 <- as.data.frame(matrix(data = NA, nrow = 15))
donor2$mean <- c(donor$F1_d, donor$BC1_d, donor$BC2_d, donor$BC3_d, donor$BC4_d)
donor2$generation <- rep(1:5, size = 5, each = 3)
donor2$gen_map <- rep(c("wt","targeted 6-times", "targeted 20-times"), size = 3, each = 1)

drag <- cbind(F1_drag, BC1_drag, BC2_drag, BC3_drag, BC4_drag)
drag <- as.data.frame(drag)
drag2 <- as.data.frame(matrix(data = NA, nrow = 15))
drag2$mean <- c(drag$F1_drag, drag$BC1_drag, drag$BC2_drag, drag$BC3_drag, drag$BC4_drag)
drag2$generation <- rep(1:5, size = 5, each = 3)
drag2$gen_map <- rep(c("wt","targeted 6-times", "targeted 20-times"), size = 3, each = 1)

group.colors = c("wt" = "#009E73", "targeted 6-times" = "#FF1493", "targeted 20-times" = "#00688B")

library(ggplot2)

#ggplot code to look at mean % of recurrent parent in progeny
recur <- ggplot(all2, aes(x=as.factor(generation), y= mean, group=gen_map)) + 
  geom_line(aes(color = gen_map)) + geom_point(size = 1) + theme_bw() + xlab("Generation") + ylab("Mean Prop. of Recurrent Parent") + ggtitle("Proportion of Recurrent Parent over 4 generations of backcrossing") +
  scale_color_manual(values=group.colors, name = "Genetic Map") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12), legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=12), plot.title = element_text(size=12), legend.position = c(0.8,0.4), legend.key.size = unit(0.3, "lines")) 

#ggplot code to look at mean % of donor parent in progeny
donor <- ggplot(donor2, aes(x=as.factor(generation), y= mean, group=gen_map)) + 
  geom_line(aes(color = gen_map)) + geom_point(size = 1) + theme_bw() + xlab("Generation") + ylab("Mean Prop. of Donor Parent") + ggtitle("Proportion of Donor Parent over 4 generations of backcrossing") +
  scale_color_manual(values=group.colors, name = "Genetic Map") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12), legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=12), plot.title = element_text(size=12), legend.position = c(0.8,0.5), legend.key.size = unit(0.3, "lines")) 

#gpplot to look at % linkage drag after each generation
drag <- ggplot(drag2, aes(x=as.factor(generation), y= mean, group=gen_map)) + 
  geom_line(aes(color = gen_map)) + geom_point(size = 1) + theme_bw() + xlab("Generation") + ylab("Mean Prop. of Linkage Drag from Donor Parent") + ggtitle("Proportion of Linkage Drag over 4 generations of backcrossing") +
  scale_color_manual(values=group.colors, name = "Genetic Map") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12), legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=12), plot.title = element_text(size=12), legend.position = c(0.8,0.7), legend.key.size = unit(0.3, "lines")) 

#binding all the variance matrices together
wtvarall <- c(wtvar2[,1], wtvar3[,1], wtvar4[,1], wtvar5[,1],
              wtvar6[,1])

targetedvarall <- c(targetedvar2[,1], targetedvar3[,1], targetedvar4[,1], targetedvar5[,1],
                targetedvar6[,1])

targeted2varall <- c(targeted2var2[,1], targeted2var3[,1], targeted2var4[,1], targeted2var5[,1],
                    targeted2var6[,1])

#putting all variance data into one data frame
allvar <- cbind(wtvarall, targetedvarall, targeted2varall)
allvar <- as.data.frame(allvar)
allvar2 <- as.data.frame(matrix(data = NA, nrow = 1500))
allvar2$var <- c(allvar$wtvarall, allvar$targetedvarall, allvar$targeted2varall)
allvar2$gen <- rep(1:5, size = 3, each = 100)
allvar2$gen_map <- rep(c("wt","targeted 6-times", "targeted 20-times"), size = 2, each = 500)

#ggplot code to plot the variance
var <- ggplot(allvar2, aes(x=as.factor(gen), y=var, fill=gen_map)) +  ggtitle("Additive Genetic Variance over 4 generations of backcrossing") +
  geom_boxplot() + theme_bw() + xlab("Generation") + ylab("Additive Genetic Variance") + scale_fill_manual(values=group.colors, name = "Genetic Map") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12), legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=12), plot.title = element_text(size=12), legend.position = c(0.8,0.65), legend.key.size = unit(0.3, "lines"))

library(ggpubr)

figure <- ggarrange(recur, donor, drag, var,
                    labels = c("A", "B", "C", "D"),
                    ncol = 2, nrow = 2)
figure
