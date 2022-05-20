library(AlphaSimR)
setwd("/home/rke27/Documents")

#generating mix of repulsion & coupling linkages between QTL
addEff_mix <- runif(300, min = -0.1, max = 0.1)
sum(addEff_mix)

#Reading in QTLs from accurate gene space
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

#SNP sites
segSites <- c(300, 237, 219, 256, 146,  92, 233, 173, 149, 184)

#Reading in genetic maps
final_map <- readRDS("final_map.RData")
ddm1_map <- readRDS("ddm1_map.RData")
zmet2_map <- readRDS("zmet2_map.RData")
recq4_map <- readRDS("recq4_map.RData")
fancm_map <- readRDS("fancm_map.RData")
ideal1_map <- readRDS("ideal1_map.RData")
ideal2_map <- readRDS("ideal2_map.RData")

#Reading in vectors of centromere positions
real_centromere <- c(147.0207, 111.57175, 66.82942, 65, 100,
                     20.490461, 77.60733, 73.698199, 67.54148, 65.189686)
real_centromere <- real_centromere/100

recq4_centromere <- c(355.4204, 285.0390, 168.4086, 162.7559, 239,
                      51.99450, 199.68464, 190.21911, 172.31881, 170.4290)
recq4_centromere <- (recq4_centromere/100)

zmet2_centromere <- c(163.7717, 113.67616, 69.41110, 67.51603, 112, 
                      20.881080, 78.00210, 73.49867, 68.812445, 67.53644)
zmet2_centromere <- zmet2_centromere/100

ideal1_centromere <- real_centromere*10

ideal2_centromere <- c(414.7823, 342.8922, 211.6182, 206.9404, 112,
                       63.02667, 257.1550, 244.20528, 216.75260, 220.07058)
ideal2_centromere <- (ideal2_centromere/100)

fancm_centromere <- c(312.2040, 248.9850, 147.2845, 142.3752, 210,
                      45.38585, 174.29785, 165.53841, 150.52966, 148.83764)
fancm_centromere <- (fancm_centromere/100)

ddm1_centromere <- c(398.6166, 340.9768, 209.0365, 204.1758, 291,
                     63.18590, 256.7602, 244.40481, 215.48164, 219.13686)
ddm1_centromere <- (ddm1_centromere/100)


#Setting up founder population
founderPop <- quickHaplo(nInd = 200, nChr = 10, inbred = TRUE, ploidy = 2L, segSites = segSites)
founderPop@genMap <- final_map
founderPop@centromere <- real_centromere
SP = SimParam$new(founderPop)
SP$addSnpChip(segSites)
SP$setTrackRec(TRUE)
SP$p = 0.15
trait_yield <- new("TraitA", nLoci = 300L, lociPerChr= c(40L, 35L, 35L, 35L, 35L, 30L, 25L, 25L, 20L, 20L),
                   lociLoc = as.integer(c(chr1_QTL, chr2_QTL, chr3_QTL, chr4_QTL, chr5_QTL,
                                          chr6_QTL, chr7_QTL, chr8_QTL, chr9_QTL, chr10_QTL)), addEff = addEff_mix, intercept = 0.1)
SP$resetPed()
SP$manAddTrait(trait_yield)

#Diverging founder population into 2 populations; "good pop" & "bad pop"
#Good pop = selecting FOR polygenic trait; no "resistance" QTLs
pop_good_sel10 <- vector(mode = "list", length = 100)
for(i in 1:100){
  pop_good <- newPop(founderPop, simParam = SP)
  pop_good <- setPheno(pop_good, h2 = 0.8, simParam = SP)
  
  pop_good1 <- randCross(pop_good, nCrosses = 10, nProgeny=10, simParam = SP)
  pop_good1 <- setPheno(pop_good1, h2 = 0.8, simParam = SP)
  
  #pop_good1_sel <- selectInd(pop_good1, nInd = 10, use = "gv", trait = 2, selectTop = FALSE, returnPop = TRUE, simParam = SP)
  pop_good_sel <- selectInd(pop_good1, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good2 <- randCross(pop_good_sel, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good2 <- setPheno(pop_good2, h2 = 0.8, simParam = SP)
  
  #pop_good2_sel2 <- selectInd(pop_good2, nInd = 10, use = "gv", trait = 2, selectTop = FALSE, returnPop = TRUE, simParam = SP)
  pop_good_sel2 <- selectInd(pop_good2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good3 <- randCross(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good3 <- setPheno(pop_good3, h2 = 0.8, simParam = SP)
  
  #pop_good2_sel3 <- selectInd(pop_good3, nInd = 10, use = "gv", trait = 2, selectTop = FALSE, returnPop = TRUE, simParam = SP)
  pop_good_sel3 <- selectInd(pop_good3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good4 <- randCross(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good4 <- setPheno(pop_good4, h2 = 0.8, simParam = SP)
  
  #pop_good2_sel4 <- selectInd(pop_good4, nInd = 10, use = "gv", trait = 2, selectTop = FALSE, returnPop = TRUE, simParam = SP)
  pop_good_sel4 <- selectInd(pop_good4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good5 <- randCross(pop_good_sel4, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good5 <- setPheno(pop_good5, h2 = 0.8, simParam = SP)
  
  #pop_good2_sel5 <- selectInd(pop_good5, nInd = 10, use = "gv", trait = 2, selectTop = FALSE, returnPop = TRUE, simParam = SP)
  pop_good_sel5 <- selectInd(pop_good5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good6 <- randCross(pop_good_sel5, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good6 <- setPheno(pop_good6, h2 = 0.8, simParam = SP)
  
  #pop_good2_sel6 <- selectInd(pop_good6, nInd = 10, use = "gv", trait = 2, selectTop = FALSE, returnPop = TRUE, simParam = SP)
  pop_good_sel6 <- selectInd(pop_good6, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good7 <- randCross(pop_good_sel6, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good7 <- setPheno(pop_good7, h2 = 0.8, simParam = SP)
  
  #pop_good2_sel7 <- selectInd(pop_good7, nInd = 10, use = "gv", trait = 2, selectTop = FALSE, returnPop = TRUE, simParam = SP)
  pop_good_sel7 <- selectInd(pop_good7, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good8 <- randCross(pop_good_sel7, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good8 <- setPheno(pop_good8, h2 = 0.8, simParam = SP)
  
  #pop_good2_sel8 <- selectInd(pop_good8, nInd = 10, use = "gv", trait = 2, selectTop = FALSE, returnPop = TRUE, simParam = SP)
  pop_good_sel8 <- selectInd(pop_good8, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good9 <- randCross(pop_good_sel8, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good9 <- setPheno(pop_good9, h2 = 0.8, simParam = SP)
  
  #pop_good2_sel9 <- selectInd(pop_good9, nInd = 10, use = "gv", trait = 2, selectTop = FALSE, returnPop = TRUE, simParam = SP)
  pop_good_sel9 <- selectInd(pop_good9, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good10 <- randCross(pop_good_sel9, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good10 <- setPheno(pop_good10, h2 = 0.8, simParam = SP)
  
  pop_good_sel10[[i]] <- selectInd(pop_good10, nInd = 5, use = "gv", trait = 1, selectop = TRUE, returnPop = TRUE, simParam = SP)
}

trait <- new("TraitAD", nLoci = 3L, lociPerChr = c(3L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
             lociLoc = c(260L, 272L, 288L), addEff = c(0.2, 0.2, 0.2), domEff = c(0.2, 0.2, 0.2), intercept = 0.1)
SP$resetPed()
SP$manAddTrait(trait)
#Bad pop = not selecting for polygenic trait even though its present; "resistance" QTLs present & selected FOR
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
#this changes every time depending on which row has highest GVs
best_bad <- 465

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

#function to decide which loci to track based on if biallelic between elite and wild parents
donor_sites <- matrix(data = NA, nrow = length(final_map[[1]]), ncol = 2)
recur_sites <- matrix(data = NA, nrow = length(final_map[[1]]), ncol = 2)

#find which wild parent alleles are present
find_donor <- function(good_trait2_geno, bad_trait1_geno){
  for(i in 1:length(good_trait1_geno)){
    if(isTRUE(good_trait1_geno[i] != bad_trait1_geno[i])){
      donor_sites[i,1] <- colnames(bad_trait1_geno)[i]
      donor_sites[i,2] <- bad_trait1_geno[i]
    }
  }
  return(donor_sites)
}
#find which elite parent alleles are present 
find_recurrent <- function(good_trait2_geno, bad_trait1_geno){
  for(i in 1:length(good_trait1_geno)){
    if(isTRUE(good_trait1_geno[i] != bad_trait1_geno[i])){
      recur_sites[i,1] <- colnames(good_trait1_geno)[i]
      recur_sites[i,2] <- good_trait1_geno[i]
    }
  }
  return(recur_sites)
}
recurrentsites <- find_recurrent(good_trait1_geno, bad_trait1_geno)
donorsites <- find_donor(good_trait1_geno, bad_trait1_geno)
sites_to_track <- find_biallelic(good_trait1_geno, bad_trait1_geno, tracked_sites)
sites_to_track <- sort(na.omit(sites_to_track))

#function to look at loci around resistance loci to find if elite or wild parents are present
#to assess linkage drag

progeny_loci <- pullQtlGeno(pop1_sel5_2, trait = 1)[,1:40]
elite_loci <- pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40]
wild_loci <- pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40]

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

range_around(progeny_loci, elite_loci, wild_loci)

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
  S4_elite_perc[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel4_cross@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_cross@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_cross@geno))
  S4_wild_perc[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_cross@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_cross@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_cross@geno))
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


#ddm1
ddm1var1 <- matrix(data = NA, ncol = 2, nrow = 100)
ddm1var2 <- matrix(data = NA, ncol = 2, nrow = 100)
ddm1var3 <- matrix(data = NA, ncol = 2, nrow = 100)
ddm1var4 <- matrix(data = NA, ncol = 2, nrow = 100)
ddm1var5 <- matrix(data = NA, ncol = 2, nrow = 100)
ddm1var6 <- matrix(data = NA, ncol = 2, nrow = 100)
ddm1var7 <- matrix(data = NA, ncol = 2, nrow = 100)
F1_elite_percddm1 <- c()
F1_wild_percddm1 <- c()
S1_elite_percddm1 <- c()
S1_wild_percddm1 <- c()
S2_elite_percddm1 <- c()
S2_wild_percddm1 <- c()
S3_elite_percddm1 <- c()
S3_wild_percddm1 <- c()
S4_elite_percddm1 <- c()
S4_wild_percddm1 <- c()
S5_elite_percddm1 <- c()
S5_wild_percddm1 <- c()
F1dragddm1 <- c()
S1dragddm1 <- c()
S2dragddm1 <- c()
S3dragddm1 <- c()
S4dragddm1 <- c()
S5dragddm1 <- c()
for(i in 1:100){
  SP$switchGenMap(ddm1_map, centromere = ddm1_centromere)
  pop <- randCross2(goodpop[best_good,], badpop[best_bad,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop <- setPheno(pop = pop, h2 = 0.8, simParam = SP)
  ddm1var1[i,] <- genicVarA(pop)
  
  pop1_sel <- selectInd(pop, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_2 <- selectInd(pop1_sel, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  F1_elite_percddm1[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel1_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel1_2@geno)) + recurrent_all(badpop[best_bad,]@geno, pop1_sel1_2@geno))
  F1_wild_percddm1[i] <- recurrent_all(badpop[best_bad,]@geno, pop1_sel1_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel1_2@geno)) + recurrent_all(badpop[best_bad,]@geno, pop1_sel1_2@geno))
  F1dragddm1[i] <- range_around(pullQtlGeno(pop1_sel1_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel1_cross <- randCross2(pop1_sel1_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel1_cross <- setPheno(pop1_sel1_cross, h2 = 0.8, simParam = SP)
  ddm1var2[i,] <- genicVarA(pop1_sel1_2)
  
  pop1_sel2 <- selectInd(pop1_sel1_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_2 <- selectInd(pop1_sel2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S1_elite_percddm1[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel2_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel2_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel2_2@geno))
  S1_wild_percddm1[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel2_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel2_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel2_2@geno))
  S1dragddm1[i] <- range_around(pullQtlGeno(pop1_sel2_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel2_cross <- randCross2(pop1_sel2_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel2_cross <- setPheno(pop1_sel2_cross, h2 = 0.8, simParam = SP)
  ddm1var3[i,] <- genicVarA(pop1_sel2_2)
  
  pop1_sel3 <- selectInd(pop1_sel2_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_2 <- selectInd(pop1_sel3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S2_elite_percddm1[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel3_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel3_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel3_2@geno))
  S2_wild_percddm1[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel3_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel3_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel3_2@geno))
  S2dragddm1[i] <- range_around(pullQtlGeno(pop1_sel3_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel3_cross <- randCross2(pop1_sel3_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel3_cross <- setPheno(pop1_sel3_cross, h2 = 0.8, simParam = SP)
  ddm1var4[i,] <- genicVarA(pop1_sel3_2)
  
  pop1_sel4 <- selectInd(pop1_sel3_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_2 <- selectInd(pop1_sel4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S3_elite_percddm1[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel4_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_2@geno))
  S3_wild_percddm1[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_2@geno))
  S3dragddm1[i] <- range_around(pullQtlGeno(pop1_sel4_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel4_cross <- randCross2(pop1_sel4_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel4_cross <- setPheno(pop1_sel4_cross, h2 = 0.8, simParam = SP)
  ddm1var5[i,] <- genicVarA(pop1_sel4_2)
  
  pop1_sel5 <- selectInd(pop1_sel4_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_2 <- selectInd(pop1_sel5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S4_elite_percddm1[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel4_cross@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_cross@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_cross@geno))
  S4_wild_percddm1[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_cross@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_cross@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_cross@geno))
  S4dragddm1[i] <- range_around(pullQtlGeno(pop1_sel5_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel5_cross <- randCross2(pop1_sel5_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel5_cross <- setPheno(pop1_sel5_cross, h2 = 0.8, simParam = SP)
  ddm1var6[i,] <- genicVarA(pop1_sel5_2)
  
  pop1_sel6 <- selectInd(pop1_sel5_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel6_2 <- selectInd(pop1_sel6, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S5_elite_percddm1[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel6_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel6_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel6_2@geno))
  S5_wild_percddm1[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel6_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel6_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel6_2@geno))
  S5dragddm1[i] <- range_around(pullQtlGeno(pop1_sel6_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
}

#zmet2
zmet2var1 <- matrix(data = NA, ncol = 2, nrow = 100)
zmet2var2 <- matrix(data = NA, ncol = 2, nrow = 100)
zmet2var3 <- matrix(data = NA, ncol = 2, nrow = 100)
zmet2var4 <- matrix(data = NA, ncol = 2, nrow = 100)
zmet2var5 <- matrix(data = NA, ncol = 2, nrow = 100)
zmet2var6 <- matrix(data = NA, ncol = 2, nrow = 100)
zmet2var7 <- matrix(data = NA, ncol = 2, nrow = 100)
F1_elite_perczmet2 <- c()
F1_wild_perczmet2 <- c()
S1_elite_perczmet2 <- c()
S1_wild_perczmet2 <- c()
S2_elite_perczmet2 <- c()
S2_wild_perczmet2 <- c()
S3_elite_perczmet2 <- c()
S3_wild_perczmet2 <- c()
S4_elite_perczmet2 <- c()
S4_wild_perczmet2 <- c()
S5_elite_perczmet2 <- c()
S5_wild_perczmet2 <- c()
F1dragzmet2 <- c()
S1dragzmet2 <- c()
S2dragzmet2 <- c()
S3dragzmet2 <- c()
S4dragzmet2 <- c()
S5dragzmet2 <- c()
for(i in 1:100){
  SP$switchGenMap(zmet2_map, centromere = zmet2_centromere)
  pop <- randCross2(goodpop[best_good,], badpop[best_bad,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop <- setPheno(pop = pop, h2 = 0.8, simParam = SP)
  zmet2var1[i,] <- genicVarA(pop)
  
  pop1_sel <- selectInd(pop, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_2 <- selectInd(pop1_sel, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  F1_elite_perczmet2[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel1_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel1_2@geno)) + recurrent_all(badpop[best_bad,]@geno, pop1_sel1_2@geno))
  F1_wild_perczmet2[i] <- recurrent_all(badpop[best_bad,]@geno, pop1_sel1_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel1_2@geno)) + recurrent_all(badpop[best_bad,]@geno, pop1_sel1_2@geno))
  F1dragzmet2[i] <- range_around(pullQtlGeno(pop1_sel1_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel1_cross <- randCross2(pop1_sel1_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel1_cross <- setPheno(pop1_sel1_cross, h2 = 0.8, simParam = SP)
  zmet2var2[i,] <- genicVarA(pop1_sel1_2)
  
  pop1_sel2 <- selectInd(pop1_sel1_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_2 <- selectInd(pop1_sel2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S1_elite_perczmet2[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel2_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel2_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel2_2@geno))
  S1_wild_perczmet2[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel2_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel2_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel2_2@geno))
  S1dragzmet2[i] <- range_around(pullQtlGeno(pop1_sel2_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel2_cross <- randCross2(pop1_sel2_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel2_cross <- setPheno(pop1_sel2_cross, h2 = 0.8, simParam = SP)
  zmet2var3[i,] <- genicVarA(pop1_sel2_2)
  
  pop1_sel3 <- selectInd(pop1_sel2_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_2 <- selectInd(pop1_sel3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S2_elite_perczmet2[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel3_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel3_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel3_2@geno))
  S2_wild_perczmet2[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel3_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel3_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel3_2@geno))
  S2dragzmet2[i] <- range_around(pullQtlGeno(pop1_sel3_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel3_cross <- randCross2(pop1_sel3_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel3_cross <- setPheno(pop1_sel3_cross, h2 = 0.8, simParam = SP)
  zmet2var4[i,] <- genicVarA(pop1_sel3_2)
  
  pop1_sel4 <- selectInd(pop1_sel3_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_2 <- selectInd(pop1_sel4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S3_elite_perczmet2[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel4_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_2@geno))
  S3_wild_perczmet2[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_2@geno))
  S3dragzmet2[i] <- range_around(pullQtlGeno(pop1_sel4_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel4_cross <- randCross2(pop1_sel4_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel4_cross <- setPheno(pop1_sel4_cross, h2 = 0.8, simParam = SP)
  zmet2var5[i,] <- genicVarA(pop1_sel4_2)
  
  pop1_sel5 <- selectInd(pop1_sel4_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_2 <- selectInd(pop1_sel5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S4_elite_perczmet2[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel4_cross@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_cross@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_cross@geno))
  S4_wild_perczmet2[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_cross@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_cross@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_cross@geno))
  S4dragzmet2[i] <- range_around(pullQtlGeno(pop1_sel5_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel5_cross <- randCross2(pop1_sel5_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel5_cross <- setPheno(pop1_sel5_cross, h2 = 0.8, simParam = SP)
  zmet2var6[i,] <- genicVarA(pop1_sel5_2)
  
  pop1_sel6 <- selectInd(pop1_sel5_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel6_2 <- selectInd(pop1_sel6, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S5_elite_perczmet2[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel6_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel6_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel6_2@geno))
  S5_wild_perczmet2[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel6_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel6_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel6_2@geno))
  S5dragzmet2[i] <- range_around(pullQtlGeno(pop1_sel6_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
}

#recq4
recq4var1 <- matrix(data = NA, ncol = 2, nrow = 100)
recq4var2 <- matrix(data = NA, ncol = 2, nrow = 100)
recq4var3 <- matrix(data = NA, ncol = 2, nrow = 100)
recq4var4 <- matrix(data = NA, ncol = 2, nrow = 100)
recq4var5 <- matrix(data = NA, ncol = 2, nrow = 100)
recq4var6 <- matrix(data = NA, ncol = 2, nrow = 100)
recq4var7 <- matrix(data = NA, ncol = 2, nrow = 100)
F1_elite_percrecq4 <- c()
F1_wild_percrecq4 <- c()
S1_elite_percrecq4 <- c()
S1_wild_percrecq4 <- c()
S2_elite_percrecq4 <- c()
S2_wild_percrecq4 <- c()
S3_elite_percrecq4 <- c()
S3_wild_percrecq4 <- c()
S4_elite_percrecq4 <- c()
S4_wild_percrecq4 <- c()
S5_elite_percrecq4 <- c()
S5_wild_percrecq4 <- c()
F1dragrecq4 <- c()
S1dragrecq4 <- c()
S2dragrecq4 <- c()
S3dragrecq4 <- c()
S4dragrecq4 <- c()
S5dragrecq4 <- c()
for(i in 1:100){
  SP$switchGenMap(recq4_map, centromere = recq4_centromere)
  pop <- randCross2(goodpop[best_good,], badpop[best_bad,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop <- setPheno(pop = pop, h2 = 0.8, simParam = SP)
  recq4var1[i,] <- genicVarA(pop)
  
  pop1_sel <- selectInd(pop, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_2 <- selectInd(pop1_sel, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  F1_elite_percrecq4[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel1_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel1_2@geno)) + recurrent_all(badpop[best_bad,]@geno, pop1_sel1_2@geno))
  F1_wild_percrecq4[i] <- recurrent_all(badpop[best_bad,]@geno, pop1_sel1_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel1_2@geno)) + recurrent_all(badpop[best_bad,]@geno, pop1_sel1_2@geno))
  F1dragrecq4[i] <- range_around(pullQtlGeno(pop1_sel1_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel1_cross <- randCross2(pop1_sel1_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel1_cross <- setPheno(pop1_sel1_cross, h2 = 0.8, simParam = SP)
  recq4var2[i,] <- genicVarA(pop1_sel1_2)
  
  pop1_sel2 <- selectInd(pop1_sel1_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_2 <- selectInd(pop1_sel2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S1_elite_percrecq4[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel2_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel2_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel2_2@geno))
  S1_wild_percrecq4[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel2_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel2_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel2_2@geno))
  S1dragrecq4[i] <- range_around(pullQtlGeno(pop1_sel2_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel2_cross <- randCross2(pop1_sel2_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel2_cross <- setPheno(pop1_sel2_cross, h2 = 0.8, simParam = SP)
  recq4var3[i,] <- genicVarA(pop1_sel2_2)
  
  pop1_sel3 <- selectInd(pop1_sel2_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_2 <- selectInd(pop1_sel3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S2_elite_percrecq4[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel3_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel3_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel3_2@geno))
  S2_wild_percrecq4[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel3_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel3_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel3_2@geno))
  S2dragrecq4[i] <- range_around(pullQtlGeno(pop1_sel3_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel3_cross <- randCross2(pop1_sel3_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel3_cross <- setPheno(pop1_sel3_cross, h2 = 0.8, simParam = SP)
  recq4var4[i,] <- genicVarA(pop1_sel3_2)
  
  pop1_sel4 <- selectInd(pop1_sel3_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_2 <- selectInd(pop1_sel4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S3_elite_percrecq4[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel4_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_2@geno))
  S3_wild_percrecq4[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_2@geno))
  S3dragrecq4[i] <- range_around(pullQtlGeno(pop1_sel4_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel4_cross <- randCross2(pop1_sel4_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel4_cross <- setPheno(pop1_sel4_cross, h2 = 0.8, simParam = SP)
  recq4var5[i,] <- genicVarA(pop1_sel4_2)
  
  pop1_sel5 <- selectInd(pop1_sel4_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_2 <- selectInd(pop1_sel5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S4_elite_percrecq4[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel4_cross@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_cross@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_cross@geno))
  S4_wild_percrecq4[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_cross@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_cross@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_cross@geno))
  S4dragrecq4[i] <- range_around(pullQtlGeno(pop1_sel5_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel5_cross <- randCross2(pop1_sel5_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel5_cross <- setPheno(pop1_sel5_cross, h2 = 0.8, simParam = SP)
  recq4var6[i,] <- genicVarA(pop1_sel5_2)
  
  pop1_sel6 <- selectInd(pop1_sel5_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel6_2 <- selectInd(pop1_sel6, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S5_elite_percrecq4[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel6_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel6_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel6_2@geno))
  S5_wild_percrecq4[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel6_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel6_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel6_2@geno))
  S5dragrecq4[i] <- range_around(pullQtlGeno(pop1_sel6_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
}

#fancm
fancmvar1 <- matrix(data = NA, ncol = 2, nrow = 100)
fancmvar2 <- matrix(data = NA, ncol = 2, nrow = 100)
fancmvar3 <- matrix(data = NA, ncol = 2, nrow = 100)
fancmvar4 <- matrix(data = NA, ncol = 2, nrow = 100)
fancmvar5 <- matrix(data = NA, ncol = 2, nrow = 100)
fancmvar6 <- matrix(data = NA, ncol = 2, nrow = 100)
fancmvar7 <- matrix(data = NA, ncol = 2, nrow = 100)
F1_elite_percfancm <- c()
F1_wild_percfancm <- c()
S1_elite_percfancm <- c()
S1_wild_percfancm <- c()
S2_elite_percfancm <- c()
S2_wild_percfancm <- c()
S3_elite_percfancm <- c()
S3_wild_percfancm <- c()
S4_elite_percfancm <- c()
S4_wild_percfancm <- c()
S5_elite_percfancm <- c()
S5_wild_percfancm <- c()
F1dragfancm <- c()
S1dragfancm <- c()
S2dragfancm <- c()
S3dragfancm <- c()
S4dragfancm <- c()
S5dragfancm <- c()
for(i in 1:100){
  SP$switchGenMap(fancm_map, centromere = fancm_centromere)
  pop <- randCross2(goodpop[best_good,], badpop[best_bad,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop <- setPheno(pop = pop, h2 = 0.8, simParam = SP)
  fancmvar1[i,] <- genicVarA(pop)
  
  pop1_sel <- selectInd(pop, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_2 <- selectInd(pop1_sel, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  F1_elite_percfancm[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel1_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel1_2@geno)) + recurrent_all(badpop[best_bad,]@geno, pop1_sel1_2@geno))
  F1_wild_percfancm[i] <- recurrent_all(badpop[best_bad,]@geno, pop1_sel1_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel1_2@geno)) + recurrent_all(badpop[best_bad,]@geno, pop1_sel1_2@geno))
  F1dragfancm[i] <- range_around(pullQtlGeno(pop1_sel1_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel1_cross <- randCross2(pop1_sel1_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel1_cross <- setPheno(pop1_sel1_cross, h2 = 0.8, simParam = SP)
  fancmvar2[i,] <- genicVarA(pop1_sel1_2)
  
  pop1_sel2 <- selectInd(pop1_sel1_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_2 <- selectInd(pop1_sel2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S1_elite_percfancm[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel2_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel2_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel2_2@geno))
  S1_wild_percfancm[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel2_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel2_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel2_2@geno))
  S1dragfancm[i] <- range_around(pullQtlGeno(pop1_sel2_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  op1_sel2_cross <- randCross2(pop1_sel2_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel2_cross <- setPheno(pop1_sel2_cross, h2 = 0.8, simParam = SP)
  fancmvar3[i,] <- genicVarA(pop1_sel2_2)
  
  pop1_sel3 <- selectInd(pop1_sel2_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_2 <- selectInd(pop1_sel3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S2_elite_percfancm[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel3_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel3_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel3_2@geno))
  S2_wild_percfancm[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel3_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel3_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel3_2@geno))
  S2dragfancm[i] <- range_around(pullQtlGeno(pop1_sel3_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel3_cross <- randCross2(pop1_sel3_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel3_cross <- setPheno(pop1_sel3_cross, h2 = 0.8, simParam = SP)
  fancmvar4[i,] <- genicVarA(pop1_sel3_2)
  
  pop1_sel4 <- selectInd(pop1_sel3_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_2 <- selectInd(pop1_sel4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S3_elite_percfancm[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel4_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_2@geno))
  S3_wild_percfancm[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_2@geno))
  S3dragfancm[i] <- range_around(pullQtlGeno(pop1_sel4_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel4_cross <- randCross2(pop1_sel4_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel4_cross <- setPheno(pop1_sel4_cross, h2 = 0.8, simParam = SP)
  fancmvar5[i,] <- genicVarA(pop1_sel4_2)
  
  pop1_sel5 <- selectInd(pop1_sel4_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_2 <- selectInd(pop1_sel5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S4_elite_percfancm[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel4_cross@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_cross@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_cross@geno))
  S4_wild_percfancm[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_cross@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_cross@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_cross@geno))
  S4dragfancm[i] <- range_around(pullQtlGeno(pop1_sel5_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel5_cross <- randCross2(pop1_sel5_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel5_cross <- setPheno(pop1_sel5_cross, h2 = 0.8, simParam = SP)
  fancmvar6[i,] <- genicVarA(pop1_sel5_2)
  
  pop1_sel6 <- selectInd(pop1_sel5_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel6_2 <- selectInd(pop1_sel6, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S5_elite_percfancm[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel6_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel6_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel6_2@geno))
  S5_wild_percfancm[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel6_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel6_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel6_2@geno))
  S5dragfancm[i] <- range_around(pullQtlGeno(pop1_sel6_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
}

#ideal1 = 10X global increase
ideal1var1 <- matrix(data = NA, ncol = 2, nrow = 100)
ideal1var2 <- matrix(data = NA, ncol = 2, nrow = 100)
ideal1var3 <- matrix(data = NA, ncol = 2, nrow = 100)
ideal1var4 <- matrix(data = NA, ncol = 2, nrow = 100)
ideal1var5 <- matrix(data = NA, ncol = 2, nrow = 100)
ideal1var6 <- matrix(data = NA, ncol = 2, nrow = 100)
ideal1var7 <- matrix(data = NA, ncol = 2, nrow = 100)
F1_elite_perc10X <- c()
F1_wild_perc10X <- c()
S1_elite_perc10X <- c()
S1_wild_perc10X <- c()
S2_elite_perc10X <- c()
S2_wild_perc10X <- c()
S3_elite_perc10X <- c()
S3_wild_perc10X <- c()
S4_elite_perc10X <- c()
S4_wild_perc10X <- c()
S5_elite_perc10X <- c()
S5_wild_perc10X <- c()
F1drag10X <- c()
S1drag10X <- c()
S2drag10X <- c()
S3drag10X <- c()
S4drag10X <- c()
S5drag10X <- c()
for(i in 1:100){
  SP$switchGenMap(ideal1_map, centromere = ideal1_centromere)
  pop <- randCross2(goodpop[best_good,], badpop[best_bad,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop <- setPheno(pop = pop, h2 = 0.8, simParam = SP)
  ideal1var1[i,] <- genicVarA(pop)
  
  pop1_sel <- selectInd(pop, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_2 <- selectInd(pop1_sel, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  F1_elite_perc10X[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel1_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel1_2@geno)) + recurrent_all(badpop[best_bad,]@geno, pop1_sel1_2@geno))
  F1_wild_perc10X[i] <- recurrent_all(badpop[best_bad,]@geno, pop1_sel1_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel1_2@geno)) + recurrent_all(badpop[best_bad,]@geno, pop1_sel1_2@geno))
  F1drag10X[i] <- range_around(pullQtlGeno(pop1_sel1_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel1_cross <- randCross2(pop1_sel1_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel1_cross <- setPheno(pop1_sel1_cross, h2 = 0.8, simParam = SP)
  ideal1var2[i,] <- genicVarA(pop1_sel1_2)
  
  pop1_sel2 <- selectInd(pop1_sel1_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_2 <- selectInd(pop1_sel2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S1_elite_perc10X[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel2_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel2_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel2_2@geno))
  S1_wild_perc10X[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel2_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel2_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel2_2@geno))
  S1drag10X[i] <- range_around(pullQtlGeno(pop1_sel2_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel2_cross <- randCross2(pop1_sel2_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel2_cross <- setPheno(pop1_sel2_cross, h2 = 0.8, simParam = SP)
  ideal1var3[i,] <- genicVarA(pop1_sel2_2)
  
  pop1_sel3 <- selectInd(pop1_sel2_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_2 <- selectInd(pop1_sel3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S2_elite_perc10X[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel3_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel3_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel3_2@geno))
  S2_wild_perc10X[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel3_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel3_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel3_2@geno))
  S2drag10X[i] <- range_around(pullQtlGeno(pop1_sel3_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel3_cross <- randCross2(pop1_sel3_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel3_cross <- setPheno(pop1_sel3_cross, h2 = 0.8, simParam = SP)
  ideal1var4[i,] <- genicVarA(pop1_sel3_2)
  
  pop1_sel4 <- selectInd(pop1_sel3_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_2 <- selectInd(pop1_sel4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S3_elite_perc10X[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel4_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_2@geno))
  S3_wild_perc10X[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_2@geno))
  S3drag10X[i] <- range_around(pullQtlGeno(pop1_sel4_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel4_cross <- randCross2(pop1_sel4_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel4_cross <- setPheno(pop1_sel4_cross, h2 = 0.8, simParam = SP)
  ideal1var5[i,] <- genicVarA(pop1_sel4_2)
  
  pop1_sel5 <- selectInd(pop1_sel4_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_2 <- selectInd(pop1_sel5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S4_elite_perc10X[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel4_cross@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_cross@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_cross@geno))
  S4_wild_perc10X[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_cross@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_cross@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_cross@geno))
  S4drag10X[i] <- range_around(pullQtlGeno(pop1_sel5_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel5_cross <- randCross2(pop1_sel5_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel5_cross <- setPheno(pop1_sel5_cross, h2 = 0.8, simParam = SP)
  ideal1var6[i,] <- genicVarA(pop1_sel5_2)
  
  pop1_sel6 <- selectInd(pop1_sel5_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel6_2 <- selectInd(pop1_sel6, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S5_elite_perc10X[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel6_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel6_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel6_2@geno))
  S5_wild_perc10X[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel6_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel6_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel6_2@geno))
  S5drag10X[i] <- range_around(pullQtlGeno(pop1_sel6_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
}

#ideal2 = ddm1/zmet2 double mutant
ideal2var1 <- matrix(data = NA, ncol = 2, nrow = 100)
ideal2var2 <- matrix(data = NA, ncol = 2, nrow = 100)
ideal2var3 <- matrix(data = NA, ncol = 2, nrow = 100)
ideal2var4 <- matrix(data = NA, ncol = 2, nrow = 100)
ideal2var5 <- matrix(data = NA, ncol = 2, nrow = 100)
ideal2var6 <- matrix(data = NA, ncol = 2, nrow = 100)
ideal2var7 <- matrix(data = NA, ncol = 2, nrow = 100)
F1_elite_percdz <- c()
F1_wild_percdz <- c()
S1_elite_percdz <- c()
S1_wild_percdz <- c()
S2_elite_percdz <- c()
S2_wild_percdz <- c()
S3_elite_percdz <- c()
S3_wild_percdz <- c()
S4_elite_percdz <- c()
S4_wild_percdz <- c()
S5_elite_percdz <- c()
S5_wild_percdz <- c()
F1dragdz <- c()
S1dragdz <- c()
S2dragdz <- c()
S3dragdz <- c()
S4dragdz <- c()
S5dragdz <- c()
for(i in 1:100){
  SP$switchGenMap(ideal2_map, centromere = ideal2_centromere)
  pop <- randCross2(goodpop[best_good,], badpop[best_bad,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop <- setPheno(pop = pop, h2 = 0.8, simParam = SP)
  ideal2var1[i,] <- genicVarA(pop)
  
  pop1_sel <- selectInd(pop, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_2 <- selectInd(pop1_sel, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  F1_elite_percdz[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel1_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel1_2@geno)) + recurrent_all(badpop[best_bad,]@geno, pop1_sel1_2@geno))
  F1_wild_percdz[i] <- recurrent_all(badpop[best_bad,]@geno, pop1_sel1_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel1_2@geno)) + recurrent_all(badpop[best_bad,]@geno, pop1_sel1_2@geno))
  F1dragdz[i] <- range_around(pullQtlGeno(pop1_sel1_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel1_cross <- randCross2(pop1_sel1_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel1_cross <- setPheno(pop1_sel1_cross, h2 = 0.8, simParam = SP)
  ideal2var2[i,] <- genicVarA(pop1_sel1_2)
  
  pop1_sel2 <- selectInd(pop1_sel1_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_2 <- selectInd(pop1_sel2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S1_elite_percdz[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel2_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel2_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel2_2@geno))
  S1_wild_percdz[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel2_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel2_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel2_2@geno))
  S1dragdz[i] <- range_around(pullQtlGeno(pop1_sel2_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel2_cross <- randCross2(pop1_sel2_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel2_cross <- setPheno(pop1_sel2_cross, h2 = 0.8, simParam = SP)
  ideal2var3[i,] <- genicVarA(pop1_sel2_2)
  
  pop1_sel3 <- selectInd(pop1_sel2_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_2 <- selectInd(pop1_sel3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S2_elite_percdz[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel3_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel3_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel3_2@geno))
  S2_wild_percdz[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel3_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel3_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel3_2@geno))
  S2dragdz[i] <- range_around(pullQtlGeno(pop1_sel3_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel3_cross <- randCross2(pop1_sel3_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel3_cross <- setPheno(pop1_sel3_cross, h2 = 0.8, simParam = SP)
  ideal2var4[i,] <- genicVarA(pop1_sel3_2)
  
  pop1_sel4 <- selectInd(pop1_sel3_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_2 <- selectInd(pop1_sel4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S3_elite_percdz[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel4_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_2@geno))
  S3_wild_percdz[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_2@geno))
  S3dragdz[i] <- range_around(pullQtlGeno(pop1_sel4_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel4_cross <- randCross2(pop1_sel4_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel4_cross <- setPheno(pop1_sel4_cross, h2 = 0.8, simParam = SP)
  ideal2var5[i,] <- genicVarA(pop1_sel4_2)
  
  pop1_sel5 <- selectInd(pop1_sel4_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_2 <- selectInd(pop1_sel5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S4_elite_percdz[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel4_cross@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_cross@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_cross@geno))
  S4_wild_percdz[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_cross@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel4_cross@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel4_cross@geno))
  S4dragdz[i] <- range_around(pullQtlGeno(pop1_sel5_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
  pop1_sel5_cross <- randCross2(pop1_sel5_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel5_cross <- setPheno(pop1_sel5_cross, h2 = 0.8, simParam = SP)
  ideal2var6[i,] <- genicVarA(pop1_sel5_2)
  
  pop1_sel6 <- selectInd(pop1_sel5_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel6_2 <- selectInd(pop1_sel6, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  S5_elite_percdz[i] <- recurrent_all(goodpop[best_good,]@geno, pop1_sel6_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel6_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel6_2@geno))
  S5_wild_percdz[i] <- recurrent_all(badpop[best_bad,]@geno,  pop1_sel6_2@geno)/((recurrent_all(goodpop[best_good,]@geno, pop1_sel6_2@geno)) + recurrent_all(badpop[best_bad,]@geno,  pop1_sel6_2@geno))
  S5dragdz[i] <- range_around(pullQtlGeno(pop1_sel1_2, trait = 1)[,1:40], pullQtlGeno(goodpop[best_good,], trait = 1)[,1:40], pullQtlGeno(badpop[best_bad,], trait = 1)[,1:40])
}

#introgression scheme analysis

#looking at recurrent parent % in progeny at each generation
F1 <- c(mean(F1_elite_perc), mean(F1_elite_percddm1), mean(F1_elite_perczmet2), mean(F1_elite_percfancm), mean(F1_elite_percrecq4), mean(F1_elite_perc10X), mean(F1_elite_percdz))

BC1 <- c(mean(S1_elite_perc), mean(S1_elite_percddm1), mean(S1_elite_perczmet2), mean(S1_elite_percfancm), mean(S1_elite_percrecq4), mean(S1_elite_perc10X), mean(S1_elite_percdz))

BC2 <- c(mean(S2_elite_perc), mean(S2_elite_percddm1), mean(S2_elite_perczmet2), mean(S2_elite_percfancm), mean(S2_elite_percrecq4), mean(S2_elite_perc10X), mean(S2_elite_percdz))

BC3 <- c(mean(S3_elite_perc), mean(S3_elite_percddm1), mean(S3_elite_perczmet2), mean(S3_elite_percfancm), mean(S3_elite_percrecq4), mean(S3_elite_perc10X), mean(S3_elite_percdz))

BC4 <- c(mean(S4_elite_perc), mean(S4_elite_percddm1), mean(S4_elite_perczmet2), mean(S4_elite_percfancm), mean(S4_elite_percrecq4), mean(S4_elite_perc10X), mean(S4_elite_percdz))

#looking at donor parent
F1_d <- c(mean(F1_wild_perc), mean(F1_wild_percddm1), mean(F1_wild_perczmet2), mean(F1_wild_percfancm), mean(F1_wild_percrecq4), mean(F1_wild_perc10X), mean(F1_wild_percdz))

BC1_d <- c(mean(S1_wild_perc), mean(S1_wild_percddm1), mean(S1_wild_perczmet2), mean(S1_wild_percfancm), mean(S1_wild_percrecq4), mean(S1_wild_perc10X), mean(S1_wild_percdz))

BC2_d <- c(mean(S2_wild_perc), mean(S2_wild_percddm1), mean(S2_wild_perczmet2), mean(S2_wild_percfancm), mean(S2_wild_percrecq4), mean(S2_wild_perc10X), mean(S2_wild_percdz))

BC3_d <- c(mean(S3_wild_perc), mean(S3_wild_percddm1), mean(S3_wild_perczmet2), mean(S3_wild_percfancm), mean(S3_wild_percrecq4), mean(S3_wild_perc10X), mean(S3_wild_percdz))

BC4_d <- c(mean(S4_wild_perc), mean(S4_wild_percddm1), mean(S4_wild_perczmet2), mean(S4_wild_percfancm), mean(S4_wild_percrecq4), mean(S4_wild_perc10X), mean(S4_wild_percdz))

#looking at linkage drag
F1_drag <- c(mean(F1drag), mean(F1dragddm1), mean(F1dragzmet2), mean(F1dragfancm), mean(F1dragrecq4), mean(F1drag10X), mean(F1dragdz))

BC1_drag <- c(mean(S1drag), mean(S1dragddm1), mean(S1dragzmet2), mean(S1dragfancm), mean(S1dragrecq4), mean(S1drag10X), mean(S1dragdz))

BC2_drag <- c(mean(S2drag), mean(S2dragddm1), mean(S2dragzmet2), mean(S2dragfancm), mean(S2dragrecq4), mean(S2drag10X), mean(S2dragdz))

BC3_drag <- c(mean(S3drag), mean(S3dragddm1), mean(S3dragzmet2), mean(S3dragfancm), mean(S3dragrecq4), mean(S3drag10X), mean(S3dragdz))

BC4_drag <- c(mean(S4drag), mean(S4dragddm1), mean(S4dragzmet2), mean(S4dragfancm), mean(S4dragrecq4), mean(S4drag10X), mean(S4dragdz))

recur_percent <- cbind(F1, BC1, BC2, BC3, BC4)
all <- as.data.frame(recur_percent)
all2 <- as.data.frame(matrix(data = NA, nrow = 35))
all2$mean <- c(all$F1, all$BC1, all$BC2, all$BC3, all$BC4)
all2$generation <- rep(1:5, size = 5, each = 7)
all2$gen_map <- rep(c("wt","ddm1", "zmet2", "recq4", "fancm", "10X", "ddm1/zmet2"), size = 7, each = 1)

donor_percent <- cbind(F1_d, BC1_d, BC2_d, BC3_d, BC4_d)
donor <- as.data.frame(donor_percent)
donor2 <- as.data.frame(matrix(data = NA, nrow = 35))
donor2$mean <- c(donor$F1_d, donor$BC1_d, donor$BC2_d, donor$BC3_d, donor$BC4_d)
donor2$generation <- rep(1:5, size = 5, each = 7)
donor2$gen_map <- rep(c("wt","ddm1", "zmet2", "recq4", "fancm", "10X", "ddm1/zmet2"), size = 7, each = 1)

drag <- cbind(F1_drag, BC1_drag, BC2_drag, BC3_drag, BC4_drag)
drag <- as.data.frame(drag)
drag2 <- as.data.frame(matrix(data = NA, nrow = 35))
drag2$mean <- c(drag$F1_drag, drag$BC1_drag, drag$BC2_drag, drag$BC3_drag, drag$BC4_drag)
drag2$generation <- rep(1:5, size = 5, each = 7)
drag2$gen_map <- rep(c("wt","ddm1", "zmet2", "recq4", "fancm", "10X", "ddm1/zmet2"), size = 7, each = 1)

#looking at introgression of resistance loci
trait2 <- cbind(wt2, ddm12)
trait2 <- as.data.frame(trait2)
trait22 <- as.data.frame(matrix(data = NA, nrow = 35))
trait22$gv <- c(trait2$wt2, trait2$ddm12, trait2$zmet22, trait2$recq42, trait2$ideal12, trait2$ideal22, trait2$fancm2)
trait22$generation <- rep(1:5, size = 7, each = 1)
trait22$gen_map <- rep(c("wt","ddm1", "zmet2", "recq4", "ideal1", "ideal2", "fancm"), size = 7, each = 1)

group.colors = c("ddm1" = "#E69F00", "recq4" = "#56B4E9", "wt" = "#009E73", "zmet2" = "#F0E442",
                 "fancm" = "#0072B2", "10X" = "#D55E00", "ddm1/zmet2" = "#CC79A7")

#ggplot code to look at mean % of recurrent parent in progeny
recur <- ggplot(all2, aes(x=as.factor(generation), y= mean, group=gen_map)) + 
  geom_line(aes(color = gen_map)) + geom_point(size = 1) + theme_bw() + xlab("Generations") + ylab("Mean % of Recurrent Parent") + ggtitle("% of Recurrent Parent after 4 gens of backcrossing") +
  scale_color_manual(values=group.colors, name = "Genetic Map", labels = c("10X", "ddm1", "ddm1/zmet2", "fancm",
                                                                           "recq4", "WT", "zmet2")) +
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=11), legend.title = element_text(size=11), #change legend title font size
        legend.text = element_text(size=11), plot.title = element_text(size=12), legend.position = c(0.8,0.4), legend.key.size = unit(0.3, "lines")) 

#ggplot code to look at mean % of donor parent in progeny
donor <- ggplot(donor2, aes(x=as.factor(generation), y= mean, group=gen_map)) + 
  geom_line(aes(color = gen_map)) + geom_point(size = 1) + theme_bw() + xlab("Generations") + ylab("Mean % of Donor Parent") + ggtitle("% of Donor Parent after 4 gens of backcrossing") +
  scale_color_manual(values=group.colors, name = "Genetic Map", labels = c("10X", "ddm1", "ddm1/zmet2", "fancm",
                                                                           "recq4", "WT", "zmet2")) +
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=11), legend.title = element_text(size=11), #change legend title font size
        legend.text = element_text(size=11), plot.title = element_text(size=12), legend.position = c(0.8,0.5), legend.key.size = unit(0.3, "lines")) 

#gpplot to look at % linkage drag after each generation
drag <- ggplot(drag2, aes(x=as.factor(generation), y= mean, group=gen_map)) + 
  geom_line(aes(color = gen_map)) + geom_point(size = 1) + theme_bw() + xlab("Generations") + ylab("Mean % of Linkage Drag from Donor Parent") + ggtitle("% of Linkage Drag over 4 generations of backcrossing") +
  scale_color_manual(values=group.colors, name = "Genetic Map", labels = c("10X", "ddm1", "ddm1/zmet2", "fancm",
                                                                           "recq4", "WT", "zmet2")) +
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=11), legend.title = element_text(size=11), #change legend title font size
        legend.text = element_text(size=11), plot.title = element_text(size=12), legend.position = c(0.8,0.7), legend.key.size = unit(0.3, "lines")) 

#binding all the variance matrices together
wtvarall <- c(wtvar1[,1], wtvar2[,1], wtvar3[,1], wtvar4[,1], wtvar5[,1],
              wtvar6[,1])

ddm1varall <- c(ddm1var1[,1], ddm1var2[,1], ddm1var3[,1], ddm1var4[,1], ddm1var5[,1],
                ddm1var6[,1])

zmet2varall <- c(zmet2var1[,1], zmet2var2[,1], zmet2var3[,1], zmet2var4[,1], zmet2var5[,1],
                 zmet2var6[,1])

recq4varall <- c(recq4var1[,1], recq4var2[,1], recq4var3[,1], recq4var4[,1], recq4var5[,1],
                 recq4var6[,1])

ideal1varall <- c(ideal1var1[,1], ideal1var2[,1], ideal1var3[,1], ideal1var4[,1], ideal1var5[,1],
                  ideal1var6[,1])

ddm1_zmet2varall <- c(ideal2var1[,1], ideal2var2[,1], ideal2var3[,1], ideal2var4[,1], ideal2var5[,1],
                      ideal2var6[,1])

fancmvarall <- c(fancmvar1[,1], fancmvar2[,1], fancmvar3[,1], fancmvar4[,1], fancmvar5[,1],
                 fancmvar6[,1])

#putting all variance data into one data frame
allvar <- cbind(wtvarall, ddm1varall, zmet2varall, recq4varall, ideal1varall, ddm1_zmet2varall, fancmvarall)
allvar <- as.data.frame(allvar)
allvar2 <- as.data.frame(matrix(data = NA, nrow = 4200))
allvar2$var <- c(allvar$wtvarall, allvar$ddm1varall, allvar$zmet2varall, allvar$recq4varall, allvar$ideal1varall, allvar$ddm1_zmet2varall, allvar$fancmvarall)
allvar2$gen <- rep(1:6, size = 7, each = 100)
allvar2$gen_map <- rep(c("wt","ddm1", "zmet2", "recq4", "10X", "ddm1/zmet2", "fancm"), size = 7, each = 600)

#ggplot code to plot the variance
var <- ggplot(allvar2, aes(x=as.factor(gen), y=var, fill=gen_map)) +  ggtitle("Additive genetic variance over 4 gens of backcrossing") +
  geom_boxplot() + theme_bw() + xlab("Generations") + ylab("Additive Genetic Variances") + scale_fill_manual(values=group.colors, name = "Genetic Map", labels = c("10X", "ddm1", "ddm1/zmet2", "fancm",
                                                                                                                                                                       "recq4", "WT", "zmet2")) +
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=11), legend.title = element_text(size=11), #change legend title font size
        legend.text = element_text(size=11), plot.title = element_text(size=12), legend.position = c(0.8,0.65), legend.key.size = unit(0.3, "lines"))

library(ggpubr)

figure <- ggarrange(recur, donor, drag, var,
                    labels = c("A", "B", "C", "D"),
                    ncol = 2, nrow = 2)
figure
