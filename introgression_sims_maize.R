library(AlphaSimR)

#generating mix of repulsion & coupling linkages between QTL
addEff_mix <- runif(300, min = -0.1, max = 0.1)
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
pop_bad_sel10 <- vector(mode = "list", length = 100)
trait <- new("TraitAD", nLoci = 3L, lociPerChr = c(3L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
             lociLoc = c(10L, 200L, 290L), addEff = c(0.1,0.1,0.1), domEff = c(0.2,0.2,0.2), intercept = 0.1)
SP$resetPed()
SP$manAddTrait(trait)
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

pop_bad_sel10 <- vector(mode = "list", length = 100)
trait <- new("TraitAD", nLoci = 3L, lociPerChr = c(0L, 0L, 0L, 3L, 0L, 0L, 0L, 0L, 0L, 0L),
             lociLoc = c(1L,130L,133L), addEff = c(0.1,0.1,0.1), domEff = c(0.2,0.2,0.2), intercept = 0.1)
SP$resetPed()
SP$manAddTrait(trait)
for(i in 1:100){
  pop_bad <- newPop(founderPop, simParam = SP)
  pop_bad <- setPheno(pop_bad, h2 = 0.8, simParam = SP)
  
  pop_bad1 <- randCross(pop_bad, nCrosses = 10, nProgeny=10, simParam = SP)
  pop_bad1 <- setPheno(pop_bad1, h2 = 0.8, simParam = SP)
  
  pop_bad2 <- selectCross(pop_bad1, nInd = 5, use = "gv", trait = 2, selectTop = TRUE, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_bad2 <- setPheno(pop_bad2, h2 = 0.8, simParam = SP)
  
  pop_bad3 <- selectCross(pop_bad2, nInd = 5, use = "gv", trait = 2, selectTop = TRUE, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_bad3 <- setPheno(pop_bad3, h2 = 0.8, simParam = SP)
  
  pop_bad4 <- selectCross(pop_bad3, nInd = 5, use = "gv", trait = 2, selectTop = TRUE, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_bad4 <- setPheno(pop_bad4, h2 = 0.8, simParam = SP)
  
  pop_bad5 <- selectCross(pop_bad4, nInd = 5, use = "gv", trait = 2, selectTop = TRUE, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_bad5 <- setPheno(pop_bad5, h2 = 0.8, simParam = SP)
  
  pop_bad6 <- selectCross(pop_bad5, nInd = 5, use = "gv", trait = 2, selectTop = TRUE, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_bad6 <- setPheno(pop_bad6, h2 = 0.8, simParam = SP)
  
  pop_bad7 <- selectCross(pop_bad6, nInd = 5, use = "gv", trait = 2, selectTop = TRUE, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_bad7 <- setPheno(pop_bad7, h2 = 0.8, simParam = SP)
  
  pop_bad8 <- selectCross(pop_bad7, nInd = 5, use = "gv", trait = 2, selectTop = TRUE, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_bad8 <- setPheno(pop_bad8, h2 = 0.8, simParam = SP)
  
  pop_bad9 <- selectCross(pop_bad8, nInd = 5, use = "gv", trait = 2, selectTop = TRUE, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_bad9 <- setPheno(pop_bad9, h2 = 0.8, simParam = SP)
  
  pop_bad10 <- selectCross(pop_bad9, nInd = 5, use = "gv", trait = 2, selectTop = TRUE, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_bad10 <- setPheno(pop_bad10, h2 = 0.8, simParam = SP)
  
  pop_bad_sel10[[i]] <- selectInd(pop_bad10, nInd = 5, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
}

goodpop = mergePops(pop_good_sel10)
badpop = mergePops(pop_bad_sel10)

##Backcrossing scheme
#selecting best individual from "good pop" and individual with resistance QTLs in "bad pop"
best_good <- which(goodpop@gv == max(goodpop@gv), arr.ind = FALSE)
best_bad <- which(badpop@gv == max(badpop@gv), arr.ind = FALSE)
#this changes every time depending on which row has highest GVs
best_bad <- 432

#wild-type first
pop1_sel2_gv <- matrix(data = NA, nrow = 100, ncol = 10)
pop1_sel3_gv <- matrix(data = NA, nrow = 100, ncol = 10)
pop1_sel4_gv <- matrix(data = NA, nrow = 100, ncol = 10)
pop1_sel5_gv <- matrix(data = NA, nrow = 100, ncol = 10)
pop1_sel6_gv <- matrix(data = NA, nrow = 100, ncol = 10)
wtvar1 <- matrix(data = NA, ncol = 2, nrow = 200)
wtvar2 <- matrix(data = NA, ncol = 2, nrow = 200)
wtvar3 <- matrix(data = NA, ncol = 2, nrow = 200)
wtvar4 <- matrix(data = NA, ncol = 2, nrow = 200)
wtvar5 <- matrix(data = NA, ncol = 2, nrow = 200)
wtvar6 <- matrix(data = NA, ncol = 2, nrow = 200)
wtvar7 <- matrix(data = NA, ncol = 2, nrow = 200)
for(i in 1:100){
  SP$switchGenMap(final_map, centromere = real_centromere)
  pop <- randCross2(badpop[best_bad,], goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop <- setPheno(pop = pop, h2 = 0.8, simParam = SP)
  wtvar1[i] <- genicVarA(pop)
  
  pop1_sel <- selectInd(pop, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_2 <- selectInd(pop1_sel, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_cross <- randCross2(pop1_sel1_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel1_cross <- setPheno(pop1_sel1_cross, h2 = 0.8, simParam = SP)
  wtvar2[i] <- genicVarA(pop1_sel1_2)
  
  pop1_sel2_gv[i,] <- gv(pop1_sel1_2)
  pop1_sel2 <- selectInd(pop1_sel1_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_2 <- selectInd(pop1_sel2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_cross <- randCross2(pop1_sel2_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel2_cross <- setPheno(pop1_sel2_cross, h2 = 0.8, simParam = SP)
  wtvar3[i] <- genicVarA(pop1_sel2_2)
  
  pop1_sel3_gv[i,] <- gv(pop1_sel2_2)
  pop1_sel3 <- selectInd(pop1_sel2_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_2 <- selectInd(pop1_sel3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_cross <- randCross2(pop1_sel3_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel3_cross <- setPheno(pop1_sel3_cross, h2 = 0.8, simParam = SP)
  wtvar4[i] <- genicVarA(pop1_sel3_2)
  
  pop1_sel4_gv[i,] <- gv(pop1_sel3_2)
  pop1_sel4 <- selectInd(pop1_sel3_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_2 <- selectInd(pop1_sel4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_cross <- randCross2(pop1_sel4_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel4_cross <- setPheno(pop1_sel4_cross, h2 = 0.8, simParam = SP)
  wtvar5[i] <- genicVarA(pop1_sel4_2)
  
  pop1_sel5_gv[i,] <- gv(pop1_sel4_2)
  pop1_sel5 <- selectInd(pop1_sel4_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_2 <- selectInd(pop1_sel5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_cross <- randCross2(pop1_sel5_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel5_cross <- setPheno(pop1_sel5_cross, h2 = 0.8, simParam = SP)
  wtvar6[i] <- genicVarA(pop1_sel5_2)
  
  pop1_sel6_gv[i,] <- gv(pop1_sel5_2)
}

#ddm1
ddm1_sel2_gv <- matrix(data = NA, nrow = 100, ncol = 10)
ddm1_sel3_gv <- matrix(data = NA, nrow = 100, ncol = 10)
ddm1_sel4_gv <- matrix(data = NA, nrow = 100, ncol = 10)
ddm1_sel5_gv <- matrix(data = NA, nrow = 100, ncol = 10)
ddm1_sel6_gv <- matrix(data = NA, nrow = 100, ncol = 10)
ddm1_sel7_gv <- matrix(data = NA, nrow = 100, ncol = 10)
ddm1var1 <- matrix(data = NA, ncol = 2, nrow = 200)
ddm1var2 <- matrix(data = NA, ncol = 2, nrow = 200)
ddm1var3 <- matrix(data = NA, ncol = 2, nrow = 200)
ddm1var4 <- matrix(data = NA, ncol = 2, nrow = 200)
ddm1var5 <- matrix(data = NA, ncol = 2, nrow = 200)
ddm1var6 <- matrix(data = NA, ncol = 2, nrow = 200)
ddm1var7 <- matrix(data = NA, ncol = 2, nrow = 200)
for(i in 1:100){
  SP$switchGenMap(ddm1_map, centromere = ddm1_centromere)
  pop <- randCross2(goodpop[best_good,], badpop[best_bad,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop <- setPheno(pop = pop, h2 = 0.8, simParam = SP)
  ddm1var1[i] <- genicVarA(pop)
  
  pop1_sel <- selectInd(pop, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_2 <- selectInd(pop1_sel, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_cross <- randCross2(pop1_sel1_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel1_cross <- setPheno(pop1_sel1_cross, h2 = 0.8, simParam = SP)
  ddm1var2[i] <- genicVarA(pop1_sel1_2)
  
  ddm1_sel2_gv[i,] <- gv(pop1_sel1_2)
  pop1_sel2 <- selectInd(pop1_sel1_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_2 <- selectInd(pop1_sel2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_cross <- randCross2(pop1_sel2_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel2_cross <- setPheno(pop1_sel2_cross, h2 = 0.8, simParam = SP)
  ddm1var3[i] <- genicVarA(pop1_sel2_2)
  
  ddm1_sel3_gv[i,] <- gv(pop1_sel2_2)
  pop1_sel3 <- selectInd(pop1_sel2_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_2 <- selectInd(pop1_sel3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_cross <- randCross2(pop1_sel3_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel3_cross <- setPheno(pop1_sel3_cross, h2 = 0.8, simParam = SP)
  ddm1var4[i] <- genicVarA(pop1_sel3_2)
  
  ddm1_sel4_gv[i,] <- gv(pop1_sel3_2)
  pop1_sel4 <- selectInd(pop1_sel3_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_2 <- selectInd(pop1_sel4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_cross <- randCross2(pop1_sel4_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel4_cross <- setPheno(pop1_sel4_cross, h2 = 0.8, simParam = SP)
  ddm1var5[i] <- genicVarA(pop1_sel4_2)
  
  ddm1_sel5_gv[i,] <- gv(pop1_sel4_2)
  pop1_sel5 <- selectInd(pop1_sel4_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_2 <- selectInd(pop1_sel5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_cross <- randCross2(pop1_sel5_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel5_cross <- setPheno(pop1_sel5_cross, h2 = 0.8, simParam = SP)
  ddm1var6[i] <- genicVarA(pop1_sel5_2)
  
  ddm1_sel6_gv[i,] <- gv(pop1_sel5_2)
}

#zmet2
zmet2_sel2_gv <- matrix(data = NA, nrow = 100, ncol = 10)
zmet2_sel3_gv <- matrix(data = NA, nrow = 100, ncol = 10)
zmet2_sel4_gv <- matrix(data = NA, nrow = 100, ncol = 10)
zmet2_sel5_gv <- matrix(data = NA, nrow = 100, ncol = 10)
zmet2_sel6_gv <- matrix(data = NA, nrow = 100, ncol = 10)
zmet2_sel7_gv <- matrix(data = NA, nrow = 100, ncol = 10)
zmet2var1 <- matrix(data = NA, ncol = 2, nrow = 200)
zmet2var2 <- matrix(data = NA, ncol = 2, nrow = 200)
zmet2var3 <- matrix(data = NA, ncol = 2, nrow = 200)
zmet2var4 <- matrix(data = NA, ncol = 2, nrow = 200)
zmet2var5 <- matrix(data = NA, ncol = 2, nrow = 200)
zmet2var6 <- matrix(data = NA, ncol = 2, nrow = 200)
zmet2var7 <- matrix(data = NA, ncol = 2, nrow = 200)
for(i in 1:100){
  SP$switchGenMap(zmet2_map, centromere = zmet2_centromere)
  pop <- randCross2(goodpop[best_good,], badpop[best_bad,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop <- setPheno(pop = pop, h2 = 0.8, simParam = SP)
  zmet2var1[i] <- genicVarA(pop)
  
  pop1_sel <- selectInd(pop, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_2 <- selectInd(pop1_sel, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_cross <- randCross2(pop1_sel1_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel1_cross <- setPheno(pop1_sel1_cross, h2 = 0.8, simParam = SP)
  zmet2var2[i] <- genicVarA(pop1_sel1_2)
  
  zmet2_sel2_gv[i,] <- gv(pop1_sel1_2)
  pop1_sel2 <- selectInd(pop1_sel1_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_2 <- selectInd(pop1_sel2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_cross <- randCross2(pop1_sel2_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel2_cross <- setPheno(pop1_sel2_cross, h2 = 0.8, simParam = SP)
  zmet2var3[i] <- genicVarA(pop1_sel2_2)
  
  zmet2_sel3_gv[i,] <- gv(pop1_sel2_2)
  pop1_sel3 <- selectInd(pop1_sel2_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_2 <- selectInd(pop1_sel3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_cross <- randCross2(pop1_sel3_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel3_cross <- setPheno(pop1_sel3_cross, h2 = 0.8, simParam = SP)
  zmet2var4[i] <- genicVarA(pop1_sel3_2)
  
  zmet2_sel4_gv[i,] <- gv(pop1_sel3_2)
  pop1_sel4 <- selectInd(pop1_sel3_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_2 <- selectInd(pop1_sel4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_cross <- randCross2(pop1_sel4_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel4_cross <- setPheno(pop1_sel4_cross, h2 = 0.8, simParam = SP)
  zmet2var5[i] <- genicVarA(pop1_sel4_2)
  
  zmet2_sel5_gv[i,] <- gv(pop1_sel4_2)
  pop1_sel5 <- selectInd(pop1_sel4_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_2 <- selectInd(pop1_sel5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_cross <- randCross2(pop1_sel5_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel5_cross <- setPheno(pop1_sel5_cross, h2 = 0.8, simParam = SP)
  zmet2var6[i] <- genicVarA(pop1_sel5_2)
  
  zmet2_sel6_gv[i,] <- gv(pop1_sel5_2)
}

#recq4
recq4_sel2_gv <- matrix(data = NA, nrow = 100, ncol = 10)
recq4_sel3_gv <- matrix(data = NA, nrow = 100, ncol = 10)
recq4_sel4_gv <- matrix(data = NA, nrow = 100, ncol = 10)
recq4_sel5_gv <- matrix(data = NA, nrow = 100, ncol = 10)
recq4_sel6_gv <- matrix(data = NA, nrow = 100, ncol = 10)
recq4_sel7_gv <- matrix(data = NA, nrow = 100, ncol = 10)
recq4var1 <- matrix(data = NA, ncol = 2, nrow = 200)
recq4var2 <- matrix(data = NA, ncol = 2, nrow = 200)
recq4var3 <- matrix(data = NA, ncol = 2, nrow = 200)
recq4var4 <- matrix(data = NA, ncol = 2, nrow = 200)
recq4var5 <- matrix(data = NA, ncol = 2, nrow = 200)
recq4var6 <- matrix(data = NA, ncol = 2, nrow = 200)
recq4var7 <- matrix(data = NA, ncol = 2, nrow = 200)
for(i in 1:100){
  SP$switchGenMap(recq4_map, centromere = recq4_centromere)
  pop <- randCross2(goodpop[best_good,], badpop[best_bad,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop <- setPheno(pop = pop, h2 = 0.8, simParam = SP)
  recq4var1[i] <- genicVarA(pop)
  
  pop1_sel <- selectInd(pop, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_2 <- selectInd(pop1_sel, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_cross <- randCross2(pop1_sel1_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel1_cross <- setPheno(pop1_sel1_cross, h2 = 0.8, simParam = SP)
  recq4var2[i] <- genicVarA(pop1_sel1_2)
  
  recq4_sel2_gv[i,] <- gv(pop1_sel1_2)
  pop1_sel2 <- selectInd(pop1_sel1_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_2 <- selectInd(pop1_sel2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_cross <- randCross2(pop1_sel2_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel2_cross <- setPheno(pop1_sel2_cross, h2 = 0.8, simParam = SP)
  recq4var3[i] <- genicVarA(pop1_sel2_2)
  
  recq4_sel3_gv[i,] <- gv(pop1_sel2_2)
  pop1_sel3 <- selectInd(pop1_sel2_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_2 <- selectInd(pop1_sel3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_cross <- randCross2(pop1_sel3_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel3_cross <- setPheno(pop1_sel3_cross, h2 = 0.8, simParam = SP)
  recq4var4[i] <- genicVarA(pop1_sel3_2)
  
  recq4_sel4_gv[i,] <- gv(pop1_sel3_2)
  pop1_sel4 <- selectInd(pop1_sel3_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_2 <- selectInd(pop1_sel4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_cross <- randCross2(pop1_sel4_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel4_cross <- setPheno(pop1_sel4_cross, h2 = 0.8, simParam = SP)
  recq4var5[i] <- genicVarA(pop1_sel4_2)
  
  recq4_sel5_gv[i,] <- gv(pop1_sel4_2)
  pop1_sel5 <- selectInd(pop1_sel4_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_2 <- selectInd(pop1_sel5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_cross <- randCross2(pop1_sel5_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel5_cross <- setPheno(pop1_sel5_cross, h2 = 0.8, simParam = SP)
  recq4var6[i] <- genicVarA(pop1_sel5_2)
  
  recq4_sel6_gv[i,] <- gv(pop1_sel5_2)
}

#fancm
fancm_sel2_gv <- matrix(data = NA, nrow = 100, ncol = 10)
fancm_sel3_gv <- matrix(data = NA, nrow = 100, ncol = 10)
fancm_sel4_gv <- matrix(data = NA, nrow = 100, ncol = 10)
fancm_sel5_gv <- matrix(data = NA, nrow = 100, ncol = 10)
fancm_sel6_gv <- matrix(data = NA, nrow = 100, ncol = 10)
fancm_sel7_gv <- matrix(data = NA, nrow = 100, ncol = 10)
fancmvar1 <- matrix(data = NA, ncol = 2, nrow = 200)
fancmvar2 <- matrix(data = NA, ncol = 2, nrow = 200)
fancmvar3 <- matrix(data = NA, ncol = 2, nrow = 200)
fancmvar4 <- matrix(data = NA, ncol = 2, nrow = 200)
fancmvar5 <- matrix(data = NA, ncol = 2, nrow = 200)
fancmvar6 <- matrix(data = NA, ncol = 2, nrow = 200)
fancmvar7 <- matrix(data = NA, ncol = 2, nrow = 200)
for(i in 1:100){
  SP$switchGenMap(fancm_map, centromere = fancm_centromere)
  pop <- randCross2(goodpop[best_good,], badpop[best_bad,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop <- setPheno(pop = pop, h2 = 0.8, simParam = SP)
  fancmvar1[i] <- genicVarA(pop)
  
  pop1_sel <- selectInd(pop, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_2 <- selectInd(pop1_sel, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_cross <- randCross2(pop1_sel1_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel1_cross <- setPheno(pop1_sel1_cross, h2 = 0.8, simParam = SP)
  fancmvar2[i] <- genicVarA(pop1_sel1_2)
  
  fancm_sel2_gv[i,] <- gv(pop1_sel1_2)
  pop1_sel2 <- selectInd(pop1_sel1_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_2 <- selectInd(pop1_sel2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_cross <- randCross2(pop1_sel2_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel2_cross <- setPheno(pop1_sel2_cross, h2 = 0.8, simParam = SP)
  fancmvar3[i] <- genicVarA(pop1_sel2_2)
  
  fancm_sel3_gv[i,] <- gv(pop1_sel2_2)
  pop1_sel3 <- selectInd(pop1_sel2_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_2 <- selectInd(pop1_sel3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_cross <- randCross2(pop1_sel3_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel3_cross <- setPheno(pop1_sel3_cross, h2 = 0.8, simParam = SP)
  fancmvar4[i] <- genicVarA(pop1_sel3_2)
  
  fancm_sel4_gv[i,] <- gv(pop1_sel3_2)
  pop1_sel4 <- selectInd(pop1_sel3_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_2 <- selectInd(pop1_sel4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_cross <- randCross2(pop1_sel4_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel4_cross <- setPheno(pop1_sel4_cross, h2 = 0.8, simParam = SP)
  fancmvar5[i] <- genicVarA(pop1_sel4_2)
  
  fancm_sel5_gv[i,] <- gv(pop1_sel4_2)
  pop1_sel5 <- selectInd(pop1_sel4_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_2 <- selectInd(pop1_sel5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_cross <- randCross2(pop1_sel5_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel5_cross <- setPheno(pop1_sel5_cross, h2 = 0.8, simParam = SP)
  fancmvar6[i] <- genicVarA(pop1_sel5_2)
  
  fancm_sel6_gv[i,] <- gv(pop1_sel5_2)
}

#ideal1 = 10X global increase
ideal1_sel2_gv <- matrix(data = NA, nrow = 100, ncol = 10)
ideal1_sel3_gv <- matrix(data = NA, nrow = 100, ncol = 10)
ideal1_sel4_gv <- matrix(data = NA, nrow = 100, ncol = 10)
ideal1_sel5_gv <- matrix(data = NA, nrow = 100, ncol = 10)
ideal1_sel6_gv <- matrix(data = NA, nrow = 100, ncol = 10)
ideal1_sel7_gv <- matrix(data = NA, nrow = 100, ncol = 10)
ideal1var1 <- matrix(data = NA, ncol = 2, nrow = 200)
ideal1var2 <- matrix(data = NA, ncol = 2, nrow = 200)
ideal1var3 <- matrix(data = NA, ncol = 2, nrow = 200)
ideal1var4 <- matrix(data = NA, ncol = 2, nrow = 200)
ideal1var5 <- matrix(data = NA, ncol = 2, nrow = 200)
ideal1var6 <- matrix(data = NA, ncol = 2, nrow = 200)
ideal1var7 <- matrix(data = NA, ncol = 2, nrow = 200)
for(i in 1:100){
  SP$switchGenMap(ideal1_map, centromere = ideal1_centromere)
  pop <- randCross2(goodpop[best_good,], badpop[best_bad,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop <- setPheno(pop = pop, h2 = 0.8, simParam = SP)
  ideal1var1[i] <- genicVarA(pop)
  
  pop1_sel <- selectInd(pop, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_2 <- selectInd(pop1_sel, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_cross <- randCross2(pop1_sel1_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel1_cross <- setPheno(pop1_sel1_cross, h2 = 0.8, simParam = SP)
  ideal1var2[i] <- genicVarA(pop1_sel1_2)
  
  ideal1_sel2_gv[i,] <- gv(pop1_sel1_2)
  pop1_sel2 <- selectInd(pop1_sel1_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_2 <- selectInd(pop1_sel2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_cross <- randCross2(pop1_sel2_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel2_cross <- setPheno(pop1_sel2_cross, h2 = 0.8, simParam = SP)
  ideal1var3[i] <- genicVarA(pop1_sel2_2)
  
  ideal1_sel3_gv[i,] <- gv(pop1_sel2_2)
  pop1_sel3 <- selectInd(pop1_sel2_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_2 <- selectInd(pop1_sel3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_cross <- randCross2(pop1_sel3_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel3_cross <- setPheno(pop1_sel3_cross, h2 = 0.8, simParam = SP)
  ideal1var4[i] <- genicVarA(pop1_sel3_2)
  
  ideal1_sel4_gv[i,] <- gv(pop1_sel3_2)
  pop1_sel4 <- selectInd(pop1_sel3_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_2 <- selectInd(pop1_sel4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_cross <- randCross2(pop1_sel4_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel4_cross <- setPheno(pop1_sel4_cross, h2 = 0.8, simParam = SP)
  ideal1var5[i] <- genicVarA(pop1_sel4_2)
  
  ideal1_sel5_gv[i,] <- gv(pop1_sel4_2)
  pop1_sel5 <- selectInd(pop1_sel4_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_2 <- selectInd(pop1_sel5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_cross <- randCross2(pop1_sel5_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel5_cross <- setPheno(pop1_sel5_cross, h2 = 0.8, simParam = SP)
  ideal1var6[i] <- genicVarA(pop1_sel5_2)
  
  ideal1_sel6_gv[i,] <- gv(pop1_sel5_2)
}

#ideal2 = ddm1/zmet2 double mutant
ideal2_sel2_gv <- matrix(data = NA, nrow = 100, ncol = 10)
ideal2_sel3_gv <- matrix(data = NA, nrow = 100, ncol = 10)
ideal2_sel4_gv <- matrix(data = NA, nrow = 100, ncol = 10)
ideal2_sel5_gv <- matrix(data = NA, nrow = 100, ncol = 10)
ideal2_sel6_gv <- matrix(data = NA, nrow = 100, ncol = 10)
ideal2_sel7_gv <- matrix(data = NA, nrow = 100, ncol = 10)
ideal2var1 <- matrix(data = NA, ncol = 2, nrow = 200)
ideal2var2 <- matrix(data = NA, ncol = 2, nrow = 200)
ideal2var3 <- matrix(data = NA, ncol = 2, nrow = 200)
ideal2var4 <- matrix(data = NA, ncol = 2, nrow = 200)
ideal2var5 <- matrix(data = NA, ncol = 2, nrow = 200)
ideal2var6 <- matrix(data = NA, ncol = 2, nrow = 200)
ideal2var7 <- matrix(data = NA, ncol = 2, nrow = 200)
for(i in 1:100){
  SP$switchGenMap(ideal2_map, centromere = ideal2_centromere)
  pop <- randCross2(goodpop[best_good,], badpop[best_bad,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop <- setPheno(pop = pop, h2 = 0.8, simParam = SP)
  ideal2var1[i] <- genicVarA(pop)
  
  pop1_sel <- selectInd(pop, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_2 <- selectInd(pop1_sel, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_cross <- randCross2(pop1_sel1_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel1_cross <- setPheno(pop1_sel1_cross, h2 = 0.8, simParam = SP)
  ideal2var2[i] <- genicVarA(pop1_sel1_2)
  
  ideal2_sel2_gv[i,] <- gv(pop1_sel1_2)
  pop1_sel2 <- selectInd(pop1_sel1_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_2 <- selectInd(pop1_sel2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_cross <- randCross2(pop1_sel2_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel2_cross <- setPheno(pop1_sel2_cross, h2 = 0.8, simParam = SP)
  ideal2var3[i] <- genicVarA(pop1_sel2_2)
  
  ideal2_sel3_gv[i,] <- gv(pop1_sel2_2)
  pop1_sel3 <- selectInd(pop1_sel2_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_2 <- selectInd(pop1_sel3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_cross <- randCross2(pop1_sel3_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel3_cross <- setPheno(pop1_sel3_cross, h2 = 0.8, simParam = SP)
  ideal2var4[i] <- genicVarA(pop1_sel3_2)
  
  ideal2_sel4_gv[i,] <- gv(pop1_sel3_2)
  pop1_sel4 <- selectInd(pop1_sel3_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_2 <- selectInd(pop1_sel4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_cross <- randCross2(pop1_sel4_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel4_cross <- setPheno(pop1_sel4_cross, h2 = 0.8, simParam = SP)
  ideal2var5[i] <- genicVarA(pop1_sel4_2)
  
  ideal2_sel5_gv[i,] <- gv(pop1_sel4_2)
  pop1_sel5 <- selectInd(pop1_sel4_cross, nInd = 10, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_2 <- selectInd(pop1_sel5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_cross <- randCross2(pop1_sel5_2, goodpop[best_good,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel5_cross <- setPheno(pop1_sel5_cross, h2 = 0.8, simParam = SP)
  ideal2var6[i] <- genicVarA(pop1_sel5_2)
  
  ideal2_sel6_gv[i,] <- gv(pop1_sel5_2)
}

#introgression scheme analysis

#putting all genetic values from all generations into one data frame
wt <- c(pop1_sel2_gv[,1:5], pop1_sel3_gv[,1:5], pop1_sel4_gv[,1:5], pop1_sel5_gv[,1:5], pop1_sel6_gv[,1:5])
wt_gain <- c(pop1_sel3_gv[,1:5]-pop1_sel2_gv[,1:5], pop1_sel4_gv[,1:5]-pop1_sel2_gv[,1:5], pop1_sel5_gv[,1:5]-pop1_sel2_gv[,1:5], pop1_sel6_gv[,1:5]-pop1_sel2_gv[,1:5])

ddm1 <- c(ddm1_sel2_gv[,1:5], ddm1_sel3_gv[,1:5], ddm1_sel4_gv[,1:5], ddm1_sel5_gv[,1:5], ddm1_sel6_gv[,1:5])
ddm1_gain <- c(ddm1_sel3_gv[,1:5] - ddm1_sel2_gv[,1:5], ddm1_sel4_gv[,1:5]- ddm1_sel2_gv[,1:5], ddm1_sel5_gv[,1:5] - ddm1_sel2_gv[,1:5], ddm1_sel6_gv[,1:5]-ddm1_sel2_gv[,1:5])

zmet2 <- c(zmet2_sel2_gv[,1:5], zmet2_sel3_gv[,1:5], zmet2_sel4_gv[,1:5], zmet2_sel5_gv[,1:5], zmet2_sel6_gv[,1:5])
zmet2_gain <- c(zmet2_sel3_gv[,1:5]-zmet2_sel2_gv[,1:5], zmet2_sel4_gv[,1:5]-zmet2_sel2_gv[,1:5], zmet2_sel5_gv[,1:5]-zmet2_sel2_gv[,1:5], zmet2_sel6_gv[,1:5]-zmet2_sel2_gv[,1:5])

ideal1 <- c(ideal1_sel2_gv[,1:5], ideal1_sel3_gv[,1:5], ideal1_sel4_gv[,1:5], ideal1_sel5_gv[,1:5], ideal1_sel6_gv[,1:5])
ideal1_gain <- c(ideal1_sel3_gv[,1:5]-ideal1_sel2_gv[,1:5], ideal1_sel4_gv[,1:5]-ideal1_sel2_gv[,1:5], ideal1_sel5_gv[,1:5]-ideal1_sel2_gv[,1:5], ideal1_sel6_gv[,1:5]-ideal1_sel2_gv[,1:5])

ideal2 <- c(ideal2_sel2_gv[,1:5], ideal2_sel3_gv[,1:5], ideal2_sel4_gv[,1:5], ideal2_sel5_gv[,1:5], ideal2_sel6_gv[,1:5])
ideal2_gain <- c(ideal2_sel3_gv[,1:5]-ideal2_sel2_gv[,1:5], ideal2_sel4_gv[,1:5]-ideal2_sel2_gv[,1:5], ideal2_sel5_gv[,1:5]-ideal2_sel2_gv[,1:5], ideal2_sel6_gv[,1:5]-ideal2_sel2_gv[,1:5])

fancm <- c(fancm_sel2_gv[,1:5], fancm_sel3_gv[,1:5], fancm_sel4_gv[,1:5], fancm_sel5_gv[,1:5], fancm_sel6_gv[,1:5])
fancm_gain <- c(fancm_sel3_gv[,1:5]-fancm_sel2_gv[,1:5], fancm_sel4_gv[,1:5]-fancm_sel2_gv[,1:5], fancm_sel5_gv[,1:5]-fancm_sel2_gv[,1:5], fancm_sel6_gv[,1:5]-fancm_sel2_gv[,1:5])

recq4_gain <- c(recq4_sel3_gv[,1:5]-recq4_sel2_gv[,1:5], recq4_sel4_gv[,1:5]-recq4_sel2_gv[,1:5], recq4_sel5_gv[,1:5]-recq4_sel2_gv[,1:5], recq4_sel6_gv[,1:5]-recq4_sel2_gv[,1:5])
recq4 <- c(recq4_sel2_gv[,1:5], recq4_sel3_gv[,1:5], recq4_sel4_gv[,1:5], recq4_sel5_gv[,1:5], recq4_sel6_gv[,1:5])

#for resistance loci
wt2 <- c(mean(pop1_sel2_gv[,6:10]), mean(pop1_sel3_gv[,6:10]), mean(pop1_sel4_gv[,6:10]), mean(pop1_sel5_gv[,6:10]), mean(pop1_sel6_gv[,6:10]))

ddm12 <- c(mean(ddm1_sel2_gv[,6:10]), mean(ddm1_sel3_gv[,6:10]), mean(ddm1_sel4_gv[,6:10]), mean(ddm1_sel5_gv[,6:10]), mean(ddm1_sel6_gv[,6:10]))

zmet22 <- c(mean(zmet2_sel2_gv[,6:10]), mean(zmet2_sel3_gv[,6:10]), mean(zmet2_sel4_gv[,6:10]), mean(zmet2_sel5_gv[,6:10]), mean(zmet2_sel6_gv[,6:10]))

recq42 <- c(mean(recq4_sel2_gv[,6:10]), mean(recq4_sel3_gv[,6:10]), mean(recq4_sel4_gv[,6:10]), mean(recq4_sel5_gv[,6:10]), mean(recq4_sel6_gv[,6:10]))

ideal12 <- c(mean(ideal1_sel2_gv[,6:10]), mean(ideal1_sel3_gv[,6:10]), mean(ideal1_sel4_gv[,6:10]), mean(ideal1_sel5_gv[,6:10]), mean(ideal1_sel6_gv[,6:10]))

ideal22 <- c(mean(ideal2_sel2_gv[,6:10]), mean(ideal2_sel3_gv[,6:10]), mean(ideal2_sel4_gv[,6:10]), mean(ideal2_sel5_gv[,6:10]), mean(ideal2_sel6_gv[,6:10]))

fancm2 <- c(mean(fancm_sel2_gv[,6:10]), mean(fancm_sel3_gv[,6:10]), mean(fancm_sel4_gv[,6:10]), mean(fancm_sel5_gv[,6:10]), mean(fancm_sel6_gv[,6:10]))


#combining all data together
all <- cbind(wt, ddm1, zmet2, recq4, fancm, ideal1, ideal2)
all <- as.data.frame(all)
all2 <- as.data.frame(matrix(data = NA, nrow = 17500))
all2$gv <- c(all$wt, all$ddm1, all$zmet2, all$recq4, all$fancm, all$ideal1, all$ideal2)
all2$generation <- rep(1:5, size = 7, each = 500)
all2$gen_map <- rep(c("wt","ddm1", "zmet2", "recq4", "fancm", "ideal1", "ideal2"), size = 7, each = 2500)

#looking at introgression of resistance loci
trait2 <- cbind(wt2, ddm12)
trait2 <- as.data.frame(trait2)
trait22 <- as.data.frame(matrix(data = NA, nrow = 35))
trait22$gv <- c(trait2$wt2, trait2$ddm12, trait2$zmet22, trait2$recq42, trait2$ideal12, trait2$ideal22, trait2$fancm2)
trait22$generation <- rep(1:5, size = 7, each = 1)
trait22$gen_map <- rep(c("wt","ddm1", "zmet2", "recq4", "ideal1", "ideal2", "fancm"), size = 7, each = 1)

group.colors = c("ddm1" = "#E69F00", "recq4" = "#56B4E9", "wt" = "#009E73", "zmet2" = "#F0E442",
                 "fancm" = "#0072B2", "ideal1" = "#D55E00", "ideal2" = "#CC79A7")

#ggplot code to look at polygenic trait
ggplot(all2, aes(x=as.factor(generation), y= gv, fill=gen_map)) + 
  geom_boxplot() + theme_bw() + xlab("Generations") + ylab("Genetic Values") + ggtitle("Genetic Values after 5 gens of backcrossing: 3 unlinked QTLs") +
  scale_fill_manual(values=group.colors, name = "Genetic Map Used", labels = c("ddm1", "recq4", "wt", "zmet2",
                                                                               "fancm", "10X", "ddm1/zmet2")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12), legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=12), plot.title = element_text(size=14), legend.position = c(0.8,0.2), legend.key.size = unit(0.5, "lines")) +
  geom_hline(yintercept=7.511954, linetype="dashed", color = "blue") + annotate('text', x = 1.5, y = 7.6, label = 'Recurrent Parent Level', color='blue', size = 4) 


#ggplot code for resistance loci
ggplot(trait22, aes(x=as.factor(generation), y= gv, color=gen_map)) + 
  geom_point() + theme_bw() + xlab("Generation") + ylab("Mean Genetic Values") + ggtitle("Genetic Values of Unlinked QTLs in arms over 5 gens of backcrossing") + 
  scale_color_manual(values=group.colors, name = "Genetic Map Used", labels = c("ddm1", "recq4", "wt", "zmet2",
                                                                               "fancm", "10X", "ddm1/zmet2")) +
  geom_jitter(width = 0.2, height = 0.01) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12), legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=12), plot.title = element_text(size=14), legend.position = c(0.8,0.2), legend.key.size = unit(0.5, "lines")) +
  geom_hline(yintercept=0.7, linetype="dashed", color = "blue") + annotate('text', x = 1.5, y = 0.71, label = 'Recurrent Parent Level', color='blue', size = 4) + ylim(0.5,0.8)

#binding all the variance matrices together
wtvarall <- c(wtvar1, wtvar2, wtvar3, wtvar4, wtvar5,
              wtvar6)

ddm1varall <- c(ddm1var1, ddm1var2, ddm1var3, ddm1var4, ddm1var5,
                ddm1var6)

zmet2varall <- c(zmet2var1, zmet2var2, zmet2var3, zmet2var4, zmet2var5,
                 zmet2var6)

recq4varall <- c(recq4var1, recq4var2, recq4var3, recq4var4, recq4var5,
                 recq4var6)

ideal1varall <- c(ideal1var1, ideal1var2, ideal1var3, ideal1var4, ideal1var5,
                  ideal1var6)

ideal2varall <- c(ideal2var1, ideal2var2, ideal2var3, ideal2var4, ideal2var5,
                  ideal2var6)

ddm1_zmet2varall <- c(ideal2var1, ideal2var2, ideal2var3, ideal2var4, ideal2var5,
                      ideal2var6)

fancmvarall <- c(fancmvar1, fancmvar2, fancmvar3, fancmvar4, fancmvar5,
                 fancmvar6)

#putting all variance data into one data frame
allvar <- cbind(wtvarall, ddm1varall, zmet2varall, recq4varall, ideal1varall, ideal2varall, ddm1_zmet2varall, fancmvarall)
allvar <- as.data.frame(allvar)
allvar2 <- as.data.frame(matrix(data = NA, nrow = 16800))
allvar2$var <- c(allvar$wtvarall, allvar$ddm1varall, allvar$zmet2varall, allvar$recq4varall, allvar$ideal1varall, allvar$ideal2varall, allvar$fancmvarall)
allvar2$gen <- rep(1:5, size = 7, each = 480)
allvar2$gen_map <- rep(c("wt","ddm1", "zmet2", "recq4", "ideal1", "ideal2", "fancm"), size = 7, each = 2400)

#ggplot code to plot the variance
ggplot(allvar2, aes(x=as.factor(gen), y=var, fill=gen_map)) +  ggtitle("Additive genetic variance over 5 gens of backcrossing: 3 centromere QTLs") +
  geom_boxplot() + theme_bw() + xlab("Generation") + ylab("Additive Genetic Variances") + scale_fill_manual(values=group.colors, name = "Genetic Map Used", labels = c("ddm1", "recq4", "wt", "zmet2",
                                                                                                                                                                       "fancm", "10X", "ddm1/zmet2")) + ggtitle("Additive Genetic Variance after 5 gens of backcrossing") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12), legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=12), plot.title = element_text(size=14), legend.position = c(0.8,0.8), legend.key.size = unit(0.5, "lines"))
