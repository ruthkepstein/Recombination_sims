library(AlphaSimR)
final_map<-readRDS("japonica_final_map.RData")
real_centromere <-readRDS("japonica_centromeres.RData")
#generating repulsion ONLY linkages between QTL
addEff_repulsion <- runif(209, min = 0, max = 0.1)
addEff_repulsion <- addEff_repulsion * c(-1,1)
sum(addEff_repulsion)

#generating mix of repulsion & coupling linkages between QTL
addEff_mix <- runif(209, min = -0.1, max = 0.1)
sum(addEff_mix)

#Setting up founder population
founderPop <- quickHaplo(nInd = 200, nChr = 10, inbred = TRUE, ploidy = 2L, segSites = 10)
founderPop@genMap <- final_map
founderPop@centromere <- real_centromere
SP = SimParam$new(founderPop)
SP$setTrackRec(TRUE)
SP$p = 0.15
trait_yield <- new("TraitA", nLoci = 200L, lociPerChr= c(35L, 24L, 24L, 21L, 22L, 16L, 18L, 18L, 17L, 14L),
                   lociLoc = c(4L,5L,6L,7L,8L,9L,10L,11L,12L,13L,14L,15L,16L,17L,18L,19L,20L,21L,22L,23L,42L,43L,45L,46L,47L,
                               48L,49L,50L,51L,52L,53L,54L,55L,56L,58L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,13L,14L,20L,30L,
                               35L,36L,37L,38L,39L,40L,41L,42L,43L,44L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,16L,25L,31L,
                               34L,35L,36L,37L,38L,39L,40L,41L,42L,43L,44L, 3L,4L,5L,9L,10L,11L,12L,17L,22L,31L,
                               36L,37L,38L,39L,40L,41L,42L,43L,44L,45L,46L, 3L,4L,5L,6L,7L,8L,9L,10L,11L,13L,25L,
                               32L,33L,34L,35L,36L,37L,38L,39L,40L,41L,42L, 3L,4L,5L,6L,7L,8L,13L,20L,
                               25L,26L,27L,28L,29L,30L,31L,32L, 2L,3L,4L,5L,6L,7L,8L,9L,17L,26L,27L,28L,29L,30L,31L,32L,33L,
                               34L, 2L,3L,4L,5L,6L,7L,8L,9L, 16L, 26L,27L,28L,29L,30L,31L,32L,33L,
                               34L, 2L,3L,4L,5L,6L,7L,8L,9L, 14L, 23L,24L,25L,26L,27L,28L,29L,30L, 
                               1L,2L,3L,4L,5L,6L,10L, 22L,23L,24L,25L,26L,27L,28L), addEff = addEff_mix, intercept = 0.1)
SP$manAddTrait(trait_yield)

#Diverging founder population into 2 populations; "good pop" & "bad pop"
#Good pop = selecting FOR polygenic trait; no "resistance" QTLs
pop_good_sel10 <- vector(mode = "list", length = 100)
for(i in 1:100){
  pop_good <- newPop(founderPop, simParam = SP)
  pop_good <- randCross(pop_good, nCrosses= 10, nProgeny=10, simParam = SP)
  pop_good <- setPheno(pop_good, h2 = 0.8, simParam = SP)
  
  pop_good_sel0 <- selectInd(pop_good, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good1 <- self(pop_good, nProgeny=20, keepParents= FALSE, simParam = SP)
  pop_good1 <- setPheno(pop_good1, h2 = 0.8, simParam = SP)
  
  pop_good_sel <- selectInd(pop_good1, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good2 <- self(pop_good_sel,  nProgeny=20, keepParents= FALSE, simParam = SP)
  pop_good2 <- setPheno(pop_good2, h2 = 0.8, simParam = SP)
  
  pop_good_sel2 <- selectInd(pop_good2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good3 <- self(pop_good_sel2,  nProgeny=20, keepParents= FALSE, simParam = SP)
  pop_good3 <- setPheno(pop_good3, h2 = 0.8, simParam = SP)
  
  pop_good_sel3 <- selectInd(pop_good3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good4 <- self(pop_good_sel2,nProgeny=20, keepParents= FALSE, simParam = SP)
  pop_good4 <- setPheno(pop_good4, h2 = 0.8, simParam = SP)
  
  pop_good_sel4 <- selectInd(pop_good4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good5 <- self(pop_good_sel4, nProgeny=20, keepParents= FALSE, simParam = SP)
  pop_good5 <- setPheno(pop_good5, h2 = 0.8, simParam = SP)
  
  pop_good_sel5 <- selectInd(pop_good5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good6 <- self(pop_good_sel5, nProgeny=20, keepParents= FALSE, simParam = SP)
  pop_good6 <- setPheno(pop_good6, h2 = 0.8, simParam = SP)
  
  pop_good_sel6 <- selectInd(pop_good6, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good7 <- self(pop_good_sel6, nProgeny=20, keepParents= FALSE, simParam = SP)
  pop_good7 <- setPheno(pop_good7, h2 = 0.8, simParam = SP)
  
  pop_good_sel7 <- selectInd(pop_good7, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good8 <- self(pop_good_sel7,  nProgeny=20, keepParents= FALSE, simParam = SP)
  pop_good8 <- setPheno(pop_good8, h2 = 0.8, simParam = SP)
  
  pop_good_sel8 <- selectInd(pop_good8, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good9 <- self(pop_good_sel8,  nProgeny=20, keepParents= FALSE, simParam = SP)
  pop_good9 <- setPheno(pop_good9, h2 = 0.8, simParam = SP)
  
  pop_good_sel9 <- selectInd(pop_good9, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good10 <- self(pop_good_sel9, nProgeny=20, keepParents= FALSE, simParam = SP)
  pop_good10 <- setPheno(pop_good10, h2 = 0.8, simParam = SP)
  
  pop_good_sel10[[i]] <- selectInd(pop_good10, nInd = 2, use = "gv", trait = 1, selectop = TRUE, returnPop = TRUE, simParam = SP)
}

#Bad pop = not selecting for polygenic trait even though its present; "resistance" QTLs present & selected FOR
pop_bad_sel10 <- vector(mode = "list", length = 100)
trait <- new("TraitAD", nLoci = 3L, lociPerChr = c(0L, 0L, 0L, 3L, 0L, 0L, 0L, 0L, 0L, 0L),
             lociLoc = c(1L,2L,3L), addEff = c(0.1,0.1,0.1), domEff = c(0.2,0.2,0.2), intercept = 0.1)
SP$resetPed()
SP$manAddTrait(trait)
for(i in 1:100){
  pop_bad <- newPop(founderPop, simParam = SP)
  pop_bad <- setPheno(pop_bad, h2 = 0.8, simParam = SP)
  
  pop_bad1 <- self(pop_bad, nProgeny=10, keepParents=FALSE, simParam = SP)
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
  
  pop_bad_sel10[[i]] <- selectInd(pop_bad10, nInd = 2, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
}

goodpop = mergePops(pop_good_sel10)
badpop = mergePops(pop_bad_sel10)

##Backcrossing scheme
#selecting best individual from "good pop" and individual with resistance QTLs in "bad pop"
which(goodpop@gv == max(goodpop@gv), arr.ind = TRUE)
which(badpop@gv == max(badpop@gv), arr.ind = TRUE)

#wild-type first
pop1_sel2_gv <- matrix(data = NA, nrow = 100, ncol = 20)
pop1_sel3_gv <- matrix(data = NA, nrow = 100, ncol = 20)
pop1_sel4_gv <- matrix(data = NA, nrow = 100, ncol = 20)
pop1_sel5_gv <- matrix(data = NA, nrow = 100, ncol = 20)
pop1_sel6_gv <- matrix(data = NA, nrow = 100, ncol = 20)
pop1_sel7_gv <- matrix(data = NA, nrow = 100, ncol = 20)
wtvar1 <- matrix(data = NA, ncol = 2, nrow = 200)
wtvar2 <- matrix(data = NA, ncol = 2, nrow = 200)
wtvar3 <- matrix(data = NA, ncol = 2, nrow = 200)
wtvar4 <- matrix(data = NA, ncol = 2, nrow = 200)
wtvar5 <- matrix(data = NA, ncol = 2, nrow = 200)
wtvar6 <- matrix(data = NA, ncol = 2, nrow = 200)
wtvar7 <- matrix(data = NA, ncol = 2, nrow = 200)
for(i in 1:100){
  SP$switchGenMap(final_map, centromere = real_centromere)
  pop <- randCross2(badpop[110,], goodpop[236,], nCrosses = 200, nProgeny = 1, simParam = SP)
  pop <- setPheno(pop = pop, h2 = 0.8, simParam = SP)
  wtvar1[i] <- varA(pop)
  
  pop1_sel <- selectInd(pop, nInd = 100, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_2 <- selectInd(pop1_sel, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_cross <- self2(pop1_sel1_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel1_cross <- setPheno(pop1_sel1_cross, h2 = 0.8, simParam = SP)
  wtvar2[i] <- varA(pop1_sel1_2)
  
  pop1_sel2_gv[i,] <- gv(pop1_sel1_2)
  pop1_sel2 <- selectInd(pop1_sel1_cross, nInd = 20, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_2 <- selectInd(pop1_sel2, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_cross <- self2(pop1_sel2_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel2_cross <- setPheno(pop1_sel2_cross, h2 = 0.8, simParam = SP)
  wtvar3[i] <- varA(pop1_sel2_2)
  
  pop1_sel3_gv[i,] <- gv(pop1_sel2_2)
  pop1_sel3 <- selectInd(pop1_sel2_cross, nInd = 20, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_2 <- selectInd(pop1_sel3, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_cross <- self2(pop1_sel3_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel3_cross <- setPheno(pop1_sel3_cross, h2 = 0.8, simParam = SP)
  wtvar4[i] <- varA(pop1_sel3_2)
  
  pop1_sel4_gv[i,] <- gv(pop1_sel3_2)
  pop1_sel4 <- selectInd(pop1_sel3_cross, nInd = 20, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_2 <- selectInd(pop1_sel4, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_cross <- self2(pop1_sel4_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel4_cross <- setPheno(pop1_sel4_cross, h2 = 0.8, simParam = SP)
  wtvar5[i] <- varA(pop1_sel4_2)
  
  pop1_sel5_gv[i,] <- gv(pop1_sel4_2)
  pop1_sel5 <- selectInd(pop1_sel4_cross, nInd = 20, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_2 <- selectInd(pop1_sel5, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_cross <- self2(pop1_sel5_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel5_cross <- setPheno(pop1_sel5_cross, h2 = 0.8, simParam = SP)
  wtvar6[i] <- varA(pop1_sel5_2)
  
  pop1_sel6_gv[i,] <- gv(pop1_sel5_2)
}

#ddm1
ddm1_sel2_gv <- matrix(data = NA, nrow = 100, ncol = 20)
ddm1_sel3_gv <- matrix(data = NA, nrow = 100, ncol = 20)
ddm1_sel4_gv <- matrix(data = NA, nrow = 100, ncol = 20)
ddm1_sel5_gv <- matrix(data = NA, nrow = 100, ncol = 20)
ddm1_sel6_gv <- matrix(data = NA, nrow = 100, ncol = 20)
ddm1_sel7_gv <- matrix(data = NA, nrow = 100, ncol = 20)
ddm1var1 <- matrix(data = NA, ncol = 2, nrow = 200)
ddm1var2 <- matrix(data = NA, ncol = 2, nrow = 200)
ddm1var3 <- matrix(data = NA, ncol = 2, nrow = 200)
ddm1var4 <- matrix(data = NA, ncol = 2, nrow = 200)
ddm1var5 <- matrix(data = NA, ncol = 2, nrow = 200)
ddm1var6 <- matrix(data = NA, ncol = 2, nrow = 200)
ddm1var7 <- matrix(data = NA, ncol = 2, nrow = 200)
for(i in 1:100){
  SP$switchGenMap(ddm1_map, centromere = ddm1_centromere)
  pop <- self2(goodpop[236,], badpop[110,], nCrosses = 200, nProgeny = 1, simParam = SP)
  pop <- setPheno(pop = pop, h2 = 0.8, simParam = SP)
  ddm1var1[i] <- varA(pop)
  
  pop1_sel <- selectInd(pop, nInd = 100, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_2 <- selectInd(pop1_sel, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_cross <- self2(pop1_sel1_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel1_cross <- setPheno(pop1_sel1_cross, h2 = 0.8, simParam = SP)
  ddm1var2[i] <- varA(pop1_sel1_2)
  
  ddm1_sel2_gv[i,] <- gv(pop1_sel1_2)
  pop1_sel2 <- selectInd(pop1_sel1_cross, nInd = 20, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_2 <- selectInd(pop1_sel2, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_cross <- self2(pop1_sel2_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel2_cross <- setPheno(pop1_sel2_cross, h2 = 0.8, simParam = SP)
  ddm1var3[i] <- varA(pop1_sel2_2)
  
  ddm1_sel3_gv[i,] <- gv(pop1_sel2_2)
  pop1_sel3 <- selectInd(pop1_sel2_cross, nInd = 20, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_2 <- selectInd(pop1_sel3, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_cross <- self2(pop1_sel3_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel3_cross <- setPheno(pop1_sel3_cross, h2 = 0.8, simParam = SP)
  ddm1var4[i] <- varA(pop1_sel3_2)
  
  ddm1_sel4_gv[i,] <- gv(pop1_sel3_2)
  pop1_sel4 <- selectInd(pop1_sel3_cross, nInd = 20, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_2 <- selectInd(pop1_sel4, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_cross <- self2(pop1_sel4_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel4_cross <- setPheno(pop1_sel4_cross, h2 = 0.8, simParam = SP)
  ddm1var5[i] <- varA(pop1_sel4_2)
  
  ddm1_sel5_gv[i,] <- gv(pop1_sel4_2)
  pop1_sel5 <- selectInd(pop1_sel4_cross, nInd = 20, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_2 <- selectInd(pop1_sel5, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_cross <- self2(pop1_sel5_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel5_cross <- setPheno(pop1_sel5_cross, h2 = 0.8, simParam = SP)
  ddm1var6[i] <- varA(pop1_sel5_2)
  
  ddm1_sel6_gv[i,] <- gv(pop1_sel5_2)
}

#zmet2
zmet2_sel2_gv <- matrix(data = NA, nrow = 100, ncol = 20)
zmet2_sel3_gv <- matrix(data = NA, nrow = 100, ncol = 20)
zmet2_sel4_gv <- matrix(data = NA, nrow = 100, ncol = 20)
zmet2_sel5_gv <- matrix(data = NA, nrow = 100, ncol = 20)
zmet2_sel6_gv <- matrix(data = NA, nrow = 100, ncol = 20)
zmet2_sel7_gv <- matrix(data = NA, nrow = 100, ncol = 20)
zmet2var1 <- matrix(data = NA, ncol = 2, nrow = 200)
zmet2var2 <- matrix(data = NA, ncol = 2, nrow = 200)
zmet2var3 <- matrix(data = NA, ncol = 2, nrow = 200)
zmet2var4 <- matrix(data = NA, ncol = 2, nrow = 200)
zmet2var5 <- matrix(data = NA, ncol = 2, nrow = 200)
zmet2var6 <- matrix(data = NA, ncol = 2, nrow = 200)
zmet2var7 <- matrix(data = NA, ncol = 2, nrow = 200)
for(i in 1:100){
  SP$switchGenMap(zmet2_map, centromere = zmet2_centromere)
  pop <- self2(goodpop[236,], badpop[110,], nCrosses = 200, nProgeny = 1, simParam = SP)
  pop <- setPheno(pop = pop, h2 = 0.8, simParam = SP)
  zmet2var1[i] <- varA(pop)
  
  pop1_sel <- selectInd(pop, nInd = 100, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_2 <- selectInd(pop1_sel, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_cross <- self2(pop1_sel1_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel1_cross <- setPheno(pop1_sel1_cross, h2 = 0.8, simParam = SP)
  zmet2var2[i] <- varA(pop1_sel1_2)
  
  zmet2_sel2_gv[i,] <- gv(pop1_sel1_2)
  pop1_sel2 <- selectInd(pop1_sel1_cross, nInd = 20, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_2 <- selectInd(pop1_sel2, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_cross <- self2(pop1_sel2_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel2_cross <- setPheno(pop1_sel2_cross, h2 = 0.8, simParam = SP)
  zmet2var3[i] <- varA(pop1_sel2_2)
  
  zmet2_sel3_gv[i,] <- gv(pop1_sel2_2)
  pop1_sel3 <- selectInd(pop1_sel2_cross, nInd = 20, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_2 <- selectInd(pop1_sel3, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_cross <- self2(pop1_sel3_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel3_cross <- setPheno(pop1_sel3_cross, h2 = 0.8, simParam = SP)
  zmet2var4[i] <- varA(pop1_sel3_2)
  
  zmet2_sel4_gv[i,] <- gv(pop1_sel3_2)
  pop1_sel4 <- selectInd(pop1_sel3_cross, nInd = 20, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_2 <- selectInd(pop1_sel4, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_cross <- self2(pop1_sel4_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel4_cross <- setPheno(pop1_sel4_cross, h2 = 0.8, simParam = SP)
  zmet2var5[i] <- varA(pop1_sel4_2)
  
  zmet2_sel5_gv[i,] <- gv(pop1_sel4_2)
  pop1_sel5 <- selectInd(pop1_sel4_cross, nInd = 20, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_2 <- selectInd(pop1_sel5, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_cross <- self2(pop1_sel5_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel5_cross <- setPheno(pop1_sel5_cross, h2 = 0.8, simParam = SP)
  zmet2var6[i] <- varA(pop1_sel5_2)
  
  zmet2_sel6_gv[i,] <- gv(pop1_sel5_2)
}

#recq4
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
for(i in 1:100){
  SP$switchGenMap(recq4_map, centromere = recq4_centromere)
  pop <- self2(goodpop[236,], badpop[110,], nCrosses = 200, nProgeny = 1, simParam = SP)
  pop <- setPheno(pop = pop, h2 = 0.8, simParam = SP)
  recq4var1[i] <- varA(pop)
  
  pop1_sel <- selectInd(pop, nInd = 100, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_2 <- selectInd(pop1_sel, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_cross <- self2(pop1_sel1_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel1_cross <- setPheno(pop1_sel1_cross, h2 = 0.8, simParam = SP)
  recq4var2[i] <- varA(pop1_sel1_2)
  
  recq4_sel2_gv[i,] <- gv(pop1_sel1_2)
  pop1_sel2 <- selectInd(pop1_sel1_cross, nInd = 20, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_2 <- selectInd(pop1_sel2, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_cross <- self2(pop1_sel2_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel2_cross <- setPheno(pop1_sel2_cross, h2 = 0.8, simParam = SP)
  recq4var3[i] <- varA(pop1_sel2_2)
  
  recq4_sel3_gv[i,] <- gv(pop1_sel2_2)
  pop1_sel3 <- selectInd(pop1_sel2_cross, nInd = 20, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_2 <- selectInd(pop1_sel3, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_cross <- self2(pop1_sel3_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel3_cross <- setPheno(pop1_sel3_cross, h2 = 0.8, simParam = SP)
  recq4var4[i] <- varA(pop1_sel3_2)
  
  recq4_sel4_gv[i,] <- gv(pop1_sel3_2)
  pop1_sel4 <- selectInd(pop1_sel3_cross, nInd = 20, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_2 <- selectInd(pop1_sel4, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_cross <- self2(pop1_sel4_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel4_cross <- setPheno(pop1_sel4_cross, h2 = 0.8, simParam = SP)
  recq4var5[i] <- varA(pop1_sel4_2)
  
  recq4_sel5_gv[i,] <- gv(pop1_sel4_2)
  pop1_sel5 <- selectInd(pop1_sel4_cross, nInd = 20, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_2 <- selectInd(pop1_sel5, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_cross <- self2(pop1_sel5_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel5_cross <- setPheno(pop1_sel5_cross, h2 = 0.8, simParam = SP)
  recq4var6[i] <- varA(pop1_sel5_2)
  
  recq4_sel6_gv[i,] <- gv(pop1_sel5_2)
}

fancm_sel2_gv <- matrix(data = NA, nrow = 100, ncol = 20)
fancm_sel3_gv <- matrix(data = NA, nrow = 100, ncol = 20)
fancm_sel4_gv <- matrix(data = NA, nrow = 100, ncol = 20)
fancm_sel5_gv <- matrix(data = NA, nrow = 100, ncol = 20)
fancm_sel6_gv <- matrix(data = NA, nrow = 100, ncol = 20)
fancm_sel7_gv <- matrix(data = NA, nrow = 100, ncol = 20)
fancmvar1 <- matrix(data = NA, ncol = 2, nrow = 200)
fancmvar2 <- matrix(data = NA, ncol = 2, nrow = 200)
fancmvar3 <- matrix(data = NA, ncol = 2, nrow = 200)
fancmvar4 <- matrix(data = NA, ncol = 2, nrow = 200)
fancmvar5 <- matrix(data = NA, ncol = 2, nrow = 200)
fancmvar6 <- matrix(data = NA, ncol = 2, nrow = 200)
fancmvar7 <- matrix(data = NA, ncol = 2, nrow = 200)
for(i in 1:100){
  SP$switchGenMap(fancm_map, centromere = fancm_centromere)
  pop <- self2(goodpop[236,], badpop[110,], nCrosses = 200, nProgeny = 1, simParam = SP)
  pop <- setPheno(pop = pop, h2 = 0.8, simParam = SP)
  fancmvar1[i] <- varA(pop)
  
  pop1_sel <- selectInd(pop, nInd = 100, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_2 <- selectInd(pop1_sel, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_cross <- self2(pop1_sel1_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel1_cross <- setPheno(pop1_sel1_cross, h2 = 0.8, simParam = SP)
  fancmvar2[i] <- varA(pop1_sel1_2)
  
  fancm_sel2_gv[i,] <- gv(pop1_sel1_2)
  pop1_sel2 <- selectInd(pop1_sel1_cross, nInd = 20, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_2 <- selectInd(pop1_sel2, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_cross <- self2(pop1_sel2_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel2_cross <- setPheno(pop1_sel2_cross, h2 = 0.8, simParam = SP)
  fancmvar3[i] <- varA(pop1_sel2_2)
  
  fancm_sel3_gv[i,] <- gv(pop1_sel2_2)
  pop1_sel3 <- selectInd(pop1_sel2_cross, nInd = 20, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_2 <- selectInd(pop1_sel3, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_cross <- self2(pop1_sel3_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel3_cross <- setPheno(pop1_sel3_cross, h2 = 0.8, simParam = SP)
  fancmvar4[i] <- varA(pop1_sel3_2)
  
  fancm_sel4_gv[i,] <- gv(pop1_sel3_2)
  pop1_sel4 <- selectInd(pop1_sel3_cross, nInd = 20, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_2 <- selectInd(pop1_sel4, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_cross <- self2(pop1_sel4_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel4_cross <- setPheno(pop1_sel4_cross, h2 = 0.8, simParam = SP)
  fancmvar5[i] <- varA(pop1_sel4_2)
  
  fancm_sel5_gv[i,] <- gv(pop1_sel4_2)
  pop1_sel5 <- selectInd(pop1_sel4_cross, nInd = 20, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_2 <- selectInd(pop1_sel5, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_cross <- self2(pop1_sel5_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel5_cross <- setPheno(pop1_sel5_cross, h2 = 0.8, simParam = SP)
  fancmvar6[i] <- varA(pop1_sel5_2)
  
  fancm_sel6_gv[i,] <- gv(pop1_sel5_2)
}

#ideal1 = 10X global increase
ideal1_sel2_gv <- matrix(data = NA, nrow = 100, ncol = 20)
ideal1_sel3_gv <- matrix(data = NA, nrow = 100, ncol = 20)
ideal1_sel4_gv <- matrix(data = NA, nrow = 100, ncol = 20)
ideal1_sel5_gv <- matrix(data = NA, nrow = 100, ncol = 20)
ideal1_sel6_gv <- matrix(data = NA, nrow = 100, ncol = 20)
ideal1_sel7_gv <- matrix(data = NA, nrow = 100, ncol = 20)
ideal1var1 <- matrix(data = NA, ncol = 2, nrow = 200)
ideal1var2 <- matrix(data = NA, ncol = 2, nrow = 200)
ideal1var3 <- matrix(data = NA, ncol = 2, nrow = 200)
ideal1var4 <- matrix(data = NA, ncol = 2, nrow = 200)
ideal1var5 <- matrix(data = NA, ncol = 2, nrow = 200)
ideal1var6 <- matrix(data = NA, ncol = 2, nrow = 200)
ideal1var7 <- matrix(data = NA, ncol = 2, nrow = 200)
for(i in 1:100){
  SP$switchGenMap(ideal1_map, centromere = ideal1_centromere)
  pop <- self2(goodpop[236,], badpop[110,], nCrosses = 200, nProgeny = 1, simParam = SP)
  pop <- setPheno(pop = pop, h2 = 0.8, simParam = SP)
  ideal1var1[i] <- varA(pop)
  
  pop1_sel <- selectInd(pop, nInd = 100, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_2 <- selectInd(pop1_sel, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_cross <- self2(pop1_sel1_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel1_cross <- setPheno(pop1_sel1_cross, h2 = 0.8, simParam = SP)
  ideal1var2[i] <- varA(pop1_sel1_2)
  
  ideal1_sel2_gv[i,] <- gv(pop1_sel1_2)
  pop1_sel2 <- selectInd(pop1_sel1_cross, nInd = 20, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_2 <- selectInd(pop1_sel2, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_cross <- self2(pop1_sel2_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel2_cross <- setPheno(pop1_sel2_cross, h2 = 0.8, simParam = SP)
  ideal1var3[i] <- varA(pop1_sel2_2)
  
  ideal1_sel3_gv[i,] <- gv(pop1_sel2_2)
  pop1_sel3 <- selectInd(pop1_sel2_cross, nInd = 20, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_2 <- selectInd(pop1_sel3, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_cross <- self2(pop1_sel3_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel3_cross <- setPheno(pop1_sel3_cross, h2 = 0.8, simParam = SP)
  ideal1var4[i] <- varA(pop1_sel3_2)
  
  ideal1_sel4_gv[i,] <- gv(pop1_sel3_2)
  pop1_sel4 <- selectInd(pop1_sel3_cross, nInd = 20, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_2 <- selectInd(pop1_sel4, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_cross <- self2(pop1_sel4_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel4_cross <- setPheno(pop1_sel4_cross, h2 = 0.8, simParam = SP)
  ideal1var5[i] <- varA(pop1_sel4_2)
  
  ideal1_sel5_gv[i,] <- gv(pop1_sel4_2)
  pop1_sel5 <- selectInd(pop1_sel4_cross, nInd = 20, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_2 <- selectInd(pop1_sel5, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_cross <- self2(pop1_sel5_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel5_cross <- setPheno(pop1_sel5_cross, h2 = 0.8, simParam = SP)
  ideal1var6[i] <- varA(pop1_sel5_2)
  
  ideal1_sel6_gv[i,] <- gv(pop1_sel5_2)
}

#ideal2 = ddm1/zmet2 double mutant
ideal2_sel2_gv <- matrix(data = NA, nrow = 100, ncol = 20)
ideal2_sel3_gv <- matrix(data = NA, nrow = 100, ncol = 20)
ideal2_sel4_gv <- matrix(data = NA, nrow = 100, ncol = 20)
ideal2_sel5_gv <- matrix(data = NA, nrow = 100, ncol = 20)
ideal2_sel6_gv <- matrix(data = NA, nrow = 100, ncol = 20)
ideal2_sel7_gv <- matrix(data = NA, nrow = 100, ncol = 20)
ideal2var1 <- matrix(data = NA, ncol = 2, nrow = 200)
ideal2var2 <- matrix(data = NA, ncol = 2, nrow = 200)
ideal2var3 <- matrix(data = NA, ncol = 2, nrow = 200)
ideal2var4 <- matrix(data = NA, ncol = 2, nrow = 200)
ideal2var5 <- matrix(data = NA, ncol = 2, nrow = 200)
ideal2var6 <- matrix(data = NA, ncol = 2, nrow = 200)
ideal2var7 <- matrix(data = NA, ncol = 2, nrow = 200)
for(i in 1:100){
  SP$switchGenMap(ideal2_map, centromere = ideal2_centromere)
  pop <- self2(goodpop[236,], badpop[110,], nCrosses = 200, nProgeny = 1, simParam = SP)
  pop <- setPheno(pop = pop, h2 = 0.8, simParam = SP)
  ideal2var1[i] <- varA(pop)
  
  pop1_sel <- selectInd(pop, nInd = 100, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_2 <- selectInd(pop1_sel, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel1_cross <- self2(pop1_sel1_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel1_cross <- setPheno(pop1_sel1_cross, h2 = 0.8, simParam = SP)
  ideal2var2[i] <- varA(pop1_sel1_2)
  
  ideal2_sel2_gv[i,] <- gv(pop1_sel1_2)
  pop1_sel2 <- selectInd(pop1_sel1_cross, nInd = 20, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_2 <- selectInd(pop1_sel2, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel2_cross <- self2(pop1_sel2_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel2_cross <- setPheno(pop1_sel2_cross, h2 = 0.8, simParam = SP)
  ideal2var3[i] <- varA(pop1_sel2_2)
  
  ideal2_sel3_gv[i,] <- gv(pop1_sel2_2)
  pop1_sel3 <- selectInd(pop1_sel2_cross, nInd = 20, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_2 <- selectInd(pop1_sel3, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel3_cross <- self2(pop1_sel3_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel3_cross <- setPheno(pop1_sel3_cross, h2 = 0.8, simParam = SP)
  ideal2var4[i] <- varA(pop1_sel3_2)
  
  ideal2_sel4_gv[i,] <- gv(pop1_sel3_2)
  pop1_sel4 <- selectInd(pop1_sel3_cross, nInd = 20, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_2 <- selectInd(pop1_sel4, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel4_cross <- self2(pop1_sel4_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel4_cross <- setPheno(pop1_sel4_cross, h2 = 0.8, simParam = SP)
  ideal2var5[i] <- varA(pop1_sel4_2)
  
  ideal2_sel5_gv[i,] <- gv(pop1_sel4_2)
  pop1_sel5 <- selectInd(pop1_sel4_cross, nInd = 20, use = "gv", trait = 2, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_2 <- selectInd(pop1_sel5, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop1_sel5_cross <- self2(pop1_sel5_2, goodpop[236,], nCrosses = 10, nProgeny = 10, simParam = SP)
  pop1_sel5_cross <- setPheno(pop1_sel5_cross, h2 = 0.8, simParam = SP)
  ideal2var6[i] <- varA(pop1_sel5_2)
  
  ideal2_sel6_gv[i,] <- gv(pop1_sel5_2)
}

#introgression scheme analysis

#putting all genetic values from all generations into one data frame
wt <- c(pop1_sel2_gv[,1:10], pop1_sel3_gv[,1:5], pop1_sel4_gv[,1:5], pop1_sel5_gv[,1:5], pop1_sel6_gv[,1:5])
wt_gain <- c(pop1_sel3_gv[,1:10]-pop1_sel2_gv[,1:10], pop1_sel4_gv[,1:10]-pop1_sel2_gv[,1:10], pop1_sel5_gv[,1:10]-pop1_sel2_gv[,1:10], pop1_sel6_gv[,1:10]-pop1_sel2_gv[,1:10])

ddm1 <- c(ddm1_sel2_gv[,1:5], ddm1_sel3_gv[,1:5], ddm1_sel4_gv[,1:5], ddm1_sel5_gv[,1:5], ddm1_sel6_gv[,1:5])
ddm1_gain <- c(ddm1_sel3_gv[,1:10] - ddm1_sel2_gv[,1:10], ddm1_sel4_gv[,1:10]- ddm1_sel2_gv[,1:10], ddm1_sel5_gv[,1:10] - ddm1_sel2_gv[,1:10], ddm1_sel6_gv[,1:10] - ddm1_sel2_gv[,1:10])

zmet2 <- c(zmet2_sel2_gv[,1:5], zmet2_sel3_gv[,1:5], zmet2_sel4_gv[,1:5], zmet2_sel5_gv[,1:5], zmet2_sel6_gv[,1:5])
zmet2_gain <- c(zmet2_sel3_gv[,1:10]-zmet2_sel2_gv[,1:10], zmet2_sel4_gv[,1:10]-zmet2_sel2_gv[,1:10], zmet2_sel5_gv[,1:10]-zmet2_sel2_gv[,1:10], zmet2_sel6_gv[,1:10] - zmet2_sel2_gv[,1:10])

ideal1 <- c(ideal1_sel2_gv[,1:5], ideal1_sel3_gv[,1:5], ideal1_sel4_gv[,1:5], ideal1_sel5_gv[,1:5], ideal1_sel6_gv[,1:5])
ideal1_gain <- c(ideal1_sel3_gv[,1:10]-ideal1_sel2_gv[,1:10], ideal1_sel4_gv[,1:10]-ideal1_sel2_gv[,1:10], ideal1_sel5_gv[,1:10]-ideal1_sel2_gv[,1:10], ideal1_sel6_gv[,1:10]-ideal1_sel2_gv[,1:10])

ideal2 <- c(ideal2_sel2_gv[,1:5], ideal2_sel3_gv[,1:5], ideal2_sel4_gv[,1:5], ideal2_sel5_gv[,1:5], ideal2_sel6_gv[,1:5])
ideal2_gain <- c(ideal2_sel3_gv[,1:10]-ideal2_sel2_gv[,1:10], ideal2_sel4_gv[,1:10]-ideal2_sel2_gv[,1:10], ideal2_sel5_gv[,1:10]-ideal2_sel2_gv[,1:10], ideal2_sel6_gv[,1:10]-ideal2_sel2_gv[,1:10])

fancm <- c(fancm_sel2_gv[,1:5], fancm_sel3_gv[,1:5], fancm_sel4_gv[,1:5], fancm_sel5_gv[,1:5], fancm_sel6_gv[,1:5])
fancm_gain <- c(fancm_sel3_gv[,1:10]-fancm_sel2_gv[,1:10], fancm_sel4_gv[,1:10]-fancm_sel2_gv[,1:10], fancm_sel5_gv[,1:10]-fancm_sel2_gv[,1:10], fancm_sel6_gv[,1:10]-fancm_sel2_gv[,1:10])

recq4_gain <- c(recq4_sel3_gv[,1:10]-recq4_sel2_gv[,1:10], recq4_sel4_gv[,1:10]-recq4_sel2_gv[,1:10], recq4_sel5_gv[,1:10]-recq4_sel2_gv[,1:10], recq4_sel6_gv[,1:10]-recq4_sel2_gv[,1:10])


wt2 <- c(mean(pop1_sel2_gv[,11:20]), mean(pop1_sel3_gv[,11:20]), mean(pop1_sel4_gv[,11:20]), mean(pop1_sel5_gv[,11:20]), mean(pop1_sel6_gv[,11:20]))

ddm12 <- c(mean(ddm1_sel2_gv[,11:20]), mean(ddm1_sel3_gv[,11:20]), mean(ddm1_sel4_gv[,11:20]), mean(ddm1_sel5_gv[,11:20]), mean(ddm1_sel6_gv[,11:20]))

zmet22 <- c(mean(zmet2_sel2_gv[,11:20]), mean(zmet2_sel3_gv[,11:20]), mean(zmet2_sel4_gv[,11:20]), mean(zmet2_sel5_gv[,11:20]), mean(zmet2_sel6_gv[,11:20]))

recq42 <- c(mean(recq4_sel2_gv[,11:20]), mean(recq4_sel3_gv[,11:20]), mean(recq4_sel4_gv[,11:20]), mean(recq4_sel5_gv[,11:20]), mean(recq4_sel6_gv[,11:20]))

ideal12 <- c(mean(ideal1_sel2_gv[,11:20]), mean(ideal1_sel3_gv[,11:20]), mean(ideal1_sel4_gv[,11:20]), mean(ideal1_sel5_gv[,11:20]), mean(ideal1_sel6_gv[,11:20]))

ideal22 <- c(mean(ideal2_sel2_gv[,11:20]), mean(ideal2_sel3_gv[,11:20]), mean(ideal2_sel4_gv[,11:20]), mean(ideal2_sel5_gv[,11:20]), mean(ideal2_sel6_gv[,11:20]))

fancm2 <- c(mean(fancm_sel2_gv[,11:20]), mean(fancm_sel3_gv[,11:20]), mean(fancm_sel4_gv[,11:20]), mean(fancm_sel5_gv[,11:20]), mean(fancm_sel6_gv[,11:20]))

#adding all genetic values
all <- cbind(wt, ddm1, zmet2, recq4)
all <- as.data.frame(all)
all2 <- as.data.frame(matrix(data = NA, nrow = 20000))
all2$gv <- c(all$wt, all$ddm1, all$zmet2, all$recq4)
all2$generation <- rep(1:5, size = 4, each = 1000)
all2$gen_map <- rep(c("wt","ddm1", "zmet2", "recq4"), size = 4, each = 5000)
gen2 <- all2[which(all2$generation == 5), ]

gain <- cbind(wt_gain, ddm1_gain, zmet2_gain, recq4_gain, fancm_gain, ideal1_gain, ideal2_gain)
gain <- as.data.frame(gain)
gain2 <- as.data.frame(matrix(data = NA, nrow = 28000))
gain2$gv <- c(gain$wt_gain, gain$ddm1_gain, gain$zmet2_gain, gain$recq4_gain, gain$fancm_gain, gain$ideal1_gain, gain$ideal2_gain)
gain2$generation <- rep(1:4, size = 7, each = 1000)
gain2$gen_map <- rep(c("wt","ddm1", "zmet2", "recq4", "fancm", "ideal1", "ideal2"), size = 7, each = 4000)
lastgen_gain <- gain2[which(gain2$generation == 4), ]

#looking at introgression of resistance loci
trait2 <- cbind(wt2, ddm12, zmet22, recq42)
trait2 <- as.data.frame(trait2)
trait22 <- as.data.frame(matrix(data = NA, nrow = 20))
trait22$gv <- c(trait2$wt2, trait2$ddm12, trait2$zmet22, trait2$recq42)
trait22$generation <- rep(1:5, size = 4, each = 1)
trait22$gen_map <- rep(c("wt","ddm1", "zmet2", "recq4"), size = 4, each = 1)
intro_lastgen <- trait22[which(trait22$generation == 5), ]

#many extreme values
ggqqplot(residuals(reduced_introgression))

ggplot(all2, aes(x=as.factor(generation), y= gv, fill=gen_map)) + 
  geom_boxplot() + theme_bw() + xlab("Generation") + ylab("GVs") + 
  scale_fill_manual(values=c("#aae4c2", "#fff87a", "#35c7ff", '#ff8b77')) + ggtitle("Introgression Genetic Values through 6 Generations")+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=14), legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14), plot.title = element_text(size=20))

colors = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")
ggplot(gain2, aes(x=as.factor(generation), y= gv, fill=gen_map)) + 
  geom_boxplot() + theme_bw() + xlab("Generations since Gen 0") + ylab("Genetic Gain from Gen 0") + 
  scale_fill_manual(values=group.colors) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=14), legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14), plot.title = element_text(size=20), legend.position = c(0.8,0.2)) +ylim(0.5,2)

group.colors = c("ddm1" = "#E69F00", "recq4" = "#56B4E9", "wt" = "#009E73", "zmet2" = "#F0E442",
                 "fancm" = "#0072B2", "ideal1" = "#D55E00", "ideal2" = "#CC79A7")
ggboxplot(lastgen_gain, x="gen_map", y= "gv", fill = "gen_map") + 
  theme_bw() + xlab("Genetic Map Used") + ylab("Genetic Gain") +
  scale_fill_manual(values = group.colors) + ggtitle("Boxplot of genetic gain with 10% selection intensity") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=14), legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14), plot.title = element_text(size=20), legend.position="non") + ylim(1,2.0)

ggplot(trait22, aes(x=as.factor(generation), y= gv, color=gen_map)) + 
  geom_point() + theme_bw() + xlab("Generations since gen 0") + ylab("GVs")
#scale_fill_manual(values=c("#aae4c2", "#fff87a", "#35c7ff", "#ff8b77"))

allvar <- cbind(wtvarall, ddm1varall, zmet2varall, recq4varall)
allvar <- as.data.frame(allvar)
allvar2 <- as.data.frame(matrix(data = NA, nrow = 9600))
allvar2$var <- c(allvar$wtvarall, allvar$ddm1varall, allvar$zmet2varall, allvar$recq4varall)
allvar2$gen <- rep(1:5, size = 4, each = 480)
allvar2$gen_map <- rep(c("wt","ddm1", "zmet", "recq4"), size = 4, each = 2400)
lastvar <- allvar2[which(allvar2$gen == 5),]

introgressionvar <- aov(var ~ as.factor(gen_map), data = lastvar)
summary(introgressionvar)
TukeyHSD(introgressionvar, conf.level = .95)

ggqqplot(residuals(introgressionvar))

ggplot(allvar2, aes(x=as.factor(gen), y=var, fill=gen_map)) + 
  geom_boxplot() + theme_bw() + xlab("Generation") + ylab("Genetic Variances") + 
  scale_fill_manual(values=c("#aae4c2", "#fff87a", "#35c7ff",'#ff8b77')) + ggtitle("Introgression Genetic Variance through 6 Generations") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=14), legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14), plot.title = element_text(size=20)) + ylim(-0.05,0.2)