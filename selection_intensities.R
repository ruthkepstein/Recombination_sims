library(AlphaSimR)
library(Rcpp)
library(ggplot2)
library(dplyr)
set.seed(420)

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

#keep overall selected indivuduals to 5, change how many u make from selected
segSites <- c(300, 237, 219, 256, 146,  92, 233, 173, 149, 184)

#burn-in to create LD amongst SNPs
burn_in_pop <- vector(mode = "list", length = 100)
for(i in 1:100){
  founderPop <- quickHaplo(nInd = 200, nChr = 10, inbred = TRUE, ploidy = 2L, segSites = segSites)
  founderPop@genMap <- final_map
  founderPop@centromere <- real_centromere
  SP = SimParam$new(founderPop)
  SP$setTrackRec(TRUE)
  SP$p = 0.15
  #polygenic trait with QTLs in accurate gene space
  trait_yield <- new("TraitA", nLoci = 300L, lociPerChr= c(40L, 35L, 35L, 35L, 35L, 30L, 25L, 25L, 20L, 20L),
                     lociLoc = as.integer(c(chr1_QTL, chr2_QTL, chr3_QTL, chr4_QTL, chr5_QTL,
                                            chr6_QTL, chr7_QTL, chr8_QTL, chr9_QTL, chr10_QTL)), addEff = addEff_mix, intercept = 0.1)
  SP$manAddTrait(trait_yield)
  pop_good <- newPop(founderPop, simParam = SP)
  pop_good <- setPheno(pop_good, h2 = 0.8, simParam = SP)
  
  pop_good1 <- randCross(pop_good, nCrosses = 5, nProgeny=20, simParam = SP)
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
  
  pop_good_sel10 <- selectInd(pop_good10, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good11 <- randCross(pop_good_sel10, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good11 <- setPheno(pop_good11, h2 = 0.8, simParam = SP)
  
  burn_in_pop[[i]] <- selectInd(pop_good10, nInd = 5, use = "gv", trait = 1, selectop = TRUE, returnPop = TRUE, simParam = SP)
}

#combine all 100 iterations together to make one population with 200 individuals
burn_in_pop <- mergePops(burn_in_pop)

##Starting intermating schemes for different genetic maps
#wild-type map
wt2var1 <- matrix(data = NA, ncol = 1, nrow = 100)
wt2var2 <- matrix(data = NA, ncol = 1, nrow = 100)
wt2var3 <- matrix(data = NA, ncol = 1, nrow = 100)
wt2var4 <- matrix(data = NA, ncol = 1, nrow = 100)
wt2var5 <- matrix(data = NA, ncol = 1, nrow = 100)
wt2var6 <- matrix(data = NA, ncol = 1, nrow = 100)
wt2var7 <- matrix(data = NA, ncol = 1, nrow = 100)
wt2var8 <- matrix(data = NA, ncol = 1, nrow = 100)
wt2var9 <- matrix(data = NA, ncol = 1, nrow = 100)
wt2var10 <- matrix(data = NA, ncol = 1, nrow = 100)
wt2var11 <- matrix(data = NA, ncol = 1, nrow = 100)
wt2var12 <- matrix(data = NA, ncol = 1, nrow = 100)
wt2var13 <- matrix(data = NA, ncol = 1, nrow = 100)
wt2var14 <- matrix(data = NA, ncol = 1, nrow = 100)
wt2var15 <- matrix(data = NA, ncol = 1, nrow = 100)
wt2var16 <- matrix(data = NA, ncol = 1, nrow = 100)
wt2var17 <- matrix(data = NA, ncol = 1, nrow = 100)
wt2gv1 <- matrix(data = NA, ncol = 2, nrow = 100)
wt2gv2 <- matrix(data = NA, ncol = 2, nrow = 100)
wt2gv3 <- matrix(data = NA, ncol = 2, nrow = 100)
wt2gv4 <- matrix(data = NA, ncol = 2, nrow = 100)
wt2gv5 <- matrix(data = NA, ncol = 2, nrow = 100)
wt2gv6 <- matrix(data = NA, ncol = 2, nrow = 100)
wt2gv7 <- matrix(data = NA, ncol = 2, nrow = 100)
wt2gv8 <- matrix(data = NA, ncol = 2, nrow = 100)
wt2gv9 <- matrix(data = NA, ncol = 2, nrow = 100)
wt2gv10 <- matrix(data = NA, ncol = 2, nrow = 100)
wt2gv11 <- matrix(data = NA, ncol = 2, nrow = 100)
wt2gv12 <- matrix(data = NA, ncol = 2, nrow = 100)
wt2gv13 <- matrix(data = NA, ncol = 2, nrow = 100)
wt2gv14 <- matrix(data = NA, ncol = 2, nrow = 100)
wt2gv15 <- matrix(data = NA, ncol = 2, nrow = 100)
for(i in 1:100){
  SP$switchGenMap(final_map, centromere = real_centromere)
  wt2var1[i,] <- genicVarA(burn_in_pop)
  
  pop_good1 <- randCross(burn_in_pop, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good1 <- setPheno(pop_good1, h2 = 0.8, simParam = SP)
  wt2var2[i,] <- genicVarA(pop_good1)
  
  pop_good_sel <- selectInd(pop_good1, nInd = 2, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good2 <- randCross(pop_good_sel, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good2 <- setPheno(pop_good2, h2 = 0.8, simParam = SP)
  wt2var3[i,] <- genicVarA(pop_good2)
  wt2gv1[i,] <- gv(pop_good_sel)
  
  pop_good_sel2 <- selectInd(pop_good2, nInd = 2, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good3 <- randCross(pop_good_sel2, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good3 <- setPheno(pop_good3, h2 = 0.8, simParam = SP)
  wt2var4[i,] <- genicVarA(pop_good3)
  wt2gv2[i,] <- gv(pop_good_sel2)
  
  pop_good_sel3 <- selectInd(pop_good3, nInd = 2, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good4 <- randCross(pop_good_sel2, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good4 <- setPheno(pop_good4, h2 = 0.8, simParam = SP)
  wt2var5[i,] <- genicVarA(pop_good4)
  wt2gv3[i,] <- gv(pop_good_sel3)
  
  pop_good_sel4 <- selectInd(pop_good4, nInd = 2, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good5 <- randCross(pop_good_sel4, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good5 <- setPheno(pop_good5, h2 = 0.8, simParam = SP)
  wt2var6[i,] <- genicVarA(pop_good5)
  wt2gv4[i,] <- gv(pop_good_sel4)
  
  pop_good_sel5 <- selectInd(pop_good5, nInd = 2, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good6 <- randCross(pop_good_sel5, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good6 <- setPheno(pop_good6, h2 = 0.8, simParam = SP)
  wt2var7[i,] <- genicVarA(pop_good6)
  wt2gv5[i,] <- gv(pop_good_sel5)
  
  pop_good_sel6 <- selectInd(pop_good6, nInd = 2, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good7 <- randCross(pop_good_sel6, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good7 <- setPheno(pop_good7, h2 = 0.8, simParam = SP)
  wt2var8[i,] <- genicVarA(pop_good7)
  wt2gv6[i,] <- gv(pop_good_sel6)
  
  pop_good_sel7 <- selectInd(pop_good7, nInd = 2, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good8 <- randCross(pop_good_sel7, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good8 <- setPheno(pop_good8, h2 = 0.8, simParam = SP)
  wt2var9[i,] <- genicVarA(pop_good8)
  wt2gv7[i,] <- gv(pop_good_sel7)
  
  pop_good_sel8 <- selectInd(pop_good8, nInd = 2, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good9 <- randCross(pop_good_sel8, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good9 <- setPheno(pop_good9, h2 = 0.8, simParam = SP)
  wt2var10[i,] <- genicVarA(pop_good9)
  wt2gv8[i,] <- gv(pop_good_sel8)
  
  pop_good_sel9 <- selectInd(pop_good9, nInd = 2, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good10 <- randCross(pop_good_sel9, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good10 <- setPheno(pop_good10, h2 = 0.8, simParam = SP)
  wt2var11[i,] <- genicVarA(pop_good10)
  wt2gv9[i,] <- gv(pop_good_sel9)
  
  pop_good_sel10 <- selectInd(pop_good10, nInd = 2, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good11 <- randCross(pop_good_sel10, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good11 <- setPheno(pop_good11, h2 = 0.8, simParam = SP)
  wt2var12[i,] <- genicVarA(pop_good11)
  wt2gv10[i,] <- gv(pop_good_sel10)
  
  pop_good_sel11 <- selectInd(pop_good11, nInd = 2, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good12 <- randCross(pop_good_sel11, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good12 <- setPheno(pop_good12, h2 = 0.8, simParam = SP)
  wt2var13[i,] <- genicVarA(pop_good12)
  wt2gv11[i,] <- gv(pop_good_sel11)
  
  pop_good_sel12 <- selectInd(pop_good12, nInd = 2, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good13 <- randCross(pop_good_sel12, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good13 <- setPheno(pop_good13, h2 = 0.8, simParam = SP)
  wt2var14[i,] <- genicVarA(pop_good13)
  wt2gv12[i,] <- gv(pop_good_sel12)
  
  pop_good_sel13 <- selectInd(pop_good13, nInd = 2, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good14 <- randCross(pop_good_sel13, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good14 <- setPheno(pop_good14, h2 = 0.8, simParam = SP)
  wt2var15[i,] <- genicVarA(pop_good14)
  wt2gv13[i,] <- gv(pop_good_sel13)
  
  pop_good_sel14 <- selectInd(pop_good14, nInd = 2, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good15 <- randCross(pop_good_sel14, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good15 <- setPheno(pop_good15, h2 = 0.8, simParam = SP)
  wt2var16[i,] <- genicVarA(pop_good15)
  wt2gv14[i,] <- gv(pop_good_sel14)
  
  pop_good_sel15 <- selectInd(pop_good15, nInd = 2, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good16 <- randCross(pop_good_sel14, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good16 <- setPheno(pop_good15, h2 = 0.8, simParam = SP)
  wt2var17[i,] <- genicVarA(pop_good16)
  wt2gv15[i,] <- gv(pop_good_sel15)
}

wt5var1 <- matrix(data = NA, ncol = 1, nrow = 100)
wt5var2 <- matrix(data = NA, ncol = 1, nrow = 100)
wt5var3 <- matrix(data = NA, ncol = 1, nrow = 100)
wt5var4 <- matrix(data = NA, ncol = 1, nrow = 100)
wt5var5 <- matrix(data = NA, ncol = 1, nrow = 100)
wt5var6 <- matrix(data = NA, ncol = 1, nrow = 100)
wt5var7 <- matrix(data = NA, ncol = 1, nrow = 100)
wt5var8 <- matrix(data = NA, ncol = 1, nrow = 100)
wt5var9 <- matrix(data = NA, ncol = 1, nrow = 100)
wt5var10 <- matrix(data = NA, ncol = 1, nrow = 100)
wt5var11 <- matrix(data = NA, ncol = 1, nrow = 100)
wt5var12 <- matrix(data = NA, ncol = 1, nrow = 100)
wt5var13 <- matrix(data = NA, ncol = 1, nrow = 100)
wt5var14 <- matrix(data = NA, ncol = 1, nrow = 100)
wt5var15 <- matrix(data = NA, ncol = 1, nrow = 100)
wt5var16 <- matrix(data = NA, ncol = 1, nrow = 100)
wt5var17 <- matrix(data = NA, ncol = 1, nrow = 100)
wt5gv1 <- matrix(data = NA, ncol = 5, nrow = 100)
wt5gv2 <- matrix(data = NA, ncol = 5, nrow = 100)
wt5gv3 <- matrix(data = NA, ncol = 5, nrow = 100)
wt5gv4 <- matrix(data = NA, ncol = 5, nrow = 100)
wt5gv5 <- matrix(data = NA, ncol = 5, nrow = 100)
wt5gv6 <- matrix(data = NA, ncol = 5, nrow = 100)
wt5gv7 <- matrix(data = NA, ncol = 5, nrow = 100)
wt5gv8 <- matrix(data = NA, ncol = 5, nrow = 100)
wt5gv9 <- matrix(data = NA, ncol = 5, nrow = 100)
wt5gv10 <- matrix(data = NA, ncol = 5, nrow = 100)
wt5gv11 <- matrix(data = NA, ncol = 5, nrow = 100)
wt5gv12 <- matrix(data = NA, ncol = 5, nrow = 100)
wt5gv13 <- matrix(data = NA, ncol = 5, nrow = 100)
wt5gv14 <- matrix(data = NA, ncol = 5, nrow = 100)
wt5gv15 <- matrix(data = NA, ncol = 5, nrow = 100)
for(i in 1:100){
  SP$switchGenMap(final_map, centromere = real_centromere)
  wt5var1[i,] <- genicVarA(burn_in_pop)
  
  pop_good1 <- randCross(burn_in_pop, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good1 <- setPheno(pop_good1, h2 = 0.8, simParam = SP)
  wt5var2[i,] <- genicVarA(pop_good1)
  
  pop_good_sel <- selectInd(pop_good1, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good2 <- randCross(pop_good_sel, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good2 <- setPheno(pop_good2, h2 = 0.8, simParam = SP)
  wt5var3[i,] <- genicVarA(pop_good2)
  wt5gv1[i,] <- gv(pop_good_sel)
  
  pop_good_sel2 <- selectInd(pop_good2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good3 <- randCross(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good3 <- setPheno(pop_good3, h2 = 0.8, simParam = SP)
  wt5var4[i,] <- genicVarA(pop_good3)
  wt5gv2[i,] <- gv(pop_good_sel2)
  
  pop_good_sel3 <- selectInd(pop_good3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good4 <- randCross(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good4 <- setPheno(pop_good4, h2 = 0.8, simParam = SP)
  wt5var5[i,] <- genicVarA(pop_good4)
  wt5gv3[i,] <- gv(pop_good_sel3)
  
  pop_good_sel4 <- selectInd(pop_good4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good5 <- randCross(pop_good_sel4, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good5 <- setPheno(pop_good5, h2 = 0.8, simParam = SP)
  wt5var6[i,] <- genicVarA(pop_good5)
  wt5gv4[i,] <- gv(pop_good_sel4)
  
  pop_good_sel5 <- selectInd(pop_good5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good6 <- randCross(pop_good_sel5, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good6 <- setPheno(pop_good6, h2 = 0.8, simParam = SP)
  wt5var7[i,] <- genicVarA(pop_good6)
  wt5gv5[i,] <- gv(pop_good_sel5)
  
  pop_good_sel6 <- selectInd(pop_good6, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good7 <- randCross(pop_good_sel6, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good7 <- setPheno(pop_good7, h2 = 0.8, simParam = SP)
  wt5var8[i,] <- genicVarA(pop_good7)
  wt5gv6[i,] <- gv(pop_good_sel6)
  
  pop_good_sel7 <- selectInd(pop_good7, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good8 <- randCross(pop_good_sel7, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good8 <- setPheno(pop_good8, h2 = 0.8, simParam = SP)
  wt5var9[i,] <- genicVarA(pop_good8)
  wt5gv7[i,] <- gv(pop_good_sel7)
  
  pop_good_sel8 <- selectInd(pop_good8, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good9 <- randCross(pop_good_sel8, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good9 <- setPheno(pop_good9, h2 = 0.8, simParam = SP)
  wt5var10[i,] <- genicVarA(pop_good9)
  wt5gv8[i,] <- gv(pop_good_sel8)
  
  pop_good_sel9 <- selectInd(pop_good9, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good10 <- randCross(pop_good_sel9, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good10 <- setPheno(pop_good10, h2 = 0.8, simParam = SP)
  wt5var11[i,] <- genicVarA(pop_good10)
  wt5gv9[i,] <- gv(pop_good_sel9)
  
  pop_good_sel10 <- selectInd(pop_good10, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good11 <- randCross(pop_good_sel10, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good11 <- setPheno(pop_good11, h2 = 0.8, simParam = SP)
  wt5var12[i,] <- genicVarA(pop_good11)
  wt5gv10[i,] <- gv(pop_good_sel10)
  
  pop_good_sel11 <- selectInd(pop_good11, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good12 <- randCross(pop_good_sel11, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good12 <- setPheno(pop_good12, h2 = 0.8, simParam = SP)
  wt5var13[i,] <- genicVarA(pop_good12)
  wt5gv11[i,] <- gv(pop_good_sel11)
  
  pop_good_sel12 <- selectInd(pop_good12, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good13 <- randCross(pop_good_sel12, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good13 <- setPheno(pop_good13, h2 = 0.8, simParam = SP)
  wt5var14[i,] <- genicVarA(pop_good13)
  wt5gv12[i,] <- gv(pop_good_sel12)
  
  pop_good_sel13 <- selectInd(pop_good13, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good14 <- randCross(pop_good_sel13, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good14 <- setPheno(pop_good14, h2 = 0.8, simParam = SP)
  wt5var15[i,] <- genicVarA(pop_good14)
  wt5gv13[i,] <- gv(pop_good_sel13)
  
  pop_good_sel14 <- selectInd(pop_good14, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good15 <- randCross(pop_good_sel14, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good15 <- setPheno(pop_good15, h2 = 0.8, simParam = SP)
  wt5var16[i,] <- genicVarA(pop_good15)
  wt5gv14[i,] <- gv(pop_good_sel14)
  
  pop_good_sel15 <- selectInd(pop_good15, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good16 <- randCross(pop_good_sel14, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good16 <- setPheno(pop_good15, h2 = 0.8, simParam = SP)
  wt5var17[i,] <- genicVarA(pop_good16)
  wt5gv15[i,] <- gv(pop_good_sel15)
}

wt10var1 <- matrix(data = NA, ncol = 1, nrow = 100)
wt10var2 <- matrix(data = NA, ncol = 1, nrow = 100)
wt10var3 <- matrix(data = NA, ncol = 1, nrow = 100)
wt10var4 <- matrix(data = NA, ncol = 1, nrow = 100)
wt10var5 <- matrix(data = NA, ncol = 1, nrow = 100)
wt10var6 <- matrix(data = NA, ncol = 1, nrow = 100)
wt10var7 <- matrix(data = NA, ncol = 1, nrow = 100)
wt10var8 <- matrix(data = NA, ncol = 1, nrow = 100)
wt10var9 <- matrix(data = NA, ncol = 1, nrow = 100)
wt10var10 <- matrix(data = NA, ncol = 1, nrow = 100)
wt10var11 <- matrix(data = NA, ncol = 1, nrow = 100)
wt10var12 <- matrix(data = NA, ncol = 1, nrow = 100)
wt10var13 <- matrix(data = NA, ncol = 1, nrow = 100)
wt10var14 <- matrix(data = NA, ncol = 1, nrow = 100)
wt10var15 <- matrix(data = NA, ncol = 1, nrow = 100)
wt10var16 <- matrix(data = NA, ncol = 1, nrow = 100)
wt10var17 <- matrix(data = NA, ncol = 1, nrow = 100)
wt10gv1 <- matrix(data = NA, ncol = 10, nrow = 100)
wt10gv2 <- matrix(data = NA, ncol = 10, nrow = 100)
wt10gv3 <- matrix(data = NA, ncol = 10, nrow = 100)
wt10gv4 <- matrix(data = NA, ncol = 10, nrow = 100)
wt10gv5 <- matrix(data = NA, ncol = 10, nrow = 100)
wt10gv6 <- matrix(data = NA, ncol = 10, nrow = 100)
wt10gv7 <- matrix(data = NA, ncol = 10, nrow = 100)
wt10gv8 <- matrix(data = NA, ncol = 10, nrow = 100)
wt10gv9 <- matrix(data = NA, ncol = 10, nrow = 100)
wt10gv10 <- matrix(data = NA, ncol = 10, nrow = 100)
wt10gv11 <- matrix(data = NA, ncol = 10, nrow = 100)
wt10gv12 <- matrix(data = NA, ncol = 10, nrow = 100)
wt10gv13 <- matrix(data = NA, ncol = 10, nrow = 100)
wt10gv14 <- matrix(data = NA, ncol = 10, nrow = 100)
wt10gv15 <- matrix(data = NA, ncol = 10, nrow = 100)
for(i in 1:100){
  SP$switchGenMap(final_map, centromere = real_centromere)
  wt10var1[i,] <- genicVarA(burn_in_pop)
  
  pop_good1 <- randCross(burn_in_pop, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good1 <- setPheno(pop_good1, h2 = 0.8, simParam = SP)
  wt10var2[i,] <- genicVarA(pop_good1)
  
  pop_good_sel <- selectInd(pop_good1, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good2 <- randCross(pop_good_sel, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good2 <- setPheno(pop_good2, h2 = 0.8, simParam = SP)
  wt10var3[i,] <- genicVarA(pop_good2)
  wt10gv1[i,] <- gv(pop_good_sel)
  
  pop_good_sel2 <- selectInd(pop_good2, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good3 <- randCross(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good3 <- setPheno(pop_good3, h2 = 0.8, simParam = SP)
  wt10var4[i,] <- genicVarA(pop_good3)
  wt10gv2[i,] <- gv(pop_good_sel2)
  
  pop_good_sel3 <- selectInd(pop_good3, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good4 <- randCross(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good4 <- setPheno(pop_good4, h2 = 0.8, simParam = SP)
  wt10var5[i,] <- genicVarA(pop_good4)
  wt10gv3[i,] <- gv(pop_good_sel3)
  
  pop_good_sel4 <- selectInd(pop_good4, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good5 <- randCross(pop_good_sel4, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good5 <- setPheno(pop_good5, h2 = 0.8, simParam = SP)
  wt10var6[i,] <- genicVarA(pop_good5)
  wt10gv4[i,] <- gv(pop_good_sel4)
  
  pop_good_sel5 <- selectInd(pop_good5, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good6 <- randCross(pop_good_sel5, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good6 <- setPheno(pop_good6, h2 = 0.8, simParam = SP)
  wt10var7[i,] <- genicVarA(pop_good6)
  wt10gv5[i,] <- gv(pop_good_sel5)
  
  pop_good_sel6 <- selectInd(pop_good6, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good7 <- randCross(pop_good_sel6, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good7 <- setPheno(pop_good7, h2 = 0.8, simParam = SP)
  wt10var8[i,] <- genicVarA(pop_good7)
  wt10gv6[i,] <- gv(pop_good_sel6)
  
  pop_good_sel7 <- selectInd(pop_good7, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good8 <- randCross(pop_good_sel7, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good8 <- setPheno(pop_good8, h2 = 0.8, simParam = SP)
  wt10var9[i,] <- genicVarA(pop_good8)
  wt10gv7[i,] <- gv(pop_good_sel7)
  
  pop_good_sel8 <- selectInd(pop_good8, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good9 <- randCross(pop_good_sel8, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good9 <- setPheno(pop_good9, h2 = 0.8, simParam = SP)
  wt10var10[i,] <- genicVarA(pop_good9)
  wt10gv8[i,] <- gv(pop_good_sel8)
  
  pop_good_sel9 <- selectInd(pop_good9, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good10 <- randCross(pop_good_sel9, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good10 <- setPheno(pop_good10, h2 = 0.8, simParam = SP)
  wt10var11[i,] <- genicVarA(pop_good10)
  wt10gv9[i,] <- gv(pop_good_sel9)
  
  pop_good_sel10 <- selectInd(pop_good10, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good11 <- randCross(pop_good_sel10, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good11 <- setPheno(pop_good11, h2 = 0.8, simParam = SP)
  wt10var12[i,] <- genicVarA(pop_good11)
  wt10gv10[i,] <- gv(pop_good_sel10)
  
  pop_good_sel11 <- selectInd(pop_good11, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good12 <- randCross(pop_good_sel11, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good12 <- setPheno(pop_good12, h2 = 0.8, simParam = SP)
  wt10var13[i,] <- genicVarA(pop_good12)
  wt10gv11[i,] <- gv(pop_good_sel11)
  
  pop_good_sel12 <- selectInd(pop_good12, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good13 <- randCross(pop_good_sel12, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good13 <- setPheno(pop_good13, h2 = 0.8, simParam = SP)
  wt10var14[i,] <- genicVarA(pop_good13)
  wt10gv12[i,] <- gv(pop_good_sel12)
  
  pop_good_sel13 <- selectInd(pop_good13, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good14 <- randCross(pop_good_sel13, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good14 <- setPheno(pop_good14, h2 = 0.8, simParam = SP)
  wt10var15[i,] <- genicVarA(pop_good14)
  wt10gv13[i,] <- gv(pop_good_sel13)
  
  pop_good_sel14 <- selectInd(pop_good14, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good15 <- randCross(pop_good_sel14, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good15 <- setPheno(pop_good15, h2 = 0.8, simParam = SP)
  wt10var16[i,] <- genicVarA(pop_good15)
  wt10gv14[i,] <- gv(pop_good_sel14)
  
  pop_good_sel15 <- selectInd(pop_good15, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good16 <- randCross(pop_good_sel14, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good16 <- setPheno(pop_good15, h2 = 0.8, simParam = SP)
  wt10var17[i,] <- genicVarA(pop_good16)
  wt10gv15[i,] <- gv(pop_good_sel15)
}

#with ddm1 genetic map
ddm12var1 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm12var2 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm12var3 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm12var4 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm12var5 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm12var6 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm12var7 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm12var8 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm12var9 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm12var10 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm12var11 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm12var12 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm12var13 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm12var14 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm12var15 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm12var16 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm12var17 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm12gv1 <- matrix(data = NA, ncol = 2, nrow = 100)
ddm12gv2 <- matrix(data = NA, ncol = 2, nrow = 100)
ddm12gv3 <- matrix(data = NA, ncol = 2, nrow = 100)
ddm12gv4 <- matrix(data = NA, ncol = 2, nrow = 100)
ddm12gv5 <- matrix(data = NA, ncol = 2, nrow = 100)
ddm12gv6 <- matrix(data = NA, ncol = 2, nrow = 100)
ddm12gv7 <- matrix(data = NA, ncol = 2, nrow = 100)
ddm12gv8 <- matrix(data = NA, ncol = 2, nrow = 100)
ddm12gv9 <- matrix(data = NA, ncol = 2, nrow = 100)
ddm12gv10 <- matrix(data = NA, ncol = 2, nrow = 100)
ddm12gv11 <- matrix(data = NA, ncol = 2, nrow = 100)
ddm12gv12 <- matrix(data = NA, ncol = 2, nrow = 100)
ddm12gv13 <- matrix(data = NA, ncol = 2, nrow = 100)
ddm12gv14 <- matrix(data = NA, ncol = 2, nrow = 100)
ddm12gv15 <- matrix(data = NA, ncol = 2, nrow = 100)
for(i in 1:100){
  SP$switchGenMap(ddm1_map, centromere = ddm1_centromere)
  ddm12var1[i,] <- genicVarA(burn_in_pop)
  
  pop_good1 <- randCross(burn_in_pop, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good1 <- setPheno(pop_good1, h2 = 0.8, simParam = SP)
  ddm12var2[i,] <- genicVarA(pop_good1)
  
  pop_good_sel <- selectInd(pop_good1, nInd = 2, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good2 <- randCross(pop_good_sel, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good2 <- setPheno(pop_good2, h2 = 0.8, simParam = SP)
  ddm12var3[i,] <- genicVarA(pop_good2)
  ddm12gv1[i,] <- gv(pop_good_sel)
  
  pop_good_sel2 <- selectInd(pop_good2, nInd = 2, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good3 <- randCross(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good3 <- setPheno(pop_good3, h2 = 0.8, simParam = SP)
  ddm12var4[i,] <- genicVarA(pop_good3)
  ddm12gv2[i,] <- gv(pop_good_sel2)
  
  pop_good_sel3 <- selectInd(pop_good3, nInd = 2, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good4 <- randCross(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good4 <- setPheno(pop_good4, h2 = 0.8, simParam = SP)
  ddm12var5[i,] <- genicVarA(pop_good4)
  ddm12gv3[i,] <- gv(pop_good_sel3)
  
  pop_good_sel4 <- selectInd(pop_good4, nInd = 2, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good5 <- randCross(pop_good_sel4, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good5 <- setPheno(pop_good5, h2 = 0.8, simParam = SP)
  ddm12var6[i,] <- genicVarA(pop_good5)
  ddm12gv4[i,] <- gv(pop_good_sel4)
  
  pop_good_sel5 <- selectInd(pop_good5, nInd = 2, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good6 <- randCross(pop_good_sel5, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good6 <- setPheno(pop_good6, h2 = 0.8, simParam = SP)
  ddm12var7[i,] <- genicVarA(pop_good6)
  ddm12gv5[i,] <- gv(pop_good_sel5)
  
  pop_good_sel6 <- selectInd(pop_good6, nInd = 2, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good7 <- randCross(pop_good_sel6, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good7 <- setPheno(pop_good7, h2 = 0.8, simParam = SP)
  ddm12var8[i,] <- genicVarA(pop_good7)
  ddm12gv6[i,] <- gv(pop_good_sel6)
  
  pop_good_sel7 <- selectInd(pop_good7, nInd = 2, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good8 <- randCross(pop_good_sel7, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good8 <- setPheno(pop_good8, h2 = 0.8, simParam = SP)
  ddm12var9[i,] <-genicVarA(pop_good8)
  ddm12gv7[i,] <- gv(pop_good_sel7)
  
  pop_good_sel8 <- selectInd(pop_good8, nInd = 2, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good9 <- randCross(pop_good_sel8, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good9 <- setPheno(pop_good9, h2 = 0.8, simParam = SP)
  ddm12var10[i,] <- genicVarA(pop_good9)
  ddm12gv8[i,] <- gv(pop_good_sel8)
  
  pop_good_sel9 <- selectInd(pop_good9, nInd = 2, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good10 <- randCross(pop_good_sel9, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good10 <- setPheno(pop_good10, h2 = 0.8, simParam = SP)
  ddm12var11[i,] <- genicVarA(pop_good10)
  ddm12gv9[i,] <- gv(pop_good_sel9)
  
  pop_good_sel10 <- selectInd(pop_good10, nInd = 2, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good11 <- randCross(pop_good_sel10, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good11 <- setPheno(pop_good11, h2 = 0.8, simParam = SP)
  ddm12var12[i,] <- genicVarA(pop_good11)
  ddm12gv10[i,] <- gv(pop_good_sel10)
  
  pop_good_sel11 <- selectInd(pop_good11, nInd = 2, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good12 <- randCross(pop_good_sel11, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good12 <- setPheno(pop_good12, h2 = 0.8, simParam = SP)
  ddm12var13[i,] <- genicVarA(pop_good12)
  ddm12gv11[i,] <- gv(pop_good_sel11)
  
  pop_good_sel12 <- selectInd(pop_good12, nInd = 2, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good13 <- randCross(pop_good_sel12, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good13 <- setPheno(pop_good13, h2 = 0.8, simParam = SP)
  ddm12var14[i,] <- genicVarA(pop_good13)
  ddm12gv12[i,] <- gv(pop_good_sel12)
  
  pop_good_sel13 <- selectInd(pop_good13, nInd = 2, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good14 <- randCross(pop_good_sel13, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good14 <- setPheno(pop_good14, h2 = 0.8, simParam = SP)
  ddm12var15[i,] <- genicVarA(pop_good14)
  ddm12gv13[i,] <- gv(pop_good_sel13)
  
  pop_good_sel14 <- selectInd(pop_good14, nInd = 2, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good15 <- randCross(pop_good_sel14, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good15 <- setPheno(pop_good15, h2 = 0.8, simParam = SP)
  ddm12var16[i,] <- genicVarA(pop_good15)
  ddm12gv14[i,] <- gv(pop_good_sel14)
  
  pop_good_sel15 <- selectInd(pop_good15, nInd = 2, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good16 <- randCross(pop_good_sel15, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good16 <- setPheno(pop_good16, h2 = 0.8, simParam = SP)
  ddm12var17[i,] <- genicVarA(pop_good16)
  ddm12gv15[i,] <- gv(pop_good_sel15)
}

ddm15var1 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm15var2 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm15var3 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm15var4 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm15var5 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm15var6 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm15var7 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm15var8 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm15var9 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm15var10 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm15var11 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm15var12 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm15var13 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm15var14 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm15var15 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm15var16 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm15var17 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm15gv1 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm15gv2 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm15gv3 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm15gv4 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm15gv5 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm15gv6 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm15gv7 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm15gv8 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm15gv9 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm15gv10 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm15gv11 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm15gv12 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm15gv13 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm15gv14 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm15gv15 <- matrix(data = NA, ncol = 5, nrow = 100)
for(i in 1:100){
  SP$switchGenMap(ddm1_map, centromere = ddm1_centromere)
  ddm15var1[i,] <- genicVarA(burn_in_pop)
  
  pop_good1 <- randCross(burn_in_pop, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good1 <- setPheno(pop_good1, h2 = 0.8, simParam = SP)
  ddm15var2[i,] <- genicVarA(pop_good1)
  
  pop_good_sel <- selectInd(pop_good1, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good2 <- randCross(pop_good_sel, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good2 <- setPheno(pop_good2, h2 = 0.8, simParam = SP)
  ddm15var3[i,] <- genicVarA(pop_good2)
  ddm15gv1[i,] <- gv(pop_good_sel)
  
  pop_good_sel2 <- selectInd(pop_good2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good3 <- randCross(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good3 <- setPheno(pop_good3, h2 = 0.8, simParam = SP)
  ddm15var4[i,] <- genicVarA(pop_good3)
  ddm15gv2[i,] <- gv(pop_good_sel2)
  
  pop_good_sel3 <- selectInd(pop_good3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good4 <- randCross(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good4 <- setPheno(pop_good4, h2 = 0.8, simParam = SP)
  ddm15var5[i,] <- genicVarA(pop_good4)
  ddm15gv3[i,] <- gv(pop_good_sel3)
  
  pop_good_sel4 <- selectInd(pop_good4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good5 <- randCross(pop_good_sel4, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good5 <- setPheno(pop_good5, h2 = 0.8, simParam = SP)
  ddm15var6[i,] <- genicVarA(pop_good5)
  ddm15gv4[i,] <- gv(pop_good_sel4)
  
  pop_good_sel5 <- selectInd(pop_good5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good6 <- randCross(pop_good_sel5, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good6 <- setPheno(pop_good6, h2 = 0.8, simParam = SP)
  ddm15var7[i,] <- genicVarA(pop_good6)
  ddm15gv5[i,] <- gv(pop_good_sel5)
  
  pop_good_sel6 <- selectInd(pop_good6, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good7 <- randCross(pop_good_sel6, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good7 <- setPheno(pop_good7, h2 = 0.8, simParam = SP)
  ddm15var8[i,] <- genicVarA(pop_good7)
  ddm15gv6[i,] <- gv(pop_good_sel6)
  
  pop_good_sel7 <- selectInd(pop_good7, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good8 <- randCross(pop_good_sel7, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good8 <- setPheno(pop_good8, h2 = 0.8, simParam = SP)
  ddm15var9[i,] <-genicVarA(pop_good8)
  ddm15gv7[i,] <- gv(pop_good_sel7)
  
  pop_good_sel8 <- selectInd(pop_good8, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good9 <- randCross(pop_good_sel8, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good9 <- setPheno(pop_good9, h2 = 0.8, simParam = SP)
  ddm15var10[i,] <- genicVarA(pop_good9)
  ddm15gv8[i,] <- gv(pop_good_sel8)
  
  pop_good_sel9 <- selectInd(pop_good9, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good10 <- randCross(pop_good_sel9, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good10 <- setPheno(pop_good10, h2 = 0.8, simParam = SP)
  ddm15var11[i,] <- genicVarA(pop_good10)
  ddm15gv9[i,] <- gv(pop_good_sel9)
  
  pop_good_sel10 <- selectInd(pop_good10, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good11 <- randCross(pop_good_sel10, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good11 <- setPheno(pop_good11, h2 = 0.8, simParam = SP)
  ddm15var12[i,] <- genicVarA(pop_good11)
  ddm15gv10[i,] <- gv(pop_good_sel10)
  
  pop_good_sel11 <- selectInd(pop_good11, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good12 <- randCross(pop_good_sel11, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good12 <- setPheno(pop_good12, h2 = 0.8, simParam = SP)
  ddm15var13[i,] <- genicVarA(pop_good12)
  ddm15gv11[i,] <- gv(pop_good_sel11)
  
  pop_good_sel12 <- selectInd(pop_good12, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good13 <- randCross(pop_good_sel12, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good13 <- setPheno(pop_good13, h2 = 0.8, simParam = SP)
  ddm15var14[i,] <- genicVarA(pop_good13)
  ddm15gv12[i,] <- gv(pop_good_sel12)
  
  pop_good_sel13 <- selectInd(pop_good13, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good14 <- randCross(pop_good_sel13, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good14 <- setPheno(pop_good14, h2 = 0.8, simParam = SP)
  ddm15var15[i,] <- genicVarA(pop_good14)
  ddm15gv13[i,] <- gv(pop_good_sel13)
  
  pop_good_sel14 <- selectInd(pop_good14, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good15 <- randCross(pop_good_sel14, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good15 <- setPheno(pop_good15, h2 = 0.8, simParam = SP)
  ddm15var16[i,] <- genicVarA(pop_good15)
  ddm15gv14[i,] <- gv(pop_good_sel14)
  
  pop_good_sel15 <- selectInd(pop_good15, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good16 <- randCross(pop_good_sel15, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good16 <- setPheno(pop_good16, h2 = 0.8, simParam = SP)
  ddm15var17[i,] <- genicVarA(pop_good16)
  ddm15gv15[i,] <- gv(pop_good_sel15)
}

ddm110var1 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm110var2 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm110var3 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm110var4 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm110var5 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm110var6 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm110var7 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm110var8 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm110var9 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm110var10 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm110var11 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm110var12 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm110var13 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm110var14 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm110var15 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm110var16 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm110var17 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm110gv1 <- matrix(data = NA, ncol = 10, nrow = 100)
ddm110gv2 <- matrix(data = NA, ncol = 10, nrow = 100)
ddm110gv3 <- matrix(data = NA, ncol = 10, nrow = 100)
ddm110gv4 <- matrix(data = NA, ncol = 10, nrow = 100)
ddm110gv5 <- matrix(data = NA, ncol = 10, nrow = 100)
ddm110gv6 <- matrix(data = NA, ncol = 10, nrow = 100)
ddm110gv7 <- matrix(data = NA, ncol = 10, nrow = 100)
ddm110gv8 <- matrix(data = NA, ncol = 10, nrow = 100)
ddm110gv9 <- matrix(data = NA, ncol = 10, nrow = 100)
ddm110gv10 <- matrix(data = NA, ncol = 10, nrow = 100)
ddm110gv11 <- matrix(data = NA, ncol = 10, nrow = 100)
ddm110gv12 <- matrix(data = NA, ncol = 10, nrow = 100)
ddm110gv13 <- matrix(data = NA, ncol = 10, nrow = 100)
ddm110gv14 <- matrix(data = NA, ncol = 10, nrow = 100)
ddm110gv15 <- matrix(data = NA, ncol = 10, nrow = 100)
for(i in 1:100){
  SP$switchGenMap(ddm1_map, centromere = ddm1_centromere)
  ddm110var1[i,] <- genicVarA(burn_in_pop)
  
  pop_good1 <- randCross(burn_in_pop, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good1 <- setPheno(pop_good1, h2 = 0.8, simParam = SP)
  ddm110var2[i,] <- genicVarA(pop_good1)
  
  pop_good_sel <- selectInd(pop_good1, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good2 <- randCross(pop_good_sel, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good2 <- setPheno(pop_good2, h2 = 0.8, simParam = SP)
  ddm110var3[i,] <- genicVarA(pop_good2)
  ddm110gv1[i,] <- gv(pop_good_sel)
  
  pop_good_sel2 <- selectInd(pop_good2, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good3 <- randCross(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good3 <- setPheno(pop_good3, h2 = 0.8, simParam = SP)
  ddm110var4[i,] <- genicVarA(pop_good3)
  ddm110gv2[i,] <- gv(pop_good_sel2)
  
  pop_good_sel3 <- selectInd(pop_good3, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good4 <- randCross(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good4 <- setPheno(pop_good4, h2 = 0.8, simParam = SP)
  ddm110var5[i,] <- genicVarA(pop_good4)
  ddm110gv3[i,] <- gv(pop_good_sel3)
  
  pop_good_sel4 <- selectInd(pop_good4, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good5 <- randCross(pop_good_sel4, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good5 <- setPheno(pop_good5, h2 = 0.8, simParam = SP)
  ddm110var6[i,] <- genicVarA(pop_good5)
  ddm110gv4[i,] <- gv(pop_good_sel4)
  
  pop_good_sel5 <- selectInd(pop_good5, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good6 <- randCross(pop_good_sel5, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good6 <- setPheno(pop_good6, h2 = 0.8, simParam = SP)
  ddm110var7[i,] <- genicVarA(pop_good6)
  ddm110gv5[i,] <- gv(pop_good_sel5)
  
  pop_good_sel6 <- selectInd(pop_good6, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good7 <- randCross(pop_good_sel6, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good7 <- setPheno(pop_good7, h2 = 0.8, simParam = SP)
  ddm110var8[i,] <- genicVarA(pop_good7)
  ddm110gv6[i,] <- gv(pop_good_sel6)
  
  pop_good_sel7 <- selectInd(pop_good7, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good8 <- randCross(pop_good_sel7, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good8 <- setPheno(pop_good8, h2 = 0.8, simParam = SP)
  ddm110var9[i,] <-genicVarA(pop_good8)
  ddm110gv7[i,] <- gv(pop_good_sel7)
  
  pop_good_sel8 <- selectInd(pop_good8, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good9 <- randCross(pop_good_sel8, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good9 <- setPheno(pop_good9, h2 = 0.8, simParam = SP)
  ddm110var10[i,] <- genicVarA(pop_good9)
  ddm110gv8[i,] <- gv(pop_good_sel8)
  
  pop_good_sel9 <- selectInd(pop_good9, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good10 <- randCross(pop_good_sel9, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good10 <- setPheno(pop_good10, h2 = 0.8, simParam = SP)
  ddm110var11[i,] <- genicVarA(pop_good10)
  ddm110gv9[i,] <- gv(pop_good_sel9)
  
  pop_good_sel10 <- selectInd(pop_good10, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good11 <- randCross(pop_good_sel10, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good11 <- setPheno(pop_good11, h2 = 0.8, simParam = SP)
  ddm110var12[i,] <- genicVarA(pop_good11)
  ddm110gv10[i,] <- gv(pop_good_sel10)
  
  pop_good_sel11 <- selectInd(pop_good11, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good12 <- randCross(pop_good_sel11, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good12 <- setPheno(pop_good12, h2 = 0.8, simParam = SP)
  ddm110var13[i,] <- genicVarA(pop_good12)
  ddm110gv11[i,] <- gv(pop_good_sel11)
  
  pop_good_sel12 <- selectInd(pop_good12, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good13 <- randCross(pop_good_sel12, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good13 <- setPheno(pop_good13, h2 = 0.8, simParam = SP)
  ddm110var14[i,] <- genicVarA(pop_good13)
  ddm110gv12[i,] <- gv(pop_good_sel12)
  
  pop_good_sel13 <- selectInd(pop_good13, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good14 <- randCross(pop_good_sel13, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good14 <- setPheno(pop_good14, h2 = 0.8, simParam = SP)
  ddm110var15[i,] <- genicVarA(pop_good14)
  ddm110gv13[i,] <- gv(pop_good_sel13)
  
  pop_good_sel14 <- selectInd(pop_good14, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good15 <- randCross(pop_good_sel14, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good15 <- setPheno(pop_good15, h2 = 0.8, simParam = SP)
  ddm110var16[i,] <- genicVarA(pop_good15)
  ddm110gv14[i,] <- gv(pop_good_sel14)
  
  pop_good_sel15 <- selectInd(pop_good15, nInd = 10, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good16 <- randCross(pop_good_sel15, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good16 <- setPheno(pop_good16, h2 = 0.8, simParam = SP)
  ddm110var17[i,] <- genicVarA(pop_good16)
  ddm110gv15[i,] <- gv(pop_good_sel15)
}


 
wt2_gain <- c(wt2gv2[,1:2] - wt2gv1[,1:2], wt2gv3[,1:2] - wt2gv1[,1:2], wt2gv4[,1:2] - wt2gv1[,1:2], wt2gv5[,1:2] - wt2gv1[,1:2], wt2gv6[,1:2] - wt2gv1[,1:2],
             wt2gv7[,1:2] - wt2gv1[,1:2], wt2gv8[,1:2] - wt2gv1[,1:2], wt2gv9[,1:2] - wt2gv1[,1:2], wt2gv10[,1:2] - wt2gv1[,1:2], wt2gv11[,1:2] - wt2gv1[,1:2],
             wt2gv12[,1:2] - wt2gv1[,1:2], wt2gv13[,1:2] - wt2gv1[,1:2], wt2gv14[,1:2] - wt2gv1[,1:2], wt2gv15[,1:2] - wt2gv1[,1:2])

wt5_gain <- c(wt5gv2[,1:5] - wt5gv1[,1:5], wt5gv3[,1:5] - wt5gv1[,1:5], wt5gv4[,1:5] - wt5gv1[,1:5], wt5gv5[,1:5] - wt5gv1[,1:5], wt5gv6[,1:5] - wt5gv1[,1:5],
             wt5gv7[,1:5] - wt5gv1[,1:5], wt5gv8[,1:5] - wt5gv1[,1:5], wt5gv9[,1:5] - wt5gv1[,1:5], wt5gv10[,1:5] - wt5gv1[,1:5], wt5gv11[,1:5] - wt5gv1[,1:5],
             wt5gv12[,1:5] - wt5gv1[,1:5], wt5gv13[,1:5] - wt5gv1[,1:5], wt5gv14[,1:5] - wt5gv1[,1:5], wt5gv15[,1:5] - wt5gv1[,1:5])

wt10_gain <- c(wt10gv2[,1:10] - wt10gv1[,1:10], wt10gv3[,1:10] - wt10gv1[,1:10], wt10gv4[,1:10] - wt10gv1[,1:10], wt10gv5[,1:10] - wt10gv1[,1:10], wt10gv6[,1:10] - wt10gv1[,1:10],
             wt10gv7[,1:10] - wt10gv1[,1:10], wt10gv8[,1:10] - wt10gv1[,1:10], wt10gv9[,1:10] - wt10gv1[,1:10], wt10gv10[,1:10] - wt10gv1[,1:10], wt10gv11[,1:10] - wt10gv1[,1:10],
             wt10gv12[,1:10] - wt10gv1[,1:10], wt10gv13[,1:10] - wt10gv1[,1:10], wt10gv14[,1:10] - wt10gv1[,1:10], wt10gv15[,1:10] - wt10gv1[,1:10])

ddm12_gain <- c(ddm12gv2[,1:2]-ddm12gv1[,1:2], ddm12gv3[,1:2]-ddm12gv1[,1:2], ddm12gv4[,1:2]-ddm12gv1[,1:2], ddm12gv5[,1:2]-ddm12gv1[,1:2], ddm12gv6[,1:2]-ddm12gv1[,1:2],
               ddm12gv7[,1:2]-ddm12gv1[,1:2], ddm12gv8[,1:2]-ddm12gv1[,1:2], ddm12gv9[,1:2]-ddm12gv1[,1:2], ddm12gv10[,1:2]-ddm12gv1[,1:2],
               ddm12gv11[,1:2]-ddm12gv1[,1:2], ddm12gv12[,1:2]-ddm12gv1[,1:2], ddm12gv13[,1:2]-ddm12gv1[,1:2], ddm12gv14[,1:2]-ddm12gv1[,1:2],ddm12gv15[,1:2]-ddm12gv1[,1:2])

ddm15_gain <- c(ddm15gv2[,1:5]-ddm15gv1[,1:5], ddm15gv3[,1:5]-ddm15gv1[,1:5], ddm15gv4[,1:5]-ddm15gv1[,1:5], ddm15gv5[,1:5]-ddm15gv1[,1:5], ddm15gv6[,1:5]-ddm15gv1[,1:5],
               ddm15gv7[,1:5]-ddm15gv1[,1:5], ddm15gv8[,1:5]-ddm15gv1[,1:5], ddm15gv9[,1:5]-ddm15gv1[,1:5], ddm15gv10[,1:5]-ddm15gv1[,1:5],
               ddm15gv11[,1:5]-ddm15gv1[,1:5], ddm15gv12[,1:5]-ddm15gv1[,1:5], ddm15gv13[,1:5]-ddm15gv1[,1:5], ddm15gv14[,1:5]-ddm15gv1[,1:5],ddm15gv15[,1:5]-ddm15gv1[,1:5])

ddm110_gain <- c(ddm110gv2[,1:10]-ddm110gv1[,1:10], ddm110gv3[,1:10]-ddm110gv1[,1:10], ddm110gv4[,1:10]-ddm110gv1[,1:10], ddm110gv5[,1:10]-ddm110gv1[,1:10], ddm110gv6[,1:10]-ddm110gv1[,1:10],
               ddm110gv7[,1:10]-ddm110gv1[,1:10], ddm110gv8[,1:10]-ddm110gv1[,1:10], ddm110gv9[,1:10]-ddm110gv1[,1:10], ddm110gv10[,1:10]-ddm110gv1[,1:10],
               ddm110gv11[,1:10]-ddm110gv1[,1:10], ddm110gv12[,1:10]-ddm110gv1[,1:10], ddm110gv13[,1:10]-ddm110gv1[,1:10], ddm110gv14[,1:10]-ddm110gv1[,1:10],ddm110gv15[,1:10]-ddm110gv1[,1:10])


wt2varall <- c(wt2var3, wt2var4, wt2var5,
              wt2var6,  wt2var7,  wt2var8,  wt2var9,  wt2var10, wt2var11,
              wt2var12, wt2var13,  wt2var14,  wt2var15,  wt2var16,  wt2var17)

wt5varall <- c(wt5var3, wt5var4, wt5var5,
              wt5var6,  wt5var7,  wt5var8,  wt5var9,  wt5var10, wt5var11,
              wt5var12, wt5var13,  wt5var14,  wt5var15,  wt5var16,  wt5var17)

wt10varall <- c(wt10var3, wt10var4, wt10var5,
              wt10var6,  wt10var7,  wt10var8,  wt10var9,  wt10var10, wt10var11,
              wt10var12, wt10var13,  wt10var14,  wt10var15,  wt10var16,  wt10var17)

ddm12varall <- c(ddm12var3, ddm12var4, ddm12var5,
                ddm12var6,  ddm12var7,  ddm12var8,  ddm12var9,  ddm12var10, 
                ddm12var11, ddm12var12, ddm12var13, ddm12var14, ddm12var15, ddm12var16, ddm12var17)

ddm15varall <- c(ddm15var3, ddm15var4, ddm15var5,
                ddm15var6,  ddm15var7,  ddm15var8,  ddm15var9,  ddm15var10, 
                ddm15var11, ddm15var12, ddm15var13, ddm15var14, ddm15var15, ddm15var16, ddm15var17)

ddm110varall <- c(ddm110var3, ddm110var4, ddm110var5,
                ddm110var6,  ddm110var7,  ddm110var8,  ddm110var9,  ddm110var10, 
                ddm110var11, ddm110var12, ddm110var13, ddm110var14, ddm110var15, ddm110var16, ddm110var17)

#combining all genetic gain data for wt & all mutants
gain_inter2 <- cbind(wt2_gain, ddm12_gain)
gain_inter5 <- cbind(wt5_gain, ddm15_gain)
gain_inter10 <- cbind(wt10_gain, ddm110_gain)

gain_inter2 <- as.data.frame(gain_inter2)
gain_inter5 <- as.data.frame(gain_inter5)
gain_inter10 <- as.data.frame(gain_inter10)

gain_inter22 <- as.data.frame(matrix(data = NA, nrow = 5600))
gain_inter25 <- as.data.frame(matrix(data = NA, nrow = 14000))
gain_inter210 <- as.data.frame(matrix(data = NA, nrow = 28000))

gain_inter22$gv <- c(gain_inter2$wt2_gain, gain_inter2$ddm12_gain)
gain_inter25$gv <- c(gain_inter5$wt5_gain, gain_inter5$ddm15_gain)
gain_inter210$gv <- c(gain_inter10$wt10_gain, gain_inter10$ddm110_gain)

gain_inter22$generation <- rep(1:14, size = 2, each = 200)
gain_inter25$generation <- rep(1:14, size = 2, each = 500)
gain_inter210$generation <- rep(1:14, size = 2, each = 1000)

gain_inter22$gen_map <- rep(c("wt 2%", "ddm1 2%"), size = 2, each = 2800)
gain_inter25$gen_map <- rep(c("wt 5%", "ddm1 5%"), size = 2, each = 7000)
gain_inter210$gen_map <- rep(c("wt 10%", "ddm1 10%"), size = 2, each = 14000)

gain_sel <- rbind(gain_inter210, gain_inter25, gain_inter22)
lastgen <- gain_sel[which(gain_sel$generation == c(1,3,5,7,9,11,13,14)),]
lastgen$gen_map <- factor(lastgen$gen_map,
                levels = paste(c("wt 2%","wt 5%","wt 10%","ddm1 2%", "ddm1 5%", "ddm1 10%")))
#firstgen <- gain_inter2[which(gain_inter2$generation == 1),]

gain <- ggplot(lastgen, aes(x=as.factor(generation), y=gv, fill=as.factor(gen_map))) + 
  geom_boxplot() + theme_bw() + xlab("Generations since generation 1") + ylab("Genetic Gain") +
  scale_fill_manual(values=colors, name = "Genetic Map & Sel Int") + ggtitle("Genetic Gain over 15 generations of Recurrent Selection")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12), legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=12), plot.title = element_text(size=12), legend.position = c(0.2,0.75), legend.key.size = unit(0.5, "lines")) +
  guides(color = guide_legend(override.aes = list(size = 2)))

colors = c("wt 5%" = "#EE6A50", "wt 2%" = "#8B3E2F", "wt 10%" = "#FF7F50", "ddm1 2%" = "#556B2F",
           "ddm1 5%" = "#6E8B3D", "ddm1 10%" = "#CAFF70")

allvar <- cbind(wt2varall, wt5varall, wt10varall, ddm12varall, ddm15varall, ddm110varall)
allvar <- as.data.frame(allvar)
allvar2 <- as.data.frame(matrix(data = NA, nrow = 9000))
allvar2$var <- c(allvar$wt2varall, allvar$wt5varall, allvar$wt10varall, allvar$ddm12varall, allvar$ddm15varall, allvar$ddm110varall)
allvar2$gen <- rep(1:15, size = 6, each = 100)
allvar2$gen_map <- rep(c("wt 2%", "wt 5%", "wt 10%", "ddm1 2%", "ddm1 5%", "ddm1 10%"), size = 6, each = 1500)
lastgenvar <- allvar2[which(allvar2$gen == c(1,3,5,7,9,11,13,15)),]
lastgenvar$gen_map <- factor(lastgenvar$gen_map,
                          levels = paste(c("wt 2%","wt 5%","wt 10%","ddm1 2%", "ddm1 5%", "ddm1 10%")))



var <- ggplot(lastgenvar, aes(x=as.factor(gen), y=var, fill=gen_map)) + 
  geom_boxplot() + theme_bw() + xlab("Generation") + ylab("Additive Genetic Variance") + 
  scale_fill_manual(values=colors, name = "Genetic Map & Sel Int") + ggtitle("Additive Genetic Variance over 15 generations of Recurrent Selection") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12), legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=12), plot.title = element_text(size=12), legend.position = c(0.8,0.75), legend.key.size = unit(0.5, "lines"))

figure <- ggarrange(gain, var,
                    labels = c("A", "B"),
                    ncol = 1, nrow = 2)
figure
