library(AlphaSimR)
set.seed(420)
final_map<-readRDS("japonica_final_map.RData")
real_centromere <-readRDS("japonica_centromeres.RData")

snp_pos<-function(lociPositions,finalpos){
  snps <-c()
  row.names(finalpos) <- NULL
  for(k in lociPositions){
    index<-finalpos[which.min(abs(k-finalpos$`SNP Start`)),]
    row<-as.integer(rownames(index))
    snps<-append(snps,row)
  }
  print(snps)
}

chr1_lociPos<-readRDS(file="chr1_locipos")
chr1_geneProb<-readRDS(file="chr1_geneProb")
chr1_lociPositions<-sample(chr1_lociPos,size=30,replace=FALSE,prob=chr1_geneProb)
chr1_finalPos<-readRDS("jap_chr1_finalpos.RData")
chr1_lociPositions2<-snp_pos(chr1_lociPositions,chr1_finalPos)
chr1_lociPositions2 <- chr1_lociPositions2[!duplicated(chr1_lociPositions2)]

chr2_lociPos<-readRDS(file="chr2_locipos")
chr2_geneProb<-readRDS(file="chr2_geneProb")
chr2_lociPositions<-sample(chr2_lociPos,size=25,replace=FALSE,prob=chr2_geneProb)
chr2_finalPos<-readRDS("jap_chr2_finalpos.RData")
chr2_lociPositions2<-snp_pos(chr2_lociPositions,chr2_finalPos)
chr2_lociPositions2 <- chr2_lociPositions2[!duplicated(chr2_lociPositions2)]

chr3_lociPos<-readRDS(file="chr3_locipos")
chr3_geneProb<-readRDS(file="chr3_geneProb")
chr3_lociPositions<-sample(chr3_lociPos,size=25,replace=FALSE,prob=chr3_geneProb)
chr3_finalPos<-readRDS("jap_chr3_finalpos.RData")
chr3_lociPositions2<-snp_pos(chr3_lociPositions,chr3_finalPos)
chr3_lociPositions2 <- chr3_lociPositions2[!duplicated(chr3_lociPositions2)]

chr4_lociPos<-readRDS(file="chr4_locipos")
chr4_geneProb<-readRDS(file="chr4_geneProb")
chr4_lociPositions<-sample(chr4_lociPos,size=25,replace=FALSE,prob=chr4_geneProb)
chr4_finalPos<-readRDS("jap_chr4_finalpos.RData")
chr4_lociPositions2<-snp_pos(chr4_lociPositions,chr4_finalPos)
chr4_lociPositions2 <- chr4_lociPositions2[!duplicated(chr4_lociPositions2)]

chr5_lociPos<-readRDS(file="chr5_locipos")
chr5_geneProb<-readRDS(file="chr5_geneProb")
chr5_lociPositions<-sample(chr5_lociPos,size=22,replace=FALSE,prob=chr5_geneProb)
chr5_finalPos<-readRDS("jap_chr5_finalpos.RData")
chr5_lociPositions2<-snp_pos(chr5_lociPositions,chr5_finalPos)
chr5_lociPositions2 <- chr5_lociPositions2[!duplicated(chr5_lociPositions2)]

chr6_lociPos<-readRDS(file="chr6_locipos")
chr6_geneProb<-readRDS(file="chr6_geneProb")
chr6_lociPositions<-sample(chr6_lociPos,size=20,replace=FALSE,prob=chr6_geneProb)
chr6_finalPos<-readRDS("jap_chr6_finalpos.RData")
chr6_lociPositions2<-snp_pos(chr6_lociPositions,chr6_finalPos)
chr6_lociPositions2 <- chr6_lociPositions2[!duplicated(chr6_lociPositions2)]

chr7_lociPos<-readRDS(file="chr7_locipos")
chr7_geneProb<-readRDS(file="chr7_geneProb")
chr7_lociPositions<-sample(chr7_lociPos,size=23,replace=FALSE,prob=chr7_geneProb)
chr7_finalPos<-readRDS("jap_chr7_finalpos.RData")
chr7_lociPositions2<-snp_pos(chr7_lociPositions,chr7_finalPos)
chr7_lociPositions2 <- chr7_lociPositions2[!duplicated(chr7_lociPositions2)]

chr8_lociPos<-readRDS(file="chr8_locipos")
chr8_geneProb<-readRDS(file="chr8_geneProb")
chr8_lociPositions<-sample(chr8_lociPos,size=18,replace=FALSE,prob=chr8_geneProb)
chr8_finalPos<-readRDS("jap_chr8_finalpos.RData")
chr8_lociPositions2<-snp_pos(chr8_lociPositions,chr8_finalPos)
chr8_lociPositions2 <- chr8_lociPositions2[!duplicated(chr8_lociPositions2)]

chr9_lociPos<-readRDS(file="chr9_locipos")
chr9_geneProb<-readRDS(file="chr9_geneProb")
chr9_lociPositions<-sample(chr9_lociPos,size=15,replace=FALSE,prob=chr9_geneProb)
chr9_finalPos<-readRDS("jap_chr9_finalpos.RData")
chr9_lociPositions2<-snp_pos(chr9_lociPositions,chr9_finalPos)
chr9_lociPositions2 <- chr9_lociPositions2[!duplicated(chr9_lociPositions2)]

chr10_lociPos<-readRDS(file="chr10_locipos")
chr10_geneProb<-readRDS(file="chr10_geneProb")
chr10_lociPositions<-sample(chr10_lociPos,size=18,replace=FALSE,prob=chr10_geneProb)
chr10_finalPos<-readRDS("jap_chr10_finalpos.RData")
chr10_lociPositions2<-snp_pos(chr10_lociPositions,chr10_finalPos)
chr10_lociPositions2 <- chr10_lociPositions2[!duplicated(chr10_lociPositions2)]

chr11_lociPos<-readRDS(file="chr11_locipos")
chr11_geneProb<-readRDS(file="chr11_geneProb")
chr11_lociPositions<-sample(chr11_lociPos,size=20,replace=FALSE,prob=chr11_geneProb)
chr11_finalPos<-readRDS("jap_chr11_finalpos.RData")
chr11_lociPositions2<-snp_pos(chr11_lociPositions,chr11_finalPos)
chr11_lociPositions2 <- chr11_lociPositions2[!duplicated(chr11_lociPositions2)]

chr12_lociPos<-readRDS(file="chr12_locipos")
chr12_geneProb<-readRDS(file="chr12_geneProb")
chr12_lociPositions<-sample(chr12_lociPos,size=20,replace=FALSE,prob=chr12_geneProb)
chr12_finalPos<-readRDS("jap_chr12_finalpos.RData")
chr12_lociPositions2<-snp_pos(chr12_lociPositions,chr12_finalPos)
chr12_lociPositions2 <- chr12_lociPositions2[!duplicated(chr12_lociPositions2)]


#generating repulsion ONLY linkages between QTL
addEff_repulsion <- runif(209, min = 0, max = 0.1)
addEff_repulsion <- addEff_repulsion * c(-1,1)
sum(addEff_repulsion)

#generating mix of repulsion & coupling linkages between QTL
addEff_mix <- runif(240, min = -0.1, max = 0.1)
sum(addEff_mix)

#burn-in to create LD amongst SNPs
burn_in_pop <- vector(mode = "list", length = 100)
segSites<-readRDS("japonica_num_SNP.RData")
for(i in 1:100){
  founderPop <- quickHaplo(nInd = 200, nChr = 10, inbred = TRUE, ploidy = 2L, segSites = segSites)
  founderPop@genMap <- final_map
  founderPop@centromere <- real_centromere
  SP = SimParam$new(founderPop)
  SP$setTrackRec(TRUE)
  SP$p = 0.15
  #polygenic trait with QTLs in accurate gene space
  trait_yield <- new("TraitA", lociPerChr= c(length(chr1_lociPositions), length(chr2_lociPositions), length(chr3_lociPositions), length(chr4_lociPositions), length(chr5_lociPositions), length(chr6_lociPositions), length(chr7_lociPositions), length(chr8_lociPositions), length(chr9_lociPositions), length(chr10_lociPositions), length(chr11_lociPositions), length(chr12_lociPositions)),
                     lociLoc = c(chr1_lociPositions2,chr2_lociPositions2,chr3_lociPositions2,chr4_lociPositions2,chr5_lociPositions2,chr6_lociPositions2,chr7_lociPositions2,chr8_lociPositions2,chr9_lociPositions2,chr10_lociPositions2,chr11_lociPositions2,chr12_lociPositions2),
                     nLoci=sum(lociPerChr),addEff = addEff_mix, intercept = 0.1)
  SP$manAddTrait(trait_yield)
  pop_good <- newPop(founderPop, simParam = SP)
  pop_good <- randCross(pop_good, nCrosses= 10, nProgeny=10, simParam = SP)
  pop_good <- setPheno(pop_good, h2 = 0.8, simParam = SP)
  
  pop_good_sel0 <- selectInd(pop_good, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good1 <- self(pop_good, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good1 <- setPheno(pop_good1, h2 = 0.8, simParam = SP)
  
  pop_good_sel <- selectInd(pop_good1, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good2 <- self(pop_good_sel, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good2 <- setPheno(pop_good2, h2 = 0.8, simParam = SP)
  
  pop_good_sel2 <- selectInd(pop_good2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good3 <- self(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good3 <- setPheno(pop_good3, h2 = 0.8, simParam = SP)
  
  pop_good_sel3 <- selectInd(pop_good3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good4 <- self(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good4 <- setPheno(pop_good4, h2 = 0.8, simParam = SP)
  
  pop_good_sel4 <- selectInd(pop_good4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good5 <- self(pop_good_sel4, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good5 <- setPheno(pop_good5, h2 = 0.8, simParam = SP)
  
  pop_good_sel5 <- selectInd(pop_good5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good6 <- self(pop_good_sel5, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good6 <- setPheno(pop_good6, h2 = 0.8, simParam = SP)
  
  pop_good_sel6 <- selectInd(pop_good6, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good7 <- self(pop_good_sel6, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good7 <- setPheno(pop_good7, h2 = 0.8, simParam = SP)
  
  pop_good_sel7 <- selectInd(pop_good7, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good8 <- self(pop_good_sel7, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good8 <- setPheno(pop_good8, h2 = 0.8, simParam = SP)
  
  pop_good_sel8 <- selectInd(pop_good8, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good9 <- self(pop_good_sel8, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good9 <- setPheno(pop_good9, h2 = 0.8, simParam = SP)
  
  pop_good_sel9 <- selectInd(pop_good9, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good10 <- self(pop_good_sel9, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good10 <- setPheno(pop_good10, h2 = 0.8, simParam = SP)
  
  pop_good_sel10 <- selectInd(pop_good10, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good11 <- self(pop_good_sel10, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good11 <- setPheno(pop_good11, h2 = 0.8, simParam = SP)
  
  burn_in_pop[[i]] <- selectInd(pop_good10, nInd = 5, use = "gv", trait = 1, selectop = TRUE, returnPop = TRUE, simParam = SP)
}

#combine all 100 iterations together to make one population with 200 individuals
burn_in_pop <- mergePops(burn_in_pop)

##Starting intermating schemes for different genetic maps
#wild-type map
wtvar1 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar2 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar3 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar4 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar5 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar6 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar7 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar8 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar9 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar10 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar11 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar12 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar13 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar14 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar15 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar16 <- matrix(data = NA, ncol = 1, nrow = 100)
wtvar17 <- matrix(data = NA, ncol = 1, nrow = 100)
wtgv1 <- matrix(data = NA, ncol = 5, nrow = 100)
wtgv2 <- matrix(data = NA, ncol = 5, nrow = 100)
wtgv3 <- matrix(data = NA, ncol = 5, nrow = 100)
wtgv4 <- matrix(data = NA, ncol = 5, nrow = 100)
wtgv5 <- matrix(data = NA, ncol = 5, nrow = 100)
wtgv6 <- matrix(data = NA, ncol = 5, nrow = 100)
wtgv7 <- matrix(data = NA, ncol = 5, nrow = 100)
wtgv8 <- matrix(data = NA, ncol = 5, nrow = 100)
wtgv9 <- matrix(data = NA, ncol = 5, nrow = 100)
wtgv10 <- matrix(data = NA, ncol = 5, nrow = 100)
wtgv11 <- matrix(data = NA, ncol = 5, nrow = 100)
wtgv12 <- matrix(data = NA, ncol = 5, nrow = 100)
wtgv13 <- matrix(data = NA, ncol = 5, nrow = 100)
wtgv14 <- matrix(data = NA, ncol = 5, nrow = 100)
wtgv15 <- matrix(data = NA, ncol = 5, nrow = 100)
for(i in 1:100){
  SP$switchGenMap(final_map, centromere = real_centromere)
  wtvar1[i,] <- genicVarA(burn_in_pop)
  
  pop_good1 <- self(burn_in_pop, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good1 <- setPheno(pop_good1, h2 = 0.8, simParam = SP)
  wtvar2[i,] <- genicVarA(pop_good1)
  
  pop_good_sel <- selectInd(pop_good1, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good2 <- self(pop_good_sel, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good2 <- setPheno(pop_good2, h2 = 0.8, simParam = SP)
  wtvar3[i,] <- genicVarA(pop_good2)
  wtgv1[i,] <- gv(pop_good_sel)
  
  pop_good_sel2 <- selectInd(pop_good2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good3 <- self(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good3 <- setPheno(pop_good3, h2 = 0.8, simParam = SP)
  wtvar4[i,] <- genicVarA(pop_good3)
  wtgv2[i,] <- gv(pop_good_sel2)
  
  pop_good_sel3 <- selectInd(pop_good3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good4 <- self(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good4 <- setPheno(pop_good4, h2 = 0.8, simParam = SP)
  wtvar5[i,] <- genicVarA(pop_good4)
  wtgv3[i,] <- gv(pop_good_sel3)
  
  pop_good_sel4 <- selectInd(pop_good4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good5 <- self(pop_good_sel4, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good5 <- setPheno(pop_good5, h2 = 0.8, simParam = SP)
  wtvar6[i,] <- genicVarA(pop_good5)
  wtgv4[i,] <- gv(pop_good_sel4)
  
  pop_good_sel5 <- selectInd(pop_good5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good6 <- self(pop_good_sel5, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good6 <- setPheno(pop_good6, h2 = 0.8, simParam = SP)
  wtvar7[i,] <- genicVarA(pop_good6)
  wtgv5[i,] <- gv(pop_good_sel5)
  
  pop_good_sel6 <- selectInd(pop_good6, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good7 <- self(pop_good_sel6, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good7 <- setPheno(pop_good7, h2 = 0.8, simParam = SP)
  wtvar8[i,] <- genicVarA(pop_good7)
  wtgv6[i,] <- gv(pop_good_sel6)
  
  pop_good_sel7 <- selectInd(pop_good7, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good8 <- self(pop_good_sel7, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good8 <- setPheno(pop_good8, h2 = 0.8, simParam = SP)
  wtvar9[i,] <- genicVarA(pop_good8)
  wtgv7[i,] <- gv(pop_good_sel7)
  
  pop_good_sel8 <- selectInd(pop_good8, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good9 <- self(pop_good_sel8, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good9 <- setPheno(pop_good9, h2 = 0.8, simParam = SP)
  wtvar10[i,] <- genicVarA(pop_good9)
  wtgv8[i,] <- gv(pop_good_sel8)
  
  pop_good_sel9 <- selectInd(pop_good9, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good10 <- self(pop_good_sel9, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good10 <- setPheno(pop_good10, h2 = 0.8, simParam = SP)
  wtvar11[i,] <- genicVarA(pop_good10)
  wtgv9[i,] <- gv(pop_good_sel9)
  
  pop_good_sel10 <- selectInd(pop_good10, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good11 <- self(pop_good_sel10, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good11 <- setPheno(pop_good11, h2 = 0.8, simParam = SP)
  wtvar12[i,] <- genicVarA(pop_good11)
  wtgv10[i,] <- gv(pop_good_sel10)
  
  pop_good_sel11 <- selectInd(pop_good11, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good12 <- self(pop_good_sel11, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good12 <- setPheno(pop_good12, h2 = 0.8, simParam = SP)
  wtvar13[i,] <- genicVarA(pop_good12)
  wtgv11[i,] <- gv(pop_good_sel11)
  
  pop_good_sel12 <- selectInd(pop_good12, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good13 <- self(pop_good_sel12, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good13 <- setPheno(pop_good13, h2 = 0.8, simParam = SP)
  wtvar14[i,] <- genicVarA(pop_good13)
  wtgv12[i,] <- gv(pop_good_sel12)
  
  pop_good_sel13 <- selectInd(pop_good13, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good14 <- self(pop_good_sel13, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good14 <- setPheno(pop_good14, h2 = 0.8, simParam = SP)
  wtvar15[i,] <- genicVarA(pop_good14)
  wtgv13[i,] <- gv(pop_good_sel13)
  
  pop_good_sel14 <- selectInd(pop_good14, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good15 <- self(pop_good_sel14, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good15 <- setPheno(pop_good15, h2 = 0.8, simParam = SP)
  wtvar16[i,] <- genicVarA(pop_good15)
  wtgv14[i,] <- gv(pop_good_sel14)
  
  pop_good_sel15 <- selectInd(pop_good15, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good16 <- self(pop_good_sel14, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good16 <- setPheno(pop_good15, h2 = 0.8, simParam = SP)
  wtvar17[i,] <- genicVarA(pop_good16)
  wtgv15[i,] <- gv(pop_good_sel15)
}

#with ddm1 genetic map
ddm1var1 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var2 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var3 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var4 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var5 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var6 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var7 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var8 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var9 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var10 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var11 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var12 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var13 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var14 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var15 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var16 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1var17 <- matrix(data = NA, ncol = 1, nrow = 100)
ddm1gv1 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm1gv2 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm1gv3 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm1gv4 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm1gv5 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm1gv6 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm1gv7 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm1gv8 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm1gv9 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm1gv10 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm1gv11 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm1gv12 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm1gv13 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm1gv14 <- matrix(data = NA, ncol = 5, nrow = 100)
ddm1gv15 <- matrix(data = NA, ncol = 5, nrow = 100)
for(i in 1:100){
  SP$switchGenMap(ddm1_map, centromere = ddm1_centromere)
  ddm1var1[i,] <- genicVarA(burn_in_pop)

  pop_good1 <- self(burn_in_pop, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good1 <- setPheno(pop_good1, h2 = 0.8, simParam = SP)
  ddm1var2[i,] <- genicVarA(pop_good1)
  
  pop_good_sel <- selectInd(pop_good1, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good2 <- self(pop_good_sel, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good2 <- setPheno(pop_good2, h2 = 0.8, simParam = SP)
  ddm1var3[i,] <- genicVarA(pop_good2)
  ddm1gv1[i,] <- gv(pop_good_sel)
  
  pop_good_sel2 <- selectInd(pop_good2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good3 <- self(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good3 <- setPheno(pop_good3, h2 = 0.8, simParam = SP)
  ddm1var4[i,] <- genicVarA(pop_good3)
  ddm1gv2[i,] <- gv(pop_good_sel2)
  
  pop_good_sel3 <- selectInd(pop_good3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good4 <- self(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good4 <- setPheno(pop_good4, h2 = 0.8, simParam = SP)
  ddm1var5[i,] <- genicVarA(pop_good4)
  ddm1gv3[i,] <- gv(pop_good_sel3)
  
  pop_good_sel4 <- selectInd(pop_good4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good5 <- self(pop_good_sel4, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good5 <- setPheno(pop_good5, h2 = 0.8, simParam = SP)
  ddm1var6[i,] <- genicVarA(pop_good5)
  ddm1gv4[i,] <- gv(pop_good_sel4)
  
  pop_good_sel5 <- selectInd(pop_good5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good6 <- self(pop_good_sel5, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good6 <- setPheno(pop_good6, h2 = 0.8, simParam = SP)
  ddm1var7[i,] <- genicVarA(pop_good6)
  ddm1gv5[i,] <- gv(pop_good_sel5)
  
  pop_good_sel6 <- selectInd(pop_good6, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good7 <- self(pop_good_sel6, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good7 <- setPheno(pop_good7, h2 = 0.8, simParam = SP)
  ddm1var8[i,] <- genicVarA(pop_good7)
  ddm1gv6[i,] <- gv(pop_good_sel6)
  
  pop_good_sel7 <- selectInd(pop_good7, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good8 <- self(pop_good_sel7, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good8 <- setPheno(pop_good8, h2 = 0.8, simParam = SP)
  ddm1var9[i,] <-genicVarA(pop_good8)
  ddm1gv7[i,] <- gv(pop_good_sel7)
  
  pop_good_sel8 <- selectInd(pop_good8, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good9 <- self(pop_good_sel8, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good9 <- setPheno(pop_good9, h2 = 0.8, simParam = SP)
  ddm1var10[i,] <- genicVarA(pop_good9)
  ddm1gv8[i,] <- gv(pop_good_sel8)
  
  pop_good_sel9 <- selectInd(pop_good9, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good10 <- self(pop_good_sel9, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good10 <- setPheno(pop_good10, h2 = 0.8, simParam = SP)
  ddm1var11[i,] <- genicVarA(pop_good10)
  ddm1gv9[i,] <- gv(pop_good_sel9)
  
  pop_good_sel10 <- selectInd(pop_good10, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good11 <- self(pop_good_sel10, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good11 <- setPheno(pop_good11, h2 = 0.8, simParam = SP)
  ddm1var12[i,] <- genicVarA(pop_good11)
  ddm1gv10[i,] <- gv(pop_good_sel10)
  
  pop_good_sel11 <- selectInd(pop_good11, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good12 <- self(pop_good_sel11, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good12 <- setPheno(pop_good12, h2 = 0.8, simParam = SP)
  ddm1var13[i,] <- genicVarA(pop_good12)
  ddm1gv11[i,] <- gv(pop_good_sel11)
  
  pop_good_sel12 <- selectInd(pop_good12, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good13 <- self(pop_good_sel12, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good13 <- setPheno(pop_good13, h2 = 0.8, simParam = SP)
  ddm1var14[i,] <- genicVarA(pop_good13)
  ddm1gv12[i,] <- gv(pop_good_sel12)
  
  pop_good_sel13 <- selectInd(pop_good13, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good14 <- self(pop_good_sel13, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good14 <- setPheno(pop_good14, h2 = 0.8, simParam = SP)
  ddm1var15[i,] <- genicVarA(pop_good14)
  ddm1gv13[i,] <- gv(pop_good_sel13)
  
  pop_good_sel14 <- selectInd(pop_good14, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good15 <- self(pop_good_sel14, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good15 <- setPheno(pop_good15, h2 = 0.8, simParam = SP)
  ddm1var16[i,] <- genicVarA(pop_good15)
  ddm1gv14[i,] <- gv(pop_good_sel14)
  
  pop_good_sel15 <- selectInd(pop_good15, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good16 <- self(pop_good_sel15, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good16 <- setPheno(pop_good16, h2 = 0.8, simParam = SP)
  ddm1var17[i,] <- genicVarA(pop_good16)
  ddm1gv15[i,] <- gv(pop_good_sel15)
}

#zmet2
zmet2var1 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var2 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var3 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var4 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var5 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var6 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var7 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var8 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var9 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var10 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var11 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var12 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var13 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var14 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var15 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var16 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2var17 <- matrix(data = NA, ncol = 1, nrow = 100)
zmet2gv1 <- matrix(data = NA, ncol = 5, nrow = 100)
zmet2gv2 <- matrix(data = NA, ncol = 5, nrow = 100)
zmet2gv3 <- matrix(data = NA, ncol = 5, nrow = 100)
zmet2gv4 <- matrix(data = NA, ncol = 5, nrow = 100)
zmet2gv5 <- matrix(data = NA, ncol = 5, nrow = 100)
zmet2gv6 <- matrix(data = NA, ncol = 5, nrow = 100)
zmet2gv7 <- matrix(data = NA, ncol = 5, nrow = 100)
zmet2gv8 <- matrix(data = NA, ncol = 5, nrow = 100)
zmet2gv9 <- matrix(data = NA, ncol = 5, nrow = 100)
zmet2gv10 <- matrix(data = NA, ncol = 5, nrow = 100)
zmet2gv11 <- matrix(data = NA, ncol = 5, nrow = 100)
zmet2gv12 <- matrix(data = NA, ncol = 5, nrow = 100)
zmet2gv13 <- matrix(data = NA, ncol = 5, nrow = 100)
zmet2gv14 <- matrix(data = NA, ncol = 5, nrow = 100)
zmet2gv15 <- matrix(data = NA, ncol = 5, nrow = 100)
for(i in 1:100){
  SP$switchGenMap(zmet2_map, centromere = zmet2_centromere)
  zmet2var1[i,] <- genicVarA(burn_in_pop)
  
  pop_good1 <- self(burn_in_pop, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good1 <- setPheno(pop_good1, h2 = 0.8, simParam = SP)
  zmet2var2[i,] <- genicVarA(pop_good1)
  
  pop_good_sel <- selectInd(pop_good1, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good2 <- self(pop_good_sel, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good2 <- setPheno(pop_good2, h2 = 0.8, simParam = SP)
  zmet2var3[i,] <- genicVarA(pop_good2)
  zmet2gv1[i,] <- gv(pop_good_sel)
  
  pop_good_sel2 <- selectInd(pop_good2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good3 <- self(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good3 <- setPheno(pop_good3, h2 = 0.8, simParam = SP)
  zmet2var4[i,] <- genicVarA(pop_good3)
  zmet2gv2[i,] <- gv(pop_good_sel2)
  
  pop_good_sel3 <- selectInd(pop_good3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good4 <- self(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good4 <- setPheno(pop_good4, h2 = 0.8, simParam = SP)
  zmet2var5[i,] <- genicVarA(pop_good4)
  zmet2gv3[i,] <- gv(pop_good_sel3)
  
  pop_good_sel4 <- selectInd(pop_good4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good5 <- self(pop_good_sel4, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good5 <- setPheno(pop_good5, h2 = 0.8, simParam = SP)
  zmet2var6[i,] <- genicVarA(pop_good5)
  zmet2gv4[i,] <- gv(pop_good_sel4)
  
  pop_good_sel5 <- selectInd(pop_good5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good6 <- self(pop_good_sel5, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good6 <- setPheno(pop_good6, h2 = 0.8, simParam = SP)
  zmet2var7[i,] <- genicVarA(pop_good6)
  zmet2gv5[i,] <- gv(pop_good_sel5)
  
  pop_good_sel6 <- selectInd(pop_good6, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good7 <- self(pop_good_sel6, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good7 <- setPheno(pop_good7, h2 = 0.8, simParam = SP)
  zmet2var8[i,] <- genicVarA(pop_good7)
  zmet2gv6[i,] <- gv(pop_good_sel6)
  
  pop_good_sel7 <- selectInd(pop_good7, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good8 <- self(pop_good_sel7, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good8 <- setPheno(pop_good8, h2 = 0.8, simParam = SP)
  zmet2var9[i,] <- genicVarA(pop_good8)
  zmet2gv7[i,] <- gv(pop_good_sel7)
  
  pop_good_sel8 <- selectInd(pop_good8, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good9 <- self(pop_good_sel8, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good9 <- setPheno(pop_good9, h2 = 0.8, simParam = SP)
  zmet2var10[i,] <- genicVarA(pop_good9)
  zmet2gv8[i,] <- gv(pop_good_sel8)
  
  pop_good_sel9 <- selectInd(pop_good9, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good10 <- self(pop_good_sel9, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good10 <- setPheno(pop_good10, h2 = 0.8, simParam = SP)
  zmet2var11[i,] <- genicVarA(pop_good10)
  zmet2gv9[i,] <- gv(pop_good_sel9)
  
  pop_good_sel10 <- selectInd(pop_good10, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good11 <- self(pop_good_sel10, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good11 <- setPheno(pop_good11, h2 = 0.8, simParam = SP)
  zmet2var12[i,] <- genicVarA(pop_good11)
  zmet2gv10[i,] <- gv(pop_good_sel10)
  
  pop_good_sel11 <- selectInd(pop_good11, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good12 <- self(pop_good_sel11, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good12 <- setPheno(pop_good12, h2 = 0.8, simParam = SP)
  zmet2var13[i,] <- genicVarA(pop_good12)
  zmet2gv11[i,] <- gv(pop_good_sel11)
  
  pop_good_sel12 <- selectInd(pop_good12, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good13 <- self(pop_good_sel12, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good13 <- setPheno(pop_good13, h2 = 0.8, simParam = SP)
  zmet2var14[i,] <- genicVarA(pop_good13)
  zmet2gv12[i,] <- gv(pop_good_sel12)
  
  pop_good_sel13 <- selectInd(pop_good13, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good14 <- self(pop_good_sel13, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good14 <- setPheno(pop_good14, h2 = 0.8, simParam = SP)
  zmet2var15[i,] <- genicVarA(pop_good14)
  zmet2gv13[i,] <- gv(pop_good_sel13)
  
  pop_good_sel14 <- selectInd(pop_good14, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good15 <- self(pop_good_sel14, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good15 <- setPheno(pop_good15, h2 = 0.8, simParam = SP)
  zmet2var16[i,] <- genicVarA(pop_good15)
  zmet2gv14[i,] <- gv(pop_good_sel14)
  
  pop_good_sel15 <- selectInd(pop_good15, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good16 <- self(pop_good_sel15, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good16 <- setPheno(pop_good16, h2 = 0.8, simParam = SP)
  zmet2var17[i,] <- genicVarA(pop_good16)
  zmet2gv15[i,] <- gv(pop_good_sel15)
}

#recq4
recq4var1 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var2 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var3 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var4 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var5 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var6 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var7 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var8 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var9 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var10 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var11 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var12 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var13 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var14 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var15 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var16 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4var17 <- matrix(data = NA, ncol = 1, nrow = 100)
recq4gv1 <- matrix(data = NA, ncol = 5, nrow = 100)
recq4gv2 <- matrix(data = NA, ncol = 5, nrow = 100)
recq4gv3 <- matrix(data = NA, ncol = 5, nrow = 100)
recq4gv4 <- matrix(data = NA, ncol = 5, nrow = 100)
recq4gv5 <- matrix(data = NA, ncol = 5, nrow = 100)
recq4gv6 <- matrix(data = NA, ncol = 5, nrow = 100)
recq4gv7 <- matrix(data = NA, ncol = 5, nrow = 100)
recq4gv8 <- matrix(data = NA, ncol = 5, nrow = 100)
recq4gv9 <- matrix(data = NA, ncol = 5, nrow = 100)
recq4gv10 <- matrix(data = NA, ncol = 5, nrow = 100)
recq4gv11 <- matrix(data = NA, ncol = 5, nrow = 100)
recq4gv12 <- matrix(data = NA, ncol = 5, nrow = 100)
recq4gv13 <- matrix(data = NA, ncol = 5, nrow = 100)
recq4gv14 <- matrix(data = NA, ncol = 5, nrow = 100)
recq4gv15 <- matrix(data = NA, ncol = 5, nrow = 100)
for(i in 1:100){
  SP$switchGenMap(recq4_map, centromere = recq4_centromere)
  recq4var1[i,] <- genicVarA(burn_in_pop)
  
  pop_good1 <- self(burn_in_pop, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good1 <- setPheno(pop_good1, h2 = 0.8, simParam = SP)
  recq4var2[i,] <- genicVarA(pop_good1)
  
  pop_good_sel <- selectInd(pop_good1, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good2 <- self(pop_good_sel, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good2 <- setPheno(pop_good2, h2 = 0.8, simParam = SP)
  recq4var3[i,] <- genicVarA(pop_good2)
  recq4gv1[i,] <- gv(pop_good_sel)
  
  pop_good_sel2 <- selectInd(pop_good2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good3 <- self(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good3 <- setPheno(pop_good3, h2 = 0.8, simParam = SP)
  recq4var4[i,] <- genicVarA(pop_good3)
  recq4gv2[i,] <- gv(pop_good_sel2)
  
  pop_good_sel3 <- selectInd(pop_good3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good4 <- self(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good4 <- setPheno(pop_good4, h2 = 0.8, simParam = SP)
  recq4var5[i,] <- genicVarA(pop_good4)
  recq4gv3[i,] <- gv(pop_good_sel3)
  
  pop_good_sel4 <- selectInd(pop_good4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good5 <- self(pop_good_sel4, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good5 <- setPheno(pop_good5, h2 = 0.8, simParam = SP)
  recq4var6[i,] <- genicVarA(pop_good5)
  recq4gv4[i,] <- gv(pop_good_sel4)
  
  pop_good_sel5 <- selectInd(pop_good5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good6 <- self(pop_good_sel5, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good6 <- setPheno(pop_good6, h2 = 0.8, simParam = SP)
  recq4var7[i,] <- genicVarA(pop_good6)
  recq4gv5[i,] <- gv(pop_good_sel5)
  
  pop_good_sel6 <- selectInd(pop_good6, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good7 <- self(pop_good_sel6, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good7 <- setPheno(pop_good7, h2 = 0.8, simParam = SP)
  recq4var8[i,] <- genicVarA(pop_good7)
  recq4gv6[i,] <- gv(pop_good_sel6)
  
  pop_good_sel7 <- selectInd(pop_good7, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good8 <- self(pop_good_sel7, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good8 <- setPheno(pop_good8, h2 = 0.8, simParam = SP)
  recq4var9[i,] <- genicVarA(pop_good8)
  recq4gv7[i,] <- gv(pop_good_sel7)
  
  pop_good_sel8 <- selectInd(pop_good8, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good9 <- self(pop_good_sel8, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good9 <- setPheno(pop_good9, h2 = 0.8, simParam = SP)
  recq4var10[i,] <- genicVarA(pop_good9)
  recq4gv8[i,] <- gv(pop_good_sel8)
  
  pop_good_sel9 <- selectInd(pop_good9, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good10 <- self(pop_good_sel9, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good10 <- setPheno(pop_good10, h2 = 0.8, simParam = SP)
  recq4var11[i,] <- genicVarA(pop_good10)
  recq4gv9[i,] <- gv(pop_good_sel9)
  
  pop_good_sel10 <- selectInd(pop_good10, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good11 <- self(pop_good_sel10, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good11 <- setPheno(pop_good11, h2 = 0.8, simParam = SP)
  recq4var12[i,] <- genicVarA(pop_good11)
  recq4gv10[i,] <- gv(pop_good_sel10)
  
  pop_good_sel11 <- selectInd(pop_good11, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good12 <- self(pop_good_sel11, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good12 <- setPheno(pop_good12, h2 = 0.8, simParam = SP)
  recq4var13[i,] <- genicVarA(pop_good12)
  recq4gv11[i,] <- gv(pop_good_sel11)
  
  pop_good_sel12 <- selectInd(pop_good12, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good13 <- self(pop_good_sel12, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good13 <- setPheno(pop_good13, h2 = 0.8, simParam = SP)
  recq4var14[i,] <- genicVarA(pop_good13)
  recq4gv12[i,] <- gv(pop_good_sel12)
  
  pop_good_sel13 <- selectInd(pop_good13, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good14 <- self(pop_good_sel13, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good14 <- setPheno(pop_good14, h2 = 0.8, simParam = SP)
  recq4var15[i,] <- genicVarA(pop_good14)
  recq4gv13[i,] <- gv(pop_good_sel13)
  
  pop_good_sel14 <- selectInd(pop_good14, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good15 <- self(pop_good_sel14, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good15 <- setPheno(pop_good15, h2 = 0.8, simParam = SP)
  recq4var16[i,] <- genicVarA(pop_good15)
  recq4gv14[i,] <- gv(pop_good_sel14)
  
  pop_good_sel15 <- selectInd(pop_good15, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good16 <- self(pop_good_sel15, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good16 <- setPheno(pop_good16, h2 = 0.8, simParam = SP)
  recq4var17[i,] <- genicVarA(pop_good16)
  recq4gv15[i,] <- gv(pop_good_sel15)
}

#fancm
fancmvar1 <- matrix(data = NA, ncol = 1, nrow = 100)
fancmvar2 <- matrix(data = NA, ncol = 1, nrow = 100)
fancmvar3 <- matrix(data = NA, ncol = 1, nrow = 100)
fancmvar4 <- matrix(data = NA, ncol = 1, nrow = 100)
fancmvar5 <- matrix(data = NA, ncol = 1, nrow = 100)
fancmvar6 <- matrix(data = NA, ncol = 1, nrow = 100)
fancmvar7 <- matrix(data = NA, ncol = 1, nrow = 100)
fancmvar8 <- matrix(data = NA, ncol = 1, nrow = 100)
fancmvar9 <- matrix(data = NA, ncol = 1, nrow = 100)
fancmvar10 <- matrix(data = NA, ncol = 1, nrow = 100)
fancmvar11 <- matrix(data = NA, ncol = 1, nrow = 100)
fancmvar12 <- matrix(data = NA, ncol = 1, nrow = 100)
fancmvar13 <- matrix(data = NA, ncol = 1, nrow = 100)
fancmvar14 <- matrix(data = NA, ncol = 1, nrow = 100)
fancmvar15 <- matrix(data = NA, ncol = 1, nrow = 100)
fancmvar16 <- matrix(data = NA, ncol = 1, nrow = 100)
fancmvar17 <- matrix(data = NA, ncol = 1, nrow = 100)
fancmgv1 <- matrix(data = NA, ncol = 5, nrow = 100)
fancmgv2 <- matrix(data = NA, ncol = 5, nrow = 100)
fancmgv3 <- matrix(data = NA, ncol = 5, nrow = 100)
fancmgv4 <- matrix(data = NA, ncol = 5, nrow = 100)
fancmgv5 <- matrix(data = NA, ncol = 5, nrow = 100)
fancmgv6 <- matrix(data = NA, ncol = 5, nrow = 100)
fancmgv7 <- matrix(data = NA, ncol = 5, nrow = 100)
fancmgv8 <- matrix(data = NA, ncol = 5, nrow = 100)
fancmgv9 <- matrix(data = NA, ncol = 5, nrow = 100)
fancmgv10 <- matrix(data = NA, ncol = 5, nrow = 100)
fancmgv11 <- matrix(data = NA, ncol = 5, nrow = 100)
fancmgv12 <- matrix(data = NA, ncol = 5, nrow = 100)
fancmgv13 <- matrix(data = NA, ncol = 5, nrow = 100)
fancmgv14 <- matrix(data = NA, ncol = 5, nrow = 100)
fancmgv15 <- matrix(data = NA, ncol = 5, nrow = 100)
for(i in 1:100){
  SP$switchGenMap(fancm_map, centromere = fancm_centromere)
  fancmvar1[i,] <- genicVarA(burn_in_pop)
  
  pop_good1 <- self(burn_in_pop, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good1 <- setPheno(pop_good1, h2 = 0.8, simParam = SP)
  fancmvar2[i,] <- genicVarA(pop_good1)
  
  pop_good_sel <- selectInd(pop_good1, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good2 <- self(pop_good_sel, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good2 <- setPheno(pop_good2, h2 = 0.8, simParam = SP)
  fancmvar3[i,] <- genicVarA(pop_good2)
  fancmgv1[i,] <- gv(pop_good_sel)
  
  pop_good_sel2 <- selectInd(pop_good2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good3 <- self(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good3 <- setPheno(pop_good3, h2 = 0.8, simParam = SP)
  fancmvar4[i,] <- genicVarA(pop_good3)
  fancmgv2[i,] <- gv(pop_good_sel2)
  
  pop_good_sel3 <- selectInd(pop_good3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good4 <- self(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good4 <- setPheno(pop_good4, h2 = 0.8, simParam = SP)
  fancmvar5[i,] <- genicVarA(pop_good4)
  fancmgv3[i,] <- gv(pop_good_sel3)
  
  pop_good_sel4 <- selectInd(pop_good4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good5 <- self(pop_good_sel4, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good5 <- setPheno(pop_good5, h2 = 0.8, simParam = SP)
  fancmvar6[i,] <- genicVarA(pop_good5)
  fancmgv4[i,] <- gv(pop_good_sel4)
  
  pop_good_sel5 <- selectInd(pop_good5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good6 <- self(pop_good_sel5, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good6 <- setPheno(pop_good6, h2 = 0.8, simParam = SP)
  fancmvar7[i,] <- genicVarA(pop_good6)
  fancmgv5[i,] <- gv(pop_good_sel5)
  
  pop_good_sel6 <- selectInd(pop_good6, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good7 <- self(pop_good_sel6, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good7 <- setPheno(pop_good7, h2 = 0.8, simParam = SP)
  fancmvar8[i,] <- genicVarA(pop_good7)
  fancmgv6[i,] <- gv(pop_good_sel6)
  
  pop_good_sel7 <- selectInd(pop_good7, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good8 <- self(pop_good_sel7, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good8 <- setPheno(pop_good8, h2 = 0.8, simParam = SP)
  fancmvar9[i,] <- genicVarA(pop_good8)
  fancmgv7[i,] <- gv(pop_good_sel7)
  
  pop_good_sel8 <- selectInd(pop_good8, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good9 <- self(pop_good_sel8, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good9 <- setPheno(pop_good9, h2 = 0.8, simParam = SP)
  fancmvar10[i,] <- genicVarA(pop_good9)
  fancmgv8[i,] <- gv(pop_good_sel8)
  
  pop_good_sel9 <- selectInd(pop_good9, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good10 <- self(pop_good_sel9, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good10 <- setPheno(pop_good10, h2 = 0.8, simParam = SP)
  fancmvar11[i,] <- genicVarA(pop_good10)
  fancmgv9[i,] <- gv(pop_good_sel9)
  
  pop_good_sel10 <- selectInd(pop_good10, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good11 <- self(pop_good_sel10, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good11 <- setPheno(pop_good11, h2 = 0.8, simParam = SP)
  fancmvar12[i,] <- genicVarA(pop_good11)
  fancmgv10[i,] <- gv(pop_good_sel10)
  
  pop_good_sel11 <- selectInd(pop_good11, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good12 <- self(pop_good_sel11, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good12 <- setPheno(pop_good12, h2 = 0.8, simParam = SP)
  fancmvar13[i,] <- genicVarA(pop_good12)
  fancmgv11[i,] <- gv(pop_good_sel11)
  
  pop_good_sel12 <- selectInd(pop_good12, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good13 <- self(pop_good_sel12, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good13 <- setPheno(pop_good13, h2 = 0.8, simParam = SP)
  fancmvar14[i,] <- genicVarA(pop_good13)
  fancmgv12[i,] <- gv(pop_good_sel12)
  
  pop_good_sel13 <- selectInd(pop_good13, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good14 <- self(pop_good_sel13, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good14 <- setPheno(pop_good14, h2 = 0.8, simParam = SP)
  fancmvar15[i,] <- genicVarA(pop_good14)
  fancmgv13[i,] <- gv(pop_good_sel13)
  
  pop_good_sel14 <- selectInd(pop_good14, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good15 <- self(pop_good_sel14, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good15 <- setPheno(pop_good15, h2 = 0.8, simParam = SP)
  fancmvar16[i,] <- genicVarA(pop_good15)
  fancmgv14[i,] <- gv(pop_good_sel14)
  
  pop_good_sel15 <- selectInd(pop_good15, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good16 <- self(pop_good_sel15, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good16 <- setPheno(pop_good16, h2 = 0.8, simParam = SP)
  fancmvar17[i,] <- genicVarA(pop_good16)
  fancmgv15[i,] <- gv(pop_good_sel15)
}

#ideal1: 10X global increase in recombination
ideal1var1 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal1var2 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal1var2 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal1var3 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal1var4 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal1var5 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal1var6 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal1var7 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal1var8 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal1var9 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal1var10 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal1var11 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal1var12 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal1var13 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal1var14 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal1var15 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal1var16 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal1var17 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal1gv1 <- matrix(data = NA, ncol = 5, nrow = 100)
ideal1gv2 <- matrix(data = NA, ncol = 5, nrow = 100)
ideal1gv3 <- matrix(data = NA, ncol = 5, nrow = 100)
ideal1gv4 <- matrix(data = NA, ncol = 5, nrow = 100)
ideal1gv5 <- matrix(data = NA, ncol = 5, nrow = 100)
ideal1gv6 <- matrix(data = NA, ncol = 5, nrow = 100)
ideal1gv7 <- matrix(data = NA, ncol = 5, nrow = 100)
ideal1gv8 <- matrix(data = NA, ncol = 5, nrow = 100)
ideal1gv9 <- matrix(data = NA, ncol = 5, nrow = 100)
ideal1gv10 <- matrix(data = NA, ncol = 5, nrow = 100)
ideal1gv11 <- matrix(data = NA, ncol = 5, nrow = 100)
ideal1gv12 <- matrix(data = NA, ncol = 5, nrow = 100)
ideal1gv13 <- matrix(data = NA, ncol = 5, nrow = 100)
ideal1gv14 <- matrix(data = NA, ncol = 5, nrow = 100)
ideal1gv15 <- matrix(data = NA, ncol = 5, nrow = 100)
for(i in 1:100){
  SP$switchGenMap(ideal1_map, centromere = ideal1_centromere)
  ideal1var1[i,] <- genicVarA(burn_in_pop)
  
  pop_good1 <- self(burn_in_pop, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good1 <- setPheno(pop_good1, h2 = 0.8, simParam = SP)
  ideal1var2[i,] <- genicVarA(pop_good1)
  
  pop_good_sel <- selectInd(pop_good1, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good2 <- self(pop_good_sel, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good2 <- setPheno(pop_good2, h2 = 0.8, simParam = SP)
  ideal1var3[i,] <- genicVarA(pop_good2)
  ideal1gv1[i,] <- gv(pop_good_sel)
  
  pop_good_sel2 <- selectInd(pop_good2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good3 <- self(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good3 <- setPheno(pop_good3, h2 = 0.8, simParam = SP)
  ideal1var4[i,] <- genicVarA(pop_good3)
  ideal1gv2[i,] <- gv(pop_good_sel2)
  
  pop_good_sel3 <- selectInd(pop_good3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good4 <- self(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good4 <- setPheno(pop_good4, h2 = 0.8, simParam = SP)
  ideal1var5[i,] <- genicVarA(pop_good4)
  ideal1gv3[i,] <- gv(pop_good_sel3)
  
  pop_good_sel4 <- selectInd(pop_good4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good5 <- self(pop_good_sel4, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good5 <- setPheno(pop_good5, h2 = 0.8, simParam = SP)
  ideal1var6[i,] <- genicVarA(pop_good5)
  ideal1gv4[i,] <- gv(pop_good_sel4)
  
  pop_good_sel5 <- selectInd(pop_good5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good6 <- self(pop_good_sel5, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good6 <- setPheno(pop_good6, h2 = 0.8, simParam = SP)
  ideal1var7[i,] <- genicVarA(pop_good6)
  ideal1gv5[i,] <- gv(pop_good_sel5)
  
  pop_good_sel6 <- selectInd(pop_good6, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good7 <- self(pop_good_sel6, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good7 <- setPheno(pop_good7, h2 = 0.8, simParam = SP)
  ideal1var8[i,] <- genicVarA(pop_good7)
  ideal1gv6[i,] <- gv(pop_good_sel6)
  
  pop_good_sel7 <- selectInd(pop_good7, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good8 <- self(pop_good_sel7, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good8 <- setPheno(pop_good8, h2 = 0.8, simParam = SP)
  ideal1var9[i,] <- genicVarA(pop_good8)
  ideal1gv7[i,] <- gv(pop_good_sel7)
  
  pop_good_sel8 <- selectInd(pop_good8, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good9 <- self(pop_good_sel8, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good9 <- setPheno(pop_good9, h2 = 0.8, simParam = SP)
  ideal1var10[i,] <- genicVarA(pop_good9)
  ideal1gv8[i,] <- gv(pop_good_sel8)
  
  pop_good_sel9 <- selectInd(pop_good9, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good10 <- self(pop_good_sel9, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good10 <- setPheno(pop_good10, h2 = 0.8, simParam = SP)
  ideal1var11[i,] <- genicVarA(pop_good10)
  ideal1gv9[i,] <- gv(pop_good_sel9)
  
  pop_good_sel10 <- selectInd(pop_good10, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good11 <- self(pop_good_sel10, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good11 <- setPheno(pop_good11, h2 = 0.8, simParam = SP)
  ideal1var12[i,] <- genicVarA(pop_good11)
  ideal1gv10[i,] <- gv(pop_good_sel10)
  
  pop_good_sel11 <- selectInd(pop_good11, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good12 <- self(pop_good_sel11, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good12 <- setPheno(pop_good12, h2 = 0.8, simParam = SP)
  ideal1var13[i,] <- genicVarA(pop_good12)
  ideal1gv11[i,] <- gv(pop_good_sel11)
  
  pop_good_sel12 <- selectInd(pop_good12, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good13 <- self(pop_good_sel12, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good13 <- setPheno(pop_good13, h2 = 0.8, simParam = SP)
  ideal1var14[i,] <- genicVarA(pop_good13)
  ideal1gv12[i,] <- gv(pop_good_sel12)
  
  pop_good_sel13 <- selectInd(pop_good13, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good14 <- self(pop_good_sel13, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good14 <- setPheno(pop_good14, h2 = 0.8, simParam = SP)
  ideal1var15[i,] <- genicVarA(pop_good14)
  ideal1gv13[i,] <- gv(pop_good_sel13)
  
  pop_good_sel14 <- selectInd(pop_good14, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good15 <- self(pop_good_sel14, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good15 <- setPheno(pop_good15, h2 = 0.8, simParam = SP)
  ideal1var16[i,] <- genicVarA(pop_good15)
  ideal1gv14[i,] <- gv(pop_good_sel14)
  
  pop_good_sel15 <- selectInd(pop_good15, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good16 <- self(pop_good_sel15, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good16 <- setPheno(pop_good16, h2 = 0.8, simParam = SP)
  ideal1var17[i,] <- genicVarA(pop_good16)
  ideal1gv15[i,] <- gv(pop_good_sel15)
}

#ideal2: ddm1/zmet2 double mutant
ideal2var1 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal2var2 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal2var3 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal2var4 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal2var5 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal2var6 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal2var7 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal2var8 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal2var9 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal2var10 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal2var11 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal2var12 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal2var13 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal2var14 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal2var15 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal2var16 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal2var17 <- matrix(data = NA, ncol = 1, nrow = 100)
ideal2gv1 <- matrix(data = NA, ncol = 5, nrow = 100)
ideal2gv2 <- matrix(data = NA, ncol = 5, nrow = 100)
ideal2gv3 <- matrix(data = NA, ncol = 5, nrow = 100)
ideal2gv4 <- matrix(data = NA, ncol = 5, nrow = 100)
ideal2gv5 <- matrix(data = NA, ncol = 5, nrow = 100)
ideal2gv6 <- matrix(data = NA, ncol = 5, nrow = 100)
ideal2gv7 <- matrix(data = NA, ncol = 5, nrow = 100)
ideal2gv8 <- matrix(data = NA, ncol = 5, nrow = 100)
ideal2gv9 <- matrix(data = NA, ncol = 5, nrow = 100)
ideal2gv10 <- matrix(data = NA, ncol = 5, nrow = 100)
ideal2gv11 <- matrix(data = NA, ncol = 5, nrow = 100)
ideal2gv12 <- matrix(data = NA, ncol = 5, nrow = 100)
ideal2gv13 <- matrix(data = NA, ncol = 5, nrow = 100)
ideal2gv14 <- matrix(data = NA, ncol = 5, nrow = 100)
ideal2gv15 <- matrix(data = NA, ncol = 5, nrow = 100)
for(i in 1:100){
  SP$switchGenMap(ideal2_map, centromere = ideal2_centromere)
  ideal2var1[i,] <- genicVarA(burn_in_pop)
  
  pop_good1 <- self(burn_in_pop, nCrosses = 5, nProgeny=20, simParam = SP)
  pop_good1 <- setPheno(pop_good1, h2 = 0.8, simParam = SP)
  ideal2var2[i,] <- genicVarA(pop_good1)
  
  pop_good_sel <- selectInd(pop_good1, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good2 <- self(pop_good_sel, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good2 <- setPheno(pop_good2, h2 = 0.8, simParam = SP)
  ideal2var3[i,] <- genicVarA(pop_good2)
  ideal2gv1[i,] <- gv(pop_good_sel)
  
  pop_good_sel2 <- selectInd(pop_good2, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good3 <- self(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good3 <- setPheno(pop_good3, h2 = 0.8, simParam = SP)
  ideal2var4[i,] <- genicVarA(pop_good3)
  ideal2gv2[i,] <- gv(pop_good_sel2)
  
  pop_good_sel3 <- selectInd(pop_good3, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good4 <- self(pop_good_sel2, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good4 <- setPheno(pop_good4, h2 = 0.8, simParam = SP)
  ideal2var5[i,] <- genicVarA(pop_good4)
  ideal2gv3[i,] <- gv(pop_good_sel3)
  
  pop_good_sel4 <- selectInd(pop_good4, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good5 <- self(pop_good_sel4, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good5 <- setPheno(pop_good5, h2 = 0.8, simParam = SP)
  ideal2var6[i,] <- genicVarA(pop_good5)
  ideal2gv4[i,] <- gv(pop_good_sel4)
  
  pop_good_sel5 <- selectInd(pop_good5, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good6 <- self(pop_good_sel5, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good6 <- setPheno(pop_good6, h2 = 0.8, simParam = SP)
  ideal2var7[i,] <- genicVarA(pop_good6)
  ideal2gv5[i,] <- gv(pop_good_sel5)
  
  pop_good_sel6 <- selectInd(pop_good6, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good7 <- self(pop_good_sel6, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good7 <- setPheno(pop_good7, h2 = 0.8, simParam = SP)
  ideal2var8[i,] <- genicVarA(pop_good7)
  ideal2gv6[i,] <- gv(pop_good_sel6)
  
  pop_good_sel7 <- selectInd(pop_good7, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good8 <- self(pop_good_sel7, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good8 <- setPheno(pop_good8, h2 = 0.8, simParam = SP)
  ideal2var9[i,] <- genicVarA(pop_good8)
  ideal2gv7[i,] <- gv(pop_good_sel7)
  
  pop_good_sel8 <- selectInd(pop_good8, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good9 <- self(pop_good_sel8, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good9 <- setPheno(pop_good9, h2 = 0.8, simParam = SP)
  ideal2var10[i,] <- genicVarA(pop_good9)
  ideal2gv8[i,] <- gv(pop_good_sel8)
  
  pop_good_sel9 <- selectInd(pop_good9, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good10 <- self(pop_good_sel9, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good10 <- setPheno(pop_good10, h2 = 0.8, simParam = SP)
  ideal2var11[i,] <- genicVarA(pop_good10)
  ideal2gv9[i,] <- gv(pop_good_sel9)
  
  pop_good_sel10 <- selectInd(pop_good10, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good11 <- self(pop_good_sel10, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good11 <- setPheno(pop_good11, h2 = 0.8, simParam = SP)
  ideal2var12[i,] <- genicVarA(pop_good11)
  ideal2gv10[i,] <- gv(pop_good_sel10)
  
  pop_good_sel11 <- selectInd(pop_good11, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good12 <- self(pop_good_sel11, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good12 <- setPheno(pop_good12, h2 = 0.8, simParam = SP)
  ideal2var13[i,] <- genicVarA(pop_good12)
  ideal2gv11[i,] <- gv(pop_good_sel11)
  
  pop_good_sel12 <- selectInd(pop_good12, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good13 <- self(pop_good_sel12, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good13 <- setPheno(pop_good13, h2 = 0.8, simParam = SP)
  ideal2var14[i,] <- genicVarA(pop_good13)
  ideal2gv12[i,] <- gv(pop_good_sel12)
  
  pop_good_sel13 <- selectInd(pop_good13, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good14 <- self(pop_good_sel13, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good14 <- setPheno(pop_good14, h2 = 0.8, simParam = SP)
  ideal2var15[i,] <- genicVarA(pop_good14)
  ideal2gv13[i,] <- gv(pop_good_sel13)
  
  pop_good_sel14 <- selectInd(pop_good14, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good15 <- self(pop_good_sel14, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good15 <- setPheno(pop_good15, h2 = 0.8, simParam = SP)
  ideal2var16[i,] <- genicVarA(pop_good15)
  ideal2gv14[i,] <- gv(pop_good_sel14)
  
  pop_good_sel15 <- selectInd(pop_good15, nInd = 5, use = "gv", trait = 1, selectTop = TRUE, returnPop = TRUE, simParam = SP)
  pop_good16 <- self(pop_good_sel15, nCrosses = 5, nProgeny = 20, simParam = SP)
  pop_good16 <- setPheno(pop_good16, h2 = 0.8, simParam = SP)
  ideal2var17[i,] <- genicVarA(pop_good16)
  ideal2gv15[i,] <- gv(pop_good_sel15)
}

#adding all generations to one data file
#one with just genetic values & one calculating genetic gain

wt <- c(wtgv1[,1:5], wtgv2[,1:5], wtgv3[,1:5], wtgv4[,1:5], wtgv5[,1:5],
        wtgv6[,1:5], wtgv7[,1:5], wtgv8[,1:5], wtgv9[,1:5], wtgv10[,1:5], wtgv11[,1:5],
        wtgv12[,1:5], wtgv13[,1:5], wtgv14[,1:5], wtgv15[,1:5])
wt_gain <- c(wtgv2[,1:5] - wtgv1[,1:5], wtgv3[,1:5] - wtgv1[,1:5], wtgv4[,1:5] - wtgv1[,1:5], wtgv5[,1:5] - wtgv1[,1:5], wtgv6[,1:5] - wtgv1[,1:5],
             wtgv7[,1:5] - wtgv1[,1:5], wtgv8[,1:5] - wtgv1[,1:5], wtgv9[,1:5] - wtgv1[,1:5], wtgv10[,1:5] - wtgv1[,1:5], wtgv11[,1:5] - wtgv1[,1:5],
             wtgv12[,1:5] - wtgv1[,1:5], wtgv13[,1:5] - wtgv1[,1:5], wtgv14[,1:5] - wtgv1[,1:5], wtgv15[,1:5] - wtgv1[,1:5])

ddm1 <- c(ddm1gv1[,1:5], ddm1gv2[,1:5], ddm1gv3[,1:5], ddm1gv4[,1:5], ddm1gv5[,1:5],
          ddm1gv6[,1:5], ddm1gv7[,1:5], ddm1gv8[,1:5], ddm1gv9[,1:5], ddm1gv10[,1:5], ddm1gv11[,1:5],
          ddm1gv12[,1:5], ddm1gv13[,1:5], ddm1gv14[,1:5], ddm1gv15[,1:5])
ddm1_gain <- c(ddm1gv2[,1:5]-ddm1gv1[,1:5], ddm1gv3[,1:5]-ddm1gv1[,1:5], ddm1gv4[,1:5]-ddm1gv1[,1:5], ddm1gv5[,1:5]-ddm1gv1[,1:5], ddm1gv6[,1:5]-ddm1gv1[,1:5],
               ddm1gv7[,1:5]-ddm1gv1[,1:5], ddm1gv8[,1:5]-ddm1gv1[,1:5], ddm1gv9[,1:5]-ddm1gv1[,1:5], ddm1gv10[,1:5]-ddm1gv1[,1:5],
               ddm1gv11[,1:5]-ddm1gv1[,1:5], ddm1gv12[,1:5]-ddm1gv1[,1:5], ddm1gv13[,1:5]-ddm1gv1[,1:5], ddm1gv14[,1:5]-ddm1gv1[,1:5],ddm1gv15[,1:5]-ddm1gv1[,1:5])

zmet2 <- c(zmet2gv1[,1:5], zmet2gv2[,1:5], zmet2gv3[,1:5], zmet2gv4[,1:5], zmet2gv5[,1:5],
           zmet2gv6[,1:5], zmet2gv7[,1:5], zmet2gv8[,1:5], zmet2gv9[,1:5], zmet2gv10[,1:5], 
           zmet2gv11[,1:5], zmet2gv12[,1:5], zmet2gv13[,1:5], zmet2gv14[,1:5], zmet2gv15[,1:5])
zmet2_gain <- c(zmet2gv2[,1:5]-zmet2gv1[,1:5], zmet2gv3[,1:5]-zmet2gv1[,1:5], zmet2gv4[,1:5]-zmet2gv1[,1:5], zmet2gv5[,1:5]-zmet2gv1[,1:5], zmet2gv6[,1:5]-zmet2gv1[,1:5],
                zmet2gv7[,1:5]-zmet2gv1[,1:5], zmet2gv8[,1:5]-zmet2gv1[,1:5], zmet2gv9[,1:5]-zmet2gv1[,1:5], zmet2gv10[,1:5]-zmet2gv1[,1:5], zmet2gv11[,1:5]-zmet2gv1[,1:5],
                zmet2gv12[,1:5]-zmet2gv1[,1:5], zmet2gv13[,1:5]-zmet2gv1[,1:5], zmet2gv14[,1:5]-zmet2gv1[,1:5], zmet2gv15[,1:5]-zmet2gv1[,1:5])

recq4 <- c(recq4gv1[,1:5], recq4gv2[,1:5], recq4gv3[,1:5], recq4gv4[,1:5], recq4gv5[,1:5],
           recq4gv6[,1:5], recq4gv7[,1:5], recq4gv8[,1:5], recq4gv9[,1:5], recq4gv10[,1:5], 
           recq4gv11[,1:5], recq4gv12[,1:5], recq4gv13[,1:5], recq4gv14[,1:5], recq4gv15[,1:5])
recq4_gain <- c(recq4gv2[,1:5]-recq4gv1[,1:5], recq4gv3[,1:5]-recq4gv1[,1:5], recq4gv4[,1:5]-recq4gv1[,1:5], recq4gv5[,1:5]-recq4gv1[,1:5],
                recq4gv6[,1:5]-recq4gv1[,1:5], recq4gv7[,1:5]-recq4gv1[,1:5], recq4gv8[,1:5]-recq4gv1[,1:5], recq4gv9[,1:5]-recq4gv1[,1:5],
                recq4gv10[,1:5]-recq4gv1[,1:5], recq4gv11[,1:5]-recq4gv1[,1:5], recq4gv12[,1:5]-recq4gv1[,1:5], recq4gv13[,1:5]-recq4gv1[,1:5],
                recq4gv14[,1:5]-recq4gv1[,1:5], recq4gv15[,1:5]-recq4gv1[,1:5])

fancm <- c(fancmgv1[,1:5], fancmgv2[,1:5], fancmgv3[,1:5], fancmgv4[,1:5], fancmgv5[,1:5],
           fancmgv6[,1:5], fancmgv7[,1:5], fancmgv8[,1:5], fancmgv9[,1:5], fancmgv10[,1:5], 
           fancmgv11[,1:5], fancmgv12[,1:5], fancmgv13[,1:5], fancmgv14[,1:5], fancmgv15[,1:5])

fancm_gain <- c(fancmgv2[,1:5]-fancmgv1[,1:5], fancmgv3[,1:5]-fancmgv1[,1:5], fancmgv4[,1:5]-fancmgv1[,1:5], fancmgv5[,1:5]-fancmgv1[,1:5],
                fancmgv6[,1:5]-fancmgv1[,1:5], fancmgv7[,1:5]-fancmgv1[,1:5], fancmgv8[,1:5]-fancmgv1[,1:5], fancmgv9[,1:5]-fancmgv1[,1:5],
                fancmgv10[,1:5]-fancmgv1[,1:5], fancmgv11[,1:5]-fancmgv1[,1:5], fancmgv12[,1:5]-fancmgv1[,1:5], fancmgv13[,1:5]-fancmgv1[,1:5],
                fancmgv14[,1:5]-fancmgv1[,1:5], fancmgv15[,1:5]-fancmgv1[,1:5])

ideal1 <- c(ideal1gv1[,1:5], ideal1gv2[,1:5], ideal1gv3[,1:5], ideal1gv4[,1:5], ideal1gv5[,1:5],
            ideal1gv6[,1:5], ideal1gv7[,1:5], ideal1gv8[,1:5], ideal1gv9[,1:5], ideal1gv10[,1:5], 
            ideal1gv11[,1:5], ideal1gv12[,1:5], ideal1gv13[,1:5], ideal1gv14[,1:5], ideal1gv15[,1:5])

ideal1_gain <- c(ideal1gv2[,1:5]-ideal1gv1[,1:5], ideal1gv3[,1:5]-ideal1gv1[,1:5], ideal1gv4[,1:5]-ideal1gv1[,1:5], ideal1gv5[,1:5]-ideal1gv1[,1:5],
                 ideal1gv6[,1:5]-ideal1gv1[,1:5], ideal1gv7[,1:5]-ideal1gv1[,1:5], ideal1gv8[,1:5]-ideal1gv1[,1:5], ideal1gv9[,1:5]-ideal1gv1[,1:5],
                 ideal1gv10[,1:5]-ideal1gv1[,1:5], ideal1gv11[,1:5]-ideal1gv1[,1:5], ideal1gv12[,1:5]-ideal1gv1[,1:5], ideal1gv13[,1:5]-ideal1gv1[,1:5],
                 ideal1gv14[,1:5]-ideal1gv1[,1:5], ideal1gv15[,1:5]-ideal1gv1[,1:5])

ideal2 <- c(ideal2gv1[,1:5], ideal2gv2[,1:5], ideal2gv3[,1:5], ideal2gv4[,1:5], ideal2gv5[,1:5],
            ideal2gv6[,1:5], ideal2gv7[,1:5], ideal2gv8[,1:5], ideal2gv9[,1:5], ideal2gv10[,1:5], 
            ideal2gv11[,1:5], ideal2gv12[,1:5], ideal2gv13[,1:5], ideal2gv14[,1:5], ideal2gv15[,1:5])

ideal2_gain <- c(ideal2gv2[,1:5]-ideal2gv1[,1:5], ideal2gv3[,1:5]-ideal2gv1[,1:5], ideal2gv4[,1:5]-ideal2gv1[,1:5], ideal2gv5[,1:5]-ideal2gv1[,1:5],
                 ideal2gv6[,1:5]-ideal2gv1[,1:5], ideal2gv7[,1:5]-ideal2gv1[,1:5], ideal2gv8[,1:5]-ideal2gv1[,1:5], ideal2gv9[,1:5]-ideal2gv1[,1:5],
                 ideal2gv10[,1:5]-ideal2gv1[,1:5], ideal2gv11[,1:5]-ideal2gv1[,1:5], ideal2gv12[,1:5]-ideal2gv1[,1:5], ideal2gv13[,1:5]-ideal2gv1[,1:5],
                 ideal2gv14[,1:5]-ideal2gv1[,1:5], ideal2gv15[,1:5]-ideal2gv1[,1:5])

#combining all genetic value data for wt & all mutants
all <- cbind(wt, zmet2, ddm1, recq4, fancm, ideal1, ideal2)
all <- as.data.frame(all)
all2 <- as.data.frame(matrix(data = NA, nrow = 52500))
all2$gv <- c(all$wt, all$ddm1, all$zmet2, all$recq4, all$fancm, all$ideal1, all$ideal2)
all2$generation <- rep(1:15, size = 7, each = 500)
all2$gen_map <- rep(c("wt","ddm1", "zmet", "recq4", "fancm", "ideal1", "ideal2"), size = 7, each = 7500)
lastgen <- all2[which(all2$generation == 15),]

#combining all genetic gain data for wt & all mutants
gain_inter <- cbind(wt_gain, zmet2_gain, ddm1_gain, recq4_gain, fancm_gain, ideal1_gain, ideal2_gain)
gain_inter <- as.data.frame(gain_inter)
gain_inter2 <- as.data.frame(matrix(data = NA, nrow = 49000))
gain_inter2$gv <- c(gain_inter$wt_gain, gain_inter$ddm1_gain, gain_inter$zmet2_gain, gain_inter$recq4_gain, gain_inter$fancm_gain, 
                    gain_inter$ideal1_gain, gain_inter$ideal2_gain)
gain_inter2$generation <- rep(1:14, size = 7, each = 500)
gain_inter2$gen_map <- rep(c("wt","ddm1", "zmet2", "recq4", "fancm", "ideal1", "ideal2"), size = 7, each = 7000)
lastgen <- gain_inter2[which(gain_inter2$generation == c(1,3,5,7,9,11,13,14)),]
firstgen <- gain_inter2[which(gain_inter2$generation == 1),]

#colorblind friendly palette with black
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
group.colors = c("ddm1" = "#E69F00", "recq4" = "#56B4E9", "wt" = "#009E73", "zmet2" = "#F0E442",
                 "fancm" = "#0072B2", "ideal1" = "#D55E00", "ideal2" = "#CC79A7")
ggplot(lastgen, aes(x=as.factor(generation), y=gv, fill=gen_map)) + 
  geom_boxplot() + theme_bw() + xlab("Generations since generation 1") + ylab("Genetic gain since generation 1") +
  scale_fill_manual(values=group.colors, name = "Genetic Map Used", labels = c("ddm1", "recq4", "WT", "zmet2",
                                                                               "fancm", "10X", "ddm1/zmet2")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12), legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=12), plot.title = element_text(size=14), legend.position = c(0.15,0.8), legend.key.size = unit(0.5, "lines")) +
  guides(color = guide_legend(override.aes = list(size = 2)))


#example of shaded SD
ggplot(gain_inter2, aes(x=generation,y=gv, colour=gen_map)) + geom_line() +
  stat_smooth(method="loess", span=0.01, se=TRUE, aes(fill=gen_map), alpha=0.1) +
  theme_bw()

#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
wtvarall <- c(wtvar1, wtvar2, wtvar3, wtvar4, wtvar5,
              wtvar6,  wtvar7,  wtvar8,  wtvar9,  wtvar10, wtvar11,
              wtvar12, wtvar13,  wtvar14,  wtvar15,  wtvar16,  wtvar17)

ddm1varall <- c(ddm1var1, ddm1var2, ddm1var3, ddm1var4, ddm1var5,
                ddm1var6,  ddm1var7,  ddm1var8,  ddm1var9,  ddm1var10, 
                ddm1var11, ddm1var12, ddm1var13, ddm1var14, ddm1var15, ddm1var16, ddm1var17)

zmet2varall <- c(zmet2var1, zmet2var2, zmet2var3, zmet2var4, zmet2var5,
                 zmet2var6,  zmet2var7,  zmet2var8,  zmet2var9,  zmet2var10,
                 zmet2var11, zmet2var12, zmet2var13, zmet2var14, zmet2var15, zmet2var16, zmet2var17)

recq4varall <- c(recq4var1, recq4var2, recq4var3, recq4var4, recq4var5,
                 recq4var6,  recq4var7,  recq4var8,  recq4var9,  recq4var10,
                 recq4var11, recq4var12, recq4var13, recq4var14, recq4var15, recq4var16, recq4var17)

ideal1varall <- c(ideal1var1, ideal1var2, ideal1var3, ideal1var4, ideal1var5,
                  ideal1var6,  ideal1var7,  ideal1var8,  ideal1var9,  ideal1var10,
                  ideal1var11, ideal1var12, ideal1var13, ideal1var14, ideal1var15, ideal1var16, ideal1var17)

ideal2varall <- c(ideal2var1, ideal2var2, ideal2var3, ideal2var4, ideal2var5,
                  ideal2var6,  ideal2var7,  ideal2var8,  ideal2var9,  ideal2var10,
                  ideal2var11, ideal2var12, ideal2var13, ideal2var14, ideal2var15, ideal2var16, ideal2var17)

ddm1_zmet2varall <- c(ideal2var1, ideal2var2, ideal2var3, ideal2var4, ideal2var5,
                  ideal2var6,  ideal2var7,  ideal2var8,  ideal2var9,  ideal2var10,
                  ideal2var11, ideal2var12, ideal2var13, ideal2var14, ideal2var15, ideal2var16, ideal2var17)

fancmvarall <- c(fancmvar1, fancmvar2, fancmvar3, fancmvar4, fancmvar5,
                 fancmvar6,  fancmvar7,  fancmvar8,  fancmvar9,  fancmvar10,
                 fancmvar11, fancmvar12, fancmvar13, fancmvar14, fancmvar15, fancmvar16, fancmvar17)

allvar <- cbind(wtvarall, ddm1varall, zmet2varall, recq4varall, fancmvarall, ideal1varall, ideal2varall)
allvar <- as.data.frame(allvar)
allvar2 <- as.data.frame(matrix(data = NA, nrow = 11900))
allvar2$var <- c(allvar$wtvarall, allvar$ddm1varall, allvar$zmet2varall, allvar$recq4varall, allvar$fancmvarall, allvar$ideal1varall, allvar$ideal2varall)
allvar2$gen <- rep(1:17, size = 7, each = 100)
allvar2$gen_map <- rep(c("wt","ddm1", "zmet2", "recq4", "fancm", "ideal1", "ideal2"), size = 7, each = 1700)
lastgen <- allvar2[which(allvar2$gen == c(1,3,5,7,9,11,13,15,16)),]

ggplot(lastgen, aes(x=as.factor(gen), y=var, fill=gen_map)) + 
  geom_boxplot() + theme_bw() + xlab("Generation") + ylab("Additive Genetic Variance") + 
  scale_fill_manual(values=group.colors, name = "Genetic Map Used", labels = c("ddm1", "recq4", "WT", "zmet2",
                                                                                 "fancm", "10X", "ddm1/zmet2")) + ggtitle("Additive Genetic Variance after 15 Generations") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12), legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=12), plot.title = element_text(size=14), legend.position = c(0.8,0.75), legend.key.size = unit(0.5, "lines"))
