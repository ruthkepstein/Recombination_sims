library(AlphaSimR)
library(Rcpp)
library(ggplot2)
library(dplyr)
#Set seed to ensure sampling is the same each time
set.seed(420)

##Read in wildtype recombination rate data and create separate dataframe for each chromosome
japonica_CO <- read.table("japonica_rec_rate.bed", header = FALSE)
colnames(japonica_CO) <- c("Chr", "CO Start", "CO End", "rate")
japonica_CO <- japonica_CO[order(japonica_CO$Chr,japonica_CO$`CO Start`),]

japonica_chr1_CO <- japonica_CO[ which(japonica_CO$Chr == "chr01"),]
japonica_chr1_CO$midpoint <- (japonica_chr1_CO$`CO Start`+ japonica_chr1_CO$`CO End`)/2
japonica_chr1_CO <- japonica_chr1_CO[order(japonica_chr1_CO$`CO Start`),]

japonica_chr2_CO <- japonica_CO[ which(japonica_CO$Chr == "chr02"),]
japonica_chr2_CO$midpoint <- (japonica_chr2_CO$`CO Start`+ japonica_chr2_CO$`CO End`)/2
japonica_chr2_CO <- japonica_chr2_CO[order(japonica_chr2_CO$`CO Start`),]

japonica_chr3_CO <- japonica_CO[ which(japonica_CO$Chr == "chr03"),]
japonica_chr3_CO$midpoint <- (japonica_chr3_CO$`CO Start`+ japonica_chr3_CO$`CO End`)/2
japonica_chr3_CO <- japonica_chr3_CO[order(japonica_chr3_CO$`CO Start`),]

japonica_chr4_CO <- japonica_CO[ which(japonica_CO$Chr == "chr04"),]
japonica_chr4_CO$midpoint <- (japonica_chr4_CO$`CO Start`+ japonica_chr4_CO$`CO End`)/2
japonica_chr4_CO <- japonica_chr4_CO[order(japonica_chr4_CO$`CO Start`),]

japonica_chr5_CO <- japonica_CO[ which(japonica_CO$Chr == "chr05"),]
japonica_chr5_CO$midpoint <- (japonica_chr5_CO$`CO Start`+ japonica_chr5_CO$`CO End`)/2
japonica_chr5_CO <- japonica_chr5_CO[order(japonica_chr5_CO$`CO Start`),]

japonica_chr6_CO <- japonica_CO[ which(japonica_CO$Chr == "chr06"),]
japonica_chr6_CO$midpoint <- (japonica_chr6_CO$`CO Start`+ japonica_chr6_CO$`CO End`)/2
japonica_chr6_CO <- japonica_chr6_CO[order(japonica_chr6_CO$`CO Start`),]

japonica_chr7_CO <- japonica_CO[ which(japonica_CO$Chr == "chr07"),]
japonica_chr7_CO$midpoint <- (japonica_chr7_CO$`CO Start`+ japonica_chr7_CO$`CO End`)/2
japonica_chr7_CO <- japonica_chr7_CO[order(japonica_chr7_CO$`CO Start`),]

japonica_chr8_CO <- japonica_CO[ which(japonica_CO$Chr == "chr08"),]
japonica_chr8_CO$midpoint <- (japonica_chr8_CO$`CO Start`+ japonica_chr8_CO$`CO End`)/2
japonica_chr8_CO <- japonica_chr8_CO[order(japonica_chr8_CO$`CO Start`),]

japonica_chr9_CO <- japonica_CO[ which(japonica_CO$Chr == "chr09"),]
japonica_chr9_CO$midpoint <- (japonica_chr9_CO$`CO Start`+ japonica_chr9_CO$`CO End`)/2
japonica_chr9_CO <- japonica_chr9_CO[order(japonica_chr9_CO$`CO Start`),]

japonica_chr10_CO <- japonica_CO[ which(japonica_CO$Chr == "chr10"),]
japonica_chr10_CO$midpoint <- (japonica_chr10_CO$`CO Start`+ japonica_chr10_CO$`CO End`)/2
japonica_chr10_CO <- japonica_chr10_CO[order(japonica_chr10_CO$`CO Start`),]

japonica_chr11_CO <- japonica_CO[ which(japonica_CO$Chr == "chr11"),]
japonica_chr11_CO$midpoint <- (japonica_chr11_CO$`CO Start`+ japonica_chr11_CO$`CO End`)/2
japonica_chr11_CO <- japonica_chr11_CO[order(japonica_chr11_CO$`CO Start`),]

japonica_chr12_CO <- japonica_CO[ which(japonica_CO$Chr == "chr12"),]
japonica_chr12_CO$midpoint <- (japonica_chr12_CO$`CO Start`+ japonica_chr12_CO$`CO End`)/2
japonica_chr12_CO <- japonica_chr12_CO[order(japonica_chr12_CO$`CO Start`),]

#Making each recombination rate intervals start at zero
japonica_chr1_CO$`CO End` <- japonica_chr1_CO$`CO End` - min(japonica_chr1_CO$`CO Start`)
japonica_chr1_CO$`CO Start` <- japonica_chr1_CO$`CO Start` - min(japonica_chr1_CO$`CO Start`)

japonica_chr2_CO$`CO End` <- japonica_chr2_CO$`CO End` - min(japonica_chr2_CO$`CO Start`)
japonica_chr2_CO$`CO Start` <- japonica_chr2_CO$`CO Start` - min(japonica_chr2_CO$`CO Start`)

japonica_chr3_CO$`CO End` <- japonica_chr3_CO$`CO End` - min(japonica_chr3_CO$`CO Start`)
japonica_chr3_CO$`CO Start` <- japonica_chr3_CO$`CO Start` - min(japonica_chr3_CO$`CO Start`)

japonica_chr4_CO$`CO End` <- japonica_chr4_CO$`CO End` - min(japonica_chr4_CO$`CO Start`)
japonica_chr4_CO$`CO Start` <- japonica_chr4_CO$`CO Start` - min(japonica_chr4_CO$`CO Start`)

japonica_chr5_CO$`CO End` <- japonica_chr5_CO$`CO End` - min(japonica_chr5_CO$`CO Start`)
japonica_chr5_CO$`CO Start` <- japonica_chr5_CO$`CO Start` - min(japonica_chr5_CO$`CO Start`)

japonica_chr6_CO$`CO End` <- japonica_chr6_CO$`CO End` - min(japonica_chr6_CO$`CO Start`)
japonica_chr6_CO$`CO Start` <- japonica_chr6_CO$`CO Start` - min(japonica_chr6_CO$`CO Start`)

japonica_chr7_CO$`CO End` <- japonica_chr7_CO$`CO End` - min(japonica_chr7_CO$`CO Start`)
japonica_chr7_CO$`CO Start` <- japonica_chr7_CO$`CO Start` - min(japonica_chr7_CO$`CO Start`)

japonica_chr8_CO$`CO End` <- japonica_chr8_CO$`CO End` - min(japonica_chr8_CO$`CO Start`)
japonica_chr8_CO$`CO Start` <- japonica_chr8_CO$`CO Start` - min(japonica_chr8_CO$`CO Start`)

japonica_chr9_CO$`CO End` <- japonica_chr9_CO$`CO End` - min(japonica_chr9_CO$`CO Start`)
japonica_chr9_CO$`CO Start` <- japonica_chr9_CO$`CO Start` - min(japonica_chr9_CO$`CO Start`)

japonica_chr10_CO$`CO End` <- japonica_chr10_CO$`CO End` - min(japonica_chr10_CO$`CO Start`)
japonica_chr10_CO$`CO Start` <- japonica_chr10_CO$`CO Start` - min(japonica_chr10_CO$`CO Start`)

japonica_chr11_CO$`CO End` <- japonica_chr11_CO$`CO End` - min(japonica_chr11_CO$`CO Start`)
japonica_chr11_CO$`CO Start` <- japonica_chr11_CO$`CO Start` - min(japonica_chr11_CO$`CO Start`)

japonica_chr12_CO$`CO End` <- japonica_chr12_CO$`CO End` - min(japonica_chr12_CO$`CO Start`)
japonica_chr12_CO$`CO Start` <- japonica_chr12_CO$`CO Start` - min(japonica_chr12_CO$`CO Start`)

library(dlookr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(OneR)

##Bin the recombination rate data into 40 equally sized bins on each chromosome
#Function to label the start SNP positions of each bin
fill_start<- function(chr_CO){
  l<-0
  for(k in 1:nrow(chr_CO)){
    if(isFALSE(is.na(chr_CO$rates[k]))){
      temp<-chr_CO$`CO Start`[k]
      chr_CO$`CO Start`[k] <-l
      l<-temp
    }
  }
  print(chr_CO)
}
#Bin size (variable: bins) calculated by dividing the total size of dataframe by 40
#Rollapply function "rolls" by interval of bin size and calculates a moving average recombination rate for each bin
#Use drop_NA function to condense table so it only contains the avg recombination rate
library(zoo)
library(tidyr)
japonica_chr1_CO_2 <- japonica_chr1_CO
bins<-as.integer(nrow(japonica_chr1_CO)/40)
japonica_chr1_CO_2$rates<- rollapply(japonica_chr1_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr1_CO_2<-fill_start(japonica_chr1_CO_2)
japonica_chr1_CO_2<- japonica_chr1_CO_2 %>% drop_na(rates)

japonica_chr2_CO_2 <- japonica_chr2_CO
bins<-as.integer(nrow(japonica_chr2_CO)/40)
japonica_chr2_CO_2$rates<- rollapply(japonica_chr2_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr2_CO_2<-fill_start(japonica_chr2_CO_2)
japonica_chr2_CO_2<- japonica_chr2_CO_2 %>% drop_na(rates)

japonica_chr3_CO_2 <- japonica_chr3_CO
bins<-as.integer(nrow(japonica_chr3_CO)/40)
japonica_chr3_CO_2$rates<- rollapply(japonica_chr3_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr3_CO_2<-fill_start(japonica_chr3_CO_2)
japonica_chr3_CO_2<- japonica_chr3_CO_2 %>% drop_na(rates)

japonica_chr4_CO_2 <- japonica_chr4_CO
bins<-as.integer(nrow(japonica_chr4_CO)/40)
japonica_chr4_CO_2$rates<- rollapply(japonica_chr4_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr4_CO_2<-fill_start(japonica_chr4_CO_2)
japonica_chr4_CO_2<- japonica_chr4_CO_2 %>% drop_na(rates)

japonica_chr5_CO_2 <- japonica_chr5_CO
bins<-as.integer(nrow(japonica_chr5_CO)/40)
japonica_chr5_CO_2$rates<- rollapply(japonica_chr5_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr5_CO_2<-fill_start(japonica_chr5_CO_2)
japonica_chr5_CO_2<- japonica_chr5_CO_2 %>% drop_na(rates)

japonica_chr6_CO_2 <- japonica_chr6_CO
bins<-as.integer(nrow(japonica_chr6_CO)/40)
japonica_chr6_CO_2$rates<- rollapply(japonica_chr6_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr6_CO_2<-fill_start(japonica_chr6_CO_2)
japonica_chr6_CO_2<- japonica_chr6_CO_2 %>% drop_na(rates)

japonica_chr7_CO_2 <- japonica_chr7_CO
bins<-as.integer(nrow(japonica_chr7_CO)/40)
japonica_chr7_CO_2$rates<- rollapply(japonica_chr7_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr7_CO_2<-fill_start(japonica_chr7_CO_2)
japonica_chr7_CO_2<- japonica_chr7_CO_2 %>% drop_na(rates)

japonica_chr8_CO_2 <- japonica_chr8_CO
bins<-as.integer(nrow(japonica_chr8_CO)/40)
japonica_chr8_CO_2$rates<- rollapply(japonica_chr8_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr8_CO_2<-fill_start(japonica_chr8_CO_2)
japonica_chr8_CO_2<- japonica_chr8_CO_2 %>% drop_na(rates)

japonica_chr9_CO_2 <- japonica_chr9_CO
bins<-as.integer(nrow(japonica_chr9_CO)/40)
japonica_chr9_CO_2$rates<- rollapply(japonica_chr9_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr9_CO_2<-fill_start(japonica_chr9_CO_2)
japonica_chr9_CO_2<- japonica_chr9_CO_2 %>% drop_na(rates)

japonica_chr10_CO_2 <- japonica_chr10_CO
bins<-as.integer(nrow(japonica_chr10_CO)/40)
japonica_chr10_CO_2$rates<- rollapply(japonica_chr10_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr10_CO_2<-fill_start(japonica_chr10_CO_2)
japonica_chr10_CO_2<- japonica_chr10_CO_2 %>% drop_na(rates)

japonica_chr11_CO_2 <- japonica_chr11_CO
bins<-as.integer(nrow(japonica_chr11_CO)/40)
japonica_chr11_CO_2$rates<- rollapply(japonica_chr11_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr11_CO_2<-fill_start(japonica_chr11_CO_2)
japonica_chr11_CO_2<- japonica_chr11_CO_2 %>% drop_na(rates)

japonica_chr12_CO_2 <- japonica_chr12_CO
bins<-as.integer(nrow(japonica_chr12_CO)/40)
japonica_chr12_CO_2$rates<- rollapply(japonica_chr12_CO$rate, width=bins, FUN=mean, by = bins, by.column = TRUE, fill = NA)
japonica_chr12_CO_2<-fill_start(japonica_chr12_CO_2)
japonica_chr12_CO_2<- japonica_chr12_CO_2 %>% drop_na(rates)

##Read in genes on each chromosomes, calculate gene density by binning
#Genome annnotation file from Nipponbare (https://rapdb.dna.affrc.go.jp/download/irgsp1.html)
ref_genes <- read.table("locus.csv", header = TRUE, sep =",")
ref_genes1 <- ref_genes[which(ref_genes$Ref == 'chr01'),]
ref_genes2 <- ref_genes[which(ref_genes$Ref == 'chr02'),]
ref_genes3 <- ref_genes[which(ref_genes$Ref == 'chr03'),]
ref_genes4 <- ref_genes[which(ref_genes$Ref == 'chr04'),]
ref_genes5 <- ref_genes[which(ref_genes$Ref == 'chr05'),]
ref_genes6 <- ref_genes[which(ref_genes$Ref == 'chr06'),]
ref_genes7 <- ref_genes[which(ref_genes$Ref == 'chr07'),]
ref_genes8 <- ref_genes[which(ref_genes$Ref == 'chr08'),]
ref_genes9 <- ref_genes[which(ref_genes$Ref == 'chr09'),]
ref_genes10 <- ref_genes[which(ref_genes$Ref == 'chr10'),]
ref_genes11 <- ref_genes[which(ref_genes$Ref == 'chr11'),]
ref_genes12 <- ref_genes[which(ref_genes$Ref == 'chr12'),]

#Bin genes on each chromosome into equal (~1 kb) intervals, create dataframe with bin info (position, freq, rate, gene density)  
#Variables foo.X1 and foo.X2 represent interval start and end positions, respectively
genes_bin1 <- as.data.frame(summary(binning(ref_genes1$X1, nbins = round(max(ref_genes1$X307041717)/100000), type = "kmeans")))
genes_bin1 <- within(genes_bin1, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin1$levels), ',', fixed=TRUE))))
genes_bin1 <- do.call(data.frame, genes_bin1)
genes_bin1 <- genes_bin1 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin1 <- genes_bin1 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin1[1,4]<-2982
genes_bin1$length <- (genes_bin1$foo.X2)-(genes_bin1$foo.X1)
genes_bin1[round(max(ref_genes1$X307041717)/100000),5] <- max(ref_genes1$X307041717)
genes_bin1$density <- (genes_bin1$freq/(genes_bin1$length/1000000))/1000
plot(genes_bin1$foo.X1/1000000, genes_bin1$density, type = "l")
saveRDS(genes_bin1$foo.X1/1000000,file="chr1_locipos")
saveRDS(genes_bin1$density,file="chr1_geneProb")

genes_bin2 <- as.data.frame(summary(binning(ref_genes2$X1, nbins = round(max(ref_genes2$X307041717)/100000), type = "kmeans")))
genes_bin2 <- within(genes_bin2, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin2$levels), ',', fixed=TRUE))))
genes_bin2 <- do.call(data.frame, genes_bin2)
genes_bin2 <- genes_bin2 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin2 <- genes_bin2 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin2[1,4]<-391
genes_bin2$length <- (genes_bin2$foo.X2-genes_bin2$foo.X1)
genes_bin2[round(max(ref_genes2$X307041717)/100000),5] <- max(ref_genes2$X307041717)
genes_bin2$density <- (genes_bin2$freq/(genes_bin2$length/1000000))/1000
plot(genes_bin2$foo.X1/1000000, genes_bin2$density,type="l")
saveRDS(genes_bin2$foo.X1/1000000,file="chr2_locipos")
saveRDS(genes_bin2$density,file="chr2_geneProb")

genes_bin3 <- as.data.frame(summary(binning(ref_genes3$X1, nbins = round(max(ref_genes3$X307041717)/100000), type = "kmeans")))
genes_bin3 <- within(genes_bin3, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin3$levels), ',', fixed=TRUE))))
genes_bin3 <- do.call(data.frame, genes_bin3)
genes_bin3 <- genes_bin3 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin3 <- genes_bin3 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin3[1,4]<-4582
genes_bin3$length <- (genes_bin3$foo.X2-genes_bin3$foo.X1)
genes_bin3[round(max(ref_genes4$X307041717)/100000),5] <- max(ref_genes3$X307041717)
genes_bin3$density <- (genes_bin3$freq/(genes_bin3$length/1000000))/1000
plot(genes_bin3$foo.X1/1000000, genes_bin3$density,type="l")
saveRDS(genes_bin3$foo.X1/1000000,file="chr3_locipos")
saveRDS(genes_bin3$density,file="chr3_geneProb")

genes_bin4 <- as.data.frame(summary(binning(ref_genes4$X1, nbins = round(max(ref_genes4$X307041717)/100000), type = "kmeans")))
genes_bin4 <- within(genes_bin4, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin4$levels), ',', fixed=TRUE))))
genes_bin4 <- do.call(data.frame, genes_bin4)
genes_bin4 <- genes_bin4 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin4 <- genes_bin4 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin4[1,4]<-58788
genes_bin4$length <- (genes_bin4$foo.X2-genes_bin4$foo.X1)
genes_bin4[round(max(ref_genes4$X307041717)/100000),5] <- max(ref_genes4$X307041717)
genes_bin4$density <- (genes_bin4$freq/(genes_bin4$length/1000000))/1000
plot(genes_bin4$foo.X1/1000000, genes_bin4$density,type="l")
saveRDS(genes_bin4$foo.X1/1000000,file="chr4_locipos")
saveRDS(genes_bin4$density,file="chr4_geneProb")

genes_bin5 <- as.data.frame(summary(binning(ref_genes5$X1, nbins = round(max(ref_genes5$X307041717)/100000), type = "kmeans")))
genes_bin5 <- within(genes_bin5, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin5$levels), ',', fixed=TRUE))))
genes_bin5 <- do.call(data.frame, genes_bin5)
genes_bin5 <- genes_bin5 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin5 <- genes_bin5 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin5[1,4]<-10946
genes_bin5$length <- (genes_bin5$foo.X2-genes_bin5$foo.X1)
genes_bin5[round(max(ref_genes5$X307041717)/100000),5] <- max(ref_genes5$X307041717)
genes_bin5$density <- (genes_bin5$freq/(genes_bin5$length/1000000))/1000
plot(genes_bin5$foo.X1/1000000, genes_bin5$density,type="l")
saveRDS(genes_bin5$foo.X1/1000000,file="chr5_locipos")
saveRDS(genes_bin5$density,file="chr5_geneProb")

genes_bin6 <- as.data.frame(summary(binning(ref_genes6$X1, nbins = round(max(ref_genes6$X307041717)/100000), type = "kmeans")))
genes_bin6 <- within(genes_bin6, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin6$levels), ',', fixed=TRUE))))
genes_bin6 <- do.call(data.frame, genes_bin6)
genes_bin6 <- genes_bin6 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin6 <- genes_bin6 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin6[1,4]<-38339
genes_bin6$length <- (genes_bin6$foo.X2-genes_bin6$foo.X1)
genes_bin6[round(max(ref_genes6$X307041717)/100000),5] <- max(ref_genes6$X307041717)
genes_bin6$density <- (genes_bin6$freq/(genes_bin6$length/1000000))/1000
plot(genes_bin6$foo.X1/1000000, genes_bin6$density,type="l")
saveRDS(genes_bin6$foo.X1/1000000,file="chr6_locipos")
saveRDS(genes_bin6$density,file="chr6_geneProb")

genes_bin7 <- as.data.frame(summary(binning(ref_genes7$X1, nbins = round(max(ref_genes7$X307041717)/100000), type = "kmeans")))
genes_bin7 <- within(genes_bin7, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin7$levels), ',', fixed=TRUE))))
genes_bin7 <- do.call(data.frame, genes_bin7)
genes_bin7 <- genes_bin7 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin7 <- genes_bin7 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin7[1,4]<-11647
genes_bin7$length <- (genes_bin7$foo.X2-genes_bin7$foo.X1)
genes_bin7[round(max(ref_genes7$X307041717)/100000),5] <- max(ref_genes7$X307041717)
genes_bin7$density <- (genes_bin7$freq/(genes_bin7$length/1000000))/1000
plot(genes_bin7$foo.X1/1000000, genes_bin7$density,type="l")
saveRDS(genes_bin7$foo.X1/1000000,file="chr7_locipos")
saveRDS(genes_bin7$density,file="chr7_geneProb")

genes_bin8 <- as.data.frame(summary(binning(ref_genes8$X1, nbins = round(max(ref_genes8$X307041717)/100000), type = "kmeans")))
genes_bin8 <- within(genes_bin8, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin8$levels), ',', fixed=TRUE))))
genes_bin8 <- do.call(data.frame, genes_bin8)
genes_bin8 <- genes_bin8 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin8 <- genes_bin8 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin8[1,4]<-17350
genes_bin8$length <- (genes_bin8$foo.X2-genes_bin8$foo.X1)
genes_bin8[round(max(ref_genes8$X307041717)/100000),5] <- max(ref_genes8$X307041717)
genes_bin8$density <- (genes_bin8$freq/(genes_bin8$length/1000000))/1000
plot(genes_bin8$foo.X1/1000000, genes_bin8$density,type="l")
saveRDS(genes_bin8$foo.X1/1000000,file="chr8_locipos")
saveRDS(genes_bin8$density,file="chr8_geneProb")

genes_bin9 <- as.data.frame(summary(binning(ref_genes9$X1, nbins = round(max(ref_genes9$X307041717)/100000), type = "kmeans")))
genes_bin9 <- within(genes_bin9, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin9$levels), ',', fixed=TRUE))))
genes_bin9 <- do.call(data.frame, genes_bin9)
genes_bin9 <- genes_bin9 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin9 <- genes_bin9 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin9[1,4]<-145180
genes_bin9$length <- (genes_bin9$foo.X2-genes_bin9$foo.X1)
genes_bin9[round(max(ref_genes9$X307041717)/100000),5] <- max(ref_genes9$X307041717)
genes_bin9$density <- (genes_bin9$freq/(genes_bin9$length/1000000))/1000
plot(genes_bin9$foo.X1/1000000, genes_bin9$density,type="l")
saveRDS(genes_bin9$foo.X1/1000000,file="chr9_locipos")
saveRDS(genes_bin9$density,file="chr9_geneProb")

genes_bin10 <- binning(ref_genes10$X1, nbins = round(max(ref_genes10$X307041717)/100000), type = "kmeans")
genes_bin10 <- as.data.frame(summary(binning(ref_genes10$X1, nbins = round(max(ref_genes10$X307041717)/100000), type = "kmeans")))
genes_bin10 <- within(genes_bin10, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin10$levels), ',', fixed=TRUE))))
genes_bin10 <- do.call(data.frame, genes_bin10)
genes_bin10 <- genes_bin10 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin10 <- genes_bin10 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin10[1,4]<-44901
genes_bin10$length <- (genes_bin10$foo.X2-genes_bin10$foo.X1)
genes_bin10[round(max(ref_genes10$X307041717)/100000),5] <- max(ref_genes10$X307041717)
genes_bin10$density <- (genes_bin10$freq/(genes_bin10$length/1000000))/1000
plot(genes_bin10$foo.X1/1000000, genes_bin10$density,type="l")
saveRDS(genes_bin10$foo.X1/1000000,file="chr10_locipos")
saveRDS(genes_bin10$density,file="chr10_geneProb")

genes_bin11 <- binning(ref_genes11$X1, nbins = round(max(ref_genes11$X307041717)/100000), type = "kmeans")
genes_bin11 <- as.data.frame(summary(binning(ref_genes11$X1, nbins = round(max(ref_genes11$X307041717)/100000), type = "kmeans")))
genes_bin11 <- within(genes_bin11, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin11$levels), ',', fixed=TRUE))))
genes_bin11 <- do.call(data.frame, genes_bin11)
genes_bin11 <- genes_bin11 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin11 <- genes_bin11 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin11[1,4]<-1179
genes_bin11$length <- (genes_bin11$foo.X2-genes_bin11$foo.X1)
genes_bin11[round(max(ref_genes11$X307041717)/100000),5] <- max(ref_genes11$X307041717)
genes_bin11$density <- (genes_bin11$freq/(genes_bin11$length/1000000))/1000
plot(genes_bin11$foo.X1/1000000, genes_bin11$density,type="l")
saveRDS(genes_bin11$foo.X1/1000000,file="chr11_locipos")
saveRDS(genes_bin11$density,file="chr11_geneProb")

genes_bin12 <- binning(ref_genes12$X1, nbins = round(max(ref_genes12$X307041717)/100000), type = "kmeans")
genes_bin12 <- as.data.frame(summary(binning(ref_genes12$X1, nbins = round(max(ref_genes12$X307041717)/100000), type = "kmeans")))
genes_bin12 <- within(genes_bin12, foo<-data.frame(do.call('rbind', strsplit(as.character(genes_bin12$levels), ',', fixed=TRUE))))
genes_bin12 <- do.call(data.frame, genes_bin12)
genes_bin12 <- genes_bin12 %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
genes_bin12 <- genes_bin12 %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
genes_bin12[1,4]<-2681
genes_bin12$length <- (genes_bin12$foo.X2-genes_bin12$foo.X1)
genes_bin12[round(max(ref_genes12$X307041717)/100000),5] <- max(ref_genes12$X307041717)
genes_bin12$density <- (genes_bin12$freq/(genes_bin12$length/1000000))/1000
plot(genes_bin12$foo.X1/1000000, genes_bin12$density,type="l")
saveRDS(genes_bin12$foo.X1/1000000,file="chr12_locipos")
saveRDS(genes_bin12$density,file="chr12_geneProb")


##Conduct Spearman's rank correlation test to calculate correlation between gene density & recombination rate
#Create dataframes with gene bins and associated gene density
chr1_corr <- genes_bin1$density
chr1_corr <- as.data.frame(chr1_corr)
colnames(chr1_corr) <- "density"
chr1_corr$start <- genes_bin1$foo.X1
chr1_corr$end <- genes_bin1$foo.X2
chr1_corr$rate <- NA

chr2_corr <- genes_bin2$density
chr2_corr <- as.data.frame(chr2_corr)
colnames(chr2_corr) <- "density"
chr2_corr$start <- genes_bin2$foo.X1
chr2_corr$end <- genes_bin2$foo.X2
chr2_corr$rate <- NA

chr3_corr <- genes_bin3$density
chr3_corr <- as.data.frame(chr3_corr)
colnames(chr3_corr) <- "density"
chr3_corr$start <- genes_bin3$foo.X1
chr3_corr$end <- genes_bin3$foo.X2
chr3_corr$rate <- NA

chr4_corr <- genes_bin4$density
chr4_corr <- as.data.frame(chr4_corr)
colnames(chr4_corr) <- "density"
chr4_corr$start <- genes_bin4$foo.X1
chr4_corr$end <- genes_bin4$foo.X2
chr4_corr$rate <- NA

chr5_corr <- genes_bin5$density
chr5_corr <- as.data.frame(chr5_corr)
colnames(chr5_corr) <- "density"
chr5_corr$start <- genes_bin5$foo.X1
chr5_corr$end <- genes_bin5$foo.X2
chr5_corr$rate <- NA

chr6_corr <- genes_bin6$density
chr6_corr <- as.data.frame(chr6_corr)
colnames(chr6_corr) <- "density"
chr6_corr$start <- genes_bin6$foo.X1
chr6_corr$end <- genes_bin6$foo.X2
chr6_corr$rate <- NA

chr7_corr <- genes_bin7$density
chr7_corr <- as.data.frame(chr7_corr)
colnames(chr7_corr) <- "density"
chr7_corr$start <- genes_bin7$foo.X1
chr7_corr$end <- genes_bin7$foo.X2
chr7_corr$rate <- NA

chr8_corr <- genes_bin8$density
chr8_corr <- as.data.frame(chr8_corr)
colnames(chr8_corr) <- "density"
chr8_corr$start <- genes_bin8$foo.X1
chr8_corr$end <- genes_bin8$foo.X2
chr8_corr$rate <- NA

chr9_corr <- genes_bin9$density
chr9_corr <- as.data.frame(chr9_corr)
colnames(chr9_corr) <- "density"
chr9_corr$start <- genes_bin9$foo.X1
chr9_corr$end <- genes_bin9$foo.X2
chr9_corr$rate <- NA

chr10_corr <- genes_bin10$density
chr10_corr <- as.data.frame(chr10_corr)
colnames(chr10_corr) <- "density"
chr10_corr$start <- genes_bin10$foo.X1
chr10_corr$end <- genes_bin10$foo.X2
chr10_corr$rate <- NA

chr11_corr <- genes_bin11$density
chr11_corr <- as.data.frame(chr11_corr)
colnames(chr11_corr) <- "density"
chr11_corr$start <- genes_bin11$foo.X1
chr11_corr$end <- genes_bin11$foo.X2
chr11_corr$rate <- NA

chr12_corr <- genes_bin12$density
chr12_corr <- as.data.frame(chr12_corr)
colnames(chr12_corr) <- "density"
chr12_corr$start <- genes_bin12$foo.X1
chr12_corr$end <- genes_bin12$foo.X2
chr12_corr$rate <- NA

#Function to fill in the recombination rates associated with the gene bin
assign_rate <- function(chr_bin, chr_corr){
  for(i in 1:nrow(chr_bin)){
    for(k in 1:nrow(chr_corr)){
      if(isTRUE((chr_bin$`CO Start`[i] <= chr_corr$start[k]) && (chr_bin$`CO End`[i] >= chr_corr$end[k]))){
        chr_corr$rate[k] <- chr_bin$rate[i]
      }
    }
  }
  return(chr_corr)
}
#Correlation test between binned recombination rates and gene density
chr1_corr_rate <- assign_rate(japonica_chr1_CO_2, chr1_corr)
chr1_corr_rate <- na.omit(chr1_corr_rate)
cor.test(chr1_corr_rate$rate, chr1_corr_rate$density,  method = "spearman", alternative = "greater")

chr2_corr_rate <- assign_rate(japonica_chr2_CO_2, chr2_corr)
chr2_corr_rate <-  na.omit(chr2_corr_rate)
cor.test(chr2_corr_rate$rate, chr2_corr_rate$density,  method = "spearman", alternative = "greater")

chr3_corr_rate <- assign_rate(japonica_chr3_CO_2, chr3_corr)
chr3_corr_rate <- na.omit(chr3_corr_rate)
cor.test(chr3_corr_rate$rate, chr3_corr_rate$density,  method = "spearman", alternative = "greater")

chr4_corr_rate <- assign_rate(japonica_chr4_CO_2, chr4_corr)
chr4_corr_rate <- na.omit(chr4_corr_rate)
cor.test(chr4_corr_rate$rate, chr4_corr_rate$density,  method = "spearman", alternative = "greater")

chr5_corr_rate <- assign_rate(japonica_chr5_CO_2, chr5_corr)
chr5_corr_rate <- na.omit(chr5_corr_rate)
cor.test(chr5_corr_rate$rate, chr5_corr_rate$density,  method = "spearman", alternative = "greater")

chr6_corr_rate <- assign_rate(japonica_chr6_CO_2, chr6_corr)
chr6_corr_rate <- na.omit(chr6_corr_rate)
cor.test(chr6_corr_rate$rate, chr6_corr_rate$density,  method = "spearman", alternative = "greater")

chr7_corr_rate <- assign_rate(japonica_chr7_CO_2, chr7_corr)
chr7_corr_rate <-na.omit(chr7_corr_rate)
cor.test(chr7_corr_rate$rate, chr7_corr_rate$density,  method = "spearman", alternative = "greater")

chr8_corr_rate <- assign_rate(japonica_chr8_CO_2, chr8_corr)
chr8_corr_rate <- na.omit(chr8_corr_rate)
cor.test(chr8_corr_rate$rate, chr8_corr_rate$density,  method = "spearman", alternative = "greater")

chr9_corr_rate <- assign_rate(japonica_chr9_CO_2, chr9_corr)
chr9_corr_rate <- na.omit(chr9_corr_rate)
cor.test(chr9_corr_rate$rate, chr9_corr_rate$density,  method = "spearman", alternative = "greater")

chr10_corr_rate <- assign_rate(japonica_chr10_CO_2, chr10_corr)
chr10_corr_rate <- na.omit(chr10_corr_rate)
cor.test(chr10_corr_rate$rate, chr10_corr_rate$density,  method = "spearman", alternative = "greater")

chr11_corr_rate <- assign_rate(japonica_chr11_CO_2, chr11_corr)
chr11_corr_rate <- na.omit(chr11_corr_rate)
cor.test(chr11_corr_rate$rate, chr11_corr_rate$density,  method = "spearman", alternative = "greater")

chr12_corr_rate <- assign_rate(japonica_chr12_CO_2, chr12_corr)
chr12_corr_rate <- na.omit(chr12_corr_rate)
cor.test(chr12_corr_rate$rate, chr12_corr_rate$density,  method = "spearman", alternative = "greater")

#Create dataframe of all chromosomes with gene bins and gene density, and associated average recombination rate 
genomewide <- rbind(chr1_corr_rate, chr2_corr_rate, chr3_corr_rate, chr4_corr_rate, chr5_corr_rate,
                    chr6_corr_rate, chr7_corr_rate, chr8_corr_rate, chr9_corr_rate, chr10_corr_rate, chr11_corr_rate, chr12_corr_rate)

#Assess the genome wide correlation between gene density and recombination rate
cor.test(genomewide$rate, genomewide$density,  method = "spearman", alternative = "greater")
