library(AlphaSimR)
library(Rcpp)
library(ggplot2)
library(dplyr)

setwd("C:/Users/16192/Documents/PNAS_Simulations")
set.seed(420)

japonica_snps <- read.table("japonica_SNPs.bed")
indica_snps <- read.table("indica_snps.bed")

