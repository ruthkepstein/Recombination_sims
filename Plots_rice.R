#making recombination distribution & genetic map plots

##this function scales graph
library(scales)
squish_trans <- function(from, to, factor) {
  
  trans <- function(x) {
    
    if (any(is.na(x))) return(x)
    
    # get indices for the relevant regions
    isq <- x > from & x < to
    ito <- x >= to
    
    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    
    return(x)
  }
  
  inv <- function(x) {
    
    if (any(is.na(x))) return(x)
    
    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor
    ito <- x >= from + (to - from)/factor
    
    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    
    return(x)
  }
  
  # return the transformation
  return(trans_new("squished", trans, inv))
}

group.colors = c("ddm1" = "#E69F00", "recq4" = "#56B4E9", "wt" = "#009E73", "zmet2" = "#F0E442",
                 "fancm" = "#0072B2", "10X" = "#D55E00", "ddm1/zmet2" = "#CC79A7")
chr1_snp2ddm1<- readRDS("chr1_snp2ddm1.RData")
chr2_snp2ddm1<- readRDS("chr2_snp2ddm1.RData")
chr3_snp2ddm1<- readRDS("chr3_snp2ddm1.RData")
chr4_snp2ddm1<- readRDS("chr4_snp2ddm1.RData")
chr5_snp2ddm1<- readRDS("chr5_snp2ddm1.RData")
chr6_snp2ddm1<- readRDS("chr6_snp2ddm1.RData")
chr7_snp2ddm1<- readRDS("chr7_snp2ddm1.RData")
chr8_snp2ddm1<- readRDS("chr8_snp2ddm1.RData")
chr9_snp2ddm1<- readRDS("chr9_snp2ddm1.RData")
chr10_snp2ddm1<- readRDS("chr10_snp2ddm1.RData")
chr11_snp2ddm1<- readRDS("chr11_snp2ddm1.RData")
chr12_snp2ddm1<- readRDS("chr12_snp2ddm1.RData")

chr1_snp2<- readRDS("chr1_snp2jap.RData")
chr2_snp2<- readRDS("chr2_snp2jap.RData")
chr3_snp2<- readRDS("chr3_snp2jap.RData")
chr4_snp2<- readRDS("chr4_snp2jap.RData")
chr5_snp2<- readRDS("chr5_snp2jap.RData")
chr6_snp2<- readRDS("chr6_snp2jap.RData")
chr7_snp2<- readRDS("chr7_snp2jap.RData")
chr8_snp2<- readRDS("chr8_snp2jap.RData")
chr9_snp2<- readRDS("chr9_snp2jap.RData")
chr10_snp2<- readRDS("chr10_snp2jap.RData")
chr11_snp2<- readRDS("chr11_snp2jap.RData")
chr12_snp2<- readRDS("chr12_snp2jap.RData")

chr1_snp2zmet2<- readRDS("chr1_snp2zmet2.RData")
chr2_snp2zmet2<- readRDS("chr2_snp2zmet2.RData")
chr3_snp2zmet2<- readRDS("chr3_snp2zmet2.RData")
chr4_snp2zmet2<- readRDS("chr4_snp2zmet2.RData")
chr5_snp2zmet2<- readRDS("chr5_snp2zmet2.RData")
chr6_snp2zmet2<- readRDS("chr6_snp2zmet2.RData")
chr7_snp2zmet2<- readRDS("chr7_snp2zmet2.RData")
chr8_snp2zmet2<- readRDS("chr8_snp2zmet2.RData")
chr9_snp2zmet2<- readRDS("chr9_snp2zmet2.RData")
chr10_snp2zmet2<- readRDS("chr10_snp2zmet2.RData")
chr11_snp2zmet2<- readRDS("chr11_snp2zmet2.RData")
chr12_snp2zmet2<- readRDS("chr12_snp2zmet2.RData")

chr1_snp2recq4<- readRDS("chr1_snp2recq4.RData")
chr2_snp2recq4<- readRDS("chr2_snp2recq4.RData")
chr3_snp2recq4<- readRDS("chr3_snp2recq4.RData")
chr4_snp2recq4<- readRDS("chr4_snp2recq4.RData")
chr5_snp2recq4<- readRDS("chr5_snp2recq4.RData")
chr6_snp2recq4<- readRDS("chr6_snp2recq4.RData")
chr7_snp2recq4<- readRDS("chr7_snp2recq4.RData")
chr8_snp2recq4<- readRDS("chr8_snp2recq4.RData")
chr9_snp2recq4<- readRDS("chr9_snp2recq4.RData")
chr10_snp2recq4<- readRDS("chr10_snp2recq4.RData")
chr11_snp2recq4<- readRDS("chr11_snp2recq4.RData")
chr12_snp2recq4<- readRDS("chr12_snp2recq4.RData")

chr1_snp2fancm<- readRDS("chr1_snp2fancm.RData")
chr2_snp2fancm<- readRDS("chr2_snp2fancm.RData")
chr3_snp2fancm<- readRDS("chr3_snp2fancm.RData")
chr4_snp2fancm<- readRDS("chr4_snp2fancm.RData")
chr5_snp2fancm<- readRDS("chr5_snp2fancm.RData")
chr6_snp2fancm<- readRDS("chr6_snp2fancm.RData")
chr7_snp2fancm<- readRDS("chr7_snp2fancm.RData")
chr8_snp2fancm<- readRDS("chr8_snp2fancm.RData")
chr9_snp2fancm<- readRDS("chr9_snp2fancm.RData")
chr10_snp2fancm<- readRDS("chr10_snp2fancm.RData")
chr11_snp2fancm<- readRDS("chr11_snp2fancm.RData")
chr12_snp2fancm<- readRDS("chr12_snp2fancm.RData")

chr1_snp2ideal1<- readRDS("chr1_snp2ideal1.RData")
chr2_snp2ideal1<- readRDS("chr2_snp2ideal1.RData")
chr3_snp2ideal1<- readRDS("chr3_snp2ideal1.RData")
chr4_snp2ideal1<- readRDS("chr4_snp2ideal1.RData")
chr5_snp2ideal1<- readRDS("chr5_snp2ideal1.RData")
chr6_snp2ideal1<- readRDS("chr6_snp2ideal1.RData")
chr7_snp2ideal1<- readRDS("chr7_snp2ideal1.RData")
chr8_snp2ideal1<- readRDS("chr8_snp2ideal1.RData")
chr9_snp2ideal1<- readRDS("chr9_snp2ideal1.RData")
chr10_snp2ideal1<- readRDS("chr10_snp2ideal1.RData")
chr11_snp2ideal1<- readRDS("chr11_snp2ideal1.RData")
chr12_snp2ideal1<- readRDS("chr12_snp2ideal1.RData")

chr1_snp2ideal2<- readRDS("chr1_snp2ideal2.RData")
chr2_snp2ideal2<- readRDS("chr2_snp2ideal2.RData")
chr3_snp2ideal2<- readRDS("chr3_snp2ideal2.RData")
chr4_snp2ideal2<- readRDS("chr4_snp2ideal2.RData")
chr5_snp2ideal2<- readRDS("chr5_snp2ideal2.RData")
chr6_snp2ideal2<- readRDS("chr6_snp2ideal2.RData")
chr7_snp2ideal2<- readRDS("chr7_snp2ideal2.RData")
chr8_snp2ideal2<- readRDS("chr8_snp2ideal2.RData")
chr9_snp2ideal2<- readRDS("chr9_snp2ideal2.RData")
chr10_snp2ideal2<- readRDS("chr10_snp2ideal2.RData")
chr11_snp2ideal2<- readRDS("chr11_snp2ideal2.RData")
chr12_snp2ideal2<- readRDS("chr12_snp2ideal2.RData")

##Plots for genetic maps for each mutant & chromosome
ggplot(chr1_snp2, aes(x = `SNP Start`, y = pos)) + theme_bw() +
  geom_line(aes(x = `SNP Start`, y = pos, color = "wt"), size = 2) +
  geom_line(data = chr1_snp2ddm1, aes(x = `SNP Start`, y = pos, color = "ddm1"), size = 2) +
  geom_line(data = chr1_snp2zmet2, aes(x = `SNP Start`, y = pos, color = "zmet2"), size = 2) +
  geom_line(data = chr1_snp2recq4, aes(x = `SNP Start`, y = pos, color = "recq4"), size = 2) +
  geom_line(data = chr1_snp2fancm, aes(x = `SNP Start`, y = pos, color = "fancm"), size = 2) +
  geom_line(data = chr1_snp2ideal1, aes(x = `SNP Start`, y = pos, color = "10X"), size = 2) +
  geom_line(data = chr1_snp2ideal2, aes(x = `SNP Start`, y = pos, color = "ddm1/zmet2"), size = 2) +
  scale_color_manual(values=c("ddm1" = "#E69F00", "recq4" = "#56B4E9", "wt" = "#009E73", "zmet2" = "#F0E442",
                              "fancm" = "#0072B2", "10X" = "#D55E00", "ddm1/zmet2" = "#CC79A7")) + labs(y = "Accumulative Genetic Pos (cM)",x = "Physical Pos (Mb)",color = "Genetic Map", title="Chromosome 1 Genetic Map") +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16), legend.title = element_text(size=11), #change legend title font size
        legend.text = element_text(size=11), plot.title = element_text(size=15, face="bold"), legend.position = c(.15,.8), legend.key.size = unit(0.3, "lines")) + annotate("rect", xmin = 16, xmax = 17, ymin = 2.5, ymax = 200, fill='blue', alpha = .2) + 
  annotate('text', x = 16.7, y = 200, label = 'Centromere', color='blue', size = 4) + scale_y_continuous(trans = squish_trans(60,215,8), breaks = c(50, 100, 150, ceiling(max(chr1_snp2ideal1$pos))))

ggplot(chr2_snp2, aes(x = `SNP Start`, y = pos)) + theme_bw() +
  geom_line(aes(x = `SNP Start`, y = pos, color = "wt"), size = 2) +
  geom_line(data = chr2_snp2ddm1, aes(x = `SNP Start`, y = pos, color = "ddm1"), size = 2) +
  geom_line(data = chr2_snp2zmet2, aes(x = `SNP Start`, y = pos, color = "zmet2"), size = 2) +
  geom_line(data = chr2_snp2recq4, aes(x = `SNP Start`, y = pos, color = "recq4"), size = 2) +
  geom_line(data = chr2_snp2fancm, aes(x = `SNP Start`, y = pos, color = "fancm"), size = 2) +
  geom_line(data = chr2_snp2ideal1, aes(x = `SNP Start`, y = pos, color = "10X"), size = 2) +
  geom_line(data = chr2_snp2ideal2, aes(x = `SNP Start`, y = pos, color = "ddm1/zmet2"), size = 2) +
  scale_color_manual(values=c("ddm1" = "#E69F00", "recq4" = "#56B4E9", "wt" = "#009E73", "zmet2" = "#F0E442",
                              "fancm" = "#0072B2", "10X" = "#D55E00", "ddm1/zmet2" = "#CC79A7")) + labs(y = "Accumulative Genetic Pos (cM)",x = "Physical Pos (Mb)",color = "Genetic Map", title="Chromosome 2 Genetic Map") +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16), legend.title = element_text(size=11), #change legend title font size
        legend.text = element_text(size=11), plot.title = element_text(size=15, face="bold"), legend.position = c(.11,.8), legend.key.size = unit(0.3, "lines")) + annotate("rect", xmin = 13, xmax = 14, ymin = 2.5, ymax = 195, fill='blue', alpha = .2) + 
  annotate('text', x = 13.6, y = 198, label = 'Centromere', color='blue', size = 4) + scale_y_continuous(trans = squish_trans(60,194,8), breaks = c(50, 100, ceiling(max(chr2_snp2ideal1$pos))))

ggplot(chr3_snp2, aes(x = `SNP Start`, y = pos)) + theme_bw() +
  geom_line(aes(x = `SNP Start`, y = pos, color = "wt"), size = 2) +
  geom_line(data = chr3_snp2ddm1, aes(x = `SNP Start`, y = pos, color = "ddm1"), size = 2) +
  geom_line(data = chr3_snp2zmet2, aes(x = `SNP Start`, y = pos, color = "zmet2"), size = 2) +
  geom_line(data = chr3_snp2recq4, aes(x = `SNP Start`, y = pos, color = "recq4"), size = 2) +
  geom_line(data = chr3_snp2fancm, aes(x = `SNP Start`, y = pos, color = "fancm"), size = 2) +
  geom_line(data = chr3_snp2ideal1, aes(x = `SNP Start`, y = pos, color = "10X"), size = 2) +
  geom_line(data = chr3_snp2ideal2, aes(x = `SNP Start`, y = pos, color = "ddm1/zmet2"), size = 2) +
  scale_color_manual(values=c("ddm1" = "#E69F00", "recq4" = "#56B4E9", "wt" = "#009E73", "zmet2" = "#F0E442",
                              "fancm" = "#0072B2", "10X" = "#D55E00", "ddm1/zmet2" = "#CC79A7")) + labs(y = "Accumulative Genetic Pos (cM)",x = "Physical Pos (Mb)",color = "Genetic Map", title="Chromosome 3 Genetic Map") +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16), legend.title = element_text(size=11), #change legend title font size
        legend.text = element_text(size=11), plot.title = element_text(size=15, face="bold"), legend.position = c(.15,.8), legend.key.size = unit(0.3, "lines")) + annotate("rect", xmin = 19, xmax = 20, ymin = 2.5, ymax =115, fill='blue', alpha = .2) + 
  annotate('text', x = 19.4, y = 130, label = 'Centromere', color='blue', size = 4) + scale_y_continuous(trans = squish_trans(50,164,8), breaks = c(50, 100, ceiling(max(chr3_snp2ideal1$pos))))

ggplot(chr4_snp2, aes(x = `SNP Start`, y = pos)) + theme_bw() +
  geom_line(aes(x = `SNP Start`, y = pos, color = "wt"), size = 2) +
  geom_line(data = chr4_snp2ddm1, aes(x = `SNP Start`, y = pos, color = "ddm1"), size = 2) +
  geom_line(data = chr4_snp2zmet2, aes(x = `SNP Start`, y = pos, color = "zmet2"), size = 2) +
  geom_line(data = chr4_snp2recq4, aes(x = `SNP Start`, y = pos, color = "recq4"), size = 2) +
  geom_line(data = chr4_snp2fancm, aes(x = `SNP Start`, y = pos, color = "fancm"), size = 2) +
  geom_line(data = chr4_snp2ideal1, aes(x = `SNP Start`, y = pos, color = "10X"), size = 2) +
  geom_line(data = chr4_snp2ideal2, aes(x = `SNP Start`, y = pos, color = "ddm1/zmet2"), size = 2) +
  scale_color_manual(values=c("ddm1" = "#E69F00", "recq4" = "#56B4E9", "wt" = "#009E73", "zmet2" = "#F0E442",
                              "fancm" = "#0072B2", "10X" = "#D55E00", "ddm1/zmet2" = "#CC79A7")) + labs(y = "Accumulative Genetic Pos (cM)",x = "Physical Pos (Mb)",color = "Genetic Map", title="Chromosome 4 Genetic Map") +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16), legend.title = element_text(size=11), #change legend title font size
        legend.text = element_text(size=11), plot.title = element_text(size=15, face="bold"), legend.position = c(.1,.82), legend.key.size = unit(0.3, "lines")) + annotate("rect", xmin = 9, xmax = 10, ymin = 2.5, ymax = 185, fill='blue', alpha = .2) + 
  annotate('text', x = 9.7, y = 200, label = 'Centromere', color='blue', size = 4) + scale_y_continuous(trans = squish_trans(60,215,8), breaks = c(50, 100, 150, ceiling(max(chr4_snp2ideal1$pos))))

ggplot(chr5_snp2, aes(x = `SNP Start`, y = pos)) + theme_bw() +
  geom_line(aes(x = `SNP Start`, y = pos, color = "wt"), size = 2) +
  geom_line(data = chr5_snp2ddm1, aes(x = `SNP Start`, y = pos, color = "ddm1"), size = 2) +
  geom_line(data = chr5_snp2zmet2, aes(x = `SNP Start`, y = pos, color = "zmet2"), size = 2) +
  geom_line(data = chr5_snp2recq4, aes(x = `SNP Start`, y = pos, color = "recq4"), size = 2) +
  geom_line(data = chr5_snp2fancm, aes(x = `SNP Start`, y = pos, color = "fancm"), size = 2) +
  geom_line(data = chr5_snp2ideal1, aes(x = `SNP Start`, y = pos, color = "10X"), size = 2) +
  geom_line(data = chr5_snp2ideal2, aes(x = `SNP Start`, y = pos, color = "ddm1/zmet2"), size = 2) +
  scale_color_manual(values=c("ddm1" = "#E69F00", "recq4" = "#56B4E9", "wt" = "#009E73", "zmet2" = "#F0E442",
                              "fancm" = "#0072B2", "10X" = "#D55E00", "ddm1/zmet2" = "#CC79A7")) + labs(y = "Accumulative Genetic Pos (cM)",x = "Physical Pos (Mb)",color = "Genetic Map", title="Chromosome 5 Genetic Map") +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16), legend.title = element_text(size=11), #change legend title font size
        legend.text = element_text(size=11), plot.title = element_text(size=15, face="bold"), legend.position = c(.15,.8), legend.key.size = unit(0.3, "lines")) + annotate("rect", xmin = 12, xmax = 13, ymin = 2.5, ymax = 130, fill='blue', alpha = .2) + 
  annotate('text', x = 12.4, y = 145, label = 'Centromere', color='blue', size = 4) + scale_y_continuous(trans = squish_trans(50,144,8), breaks = c(50, 100,  ceiling(max(chr5_snp2ideal1$pos))))

ggplot(chr6_snp2, aes(x = `SNP Start`, y = pos)) + theme_bw() +
  geom_line(aes(x = `SNP Start`, y = pos, color = "wt"), size = 2) +
  geom_line(data = chr6_snp2ddm1, aes(x = `SNP Start`, y = pos, color = "ddm1"), size = 2) +
  geom_line(data = chr6_snp2zmet2, aes(x = `SNP Start`, y = pos, color = "zmet2"), size = 2) +
  geom_line(data = chr6_snp2recq4, aes(x = `SNP Start`, y = pos, color = "recq4"), size = 2) +
  geom_line(data = chr6_snp2fancm, aes(x = `SNP Start`, y = pos, color = "fancm"), size = 2) +
  geom_line(data = chr6_snp2ideal1, aes(x = `SNP Start`, y = pos, color = "10X"), size = 2) +
  geom_line(data = chr6_snp2ideal2, aes(x = `SNP Start`, y = pos, color = "ddm1/zmet2"), size = 2) +
  scale_color_manual(values=c("ddm1" = "#E69F00", "recq4" = "#56B4E9", "wt" = "#009E73", "zmet2" = "#F0E442",
                              "fancm" = "#0072B2", "10X" = "#D55E00", "ddm1/zmet2" = "#CC79A7")) + labs(y = "Accumulative Genetic Pos (cM)",x = "Physical Pos (Mb)",color = "Genetic Map", title="Chromosome 6 Genetic Map") +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16), legend.title = element_text(size=11), #change legend title font size
        legend.text = element_text(size=11), plot.title = element_text(size=15, face="bold"), legend.position = c(.1,.8), legend.key.size = unit(0.3, "lines")) + annotate("rect", xmin = 15, xmax = 16, ymin = 2.5, ymax = 165, fill='blue', alpha = .2) + 
  annotate('text', x = 15.3, y = 175, label = 'Centromere', color='blue', size = 4) + scale_y_continuous(trans = squish_trans(50,177,8), breaks = c(50, 100, ceiling(max(chr6_snp2ideal1$pos))))

ggplot(chr7_snp2, aes(x = `SNP Start`, y = pos)) + theme_bw() +
  geom_line(aes(x = `SNP Start`, y = pos, color = "wt"), size = 2) +
  geom_line(data = chr7_snp2ddm1, aes(x = `SNP Start`, y = pos, color = "ddm1"), size = 2) +
  geom_line(data = chr7_snp2zmet2, aes(x = `SNP Start`, y = pos, color = "zmet2"), size = 2) +
  geom_line(data = chr7_snp2recq4, aes(x = `SNP Start`, y = pos, color = "recq4"), size = 2) +
  geom_line(data = chr7_snp2fancm, aes(x = `SNP Start`, y = pos, color = "fancm"), size = 2) +
  geom_line(data = chr7_snp2ideal1, aes(x = `SNP Start`, y = pos, color = "10X"), size = 2) +
  geom_line(data = chr7_snp2ideal2, aes(x = `SNP Start`, y = pos, color = "ddm1/zmet2"), size = 2) +
  scale_color_manual(values=c("ddm1" = "#E69F00", "recq4" = "#56B4E9", "wt" = "#009E73", "zmet2" = "#F0E442",
                              "fancm" = "#0072B2", "10X" = "#D55E00", "ddm1/zmet2" = "#CC79A7")) + labs(y = "Accumulative Genetic Pos (cM)",x = "Physical Pos (Mb)",color = "Genetic Map", title="Chromosome 7 Genetic Map") +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16), legend.title = element_text(size=11), #change legend title font size
        legend.text = element_text(size=11), plot.title = element_text(size=15, face="bold"), legend.position = c(.1,.8), legend.key.size = unit(0.3, "lines")) + annotate("rect", xmin = 12, xmax = 13, ymin = 2.5, ymax = 170, fill='blue', alpha = .2) + 
  annotate('text', x = 12.1, y = 178, label = 'Centromere', color='blue', size = 4) + scale_y_continuous(trans = squish_trans(50,175,8), breaks = c(50, 100, ceiling(max(chr7_snp2ideal1$pos))))

ggplot(chr8_snp2, aes(x = `SNP Start`, y = pos)) + theme_bw() +
  geom_line(aes(x = `SNP Start`, y = pos, color = "wt"), size = 2) +
  geom_line(data = chr8_snp2ddm1, aes(x = `SNP Start`, y = pos, color = "ddm1"), size = 2) +
  geom_line(data = chr8_snp2zmet2, aes(x = `SNP Start`, y = pos, color = "zmet2"), size = 2) +
  geom_line(data = chr8_snp2recq4, aes(x = `SNP Start`, y = pos, color = "recq4"), size = 2) +
  geom_line(data = chr8_snp2fancm, aes(x = `SNP Start`, y = pos, color = "fancm"), size = 2) +
  geom_line(data = chr8_snp2ideal1, aes(x = `SNP Start`, y = pos, color = "10X"), size = 2) +
  geom_line(data = chr8_snp2ideal2, aes(x = `SNP Start`, y = pos, color = "ddm1/zmet2"), size = 2) +
  scale_color_manual(values=c("ddm1" = "#E69F00", "recq4" = "#56B4E9", "wt" = "#009E73", "zmet2" = "#F0E442",
                              "fancm" = "#0072B2", "10X" = "#D55E00", "ddm1/zmet2" = "#CC79A7")) + labs(y = "Accumulative Genetic Pos (cM)",x = "Physical Pos (Mb)",color = "Genetic Map", title="Chromosome 8 Genetic Map") +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16), legend.title = element_text(size=11), #change legend title font size
        legend.text = element_text(size=11), plot.title = element_text(size=15, face="bold"), legend.position = c(.1,.8), legend.key.size = unit(0.3, "lines")) + annotate("rect", xmin = 12, xmax = 13, ymin = 2.5, ymax = 140, fill='blue', alpha = .2) + 
  annotate('text', x = 12.9, y = 150, label = 'Centromere', color='blue', size = 4) + scale_y_continuous(trans = squish_trans(50,150,8), breaks = c(50, 100, ceiling(max(chr8_snp2ideal1$pos))))

ggplot(chr9_snp2, aes(x = `SNP Start`, y = pos)) + theme_bw() +
  geom_line(aes(x = `SNP Start`, y = pos, color = "wt"), size = 2) +
  geom_line(data = chr9_snp2ddm1, aes(x = `SNP Start`, y = pos, color = "ddm1"), size = 2) +
  geom_line(data = chr9_snp2zmet2, aes(x = `SNP Start`, y = pos, color = "zmet2"), size = 2) +
  geom_line(data = chr9_snp2recq4, aes(x = `SNP Start`, y = pos, color = "recq4"), size = 2) +
  geom_line(data = chr9_snp2fancm, aes(x = `SNP Start`, y = pos, color = "fancm"), size = 2) +
  geom_line(data = chr9_snp2ideal1, aes(x = `SNP Start`, y = pos, color = "10X"), size = 2) +
  geom_line(data = chr9_snp2ideal2, aes(x = `SNP Start`, y = pos, color = "ddm1/zmet2"), size = 2) +
  scale_color_manual(values=c("ddm1" = "#E69F00", "recq4" = "#56B4E9", "wt" = "#009E73", "zmet2" = "#F0E442",
                              "fancm" = "#0072B2", "10X" = "#D55E00", "ddm1/zmet2" = "#CC79A7")) + labs(y = "Accumulative Genetic Pos (cM)",x = "Physical Pos (Mb)",color = "Genetic Map", title="Chromosome 9 Genetic Map") +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16), legend.title = element_text(size=11), #change legend title font size
        legend.text = element_text(size=11), plot.title = element_text(size=15, face="bold"), legend.position = c(.8,.7), legend.key.size = unit(0.3, "lines")) + annotate("rect", xmin = 2, xmax = 3, ymin = 2.5, ymax = 140, fill='blue', alpha = .2) + 
  annotate('text', x = 2.8, y = 150, label = 'Centromere', color='blue', size = 4) + scale_y_continuous(trans = squish_trans(50,151,8), breaks = c(50, 100, ceiling(max(chr9_snp2ideal1$pos))))

ggplot(chr10_snp2, aes(x = `SNP Start`, y = pos)) + theme_bw() +
  geom_line(aes(x = `SNP Start`, y = pos, color = "wt"), size = 2) +
  geom_line(data = chr10_snp2ddm1, aes(x = `SNP Start`, y = pos, color = "ddm1"), size = 2) +
  geom_line(data = chr10_snp2zmet2, aes(x = `SNP Start`, y = pos, color = "zmet2"), size = 2) +
  geom_line(data = chr10_snp2recq4, aes(x = `SNP Start`, y = pos, color = "recq4"), size = 2) +
  geom_line(data = chr10_snp2fancm, aes(x = `SNP Start`, y = pos, color = "fancm"), size = 2) +
  geom_line(data = chr10_snp2ideal1, aes(x = `SNP Start`, y = pos, color = "10X"), size = 2) +
  geom_line(data = chr10_snp2ideal2, aes(x = `SNP Start`, y = pos, color = "ddm1/zmet2"), size = 2) +
  scale_color_manual(values=c("ddm1" = "#E69F00", "recq4" = "#56B4E9", "wt" = "#009E73", "zmet2" = "#F0E442",
                              "fancm" = "#0072B2", "10X" = "#D55E00", "ddm1/zmet2" = "#CC79A7")) + labs(y = "Accumulative Genetic Pos (cM)",x = "Physical Pos (Mb)",color = "Genetic Map", title="Chromosome 10 Genetic Map") +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16), legend.title = element_text(size=11), #change legend title font size
        legend.text = element_text(size=11), plot.title = element_text(size=15, face="bold"), legend.position = c(.1,.8), legend.key.size = unit(0.3, "lines")) + annotate("rect", xmin =8, xmax = 9, ymin = 2.5, ymax = 115, fill='blue', alpha = .2) + 
  annotate('text', x = 8.2, y = 128, label = 'Centromere', color='blue', size = 4) + scale_y_continuous(trans = squish_trans(40,128,8), breaks = c(25, 50, ceiling(max(chr10_snp2ideal1$pos))))

ggplot(chr11_snp2, aes(x = `SNP Start`, y = pos)) + theme_bw() +
  geom_line(aes(x = `SNP Start`, y = pos, color = "wt"), size = 2) +
  geom_line(data = chr11_snp2ddm1, aes(x = `SNP Start`, y = pos, color = "ddm1"), size = 2) +
  geom_line(data = chr11_snp2zmet2, aes(x = `SNP Start`, y = pos, color = "zmet2"), size = 2) +
  geom_line(data = chr11_snp2recq4, aes(x = `SNP Start`, y = pos, color = "recq4"), size = 2) +
  geom_line(data = chr11_snp2fancm, aes(x = `SNP Start`, y = pos, color = "fancm"), size = 2) +
  geom_line(data = chr11_snp2ideal1, aes(x = `SNP Start`, y = pos, color = "10X"), size = 2) +
  geom_line(data = chr11_snp2ideal2, aes(x = `SNP Start`, y = pos, color = "ddm1/zmet2"), size = 2) +
  scale_color_manual(values=c("ddm1" = "#E69F00", "recq4" = "#56B4E9", "wt" = "#009E73", "zmet2" = "#F0E442",
                              "fancm" = "#0072B2", "10X" = "#D55E00", "ddm1/zmet2" = "#CC79A7")) + labs(y = "Accumulative Genetic Pos (cM)",x = "Physical Pos (Mb)",color = "Genetic Map", title="Chromosome 11 Genetic Map") +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16), legend.title = element_text(size=11), #change legend title font size
        legend.text = element_text(size=11), plot.title = element_text(size=15, face="bold"), legend.position = c(.1,.8), legend.key.size = unit(0.3, "lines")) + annotate("rect", xmin = 11.5, xmax = 12.5, ymin = 2.5, ymax = 175, fill='blue', alpha = .2) + 
  annotate('text', x = 12, y = 187, label = 'Centromere', color='blue', size = 4) + scale_y_continuous(trans = squish_trans(55,187,8), breaks = c(50, 100, ceiling(max(chr11_snp2ideal1$pos))))

ggplot(chr12_snp2, aes(x = `SNP Start`, y = pos)) + theme_bw() +
  geom_line(aes(x = `SNP Start`, y = pos, color = "wt"), size = 2) +
  geom_line(data = chr12_snp2ddm1, aes(x = `SNP Start`, y = pos, color = "ddm1"), size = 2) +
  geom_line(data = chr12_snp2zmet2, aes(x = `SNP Start`, y = pos, color = "zmet2"), size = 2) +
  geom_line(data = chr12_snp2recq4, aes(x = `SNP Start`, y = pos, color = "recq4"), size = 2) +
  geom_line(data = chr12_snp2fancm, aes(x = `SNP Start`, y = pos, color = "fancm"), size = 2) +
  geom_line(data = chr12_snp2ideal1, aes(x = `SNP Start`, y = pos, color = "10X"), size = 2) +
  geom_line(data = chr12_snp2ideal2, aes(x = `SNP Start`, y = pos, color = "ddm1/zmet2"), size = 2) +
  scale_color_manual(values=c("ddm1" = "#E69F00", "recq4" = "#56B4E9", "wt" = "#009E73", "zmet2" = "#F0E442",
                              "fancm" = "#0072B2", "10X" = "#D55E00", "ddm1/zmet2" = "#CC79A7")) + labs(y = "Accumulative Genetic Pos (cM)",x = "Physical Pos (Mb)",color = "Genetic Map", title="Chromosome 12 Genetic Map") +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16), legend.title = element_text(size=11), #change legend title font size
        legend.text = element_text(size=11), plot.title = element_text(size=15, face="bold"), legend.position = c(.1,.8), legend.key.size = unit(0.3, "lines")) + annotate("rect", xmin = 11, xmax = 12, ymin = 2.5, ymax = 165, fill='blue', alpha = .2) + 
  annotate('text', x = 11.9, y = 174, label = 'Centromere', color='blue', size = 4) + scale_y_continuous(trans = squish_trans(50,174,8), breaks = c(50, 100, ceiling(max(chr12_snp2ideal1$pos))))

##This plot makes recombination rate graph
ggplot(chr1_snp2, aes(x = `SNP Start`, y = pos2/`SNP Start`, color = "wt")) + geom_line(size = 2) + theme_bw() +
  geom_line(data = chr1_snp2ddm1, aes(x = `SNP Start`, y = pos2/`SNP Start`, color = "ddm1"), size = 2) +
  geom_line(data = chr1_snp2zmet2, aes(x = `SNP Start`, y = pos2/`SNP Start`, color = "zmet2"), size = 2) +
  geom_line(data = chr1_snp2recq4, aes(x = `SNP Start`, y = pos2/`SNP Start`, color = "recq4"), size = 2) +
  geom_line(data = chr1_snp2fancm, aes(x = `SNP Start`, y = pos2/`SNP Start`, color = "fancm"), size = 2) +
  geom_line(data = chr1_snp2ideal1, aes(x = `SNP Start`, y = pos2/`SNP Start`, color = "10X"), size = 2) +
  geom_line(data = chr1_snp2ideal2, aes(x = `SNP Start`, y = pos2/`SNP Start`, color = "ddm1/zmet2"), size = 2) +
  scale_fill_manual(values=group.colors) + labs(y = "Recombination Rate (cM/Mb)",x = "Physical Pos (Mb)",color = "Genetic Map") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=14), legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14), plot.title = element_text(size=20), legend.position = c(0.2,0.8)) + annotate("rect", xmin = 121, xmax = 133, ymin = 0, ymax = 15, fill='blue', alpha = .2) +
  scale_y_continuous(trans = squish_trans(40,140,100), breaks = c(10,40)) + annotate('text', x = 126, y = 15, 
                                                                                     label = 'Centromere', color='blue',
                                                                                     size = 4) + xlim(0,299)