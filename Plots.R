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

##Plots for genetic maps for each mutant & chromosome
ggplot(chr1_snp2, aes(x = `SNP Start`, y = log(pos))) + theme_bw() +
  geom_line(aes(x = `SNP Start`, y = log(pos), color = "wt"), size = 2) +
  geom_line(data = chr1_snp2ddm1, aes(x = `SNP Start`, y = log(pos), color = "ddm1"), size = 2) +
  geom_line(data = chr1_snp2zmet2, aes(x = `SNP Start`, y = log(pos), color = "zmet2"), size = 2) +
  geom_line(data = chr1_snp2recq4, aes(x = `SNP Start`, y = log(pos), color = "recq4"), size = 2) +
  geom_line(data = chr1_snp2fancm, aes(x = `SNP Start`, y = log(pos), color = "fancm"), size = 2) +
  geom_line(data = chr1_snp2ideal1, aes(x = `SNP Start`, y = log(pos), color = "10X"), size = 2) +
  geom_line(data = chr1_snp2ideal2, aes(x = `SNP Start`, y = log(pos), color = "ddm1/zmet2"), size = 2) +
  scale_color_manual(values=c("ddm1" = "#E69F00", "recq4" = "#56B4E9", "wt" = "#009E73", "zmet2" = "#F0E442",
                              "fancm" = "#0072B2", "10X" = "#D55E00", "ddm1/zmet2" = "#CC79A7")) + labs(y = "Log Transformation of Accumulative Genetic Pos (cM)",x = "Physical Pos (Mb)",color = "Genetic Map") +
  theme(axis.text=element_text(size=16),
        
        axis.title=element_text(size=16), legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=12), plot.title = element_text(size=12), legend.position = c(0.8,0.35), legend.key.size = unit(0.5, "lines")) + ylim(2.5,8) + annotate("rect", xmin = 121, xmax = 133, ymin = 2.5, ymax = 7.8, fill='blue', alpha = .2) + 
  annotate('text', x = 126, y = 7.9, label = 'Centromere', color='blue', size = 4) + guides(color = guide_legend(override.aes = list(size = 2)))

ggplot(chr2_snp2, aes(x = `SNP Start`, y = pos, color = "wt")) + geom_line(size = 2) + theme_bw() +
  geom_line(data = chr2_snp2ddm1, aes(x = `SNP Start`, y = pos, color = "ddm1"), size = 2) +
  geom_line(data = chr2_snp2zmet2, aes(x = `SNP Start`, y = pos, color = "zmet2"), size = 2) +
  geom_line(data = chr2_snp2recq4, aes(x = `SNP Start`, y = pos, color = "recq4"), size = 2) +
  geom_line(data = chr2_snp2fancm, aes(x = `SNP Start`, y = pos, color = "fancm"), size = 2) +
  geom_line(data = chr2_snp2ideal1, aes(x = `SNP Start`, y = pos, color = "10X"), size = 2) +
  geom_line(data = chr2_snp2ideal2, aes(x = `SNP Start`, y = pos, color = "ddm1/zmet2"), size = 2) +
  scale_color_manual(values=group.colors) + labs(y = "Accumulative Genetic Pos (cM)",x = "Physical Pos (Mb)",color = "Gen Map Used") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=14), legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14), plot.title = element_text(size=20), legend.position = c(0.2,0.8)) + annotate("rect", xmin = 92, xmax = 102, ymin = 0, ymax = 1200,
                                                                                                                          alpha = .2, fill = "blue")  + ylim(0, 1200) +
  annotate('text', x = 98, y = 980, 
           label = 'Centromere', color='blue',
           size = 4) +
  scale_y_continuous(trans = squish_trans(700,1800,8))

ggplot(chr3_snp2, aes(x = `SNP Start`, y = pos, color = "wt")) + geom_line(size = 2) + theme_bw() +
  geom_line(data = chr3_snp2ddm1, aes(x = `SNP Start`, y = pos, color = "ddm1"), size = 2) +
  geom_line(data = chr3_snp2zmet2, aes(x = `SNP Start`, y = pos, color = "zmet2"), size = 2) +
  geom_line(data = chr3_snp2recq4, aes(x = `SNP Start`, y = pos, color = "recq4"), size = 2) +
  geom_line(data = chr3_snp2fancm, aes(x = `SNP Start`, y = pos, color = "fancm"), size = 2) +
  geom_line(data = chr3_snp2ideal1, aes(x = `SNP Start`, y = pos, color = "10X"), size = 2) +
  geom_line(data = chr3_snp2ideal2, aes(x = `SNP Start`, y = pos, color = "ddm1/zmet2"), size = 2) +
  scale_fill_manual(values=group.colors) + labs(y = "Accumulative Genetic Pos (cM)",x = "Physical Pos (Mb)",color = "Gen Map Used") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=14), legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14), plot.title = element_text(size=20), legend.position = c(0.2,0.8)) + annotate("rect", xmin = 87, xmax = 89.5, ymin = 0, ymax = 3650,
                                                                                                                          alpha = .2, fill = "blue") +
  annotate('text', x = 85, y = 3680, label = 'Centromere', size = 4, color = 'blue') + scale_y_continuous(trans = squish_trans(800,3700,8), breaks = c(300, 800, 2000, 3000))

ggplot(chr4_snp2, aes(x = `SNP Start`, y = pos, color = "wt")) + geom_line(size = 2) + theme_bw() +
  geom_line(data = chr4_snp2ddm1, aes(x = `SNP Start`, y = pos, color = "ddm1"), size = 2) +
  geom_line(data = chr4_snp2zmet2, aes(x = `SNP Start`, y = pos, color = "zmet2"), size = 2) +
  geom_line(data = chr4_snp2recq4, aes(x = `SNP Start`, y = pos, color = "recq4"), size = 2) +
  geom_line(data = chr4_snp2fancm, aes(x = `SNP Start`, y = pos, color = "fancm"), size = 2) +
  geom_line(data = chr4_snp2ideal1, aes(x = `SNP Start`, y = pos, color = "10X"), size = 2) +
  geom_line(data = chr4_snp2ideal2, aes(x = `SNP Start`, y = pos, color = "ddm1/zmet2"), size = 2) +
  scale_fill_manual(values=group.colors) + labs(y = "Accumulative Genetic Pos (cM)",x = "Physical Pos (Mb)",color = "Gen Map Used") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=14), legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14), plot.title = element_text(size=20), legend.position = c(0.2,0.8)) + annotate("rect", xmin = 72, xmax = 93, ymin = 0, ymax = 2850,
                                                                                                                          alpha = .2, fill = "blue") +
  annotate('text', x = 83, y = 2890, label = 'Centromere', size = 4, color = 'blue') + scale_y_continuous(trans = squish_trans(700,2900,8), breaks = c(300, 700, 1000, 2500))


ggplot(chr5_snp2, aes(x = `SNP Start`, y = pos, color = "wt")) + geom_line(size = 2) + theme_bw() +
  geom_line(data = chr5_snp2ddm1, aes(x = `SNP Start`, y = pos, color = "ddm1"), size = 2) +
  geom_line(data = chr5_snp2zmet2, aes(x = `SNP Start`, y = pos, color = "zmet2"), size = 2) +
  geom_line(data = chr5_snp2recq4, aes(x = `SNP Start`, y = pos, color = "recq4"), size = 2) +
  geom_line(data = chr5_snp2fancm, aes(x = `SNP Start`, y = pos, color = "fancm"), size = 2) +
  geom_line(data = chr5_snp2ideal1, aes(x = `SNP Start`, y = pos, color = "10X"), size = 2) +
  geom_line(data = chr5_snp2ideal2, aes(x = `SNP Start`, y = pos, color = "ddm1/zmet2"), size = 2) +
  scale_fill_manual(values=group.colors) + labs(y = "Accumulative Genetic Pos (cM)",x = "Physical Pos (Mb)",color = "Genetic Map") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=14), legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14), plot.title = element_text(size=20), legend.position = c(0.2,0.8)) + annotate("rect", xmin = 94, xmax = 118, ymin = 0, ymax = 3150,
                                                                                                                          alpha = .2, fill = "blue") +
  annotate('text', x = 106, y = 3198, label = 'Centromere', size = 4, color = 'blue') + scale_y_continuous(trans = squish_trans(800,3200,8), breaks = c(300, 700, 1000, 2500))

ggplot(chr6_snp2, aes(x = `SNP Start`, y = pos, color = "wt")) + geom_line(size = 2) + theme_bw() +
  geom_line(data = chr6_snp2ddm1, aes(x = `SNP Start`, y = pos, color = "ddm1"), size = 2) +
  geom_line(data = chr6_snp2zmet2, aes(x = `SNP Start`, y = pos, color = "zmet2"), size = 2) +
  geom_line(data = chr6_snp2recq4, aes(x = `SNP Start`, y = pos, color = "recq4"), size = 2) +
  geom_line(data = chr6_snp2fancm, aes(x = `SNP Start`, y = pos, color = "fancm"), size = 2) +
  geom_line(data = chr6_snp2ideal1, aes(x = `SNP Start`, y = pos, color = "10X"), size = 2) +
  geom_line(data = chr6_snp2ideal2, aes(x = `SNP Start`, y = pos, color = "ddm1/zmet2"), size = 2) +
  scale_fill_manual(values=group.colors) + labs(y = "Accumulative Genetic Pos (cM)",x = "Physical Pos (Mb)",color = "Genetic Map") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=14), legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14), plot.title = element_text(size=20), legend.position = c(0.6,0.55)) + annotate("rect", xmin = 31, xmax = 33, ymin = 0, ymax = 3020,
                                                                                                                           alpha = .2, fill = "blue") +
  annotate('text', x = 32, y = 3050, label = 'Centromere', size = 4, color = 'blue') + scale_y_continuous(trans = squish_trans(850,3050,8), breaks = c(300, 700, 1000, 2500))

ggplot(chr7_snp2, aes(x = `SNP Start`, y = pos, color = "wt")) + geom_line(size = 2) + theme_bw() +
  geom_line(data = chr7_snp2ddm1, aes(x = `SNP Start`, y = pos, color = "ddm1"), size = 2) +
  geom_line(data = chr7_snp2zmet2, aes(x = `SNP Start`, y = pos, color = "zmet2"), size = 2) +
  geom_line(data = chr7_snp2recq4, aes(x = `SNP Start`, y = pos, color = "recq4"), size = 2) +
  geom_line(data = chr7_snp2fancm, aes(x = `SNP Start`, y = pos, color = "fancm"), size = 2) +
  geom_line(data = chr7_snp2ideal1, aes(x = `SNP Start`, y = pos, color = "10X"), size = 2) +
  geom_line(data = chr7_snp2ideal2, aes(x = `SNP Start`, y = pos, color = "ddm1/zmet2"), size = 2) +
  scale_fill_manual(values=group.colors) + labs(y = "Accumulative Genetic Pos (cM)",x = "Physical Pos (Mb)",color = "Genetic Map") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=14), legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14), plot.title = element_text(size=20), legend.position = c(0.15,0.7)) + annotate("rect", xmin = 45, xmax = 54, ymin = 0, ymax = 6550,
                                                                                                                           alpha = .2, fill = "blue") +
  annotate('text', x = 50, y = 6590, label = 'Centromere', size = 4, color = 'blue') + scale_y_continuous(trans = squish_trans(1800,6600,8), breaks = c(300, 1500, 4000, 6000))

ggplot(chr8_snp2, aes(x = `SNP Start`, y = pos, color = "wt")) + geom_line(size = 2) + theme_bw() +
  geom_line(data = chr8_snp2ddm1, aes(x = `SNP Start`, y = pos, color = "ddm1"), size = 2) +
  geom_line(data = chr8_snp2zmet2, aes(x = `SNP Start`, y = pos, color = "zmet2"), size = 2) +
  geom_line(data = chr8_snp2recq4, aes(x = `SNP Start`, y = pos, color = "recq4"), size = 2) +
  geom_line(data = chr8_snp2fancm, aes(x = `SNP Start`, y = pos, color = "fancm"), size = 2) +
  geom_line(data = chr8_snp2ideal1, aes(x = `SNP Start`, y = pos, color = "10X"), size = 2) +
  geom_line(data = chr8_snp2ideal2, aes(x = `SNP Start`, y = pos, color = "ddm1/zmet2"), size = 2) +
  scale_fill_manual(values=group.colors) + labs(y = "Accumulative Genetic Pos (cM)",x = "Physical Pos (Mb)",color = "Genetic Map") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=14), legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14), plot.title = element_text(size=20), legend.position = c(0.2,0.8)) + annotate("rect", xmin = 56, xmax = 80, ymin = 0, ymax = 4950,
                                                                                                                          alpha = .2, fill = "blue") +
  annotate('text', x = 68, y = 4998, label = 'Centromere', size = 4, color = 'blue') + scale_y_continuous(trans = squish_trans(1400,5000,8), breaks = c(500, 1400, 2000, 4000))

ggplot(chr9_snp2, aes(x = `SNP Start`, y = pos, color = "wt")) + geom_line(size = 2) + theme_bw() +
  geom_line(data = chr9_snp2ddm1, aes(x = `SNP Start`, y = pos, color = "ddm1"), size = 2) +
  geom_line(data = chr9_snp2zmet2, aes(x = `SNP Start`, y = pos, color = "zmet2"), size = 2) +
  geom_line(data = chr9_snp2recq4, aes(x = `SNP Start`, y = pos, color = "recq4"), size = 2) +
  geom_line(data = chr9_snp2fancm, aes(x = `SNP Start`, y = pos, color = "fancm"), size = 2) +
  geom_line(data = chr9_snp2ideal1, aes(x = `SNP Start`, y = pos, color = "10X"), size = 2) +
  geom_line(data = chr9_snp2ideal2, aes(x = `SNP Start`, y = pos, color = "ddm1/zmet2"), size = 2) +
  scale_fill_manual(values=group.colors) + labs(y = "Accumulative Genetic Pos (cM)",x = "Physical Pos (Mb)",color = "Genetic Map") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=14), legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14), plot.title = element_text(size=20), legend.position = c(0.4,0.8)) + annotate("rect", xmin = 34, xmax = 34.2, ymin = 0, ymax = 4850,
                                                                                                                          alpha = .2, fill = "blue") +
  annotate('text', x = 34, y = 4890, label = 'Centromere', size = 4, color = 'blue') + scale_y_continuous(trans = squish_trans(1450,4900,8), breaks = c(500, 1400, 2000, 4000))

ggplot(chr10_snp2, aes(x = `SNP Start`, y = pos, color = "wt")) + geom_line(size = 2) + theme_bw() +
  geom_line(data = chr10_snp2ddm1, aes(x = `SNP Start`, y = pos, color = "ddm1"), size = 2) +
  geom_line(data = chr10_snp2zmet2, aes(x = `SNP Start`, y = pos, color = "zmet2"), size = 2) +
  geom_line(data = chr10_snp2recq4, aes(x = `SNP Start`, y = pos, color = "recq4"), size = 2) +
  geom_line(data = chr10_snp2fancm, aes(x = `SNP Start`, y = pos, color = "fancm"), size = 2) +
  geom_line(data = chr10_snp2ideal1, aes(x = `SNP Start`, y = pos, color = "10X"), size = 2) +
  geom_line(data = chr10_snp2ideal2, aes(x = `SNP Start`, y = pos, color = "ddm1/zmet2"), size = 2) +
  scale_fill_manual(values=group.colors) + labs(y = "Accumulative Genetic Pos (cM)",x = "Physical Pos (Mb)",color = "Genetic Map") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=14), legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14), plot.title = element_text(size=20), legend.position = c(0.15,0.8)) + annotate("rect", xmin = 37, xmax = 46, ymin = 0, ymax = 6150,
                                                                                                                           alpha = .2, fill = "blue") +
  annotate('text', x = 42, y = 6190, label = 'Centromere', size = 4, color = 'blue') + scale_y_continuous(trans = squish_trans(1900,6200,8), breaks = c(500, 1800, 3000, 6000))

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