library(tidyverse)
library(plyr)
library(reshape2)
library(gridExtra)

### Set working directories
setwd("/Users/derekwong/Google Drive/Post-Doc/Yip_Exomes/Figures/Figure 2A")
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/Yip_Exomes/somatic"
samples <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/Yip_Exomes/yip_exome_sample_list.txt"

### Import data
data_muts <- read.delim(file.path(path, "Yip_exomes_snps.txt"))
samples <- read.delim(samples)

### Only keep recurrances
keep <- c("T1", "T9", "T5", "T10", "T12", "T13", "T15", "T16")
samples <- samples[samples$STY_code %in% keep, ]

### Seperate out Mutations
keep <- c("Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "vaf")
data_muts <- data_muts[ , colnames(data_muts) %in% keep]
data_muts$Tumor_Sample_Barcode <- gsub("_.*", "", data_muts$Tumor_Sample_Barcode)

T1 <- data_muts[data_muts$Tumor_Sample_Barcode == "T1", ]
T9 <- data_muts[data_muts$Tumor_Sample_Barcode == "T9", ]

T5 <- data_muts[data_muts$Tumor_Sample_Barcode == "T5", ]
T10 <- data_muts[data_muts$Tumor_Sample_Barcode == "T10", ]

T12 <- data_muts[data_muts$Tumor_Sample_Barcode == "T12", ]
T13 <- data_muts[data_muts$Tumor_Sample_Barcode == "T13", ]

T15 <- data_muts[data_muts$Tumor_Sample_Barcode == "T15", ]
T16 <- data_muts[data_muts$Tumor_Sample_Barcode == "T16", ]

### Merge pairs
T1_9 <- merge(T1, T9, by = c("Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"), all = TRUE)
T1_9[is.na(T1_9)] <- 0

T5_10 <- merge(T5, T10, by = c("Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"), all = TRUE)
T5_10[is.na(T5_10)] <- 0

T12_13 <- merge(T12, T13, by = c("Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"), all = TRUE)
T12_13[is.na(T12_13)] <- 0

T15_16 <- merge(T15, T16, by = c("Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"), all = TRUE)
T15_16[is.na(T15_16)] <- 0

### Set colors
T1_9$color <- ifelse(T1_9$vaf.x > 0 & T1_9$vaf.y > 0, "yes", "no")
T5_10$color <- ifelse(T5_10$vaf.x > 0 & T5_10$vaf.y > 0, "yes", "no")
T12_13$color <- ifelse(T12_13$vaf.x > 0 & T12_13$vaf.y > 0, "yes", "no")
T15_16$color <- ifelse(T15_16$vaf.x > 0 & T15_16$vaf.y > 0, "yes", "no")

### Make annotations
T1 <- nrow(T1)
T9 <- nrow(T9)
T19 <- abs(nrow(T1_9) - T1 - T9)

T5 <- nrow(T5)
T10 <- nrow(T10)
T510 <- abs(nrow(T5_10) - T5 - T10)

T12 <- nrow(T12)
T13 <- nrow(T13)
T1213 <- abs(nrow(T12_13) - T12 - T13)

T15 <- nrow(T15)
T16 <- nrow(T16)
T1516 <- abs(nrow(T15_16) - T15 - T16)

### Set Theme
theme <- theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
               axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.position="none",
               legend.title = element_text(size = 15), 
               legend.text = element_text(size = 12), 
               axis.text.x = element_text(size = 15),
               axis.text.y = element_text(size = 15),
               axis.title = element_text(size = 15),
               plot.margin = margin(t = 5,  # Top margin
                                    r = 20,  # Right margin
                                    b = 5,  # Bottom margin
                                    l = 5)) # Left mar)

### Plot Mutation Counts
T1_9_plot <- ggplot(T1_9, aes(x = vaf.y, y = vaf.x, fill = color)) + 
  geom_point(pch = 21, alpha = 0.5, size = 3) +
  scale_fill_manual(values = c("grey", "red")) +
  xlab("Primary (T9) VAF") + 
  ylab("Recurrance (T1) VAF") +
  ggtitle("T1-T9") + 
  theme + 
  scale_y_continuous(limits=c(-10, 100), expand = c(0,0)) + 
  scale_x_continuous(limits=c(-10, 100), expand = c(0,0)) + 
  annotate("text", x = 40, y = 85, label = paste0("Primary (T9) = ", T9, 
                                                  "\nRecurrance (T1) = ", T1, 
                                                  "\nOverlap = ", T19),
           size = 6)
T1_9_plot

T5_10_plot <- ggplot(T5_10, aes(x = vaf.y, y = vaf.x, fill = color)) + 
  geom_point(pch = 21, alpha = 0.5, size = 3) +
  scale_fill_manual(values = c("grey", "red")) +
  xlab("Primary (T10) VAF") + 
  ylab("Recurrance (T5) VAF") +
  ggtitle("T5-T10") + 
  theme + 
  scale_y_continuous(limits=c(-10, 100), expand = c(0,0)) + 
  scale_x_continuous(limits=c(-10, 100), expand = c(0,0)) + 
  annotate("text", x = 40, y = 85, label = paste0("Primary (T10) = ", T10, 
                                                  "\nRecurrance (T5) = ", T5, 
                                                  "\nOverlap = ", T510),
           size = 6)
T5_10_plot

T12_13_plot <- ggplot(T12_13, aes(x = vaf.y, y = vaf.x, fill = color)) + 
  geom_point(pch = 21, alpha = 0.5, size = 3) +
  scale_fill_manual(values = c("grey", "red")) +
  xlab("Primary (T13) VAF") + 
  ylab("Recurrance (T12) VAF") +
  ggtitle("T12-T13") + 
  theme + 
  scale_y_continuous(limits=c(-10, 100), expand = c(0,0)) + 
  scale_x_continuous(limits=c(-10, 100), expand = c(0,0)) + 
  annotate("text", x = 40, y = 85, label = paste0("Primary (T13) = ", T13, 
                                                  "\nRecurrance (T12) = ", T12, 
                                                  "\nOverlap = ", T1213),
           size = 6)
T12_13_plot

T15_16_plot <- ggplot(T15_16, aes(x = vaf.y, y = vaf.x, fill = color)) + 
  geom_point(pch = 21, alpha = 0.5, size = 3) +
  scale_fill_manual(values = c("grey", "red")) +
  xlab("Primary (T16) VAF") + 
  ylab("Recurrance (T15) VAF") +
  ggtitle("T15-T16") + 
  theme + 
  scale_y_continuous(limits=c(-10, 100), expand = c(0,0)) + 
  scale_x_continuous(limits=c(-10, 100), expand = c(0,0)) + 
  annotate("text", x = 40, y = 85, label = paste0("Primary (T16) = ", T16, 
                                                  "\nRecurrance (T15) = ", T15, 
                                                  "\nOverlap = ", T1516),
           size = 6)
T15_16_plot

### Ensemble the plots
plot <- ggarrange(T1_9_plot, T5_10_plot, T12_13_plot, T15_16_plot, ncol = 2, nrow = 2)

pdf("tumour_pairs.pdf", height = 9, width = 9)
annotate_figure(plot, top = text_grob("SNV Comparison (Tumour Pairs)", color = "black", face = "bold", size = 20))
dev.off()
