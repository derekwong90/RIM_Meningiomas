library(tidyverse)
library(plyr)
library(reshape2)

### Set working directories
setwd("/Users/derekwong/Google Drive/Post-Doc/Yip_Exomes/Figures/mutation_type")
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/Yip_Exomes/somatic"
samples <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/Yip_Exomes/yip_exome_sample_list.txt"

### Import data
data_muts <- read.delim(file.path(path, "Yip_exomes_snps.txt"))
samples <- read.delim(samples)

### Format Mutations
keep <- c("Variant_Classification", "Tumor_Sample_Barcode")
data_muts <- data_muts[ , keep]
data_muts$Tumor_Sample_Barcode <- gsub("_.*", "", data_muts$Tumor_Sample_Barcode)

variants <- ddply(data_muts, .(data_muts$Variant_Classification, data_muts$Tumor_Sample_Barcode), nrow)
names(variants) <- c("type", "sample", "freq")

totals <- dcast(variants, sample ~ type)
totals[is.na(totals)] <- 0
row.names(totals) <- totals$sample
totals <- totals[ , -1]
totals$other <- rowSums(totals[ , c("3'Flank", "5'Flank", "IGR", "Nonstop_Mutation", "RNA", "Splice_Region")], na.rm=TRUE)
totals$UTR <- rowSums(totals[ , c("3'UTR", "5'UTR")], na.rm=TRUE)
totals <- totals[ , !(colnames(totals) %in% c("3'Flank", "5'Flank", "IGR", "Nonstop_Mutation", "RNA", "Splice_Region", "3'UTR", "5'UTR"))]
totals$total <- rowSums(totals)

### Melt data
variants <- totals[ , !(colnames(totals) == c("total"))]
variants$sample <- rownames(variants)
variants <- melt(variants)
variants$variable <- factor(variants$variable, levels = c("Intron", "Missense_Mutation", "Silent", "other", "UTR", "Nonsense_Mutation", "Splice_Site"),
                            labels = c("Intron", "Missense", "Synonymous", "Other", "Untranslated Regions", "Nonsense", "Splice Site"))

### Calculate TMB
mutations <- as.matrix(rowSums(totals[ , c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site")], na.rm=TRUE))
mutations <- mutations/39 # Total size of the exome target regions is 39Mb
mutations <- as.data.frame(mutations)
mutations$sample <- row.names(mutations)

### Order by decreasing
totals <- totals[order(totals$total, decreasing = TRUE), ]
variants$sample <- factor(variants$sample, levels = c(row.names(totals)))
mutations$sample <- factor(mutations$sample, levels = c(row.names(totals)))

### Plot Mutation Counts
var_plot <- ggplot() + 
  geom_bar(data = variants, mapping = aes(x = sample, y = value, fill = variable),
           position = position_stack(reverse = TRUE), stat = "identity") +
  geom_line(data = mutations, mapping = aes(x = sample, y = V1*100, group = 1), size = 1, color = "black") +
  xlab("Sample") + 
  labs(fill = "Variant Type") + 
  ggtitle("SNVs") + 
  scale_fill_manual(values = c("Intron" = "#A6CEE3", "Missense" = "#B2DF8A", "Synonymous" = "grey", "Other" = "#FDBF6F", 
                               "Untranslated Regions" = "#FFFF99", "Nonsense" = "#FB9A99", "Splice Site" = "#CAB2D6")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        #legend.position=c(0.75, 0.75),
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 12), 
        axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15)) + 
  scale_y_continuous(name = "Count", limits=c(-40, 400), expand = c(0,0),
                     sec.axis = sec_axis(trans = ~./100, name = "TMB")) + 
  annotate(geom = "curve", x = 12, y = -1, xend = 8, yend = -1, curvature = -0.4, arrow = arrow(length = unit(2, "mm"))) +
  annotate(geom = "curve", x = 13, y = -1, xend = 9, yend = -1, curvature = -0.4, arrow = arrow(length = unit(2, "mm"))) +
  annotate(geom = "curve", x = 3, y = -1, xend = 5, yend = -1, curvature = 0.6, arrow = arrow(length = unit(2, "mm"))) +
  annotate(geom = "curve", x = 6, y = -1, xend = 4, yend = -1, curvature = -0.6, arrow = arrow(length = unit(2, "mm")))

var_plot
ggsave("mutation_count.pdf", var_plot, device = "pdf", width = 7, height = 4.5, units = "in")

write.table(mutations, "TMB.txt", sep = "\t", row.names = FALSE)
write.table(totals, "Mutation_counts.txt", sep = "\t", row.names = TRUE)
