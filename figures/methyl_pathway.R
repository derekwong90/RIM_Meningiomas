library(tidyverse)
library(dplyr)
library(limma)
library(missMethyl)

### Set working directories
setwd("/Users/derekwong/Google Drive/Post-Doc/Yip_Exomes/Figures/methyl_pathway")
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/Yip_Exomes/methylation"

### Read in files
dat <- read_tsv(file.path(path, "1-4_betas-norm-filt.tsv"))
dat <- as.data.frame(dat)
master <- read.csv(file.path(path, "1-4_pdat-filt.csv"))
master <- master[complete.cases(master), ]

### Order data
row.names(dat) <- dat$probe
dat <- dat[ , colnames(dat) %in% master$x850k_file]
dat <- dat[ , master$x850k_file]
dat <- as.matrix(dat)

### Set variables
master$group1 <- ifelse(master$group == 1, 1, 0)
master$group2 <- ifelse(master$group == 2, 1, 0)

group1 <- factor(master$group1)
group2 <- factor(master$group2)
design <- model.matrix(~0 + group1 + group2, data = master)
colnames(design) <- c("group1", "group2")
fit <- lmFit(dat, design)
contMatrix <- makeContrasts(group1-group2,
                            levels = design)

### Find Differentially Methylated probes
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
summary <- summary(decideTests(fit2))
write.table(summary, file.path(path, "differentially_methylated_summary.txt"), sep = "\t", row.names = FALSE)

### Extract list of DMPs
DMPs <- topTable(fit2, num=Inf, coef=1)
sigDMPs <- DMPs[DMPs$adj.P.Val <= 0.05, ]
sigDMPs_up <- sigDMPs[sigDMPs$logFC > 0, ]
sigDMPs_down <- sigDMPs[sigDMPs$logFC < -0, ]

### Do functional pathway analysis
all <- row.names(DMPs)
sigCpGs <- row.names(sigDMPs)
sigCpGs_up <- row.names(sigDMPs_up)
sigCpGs_down <- row.names(sigDMPs_down)

par(mfrow=c(1,1))
gst <- gometh(sig.cpg=sigCpGs, all.cpg=all, plot.bias=TRUE)
gst_up <- gometh(sig.cpg=sigCpGs_up, all.cpg=all, plot.bias=TRUE)
gst_down <- gometh(sig.cpg=sigCpGs_down, all.cpg=all, plot.bias=TRUE)

### Subset significant pathways
gst_up <- gst_up[gst_up$FDR < 0.05, ]
gst_up <- gst_up[order(gst_up$FDR), ]
gst_down <- gst_down[gst_down$FDR < 0.05, ]
gst_down <- gst_down[order(gst_down$FDR), ]

### Write tables
DMPs$probe <- row.names(DMPs)
DMPs <- DMPs %>% relocate(probe, .before = logFC)
write_tsv(DMPs, file.path(path, "differentially_methylated_probes.tsv"))

write.table(gst_up, file.path(path, "methylated_pathways_up.txt"), sep = "\t")
write.table(gst_down, file.path(path, "methylated_pathways_down.txt"), sep = "\t")

### Plot dysregulated pathways
gst_down <- gst_down[1:5, ]
gst_down$dir <- "down"
gst_up$dir <- "up"
gst_dys <- rbind(gst_up, gst_down)
gst_dys$logFDR <- -log10(gst_dys$FDR)
gst_dys$dir <- factor(gst_dys$dir, levels = c("up", "down"),
                      labels = c("hyper", "hypo"))
gst_dys <- gst_dys[order(gst_dys$logFDR), ]
gst_dys$TERM <- factor(gst_dys$TERM, levels = c(gst_dys$TERM))

gst_plot <- ggplot(gst_dys, aes(x = TERM, y = logFDR, fill = dir)) + 
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(dir),
             scales = "free",
             space = "free") +
  xlab("") + 
  ylab("FDR (-log10)") +
  ggtitle("Cluster 1 vs Cluster 2") + 
  scale_fill_manual(values = c("#FB9A99", "#A6CEE3")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 15),
        legend.position = "none",
        axis.line.x = element_line(),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        axis.title = element_text(size = 15)) +
  coord_flip()
gst_plot
ggsave("gene_ontology.pdf", gst_plot, device = "pdf", width = 8, height = 3, units = "in")
