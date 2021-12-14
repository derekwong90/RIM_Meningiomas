library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(matrixStats)

### Set working variables ###
setwd("/Users/derekwong/Google Drive/Post-Doc/Yip_Exomes/Figures/Figure 1D")
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/Yip_Exomes/copy_number"
samples <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/Yip_Exomes/yip_exome_sample_list.txt"
project <- "Yip_exomes"

#sequenza <- list.files(path, "sequenza", full.names = TRUE)
gatk <- list.files(path, "gatk", full.names = TRUE)

### Read in CNAs (Sequenza)
#sequenza <- read.delim(sequenza)
gatk <- read.delim(gatk)
samples <- read.delim(samples)

### Order samples
samples <- samples[order(factor(samples$Type, levels = c("Primary", "Recurrance")),
                         factor(samples$Grade, levels = c("I", "II", "III"))), ]

### Sort and format Sequenza
#sequenza$Chromosome <- factor(sequenza$Chromosome, levels = paste0('chr',c(1:22,'X','Y')))
#sequenza <- sequenza[order(sequenza$Chromosome, sequenza$Start),]
#sequenza <- sequenza[which(sequenza$Chromosome != 'chrY'),]
#rownames(sequenza) <- paste0(sequenza$Chromosome, "_", sequenza$Start, "_", sequenza$End)
#sequenza <- sequenza[ , !(colnames(sequenza) %in% c("Chromosome", "Start", "End", "GeneID", "Symbol", "Target"))]
#colnames(sequenza) <- gsub("_.*", "", colnames(sequenza))
#sequenza <- sequenza[ , colnames(sequenza) %in% samples$STY_code]
#samples2 <- samples[samples$STY_code %in% colnames(sequenza), ]
#sequenza <- sequenza[ , samples2$STY_code]
#sequenza <- as.matrix(t(sequenza))

### Sort and format GATK
gatk$Chromosome <- factor(gatk$Chromosome, levels = paste0('chr',c(1:22,'X','Y')))
gatk <- gatk[order(gatk$Chromosome, gatk$Start),]
gatk <- gatk[which(gatk$Chromosome != 'chrY'),]
rownames(gatk) <- paste0(gatk$Chromosome, "_", gatk$Start, "_", gatk$End)
gatk <- gatk[ , !(colnames(gatk) %in% c("Chromosome", "Start", "End", "GeneID", "Symbol", "Target"))]
colnames(gatk) <- gsub("_.*", "", colnames(gatk))
gatk <- gatk[ , colnames(gatk) %in% samples$STY_code]
gatk <- gatk[ , samples$STY_code]
gatk <- as.matrix(gatk)

### Set row split names
row_split <- as.matrix(rownames(gatk))
row_split <- gsub("_.*", "", row_split)
row_split <- gsub("chr", "", row_split)
row_split <- factor(row_split, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"),
                    labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21     ", "22", "X"))
dim(row_split) <- c(length(row_split), 1)
row.names(row_split) <- row.names(gatk)

### Calculate Frequency
gatk_gain <- gatk
gatk_gain[gatk_gain >= 0.25] <- 1
gatk_gain[gatk_gain < 0.25] <- 0
data_gain <- rowSums(gatk_gain/ncol(gatk_gain))

gatk_loss <- gatk
gatk_loss[gatk_loss > -0.25] <- 0
gatk_loss[gatk_loss <= -0.25] <- -1
data_loss <- rowSums(gatk_loss/ncol(gatk_loss))

### Calculate % Genome altered
gain <- colSums(gatk_gain)
loss <- colSums(abs(gatk_loss))
data_genome <- as.matrix((gain + loss)/nrow(gatk_loss))

### Set clinical information
samples$Sex <- factor(samples$Sex, levels = c("male", "female"),
                      labels = c("Male", "Female"))
data_sex <- as.matrix(samples$Sex)
row.names(data_sex) <- samples$STY_code

samples$Type <- factor(samples$Type, levels = c("Primary", "Recurrance"))
data_type <- as.matrix(samples$Type)
row.names(data_type) <- samples$STY_code

samples$Grade <- factor(samples$Grade, levels = c("I", "II", "III"))
data_grade <- as.matrix(samples$Grade)
row.names(data_grade) <- samples$STY_code

samples$Subtype <- factor(samples$Subtype, levels = c("transitional", "fibrous", "chordoid", "meningothelial"),
                          labels = c("Transitional", "Fibrous", "Chordoid", "Meningothelial"))
data_subtype <- as.matrix(samples$Subtype)
row.names(data_subtype) <- samples$STY_code

### Set colors
col_fun <- colorRamp2(c(-1, -0.1, 0, 0.1, 1), 
                      c("#1f78b4", "white", "white", "white", "#e31a1c"))
col_sex <- c(Male = "#a6cee3", Female = "#fb9a99")
col_type <- c(Primary = "#FEE0D2", Recurrance = "#fb9a99")
col_grade <- c(I = "#FED976", II = "#FEB24C", III = "#FD8D3C")
col_subtype <- c(Chordoid = "#A1D99B", Fibrous = "#74C476", Meningothelial = "#41AB5D", Transitional = "#238B45")

### Set Annotations
top_annotation <- HeatmapAnnotation(Sex = data_sex,
                                    Subtype = data_subtype,
                                    "WHO Grade" = data_grade,
                                    Sample = data_type,
                                    "% Genome Altered" = anno_barplot(data_genome,
                                                                      axis_param = list(at = c(0, 0.3),
                                                                                        gp = gpar(fontsize = 15)),
                                                                      height = unit(1.25, "cm")),
                                    col = list(Sex = col_sex,
                                               Sample = col_type,
                                               "WHO Grade" = col_grade,
                                               Subtype = col_subtype),
                                    annotation_legend_param = list(title_gp = gpar(fontsize = 15),
                                                                   labels_gp = gpar(fontsize = 13),
                                                                   Subtype = list(labels = c("Transitional", "Fibrous", "Chordoid", "Meningothelial", "NA"))),
                                    annotation_name_gp = gpar(fontsize = 15),
                                    border = TRUE,
                                    simple_anno_size = unit(0.45, "cm"))

right_annotation <- rowAnnotation("CNV Frequency" = anno_barplot(data_loss,
                                                                 bar_width = 1,
                                                                 gp = gpar(fill = "#1f78b4",
                                                                           col = "#1f78b4"),
                                                                 axis_param = list(at = c(-1, 0),
                                                                                   labels = c(1, 0),
                                                                                   labels_rot = 0,
                                                                                   gp = gpar(fontsize = 13)),
                                                                 ylim = c(-1, 0)),
                                  " " = anno_barplot(data_gain,
                                                     bar_width = 1,
                                                     gp = gpar(fill = "#e31a1c",
                                                               col = "#e31a1c"),
                                                     axis_param = list(at = c(1),
                                                                       labels_rot = 0,
                                                                       gp = gpar(fontsize = 13)),
                                                     ylim = c(0, 1)),
                                  "  " = anno_mark(at = c(18771),
                                                   labels = c("NF2"),
                                                   which = "row",
                                                   side = "right",
                                                   labels_rot = 0,
                                                   labels_gp = gpar(fontsize = 13),
                                                   padding = unit(5, "mm")),
                                  annotation_name_gp = gpar(fontsize = 15),
                                  width = unit(6, "cm"))

### Set sample orders and splits
col_order <- samples$STY_code
row_order <- rownames(gatk)

### Set Legend Parameters
heatmap_legend_param = list(title = "CNV", 
                            border = FALSE,
                            at = c(-1, 0, 1),
                            title_gp = gpar(fontsize = 15),
                            labels_gp = gpar(fontsize = 13))

### Plot CNV heatmap
cnv <- Heatmap(gatk,
        col = col_fun,
        top_annotation = top_annotation,
        right_annotation = right_annotation,
        heatmap_legend_param = heatmap_legend_param,
        row_order = row_order,
        column_order = col_order,
        show_column_dend = FALSE,
        show_column_names = FALSE,
        show_row_dend = FALSE,
        show_row_names = FALSE,
        row_split = row_split,
        row_title_rot = 0,
        row_title_gp = gpar(fontsize = 15),
        border = TRUE)
cnv

pdf("copy_number.pdf", height = 9, width = 7)
draw(cnv, heatmap_legend_side = "right", merge_legend = TRUE)
dev.off()

