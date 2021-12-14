library(tidyverse)
library(ComplexHeatmap)
library(dplyr)
library(stringr)

### Set working directories
setwd("/Users/derekwong/Google Drive/Post-Doc/Yip_Exomes/Figures/Figure 1B")
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/Yip_Exomes/somatic"
samples <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/Yip_Exomes/yip_exome_sample_list.txt"

### Import data
data_onco <- read.delim(file.path(path, "Yip_oncoplot.txt"))
seq <- read.delim(file.path(path, "Yip_sequencing.txt"))
samples <- read.delim(samples)

### Order data
data_onco <- data_onco[order(data_onco$NF2_mutation,
                             data_onco$NF2_SV_RNA, decreasing = TRUE), ]
row.names(data_onco) <- data_onco$STY_code

seq <- seq[seq$STY_code %in% data_onco$STY_code, ]
row.names(seq) <- seq$STY_code
seq <- seq[data_onco$STY_code, ]

samples <- samples[samples$STY_code %in% data_onco$STY_code, ]
row.names(samples) <- samples$STY_code
samples <- samples[data_onco$STY_code, ]

## Set annotation variables
samples$Sex <- factor(samples$Sex, levels = c("male", "female"),
                      labels = c("Male", "Female"))
data_sex <- as.matrix(samples$Sex)
row.names(data_sex) <- data_onco$STY_code

samples$Type <- factor(samples$Type, levels = c("Primary", "Recurrence"))
data_type <- as.matrix(samples$Type)
row.names(data_type) <- data_onco$STY_code

samples$Subtype <- factor(samples$Subtype, levels = c("transitional", "fibrous", "chordoid", "meningothelial"),
                          labels = c("Transitional", "Fibrous", "Chordoid", "Meningothelial"))
data_subtype <- as.matrix(samples$Subtype)
row.names(data_subtype) <- data_onco$STY_code

samples$Grade <- factor(samples$Grade, levels = c("I", "II", "III"))
data_grade <- as.matrix(samples$Grade)
row.names(data_grade) <- data_onco$STY_code

samples$Recurrence <- factor(samples$Recurrence, levels = c("no", "yes"),
                             labels = c("Did Not Recur", "Recurred"))
data_recur <- as.matrix(samples$Recurrence)
row.names(data_recur) <- data_onco$STY_code

samples$Initial.cancer.diagnosis <- factor(samples$Initial.cancer.diagnosis, levels = c("ALL", "infant ALL", "lymphoblastic lynphoma", "neuroblastoma"),
                                           labels = c("ALL", "Infant ALL", "LBL", "Neuroblastoma"))
data_initial <- as.matrix(samples$Initial.cancer.diagnosis)
row.names(data_initial) <- data_onco$STY_code

seq$exome <- factor(seq$exome, levels = c("pass", "fail"),
                        labels = c("Analyzed", "No Germline"))
data_exome <- as.matrix(seq$exome)
row.names(data_exome) <- data_onco$STY_code

seq$rna <- factor(seq$rna, levels = c("pass", "fail"),
                    labels = c("Analyzed", "Failed"))
data_rna <- as.matrix(seq$rna)
row.names(data_rna) <- data_onco$STY_code

## Format heatmap
data_onco[is.na(data_onco)] <- ""
data_onco <- data_onco[, c(2:10)]
data_onco <- as.data.frame(t(data_onco))
row.names(data_onco) <- c("Chr 1p", "Chr 22q", "NF2 (DNA)", "NF2 (RNA)", "AKT1", "SMO", "KLF4", "PIK3CA", "TRAF7")
data_onco <- as.matrix(data_onco)

## Set colours
col <- c(missense = "#33a02c", frameshift = "black", nonsense = "#E31A1C",
         splice_site = "#E31A1C", loss = "#1f78b4", fusion = "#6A3D9A",
         structural_variant = "#FF7F00", inframe_deletion = "#FFFF99")
col_sex <- c(Male = "#a6cee3", Female = "#fb9a99")
col_type <- c(Primary = "#FEE0D2", Recurrence = "#fb9a99")
col_grade <- c(I = "#FED976", II = "#FEB24C", III = "#FD8D3C")
col_subtype <- c(Chordoid = "#A1D99B", Fibrous = "#74C476", Meningothelial = "#41AB5D", Transitional = "#238B45")
col_recur <- c("Did Not Recur" = "grey", "Recurred" = "#fb9a99")
col_initial <- c(ALL = "#A6CEE3", "Infant ALL" = "#B2DF8A", LBL = "#CAB2D6", Neuroblastoma = "#FDBF6F")
col_exome <- c("Analyzed" = "grey", "No Germline" = "#fb9a99")
col_rna <- c("Analyzed" = "grey", "Failed" = "#fb9a99")

## Set variables
alter_fun = function(x, y, w, h, v) {
  # background
  grid.rect(x, y, w*0.75, h*0.9, gp = gpar(fill = "grey", col = NA))
  # alterations
  n = sum(v)  # how many alterations for current gene in current sample
  h = h
  w = w
  # use `names(which(v))` to correctly map between `v` and `col`
  if(n) grid.rect(x, y - h*0.5 + 1:n/n*h, w, 1/n*h, 
                  gp = gpar(fill = col[names(which(v))], col = "white"), just = "top")
}

# Set annotation layers
top_annotation <- HeatmapAnnotation(Sex = data_sex,
                                    "Initial Diagnosis" = data_initial,
                                    Subtype = data_subtype,
                                    "WHO Grade" = data_grade,
                                    Recurrence = data_recur,
                                    Timepoint = data_type,
                                    "Exome Sequencing" = data_exome,
                                    "RNA Sequencing" = data_rna,
                                    col = list(Sex = col_sex,
                                               "Initial Diagnosis" = col_initial,
                                               Subtype = col_subtype,
                                               "WHO Grade" = col_grade,
                                               Recurrence = col_recur,
                                               Timepoint = col_type,
                                               "Exome Sequencing" = col_exome,
                                               "RNA Sequencing" = col_rna),
                                    annotation_legend_param = list(title_gp = gpar(fontsize = 15),
                                                                   labels_gp = gpar(fontsize = 13),
                                                                   Subtype = list(labels = c("Transitional", "Fibrous", "Chordoid", "Meningothelial", "NA"))),
                                    annotation_name_gp = gpar(fontsize = 15),
                                    annotation_name_side = "left",
                                    border = TRUE,
                                    simple_anno_size = unit(0.6, "cm"))

right_annotation <- rowAnnotation(row_barplot = anno_oncoprint_barplot(border = TRUE, height = unit(4, "cm"), 
                                                                       axis_param = list(side = "bottom", 
                                                                                         labels_rot = 0,
                                                                                         gp = gpar(fontsize = 10))))
  

## Set labels
heatmap_legend_param = list(title = "Alteration", 
                            at = c("loss", "structural_variant", "nonsense",  "frameshift", "fusion", "missense", "splice_site", "inframe_deletion"), 
                            labels = c(loss = "Deletion", structural_variant = "Structural Variant", nonsense = "Nonsense", frameshift =  "Frameshift", 
                                       fusion = "Fusion", missense = "Missense", splice_site = "Splice Site", inframe_deletion = "Inframe Deletion"),
                            background = "white",
                            border = FALSE,
                            title_gp = gpar(fontsize = 15),
                            labels_gp = gpar(fontsize = 13))

## Set orders
col_order <- colnames(data_onco)
row_order <- row.names(data_onco)

## Generate oncoprint
oncoPrint <- oncoPrint(data_onco,
                       alter_fun = alter_fun, 
                       col = col,
                       row_order = row_order,
                       column_order = col_order,
                       row_names_side = "left", pct_side = "right",
                       top_annotation = top_annotation,
                       right_annotation = right_annotation,
                       heatmap_legend_param = heatmap_legend_param,
                       row_names_gp = gpar(fontsize = 15),
                       border = TRUE,
                       border_gp = gpar(col = "black"))
oncoPrint

pdf("Oncoprint.pdf", width = 9, height = 4.5)
draw(oncoPrint, heatmap_legend_side = "right", merge_legend = TRUE)
dev.off()
