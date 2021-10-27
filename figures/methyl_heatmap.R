library(tidyverse)
library(dplyr)
library(tidyr)
library(circlize)
library(ComplexHeatmap)

### Set working directories
setwd("/Users/derekwong/Google Drive/Post-Doc/Yip_Exomes/Figures/methyl_heatmap")
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/Yip_Exomes/methylation"
samples <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/Yip_Exomes/yip_exome_sample_list.txt"
NF2 <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/Yip_Exomes/somatic/Yip_oncoplot.txt"

### Import data
DKFZ <- read.csv(file.path(path, "DKFZ_outputs.csv"))
NF2 <- read.delim(NF2)
dat <- read_tsv(file.path(path, "1-4_betas-norm-filt.tsv"))
dat <- as.data.frame(dat)
master <- read.csv(file.path(path, "1-4_pdat-filt.csv"))
samples <- read.delim(samples)

### Remove failed samples
master <- master[master$classifier == "meningioma", ]
master <- master[master$sty_internal_tumour_code %in% samples$STY_code, ]

row.names(dat) <- dat$probes
dat <- dat[ , colnames(dat) %in% master$x850k_file]

samples <- samples[samples$STY_code %in% master$sty_internal_tumour_code, ]
NF2 <- NF2[NF2$STY_code %in% master$sty_internal_tumour_code, ]
DKFZ <- DKFZ[DKFZ$sample %in% master$sty_internal_tumour_code, ]

### Order samples
samples <- samples[order(samples$STY_code), ]
master <- master[order(master$sty_internal_tumour_code), ]
NF2 <- NF2[order(NF2$STY_code), ]
dat <- dat[ , master$x850k_file]
DKFZ <- DKFZ[order(DKFZ$sample), ]

### Format NF2 mutation data
NF2$mutation <- ifelse(NF2$NF2_mutation == "frameshift", "frameshift",
                       ifelse(NF2$NF2_mutation == "structural_variant", "rearrangement",
                              ifelse(NF2$NF2_SV_RNA == "fusion", "rearrangement", NA)))
NF2$mutation[is.na(NF2$mutation)] <- ""

#Log transform data
dat_log <- as.data.frame(lapply(dat, function(x) log2(as.numeric(as.character(x))) ))
colnames(dat_log) <- gsub("X", "", colnames(dat_log))
row.names(dat_log) <- row.names(dat)

#scale data average = 0 variant =1
sprDat <- t(scale(t(dat_log))) %>% as.data.frame()
str(sprDat, max.level = 0, give.attr = FALSE)

round(data.frame(avgBefore = rowMeans(head(dat_log)), avgAfter = rowMeans(head(sprDat)),
                 varBefore = apply(head(dat_log), 1, var), varAfter = apply(head(sprDat), 1, var)), 2)

dat_scaled <- sprDat %>% as.data.frame()

#Take the top 5000 most variable by MAD
mads <- apply(dat_scaled, 1, mad)
dat_top <- dat_scaled[rev(order(mads))[1:5000],] %>% as.data.frame()

#set max and min for clipped heatmap
dat_top[] <- lapply(dat_top, function(x) ifelse(x > 3, 3, x))
dat_top[] <- lapply(dat_top, function(x) ifelse(x < -3, -3, x))

# Convert to matrix
dat_top <- as.matrix(dat_top)

### Remove extra large data
rm(dat, sprDat, dat_log)

# Set annotation variables
data_labels <- as.matrix(master$sty_internal_tumour_code)
row.names(data_labels) <- master$x850k_file

samples$Sex <- factor(samples$Sex, levels = c("male", "female"),
                      labels = c("Male", "Female"))
data_sex <- as.matrix(samples$Sex)
row.names(data_sex) <- master$x850k_file

samples$Type <- factor(samples$Type, levels = c("Primary", "Recurrance"))
data_type <- as.matrix(samples$Type)
row.names(data_type) <- master$x850k_file

samples$Subtype <- factor(samples$Subtype, levels = c("transitional", "fibrous", "chordoid", "meningothelial"),
                          labels = c("Transitional", "Fibrous", "Chordoid", "Meningothelial"))
data_subtype <- as.matrix(samples$Subtype)
row.names(data_subtype) <- master$x850k_file

samples$Grade <- factor(samples$Grade, levels = c("I", "II", "III"))
data_grade <- as.matrix(samples$Grade)
row.names(data_grade) <- master$x850k_file

samples$Recurrance <- factor(samples$Recurrance, levels = c("no", "yes"),
                             labels = c("Did Not Recur", "Recurred"))
data_recur <- as.matrix(samples$Recurrance)
row.names(data_recur) <- master$x850k_file

samples$Initial.cancer.diagnosis <- factor(samples$Initial.cancer.diagnosis, levels = c("ALL", "infant ALL", "lymphoblastic lynphoma", "neuroblastoma"),
                                           labels = c("ALL", "Infant ALL", "LBL", "Neuroblastoma"))
data_initial <- as.matrix(samples$Initial.cancer.diagnosis)
row.names(data_initial) <- master$x850k_file

NF2$mutation <- factor(NF2$mutation, levels = c("frameshift", "rearrangement", ""),
                       labels = c("Frameshift", "Rearrangement", ""))
data_NF2 <- as.matrix(NF2$mutation)
row.names(data_NF2) <- master$x850k_file

# Set colors
col_fun <- colorRamp2(c(-3, 0, 3), 
                      c("#313695", "white", "#A50026"))
col_sex <- c(Male = "#a6cee3", Female = "#fb9a99")
col_type <- c(Primary = "#FEE0D2", Recurrance = "#fb9a99")
col_grade <- c(I = "#FED976", II = "#FEB24C", III = "#FD8D3C")
col_subtype <- c(Chordoid = "#A1D99B", Fibrous = "#74C476", Meningothelial = "#41AB5D", Transitional = "#238B45")
col_recur <- c("Did Not Recur" = "grey", "Recurred" = "#fb9a99")
col_initial <- c(ALL = "#A6CEE3", "Infant ALL" = "#B2DF8A", LBL = "#CAB2D6", Neuroblastoma = "#FDBF6F")
col_NF2 <- c(Frameshift = "#238B45", Rearrangement = "#74C476", " " = "grey")

# Set annotation layers
top_annotation <- HeatmapAnnotation(Sex = data_sex,
                                    "Initial Diagnosis" = data_initial,
                                    Subtype = data_subtype,
                                    "WHO Grade" = data_grade,
                                    Recurrance = data_recur,
                                    Timepoint = data_type,
                                    "NF2 Status" = data_NF2,
                                    col = list(Sex = col_sex,
                                               "Initial Diagnosis" = col_initial,
                                               Subtype = col_subtype,
                                               "WHO Grade" = col_grade,
                                               Recurrance = col_recur,
                                               Timepoint = col_type,
                                               "NF2 Status" = col_NF2),
                                    annotation_legend_param = list(title_gp = gpar(fontsize = 15),
                                                                   labels_gp = gpar(fontsize = 13),
                                                                   Subtype = list(labels = c("Transitional", "Fibrous", "Chordoid", "Meningothelial", "NA"))),
                                    annotation_name_gp = gpar(fontsize = 15),
                                    border = TRUE,
                                    simple_anno_size = unit(0.45, "cm"))

### Set Legend Parameters
heatmap_legend_param = list(title = "Log2beta", 
                            border = TRUE,
                            at = c(-3, 0, 3),
                            title_gp = gpar(fontsize = 15),
                            labels_gp = gpar(fontsize = 13),
                            direction = "horizontal",
                            legend_width = unit(4, "cm"))

### Plot and save heatmap
heatmap <- Heatmap(dat_top,
                   col = col_fun,
                   top_annotation = top_annotation,
                   heatmap_legend_param = heatmap_legend_param,
                   show_column_dend = TRUE,
                   show_column_names = FALSE,
                   show_row_dend = FALSE,
                   show_row_names = FALSE,
                   row_title_rot = 0,
                   row_title_gp = gpar(fontsize = 15),
                   border = TRUE,
                   height = unit(4, "inches"))
heatmap

### Get heatmap sample order
heatmap_dend <- column_dend(heatmap)

### Format DKFZ classifier values
row.names(DKFZ) <- DKFZ$sample
DKFZ <- DKFZ[ , 2:7]
DKFZ <- t(DKFZ)

### Set DKFZ scaling colours
col_DKFZ <- colorRamp2(c(0, 1), 
                      c("white", "#A50026"))

### Set DKFZ Row names
row.names(DKFZ) <- c("Malignant", "Benign - 1", "Benign - 2", "Benign - 3", "Intermediate - A", "Intermediate - B")

### Heatmap parameters
heatmap_legend_DKFZ = list(title = "Score", 
                           border = TRUE,
                           title_gp = gpar(fontsize = 15),
                           labels_gp = gpar(fontsize = 13),
                           legend_width = unit(4, "cm"),
                           direction = "horizontal")

### Plot DKFZ heatmap
DKFZ_heatmap <- Heatmap(DKFZ,
                        col = col_DKFZ,
                        cluster_columns = heatmap_dend, 
                        cluster_rows = FALSE,
                        row_title = "",
                        heatmap_legend_param = heatmap_legend_DKFZ,
                        row_names_gp = gpar(fontsize = 15),
                        border = TRUE,
                        height = unit(1.25, "inches"))
DKFZ_heatmap

pdf("methyl_heatmap.pdf", height = 8, width = 7)
draw(heatmap %v% DKFZ_heatmap, heatmap_legend_side = "right", merge_legend = TRUE)
dev.off()

### Save top 5000 differential probes
dat_top <- as.data.frame(dat_top)
dat_top$probe <- row.names(dat_top)
write_tsv(dat_top, file.path(path, "top_5000_probes.tsv"))
