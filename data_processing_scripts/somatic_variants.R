library(tidyverse)
library(dplyr)
library(janitor)
library(purrr)
library(data.table)

### Set working variables ###
path <- "/Users/derekwong/Desktop/H4H/projects/yip_exomes/pipeline_output/Report/plots"
bed <- "/Users/derekwong/Desktop/H4H/projects/yip_exomes/intervals/IDT_xGen_exome_v1_targets_hg38_padding100bp.bed"
outdir <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/Yip_Exomes/somatic"
project <- "Yip_exomes"

somatic_dir <- file.path(path, "ensemble_mutation_data.tsv")

### Make directory files for outputs
dir.create(outdir, showWarnings = FALSE)

### Read in Files (Ensemble) ###
bed <- read.delim(bed, header = FALSE)
somatics <- read_tsv(somatic_dir)
somatic <- somatics

### Remove failed samples and rename sample T7
somatic <- somatic[!(somatic$Tumor_Sample_Barcode %in% c("T8_primary", "T14_primary", "T16b_primary")), ]
somatic$Tumor_Sample_Barcode <- gsub("T7_recurrance", "T7_primary", somatic$Tumor_Sample_Barcode)

### Rough filtering
somatic <- remove_empty(somatic, which = c("cols"), quiet = TRUE)
somatic <- somatic[!(somatic$Chromosome == "chrY"), ]
somatic <- somatic[somatic$FLAG.low_vaf == "FALSE", ]
somatic <- somatic[somatic$FLAG.low_coverage == "FALSE", ]

### Filter for Variants in bed file
colnames(bed) <- c("chromosome", "start", "end")
mutations <- somatic[ , 5:7]
mutations <- mutations[complete.cases(mutations), ]
colnames(mutations) <- c("chromosome", "start", "end")

setDT(mutations)
setDT(bed)
setkey(bed)
mutations_bed <- foverlaps(mutations, bed, type="within", nomatch=0L)

somatic <- somatic[row.names(somatic) %in% row.names(mutations_bed), ]
rm(mutations, mutations_bed)

### Filter out variants > 0.01 population AF
somatic[, c("AF", "AFR_AF", "AMR_AF", "EAS_AF", 
            "EUR_AF", "SAS_AF", "AA_AF", "EA_AF")][is.na(somatic[ , c("AF", "AFR_AF", "AMR_AF", "EAS_AF", 
                                                                      "EUR_AF", "SAS_AF", "AA_AF", "EA_AF")])] <- 0
somatic <- somatic[somatic$AF <= 0.01 &
                     somatic$AFR_AF <= 0.01 &
                     somatic$AMR_AF <= 0.01 &
                     somatic$EAS_AF <= 0.01 &
                     somatic$EUR_AF <= 0.01 &
                     somatic$SAS_AF <= 0.01 &
                     somatic$AA_AF <= 0.01 &
                     somatic$EA_AF <= 0.01, ]
somatic[, c("gnomAD_AF", "gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_ASJ_AF", "gnomAD_EAS_AF", "gnomAD_FIN_AF", 
            "gnomAD_NFE_AF", "gnomAD_OTH_AF", "gnomAD_SAS_AF")][is.na(somatic[ , c("gnomAD_AF", "gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_ASJ_AF", "gnomAD_EAS_AF", 
                                                                                   "gnomAD_FIN_AF", "gnomAD_NFE_AF", "gnomAD_OTH_AF", "gnomAD_SAS_AF")])] <- 0
somatic <- somatic[somatic$gnomAD_AF <= 0.01 &
                     somatic$gnomAD_AFR_AF <= 0.01 &
                     somatic$gnomAD_AMR_AF <= 0.01 &
                     somatic$gnomAD_ASJ_AF <= 0.01 &
                     somatic$gnomAD_EAS_AF <= 0.01 &
                     somatic$gnomAD_FIN_AF <= 0.01 &
                     somatic$gnomAD_NFE_AF <= 0.01 &
                     somatic$gnomAD_OTH_AF <= 0.01 &
                     somatic$gnomAD_SAS_AF <= 0.01, ]
somatic <- somatic[somatic$FLAG.high_pop == "FALSE", ]

### Filter out low impact and benign
somatic <- somatic[!grepl("benign", somatic$CLIN_SIG),]
somatic <- somatic[!grepl("benign", somatic$PolyPhen),]
somatic <- somatic[!grepl("tolerated", somatic$SIFT),]

### Calculate VAFs
somatic$vaf <- somatic$t_alt_count/somatic$t_depth*100
somatic <- somatic %>% relocate(vaf, .after = n_alt_count)

### Remove low coverage (>40) and low vaf (>10)
somatic <- somatic[somatic$t_depth >= 40, ]
somatic <- somatic[!((somatic$Variant_Type == "INS" | somatic$Variant_Type == "DEL") &
                     somatic$vaf < 10), ]
somatic <- somatic[somatic$n_alt_count < 1, ]
somatic <- somatic[somatic$t_alt_count >= 5, ]

### Seperate out SNPs
snps <- somatic[somatic$Variant_Type == "SNP", ]

### Keep mutations of interest
keep <- c("Missense_Mutation", "Splice_Site", "In_Frame_Ins", "In_Frame_Del", "Frame_Shift_Ins", 
          "Frame_Shift_Del", "Nonstop_Mutation", "Nonsense_Mutation", "Translation_Start_Site")
somatic <- somatic[somatic$Variant_Classification %in% keep, ]
somatic <- somatic[!(somatic$Variant_Type == "SNP"), ]
somatic <- rbind(somatic, snps)

### NF2 Mutations
somatic_NF2 <- somatic[somatic$Hugo_Symbol == "NF2", ]

### Save files
write.table(somatic, file.path(outdir, paste0(project, "_all_somatic.txt")), sep = "\t", row.names = FALSE)
write.table(snps, file.path(outdir, paste0(project, "_snps.txt")), sep = "\t", row.names = FALSE)
write.table(somatic_NF2, file.path(outdir, paste0(project, "_NF2.txt")), sep = "\t", row.names = FALSE)
