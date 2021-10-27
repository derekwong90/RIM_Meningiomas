library(tidyverse)
library(dplyr)
library(janitor)
library(purrr)

### Set working variables ###
path <- "/Users/derekwong/Desktop/H4H/projects/yip_exomes/pipeline_output/HaplotypeCaller/CPSR"
outdir <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/Yip_Exomes/germline"
samples <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/Yip_Exomes/yip_exome_sample_list.txt"
project <- "Yip_exomes"

### Make directory files for outputs
dir.create(outdir, showWarnings = FALSE)

### Germline Mutations ###
samples <- read.delim(samples)
failed <- c("T14_primary", "T8_primary")
germline <- list.files(path, "*mutations_for_cbioportal.tsv", full.names = TRUE)
germline <- read_tsv(germline)
germline <- germline[!(germline$Tumor_Sample_Barcode %in% failed), ]
germline$Tumor_Sample_Barcode <- gsub("_.*", "", germline$Tumor_Sample_Barcode)
germline <- germline[germline$Tumor_Sample_Barcode %in% samples$STY_code, ]
germline <- remove_empty(germline, which = c("cols"), quiet = TRUE)
keep <- c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Variant_Classification", "Variant_Type", "Reference_Allele",
          "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "dbSNP_RS", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "Match_Norm_Seq_Allele1",
          "Match_Norm_Seq_Allele2", "HGVSc", "HGVSp_Short", "t_depth", "t_ref_count", "t_alt_count", "n_depth", "n_ref_count", "n_alt_count",
          "SIFT", "PolyPhen", "CLIN_SIG")
germline <- germline[, colnames(germline) %in% keep]
germline <- germline[!grepl("benign", germline$SIFT) &
                     !grepl("benign", germline$PolyPhen) &
                     !grepl("benign", germline$CLIN_SIG), ]
germline$t_vaf <- germline$t_alt_count/germline$t_depth*100
germline$n_vaf <- germline$n_alt_count/germline$n_depth*100
germline <- germline %>% relocate(t_vaf, .after = n_alt_count)
germline <- germline %>% relocate(n_vaf, .after = t_vaf)
write.table(germline, file.path(outdir, paste0(project, "_haplotypecaller.txt")), sep = "\t", row.names = FALSE)
