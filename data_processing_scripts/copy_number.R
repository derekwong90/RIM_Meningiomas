library(tidyverse)

### Set working variables ###
path <- "/Users/derekwong/Desktop/H4H/projects/yip_exomes/pipeline_output/"
outdir <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/Yip_Exomes/copy_number"
samples <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/Yip_Exomes/yip_exome_sample_list.txt"
project <- "Yip_exomes"

sequenza <- list.files(file.path(path, "VarScan"), "Sequenza_ratio_gene_matrix.tsv", full.names = TRUE)
gatk <- list.files(file.path(path, "GATK_CNV"), "gatk_ratio_gene_matrix.tsv", full.names = TRUE)

### Make directory files for outputs
dir.create(outdir, showWarnings = FALSE)

### Read in CNAs (Sequenza)
sequenza <- read.delim(sequenza)
gatk <- read.delim(gatk)

### Write tables
write.table(sequenza, file.path(outdir, paste0(project, "_sequenza_cna.txt")), row.names = FALSE, sep = "\t")
write.table(gatk, file.path(outdir, paste0(project, "_gatk_cna.txt")), row.names = FALSE, sep = "\t")