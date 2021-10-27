library(tidyverse)
library(plyr)

# Set working variables
path <- "/Users/derekwong/Desktop/H4H/projects/yip_exomes/pipeline_output/alignment/BAMQC"
outdir <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/Yip_Exomes/sequencing_qc"
samples <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/Yip_Exomes/yip_exome_sample_list.txt"
project <- "Yip_exomes"

# Set sequencing QC file paths and files
dir.create(outdir, showWarnings = FALSE)
cov_path <- file.path(path, "Coverage")
seq_path <- file.path(path, "SequenceMetrics")
coverage <- list.files(cov_path, "Coverage_summary.tsv", full.names = TRUE)
alignment <- list.files(seq_path, "AlignmentMetrics.tsv", full.names = TRUE)

# Read in files
samples <- read.delim(samples)
coverage <- read.delim(coverage)
alignment <- read.delim(alignment)

# Format files
coverage$sample <- gsub("_.*", "", coverage$X)
coverage <- coverage[coverage$sample %in% samples$STY_code, ]
coverage <- coverage[ , 1:7]

alignment <- alignment[alignment$LEVEL == "ALL_READS" & alignment$CATEGORY == "PAIR",]
alignment$sample <- gsub("_.*", "", alignment$SAMPLE)
alignment <- alignment[alignment$sample %in% samples$STY_code, ]
alignment <- alignment[ , 1:26]

# Write tables
write.table(coverage, file.path(outdir, paste0(project, "_coverage.txt")), row.names = FALSE, sep = "\t")
write.table(alignment, file.path(outdir, paste0(project, "_alignment.txt")), row.names = FALSE, sep = "\t")
