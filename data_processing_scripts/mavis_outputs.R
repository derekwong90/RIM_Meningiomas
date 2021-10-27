library(tidyverse)

### Set working variables ###
path <- "/Users/derekwong/Desktop/H4H/projects/yip_exomes/pipeline_output/Mavis"
outdir <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/Yip_Exomes/mavis_sv"
samples <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/Yip_Exomes/yip_exome_sample_list.txt"
project <- "Yip_exomes"

### Make directory files for outputs
dir.create(outdir, showWarnings = FALSE)

### Read in Structural Variants (Mavis)
samples <- read.delim(samples)
mavis <- list.files(path = file.path(path),  pattern = "*mavis_output.tsv", full.names = TRUE)
mavis <- read_tsv(mavis)
mavis <- mavis[ , c(1:18)]

### Remove failed samples
mavis$sample <- gsub("-.*", "", mavis$library)
mavis <- mavis[mavis$sample %in% samples$STY_code, ]
mavis <- mavis[ , 1:18]

### Remove multi-supported SVs
mavis_multi <- mavis[grepl(";", mavis$break1_split_reads) |
                       grepl(";", mavis$break2_split_reads) | 
                       grepl(";", mavis$spanning_reads) |
                       grepl(";", mavis$flanking_pairs), ]
mavis <- mavis[!(grepl(";", mavis$break1_split_reads) |
                   grepl(";", mavis$break2_split_reads) | 
                   grepl(";", mavis$spanning_reads) |
                   grepl(";", mavis$flanking_pairs)), ]
mavis_multi$break1_split_reads <- sapply(strsplit(as.character(mavis_multi$break1_split_reads), '[;"]'), function(x) 
  if(all(is.na(as.numeric(x)))) NA else max(as.numeric(x), na.rm = TRUE))
mavis_multi$break2_split_reads <- sapply(strsplit(as.character(mavis_multi$break2_split_reads), '[;"]'), function(x) 
  if(all(is.na(as.numeric(x)))) NA else max(as.numeric(x), na.rm = TRUE))
mavis_multi$flanking_pairs <- sapply(strsplit(as.character(mavis_multi$flanking_pairs), '[;"]'), function(x) 
  if(all(is.na(as.numeric(x)))) NA else max(as.numeric(x), na.rm = TRUE))
mavis <- rbind(mavis, mavis_multi)

### Format numbers
mavis$break1_split_reads <- as.numeric(mavis$break1_split_reads)
mavis$break2_split_reads <- as.numeric(mavis$break2_split_reads)
mavis$flanking_pairs <- as.numeric(mavis$flanking_pairs)

### Seperate out NF2
mavis_NF2 <- mavis[mavis$gene1_aliases == "NF2" |
                     mavis$gene2_aliases == "NF2", ]

### Filter out germline SVs
mavis_germline <- mavis[grepl("normal", mavis$library), ]
mavis <- mavis[!(mavis$break1_chromosome %in% mavis_germline$break1_chromosome &
                   mavis$break1_position_start %in% mavis_germline$break1_position_start &
                   mavis$break1_position_end %in% mavis_germline$break1_position_end &
                   mavis$break2_chromosome %in% mavis_germline$break2_chromosome &
                   mavis$break2_position_start %in% mavis_germline$break2_position_start &
                   mavis$break2_position_end %in% mavis_germline$break2_position_end), ]

### Remove Recurrent artifacts
mavis_artifacts <- mavis[duplicated(mavis[5:10]),]
mavis <- mavis[!(mavis$break1_chromosome %in% mavis_artifacts$break1_chromosome &
                   mavis$break1_position_start %in% mavis_artifacts$break1_position_start &
                   mavis$break1_position_end %in% mavis_artifacts$break1_position_end &
                   mavis$break2_chromosome %in% mavis_artifacts$break2_chromosome &
                   mavis$break2_position_start %in% mavis_artifacts$break2_position_start &
                   mavis$break2_position_end %in% mavis_artifacts$break2_position_end), ]

### Filter out low confidence and uninteresting
mavis <- mavis[!(mavis$break1_split_reads < 1 |
                   mavis$break1_split_reads < 1), ]
mavis <- mavis[mavis$flanking_pairs >= 10, ]
mavis <- mavis[mavis$linking_split_reads >= 5, ]
mavis <- mavis[!(mavis$event_type %in% c("insertion", "deletion", "duplication")), ]
mavis <- mavis[!(mavis$fusion_splicing_pattern %in% c("None", "retained intron", "normal", "retained multiple introns")), ]

write.table(mavis, file.path(outdir, paste0(project, "_Mavis_all.txt")), sep = "\t", row.names = FALSE)
write.table(mavis_NF2, file.path(outdir, paste0(project, "_Mavis_NF2.txt")), sep = "\t", row.names = FALSE)

