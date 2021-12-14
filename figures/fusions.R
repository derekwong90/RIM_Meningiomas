library(circlize)

### Set working variables ###
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/Yip_Exomes/fusions/fusion_summary.txt"
outdir <- "/Users/derekwong/Google Drive/Post-Doc/Yip_Exomes/Figures/fusions"

### Read in fusions
fusions <- read.delim(path)

### Subset fusions
NF2 <- fusions[grep("NF2", fusions$fusion), ]
other <- fusions[!(fusions$fusion %in% NF2$fusion), ]

### Make bed files for fusion links
NF2_bed1 <- NF2[ , c("donor_chr", "donor_breakpoint", "donor_breakpoint")]
colnames(NF2_bed1) <- c("chromosome", "start", "end")
NF2_bed1$end <- NF2_bed1$end - 10000

NF2_bed2 <- NF2[ , c("acceptor_chr", "acceptor_breakpoint", "acceptor_breakpoint")]
colnames(NF2_bed2) <- c("chromosome", "start", "end")
NF2_bed2$end <- NF2_bed2$end - 10000

other_bed1 <- other[ , c("donor_chr", "donor_breakpoint", "donor_breakpoint")]
colnames(other_bed1) <- c("chromosome", "start", "end")
other_bed1$end <- other_bed1$end - 10000

other_bed2 <- other[ , c("acceptor_chr", "acceptor_breakpoint", "acceptor_breakpoint")]
colnames(other_bed2) <- c("chromosome", "start", "end")
other_bed2$end <- other_bed2$end - 10000

### Set colors
col_NF2 <- "red"
col_other <- "black"

### Create Circos plot
pdf(file.path(outdir, "fusions.pdf"), height = 3.5, width = 3.5)
circos.initializeWithIdeogram(plotType = c("ideogram", "labels"))
#text(-0.25, -0.25, "NF2 Fusions", col = col_NF2, cex = 1)
#text(-0.25, -0.35, "Other Fusions", col = col_other, cex = 1)
circos.genomicLink(other_bed1, other_bed2, col = col_other, border = col_other)
circos.genomicLink(NF2_bed1, NF2_bed2, col = col_NF2, border = col_NF2)
text(-0.32, 0.5, "NF2::ASPG\nNF2:CEP170B", col = col_NF2, cex = 0.75)
text(-0.3, 0.1, "NF2::BAG3", col = col_NF2, cex = 0.75)
text(0.25, -0.35, "TTC32::NF2", col = col_NF2, cex = 0.75)
dev.off()

