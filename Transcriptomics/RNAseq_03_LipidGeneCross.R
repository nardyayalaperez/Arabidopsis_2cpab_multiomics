###############################################################################
# Transcriptomic Analysis - Lipid Metabolism Gene Cross-Reference
#
# Description: To identify lipid metabolism-related genes among the
# differentially expressed genes (DEGs) identified in the 2cpab vs WT
# comparison (DESeq2), the RNA-seq gene lists were crossed with a curated
# list of Arabidopsis thaliana genes, organized into nine functional
# categories: plastidial fatty acid synthesis, prokaryotic glycerolipid
# synthesis, eukaryotic galactolipid and sulfolipid synthesis, galactolipid
# degradation, eukaryotic phospholipid synthesis, TAG synthesis, TAG
# catabolism, lipases, and lipid signaling.
#
# Organism: Arabidopsis thaliana
#
# Input:
#   tables/activated_genes_deseq2.txt : list of activated DEGs (DESeq2)
#   tables/repressed_genes_deseq2.txt : list of repressed DEGs (DESeq2)
#   Lipid metabolism gene.xlsx        : curated lipid metabolism gene 
#                                       one sheet per pathway category
#
# Output:
#   tables/Lipid_DEGs_cross.xlsx    : full results (activated/repressed
#                                     lipid DEGs and complete lipid gene list)
#   images/Venn_activated_lipid.png : overlap between activated DEGs
#                                     and the curated lipid gene list
#   images/Venn_repressed_lipid.png : overlap between repressed DEGs
#                                     and the curated lipid gene list
#
# Author: Nardy Celeste Ayala Pérez
# Date: 2026
###############################################################################

#-------------------------------------------------------------------------------
# PACKAGE LOADING
#-------------------------------------------------------------------------------
library(readxl)       # Reading the lipid metabolism gene excel database
library(dplyr)        # Data manipulation (group_by, summarise)
library(writexl)      # Exporting results to Excel
library(VennDiagram)  # Venn diagrams
library(grid)         # grid.newpage(), required by VennDiagram

dir.create("tables", showWarnings = FALSE)
dir.create("images", showWarnings = FALSE)

#-------------------------------------------------------------------------------
# LOAD DEG LISTS FROM 01_RNAseq_DEG_analysis.R OUTPUT
#-------------------------------------------------------------------------------

activated.genes.deseq2 <- read.table("tables/activated_genes_deseq2.txt")$V1
repressed.genes.deseq2 <- read.table("tables/repressed_genes_deseq2.txt")$V1

cat("Activated DEGs:", length(activated.genes.deseq2))
cat("Repressed DEGs:", length(repressed.genes.deseq2))

#-------------------------------------------------------------------------------
# LOAD THE CURATED LIPID METABOLISM GENE DATABASE
#-------------------------------------------------------------------------------

# The Excel file contains one sheet per lipid metabolic pathway category.
# Only rows with a valid Arabidopsis locus identifier (e.g., AT1G01010)
# are retained. The pathway category is recorded for each gene.

lipid.gene.file <- "Lipid metabolism gene.xlsx"

lipid.gene.sheets <- c("Plastidial FA Synt",
                   "Prokaryotic Glycerolipid Synt",
                   "Eukaryotic Galacto Sulfo synt ",
                   "Galactolipid degradation",
                   "Eukaryotic Phospholipid Synt",
                   "TAG synthesis",
                   "TAG catabolism",
                   "Lipase",
                   "Lipid signaling")

lipid.genes.all <- data.frame()

for (sheet in lipid.gene.sheets) {
  sheet.data <- read_excel(lipid.gene.file, sheet = sheet)
  colnames(sheet.data)[1] <- "GENE"
  
  # Keep only rows with a valid AGI locus ID (e.g., AT1G01010), removing
  # header rows, notes, or empty cells present in the original spreadsheet
  
  sheet.data <- sheet.data[!is.na(sheet.data$GENE) &
                             grepl("^AT[0-9CM]G[0-9]{5}$", sheet.data$GENE), ]
  
  if (nrow(sheet.data) > 0) {
    sheet.data$Sheet <- sheet
    lipid.genes.all <- rbind(lipid.genes.all, sheet.data[, c("GENE", "Sheet")])
  }
}

# Remove duplicate gene entries that appear in more than one pathway category,
# keeping only the first occurrence
lipid.genes.all <- lipid.genes.all[!duplicated(lipid.genes.all$GENE), ]

cat("Total lipid metabolism genes loaded:", nrow(lipid.genes.all))

#-------------------------------------------------------------------------------
# CROSS-REFERENCE LIPID METABOLISM GENES WITH DEGs
#-------------------------------------------------------------------------------

# Lipid metabolism genes among the activated DEGs
act.lipid <- lipid.genes.all[lipid.genes.all$GENE %in% activated.genes.deseq2, ]
act.lipid$Direction <- "Activated"

# Lipid metabolism genes among the repressed DEGs
rep.lipid <- lipid.genes.all[lipid.genes.all$GENE %in% repressed.genes.deseq2, ]
rep.lipid$Direction <- "Repressed"

# Display results
cat("Lipid metabolism genes activated in 2cpab")
print(act.lipid)

cat("Lipid metabolism genes repressed in 2cpab")
print(rep.lipid)

#-------------------------------------------------------------------------------
# VENN DIAGRAMS: OVERLAP BETWEEN DEGs AND LIPID METABOLISM GENES
#-------------------------------------------------------------------------------

# Visualizes, for each direction (activated/repressed), how many DEGs
# belong to the curated lipid metabolism gene set, out of the total
# number of DEGs and the total number of lipid metabolism genes.

# Activated DEGs vs lipid metabolism genes
png("images/Venn_activated_lipid.png", width = 800, height = 600, res = 150)
grid.newpage()
draw.pairwise.venn(
  area1      = length(activated.genes.deseq2),
  area2      = length(lipid.genes.all$GENE),
  cross.area = length(intersect(activated.genes.deseq2, lipid.genes.all$GENE)),
  lwd        = 3,
  category   = c("Upregulated DEGs", "Lipid metabolism genes"),
  euler.d    = FALSE,
  scaled     = FALSE,
  col        = c("firebrick2", "forestgreen"),
  fill       = c("firebrick2", "forestgreen"),
  alpha      = 0.3,
  cex        = 1.5,
  cat.cex    = 1.5,
  cat.pos    = c(180, 0)
)
dev.off()

# Repressed DEGs vs lipid metabolism genes
png("images/Venn_repressed_lipid.png", width = 800, height = 600, res = 150)
grid.newpage()
draw.pairwise.venn(
  area1      = length(repressed.genes.deseq2),
  area2      = length(lipid.genes.all$GENE),
  cross.area = length(intersect(repressed.genes.deseq2, lipid.genes.all$GENE)),
  lwd        = 3,
  category   = c("Downregulated DEGs", "Lipid metabolism genes"),
  euler.d    = FALSE,
  scaled     = FALSE,
  col        = c("dodgerblue", "forestgreen"),
  fill       = c("dodgerblue", "forestgreen"),
  alpha      = 0.3,
  cex        = 1.5,
  cat.cex    = 1.5,
  cat.pos    = c(180, 0)
)
dev.off()


#-------------------------------------------------------------------------------
# SAVE RESULTS TO A SINGLE EXCEL FILE
#-------------------------------------------------------------------------------

results <- list("Activated_lipid_genes" = act.lipid, 
                "Repressed_lipid_genes" = rep.lipid, 
                "Al_lipid_genes"      = lipid.genes.all)

write_xlsx(results, "tables/Lipid_DEGs_cross.xlsx")

# End of script