###############################################################################
# Transcriptomic Analysis - Functional Enrichment Analysis (GO + KEGG)
#
# Description: Functional enrichment analysis of differentially expressed
# genes (DEGs) identified in the 2cpab vs WT comparison (DESeq2). Includes
# Gene Ontology (Biological Process) enrichment with gene-concept network
# (cnetplot) visualization, export of GO terms for REVIGO-based grouping by
# semantic similarity, and KEGG pathway enrichment for activated and
# repressed gene sets.
#
# Organism: Arabidopsis thaliana
#
# Input:
#   tables/activated_genes_deseq2.txt   : list of activated DEGs (DESeq2)
#   tables/repressed_genes_deseq2.txt   : list of repressed DEGs (DESeq2)
#   tables/full_results_deseq2.tsv      : complete DESeq2 results table
#
# Output:
#   images/cnetplot_activated.png     : GO BP gene-concept network, activated
#   images/cnetplot_repressed.png     : GO BP gene-concept network, repressed
#   images/KEGG_activated_dotplot.png : KEGG pathway dotplot, activated genes
#   images/KEGG_repressed_dotplot.png : KEGG pathway dotplot, repressed genes
#   tables/GO_activated_deseq2.tsv    : full GO enrichment table, activated
#   tables/GO_repressed_deseq2.tsv    : full GO enrichment table, repressed
#   tables/KEGG_activated_deseq2.tsv  : full KEGG enrichment table, activated
#   tables/KEGG_repressed_deseq2.tsv  : full KEGG enrichment table, repressed
#   tables/GO_activated_revigo.txt    : GO IDs + p.adjust, formatted for REVIGO
#   tables/GO_repressed_revigo.txt    : GO IDs + p.adjust, formatted for REVIGO
#
# Author: Nardy Celeste Ayala Pérez
# Date: 2026
###############################################################################

#-------------------------------------------------------------------------------
# PACKAGE LOADING
#-------------------------------------------------------------------------------
library(clusterProfiler)  # enrichGO, enrichKEGG
library(org.At.tair.db)   # Arabidopsis thaliana genome annotation
library(enrichplot)       # cnetplot, dotplot visualization
library(ggplot2)          # Required by ggsave() for saving enrichplot figures

dir.create("tables", showWarnings = FALSE)
dir.create("images", showWarnings = FALSE)

#-------------------------------------------------------------------------------
# LOAD DEG LISTS AND BACKGROUND UNIVERSE FROM RNAseq_01_DEGAnalysis.R OUTPUT
#-------------------------------------------------------------------------------

# Load DEG lists generated in RNAseq_01_DEGAnalysis.R
activated.genes.deseq2 <- read.table("tables/activated_genes_deseq2.txt")$V1
repressed.genes.deseq2 <- read.table("tables/repressed_genes_deseq2.txt")$V1

# Load full DESeq2 results table, needed to define gene.ids.deseq2 as the
# background gene universe for GO enrichment
res.deseq2.df   <- read.table("tables/full_results_deseq2.tsv", 
                              sep = "\t", header = TRUE)
gene.ids.deseq2 <- rownames(res.deseq2.df)

## -----------------------------------------------------------------------------
## GO ENRICHMENT — ACTIVATED GENES (DESeq2)
## -----------------------------------------------------------------------------

# The background universe is set to all genes tested by DESeq2, ensuring that 
# enrichment significance is calculated relative to genes that were actually 
# detected and tested in this experiment

universe <- gene.ids.deseq2  

activated.enrich.go <- enrichGO(
  gene          = activated.genes.deseq2,
  universe      = universe,
  OrgDb         = org.At.tair.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE,
  keyType       = "TAIR"
)

df.act.go <- as.data.frame(activated.enrich.go)
head(df.act.go)

# Gene-concept network: connects enriched GO terms to their associated genes,
# revealing shared functional modules among the activated DEGs

cnetplot(activated.enrich.go, showCategory = 20)
ggsave("images/cnetplot_activated.png", cnetplot(activated.enrich.go, 
          showCategory = 20), width = 12, height = 10, dpi = 600)

## -----------------------------------------------------------------------------
## GO ENRICHMENT — REPRESSED GENES (DESeq2)
## -----------------------------------------------------------------------------

repressed.enrich.go <- enrichGO(
  gene          = repressed.genes.deseq2,
  universe      = universe,
  OrgDb         = org.At.tair.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE,
  keyType       = "TAIR"
)
df.rep.go <- as.data.frame(repressed.enrich.go)
head(df.rep.go)

cnetplot(repressed.enrich.go, showCategory = 20)
ggsave("images/cnetplot_repressed.png",
       cnetplot(repressed.enrich.go, showCategory = 20),
       width = 12, height = 10, dpi = 600)

## -----------------------------------------------------------------------------
## KEGG ENRICHMENT — (DESeq2)
## KEGG pathway enrichment using the Arabidopsis thaliana organism code "ath".
## -----------------------------------------------------------------------------

# Activated pathways

activated.enrich.kegg <- enrichKEGG(
  gene          = activated.genes.deseq2,
  organism      = "ath",
  keyType       = "kegg",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)
df.act.kegg <- as.data.frame(activated.enrich.kegg)

if (!is.null(activated.enrich.kegg) && nrow(df.act.kegg) > 0) {
  ggsave("images/KEGG_activated_dotplot.png",
         dotplot(activated.enrich.kegg, showCategory = 20,
                 title = "KEGG — Activated genes (DESeq2)"),
         width = 8, height = 6, dpi = 300)
} else {
  message("No enriched KEGG pathways found for activated genes.")
}

# Repressed pathways

repressed.enrich.kegg <- enrichKEGG(
  gene          = repressed.genes.deseq2,
  organism      = "ath",
  keyType       = "kegg",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)
df.rep.kegg <- as.data.frame(repressed.enrich.kegg)

if (!is.null(repressed.enrich.kegg) && nrow(df.rep.kegg) > 0) {
  ggsave("images/KEGG_repressed_dotplot.png",
         dotplot(repressed.enrich.kegg, showCategory = 20,
                 title = "KEGG — Repressed genes (DESeq2)"),
         width = 8, height = 6, dpi = 300)
} else {
  message("No enriched KEGG pathways found for repressed genes.")
}

dotplot(activated.enrich.kegg, showCategory = 20,
        title = "KEGG — Activated genes (DESeq2)")
dotplot(repressed.enrich.kegg, showCategory = 20,
        title = "KEGG — Repressed genes (DESeq2)")

## -----------------------------------------------------------------------------
## TOP 10 ACTIVATED AND REPRESSED GENES BY LOG2 FOLD CHANGE
## -----------------------------------------------------------------------------

log.fc.deseq2 <- res.deseq2.df$log2FoldChange
adj.p.deseq2  <- res.deseq2.df$padj
names(log.fc.deseq2) <- gene.ids.deseq2
names(adj.p.deseq2)  <- gene.ids.deseq2

# Top 10 activated genes (highest expression in 2cpab), ordered by log2FC
top10.act <- head(activated.genes.deseq2[
  order(log.fc.deseq2[activated.genes.deseq2], decreasing = TRUE)], 10)

# Top 10 repressed genes (lowest expression in 2cpab), ordered by log2FC
top10.rep <- head(repressed.genes.deseq2[
  order(log.fc.deseq2[repressed.genes.deseq2], decreasing = FALSE)], 10)

cat("TOP 10 ACTIVATED IN 2cpab")
data.frame(
  Gene   = top10.act,
  log2FC = round(log.fc.deseq2[top10.act], 3),
  adj.p  = round(adj.p.deseq2[top10.act], 6)
)

cat("TOP 10 REPRESSED IN 2cpab")
data.frame(
  Gene   = top10.rep,
  log2FC = round(log.fc.deseq2[top10.rep], 3),
  adj.p  = round(adj.p.deseq2[top10.rep], 6)
)

## -----------------------------------------------------------------------------
## EXPORT GO TERMS FOR REVIGO
## GO IDs and adjusted p-values are exported in the space-separated format
## required by REVIGO (http://revigo.irb.hr/), used to group enriched terms
## by semantic similarity and generate a more structured visualization.
## -----------------------------------------------------------------------------

go.revigo.activated <- as.data.frame(activated.enrich.go)[, c("ID", "p.adjust")]

go.revigo.repressed <- as.data.frame(repressed.enrich.go)[, c("ID", "p.adjust")]

## -----------------------------------------------------------------------------
## SAVE ALL RESULTS TABLES
## -----------------------------------------------------------------------------

write.table(df.act.go, file = "tables/GO_activated_deseq2.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(df.rep.go, file = "tables/GO_repressed_deseq2.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(df.act.kegg, file = "tables/KEGG_activated_deseq2.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(df.rep.kegg, file = "tables/KEGG_repressed_deseq2.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(go.revigo.activated, "tables/GO_activated_revigo.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(go.revigo.repressed, "tables/GO_repressed_revigo.txt", sep = " ",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# End of script