#!/usr/bin/env Rscript
# R script for DE analysis, GO and KEGG analyses with the R package 'clusterProfiler', visualizations for the most expressed genes (Volcano graph, GO and KEGG graphs)
# (download and) open file <OUTPUT_DIR_OF_RNA-SEQ.SBATCH>/star_salmon/deseq2_qc/deseq2.dds.RData

library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(tibble)
library(dplyr)
library(EnhancedVolcano)
library(clusterProfiler)
library(enrichplot)
library(data.table)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(DOSE)

# define design (here it's 2 factors being accounted for: cell line and treatment groups (C for Control samples and T for treated samples))
# the current design(dds) function is made to answer: What genes change due to Treatment (T vs C), after adjusting for differences between cell lines?
dds$CellLine  <- factor(dds$CellLine)
dds$Treatment <- factor(dds$Treatment, levels = c("C", "T"))
design(dds) <- ~ CellLine + Treatment

# Run DESeq2
dds <- DESeq(dds)

# Extract results (Treated vs Control)
res_tvk <- results(dds, contrast = c("Treatment", "T", "C"))

res_tvk_df <- as.data.frame(res_tvk) %>%
  rownames_to_column(var = "ENSEMBL")

# Clean ENSEMBL IDs
res_tvk_df$ENSEMBL_clean <- sub("\\..*", "", res_tvk_df$ENSEMBL)

# Add gene symbols and Entrez IDs
res_tvk_df$SYMBOL <- mapIds(
  org.Hs.eg.db,
  keys = res_tvk_df$ENSEMBL_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
res_tvk_df$ENTREZID <- mapIds(
  org.Hs.eg.db,
  keys = res_tvk_df$ENSEMBL_clean,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

sig_all <- res_tvk_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1.0)

# Save annotated table to .csv file
write.csv(sig_all, "Treatment_vs_Control_sig_annotated.csv", row.names = FALSE)

# -------------------------------------------------
# Graphs (MA, PCA, Volcano)

# MA plot ylim = c(-3, 3))
plotMA(res_tvk,
       ylim = c(-5, 5),
       main = "MA diagram: Differential gene expression",
       xlab = "Average normalized expression (log2)",
       ylab = "Change in log2 expression")

# PCA
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = c("CellLine", "Treatment")) +
  ggtitle("PCA diagram: Differences between sample groups") +
  labs(x="Principal Component 1", y="Principal Component 2", color="Group") +
  theme_bw() +
  theme(plot.title = element_text(face="bold", size=16))

# Volcano plot
res_sorted <- sig_all %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1.0) %>%
  arrange(desc(abs(log2FoldChange)))
  
top5_rna <- res_sorted %>%
     head(5)

top5_symbols <- top5_rna %>%
  pull(SYMBOL)

top5_ids <- top5_rna$ENSEMBL_clean

volcano_data <- res_tvk_df %>%
  filter(!is.na(padj), !is.na(log2FoldChange)) %>%
  mutate(
    minuslog10padj = -log10(padj),
    gene_category = case_when(
      ENSEMBL_clean %in% top5_ids ~ "TOP 5 RNA",
      padj < 0.05 & abs(log2FoldChange) > 1.0 ~ "Significant",
      TRUE ~ "Not significant"
    ),
    label = ifelse(ENSEMBL_clean %in% top5_ids,
                  ifelse(!is.na(SYMBOL), SYMBOL, ENSEMBL_clean),
                  "")
  )

# Order for plotting (TOP 5 on top)
volcano_data$gene_category <- factor(volcano_data$gene_category,
                                     levels = c("Not significant", "Significant",
                                               "TOP 5 RNA"))

pdf("Volcano_TOP5_RNA.pdf", width = 12, height = 10)
ggplot(volcano_data, aes(x = log2FoldChange, y = minuslog10padj)) +
    geom_point(aes(color = gene_category, size = gene_category, alpha = gene_category)) +
    scale_color_manual(values = c("Not significant" = "gray70",
                                  "Significant" = "orange",
                                  "TOP 5 RNA" = "red")) +
    scale_size_manual(values = c("Not significant" = 0.5,
                                 "Significant" = 1,
                                 "TOP 5 RNA" = 4)) +
    scale_alpha_manual(values = c("Not significant" = 0.3,
                                  "Significant" = 0.5,
                                  "TOP 5 RNA" = 1)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", linewidth = 0.5) +
    geom_vline(xintercept = c(-1.0, 1.0), linetype = "dashed", color = "blue", linewidth = 0.5) +
    geom_text_repel(aes(label = label),
                    size = 5,
                    fontface = "bold",
                    box.padding = 1,
                    point.padding = 0.5,
                    max.overlaps = Inf,
                    check_overlap = FALSE,
                    force = 5,
                    min.segment.length = 0,
                    segment.color = "red",
                    segment.size = 0.5) +
    labs(title = "Volcano Plot: Differential Expression Analysis\n(TOP 5 lncRNAs Highlighted in Red)",
         subtitle = paste0("Total significant genes: ", nrow(sig_all)),
         x = "Log2 Fold Change (TMZ vs Control)",
         y = "-Log10 Adjusted P-value",
         color = "Gene Category") +
    guides(size = "none", alpha = "none") +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          legend.position = "right",
          legend.title = element_text(face = "bold"),
          panel.grid.major = element_line(color = "gray90"))
dev.off()

# --------------------------------------------------
# GO Enrichment

sig_genes <- res_tvk_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1.0) %>%
  pull(ENTREZID) %>%
  na.omit()

ego <- enrichGO(
  gene          = sig_all$SYMBOL,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

ego_sig <- ego@result %>%
    filter(p.adjust < 0.05 & qvalue < 0.05) %>%
    arrange(desc(abs(Count)))

# save top5 terms based on gene count for simplification of graphs
top5_terms <- head(ego_sig$Description, 5)

# save .csv file of GO significant pathways
fwrite(as.data.frame(ego_sig), "GO_BP_enrichment.csv")
 
# Dotplot - TOP 5
pdf("GO_BP_dot_RNA.pdf", width = 10, height = 6)
print(dotplot(ego, showCategory = top5_terms, font.size = 12) +
     ggtitle("GO Biological Process Enrichment\n(TOP 5 Pathways)") +
     theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16)))
dev.off()
 
# Cnetplot - TOP 5
pdf("GO_BP_cnet_RNA.pdf", width = 12, height = 10)
print(cnetplot(ego, showCategory = top5_terms,
             node_label = "all",
             cex_label_gene = 0.7,
             cex_label_category = 1.2) +
      ggtitle("GO BP Gene-Concept Network\n(TOP 5 Pathways)") +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16)))
dev.off()

# --------------------------------------------------
# KEGG Enrichment

ekegg <- enrichKEGG(
  gene          = sig_all$ENTREZID,
  organism      = "hsa",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

ekegg_sig <- ekegg@result %>%
    filter(p.adjust < 0.05 & qvalue < 0.05) %>%
    arrange(desc(abs(Count)))

top5_kegg <- head(ekegg_sig$Description, 5)

fwrite(as.data.frame(ekegg_sig), "KEGG_RNA_enrichment.csv")

# Dotplot - TOP 5
pdf("KEGG_dot_RNA.pdf", width = 10, height = 6)
print(dotplot(ekegg, showCategory = top5_kegg, font.size = 12) +
      ggtitle("KEGG Pathway Enrichment\n(TOP 5 Pathways)") +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16)))
dev.off()
 
# Cnetplot - TOP 5
pdf("KEGG_cnet_RNA.pdf", width = 12, height = 10)
print(cnetplot(ekegg, showCategory = top5_kegg,
               node_label = "all",
               cex_label_gene = 0.7,
               cex_label_category = 1.2) +
      ggtitle("KEGG Gene-Concept Network\n(TOP 5 Pathways)") +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16)))
dev.off()
