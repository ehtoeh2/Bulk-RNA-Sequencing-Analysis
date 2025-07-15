
#-------------------------------ORA analysis-----------------------------
library(clusterProfiler)
library(org.Mm.eg.db)  
library(dplyr)
library(ggplot2)


# Extract significantly upregulated genes (log2FC > 1, FDR < 0.05)
deg_up <- res_df1 %>%
  filter(log2FoldChange > 1 & padj < 0.05) %>%
  pull(gene_symbol)

# Extract significantly downregulated genes (log2FC < -1, FDR < 0.05)
deg_down <- res_df1 %>%
  filter(log2FoldChange < -1 & padj < 0.05) %>%
  pull(gene_symbol)


# GO enrichment analysis (Biological Process) for upregulated genes
go_up <- enrichGO(gene = deg_up,
                  OrgDb = org.Mm.eg.db,
                  keyType = "SYMBOL",
                  ont = "BP",  # or MF, CC
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2)

# Filter significant GO terms by adjusted p-value
sig_go_up <- go_up@result %>% 
  filter(p.adjust < 0.05)



# GO enrichment analysis for downregulated genes
go_down <- enrichGO(gene = deg_down,
                    OrgDb = org.Mm.eg.db,
                    keyType = "SYMBOL",
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2)

# Filter significant GO terms by adjusted p-value
sig_go_down <- go_down@result %>% 
  filter(p.adjust < 0.05)



# Global GO dotplot visualization
dotplot(go_up, showCategory = 15, font.size = 12, title = "GO terms (Up)") + theme_minimal()
dotplot(go_down, showCategory = 15, font.size = 12, title = "GO terms (Down)") + theme_minimal()

# Barplot visualization
barplot(go_up, showCategory = 15, font.size = 12, title = "GO terms (Up)")
barplot(go_down, showCategory = 15, font.size = 12, title = "GO terms (Down)")


# --- Subset for selected GO terms (optional) ---

# Define specific GO term IDs you are interested in
wanted_ID_up <- c(
  "GO:0000000", # ????
  "GO:0000000",
  "GO:0000000",
  "GO:0000000",
  "GO:0000000",
  "GO:0000000",
  "GO:0000000",
  "GO:0000000",
  "GO:0000000",
  "GO:0000000"
)

# Copy the full enrichment result object
go_up_sel <- go_up

# Filter the @result slot to only keep selected GO terms
go_up_sel@result <- go_up@result[
  go_up@result$ID %in% wanted_ID_up,
]

# Dotplot for selected GO terms
dotplot(
  go_up_sel,
  showCategory = length(wanted_ID_up),
  font.size    = 12,
  title        = "GO BP (????)"
) + theme_minimal()


# Create Barplot for selected GO terms
barplot(go_up_sel, showCategory = 7, font.size = 10 ,title = "GO BP UP (???)") + 
  scale_fill_gradient(low = "red3", high = "#FADBD8")  





# 1) 원하는 15개 GO term ID (down)
wanted_ID_down <- c(
  "GO:0000000", # ????
  "GO:0000000",
  "GO:0000000",
  "GO:0000000",
  "GO:0000000",
  "GO:0000000",
  "GO:0000000",
  "GO:0000000",
  "GO:0000000",
  "GO:0000000"
)

# Copy the full enrichment result object
go_down_sel <- go_down

# Filter the @result slot to only keep selected GO terms
go_down_sel@result <- go_down@result[
  go_down@result$ID %in% wanted_ID_down,
]


# Dotplot for selected GO terms
dotplot(
  go_down_sel,
  showCategory = length(wanted_ID_down),
  font.size    = 12,
  title        = "GO BP (Muscle Atrophy)"
) + theme_minimal()


# Create Barplot for selected GO terms
barplot(go_down_sel, showCategory = 10, font.size = 10, title = "GO BP DOWN (???)") + 
  scale_fill_gradient(low = "blue3", high = "#D6EAF8")  


