
#-------------------------------GSEA analysis-----------------------------
library(enrichplot)

# Extract log2 fold change values
original_gene_list <- res_df1$log2FoldChange

names(original_gene_list) <- res_df1$gene_symbol

gene_list<-na.omit(original_gene_list)

# Sort the gene list in decreasing order (required for GSEA)
gene_list = sort(gene_list, decreasing = TRUE)


# Define organism database (Mus musculus)
organism = "org.Mm.eg.db" 

# Run GSEA using Gene Ontology Biological Process terms
gse <- gseGO(geneList      = gene_list, 
             ont          = "BP", 
             keyType      = "SYMBOL", 
             minGSSize    = 3, 
             maxGSSize    = 800, 
             pvalueCutoff = 0.05, 
             verbose      = TRUE, 
             OrgDb        = organism, )


# Dotplot of top enriched GO terms with color scale by NES
dotplot(gse, showCategory=10, split=".sign",) + 
  scale_color_gradient2(
    low      = "blue",           
    mid      = "white",
    high     = "red",            
    midpoint = 0
  ) + facet_grid(.~.sign)

# Ridge plot showing enrichment score distributions
ridgeplot(gse_sel, showCategory = 7, fill = "NES") +
  scale_fill_gradient2(
    low      = "blue3",    
    mid      = "white",   
    high     = "#bb0c00", 
  ) +  labs(x = "enrichment distribution") +
  theme_minimal(base_size = 10)+
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()   
  )



# --- Subset for selected GSE results (optional) ---

# Define a vector of GO term IDs to visualize
wanted_ID_gse <- c(
  "GO:000000",
  "GO:000000",
  "GO:000000",
  "GO:000000",
  "GO:000000"
)


# Copy GSEA result and filter by selected GO terms
gse_sel <- gse
gse_sel@result <- gse@result[gse@result$ID %in% wanted_ID_gse, ]


# Dotplot visualization of selected GO terms
dotplot(gse_sel, showCategory=5, split=".sign", font.size = 10) + facet_grid(.~.sign) +
  labs(title = "GSEA GO (???)") +
  guides(color = guide_colorbar(title = "p.adjust"))


# Ridgeplot visualization with NES coloring
ridgeplot(gse_sel, showCategory = 7, fill = "NES") +
  scale_fill_gradient2(
    low      = "blue3",    
    mid      = "white",   
    high     = "#bb0c00", 
  ) +  labs(x = "enrichment distribution") +
  theme_minimal(base_size = 10)+
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()   
  )




# ------------------ Single Pathway Enrichment Plot ------------------


# Specify the GO term of interest
my_pathway <- "GO:0000000"

# Get the index of the pathway in the GSEA result
idx <- which(gse_sel@result$ID == my_pathway)

# Extract NES value and round to 2 decimal places
NES_val <- round(gse_sel@result$NES[idx], 2)

#Extract adjusted p-value
p.adjust <- round(gse_sel@result$p.adjust[idx])


# 4) enrichment plot 그리기
gp <- gseaplot2(
  gse_sel,
  geneSetID = idx,
  title     = paste0(
    gse_sel@result$Description[idx], 
    " (", my_pathway, ")"
  ),
  color     = "green"       
)

gp





