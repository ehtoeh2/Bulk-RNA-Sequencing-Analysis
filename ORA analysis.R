
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



#GSEA
# we want the log2 fold change 
res_df1
original_gene_list <- res_df1$log2FoldChange

original_gene_list

# name the vector
names(original_gene_list) <- res_df1$gene_symbol

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order
gene_list = sort(gene_list, decreasing = TRUE)

gene_list

organism = "org.Mm.eg.db" 

gse <- gseGO(geneList      = gene_list, 
             ont          = "BP", 
             keyType      = "SYMBOL", 
             minGSSize    = 3, 
             maxGSSize    = 800, 
             pvalueCutoff = 0.05, 
             verbose      = TRUE, 
             OrgDb        = organism, )

gse
View(gse@result)

dotplot(gse, showCategory=20)

dotplot(gse, showCategory=10, split=".sign",) + 
  scale_color_gradient2(
    low      = "blue",           # suppressed 쪽(음수)
    mid      = "white",
    high     = "red",            # activated 쪽(양수)
    midpoint = 0
  ) + facet_grid(.~.sign)
ridgeplot(gse) + labs(x = "enrichment distribution")


write.xlsx(gse, "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/Gse.xlsx", rownames = TRUE)





#내가 원하는 ID의 GSE 결과만 남겨서 barplot 보여주기 

# 1) 원하는 GO term ID 벡터 정의
wanted_ID_gse <- c(
  "GO:0010498",  # proteasomal protein catabolic process
  "GO:0043161",  # proteasome-mediated ubiquitin-dependent protein catabolic process
  "GO:0006511",  # ubiquitin-dependent protein catabolic process
  "GO:0051603",  # proteolysis involved in protein catabolic process
  "GO:0031329",  # regulation of cellular catabolic process
  "GO:0097190",  # apoptotic signaling pathway
  "GO:0034976",  # response to endoplasmic reticulum stress
  "GO:0000302",  # response to reactive oxygen species
  "GO:0042594",  # response to starvation
  "GO:0016567",   # protein ubiquitination
  "GO:0030198",  # extracellular matrix organization
  "GO:0043062",  # extracellular structure organization
  "GO:0031589",  # cell-substrate adhesion
  "GO:0098742",  # cell-cell adhesion via plasma-membrane adhesion molecules
  "GO:0043408",  # regulation of MAPK cascade
  "GO:0043410",  # positive regulation of MAPK cascade
  "GO:0010720",  # positive regulation of cell development
  "GO:0048589",  # developmental growth
  "GO:0060485",  # mesenchyme development
  "GO:1901890"   # positive regulation of cell junction assembly
)


wanted_ID_gse <- c(
  "GO:0010498",  # proteasomal protein catabolic process
  "GO:0043161",  # proteasome-mediated ubiquitin-dependent protein catabolic process
  "GO:0006511",  # ubiquitin-dependent protein catabolic process
  "GO:0051603",  # proteolysis involved in protein catabolic process
  "GO:0031329",  # regulation of cellular catabolic process
  "GO:0097190"  # apoptotic signaling pathway
)


# 2) gseGO 결과 객체 복사 및 @result 필터링
gse_sel <- gse
gse_sel@result <- gse@result[gse@result$ID %in% wanted_ID_gse, ]



View(gse_sel@result)
# 3) dotplot 시각화

dotplot(gse_sel, showCategory=5, split=".sign", font.size = 10) + facet_grid(.~.sign) +
  labs(title = "GSEA GO (Muscle Homeostasis)") +
  guides(color = guide_colorbar(title = "p.adjust"))



library(enrichplot)
library(ggplot2)
library(scales)






ridgeplot(gse_sel, showCategory = 7, fill = "NES") +
  scale_fill_gradient2(
    low      = "blue3",    # NES가 낮을 때 색
    mid      = "white",   # NES = 0 일 때 색
    high     = "#bb0c00", # NES가 높을 때 색
  ) +  labs(x = "enrichment distribution") +
  theme_minimal(base_size = 10)+
  theme(
    panel.grid.major = element_blank(),  # 주 격자선 제거
    panel.grid.minor = element_blank()   # 보조 격자선 제거
  )





library(enrichplot)

# 1) 그리고 싶은 pathway ID (예: proteasomal protein catabolic process)
my_pathway <- "GO:0097190"

# 2) gseResult 객체에서 해당 pathway의 숫자 인덱스 찾기
idx <- which(gse_sel@result$ID == my_pathway)

# 3) NES 값을 미리 뽑아서 변수에 저장
NES_val <- round(gse_sel@result$NES[idx], 2) #소수점 둘째까지만
NES_val

p.adjust <- round(gse_sel@result$p.adjust[idx])

p.adjust

# 4) enrichment plot 그리기
gp <- gseaplot2(
  gse_sel,
  geneSetID = idx,
  title     = paste0(
    gse_sel@result$Description[idx], 
    " (", my_pathway, ")"
  ),
  color     = "green"       # curve 색
)

gp





