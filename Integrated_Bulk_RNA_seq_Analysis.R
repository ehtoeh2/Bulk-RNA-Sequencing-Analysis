
library(dplyr)
library(ggrepel)
library(tximport)
library(DESeq2)
library(pheatmap)
library(openxlsx)
library(sva)
library(enrichplot)
set.seed(1234)  

res_df1 <- readRDS(file = "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/res_df1.rds")



TXNAME <- read.delim('/Users/daehwankim/Desktop/sequencing data/Bulk_seq_practice/TXNAME',header = F)
SYMBOL <-read.delim('/Users/daehwankim/Desktop/sequencing data/Bulk_seq_practice/SYMBOL',header = F)
tx2gene <- data.frame(TXNAME, SYMBOL)
colnames(tx2gene) <- c("TXNAME", "SYMBOL")

sample1 <- c('JCI(Joshua)_Cont1.tsv','JCI(Joshua)_Cont2.tsv','JCI(Joshua)_Cont3.tsv','JCI(Joshua)_Cont4.tsv',
             "JCI(Joshua)_Sham1.tsv","JCI(Joshua)_Sham2.tsv","JCI(Joshua)_Sham3.tsv","JCI(Joshua)_Sham4.tsv",
             "JNCI(Tseng)_Cont1.tsv","JNCI(Tseng)_Cont2.tsv","JNCI(Tseng)_Cont3.tsv",
             "JEM(Rupert)_Sham1.tsv", "JEM(Rupert)_Sham2.tsv", "JEM(Rupert)_Sham3.tsv",
             "Embo(sophia)_Cont1.tsv","Embo(sophia)_Cont2.tsv","Embo(sophia)_Cont3.tsv","Embo(sophia)_Cont4.tsv","Embo(sophia)_Cont5.tsv","Embo(sophia)_Cont6.tsv",
             
             "JCI(Joshua)_SC1.tsv","JCI(Joshua)_SC2.tsv","JCI(Joshua)_SC3.tsv","JCI(Joshua)_SC4.tsv",
             "JCI(Joshua)_SPC1.tsv","JCI(Joshua)_SPC2.tsv","JCI(Joshua)_SPC3.tsv","JCI(Joshua)_SPC4.tsv",
             "JNCI(Tseng)_C261.tsv","JNCI(Tseng)_C262.tsv","JNCI(Tseng)_C263.tsv",
             "JEM(Rupert)_KPC1.tsv","JEM(Rupert)_KPC2.tsv","JEM(Rupert)_KPC3.tsv","JEM(Rupert)_KPC4.tsv",
             "Embo(sophia)_CX1.tsv","Embo(sophia)_CX2.tsv","Embo(sophia)_CX3.tsv", "Embo(sophia)_CX4.tsv")
             

sample1 <- as.data.frame(sample1)


files1 <- file.path('/Users/daehwankim/Desktop/sequencing data/Integrated data', sample1$sample1)
names(files1) <- c('JCI(Joshua)_Cont1.tsv','JCI(Joshua)_Cont2.tsv','JCI(Joshua)_Cont3.tsv','JCI(Joshua)_Cont4.tsv',
                   "JCI(Joshua)_Sham1.tsv","JCI(Joshua)_Sham2.tsv","JCI(Joshua)_Sham3.tsv","JCI(Joshua)_Sham4.tsv",
                   "JNCI(Tseng)_Cont1.tsv","JNCI(Tseng)_Cont2.tsv","JNCI(Tseng)_Cont3.tsv",
                   "JEM(Rupert)_Sham1.tsv", "JEM(Rupert)_Sham2.tsv", "JEM(Rupert)_Sham3.tsv",
                   "Embo(sophia)_Cont1.tsv","Embo(sophia)_Cont2.tsv","Embo(sophia)_Cont3.tsv","Embo(sophia)_Cont4.tsv","Embo(sophia)_Cont5.tsv","Embo(sophia)_Cont6.tsv",
                   
                   "JCI(Joshua)_SC1.tsv","JCI(Joshua)_SC2.tsv","JCI(Joshua)_SC3.tsv","JCI(Joshua)_SC4.tsv",
                   "JCI(Joshua)_SPC1.tsv","JCI(Joshua)_SPC2.tsv","JCI(Joshua)_SPC3.tsv","JCI(Joshua)_SPC4.tsv",
                   "JNCI(Tseng)_C261.tsv","JNCI(Tseng)_C262.tsv","JNCI(Tseng)_C263.tsv",
                   "JEM(Rupert)_KPC1.tsv","JEM(Rupert)_KPC2.tsv","JEM(Rupert)_KPC3.tsv","JEM(Rupert)_KPC4.tsv",
                   "Embo(sophia)_CX1.tsv","Embo(sophia)_CX2.tsv","Embo(sophia)_CX3.tsv", "Embo(sophia)_CX4.tsv")


txi.kallisto.tsv1 <- tximport(files1, type = 'kallisto', tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion = TRUE)

sampleTable1 <- data.frame(condition=factor(rep(c("Control", "Cachexia"), times = c(20, 19))))

batch <- c(
  #Control
  rep("JCI_Joshua", 4),     
  rep("JCI_Joshua2", 4),
  rep("JNCI_Tseng", 3),
  rep("JEM_Rupert", 3),
  rep("Embo_sophia", 6),
  
  #Cachexia
rep("JCI_Joshua", 4),     
rep("JCI_Joshua2", 4),
rep("JNCI_Tseng", 3),
rep("JEM_Rupert", 4),
rep("Embo_sophia", 4))

length(batch)
sampleTable1$batch <- factor(batch)
sampleTable1 

rownames(sampleTable1) <- colnames(txi.kallisto.tsv1$counts)
dds1 <- DESeqDataSetFromTximport(txi.kallisto.tsv1, sampleTable1, ~condition)

dds1

#1
vst_dds <- vst(dds1, blind = TRUE)

View(vst_dds@assays@data)

# 4. 교정 전 PCA 확인 (옵션)
pca_dat_raw <- plotPCA(vst_dds, intgroup = c("batch","condition"), returnData=TRUE)
ggplot(pca_dat_raw, aes(PC1, PC2, color=condition, shape=condition)) +
  geom_point(size=3) +
  ggtitle("Before ComBat")

# 5. ComBat 으로 배치 교정
mat       <- assay(vst_dds)                   # 행렬 추출
batch_vec <- sampleTable1$batch               # 배치 정보
mod       <- model.matrix(~ condition, sampleTable1)  # condition 보존 모델

mat_combat <- ComBat(dat   = mat,
                     batch = batch_vec,
                     mod   = mod)

# 6. 교정된 값으로 덮어쓰기
assay(vst_dds) <- mat_combat

assay(vst_dds)


# 7. 교정 후 PCA 확인
pca_dat_bc <- plotPCA(vst_dds, intgroup = c("batch","condition"), returnData=TRUE)
ggplot(pca_dat_bc, aes(PC1, PC2, color=condition, shape=batch)) +
  geom_point(size=3) +
  ggtitle("After ComBat")



#PCA plot
pca_dat_bc <- plotPCA(vst_dds, intgroup = c("batch","condition"), returnData=TRUE) #returnData를 ture로 설정하면 그림이 아닌 데이터프레임으로 결과를 받음
pca_dat_bc

percentVar <- round(100 * attr(pca_dat_bc, "percentVar")) #분산 기여도를 %로 나타내서 보여줌
percentVar


library(ggplot2)

ggplot(pca_dat_bc, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(aes(shape = batch), size = 2) +
  scale_shape_manual(values = c(15, 16, 17, 0, 1, 2, 3, 6)) +
  scale_color_manual(values = c("Control" = "#077297", "Cachexia" = "#B24745")) +  # condition 별 색
  
  # 95% 신뢰 타원(ellipse) 추가
  stat_ellipse(aes(group = condition), 
               type  = "norm",    # 모양: 정규분포 기반 타원
               level = 0.95,      # 95% 신뢰구간
               linetype = "dashed", 
               size = 0.5) +
  
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("After ComBat with Ellipses") +
  xlim(-40, 70) +
  ylim(-40, 40) +
  
  theme(
    panel.background    = element_blank(),
    panel.grid.major    = element_blank(),
    panel.grid.minor    = element_blank(),
    axis.line           = element_line(color = "black")
  )






#DESeq2돌리기 위해 다시 batch effect 줄이기 = 이거는 DEseq2를 위한 보정이고, 아까 combat은 시각화용 데이터임

dds_batch <- DESeqDataSetFromTximport(txi.kallisto.tsv1, sampleTable1, 
                                      design = ~ batch + condition)



# (필요시) low count 필터링
keep <- rowSums(counts(dds_batch) >= 10) >= 2
dds_batch  <- dds_batch[keep,]






#1
deseq2.dds1 <- DESeq(dds_batch)										
deseq2.res1 <- results(deseq2.dds1)			
deseq2.res1 <- deseq2.res1[order(rownames(deseq2.res1)),]



# Select top 50 differentially expressed genes

res1 <- deseq2.res1
res1 <- na.omit(res1) #결측값 제거
res_ordered1 <- res1[order(res1$padj),] #padj값 오름차순으로 정렬
top_genes1 <- row.names(res_ordered1)[1:50] 


# Extract counts and normalize
counts1 <- counts(deseq2.dds1, normalized = T)
counts_top1 <- counts1[top_genes1,]

log_counts_top1 <- log2(counts_top1 + 1)

# Generate heatmap

#1
annotation_col <- data.frame(
  Condition = c(rep("Control", 20), rep("Cachexia", 19))
)

rownames(annotation_col)
colnames(log_counts_top1)

# 샘플 이름과 맞추기
rownames(annotation_col) <- colnames(log_counts_top1)

# Condition별 색상 지정

colors()

ann_colors <- list(
  Condition = c(Control = "cyan2", Cachexia = "darkred")  # 원하는 색상으로 변경 가능
)


ann_colors <- list(
  Condition = c(Control = "#077297", Cachexia = "#B24745")  # 원하는 색상으로 변경 가능
)

heatmap_colors <- colorRampPalette(c("blue3", "white", "red3"))(100)
heatmap_colors <- colorRampPalette(c("#4F7090", "white", "#FF7171"))(100)


pastel_colors <- colorRampPalette(c("#FFE6B3", "#FFFFFF", "#B3FFCC"))(100)

cyborg_colors <- colorRampPalette(c("#00FFFF", "#000000", "#FF00FF"))(100)

heatmap_colors <- colorRampPalette(c("blue", "green", "yellow"))(100)


cyborg_light <- colorRampPalette(
  c("#66FFFF",  # 옅은 시안
    "#FFFFFF",  # 화이트
    "#FF66FF"   # 옅은 마젠타
  )
)(100)

volcanoplot_color <- colorRampPalette(c("#173379", "white", "#bb0c00"))(100)

colo <- viridis::viridis(101)


# 🔹 관심 있는 유전자 리스트
genes_of_interest <- c("Tgfb1", "Col1a1", "Mmp2", "Apod","Mmp3")

# 🔹 해당 유전자들의 count 데이터 추출 (정규화된 값에서)
counts_subset <- counts1[rownames(counts1) %in% genes_of_interest, ]


# Heatmap 그리기
pheatmap(log_counts_top1, 
         scale = "row", 
         angle_col = 315,
         annotation_col = annotation_col,  # 그룹화된 컬러 표시
         annotation_colors = ann_colors,  # 그룹색상 지정
         cluster_cols = FALSE,  # 샘플 순서 유지
         show_colnames = FALSE, # 샘플 이름 숨기기
         color = heatmap_colors)  






#Volcanoplot 그리기
volcanoplot_color <- colorRampPalette(c("#00AFBB", "grey", "#bb0c00"))(100)

## Volcano plot
# Prepare the data for plotting
res_df1 <- as.data.frame(deseq2.res1)
res_df1 <- na.omit(res_df1)
res_df1$gene <- row.names(res_df1)
res_df1$gene_symbol <- row.names(res_df1)
head(res_df1)


# Create volcano plot
head(res_df1 %>% filter(gene_symbol == "Apod"))
head(res_df1 %>% filter(gene_symbol == "Trim63"))
head(res_df1 %>% filter(gene_symbol == "Fbxo32"))
head(res_df1 %>% filter(gene_symbol == "Apln"))

write.xlsx(upgene, "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/upgene2.xlsx", rownames = TRUE)
write.xlsx(downgene, "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/downgene.xlsx", rownames = TRUE)

upgenes <- c(
  "Gpx3",  "Scgb3a1", "Inhbb",  "Apod",   "Il4ra",
  "Hmgb2", "Adamts4","Ctsl",   "Il6ra",  "Cxcl13",
  "Lrg1",  "Vegfd",  "Lcn2",   "Serpina3n","Lox",
  "Clec18a","Itih4", "Chi3l1", "Adamts9","Wfdc17","Gapdh"
)

upgene <- tibble(
  res_df1 %>%
    filter(gene_symbol %in% upgenes)
)

down_genes <- c(
  "Ighm",   "Grem2",  "Postn",   "Igfbp5",  "Mfap4",
  "Gdf11",  "Lingo3", "Col1a1",  "Kera",    "Col1a2",
  "Islr",   "Angptl2","Col6a2",  "Tfrc",    "Col6a1",
  "Sparcl1","Ssc5d", "Col3a1",  "Pthlh",   "Sema3a",
  "Metrn",  "Ecrg4",  "Wnt5a",   "Col11a2", "Tnc",
  "Cdh13",  "Loxl2",  "Fn1",     "Frzb",    "Col4a5",
  "Xylt1",  "Gkn3",   "Apln",    "Smpdl3b", "Col11a1",
  "Col4a3", "Cpe",    "Ostn",    "Ntf5",    "Wnt10b",
  "Fbln7",  "Col5a2", "Ache",    "Mmp15",   "Serpinb6b",
  "Col14a1"
)


downgene <- tibble(
  res_df1 %>%
    filter(gene_symbol %in% down_genes)
)



#- 곱해서 부호 바꿔주기 (이거는 상황에 따라서 맞게 사용해야함)
res_df1 <- res_df1 %>%
  mutate(log2FoldChange = -log2FoldChange) 



res_df1$gene <- "Non significant"
res_df1$gene[res_df1$log2FoldChange > 1 & res_df1$padj < 0.05] <- "Upregulated"
res_df1$gene[res_df1$log2FoldChange < -1 & res_df1$padj < 0.05] <- "Downregulated"
res_df1 <- res_df1[order(res_df1$padj),]


mycolors <- c("Downregulated" = "#00AFBB", "Upregulated" = "#bb0c00", "Non significant" = "grey")
mycolors2 <- c("Downregulated" = "#173379", "Upregulated" = "red3", "Non significant" = "grey")
scale_color_manual(values = mycolors, guide = "none")


res_df1$delabel <- NA
res_df1$delabel[res_df1$gene != "Non significant"] <- res_df1$gene_symbol[res_df1$gene != "Non significant"]
head(res_df1, 10)




annotation1 <- res_df1 %>% filter(gene_symbol %in% c("Eda2r", "Fbxo32", "Trim63", "Apod"))
head(annotation1, 10)

annotation2 <- res_df1 %>% filter(gene_symbol %in% c("Gpx3", "Inhbb","Apod","Hmgb2","Scgb3a1","Tnc", "Angptl2","Ostn","Apln","Gkn3","Eda2r", "Fbxo32","Trim63"))






volcano_plot1 <- ggplot(res_df1, aes(x = log2FoldChange, y = -log10(padj), col = gene, label = delabel)) + 
  geom_point(alpha = 0.6, size = 1) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_color_manual(values = mycolors2) + 
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right") + 
  #geom_vline(xintercept = c(-1, 1), col = "grey") + <-이거는 optional
  geom_hline(yintercept = -log10(0.09), col = "grey") + 
  
  geom_label_repel(                  #geom_label_repel을 사용하면 박스로 표시가 됨
    data = annotation2, 
    box.padding = 0.5,           # 텍스트 주변의 여백 조정
    segment.color = 'black', 
    segment.size = 0.4, 
    size = 3.0, 
    max.overlaps = 100,        # 겹침을 최소화하기 위해 overlaps 수 줄임
    key_glyph = draw_key_point, 
    force = 3,                   # 텍스트 밀어내기 힘 설정
    nudge_y = 2                  # y축으로 텍스트를 더 높이 이동
  ) + 
  ylim(c(0, 87)) +              # y축 데이터를 0부터 45까지 보여줌 
  xlim(c(-9, 9))                # X축 데이터를 -5부터 5까지 보여줌

print(volcano_plot1)

saveRDS(res_df1, file = "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/res_df1.rds")


# ORA analysis
library(clusterProfiler)
library(org.Mm.eg.db)  # 마우스일 경우
library(dplyr)
library(ggplot2)



# upregulated
deg_up <- res_df1 %>%
  filter(log2FoldChange > 1 & padj < 0.05) %>%
  pull(gene_symbol)

# downregulated
deg_down <- res_df1 %>%
  filter(log2FoldChange < -1 & padj < 0.05) %>%
  pull(gene_symbol)



# Upregulated ORA
go_up <- enrichGO(gene = deg_up,
                  OrgDb = org.Mm.eg.db,
                  keyType = "SYMBOL",
                  ont = "BP",  # 또는 MF, CC
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2)

# FDR (q-value, 즉 p.adjust) 기준으로 유의한 결과만 추출
sig_go_up <- go_up@result %>% 
  filter(p.adjust < 0.05)



# Downregulated ORA
go_down <- enrichGO(gene = deg_down,
                    OrgDb = org.Mm.eg.db,
                    keyType = "SYMBOL",
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2)

# FDR (q-value, 즉 p.adjust) 기준으로 유의한 결과만 추출
sig_go_down <- go_down@result %>% 
  filter(p.adjust < 0.05)


View(sig_go_up)
View(sig_go_down)

write.xlsx(sig_go_up, "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/Go_up.xlsx", rownames = TRUE)
write.xlsx(sig_go_down, "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/Go_down.xlsx", rownames = TRUE)
write.xlsx(top50, "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/top50.xlsx", rownames = TRUE)
top50 <- data.frame(head(res_df1$gene_symbol, 50)) 

#전체 BP bar plot
dotplot(go_up, showCategory = 15, font.size = 12, title = "GO BP Enrichment (Upregulated)") +
  theme_minimal()

dotplot(go_down, showCategory = 15, font.size = 12, title = "GO BP Enrichment (downregulated)") +
  theme_minimal()


barplot(go_up, showCategory = 15, font.size = 12, title = "GO Analysis (Up)")
barplot(go_down, showCategory = 15, font.size = 12, title = "GO Analysis (Down)")




# 1) 원하는 ID으로 subset 

# 1) 원하는 15개 GO term ID (up)
wanted_ID_up <- c(
  "GO:0014889", # muscle atrophy
  "GO:0014741", # negative regulation of muscle hypertrophy
  "GO:0014745", # negative regulation of muscle adaptation
  "GO:0010508", # positive regulation of autophagy
  "GO:0016236", # macroautophagy
  "GO:0010657", # muscle cell apoptotic process
  "GO:0010661", # positive regulation of muscle cell apoptotic process
  "GO:0010663", # positive regulation of striated muscle cell apoptotic process
  "GO:0097193", # intrinsic apoptotic signaling pathway
  "GO:0097191", # extrinsic apoptotic signaling pathway
  "GO:0006919", # activation of cysteine-type endopeptidase activity involved in apoptotic process
  "GO:0030968", # endoplasmic reticulum unfolded protein response
  "GO:0035966", # response to topologically incorrect protein
  "GO:0031329", # regulation of cellular catabolic process
  "GO:0031396"  # regulation of protein ubiquitination
)

# 1) 원본 복사
go_up_sel <- go_up

# 2) @result 데이터프레임에서 ID 필터링
go_up_sel@result <- go_up@result[
  go_up@result$ID %in% wanted_ID_up,
]


# 3) dotplot 그리기
dotplot(
  go_up_sel,
  showCategory = length(wanted_ID_up),
  font.size    = 12,
  title        = "GO BP (Muscle Atrophy)"
) + theme_minimal()


#4 Barplot 그리기 
barplot(go_up_sel, showCategory = 7, font.size = 10 ,title = "GO BP UP (Muscle Atrophy)") + 
  scale_fill_gradient(low = "red3", high = "#FADBD8")  # 값이 작을수록 진한 파랑
help(barplot)

help(scale_fill_gradient)




# 1) 원하는 15개 GO term ID (down)
wanted_ID_down <- c(
  "GO:0007517",
  "GO:0055001",
  "GO:0055002",
  "GO:0033002",
  "GO:0051146",
  "GO:0048741",
  "GO:0014902",
  "GO:0030239",
  "GO:0045214",
  "GO:0006936",
  "GO:0003012",
  "GO:0006942",
  "GO:0042692",
  "GO:0060538",
  "GO:0007519"
)

# 1) 원본 복사
go_down_sel <- go_down

# 2) @result 데이터프레임에서 ID 필터링
go_down_sel@result <- go_down@result[
  go_down@result$ID %in% wanted_ID_down,
]


# 3) dotplot 그리기
dotplot(
  go_down_sel,
  showCategory = length(wanted_ID_down),
  font.size    = 12,
  title        = "GO BP (Muscle Atrophy)"
) + theme_minimal()


#4 Barplot 그리기 
barplot(go_down_sel, showCategory = 10, font.size = 10, title = "GO BP DOWN (Muscle Maintenance)") + 
  scale_fill_gradient(low = "blue3", high = "#D6EAF8")  # 값이 작을수록 진한 파랑



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




