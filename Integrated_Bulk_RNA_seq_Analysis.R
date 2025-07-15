# Road required library
library(dplyr)
library(ggrepel)
library(tximport)
library(DESeq2)
library(pheatmap)
library(openxlsx)
library(sva)
library(enrichplot)
library(ggplot2)

#Set the seed for reproducibility 
set.seed(1234)  


# Read in transcript-to-gene mapping files 
TXNAME <- read.delim('/Users/daehwankim/Desktop/sequencing data/Bulk_seq_practice/TXNAME',header = F)
SYMBOL <-read.delim('/Users/daehwankim/Desktop/sequencing data/Bulk_seq_practice/SYMBOL',header = F)
tx2gene <- data.frame(TXNAME, SYMBOL)
colnames(tx2gene) <- c("TXNAME", "SYMBOL")


sample1 <- c(
  "Sample1_Cont1.tsv", "Sample1_Cont2.tsv", "Sample1_Cont3.tsv", "Sample1_Cont4.tsv",
  "Sample1_Sham1.tsv", "Sample1_Sham2.tsv", "Sample1_Sham3.tsv", "Sample1_Sham4.tsv",
  "Sample2_Cont1.tsv", "Sample2_Cont2.tsv", "Sample2_Cont3.tsv",
  "Sample3_Sham1.tsv", "Sample3_Sham2.tsv", "Sample3_Sham3.tsv",
  "Sample4_Cont1.tsv", "Sample4_Cont2.tsv", "Sample4_Cont3.tsv", "Sample4_Cont4.tsv", "Sample4_Cont5.tsv", "Sample4_Cont6.tsv",
  
  "Sample1_SC1.tsv", "Sample1_SC2.tsv", "Sample1_SC3.tsv", "Sample1_SC4.tsv",
  "Sample1_SPC1.tsv", "Sample1_SPC2.tsv", "Sample1_SPC3.tsv", "Sample1_SPC4.tsv",
  "Sample2_C261.tsv", "Sample2_C262.tsv", "Sample2_C263.tsv",
  "Sample3_KPC1.tsv", "Sample3_KPC2.tsv", "Sample3_KPC3.tsv", "Sample3_KPC4.tsv",
  "Sample4_CX1.tsv", "Sample4_CX2.tsv", "Sample4_CX3.tsv", "Sample4_CX4.tsv"
)

             

sample1 <- as.data.frame(sample1)


files1 <- file.path('/My path to/Integrated data', sample1$sample1)
names(files1) <- c(
  "Sample1_Cont1.tsv", "Sample1_Cont2.tsv", "Sample1_Cont3.tsv", "Sample1_Cont4.tsv",
  "Sample1_Sham1.tsv", "Sample1_Sham2.tsv", "Sample1_Sham3.tsv", "Sample1_Sham4.tsv",
  "Sample2_Cont1.tsv", "Sample2_Cont2.tsv", "Sample2_Cont3.tsv",
  "Sample3_Sham1.tsv", "Sample3_Sham2.tsv", "Sample3_Sham3.tsv",
  "Sample4_Cont1.tsv", "Sample4_Cont2.tsv", "Sample4_Cont3.tsv", "Sample4_Cont4.tsv", "Sample4_Cont5.tsv", "Sample4_Cont6.tsv",
  
  "Sample1_SC1.tsv", "Sample1_SC2.tsv", "Sample1_SC3.tsv", "Sample1_SC4.tsv",
  "Sample1_SPC1.tsv", "Sample1_SPC2.tsv", "Sample1_SPC3.tsv", "Sample1_SPC4.tsv",
  "Sample2_C261.tsv", "Sample2_C262.tsv", "Sample2_C263.tsv",
  "Sample3_KPC1.tsv", "Sample3_KPC2.tsv", "Sample3_KPC3.tsv", "Sample3_KPC4.tsv",
  "Sample4_CX1.tsv", "Sample4_CX2.tsv", "Sample4_CX3.tsv", "Sample4_CX4.tsv"
)



# Import transcript-level quantification results from Kallisto and summarize to gene-level
txi.kallisto.tsv1 <- tximport(files1, type = 'kallisto', tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion = TRUE)

sampleTable1 <- data.frame(condition=factor(rep(c("Control", "Experimental"), times = c(20, 19))))


# Define Batch information for each sample
batch <- c(
  #Control
  rep("Sample1", 4),     
  rep("Sample2", 4),
  rep("Sample3", 3),
  rep("Sample4", 3),
  rep("Sample5", 6),
  
  #Experimental
  rep("Sample1", 4),     
  rep("Sample2", 4),
  rep("Sample3", 3),
  rep("Sample4", 4),
  rep("Sample5", 4))


sampleTable1$batch <- factor(batch)

# Set the row names of the sample table 
rownames(sampleTable1) <- colnames(txi.kallisto.tsv1$counts)


# Create a DESeq2 dataset object 
dds1 <- DESeqDataSetFromTximport(txi.kallisto.tsv1, sampleTable1, ~condition)


#-------------------------------PCA Plot --------------------------------------------------
# Apply variance-stabilizing transformation (VST)
vst_dds <- vst(dds1, blind = TRUE)


# PCA before batch correction (optional)
pca_dat_raw <- plotPCA(vst_dds, intgroup = c("batch","condition"), returnData=TRUE)
ggplot(pca_dat_raw, aes(PC1, PC2, color=condition, shape=condition)) +
  geom_point(size=3) +
  ggtitle("Before ComBat")

# Apply ComBat for batch correction (For visualization)
mat       <- assay(vst_dds)                   
batch_vec <- sampleTable1$batch               
mod       <- model.matrix(~ condition, sampleTable1)  

mat_combat <- ComBat(dat   = mat,
                     batch = batch_vec,
                     mod   = mod)

# Overwrite vst values with ComBat-corrected matrix
assay(vst_dds) <- mat_combat


# PCA after batch correction
pca_dat_bc <- plotPCA(vst_dds, intgroup = c("batch","condition"), returnData=TRUE)
ggplot(pca_dat_bc, aes(PC1, PC2, color=condition, shape=batch)) +
  geom_point(size=3) +
  ggtitle("After ComBat")



#Create PCA plot
pca_dat_bc <- plotPCA(vst_dds, intgroup = c("batch","condition"), returnData=TRUE) 
percentVar <- round(100 * attr(pca_dat_bc, "percentVar")) 


ggplot(pca_dat_bc, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(aes(shape = batch), size = 2) +
  scale_shape_manual(values = c(15, 16, 17, 0, 1, 2, 3, 6)) +
  scale_color_manual(values = c("Control" = "#077297", "Experimental" = "#B24745")) +  
  
  stat_ellipse(aes(group = condition), 
               type  = "norm",    
               level = 0.95,      
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



#-------------------------------Heatmap Plot --------------------------------------------------

# Rebuild DESeq2 object with batch effect 
# (This is for statistical testing, unlike ComBat which was used for visualization only)

dds_batch <- DESeqDataSetFromTximport(txi.kallisto.tsv1, sampleTable1, 
                                      design = ~ batch + condition)



#Run DESeq2 differential expression analysis
deseq2.dds1 <- DESeq(dds_batch)										
deseq2.res1 <- results(deseq2.dds1)			
deseq2.res1 <- deseq2.res1[order(rownames(deseq2.res1)),]


# Get top 50 most significant genes (lowest padj)
res1 <- deseq2.res1
res1 <- na.omit(res1) 
res_ordered1 <- res1[order(res1$padj),] 
top_genes1 <- row.names(res_ordered1)[1:50] 


# Extract counts and normalize
counts_norm <- counts(deseq2.dds1, normalized = TRUE)
counts_top <- counts_norm[top_genes,]

# Apply log2 transformation to stabilize variance
log_counts_top1 <- log2(counts_top1 + 1)



# Annotation for sample groups
annotation_col <- data.frame(
  Condition = c(rep("Control", 20), rep("Experimental", 19))
)

# Match row names of the annotation to column names of the count matrix
rownames(annotation_col) <- colnames(log_counts_top1)

# Define colors for each condition
ann_colors <- list(Condition = c(Control = "#077297", Experimental = "#B24745"))
heatmap_colors <- colorRampPalette(c("blue3", "white", "red3"))(100)


# Create a heatmap plot
pheatmap(log_counts_top1, 
         scale = "row", 
         angle_col = 315,
         annotation_col = annotation_col,  
         annotation_colors = ann_colors,  
         cluster_cols = FALSE,  
         show_colnames = FALSE, 
         color = heatmap_colors)  




#-------------------------------Volcano Plot --------------------------------------------------


# Prepare the data for plotting
res_df1 <- as.data.frame(deseq2.res1)
res_df1 <- na.omit(res_df1)
res_df1$gene <- row.names(res_df1)
res_df1$gene_symbol <- row.names(res_df1)


#Flip sign of log2FoldChange based on marker gene direction (Optional)
head(res_df1 %>% filter(gene_symbol == "Marker gene1"))
res_df1 <- res_df1 %>%
  mutate(log2FoldChange = -log2FoldChange) # <-- Use only if necessary


# Label genes by significance and fold change
res_df1$gene <- "Non significant"
res_df1$gene[res_df1$log2FoldChange > 1 & res_df1$padj < 0.05] <- "Upregulated"
res_df1$gene[res_df1$log2FoldChange < -1 & res_df1$padj < 0.05] <- "Downregulated"

res_df1 <- res_df1[order(res_df1$padj),]

res_df1$delabel <- NA
res_df1$delabel[res_df1$gene != "Non significant"] <- res_df1$gene_symbol[res_df1$gene != "Non significant"]



#Create Volcano plot
mycolors <- c("Downregulated" = "#173379", "Upregulated" = "red3", "Non significant" = "grey")
annotation1 <- res_df1 %>% filter(gene_symbol %in% c("Gene1", "Gene2", "Gene3", "Gene4"))


volcano_plot1 <- ggplot(res_df1, aes(x = log2FoldChange, y = -log10(padj), col = gene, label = delabel)) + 
  geom_point(alpha = 0.6, size = 1) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_color_manual(values = mycolors) + 
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right") + 
  #geom_vline(xintercept = c(-1, 1), col = "grey") <- logFC cutoff lines (Optional)
  geom_hline(yintercept = -log10(0.09), col = "grey") + 
  
  # Label selected genes using repelling boxes
  geom_label_repel(                  
    data = annotation1, 
    box.padding = 0.5,          
    segment.color = 'black', 
    segment.size = 0.4, 
    size = 3.0, 
    max.overlaps = 100,        
    key_glyph = draw_key_point, 
    force = 3,                   
    nudge_y = 2                  
  ) + 
  ylim(c(0, 87)) +              
  xlim(c(-9, 9))             

print(volcano_plot1)

saveRDS(res_df1, file = "/My path to/res_df1.rds")






