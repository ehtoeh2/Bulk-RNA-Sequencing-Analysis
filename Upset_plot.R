

library(UpSetR)
library(readxl)
library(ComplexUpset)
library(ComplexHeatmap)


# 엑셀 파일 불러오기
GSE65936 <- c(read_excel("/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Bulk_seq_analysis_practice/JNCI_bulk (Tseng)/DEG_results(Tseng).xlsx", sheet = "Upregulated")[[7]],
              read_excel("/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Bulk_seq_analysis_practice/JNCI_bulk (Tseng)/DEG_results(Tseng).xlsx", sheet = "Downregulated")[[7]])

GSE123310 <- c(read_excel("/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Bulk_seq_analysis_practice/JEM_bulk(Rupert)/DEG_results(Rupert).xlsx", sheet = "Upregulated")[[7]],
               read_excel("/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Bulk_seq_analysis_practice/JEM_bulk(Rupert)/DEG_results(Rupert).xlsx", sheet = "Downregulated")[[7]])

GSE142455_SC = c(read_excel("/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Bulk_seq_analysis_practice/250304_JCI(joshua)/DEG_results(joshua_SC).xlsx", sheet = "Upregulated")[[7]], 
                 read_excel("/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Bulk_seq_analysis_practice/250304_JCI(joshua)/DEG_results(joshua_SC).xlsx", sheet = "Downregulated")[[7]])

GSE142455_Orth = c(read_excel("/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Bulk_seq_analysis_practice/250304_JCI(joshua)/DEG_results(joshua).xlsx", sheet = "Upregulated")[[7]], 
                   read_excel("/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Bulk_seq_analysis_practice/250304_JCI(joshua)/DEG_results(joshua).xlsx", sheet = "Downregulated")[[7]])

GSE138464 = c(read_excel("/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Bulk_seq_analysis_practice/250303_EMBO(sophia)/DEG_results(EMBO).xlsx", sheet = "Upregulated")[[7]], 
              read_excel("/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Bulk_seq_analysis_practice/250303_EMBO(sophia)/DEG_results(EMBO).xlsx", sheet = "Downregulated")[[7]])



deg_list <- list(
  GSE65936 = GSE65936,
  GSE123310 = GSE123310,
  GSE142455_SC = GSE142455_SC,
  GSE142455_Orth = GSE142455_Orth,
  GSE138464 = GSE138464
)



# 모든 유전자 리스트 통합
all_genes <- unique(unlist(deg_list))

# 바이너리 매트릭스 생성
deg_matrix <- data.frame(Gene = all_genes)

for (name in names(deg_list)) {
  deg_matrix[[name]] <- as.integer(deg_matrix$Gene %in% deg_list[[name]])
}

# Gene 열을 rownames으로 변경
rownames(deg_matrix) <- deg_matrix$Gene
deg_matrix$Gene <- NULL



m = make_comb_mat(deg_matrix, mode = 'intersect')

m2 <- m[comb_degree(m) %in% c(1, 2, 3, 4, 5)]
ss2 = set_size(m2)
cs2 = comb_size(m2)

ss2
cs2
comb_colors <- ifelse(comb_degree(m2) == 5, "black", "black") 

ht = UpSet(m2, 
           pt_size = unit(5, "mm"),  # 점 크기 키우기
           lwd = 4,
           comb_col = comb_colors,
           set_order = order(ss2),
           comb_order = order(comb_degree(m2), -cs2),
           top_annotation = HeatmapAnnotation(
             "DEG Intersections" = anno_barplot(cs2, 
                                                ylim = c(0, max(cs2)*1.1),
                                                border = FALSE, 
                                                gp = gpar(fill = "black"), 
                                                height = unit(8, "cm")
             ), 
             annotation_name_side = "left", 
             annotation_name_rot = 90),
           left_annotation = rowAnnotation(
             "DEGs Per data set " = anno_barplot(-ss2, 
                                                 baseline = 0,
                                                 bar_width = 0.5,
                                                 axis_param = list(
                                                   at = c(0, -2000, -4000, -6000),
                                                   labels = c(0, 2000, 4000, 6000),
                                                   labels_rot = 0),
                                                 border = FALSE, 
                                                 gp = gpar(fill = "black"),
                                                 labels_gp = gpar(fontsize = 10),
                                                 width = unit(7, "cm")
             ),
             set_name = anno_text(set_name(m2), 
                                  location = 0.5, 
                                  just = "center",
                                  width = max_text_width(set_name(m2)) + unit(4, "mm"))
           ), 
           right_annotation = NULL,
           show_row_names = FALSE)
ht = draw(ht)
od = column_order(ht)

decorate_annotation("DEG Intersections", {
  grid.text(cs2[od], x = seq_along(cs2), y = unit(cs2[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 6, col = "#404040"), rot = 45)
})














