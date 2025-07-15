
# Load required libraries
library(UpSetR)
library(readxl)
library(ComplexUpset)
library(ComplexHeatmap)


# Read DEG gene symbols from each Excel file 
GSE1 <- c(read_excel("/My path to/sample1.xlsx", sheet = "Upregulated")[[7]],
              read_excel("/My path to/sample1.xlsx", sheet = "Downregulated")[[7]])

GSE2 <- c(read_excel("/My path to/sample2.xlsx", sheet = "Upregulated")[[7]],
               read_excel("/My path to/sample2.xlsx", sheet = "Downregulated")[[7]])

GSE3 <- c(read_excel("/My path to/sample3.xlsx", sheet = "Upregulated")[[7]],
               read_excel("/My path to/sample3.xlsx", sheet = "Downregulated")[[7]])

GSE4 <- c(read_excel("/My path to/sample4.xlsx", sheet = "Upregulated")[[7]],
               read_excel("/My path to/sample4.xlsx", sheet = "Downregulated")[[7]])

GSE5 <- c(read_excel("/My path to/sample5.xlsx", sheet = "Upregulated")[[7]],
               read_excel("/My path to/sample5.xlsx", sheet = "Downregulated")[[7]])


deg_list <- list(
  Sample1 = GSE1,
  Sample2 = GSE2,
  Sample3 = GSE3,
  Sample4 = GSE4,
  Sample5 = GSE5
)



# Get all unique gene symbols across datasets
all_genes <- unique(unlist(deg_list))

# Fill binary 1/0 indicating whether a gene appears in each dataset
deg_matrix <- data.frame(Gene = all_genes)

for (name in names(deg_list)) {
  deg_matrix[[name]] <- as.integer(deg_matrix$Gene %in% deg_list[[name]])
}

# Set row names as gene names, then remove gene column
rownames(deg_matrix) <- deg_matrix$Gene
deg_matrix$Gene <- NULL



#----------------------------------------------------------------------

# Create intersection matrix using ComplexHeatmap

m = make_comb_mat(deg_matrix, mode = 'intersect')

m2 <- m[comb_degree(m) %in% c(1, 2, 3, 4, 5)]
ss2 = set_size(m2)
cs2 = comb_size(m2)

comb_colors <- ifelse(comb_degree(m2) == 5, "black", "black") 





#----------------------------------------------------------------------

# Create the UpSet plot with customized annotations
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



