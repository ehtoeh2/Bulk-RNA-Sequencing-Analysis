# Load required libraries
library(ggplot2)
library(readxl)
library(dplyr)

# Define a list of DEG Excel files (each with Up/Downregulated sheets)
file_list <- c(
  "/My path to/Sample1.xlsx", 
  "/My path to/Sample2.xlsx", 
  "/My path to/Sample3.xlsx", 
  "/My path to/Sample4.xlsx", 
  "/My path to/Sample5.xlsx"
)


deg_data <- data.frame()


# Read and combine data from all files
for (file in file_list) {
  condition <- tools::file_path_sans_ext(basename(file))  # 파일명을 condition으로 사용
  
  # Read upregulated genes
  up_data <- read_excel(file, sheet = "Upregulated") %>%
    mutate(condition = condition, direction = "Up")
  
  # Read downregulated genes
  down_data <- read_excel(file, sheet = "Downregulated") %>%
    mutate(condition = condition, direction = "Down")
  
  # Combine both up and down into one table
  deg_data <- bind_rows(deg_data, up_data, down_data)
}



# Calculate signed -log10(padj) values for Miami plot
deg_data <- deg_data %>%
  mutate(negLogP = ifelse(direction == "Up", -log10(padj), log10(padj)))


# Create Miami plot (fold-change based)
ggplot(deg_data, aes(x = condition, y = log2FoldChange, color = direction)) +
  geom_vline(xintercept = 0.5, linetype = "solid", color = "black", linewidth = 0.25) +  
  geom_jitter(alpha = 1.0, width = 0.25, size = 0.2) +  
  scale_color_manual(values = c("Up" = "red3",  "Down" = "blue4")) +
  labs(x = "Dataset", y = "log2FoldChange", title = "Miami Plot for Multiple DEG Analyses") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank())   



