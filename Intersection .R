# Load required libraries
library(readxl)
library(writexl)
library(dplyr)    
library(openxlsx)  
library(org.Mm.eg.db) 
library(AnnotationDbi)

# Read Excel files containing upregulated and downregulated gene lists
df_up <- read_excel("/My path to/up_file.xlsx")
df_down <- read_excel("/My path to/down_file.xlsx")


# Extract common genes across all samples (columns)
common_values_up <- Reduce(intersect, as.list(df_up[, c(
  "smaple1", 
  "smaple2", 
  "smaple3", 
  "smaple4", 
  "smaple5"
)]))

print(common_values_up)


common_values_down <- Reduce(intersect, as.list(df_down[, c(
  "smaple1", 
  "smaple2", 
  "smaple3", 
  "smaple4", 
  "smaple5"
)]))

print(common_values_down)


# Save both common genes lists in a single Excel file with separate sheets
df_common_up <- data.frame(Common_Values = common_values_up)  
df_common_down <- data.frame(Common_Values = common_values_down)  

sheets <- list(
  Up   = data.frame(Common_Values = common_values_up),
  Down = data.frame(Common_Values = common_values_down)
)

out_file <- "/My path to/common_value.xlsx"
write_xlsx(sheets, path = out_file)





# ----------- Extract genes associated with a specific GO term -----


go_id <- "GO:0000000"
eg_in_go <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys    = go_id,
  keytype = "GO",         
  columns = c("ENTREZID", "SYMBOL")
)

# Remove duplicates and select only necessary columns
eg_in_go <- unique(eg_in_go[, c("ENTREZID", "SYMBOL")])
cat("Extracted", nrow(eg_in_go), "mouse genes associated with GO term\n")
View(eg_in_go)

# Save GO term genes to Excel
output_file <- "/My path to/GO_genes.xlsx"
write.xlsx(
  eg_in_go,
  file      = output_file,
  sheetName = "GO term",
  rowNames  = FALSE
)


#Load common gene lists and GO genes to intersect
common_file <- "/My path to/common_value.xlsx"
Go_genes_file <- "My path to/GO_genes.xlsx"

# Load up and downregulated common gene lists
up   <- read_excel(common_file, sheet = "Up")   %>% pull(1) %>% unique()
down <- read_excel(common_file, sheet = "Down") %>% pull(1) %>% unique()


# Load GO term gene list and extract gene symbols
protein_secretion <- read_excel(Protein_secretion_file, sheet = 1) %>% 
  rename_with(tolower) %>%              
  pull(symbol) %>%                      
  unique()

#Find intersection between up/downregulated genes and GO term genes
up_in_protein_secretion   <- intersect(up, protein_secretion)
down_in_protein_secretion <- intersect(down, protein_secretion)


cat("Upregulated ∩ Extracellular space: ", length(up_in_protein_secretion), " genes\n")
cat("Downregulated ∩ Extracellular space: ", length(down_in_protein_secretion), " genes\n")

#Combine intersection results into a single table and save to Excel

# find the longest length
n <- max(length(up_in_extracellular_space), length(down_in_extracellular_space))

# pad each vector to length n
up_padded   <- c(up_in_extracellular_space,   rep(NA, n - length(up_in_extracellular_space)))
down_padded <- c(down_in_extracellular_space, rep(NA, n - length(down_in_extracellular_space)))


# Combine into a data frame
out_df <- data.frame(
  up_in_extracellular_space   = up_padded,
  down_in_extracellular_space = down_padded
)

# Save intersection table to Excel
write.xlsx(out_df,
           file      = "/My path to/common_vs_Goterm_intersection.xlsx",
           sheetName = "intersection",
           rowNames  = FALSE)


