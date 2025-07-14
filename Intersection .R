library(readxl)
library(writexl)



# 엑셀 파일 불러오기 (파일명을 실제 파일명으로 변경)
df <- read_excel("/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/Bulk DEG 통합본.xlsx")
df_up <- read_excel("/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/BULK_upregulated.xlsx")
df_down <- read_excel("/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/BULK_downregulated.xlsx")


common_values_up <- Reduce(intersect, as.list(df_up[, c(
  "EMBO (Sophia)", 
  "JCI (joshua)...2", 
  "JCI (joshua)...3", 
  "JEM (Rupert)", 
  "JNCI (Tseng)"
)]))


print(common_values_up)


common_values_down <- Reduce(intersect, as.list(df_down[, c(
  "EMBO (Sophia)", 
  "JCI (joshua)...2", 
  "JCI (joshua)...3", 
  "JEM (Rupert)", 
  "JNCI (Tseng)"
)]))


print(common_values_down)

df_common_up <- data.frame(Common_Values = common_values_up)  # 데이터프레임 변환
write_xlsx(df_common_up, "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/common_values_up_wo_talbert.xlsx")  # 엑셀 파일로 저장

df_common_down <- data.frame(Common_Values = common_values_down)  # 데이터프레임 변환
write_xlsx(df_common_down, "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/common_values_down_wo_talbert.xlsx")  # 엑셀 파일로 저장


# 2) 리스트로 시트별 데이터 준비
sheets <- list(
  Up   = data.frame(Common_Values = common_values_up),
  Down = data.frame(Common_Values = common_values_down)
)

# 3) 한 파일에 쓰기
out_file <- "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/common_values_updown_wo_talbert.xlsx"
write_xlsx(sheets, path = out_file)





# 1) 패키지 로드
library(readxl)    # 엑셀 불러오기
library(dplyr)     # 데이터 처리
library(openxlsx)  # 엑셀 쓰기

# 2) 파일 경로 지정
common_file <- "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/common_values_updown_wo_talbert.xlsx"
Protein_secretion_file <- "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/mouse_protein_secretion_genes.xlsx"

# 3) Up / Down 리스트 불러오기
up   <- read_excel(common_file, sheet = "Up")   %>% pull(1) %>% unique()
down <- read_excel(common_file, sheet = "Down") %>% pull(1) %>% unique()

#    — 만약 시트명이 정확하지 않으면 sheet = 1, sheet = 2 처럼 숫자로 지정하셔도 됩니다.

# 4) Extracellular space 유전자 불러오기
protein_secretion <- read_excel(Protein_secretion_file, sheet = 1) %>% 
  rename_with(tolower) %>%               # 컬럼명 소문자 통일
  pull(symbol) %>%                      # SYMBOL 컬럼만 벡터로
  unique()

# 5) 교집합 구하기
up_in_protein_secretion   <- intersect(up, protein_secretion)
down_in_protein_secretion <- intersect(down, protein_secretion)

up_in_protein_secretion
down_in_protein_secretion

# 6) 결과 확인
cat("Upregulated ∩ Extracellular space: ", length(up_in_protein_secretion), " genes\n")
cat("Downregulated ∩ Extracellular space: ", length(down_in_protein_secretion), " genes\n")

# 7) 결과를 엑셀로 저장 (옵션)
out_df <- data.frame(
  up_in_protein_secretion   = c(up_in_protein_secretion,   rep(NA, max(0, length(up_in_protein_secretion) - length(up_in_protein_secretion)))),
  down_in_protein_secretion = c(down_in_protein_secretion, rep(NA, max(0, length(down_in_protein_secretion)   - length(down_in_protein_secretion))))
)

out_df

write.xlsx(out_df,
           file      = "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/common_vs_proteinsecretion_intersection.xlsx",
           sheetName = "intersection",
           rowNames  = FALSE)






# 1) 필요한 패키지 설치 (한 번만)


# 2) 패키지 로드
library(org.Mm.eg.db)       # 마우스 유전자 어노테이션
library(AnnotationDbi)
library(openxlsx)           # 엑셀 쓰기용

# 3) GO:0005615에 속한 Entrez Gene ID 및 기호(Symbol) 추출
go_id <- "GO:0005615"
eg_in_go <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys    = go_id,
  keytype = "GO",         # GO term을 키로 사용
  columns = c("ENTREZID", "SYMBOL")
)

# 4) 중복 제거 및 컬럼 정리
eg_in_go <- unique(eg_in_go[, c("ENTREZID", "SYMBOL")])
cat("총", nrow(eg_in_go), "개의 마우스 유전자 추출됨\n")

View(eg_in_go)

# 5) 엑셀 파일로 저장
output_file <- "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/mouse_extracellular_space_genes.xlsx"
write.xlsx(
  eg_in_go,
  file      = output_file,
  sheetName = "extracellular_space",
  rowNames  = FALSE
)




# 3) GO:0005615에 속한 Entrez Gene ID 및 기호(Symbol) 추출
go_id_protein_secretion <- "GO:0046903"
eg_in_go <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys    = go_id_protein_secretion,
  keytype = "GO",         # GO term을 키로 사용
  columns = c("ENTREZID", "SYMBOL")
)

eg_in_go

# 4) 중복 제거 및 컬럼 정리
eg_in_go <- unique(eg_in_go[, c("ENTREZID", "SYMBOL")])
cat("총", nrow(eg_in_go), "개의 마우스 유전자 추출됨\n")

View(eg_in_go)

# 5) 엑셀 파일로 저장
output_file <- "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/mouse_protein_secretion_genes.xlsx"
write.xlsx(
  eg_in_go,
  file      = output_file,
  sheetName = "Protein_Secretion",
  rowNames  = FALSE
)



#GO:0046903 – secretion
#GO:0009306 – protein secretion
#GO:0070254 – mucus secretion


# 필요한 패키지 로드
library(AnnotationDbi)
library(org.Mm.eg.db)
library(openxlsx)

# 필요한 패키지 로드
library(AnnotationDbi)
library(org.Mm.eg.db)
library(openxlsx)

# 1) 처리할 GO term 목록과 시트 이름
go_terms <- c(
  "GO:0005615", # extracellular space
  "GO:0005576", # extracellular region
  "GO:0044421", # extracellular region part
  "GO:0070062", # extracellular exosome
  "GO:0046903", # secretion
  "GO:0050658", # secretion by cell
  "GO:0009306", # protein secretion
  "GO:0018250"  # peptide secretion
)
sheet_names <- c(
  "extracellular_space",
  "extracellular_region",
  "extracellular_region_part",
  "extracellular_exosome",
  "secretion",
  "secretion_by_cell",
  "protein_secretion",
  "peptide_secretion"
)

# 2) org.Mm.eg.db 에서 유효한 GO ID 목록 가져오기
valid_go_ids <- keys(org.Mm.eg.db, keytype = "GO")

# 3) 결과를 저장할 폴더
output_dir <- "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/"

# 4) 반복 처리
for (i in seq_along(go_terms)) {
  go_id <- go_terms[i]
  sheet <- sheet_names[i]
  
  # GO ID 유효성 체크
  if (!go_id %in% valid_go_ids) {
    warning(sprintf("'%s'는 org.Mm.eg.db에 존재하지 않는 GO ID입니다. 건너뜁니다.", go_id))
    next
  }
  
  # EntrezID & SYMBOL 추출
  eg_df <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys    = go_id,
    keytype = "GO",
    columns = c("ENTREZID", "SYMBOL")
  )
  
  # 중복 제거 & 컬럼 정리
  eg_df <- unique(eg_df[, c("ENTREZID", "SYMBOL")])
  message(sprintf("%s (%s): %d genes", sheet, go_id, nrow(eg_df)))
  
  # 엑셀 파일로 저장
  out_file <- file.path(output_dir, sprintf("mouse_%s_genes.xlsx", sheet))
  write.xlsx(
    eg_df,
    file      = out_file,
    sheetName = sheet,
    rowNames  = FALSE
  )
}




library(AnnotationDbi)
library(org.Mm.eg.db)
library(openxlsx)

# 1) 처리할 GO term 목록과 시트 이름
go_terms <- c(
  "GO:0005615", # extracellular space
  "GO:0005576", # extracellular region
  "GO:0044421", # extracellular region part
  "GO:0070062", # extracellular exosome
  "GO:0046903", # secretion
  "GO:0050658", # secretion by cell
  "GO:0009306", # protein secretion
  )
sheet_names <- c(
  "extracellular_space",
  "extracellular_region",
  "extracellular_region_part",
  "extracellular_exosome",
  "secretion",
  "secretion_by_cell",
  "protein_secretion",
  "peptide_secretion"
)

# 2) org.Mm.eg.db 에서 유효한 GO ID 목록 가져오기
valid_go_ids <- keys(org.Mm.eg.db, keytype = "GO")

# 3) 워크북 생성
wb <- createWorkbook()

# 4) 각 GO term 별로 시트 추가 및 데이터 쓰기
for (i in seq_along(go_terms)) {
  go_id <- go_terms[i]
  sheet <- sheet_names[i]
  
  # GO ID 유효성 체크
  if (!go_id %in% valid_go_ids) {
    warning(sprintf("'%s'는 org.Mm.eg.db에 매핑된 유전자가 없어 건너뜁니다.", go_id))
    next
  }
  
  # EntrezID & SYMBOL 추출
  eg_df <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys    = go_id,
    keytype = "GO",
    columns = c("ENTREZID", "SYMBOL")
  )
  
  # 중복 제거 & 컬럼 정리
  eg_df <- unique(eg_df[, c("ENTREZID", "SYMBOL")])
  message(sprintf("%s (%s): %d genes", sheet, go_id, nrow(eg_df)))
  
  # 시트 추가
  addWorksheet(wb, sheetName = sheet)
  # 데이터 기록
  writeData(wb, sheet = sheet, eg_df)
}

# 5) 파일로 저장
output_file <- file.path(
  "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/",
  "mouse_extracellular_secretion_GO_genes.xlsx"
)
saveWorkbook(wb, file = output_file, overwrite = TRUE)














library(readxl)
library(dplyr)
library(openxlsx)

# 1) 공통 리스트 불러오기
common_file <- "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/common_values_updown_wo_talbert.xlsx"
up   <- read_excel(common_file, sheet = "Up")   %>% pull(1) %>% unique()
down <- read_excel(common_file, sheet = "Down") %>% pull(1) %>% unique()

# 2) 처리할 파일 경로 목록
go_files <- list(
  extracellular_space   = "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/mouse_extracellular_space_genes.xlsx",
  extracellular_region  = "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/mouse_extracellular_region_genes.xlsx",
  extracellular_exosome  = "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/mouse_extracellular_exosome_genes.xlsx",
  protein_secretion     = "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/mouse_protein_secretion_genes.xlsx",
  secretion_by_cell     = "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/mouse_secretion_by_cell_genes.xlsx",
  secretion             = "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/mouse_secretion_genes.xlsx",
  peptide_secretion     = "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/mouse_peptide_secretion_genes.xlsx"
)

# 3) 결과 워크북 생성
wb <- createWorkbook()

# 4) 각 파일별로 교집합 계산 후 시트에 기록
for (term in names(go_files)) {
  infile <- go_files[[term]]
  
  # 파일에서 SYMBOL 컬럼만 뽑아 유니크 벡터 생성
  genes <- read_excel(infile, sheet = 1) %>%
    rename_with(tolower) %>%
    pull(symbol) %>%
    unique()
  
  # Up/Down 교집합
  up_int   <- intersect(up,   genes)
  down_int <- intersect(down, genes)
  
  # 같은 길이로 맞추기
  nmax <- max(length(up_int), length(down_int))
  df   <- data.frame(
    up   = c(up_int,   rep(NA, nmax - length(up_int))),
    down = c(down_int, rep(NA, nmax - length(down_int)))
  )
  
  # 워크북에 시트 추가 & 데이터 쓰기
  addWorksheet(wb, sheetName = term)
  writeData(wb, sheet = term, df)
}

# 5) 최종 파일로 저장
output_file <- "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/common_vs_GO_intersections.xlsx"
saveWorkbook(wb, file = output_file, overwrite = TRUE)

message("✅ 모든 교집합 결과가 'common_vs_GO_intersections.xlsx' 파일에 시트별로 저장되었습니다.")






if (!requireNamespace("GO.db", quietly=TRUE))
  BiocManager::install("GO.db")

BiocManager::install("GO.db")
a

library(GO.db)
data(GOBPOFFSPRING)

# 패키지 로드
library(GO.db)
library(AnnotationDbi)
library(org.Mm.eg.db)

# 1) 부모 GO term 설정
go_parent <- "GO:0046903"  # secretion
go_parent <- "GO:0005615" #extra cellular space
# 2) BP ontology에서 이 term의 모든 후손(자식, 손자, ...) GO ID 가져오기
bp_offspring <- GOBPOFFSPRING[[go_parent]]
all_terms    <- c(go_parent, bp_offspring)


# 3) org.Mm.eg.db 에서 이 term들에 매핑된 유전자 추출
eg_all <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys    = all_terms,
  keytype = "GO",
  columns = c("ENTREZID", "SYMBOL")
)

# 4) 중복 제거
eg_all <- unique(eg_all[, c("ENTREZID", "SYMBOL")])

# 5) 결과 확인
cat("총", nrow(eg_all), "개의 유전자 추출됨\n")
head(eg_all)


write.xlsx(
  eg_all,
  file      = "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/extracellular_space_gene.xlsx",
  rowNames  = TRUE
)



# 2) 파일 경로 지정
common_file <- "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/common_values_updown_wo_talbert.xlsx"
secretion_file <- "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/secretion_gene.xlsx"
extracellular_space_file <- "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/extracellular_space_gene.xlsx"

# 3) Up / Down 리스트 불러오기
up   <- read_excel(common_file, sheet = "Up")   %>% pull(1) %>% unique()
down <- read_excel(common_file, sheet = "Down") %>% pull(1) %>% unique()

#    — 만약 시트명이 정확하지 않으면 sheet = 1, sheet = 2 처럼 숫자로 지정하셔도 됩니다.

# 4) Extracellular space 유전자 불러오기
extracellular_space <- read_excel(extracellular_space_file, sheet = 1) %>% 
  rename_with(tolower) %>%               # 컬럼명 소문자 통일
  pull(symbol) %>%                      # SYMBOL 컬럼만 벡터로
  unique()

# 5) 교집합 구하기
up_in_extracellular_space   <- intersect(up, extracellular_space)
down_in_extracellular_space <- intersect(down, extracellular_space)

up_in_extracellular_space
down_in_extracellular_space

# 6) 결과 확인
cat("Upregulated ∩ extracellular_space: ", length(up_in_extracellular_space), " genes\n")
cat("Downregulated ∩ extracellular_space: ", length(down_in_extracellular_space), " genes\n")

# 7) 결과를 엑셀로 저장 (옵션)

# 1) find the longest length
n <- max(length(up_in_extracellular_space), length(down_in_extracellular_space))

# 2) pad each vector to length n
up_padded   <- c(up_in_extracellular_space,   rep(NA, n - length(up_in_extracellular_space)))
down_padded <- c(down_in_extracellular_space, rep(NA, n - length(down_in_extracellular_space)))

# 3) build your data.frame
out_df <- data.frame(
  up_in_extracellular_space   = up_padded,
  down_in_extracellular_space = down_padded
)


out_df


write.xlsx(out_df,
           file      = "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Cachexia (통합본)/common_vs_extracellular_space_intersection.xlsx",
           sheetName = "intersection",
           rowNames  = FALSE)




