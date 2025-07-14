# 0) 패키지 로드

BiocManager::install("ggVennDiagram")

library(readxl)
library(ggVennDiagram)
library(sf)

# 파일 목록 (예: 8개의 엑셀 파일)
file_list <- c(
  "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Bulk_seq_analysis_practice/250303_EMBO(sophia)/DEG_results(EMBO).xlsx", 
  "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Bulk_seq_analysis_practice/250304_JCI(joshua)/DEG_results(joshua_SC).xlsx", 
  "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Bulk_seq_analysis_practice/250304_JCI(joshua)/DEG_results(joshua).xlsx",
  "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Bulk_seq_analysis_practice/JEM_bulk(Rupert)/DEG_results(Rupert).xlsx",
  "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Bulk_seq_analysis_practice/JNCI_bulk (Tseng)/DEG_results(Tseng).xlsx"
)

names(file_list) <- c("GSE138464","GSE142455_(SC)","GSE142455_(Orth)","GSE123310","GSE65936")
names(file_list)

# 2) Up / Down 리스트 생성
#    – 시트 이름이 각각 "upregulated" / "downregulated" 라고 가정
deg_lists_up <- lapply(file_list, function(path) {
  df_up <- read_xlsx(path, sheet = "Upregulated")
  unique(df_up$genesymbol)
})
deg_lists_down <- lapply(file_list, function(path) {
  df_dn <- read_xlsx(path, sheet = "Downregulated")
  unique(df_dn$genesymbol)
})

names(deg_lists_up)
names(deg_lists_down)

# 4) ggVennDiagram으로 5세트 Venn 그리기
ggVennDiagram(deg_lists_up, label = "count", label_alpha = 0, 
              label_size = 3, set_color = "black", set_size = 2, 
              edge_size = 0.5,  mapping = aes(fill = degree)) +
  scale_fill_gradient(low = "white", high = "red") +
  theme(legend.position = "none") +
  labs(title = "DEGs Upregulated")

help(scale_fill_gradient)


ggVennDiagram(deg_lists_down, label = "count", label_alpha = 0, label_size = 3) +
  scale_fill_gradient(low = "white", high = "white") +
  theme(legend.position = "none") +
  labs(title = "DEGs Downregulated")


help("ggVennDiagram")


c("Up" = "#d63c3c",  "Down" = "#3c76d6")


#force_upset | FALSE | TRUE로 설정하면 Venn이 아닌 UpSet 스타일로 강제 변환합니다. |



# 1) 리스트 → Venn 객체 생성
venn_up <- Venn(deg_lists_up)
vd_up   <- process_data(venn_up)

# 2) 교집합 정보 테이블
vd_up$regionLabel$id
up_ids <- vd_up$regionLabel$id

sets_list <- strsplit(up_ids, "/", fixed = TRUE) # 1) 문자열을 "/" 기준으로 분리하여 리스트로

# 2) 각 요소 길이를 재서 degree 벡터로
degree <- lengths(sets_list)

data.frame(
  id     = up_ids,
  degree = degree
)

vd_up$regionLabel$degree <- degree
vd_up$regionLabel


vd_up$regionData$degree <- degree

View(vd_up$regionData)





