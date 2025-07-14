library(ggplot2)
library(readxl)
library(dplyr)

# 파일 목록 (예: 8개의 엑셀 파일)
file_list <- c(
  "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Bulk_seq_analysis_practice/250303_EMBO(sophia)/DEG_results(EMBO).xlsx", 
  "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Bulk_seq_analysis_practice/250304_JCI(joshua)/DEG_results(joshua_SC).xlsx", 
  "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Bulk_seq_analysis_practice/250304_JCI(joshua)/DEG_results(joshua).xlsx",
  "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Bulk_seq_analysis_practice/JEM_bulk(Rupert)/DEG_results(Rupert).xlsx",
  "/Users/daehwankim/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Seqencing_Practicing/Bulk_seq_analysis_practice/JNCI_bulk (Tseng)/DEG_results(Tseng).xlsx"
)


# 데이터프레임 초기화
deg_data <- data.frame()


#모든 파일 읽기
for (file in file_list) {
  condition <- tools::file_path_sans_ext(basename(file))  # 파일명을 condition으로 사용
  
  # Upregulated 시트 읽기
  up_data <- read_excel(file, sheet = "Upregulated") %>%
    mutate(condition = condition, direction = "Up")
  
  # Downregulated 시트 읽기
  down_data <- read_excel(file, sheet = "Downregulated") %>%
    mutate(condition = condition, direction = "Down")
  
  # 통합
  deg_data <- bind_rows(deg_data, up_data, down_data)
}

head(deg_data)

# p-value 변환 및 Miami Plot을 위한 y값 조정
deg_data <- deg_data %>%
  mutate(negLogP = ifelse(direction == "Up", -log10(padj), log10(padj)))


# Miami plot 생성
ggplot(deg_data, aes(x = condition, y = log2FoldChange, color = direction)) +
  geom_vline(xintercept = 0.5, linetype = "solid", color = "black", linewidth = 0.25) +  # y축 왼쪽 세로선 추가
  geom_jitter(alpha = 1.0, width = 0.25, size = 0.2) +  # 점을 약간 퍼트려서 겹침 방지
  scale_color_manual(values = c("Up" = "red3",  "Down" = "blue4")) +
  labs(x = "Dataset", y = "log2FoldChange", title = "Miami Plot for Multiple DEG Analyses") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # x축 글자 회전 글자 없에고 싶으면, element_blank()
    panel.grid.major = element_blank(),  # 주요 격자선 제거
    panel.grid.minor = element_blank())   # 보조 격자선 제거

#ㅇ여러 색 조합
c("Up" = "#e64949",  "Down" = "#4989e6")
c("Up" = "#d63c3c",  "Down" = "#3c76d6")

element_text()

mycolors <- c("Downregulated" = "#00AFBB", "Upregulated" = "#bb0c00", "Non significant" = "grey")