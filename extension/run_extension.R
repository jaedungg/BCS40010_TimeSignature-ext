#####################################################################
# Master driver for the TimeSignature gene-selection extension
#
# BCS40010 project (Kim Na-hye / Lee Jae-hyun)
#
# Runs both experiments and assembles the headline deliverable:
# a single gene-count vs accuracy trade-off curve (Exp 1 JTK pool
# sweep) with the published gene sets (Exp 2) overlaid as points.
#
# Usage (working directory = TimeSignatR/extension/):
#   Rscript run_extension.R
#####################################################################

library(ggplot2)
library(ggrepel)

source("00_setup.R")
source("Exp1_JTK_prefilter.R")
source("Exp2_geneset_substitution.R")

#--- combined trade-off figure ---------------------------------------
e1 <- exp1.pooled[order(exp1.pooled$nGene), ]
e2 <- exp2.pooled[exp2.pooled$set != "fullPool", ]

pdf("output/Fig_tradeoff_combined.pdf", width = 8, height = 5.5)

p <- ggplot() +
  # 1. Clinical floor (nAUC = 0.80) 가이드라인
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "firebrick", alpha = 0.7) +
  annotate("text", x = 4, y = 0.802, label = "clinical floor nAUC = 0.80", 
           color = "firebrick", size = 3.5, hjust = 0, fontface = "italic") +
  
  # 2. Exp 1: JTK-filtered top-K pool (선과 점)
  geom_line(data = e1, aes(x = nGene, y = nAUC), color = "navy", alpha = 0.6) +
  geom_point(data = e1, aes(x = nGene, y = nAUC, shape = "Exp 1: JTK-filtered top-K pool"), 
             color = "navy", size = 2.5) +
  
  # 3. Exp 2: Published gene sets (포인트)
  geom_point(data = e2, aes(x = nGene, y = nAUC, color = set, shape = "Exp 2: published gene sets"), 
             size = 3.5) +
  
  # 4. Exp 2: 겹치지 않는 라벨 (ggrepel 사용)
  geom_text_repel(data = e2, aes(x = nGene, y = nAUC, label = set, color = set),
                  size = 3.5, fontface = "bold", box.padding = 0.6, 
                  point.padding = 0.5, show.legend = FALSE) +
  
  # 5. 스케일 및 축 범위 설정 (원래 범위 완벽 재현)
  scale_x_log10(breaks = c(10, 50, 100, 500, 1000, 5000)) +
  scale_shape_manual(name = NULL, values = c("Exp 1: JTK-filtered top-K pool" = 16, 
                                             "Exp 2: published gene sets" = 17)) +
  scale_color_brewer(palette = "Set1") + 
  guides(color = "none") + 
  
  # 6. 테마 설정 (그리드 제거 및 범례 내부 배치 반영)
  labs(
    title = "TimeSignature gene-count vs accuracy trade-off",
    x = "Number of genes in candidate pool (log scale)",
    y = "nAUC (pooled validation: GSE48113/56931/113883)"
  ) +
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(),  # 주 격자선 제거
    panel.grid.minor = element_blank(),  # 보조 격자선 제거
    legend.position = c(0.98, 0.02),      # 그래프 내부 우측 하단 배치 (상대좌표 x=0.98, y=0.02)
    legend.justification = c("right", "bottom"),
    legend.background = element_blank(),  # 범례 박스 배경 투명화 (원래 코드의 bty = "n")
    legend.key = element_blank(),         # 범례 아이콘 심볼 배경 투명화

    plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 15)), # 메인 타이틀 아래쪽(b) 여백 추가
    axis.title.x = element_text(margin = margin(t = 15)), # x축 제목 위쪽(t) 여백 추가 (그림과 멀어짐)
    axis.title.y = element_text(margin = margin(r = 15))  # y축 제목 오른쪽(r) 여백 추가 (그림과 멀어짐)
  )

print(p)
dev.off()

#--- master summary table --------------------------------------------
master <- rbind(
	data.frame(experiment = "Exp1-JTKpool",
	           label = paste0("top", e1$K),
	           nGene = e1$nGene, nSel = e1$nSel,
	           MAE = e1$MAE, nAUC = e1$nAUC, pct2h = e1$pct2h),
	data.frame(experiment = "Exp2-geneset",
	           label = exp2.pooled$set,
	           nGene = exp2.pooled$nGene, nSel = exp2.pooled$nSel,
	           MAE = exp2.pooled$MAE, nAUC = exp2.pooled$nAUC,
	           pct2h = exp2.pooled$pct2h)
)
write.csv(master, "output/master_summary.csv", row.names = FALSE)

cat("\n=====================================================\n")
cat("Extension complete. Outputs in extension/output/:\n")
cat("  Exp1_JTK_prefilter_results.csv\n")
cat("  Exp1_tradeoff.pdf\n")
cat("  Exp2_geneset_substitution_results.csv\n")
cat("  Exp2_geneset_nAUC.pdf\n")
cat("  Fig_tradeoff_combined.pdf   <- headline deliverable\n")
cat("  master_summary.csv\n")
cat("=====================================================\n")
