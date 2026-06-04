#####################################################################
# Experiment 2 -- direct substitution of published gene sets
#
# BCS40010 project (Kim Na-hye / Lee Jae-hyun)
#
# Question: which published circadian gene set transfers best when
#           dropped into the TimeSignature elastic-net framework?
#
# Sets compared (all trained on GSE39445, validated on GSE48113 /
# GSE56931 / GSE113883 with two blood draws per subject):
#   * TScore18      -- Braun 2018 core-18 predictors (Table 1)
#   * TimeMachine37 -- Huang & Braun 2024 JTK panel (37 genes)
#   * tauFisher     -- Duan & Ngo 2024 core-clock + top-10 rhythmic
#   * clockOnly     -- canonical core-clock genes only
#   * intersection  -- genes shared by the three published sets
# A full-pool model (all 7,615 genes) is included as the reference.
#
# Assumes 00_setup.R has been sourced.
#####################################################################

if (!exists("jtkRanked")) source("00_setup.R")

#--- assemble the gene sets (literature + JTK-derived panels) --------
sets <- buildGeneSets(jtk, available = rownames(all.expr),
                      nTimeMachine = 37, nTauTop = 10)
sets$fullPool <- rownames(all.expr)         # original TimeSignature reference

cat("\n--- Experiment 2 gene sets (size, members) ---\n")
for (nm in names(sets)) {
	g <- sets[[nm]]
	cat(sprintf("  %-13s n=%4d", nm, length(g)))
	if (length(g) <= 20) cat("  : ", paste(g, collapse = ", "))
	cat("\n")
}

#--- train + validate each set ---------------------------------------
exp2 <- do.call(rbind, lapply(names(sets), function(nm) {
	ts  <- trainTSpool(trainWSN, trainSubjs, trainTimes,
	                   genes = sets[[nm]], alpha = 0.5, s = NULL, seed = 194)
	res <- evalAcrossStudies(ts, calibAll, all.meta,
	                         studies = c("V1", "V2", "V3"))
	res$set <- nm
	cat(sprintf("  %-13s nGene=%4d selected=%2d  Valid.all: nAUC=%.3f MAE=%.2f %%2h=%.0f\n",
	            nm, length(sets[[nm]]), ts$nSelected,
	            res$nAUC[res$study == "Valid.all"],
	            res$MAE[res$study == "Valid.all"],
	            res$pct2h[res$study == "Valid.all"]))
	res
}))

write.csv(exp2, "output/Exp2_geneset_substitution_results.csv", row.names = FALSE)

#--- summary table (pooled validation) -------------------------------
library(ggplot2)

# 1. x축 라벨을 '이름\n(n=유전자수)' 형태로 미리 생성
exp2.pooled$label <- sprintf("%s\n(n=%d)", exp2.pooled$set, exp2.pooled$nGene)

# 2. x축 막대 순서가 알파벳순으로 꼬이지 않도록 기존 데이터프레임 순서로 고정(Factor 화)
exp2.pooled$label <- factor(exp2.pooled$label, levels = exp2.pooled$label)

pdf("output/Exp2_geneset_nAUC.pdf", width = 7.5, height = 5)

p2 <- ggplot(exp2.pooled, aes(x = label, y = nAUC)) +
  # 바 차트 생성 (steelblue 색상, 두께와 투명도 조절로 세련미 추가)
  geom_col(fill = "steelblue", width = 0.55, alpha = 0.9) +
  
  # Clinical floor (nAUC = 0.80) 가이드라인 및 텍스트
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "firebrick", linewidth = 0.7) +
  annotate("text", x = 5, y = 0.80 - 0.03, 
           label = "clinical floor nAUC = 0.80", 
           color = "firebrick", size = 3.5, fontface = "italic", hjust = 0) +
  
  # 막대 위에 nAUC 수치 표시 (소수점 3자리까지 표시하여 정밀도 향상)
  geom_text(aes(label = sprintf("%.3f", nAUC)), vjust = -0.8, size = 4, fontface = "bold", color = "black") +
  
  # y축 범위 설정 (0부터 1까지 꽉 채우되, 글자가 잘리지 않도록 상단 여유 1.05 부여)
  scale_y_continuous(limits = c(0, 1.05), expand = c(0, 0)) +
  
  # 테마 및 라벨링
  labs(
    title = "Exp 2: Published gene sets in the TimeSignature framework",
    x = NULL, # x축 라벨에 이미 정보가 충분하므로 축 제목은 생략
    y = "nAUC (pooled validation)"
  ) +
  theme_classic() + # 막대 그래프에 최적화된 깔끔한 클래식 테마
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 15)),
    axis.text.x = element_text(size = 10, color = "black", margin = margin(t = 10)),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title.y = element_text(margin = margin(r = 10), color = "black"),
    axis.ticks.x = element_blank() # x축의 자잘한 눈금선 제거로 깔끔함 극대화
  )

print(p2)
dev.off()

#--- summary text output ---------------------------------------------
best <- exp2.pooled[which.max(exp2.pooled$nAUC), ]
cat(sprintf("\nBest-transferring set: %s (n=%d genes, nAUC=%.3f, MAE=%.2f h)\n",
            best$set, best$nGene, best$nAUC, best$MAE))
cat("Exp 2 done -> output/Exp2_geneset_substitution_results.csv, output/Exp2_geneset_nAUC.pdf\n")