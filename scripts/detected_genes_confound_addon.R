###############################################################################
# 补算 per-donor "detected genes"(median nFeature)作为 confound 平衡检查
#
# 背景:图审核发现 Results §3.1 曾列 detected genes(P=0.72),但
#       step3_confound_check.csv 无此列、fig5 未画 → 已暂从正文删。
#       本脚本从 h5ad 真算一遍,据结果决定补回 or 维持删。
#
# 维度区分(两个不同的技术维度):
#   · median_depth   = 每核 UMI 数(测序量,CSV 已有)
#   · detected genes = 每核检出的不同基因数(nFeature,文库复杂度,本次补)
#
# ★ 决策(跑完看打印的 P 与 Cliff's δ):
#   · ns(P>0.05)              → 平衡,可补回 fig5/正文
#   · 显著但 δ<0(HFpEF 更低)  → conservative(同 depth),可用,正文点明方向
#   · 显著且 δ>0(HFpEF 更高)  → 不利,弃(维持已删)
#
# 数据层:raw_count_cellranger(干净整数 counts)
# 安全:不覆盖原 CSV,结果先写 step3_confound_check_withGenes.csv 供你核
###############################################################################

source("/Users/franxy/Documents/博士课题 6.15/FLRT2&UNC5C/fig_common.R")  # two_group_summary, GRP_LEVELS
suppressPackageStartupMessages({ library(Matrix); library(dplyr) })

CONF  <- "/Users/franxy/Documents/博士课题 6.15/严谨/results/step3_confound_check.csv"
H5AD  <- "/Users/franxy/Documents/博士课题 6.15/严谨/data/Hahn 人 HFpEF 单核 snRNA/HFpEF_snRNAseq.h5ad"
LAYER <- "raw_count_cellranger"

conf <- read.csv(CONF, stringsAsFactors = FALSE)
conf$donor <- as.character(conf$donor)

## === 取每核 counts(基因×核, dgCMatrix)+ 每核所属 donor ====================
## 路径 A(优先,最省时):若你生成 step3 的对象还在内存(SCE 叫 sce / Seurat 叫
##   seu),直接用它的 raw_count_cellranger,跳过重读 h5ad——把下面 B 块注释掉,
##   自己给 m(基因×核 dgCMatrix)和 cd(含 donor 的 data.frame)即可。
##
## 路径 B(自包含):zellkonverter 从 h5ad 重读。首次会装 basilisk 环境,稍慢。
sce <- zellkonverter::readH5AD(H5AD, use_hdf5 = FALSE)   # 48866 核 raw counts 进内存 ~3-4GB
stopifnot(LAYER %in% SummarizedExperiment::assayNames(sce))
m  <- SummarizedExperiment::assay(sce, LAYER)            # 基因 × 核
cd <- as.data.frame(SummarizedExperiment::colData(sce))

## === 每核 detected genes(非零基因数;sparse 列指针差分,不 densify)========
m <- as(m, "CsparseMatrix")                              # 确保 CSC
nfeat_per_nucleus <- diff(m@p)                            # 每列非零元素数 = nFeature
stopifnot(length(nfeat_per_nucleus) == ncol(m))

## === 自动识别 donor 列(与 CSV donor 重合度最高那列,免猜列名)=============
score <- sapply(cd, function(x) mean(as.character(x) %in% conf$donor))
donor_col <- names(which.max(score))
cat(sprintf("[donor 列] 自动识别 = '%s'(重合度 %.0f%%)\n", donor_col, 100 * max(score)))
if (max(score) <= 0.5) {                                  # scalar 条件用 if 不用 ifelse
  cat("  colData 各列:\n"); print(colnames(cd))
  stop("没匹配上 donor 列(重合度<50%);请手动设 donor_col 后重跑")
}

## === per-donor median nFeature → merge 进 confound 表 ======================
per_donor <- data.frame(donor = as.character(cd[[donor_col]]),
                        nfeat = nfeat_per_nucleus) %>%
  group_by(donor) %>% summarise(detected_genes = median(nfeat), .groups = "drop")

out <- conf %>% left_join(per_donor, by = "donor")
stopifnot(!anyNA(out$detected_genes))                     # 每个 donor 都配上

## === 同口径 Wilcoxon(与 fig5 一致:two_group_summary, exact)===============
s <- two_group_summary(out$detected_genes, factor(out$group, levels = GRP_LEVELS))
cat("\n==== detected genes / donor:HFpEF vs control ====\n")
cat(sprintf("  median  control = %.0f    HFpEF = %.0f\n", s$median_control, s$median_HFpEF))
cat(sprintf("  Wilcoxon P = %.3g     Cliff's δ = %+.2f\n", s$p_wilcox, s$cliffs_d))
if (s$p_wilcox > 0.05) {
  cat("  → ns,平衡。可补回 fig5/正文。\n")
} else if (s$cliffs_d < 0) {
  cat("  → 显著但 HFpEF 更低(conservative,同 depth)。可用,正文点明方向。\n")
} else {
  cat("  → 显著且 HFpEF 更高(不利)。建议维持已删,不补。\n")
}

## === 不覆盖原文件,写副本供你核 ============================================
OUT <- sub("step3_confound_check.csv", "step3_confound_check_withGenes.csv", CONF)
write.csv(out, OUT, row.names = FALSE)
cat(sprintf("\n已写:%s\n把上面 median / P / δ / 方向那几行贴回给我:好就更新 fig5 + 正文。\n", OUT))
