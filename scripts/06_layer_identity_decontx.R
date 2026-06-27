###############################################################################
# Step 2a-2:小数 layer 身份终极判别 + DecontX 可行性
#
# 背景:Step 2a(脚本05)已测到 raw_count_cellbender(=X)是小数,且 22.9% 位置
#       cb > cr(最大 2.88)——这不符合 CellBender 去污染"只减不增"的指纹,
#       说明它不是 CellBender 原始去污染counts,而更可能是被归一化/变换过。
#       本脚本用专家建议的 colSums 指纹 + log-norm 指纹把它的身份彻底钉死,
#       并探查 DecontX(de-novo 去污染,不需空液滴)的可行性。
#
# 判别逻辑:
#   指纹1 colSums:小数layer的每细胞总量 vs cellranger
#     - 系统性偏低且各细胞不同 → 像被去噪的 CellBender counts
#     - 被拉到常数附近(各细胞几乎相同) → 是 library-size 归一化过的
#   指纹2 log-norm:小数layer 是否 = log1p(1e4 * cellranger / colSums(cellranger))
#     - 吻合 → 它就是标准 log-normalized,不能当 CellBender counts 喂 CellChat
#   指纹3 DecontX:celda 包在不在,为 de-novo 去污染探路
###############################################################################

suppressPackageStartupMessages({
  library(zellkonverter)
  library(SingleCellExperiment)
  library(Matrix)
})

h5 <- "/Users/franxy/Documents/博士课题 6.15/严谨/data/Hahn 人 HFpEF 单核 snRNA/HFpEF_snRNAseq.h5ad"

message("== 读入 h5ad…")
sce <- readH5AD(h5, use_hdf5 = TRUE, reader = "R")

# 取前 3000 个细胞做判别(够了,省内存)
sub <- sce[, 1:3000]
cb <- as(assay(sub, "raw_count_cellbender"),  "CsparseMatrix")  # 小数 layer
cr <- as(assay(sub, "raw_count_cellranger"),  "CsparseMatrix")  # 整数 counts

## ===========================================================================
## 指纹 1:colSums(每细胞总量)比较
## ===========================================================================
cat("\n========== 指纹 1:colSums(每细胞总表达量)==========\n")
cs_cb <- Matrix::colSums(cb)
cs_cr <- Matrix::colSums(cr)

cat(sprintf("cellranger colSums:  范围 [%.1f, %.1f],中位 %.1f,变异系数 %.3f\n",
            min(cs_cr), max(cs_cr), median(cs_cr), sd(cs_cr)/mean(cs_cr)))
cat(sprintf("小数layer  colSums:  范围 [%.2f, %.2f],中位 %.2f,变异系数 %.3f\n",
            min(cs_cb), max(cs_cb), median(cs_cb), sd(cs_cb)/mean(cs_cb)))

cat(sprintf("\n小数layer/cellranger 的 colSums 比值: 中位 %.4f,范围 [%.4f, %.4f]\n",
            median(cs_cb/cs_cr), min(cs_cb/cs_cr), max(cs_cb/cs_cr)))

cat("\n--- 指纹1 判读 ---\n")
cv_cb <- sd(cs_cb)/mean(cs_cb)
if (cv_cb < 0.05) {
  cat(">> 小数layer 的 colSums 变异系数极小(各细胞总量几乎相同)\n")
  cat("   → 强烈提示它被 library-size 归一化过(每细胞拉到相同总量)\n")
  cat("   → 它【不是】CellBender 去污染counts,不能直接喂 CellChat(会双重归一化)\n")
} else {
  cat(">> 小数layer 的 colSums 各细胞差异较大(变异系数 > 0.05)\n")
  cat("   → 不像被简单 library-size 归一化;需结合指纹2判断\n")
}

## ===========================================================================
## 指纹 2:小数 layer 是否 = log1p(1e4-normalized cellranger)
## ===========================================================================
cat("\n========== 指纹 2:是否为标准 log-normalize ==========\n")
# 标准 scanpy log-normalize:normalize_total(target_sum=1e4) + log1p
cr_norm <- cr
cr_norm@x <- cr_norm@x  # copy
# 每列除以该列总和、乘 1e4
scaling <- 1e4 / cs_cr
cr_lognorm <- cr
cr_lognorm <- as(cr_lognorm, "CsparseMatrix")
# 用对角缩放实现按列归一化
D <- Diagonal(x = scaling)
cr_lognorm <- cr_lognorm %*% D            # 每列 *scaling
cr_lognorm@x <- log1p(cr_lognorm@x)       # log1p

# 对齐非零位置比较小数layer 和 重算的 log-norm
cb_tri <- summary(cb)
ln_tri <- summary(cr_lognorm)
key_cb <- paste(cb_tri$i, cb_tri$j)
key_ln <- paste(ln_tri$i, ln_tri$j)
common <- intersect(key_cb, key_ln)
icb <- match(common, key_cb); iln <- match(common, key_ln)

diff <- abs(cb_tri$x[icb] - ln_tri$x[iln])
corr <- cor(cb_tri$x[icb], ln_tri$x[iln])

cat(sprintf("对齐位置数: %d\n", length(common)))
cat(sprintf("小数layer vs 重算log-norm: 相关 r = %.4f\n", corr))
cat(sprintf("绝对差: 中位 %.4f,均值 %.4f,最大 %.4f\n",
            median(diff), mean(diff), max(diff)))
cat("\n小数layer 与 重算log-norm 的并排示例(前8个非零位置):\n")
show_n <- min(8, length(common))
print(data.frame(
  小数layer    = round(cb_tri$x[icb][1:show_n], 4),
  重算lognorm  = round(ln_tri$x[iln][1:show_n], 4)
))

cat("\n--- 指纹2 判读 ---\n")
if (corr > 0.99 && median(diff) < 0.05) {
  cat(">> 小数layer 与标准 log-normalize 高度吻合(r>0.99,差异极小)\n")
  cat("   → 【钉死】小数layer = log1p(1e4-normalized CellRanger),是 log-normalized 表达\n")
  cat("   → 之前 CellChat 用 raw_count_cellranger 当输入是对的;\n")
  cat("      这个小数layer 不能当 CellBender counts(它是归一化数据,非counts)\n")
} else {
  cat(sprintf(">> 与标准 log-normalize 不完全吻合(r=%.3f)\n", corr))
  cat("   → 它可能是 CellBender-corrected 后再做了某种变换,或其他处理;\n")
  cat("      需进一步判断(把上面示例值贴给我)\n")
}

## ===========================================================================
## 指纹 3:DecontX(celda)可行性 —— de-novo 去污染,不需空液滴
## ===========================================================================
cat("\n========== 指纹 3:DecontX 可行性 ==========\n")
has_celda <- requireNamespace("celda", quietly = TRUE)
cat(sprintf("celda 包(含 DecontX)是否已安装: %s\n", has_celda))
if (has_celda) {
  cat(sprintf("  celda 版本: %s\n", as.character(packageVersion("celda"))))
  cat("  → 可直接用 raw_count_cellranger + 细胞类型标签跑 DecontX 做 de-novo 去污染\n")
} else {
  cat("  → 未安装。下一步可装:BiocManager::install('celda')\n")
  cat("    DecontX 从 filtered matrix + cluster labels 估计污染,不需空液滴,\n")
  cat("    是'只有处理后矩阵'场景做 de-novo 去污染的对口工具(Yang et al. Genome Biol 2020)\n")
}

## ===========================================================================
## 结论
## ===========================================================================
cat("\n========== 结论 ==========\n")
cat("把以上三个指纹的输出贴回对话,据此确定:\n")
cat("  1. 小数layer 到底是什么(log-norm? 还是 post-CellBender?)\n")
cat("  2. 能否用它当 CellBender counts 重跑 CellChat(若是 log-norm 则不能)\n")
cat("  3. DecontX 这条 de-novo 去污染路要不要走(大概率要,作为主稳健性证据)\n")
