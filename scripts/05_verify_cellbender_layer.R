###############################################################################
# Step 2a(承重闸门):坐实 raw_count_cellbender 确实是 CellBender 去环境RNA后的counts
#
# 为什么必须先做(专家意见):
#   把一个来源未坐实的 layer 喂进 CellChat 的 normalizeData,如果它其实是
#   normalized/imputed 输出,就会"在归一化上再叠归一化",引入新 artifact。
#   所以用它当 CellChat 输入前,必须先证实它是 CellBender 原始(去污染)counts。
#
# 判断标准(三问全为"是"才算闸门通过):
#   问题1:layer 命名 + 数值特征  —— cellranger 应为整数,cellbender 应为小数
#   问题2:两者关系               —— 对非零位置,cellbender/cellranger 应落在 (0,1],
#                                    即 CellBender ≤ CellRanger(去污染只减不增),
#                                    符合"去污染期望counts = 原始 × P(signal)"
#   问题3:原文是否跑了 CellBender —— 已从 Hahn 原文证实(见下方注释),作为命名佐证
#
# 【已证实的原文事实(来自你上传的 Circulation Research 2025 正文 + 补充材料)】:
#   - 正文 Methods:CellRanger v4.0.0 + CellBender remove-background v0.2.0,FPR=0.01
#   - 补充材料:基因需在 CellBender 与 CellRanger 两套counts里都显著且方向一致
#   => Hahn 确实跑了 CellBender;h5ad 里名为 raw_count_cellbender 的 layer
#      命名与原文处理一致。本脚本再用数值关系做最后坐实。
###############################################################################

suppressPackageStartupMessages({
  library(zellkonverter)
  library(SingleCellExperiment)
  library(Matrix)
})

h5 <- "/Users/franxy/Documents/博士课题 6.15/严谨/data/Hahn 人 HFpEF 单核 snRNA/HFpEF_snRNAseq.h5ad"

message("== 读入 h5ad(磁盘模式)…")
sce <- readH5AD(h5, use_hdf5 = TRUE, reader = "R")

## ---- 问题 1:所有 assay 的命名 + 数值特征 ----------------------------------
cat("\n========== 问题 1:assay 命名与数值特征 ==========\n")
cat("assay 名:\n"); print(assayNames(sce))

# 抽一个细胞类型的小子集来快速检查数值特征(省内存)
# 用前 2000 个细胞即可判断整数/小数/范围
sub <- sce[, 1:2000]
for (a in assayNames(sub)) {
  m <- as(assay(sub, a), "CsparseMatrix")
  vals <- m@x                       # 非零值
  is_int <- all(vals == round(vals))
  cat(sprintf("\n--- assay '%s' ---\n", a))
  cat(sprintf("  是否全整数: %s\n", is_int))
  cat(sprintf("  非零值范围: [%.4f, %.2f]\n", min(vals), max(vals)))
  cat(sprintf("  非零值示例(前5): %s\n", paste(round(head(vals,5),4), collapse=", ")))
  cat(sprintf("  稀疏度(非零占比): %.3f%%\n", 100*length(vals)/(nrow(m)*ncol(m))))
}

## ---- 问题 2:cellbender 与 cellranger 的关系(核心)-------------------------
cat("\n========== 问题 2:cellbender vs cellranger 关系 ==========\n")
cb <- as(assay(sub, "raw_count_cellbender"), "CsparseMatrix")
cr <- as(assay(sub, "raw_count_cellranger"), "CsparseMatrix")

# 找两者都非零的位置,算比值 cb/cr
# 转成三元组对齐
cb_tri <- summary(cb)   # i, j, x
cr_tri <- summary(cr)
key_cb <- paste(cb_tri$i, cb_tri$j)
key_cr <- paste(cr_tri$i, cr_tri$j)
common <- intersect(key_cb, key_cr)
cat(sprintf("两 layer 都非零的位置数: %d\n", length(common)))

idx_cb <- match(common, key_cb)
idx_cr <- match(common, key_cr)
ratio <- cb_tri$x[idx_cb] / cr_tri$x[idx_cr]

cat(sprintf("\ncb/cr 比值统计(对两者都非零的位置):\n"))
cat(sprintf("  最小值: %.4f\n", min(ratio)))
cat(sprintf("  最大值: %.4f\n", max(ratio)))
cat(sprintf("  中位数: %.4f\n", median(ratio)))
cat(sprintf("  比值 > 1 的占比: %.3f%%(应接近 0%%——去污染只减不增)\n",
            100*mean(ratio > 1.0001)))
cat(sprintf("  比值 ≤ 1 的占比: %.3f%%\n", 100*mean(ratio <= 1.0001)))

# 关键判断
cat("\n--- 关系判断 ---\n")
n_cb_only <- length(setdiff(key_cb, key_cr))  # cellbender 有、cellranger 没有的位置
cat(sprintf("cellbender 非零但 cellranger 为零的位置数: %d\n", n_cb_only))
cat("  (去污染不应凭空产生新计数,此数应为 0 或极少)\n")

if (max(ratio) <= 1.0001 && n_cb_only == 0) {
  cat("\n>> 判断:cb ≤ cr 且无凭空新增 → 符合 CellBender 去污染特征 ✓\n")
} else {
  cat("\n>> 警告:出现 cb > cr 或凭空新增计数 → 这个 layer 可能不是简单的\n")
  cat("   CellBender 去污染counts(也许是 normalized/imputed),需停下来核查!\n")
}

## ---- 问题 3:原文事实(已证实,见脚本头注释)-------------------------------
cat("\n========== 问题 3:Hahn 原文是否跑 CellBender ==========\n")
cat("已从原文证实(见脚本头部注释):\n")
cat("  - 正文 Methods:CellRanger v4.0.0 + CellBender remove-background v0.2.0\n")
cat("  - 补充材料:CellBender 与 CellRanger 双层一致性过滤\n")
cat("  => layer 命名 'raw_count_cellbender' 与原文处理一致。\n")

## ---- 闸门结论 --------------------------------------------------------------
cat("\n========== 闸门结论 ==========\n")
cat("若问题1(cellranger整数/cellbender小数)、问题2(cb≤cr无新增)、问题3(原文跑了CellBender)\n")
cat("三者皆成立 → raw_count_cellbender 坐实为 CellBender 去污染counts,可作 CellChat 输入。\n")
cat("把上面输出贴回对话,据此决定 Step 2c 用哪个 layer 重跑 CellChat。\n")
