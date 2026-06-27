###############################################################################
# Step 2a-2(第3版):小数 layer 身份终极判别 —— 修正逻辑矛盾 + 反解缩放因子
#
# 相比 06b 的修正(按代码审查):
#   1.【逻辑矛盾,必须改】struct_match=TRUE 在逻辑上排除 normalize(CellBender):
#       归一化不改稀疏结构,CellBender 去污染会把 ambient 项削成 0(减少非零)。
#       所以"结构一致"=无任何 CellBender 削除痕迹=与独立CellBender无关。
#       - branch 2(结构一致但三候选不吻合)正确解释 = CellRanger 的某种未测归一化
#       - normalize(CellBender) 假设只能进 branch 3,签名 n_cr_only>0 且 n_cb_only≈0
#   2. med_target 在【全量】上算(scanpy normalize_total 默认用全数据中位数)
#   3. 加【反解缩放因子】判据:逐细胞反推隐含 target,免疫于 target 猜错
#   4. 去掉 <<- 全局变量,缩放向量当参数传;删重复 set.seed
###############################################################################

suppressPackageStartupMessages({
  library(zellkonverter)
  library(SingleCellExperiment)
  library(Matrix)
})

h5 <- "/Users/franxy/Documents/博士课题 6.15/严谨/data/Hahn 人 HFpEF 单核 snRNA/HFpEF_snRNAseq.h5ad"

message("== 读入 h5ad…")
sce <- readH5AD(h5, use_hdf5 = TRUE, reader = "R")

# 先在全量上算 cellranger 每细胞总量(供 scanpy 默认 median target)
message("== 全量计算 cellranger colSums(供 median target,分块流式)…")
cs_cr_full <- Matrix::colSums(assay(sce, "raw_count_cellranger"))
med_target_full <- median(cs_cr_full)
cat(sprintf("全量 cellranger 每细胞总量中位数(scanpy 默认 target)= %.2f\n", med_target_full))

# 随机抽 3000 细胞做判别
set.seed(1)
sel <- sort(sample(ncol(sce), 3000))
sub <- sce[, sel]
cb <- as(assay(sub, "raw_count_cellbender"), "CsparseMatrix")
cr <- as(assay(sub, "raw_count_cellranger"), "CsparseMatrix")

## ---- 方向断言 -------------------------------------------------------------
cat("\n========== 方向断言 ==========\n")
cat(sprintf("dim(cb) = %d × %d(应为 ~36601 基因 × 3000 细胞)\n", nrow(cb), ncol(cb)))
stopifnot(nrow(cb) > ncol(cb))
cat(">> 方向正确:行=基因,列=细胞 ✓\n")

## ---- 取值范围 -------------------------------------------------------------
cat("\n========== 取值范围(决定性线索)==========\n")
cat(sprintf("小数layer 非零值范围: [%.4f, %.4f]\n", min(cb@x), max(cb@x)))
cat(sprintf("cellranger 非零值范围: [%.0f, %.0f]\n", min(cr@x), max(cr@x)))
cat("判读:max≈9 上下 → log1p 形态;max 到几百~上万 → 线性 CP10K(未log)\n")

## ---- 结构断言(双向 setdiff)----------------------------------------------
cat("\n========== 结构断言 ==========\n")
cb_tri <- summary(cb); cr_tri <- summary(cr)
key_cb <- paste(cb_tri$i, cb_tri$j)
key_cr <- paste(cr_tri$i, cr_tri$j)
common <- intersect(key_cb, key_cr)
n_cb_only <- length(setdiff(key_cb, key_cr))
n_cr_only <- length(setdiff(key_cr, key_cb))
cat(sprintf("cb 非零: %d   cr 非零: %d   交集: %d\n", length(key_cb), length(key_cr), length(common)))
cat(sprintf("cb 独有(cr为0): %d   cr 独有(cb为0): %d\n", n_cb_only, n_cr_only))
struct_match <- (length(key_cb) == length(key_cr)) && (n_cb_only == 0) && (n_cr_only == 0)
cat(sprintf(">> struct_match = %s\n", struct_match))

## ---- 三候选归一化比对(缩放向量当参数传,无全局变量)----------------------
cat("\n========== 三候选归一化比对 ==========\n")
cs_cr <- Matrix::colSums(cr)
cs_cb <- Matrix::colSums(cb)
cat(sprintf("cellranger colSums(子集): 中位 %.1f,变异系数 %.3f\n", median(cs_cr), sd(cs_cr)/mean(cs_cr)))
cat(sprintf("小数layer  colSums:        中位 %.2f,变异系数 %.3f\n\n", median(cs_cb), sd(cs_cb)/mean(cs_cb)))

icb <- match(common, key_cb)
cb_vals <- cb_tri$x[icb]

# 给定 target 与变换,返回 cr 归一化后在 common 位置的值
cr_normed_at_common <- function(target, transform) {
  scale_vec <- ifelse(cs_cr > 0, target / cs_cr, 0)
  m <- cr %*% Diagonal(x = scale_vec)
  m@x <- transform(m@x)
  m <- as(m, "CsparseMatrix")
  mt <- summary(m)
  key_m <- paste(mt$i, mt$j)
  mt$x[match(common, key_m)]
}

candidates <- list(
  "线性 CP10K (target=1e4, 不log)"  = list(target=1e4,             fn=function(x) x),
  "log1p CP10K (target=1e4)"        = list(target=1e4,             fn=function(x) log1p(x)),
  "log1p (target=全量median)"       = list(target=med_target_full, fn=function(x) log1p(x))
)

cat(sprintf("%-36s %8s %12s\n", "候选", "相关r", "中位绝对差"))
cat(strrep("-", 58), "\n")
best_name <- ""; best_r <- -1; best_md <- NA
for (nm in names(candidates)) {
  vals <- cr_normed_at_common(candidates[[nm]]$target, candidates[[nm]]$fn)
  r  <- cor(cb_vals, vals)
  md <- median(abs(cb_vals - vals))
  cat(sprintf("%-36s %8.4f %12.4f\n", nm, r, md))
  if (r > best_r) { best_r <- r; best_md <- md; best_name <- nm }
}

## ---- 反解缩放因子(免疫于 target 猜错)------------------------------------
cat("\n========== 反解每细胞缩放因子 ==========\n")
# 取一个有代表性的细胞(非零基因较多者),反推它的缩放因子
# 线性假设:cb = k_j * cr  → k_j = cb/cr 在该细胞所有基因上应为常数
# log 假设: cb = log1p(k_j * cr) → (exp(cb)-1)/cr = k_j 应为常数
pick_cell <- which.max(Matrix::colSums(cr > 0))  # 非零基因最多的细胞
rows_lin <- cb[, pick_cell] / cr[, pick_cell]              # 线性下应为常数
rows_log <- (expm1(cb[, pick_cell])) / cr[, pick_cell]     # log 下应为常数
rows_lin <- rows_lin[is.finite(rows_lin) & cr[, pick_cell] > 0]
rows_log <- rows_log[is.finite(rows_log) & cr[, pick_cell] > 0]

cat(sprintf("选细胞 #%d(非零基因 %d 个),它的 cellranger 总量 = %.0f\n",
            pick_cell, sum(cr[, pick_cell] > 0), cs_cr[pick_cell]))
cat(sprintf("  线性假设 cb/cr:        变异系数 %.4f,中位 %.6f\n",
            sd(rows_lin)/mean(rows_lin), median(rows_lin)))
cat(sprintf("  log假设 (exp(cb)-1)/cr: 变异系数 %.4f,中位 %.6f\n",
            sd(rows_log)/mean(rows_log), median(rows_log)))
cat("判读:哪个假设下变异系数≈0(常数),就是哪种归一化。\n")
# 反推隐含 target
if (sd(rows_log)/mean(rows_log) < sd(rows_lin)/mean(rows_lin)) {
  implied_target <- median(rows_log) * cs_cr[pick_cell]
  cat(sprintf("  → log 形态更像;隐含 target ≈ %.1f(对比 1e4 / 全量median %.1f)\n",
              implied_target, med_target_full))
} else {
  implied_target <- median(rows_lin) * cs_cr[pick_cell]
  cat(sprintf("  → 线性形态更像;隐含 target ≈ %.1f(对比 1e4 / 全量median %.1f)\n",
              implied_target, med_target_full))
}

## ===========================================================================
## 自动结论(修正后的 branch 逻辑)
## ===========================================================================
cat("\n========== 自动结论 ==========\n")
exact_match <- (best_r > 0.999 && best_md < 0.01)
cat(sprintf("最佳候选: %s(r=%.4f, 中位差=%.4f)\n", best_name, best_r, best_md))

if (struct_match && exact_match) {
  # branch 1
  cat("\n>> 【钉死·branch1】小数layer = ", best_name, " of CellRanger\n", sep="")
  cat("   结构一致 + 精确吻合 → 它是归一化的 CellRanger,与 CellBender 无关。\n")
  cat("   'raw_count_cellbender' 名字误导;h5ad 无可用 CellBender 衍生数据。\n")
  cat("   决定:不能当 counts 喂 CellChat;之前用 raw_count_cellranger 是对的;去污染走 DecontX。\n")
} else if (struct_match && !exact_match) {
  # branch 2(已修正:结构一致 = 无 CellBender 成分)
  cat("\n>> 【branch2】结构一致,但三候选不精确吻合(最佳 r=", round(best_r,4), ")\n", sep="")
  cat("   重要:结构一致已在逻辑上排除 normalize(CellBender)——\n")
  cat("        归一化不改稀疏结构,而 CellBender 去污染会把 ambient 项削成0(减少非零)。\n")
  cat("   所以它是 CellRanger 的【某种我没测到的归一化】(full-data median/log2/sqrt/其他target),\n")
  cat("   仍与独立 CellBender 无关、仍不能当 counts。\n")
  cat("   决定:去污染走 DecontX。(看上面'反解'结果可进一步定型;需要可把 r 值贴我)\n")
} else {
  # branch 3(结构不一致:这里才可能有 CellBender 成分)
  cat("\n>> 【branch3】结构不一致 → 按方向细分:\n")
  if (n_cr_only > 0 && n_cb_only <= 0.001 * length(key_cr)) {
    cat("   n_cr_only>0 且 n_cb_only≈0 → cb 很可能源自 CellBender 去污染counts\n")
    cat("   (被清零的正是 ambient 项)。这是'文件确含 CellBender 成分、只是被归一化抹成不可逆'的实锤。\n")
    cat("   对写 Data availability 有用,但仍不能当 counts 喂 CellChat(归一化不可逆)。\n")
  } else if (n_cb_only > 0) {
    cat("   n_cb_only>0(cb 有 cr 没有的非零)→ 归一化和去污染都造不出新非零,\n")
    cat("   来源更怪(imputation/另一个矩阵),需细查。\n")
  }
  cat("   决定:无论如何去污染走 DecontX。把输出贴我细判。\n")
}

## ---- DecontX 可行性 -------------------------------------------------------
cat("\n========== DecontX 可行性 ==========\n")
has_celda <- requireNamespace("celda", quietly = TRUE)
cat(sprintf("celda 包已安装: %s\n", has_celda))
if (has_celda) cat(sprintf("  版本: %s\n", as.character(packageVersion("celda")))) else
  cat("  未安装,下一步:BiocManager::install('celda')\n")
cat("\nDecontX 真跑设置:输入 raw_count_cellranger;z=9个细胞类型标签;batch=biosample_id(8池);seed 固定\n")

cat("\n========== sessionInfo ==========\n")
print(sessionInfo())
