###############################################################################
# Step 2a-2(重写版):小数 layer 身份终极判别 —— 自动钉死,不用人肉贴值
#
# 按代码审查意见重写:
#   - 一次性测三种归一化候选:线性CP10K / log1p@1e4 / log1p@median
#   - 取值范围 range(cb@x):几乎单独定生死(max≈9+colSums各异→log;max几百~万+colSums≈常数→线性)
#   - 结构断言:normalize 不改稀疏结构,cb 与 cr 非零位置必须逐一相同
#   - 方向断言:防 reader="R" 把 features×cells 方向搞反导致 colSums 静默全错
#   - 顺带解悖论:cb==normalize(cellranger)? → 判断这层到底有没有 CellBender 成分
#   - 随机抽 3000 细胞(防 h5ad 按 donor 排序导致 subset 有偏)
#
# 背景:script-05 已测到 raw_count_cellbender 是小数、22.9% 位置 cb>cr(最大2.88),
#       不符合 CellBender"只减不增" → 它不是 CellBender 原始counts,是归一化/变换产物。
#       本脚本钉死它到底是哪种归一化、以及有没有 CellBender 成分。
###############################################################################

suppressPackageStartupMessages({
  library(zellkonverter)
  library(SingleCellExperiment)
  library(Matrix)
})

set.seed(1)
h5 <- "/Users/franxy/Documents/博士课题 6.15/严谨/data/Hahn 人 HFpEF 单核 snRNA/HFpEF_snRNAseq.h5ad"

message("== 读入 h5ad…")
sce <- readH5AD(h5, use_hdf5 = TRUE, reader = "R")

# 随机抽 3000 细胞(防按 donor 排序的偏倚)
set.seed(1)
sel <- sort(sample(ncol(sce), 3000))
sub <- sce[, sel]
cb <- as(assay(sub, "raw_count_cellbender"), "CsparseMatrix")  # 小数 layer
cr <- as(assay(sub, "raw_count_cellranger"), "CsparseMatrix")  # 整数 counts

## ===========================================================================
## 方向断言(最该防的静默错误):SCE 应是 features × cells
## ===========================================================================
cat("========== 方向断言 ==========\n")
cat(sprintf("dim(cb) = %d × %d(应为 ~36601 基因 × 3000 细胞)\n", nrow(cb), ncol(cb)))
stopifnot(nrow(cb) > ncol(cb))   # 基因数 > 细胞数,否则方向翻了,后面 colSums 全错
cat(">> 方向正确:行=基因,列=细胞,colSums = 每细胞总量 ✓\n")

## ===========================================================================
## 关键事实 0:小数 layer 的取值范围(几乎单独定生死)
## ===========================================================================
cat("\n========== 取值范围(决定性线索)==========\n")
cat(sprintf("小数layer 非零值范围: [%.4f, %.4f]\n", min(cb@x), max(cb@x)))
cat(sprintf("cellranger 非零值范围: [%.0f, %.0f]\n", min(cr@x), max(cr@x)))
cat("判读:max≈9 上下 → log1p 形态;max 到几百~上万 → 线性 CP10K(未log)\n")

## ===========================================================================
## 结构断言:normalize 不改变稀疏结构 → cb 与 cr 非零位置必须逐一相同
## ===========================================================================
cat("\n========== 结构断言(非零位置是否一致)==========\n")
cb_tri <- summary(cb); cr_tri <- summary(cr)
key_cb <- paste(cb_tri$i, cb_tri$j)
key_cr <- paste(cr_tri$i, cr_tri$j)
cat(sprintf("cb 非零位置数: %d\n", length(key_cb)))
cat(sprintf("cr 非零位置数: %d\n", length(key_cr)))
common <- intersect(key_cb, key_cr)
cat(sprintf("交集位置数:   %d\n", length(common)))
n_cb_only <- length(setdiff(key_cb, key_cr))
n_cr_only <- length(setdiff(key_cr, key_cb))
cat(sprintf("cb 独有(cr为0): %d   cr 独有(cb为0): %d\n", n_cb_only, n_cr_only))

struct_match <- (length(key_cb) == length(key_cr)) && (n_cb_only == 0) && (n_cr_only == 0)
cat("\n--- 结构判读 ---\n")
if (struct_match) {
  cat(">> 非零结构完全一致 → 符合'cb 是 cr 的某种逐元素归一化'(不改稀疏结构)✓\n")
} else {
  cat(">> 非零结构【不一致】→ cb 不是 cr 的简单归一化!\n")
  cat("   (cb 独有或 cr 独有位置非零,说明两者非零集不同,排除 CP10K/log1p(cr))\n")
  cat("   → cb 可能是 normalize(CellBender) 等其他来源,或经过改变结构的处理\n")
}

## ===========================================================================
## 一次性测三种归一化候选:cb 是否 = 某种 normalize(cr)
## ===========================================================================
cat("\n========== 三种归一化候选比对 ==========\n")
cs_cr <- Matrix::colSums(cr)
cs_cb <- Matrix::colSums(cb)
cat(sprintf("cellranger colSums: 中位 %.1f,变异系数 %.3f\n", median(cs_cr), sd(cs_cr)/mean(cs_cr)))
cat(sprintf("小数layer  colSums: 中位 %.2f,变异系数 %.3f\n", median(cs_cb), sd(cs_cb)/mean(cs_cb)))
cat("(若小数layer colSums 变异系数≈0 → 被归一化到常数;线性CP10K会这样)\n\n")

# 防空细胞 Inf
safe_scale <- function(target) ifelse(cs_cr > 0, target / cs_cr, 0)

# 在交集位置上比较(结构一致时交集=全部)
icb <- match(common, key_cb)
get_cr_normed <- function(transform) {
  # transform: 函数,作用于按列缩放后的 cr
  m <- cr %*% Diagonal(x = safe_scale_target)
  m@x <- transform(m@x)
  mt <- summary(as(m, "CsparseMatrix"))
  key_m <- paste(mt$i, mt$j)
  im <- match(common, key_m)
  mt$x[im]
}

cb_vals <- cb_tri$x[icb]
med_target <- median(cs_cr)

candidates <- list(
  "线性 CP10K (target=1e4, 不log)"      = list(target=1e4,        fn=function(x) x),
  "log1p CP10K (target=1e4)"            = list(target=1e4,        fn=function(x) log1p(x)),
  "log1p (target=median counts)"        = list(target=med_target, fn=function(x) log1p(x))
)

cat(sprintf("%-38s %8s %10s\n", "候选归一化", "相关r", "中位绝对差"))
cat(strrep("-", 60), "\n")
best_name <- ""; best_r <- -1
for (nm in names(candidates)) {
  tg <- candidates[[nm]]$target
  fn <- candidates[[nm]]$fn
  safe_scale_target <<- safe_scale(tg)   # 供 get_cr_normed 用
  vals <- get_cr_normed(fn)
  r <- cor(cb_vals, vals)
  md <- median(abs(cb_vals - vals))
  cat(sprintf("%-38s %8.4f %10.4f\n", nm, r, md))
  if (r > best_r) { best_r <- r; best_name <- nm }
}

## ===========================================================================
## 自动结论
## ===========================================================================
cat("\n========== 自动结论 ==========\n")
# 找最佳候选并判断是否精确吻合
safe_scale_target <<- safe_scale(candidates[[best_name]]$target)
best_vals <- get_cr_normed(candidates[[best_name]]$fn)
best_md <- median(abs(cb_vals - best_vals))
exact_match <- (best_r > 0.999 && best_md < 0.01)

cat(sprintf("最佳匹配候选: %s(r=%.4f)\n", best_name, best_r))

if (struct_match && exact_match) {
  cat("\n>> 【钉死】小数layer = ", best_name, "of CellRanger counts\n", sep="")
  cat("   含义:这层只是归一化的 CellRanger,与 CellBender 无关,\n")
  cat("        'raw_count_cellbender' 这个名字是误导;h5ad 里没有可用的 CellBender 衍生数据。\n")
  cat("   决定:① 它不能当 counts 喂 CellChat(它是归一化数据);\n")
  cat("        ② 之前 CellChat 用 raw_count_cellranger 当输入是对的;\n")
  cat("        ③ 去污染只能走 de-novo DecontX(下一步)。\n")
} else if (struct_match && !exact_match) {
  cat("\n>> 结构一致但三种 normalize(cellranger) 都不精确吻合(最佳 r=", round(best_r,4), ")\n", sep="")
  cat("   含义:它可能是 normalize(CellBender counts)——即去污染counts垫在底下、\n")
  cat("        被归一化抹掉、不可逆。或用了其他 target/方法。\n")
  cat("   决定:仍不能当 counts 喂 CellChat;去污染走 DecontX。\n")
  cat("   (把上面三行 r 值贴给我,可进一步判断是否含 CellBender 成分)\n")
} else {
  cat("\n>> 非零结构不一致 → 它不是 cr 的简单归一化,来源更复杂。\n")
  cat("   决定:不能当 counts 喂 CellChat;去污染走 DecontX。把输出贴我细判。\n")
}

## ===========================================================================
## DecontX 可行性(探包)
## ===========================================================================
cat("\n========== DecontX 可行性 ==========\n")
has_celda <- requireNamespace("celda", quietly = TRUE)
cat(sprintf("celda 包(含 DecontX)已安装: %s\n", has_celda))
if (has_celda) {
  cat(sprintf("  版本: %s → 可直接跑 DecontX\n", as.character(packageVersion("celda"))))
} else {
  cat("  未安装。下一步:BiocManager::install('celda')\n")
}
cat("\nDecontX 真跑时的关键设置(记下):\n")
cat("  - 输入:raw_count_cellranger(原始整数)\n")
cat("  - z = 9 个已发表细胞类型标签(不让它自己重聚类,可复现+用原文生物学)\n")
cat("  - batch = biosample_id(8个池,ambient profile 逐池不同,分池估计更准)\n")
cat("  - seed 固定(DecontX 是变分推断,带随机初始化)\n")

## ---- 存档(可复现性纪律)----
cat("\n========== sessionInfo ==========\n")
print(sessionInfo())
