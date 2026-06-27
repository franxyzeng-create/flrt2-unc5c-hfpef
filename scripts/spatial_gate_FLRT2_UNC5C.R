###############################################################################
# spatial_gate_FLRT2_UNC5C.R
# 可行性 gate:查 FLRT2 / UNC5C 在 fibroblast-rich Visium spot 的检出率
#
# 目的:决定"公开 Visium 空间共定位"这条路走不走得通(低丰度 dropout?)。
#       这是纯生信、零湿实验花费;FLRT2/UNC5C 是低丰度 guidance 基因,
#       Visium 可能 dropout,所以先用一个切片做 gate,够了再正式做共定位图。
#
# 环境:R 4.6 / Bioconductor 3.23;.h5ad 用 zellkonverter 读(你装过)。
# 数据:CELLxGENE / heartcellatlas 下载的 Kanemaru(正常)/ Kuppe(MI) Visium .h5ad
# 用法:改下面 H5AD 路径 → source 本脚本 → 先看 Part1 结构,再看 Part3 检出率。
#
# 判读(Part3 末):fib-rich 检出率 >~15-20% 才有作图价值;<~5% = 留白。
#
# R 易错点遵循:标量条件用 if/else(不用 ifelse);辅助函数先定义后用。
###############################################################################

suppressPackageStartupMessages({
  library(zellkonverter); library(SingleCellExperiment); library(Matrix)
})

## ====== 改这里:你下载的 .h5ad 路径(先用一个切片) ======
H5AD <- "/Users/franxy/Downloads/your_visium.h5ad"

## use_hdf5=TRUE → 不把矩阵全载进内存(M1 16GB 友好)
sce <- readH5AD(H5AD, use_hdf5 = TRUE)

## ---- 辅助函数(先定义后用) ----
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# 基因 symbol 向量(优先 rowData$feature_name,否则 rownames)
gene_symbols <- function(sce) {
  rd <- rowData(sce)
  if ("feature_name" %in% colnames(rd)) as.character(rd$feature_name) else rownames(sce)
}
# 按 symbol 或 Ensembl ID 找行号(兼容两种 rownames)
find_gene <- function(sce, sym, ens) {
  gs  <- gene_symbols(sce)
  hit <- which(gs == sym)
  if (length(hit) == 0) hit <- which(rownames(sce) == ens)
  if (length(hit) == 0) hit <- which(rownames(sce) == sym)
  hit
}

## ====== Part 1:探查结构(先看清楚再分析) ======
cat("========== 结构探查 ==========\n")
cat("dim (genes x spots):", nrow(sce), "x", ncol(sce), "\n")
cat("assayNames:", paste(assayNames(sce), collapse = ", "), "\n")
cat("\ncolData(obs) 列名:\n"); print(colnames(colData(sce)))
cat("\nrowData(var) 列名:\n"); print(colnames(rowData(sce)))
cat("\nrownames 头 3:\n"); print(head(rownames(sce), 3))

iF <- find_gene(sce, "FLRT2", "ENSG00000185070")
iU <- find_gene(sce, "UNC5C", "ENSG00000182168")
cat("\nFLRT2 行号:", iF %||% NA, "   UNC5C 行号:", iU %||% NA, "\n")
if (length(iF) == 0 || length(iU) == 0) {
  cat("[!] 有基因没找到 → 看上面 rowData 列名/rownames,确认是 symbol 还是 Ensembl,必要时改 find_gene\n")
}

## 选检出率用的 assay(检出率=是否>0,normalized 或 counts 都等价)
counts_assay <- if ("counts" %in% assayNames(sce)) "counts" else assayNames(sce)[1]
A <- assay(sce, counts_assay)
cat("\n检出率使用 assay:", counts_assay, "\n")

## ====== Part 2:定 fibroblast-rich spot ======
gs <- gene_symbols(sce)

## (A) 先找 obs 里现成的 fibroblast 注释/deconvolution 列(列名含 fibro 且为数值)
fib_cols <- grep("fibro", colnames(colData(sce)), ignore.case = TRUE, value = TRUE)
cat("\nobs 里疑似 fibroblast 列:", paste(fib_cols, collapse = ", ") %||% "(无)", "\n")

use_decon <- FALSE
if (length(fib_cols) >= 1) {
  if (is.numeric(colData(sce)[[fib_cols[1]]])) use_decon <- TRUE
}

if (use_decon) {
  fp <- colData(sce)[[fib_cols[1]]]
  fib_rich <- fp >= quantile(fp, 0.75, na.rm = TRUE)
  cat("→ 用 deconvolution 列:", fib_cols[1], "(top 25% 作 fibroblast-rich)\n")
} else {
  ## (B) 兜底:fibroblast marker score(spot 上 marker 总表达,取 top 25%)
  fib_markers <- c("PDGFRA","DCN","COL1A1","COL1A2","GSN","LUM","FBLN1")
  fib_idx <- which(gs %in% fib_markers)
  cat("→ 用 marker score 兜底,找到 markers:", paste(gs[fib_idx], collapse = ", "), "\n")
  fib_score <- Matrix::colSums(A[fib_idx, , drop = FALSE])
  fib_rich  <- fib_score >= quantile(fib_score, 0.75, na.rm = TRUE)
}
cat("fibroblast-rich spot:", sum(fib_rich, na.rm = TRUE), "/", length(fib_rich), "\n")

## ====== Part 3:FLRT2 / UNC5C 检出率 ======
detrate <- function(A, gi, cells) {
  if (length(gi) == 0) return(NA_real_)
  100 * mean(A[gi[1], cells, drop = FALSE] > 0, na.rm = TRUE)
}

cat("\n========== 检出率 (% spot with expr > 0) ==========\n")
cat(sprintf("FLRT2   fib-rich %5.1f%%  |  其他 %5.1f%%\n", detrate(A, iF, fib_rich), detrate(A, iF, !fib_rich)))
cat(sprintf("UNC5C   fib-rich %5.1f%%  |  其他 %5.1f%%\n", detrate(A, iU, fib_rich), detrate(A, iU, !fib_rich)))
if (length(iF) > 0 && length(iU) > 0) {
  co <- (A[iF[1], ] > 0) & (A[iU[1], ] > 0)
  cat(sprintf("共检出  fib-rich %5.1f%%  |  其他 %5.1f%%\n",
              100 * mean(co[fib_rich], na.rm = TRUE), 100 * mean(co[!fib_rich], na.rm = TRUE)))
}

cat("\n========== gate 判读 ==========\n")
cat("fib-rich 里检出率 >~15-20% → 空间共定位可作图(进正式分析);\n")
cat("~5-15%             → 勉强,可考虑合并多切片 / 邻域平滑 / 只做疾病(Kuppe)区;\n")
cat("<~5%               → dropout 太重,空间这块诚实留白,只用 HPA/GTEx(FLRT2)+ snRNA。\n")
cat("\n提示:UNC5C 正常组织本就低(见 HPA/GTEx),Kanemaru(正常)里 UNC5C 低是预期;\n")
cat("     重点看 Kuppe(疾病/纤维化区)里 UNC5C 在 fib-rich spot 是否抬头(疾病诱导)。\n")
