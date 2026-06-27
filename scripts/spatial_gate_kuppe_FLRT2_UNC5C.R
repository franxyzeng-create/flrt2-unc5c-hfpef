###############################################################################
# spatial_gate_kuppe_FLRT2_UNC5C.R
# RStudio 本地复算:FLRT2 / UNC5C 在 Kuppe Visium fibroblast spot 的检出率
#
# 目的:可行性 gate —— 决定"公开 Visium 空间共定位"走不走得通(低丰度 dropout)。
# 数据:Kuppe 2022 CELLxGENE Visium 4 个 .h5ad(control_P1/P7、FZ_P18/P20)。
# 环境:R 4.6 / Bioconductor 3.23;读 .h5ad 用 zellkonverter。
# 注:首次跑 zellkonverter 会用 basilisk 自动配一个 python 环境(几分钟,只一次)。
#
# 用法:改下面 DATA_DIR 为你的 Kuppe 文件夹 → 全选 source。
# 预期(应与之前沙箱结果一致):
#   P1(normal) UMI~3264  FLRT2 3.4%  UNC5C 5.1%  共检出 0.0%
#   P7(normal) UMI~2821  FLRT2 2.3%  UNC5C 1.9%  共检出 0.3%
#   P18(MI)    UMI~1365  FLRT2 7.8%  UNC5C 3.0%  共检出 0.4%
#   P20(MI)    UMI~1349  FLRT2 6.6%  UNC5C 3.3%  共检出 0.2%
###############################################################################

suppressPackageStartupMessages({
  library(zellkonverter); library(SingleCellExperiment); library(Matrix)
})

## ====== 自动定位 Kuppe 文件夹(无需手改,直接运行) ======
ROOT <- "/Users/franxy/Documents/博士课题 6.15/严谨/data"
.dirs    <- list.dirs(ROOT, recursive = FALSE)
DATA_DIR <- .dirs[grepl("Kuppe", basename(.dirs))][1]
stopifnot(!is.na(DATA_DIR), dir.exists(DATA_DIR))
cat("使用 Kuppe 文件夹:", DATA_DIR, "\n")

## cell2location 反卷积的细胞类型丰度列(obs 里现成)
CT_ALL <- c("Adipocyte","Cardiomyocyte","Endothelial","Fibroblast","Lymphoid",
            "Mast","Myeloid","Neuronal","Pericyte","Cycling.cells","vSMCs")

## ---- 辅助函数(先定义后用;标量条件用 if/else,不用 ifelse) ----
gene_row <- function(sce, sym) {                 # 按 symbol 找行号(feature_name 列)
  fn <- as.character(rowData(sce)$feature_name)
  w  <- which(fn == sym)
  if (length(w) == 0) w <- which(rownames(sce) == sym)
  w
}
det_rate <- function(vec, mask) {                # 检出率 = % spot 表达>0
  if (sum(mask, na.rm = TRUE) == 0) return(NA_real_)
  100 * mean(vec[mask] > 0, na.rm = TRUE)
}

## ---- 主循环:文件夹里所有 .h5ad ----
h5s <- list.files(DATA_DIR, pattern = "\\.h5ad$", full.names = TRUE)
stopifnot(length(h5s) > 0)

res <- data.frame()
for (fp in h5s) {
  sce <- readH5AD(fp, verbose = FALSE)
  cd  <- colData(sce)

  ## 自动标签:donor + disease(不依赖文件名)
  donor <- if ("donor_id" %in% colnames(cd)) as.character(cd$donor_id[1]) else "?"
  dis   <- if ("disease"  %in% colnames(cd)) as.character(cd$disease[1])  else "?"
  lab   <- paste0(donor, " (", dis, ")")

  ## X(normalized);检出率看 >0,与 counts>0 等价
  X <- assay(sce, 1)

  ## cell2location 丰度 → fibroblast-dominant spot(fibroblast 为最高丰度的 spot)
  CT  <- intersect(CT_ALL, colnames(cd))
  M   <- as.matrix(cd[, CT]); storage.mode(M) <- "numeric"
  fib <- M[, "Fibroblast"]
  valid <- !is.na(fib)                           # 部分 spot 没反卷积 → NaN
  Mz  <- M; Mz[is.na(Mz)] <- -1
  fib_dom <- valid & (max.col(Mz, ties.method = "first") == which(CT == "Fibroblast"))

  ## FLRT2 / UNC5C 表达向量
  iF <- gene_row(sce, "FLRT2"); iU <- gene_row(sce, "UNC5C")
  xF <- as.numeric(X[iF[1], ]); xU <- as.numeric(X[iU[1], ])
  co <- (xF > 0) & (xU > 0)                       # 同一 spot 两个都检出

  med_umi <- if ("n_counts" %in% colnames(cd)) median(cd$n_counts, na.rm = TRUE) else NA_real_

  res <- rbind(res, data.frame(
    slice         = lab,
    median_UMI    = round(med_umi),
    fib_dom_spots = sum(fib_dom),
    FLRT2_pct     = round(det_rate(xF, fib_dom), 1),
    UNC5C_pct     = round(det_rate(xU, fib_dom), 1),
    codetect_pct  = round(det_rate(co, fib_dom), 1)
  ))
  rm(sce, X); gc()
}

cat("\n== FLRT2/UNC5C 在 fibroblast-dominant spot 的检出率(% spot expr>0) ==\n")
print(res, row.names = FALSE)
cat("\n判读:codetect_pct(共检出)< 1% → 空间共定位作图不可行(低丰度 dropout)→ 空间留白。\n")
cat("       FLRT2 在 MI(FZ)fib spot 检出 > 正常 = 弱正面趋势(疾病升高)。\n")
