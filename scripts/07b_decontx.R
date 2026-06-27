###############################################################################
# Step 2b-1(内存安全版):DecontX de-novo 去污染
#
# 修正:07 版用 as.matrix(counts) 把 4.7万×3.6万 转稠密(~13GB)爆了 16GB 上限。
#       本版改为传 SingleCellExperiment 对象给 decontX,内部走稀疏,不实化稠密。
#
# 专家设置不变:
#   - 输入 raw_count_cellranger(原始整数 counts)
#   - z = 9 个细胞类型标签(不重聚类)
#   - batch = biosample_id(8 池,ambient profile 逐池估计)
#   - seed = 1 固定
###############################################################################

suppressPackageStartupMessages({
  library(zellkonverter)
  library(SingleCellExperiment)
  library(Matrix)
  library(celda)
})

set.seed(1)
h5      <- "/Users/franxy/Documents/博士课题 6.15/严谨/data/Hahn 人 HFpEF 单核 snRNA/HFpEF_snRNAseq.h5ad"
out_dir <- "/Users/franxy/Documents/博士课题 6.15/严谨/results"

KEEP <- c(
  "CL:0000746" = "Cardiomyocyte", "CL:0002548" = "Fibroblast",
  "CL:0010008" = "Endothelial1",  "CL:4033076" = "Endothelial2",
  "CL:0000669" = "Pericyte",      "CL:0000359" = "VSMC",
  "CL:0000235" = "Macrophage",    "CL:0000542" = "Lymphocyte",
  "CL:0002350" = "Endocardial"
)

message("== 读入 h5ad…")
sce0 <- readH5AD(h5, use_hdf5 = TRUE, reader = "R")

idx <- which(as.character(sce0$cell_type) %in% names(KEEP))
sce0 <- sce0[, idx]
celltype <- factor(KEEP[as.character(sce0$cell_type)], levels = unname(KEEP))
batch    <- factor(sce0$biosample_id)
group    <- ifelse(sce0$disease == "MONDO_0005252", "HFpEF", "control")

cat("保留后维度:", dim(sce0)[1], "基因 ×", dim(sce0)[2], "细胞\n")
cat("细胞类型分布:\n"); print(table(celltype))

# 取原始整数 counts 为稀疏矩阵(genes × cells),不转稠密
counts <- as(assay(sce0, "raw_count_cellranger"), "CsparseMatrix")
gene_names <- rownames(sce0)
cell_names <- colnames(sce0)

# 构造 SCE 时:先清掉 assay 矩阵的 dimnames,再由 SCE 统一赋名
# (否则 assay 自带行列名与 SCE 命名规则冲突,报 rownames/colnames must be NULL or identical)
dimnames(counts) <- NULL
sce_in <- SingleCellExperiment(
  assays  = list(counts = counts),
  colData = DataFrame(celltype = celltype, batch = batch, group = group)
)
rownames(sce_in) <- gene_names
colnames(sce_in) <- cell_names

## ---- 运行 DecontX(传 SCE 对象)-------------------------------------------
message("\n== 运行 DecontX(SCE 输入, z=细胞类型, batch=池, seed=1)…")
t0 <- Sys.time()
sce_dx <- decontX(
  x     = sce_in,
  z     = celltype,
  batch = batch,
  seed  = 1
)
cat(sprintf("DecontX 完成,用时 %.1f 分钟\n", as.numeric(difftime(Sys.time(), t0, units="mins"))))

## ---- 提取结果 ------------------------------------------------------------
# decontX 把去污染 counts 放进 assay 'decontXcounts',污染比例放 colData$decontX_contamination
decontaminated <- assay(sce_dx, "decontXcounts")          # genes × cells
decontaminated <- as(decontaminated, "CsparseMatrix")
contamination  <- sce_dx$decontX_contamination

cat("\n== 去污染结果概览 ==\n")
cat(sprintf("每细胞污染比例: 中位 %.3f,范围 [%.3f, %.3f]\n",
            median(contamination), min(contamination), max(contamination)))
cat("按细胞类型的中位污染比例:\n")
print(round(tapply(contamination, celltype, median), 3))

# 去污染 counts 为小数(期望值),取整后给 CellChat
decont_int <- decontaminated
decont_int@x <- round(decont_int@x)
decont_int <- drop0(decont_int)

## ---- 保存 ----------------------------------------------------------------
saveRDS(list(
  counts_decont     = decont_int,
  counts_decont_raw = decontaminated,
  contamination     = contamination,
  celltype          = celltype,
  batch             = batch,
  group             = group,
  cellnames         = cell_names
), file.path(out_dir, "decontx_counts.rds"))

cat(sprintf("\n== 完成。已存 %s\n", file.path(out_dir, "decontx_counts.rds")))
cat("下一步:Step 2b-2 三臂 CellChat(raw / CellBender layer / DecontX),参数全锁死。\n")

cat("\n========== sessionInfo ==========\n")
print(sessionInfo())
