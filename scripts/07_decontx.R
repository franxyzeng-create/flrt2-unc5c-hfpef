###############################################################################
# Step 2b-1:DecontX de-novo 去污染(为三臂稳健性的 (b) 臂准备输入)
#
# 专家拍板:(c) 三臂稳健性,(a) CellBender layer 为主,(b) DecontX 为方法独立交叉验证。
# DecontX 价值:cluster-based、不用空液滴,模型假设与 CellBender 完全不重叠。
#
# 关键设置(专家记下):
#   - 输入:raw_count_cellranger(原始整数 counts)
#   - z = 9 个已发表细胞类型标签(不让它自己重聚类 → 可复现 + 用原文生物学)
#   - batch = biosample_id(8 个池;ambient profile 逐池不同,分池估计更准)
#   - seed 固定(DecontX 是变分推断,带随机初始化)
#
# 产出:decontx_counts.rds(去污染后的 counts 矩阵 + 每细胞污染比例),
#       供 Step 2b-2 三臂 CellChat 使用,不用重跑 DecontX。
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

# 9 个可分析细胞类型(与 CellChat 探查/出版版完全一致)
KEEP <- c(
  "CL:0000746" = "Cardiomyocyte", "CL:0002548" = "Fibroblast",
  "CL:0010008" = "Endothelial1",  "CL:4033076" = "Endothelial2",
  "CL:0000669" = "Pericyte",      "CL:0000359" = "VSMC",
  "CL:0000235" = "Macrophage",    "CL:0000542" = "Lymphocyte",
  "CL:0002350" = "Endocardial"
)

message("== 读入 h5ad…")
sce <- readH5AD(h5, use_hdf5 = TRUE, reader = "R")

# 只保留 9 个细胞类型
idx <- which(as.character(sce$cell_type) %in% names(KEEP))
sce <- sce[, idx]
celltype <- factor(KEEP[as.character(sce$cell_type)], levels = unname(KEEP))
batch    <- factor(sce$biosample_id)   # 8 个池
group    <- ifelse(sce$disease == "MONDO_0005252", "HFpEF", "control")

cat("保留后维度:", dim(sce)[1], "基因 ×", dim(sce)[2], "细胞\n")
cat("细胞类型分布:\n"); print(table(celltype))
cat("池(batch)分布:\n"); print(table(batch))

# 取原始整数 counts(genes × cells)
counts <- as(assay(sce, "raw_count_cellranger"), "CsparseMatrix")
rownames(counts) <- rownames(sce)
colnames(counts) <- colnames(sce)

## ---- 运行 DecontX --------------------------------------------------------
# z = 细胞类型标签(已知,不重聚类);batch = 池;seed 固定
message("\n== 运行 DecontX(z=细胞类型, batch=池, seed=1)…这步稍慢,请耐心")
t0 <- Sys.time()
dx <- decontX(
  x     = as.matrix(counts),    # DecontX 接受 matrix 或 SCE;大矩阵用稀疏更省内存,但 decontX 内部会处理
  z     = celltype,             # 已知细胞类型,不自己聚类
  batch = batch,                # 分池估计 ambient profile
  seed  = 1
)
cat(sprintf("DecontX 完成,用时 %.1f 分钟\n", as.numeric(difftime(Sys.time(), t0, units="mins"))))

## ---- 提取去污染 counts + 污染比例 ----------------------------------------
decontaminated <- dx$decontXcounts        # 去污染后的 counts(genes × cells)
contamination  <- dx$contamination        # 每细胞的污染比例估计

cat("\n== 去污染结果概览 ==\n")
cat(sprintf("每细胞污染比例: 中位 %.3f,范围 [%.3f, %.3f]\n",
            median(contamination), min(contamination), max(contamination)))
cat("按细胞类型的中位污染比例:\n")
print(round(tapply(contamination, celltype, median), 3))
cat("\n(预期:高丰度细胞如心肌/成纤维贡献 ambient,各类型污染比例可能不同)\n")

# 去污染 counts 可能是小数(期望值),CellChat normalizeData 能处理;但为稳妥取整
decont_int <- decontaminated
decont_int@x <- round(decont_int@x)
decont_int <- drop0(decont_int)   # 去掉变成 0 的

## ---- 保存(供 Step 2b-2 三臂 CellChat 用)---------------------------------
saveRDS(list(
  counts_decont     = decont_int,        # 去污染后取整 counts(给 CellChat)
  counts_decont_raw = decontaminated,    # 去污染原始(小数)counts
  contamination     = contamination,
  celltype          = celltype,
  batch             = batch,
  group             = group,
  cellnames         = colnames(counts)
), file.path(out_dir, "decontx_counts.rds"))

cat(sprintf("\n== 完成。已存 %s\n", file.path(out_dir, "decontx_counts.rds")))
cat("下一步:Step 2b-2 三臂 CellChat(raw / CellBender layer / DecontX),参数全锁死。\n")

cat("\n========== sessionInfo ==========\n")
print(sessionInfo())
