###############################################################################
# 第 1 步(探查):CellChat 全局细胞通讯 —— HFpEF vs Control
#
# 目的:看 9 个可分析细胞类型之间的细胞通讯,在 HFpEF 与对照之间
#       哪些信号通路 / 配受体对显著改变;血管龛(内皮-周细胞-VSMC-心内膜)
#       的通讯是否突出。用数据帮我们决定 β(CellChat)是否当主线。
#
# 关键设计:
#   - 输入用 raw_count_cellranger(整数),CellChat 标准 normalizeData 自己归一化
#     (不用 X,因为 X 是 CellBender 的小数 counts,不是 log-normalized)
#   - 只保留原文做 DE 的 9 个细胞类型,用可读名
#   - HFpEF / control 各建一个 CellChat 对象,同库同预处理,最后 mergeCellChat 比较
#   - 所有对象存 RDS,后续画图/深挖不重跑
#
# 预期耗时:M1 16G 上约 30-60 分钟(computeCommunProb 最慢)。
###############################################################################

suppressPackageStartupMessages({
  library(zellkonverter)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(CellChat)
  library(patchwork)
})

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 4 * 1024^3)  # 抬高并行对象阈值到 4 GiB(默认 500MiB 会报错)
future::plan("multisession", workers = 4)   # 并行加速;若仍报错改成 future::plan("sequential")

## ---- 0. 配置 --------------------------------------------------------------
h5      <- "/Users/franxy/Documents/博士课题 6.15/严谨/data/Hahn 人 HFpEF 单核 snRNA/HFpEF_snRNAseq.h5ad"
out_dir <- "/Users/franxy/Documents/博士课题 6.15/严谨/results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 保留原文做 DE 的 9 个细胞类型(CL 编码 → 可读名)
KEEP <- c(
  "CL:0000746" = "Cardiomyocyte",
  "CL:0002548" = "Fibroblast",
  "CL:0010008" = "Endothelial1",
  "CL:4033076" = "Endothelial2",      # 注:这是较小的 EC2;若与原文命名有出入下一步再核
  "CL:0000669" = "Pericyte",
  "CL:0000359" = "VSMC",
  "CL:0000235" = "Macrophage",
  "CL:0000542" = "Lymphocyte",
  "CL:0002350" = "Endocardial"
)

## ---- 1. 读入 h5ad,取 CellRanger 整数 counts ------------------------------
message("== 读入 h5ad…")
sce <- readH5AD(h5, use_hdf5 = TRUE, reader = "R")

# 只保留 9 个目标细胞类型的核
idx <- which(as.character(sce$cell_type) %in% names(KEEP))
sce <- sce[, idx]
sce$celltype <- factor(KEEP[as.character(sce$cell_type)], levels = unname(KEEP))
sce$group    <- ifelse(sce$disease == "MONDO_0005252", "HFpEF", "control")

message("== 保留后各细胞类型核数(按组):")
print(table(sce$celltype, sce$group))

# 取 CellRanger 整数 counts,实化为稀疏矩阵(genes × cells)
counts <- as(assay(sce, "raw_count_cellranger"), "CsparseMatrix")
rownames(counts) <- rownames(sce)
colnames(counts) <- colnames(sce)

## ---- 2. 构建函数:对单组建 CellChat 并跑完整流程 -------------------------
run_cellchat_one <- function(counts_sub, labels) {
  # CellChat 标准 normalizeData:library-size 归一 + log1p
  data.norm <- normalizeData(counts_sub)
  meta <- data.frame(labels = labels,
                     samples = "sample1",          # 消除 createCellChat 的 samples 列 warning
                     row.names = colnames(counts_sub))

  cc <- createCellChat(object = data.norm, meta = meta, group.by = "labels")
  cc@DB <- CellChatDB.human                       # 人配受体数据库(包自带,无需下载)

  cc <- subsetData(cc)                            # 只保留信号基因,省内存
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)

  # 最慢的一步:置换检验估通讯概率
  cc <- computeCommunProb(cc, type = "truncatedMean", trim = 0.1,
                          population.size = TRUE)
  cc <- filterCommunication(cc, min.cells = 10)   # 少于10核的细胞对不算
  cc <- computeCommunProbPathway(cc)              # 汇总到信号通路层
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc)         # 计算各细胞的发送/接收中心性
  cc
}

## ---- 3. 分别对 HFpEF 和 control 跑 ----------------------------------------
groups <- c("control", "HFpEF")
cc_list <- list()
for (g in groups) {
  message(sprintf("\n========== 跑 CellChat: %s ==========", g))
  cells_g <- colnames(sce)[sce$group == g]
  t0 <- Sys.time()
  cc_list[[g]] <- run_cellchat_one(counts[, cells_g, drop = FALSE],
                                   droplevels(sce$celltype[sce$group == g]))
  message(sprintf("   %s 完成,用时 %.1f 分钟",
                  g, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  saveRDS(cc_list[[g]], file.path(out_dir, sprintf("cellchat_%s.rds", g)))
}

## ---- 4. 合并两组,做差异比较 ----------------------------------------------
message("\n== 合并两组并比较…")
cc_list <- lapply(cc_list, function(x) { x@meta$labels <- droplevels(x@meta$labels); x })
cellchat <- mergeCellChat(cc_list, add.names = names(cc_list))
saveRDS(cellchat, file.path(out_dir, "cellchat_merged.rds"))

## ---- 5. 关键结果打印 + 出图 -----------------------------------------------
# (1) 总体:交互数量与强度,HFpEF vs control
message("\n== 总体交互数量/强度对比(HFpEF vs control):")
gg1 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1,2),
                           measure = "weight")
ggsave(file.path(out_dir, "cc_overall_interactions.png"),
       gg1 + gg2, width = 8, height = 4, dpi = 150)

# (2) 差异网络图:哪些细胞对之间通讯增强(红)/减弱(蓝)
png(file.path(out_dir, "cc_diff_network.png"), width = 1000, height = 500)
par(mfrow = c(1,2), xpd = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")
dev.off()

# (3) 各信号通路在两组中的总强度排名(最有信息量的一张)
message("\n== 各信号通路在两组中的相对强度(看哪些通路 HFpEF 偏强/偏弱):")
gg_rank <- rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE)
print(gg_rank$data[order(-abs(gg_rank$data$contribution.scaled)), ][1:30, ])
ggsave(file.path(out_dir, "cc_pathway_rank.png"), gg_rank, width = 6, height = 9, dpi = 150)

# (4) 每组各细胞类型作为"发送方/接收方"的信号强度(2D 散点)
for (g in groups) {
  p <- netAnalysis_signalingRole_scatter(cc_list[[g]], title = g)
  ggsave(file.path(out_dir, sprintf("cc_signalingRole_%s.png", g)),
         p, width = 5, height = 5, dpi = 150)
}

message(sprintf("\n== 完成。图与对象已存到 %s", out_dir))
message("   关键看:cc_pathway_rank.png(哪些通路 HFpEF 偏强/弱)")
message("            cc_diff_network.png(哪些细胞对通讯增强/减弱)")
message("            cellchat_merged.rds(后续深挖,不用重跑)")

future::plan("sequential")   # 收尾,释放并行 worker
