###############################################################################
# 第 1 步(出版级重跑):CellChat HFpEF vs Control —— 单线程 + 固定种子
#
# 与探查版的区别(全部为了可复现 + 干净 p 值):
#   - future::plan("sequential")          关闭并行,根除并行随机数问题
#   - computeCommunProb(seed.use = 1)      固定置换检验种子
#   - nboot = 100 显式写出                  置换次数可报告
#   - set.seed(1) 全程兜底
#   - 结果存为 *_final.rds,不覆盖探查版
#
# 诚实预期:单线程下 computeCommunProb 很慢,两组可能 8-12 小时。
#           建议睡前启动 + 电脑接电源 + 关闭自动睡眠(否则进程会被挂起)。
#
# 防睡眠建议:另开一个终端跑  caffeinate -i  (在脚本运行期间保持系统不休眠)
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
future::plan("sequential")          # 关键:单线程,根除并行随机数问题
set.seed(1)                         # 全局种子兜底

## ---- 0. 配置 --------------------------------------------------------------
h5      <- "/Users/franxy/Documents/博士课题 6.15/严谨/data/Hahn 人 HFpEF 单核 snRNA/HFpEF_snRNAseq.h5ad"
out_dir <- "/Users/franxy/Documents/博士课题 6.15/严谨/results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

SEED  <- 1      # computeCommunProb 的种子
NBOOT <- 100    # 置换次数(默认 100;显式写出以便方法学报告)

KEEP <- c(
  "CL:0000746" = "Cardiomyocyte",
  "CL:0002548" = "Fibroblast",
  "CL:0010008" = "Endothelial1",
  "CL:4033076" = "Endothelial2",
  "CL:0000669" = "Pericyte",
  "CL:0000359" = "VSMC",
  "CL:0000235" = "Macrophage",
  "CL:0000542" = "Lymphocyte",
  "CL:0002350" = "Endocardial"
)

## ---- 1. 读入 + 取 CellRanger 整数 counts ----------------------------------
message("== 读入 h5ad…")
sce <- readH5AD(h5, use_hdf5 = TRUE, reader = "R")
idx <- which(as.character(sce$cell_type) %in% names(KEEP))
sce <- sce[, idx]
sce$celltype <- factor(KEEP[as.character(sce$cell_type)], levels = unname(KEEP))
sce$group    <- ifelse(sce$disease == "MONDO_0005252", "HFpEF", "control")

counts <- as(assay(sce, "raw_count_cellranger"), "CsparseMatrix")
rownames(counts) <- rownames(sce); colnames(counts) <- colnames(sce)

## ---- 2. 单组建 CellChat(出版级:固定种子) -------------------------------
run_cellchat_one <- function(counts_sub, labels) {
  set.seed(SEED)                                       # 每组前重置种子
  data.norm <- normalizeData(counts_sub)
  meta <- data.frame(labels  = factor(labels),
                     samples = factor("sample1"),
                     row.names = colnames(counts_sub))

  cc <- createCellChat(object = data.norm, meta = meta, group.by = "labels")
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)

  # 关键:固定置换检验种子 + 显式 nboot;单线程下可复现
  cc <- computeCommunProb(cc, type = "truncatedMean", trim = 0.1,
                          population.size = TRUE,
                          seed.use = SEED, nboot = NBOOT)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc)
  cc
}

## ---- 3. 两组分别跑(单线程,慢) -------------------------------------------
groups <- c("control", "HFpEF")
cc_list <- list()
for (g in groups) {
  message(sprintf("\n========== 跑 CellChat(出版级): %s ==========", g))
  cells_g <- colnames(sce)[sce$group == g]
  t0 <- Sys.time()
  cc_list[[g]] <- run_cellchat_one(counts[, cells_g, drop = FALSE],
                                   droplevels(sce$celltype[sce$group == g]))
  message(sprintf("   %s 完成,用时 %.1f 分钟",
                  g, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  saveRDS(cc_list[[g]], file.path(out_dir, sprintf("cellchat_%s_final.rds", g)))
}

## ---- 4. 合并 + 比较 -------------------------------------------------------
message("\n== 合并两组并比较…")
cc_list <- lapply(cc_list, function(x){ x@meta$labels <- droplevels(x@meta$labels); x })
cellchat <- mergeCellChat(cc_list, add.names = names(cc_list))
saveRDS(cellchat, file.path(out_dir, "cellchat_merged_final.rds"))

## ---- 5. 结果打印 + 出图 ---------------------------------------------------
message("\n== 通路相对强度(出版级 p 值):")
gg_rank <- rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE)
rank_tab <- gg_rank$data[order(-abs(gg_rank$data$contribution.scaled)), ]
print(rank_tab[1:40, ])
write.csv(rank_tab, file.path(out_dir, "cc_pathway_rank_final.csv"), row.names = FALSE)
ggsave(file.path(out_dir, "cc_pathway_rank_final.png"), gg_rank, width = 6, height = 9, dpi = 150)

gg1 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1,2), measure = "weight")
ggsave(file.path(out_dir, "cc_overall_interactions_final.png"), gg1 + gg2,
       width = 8, height = 4, dpi = 150)

png(file.path(out_dir, "cc_diff_network_final.png"), width = 1000, height = 500)
par(mfrow = c(1,2), xpd = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")
dev.off()

for (g in groups) {
  p <- netAnalysis_signalingRole_scatter(cc_list[[g]], title = g)
  ggsave(file.path(out_dir, sprintf("cc_signalingRole_%s_final.png", g)),
         p, width = 5, height = 5, dpi = 150)
}

message(sprintf("\n== 完成(出版级)。结果存到 %s,文件名带 _final", out_dir))
message("   cc_pathway_rank_final.csv / .png(可复现 p 值的通路排名)")
message("   cellchat_merged_final.rds(后续深挖)")

print(sessionInfo())
