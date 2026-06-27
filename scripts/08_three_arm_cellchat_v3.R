###############################################################################
# Step 2b-2(v3):三臂 CellChat 环境RNA稳健性检验 —— 出版级
#
# v3 相比 v2 的修改(全部按专家第二轮审查):
#   ① 判定规则补"方向":robust ⟺ 三臂都 BH<0.05 且三臂都 flow_HFpEF > flow_control
#      (rankNet Wilcoxon 只测"有无差异"不测方向;漏方向会把 HFpEF↓的通路误判 robust)
#   ② raw 臂【本脚本重建】(不再 readRDS 复用 merged_final),三臂同一段代码同参数,
#      零管线漂移(消除"merged_final 当初是否同管线"的隐性混杂)
#   ③ 三臂 cell set 一致性断言(setequal)
#   ④ HEADLINE 通路名 vs CellChatDB 逐字核对(打印所有出现的 name,防静默 NA 误判)
#   ⑤ 防御性断言:每通路一个 p、decont 行列对齐、合并前两组 idents levels 一致
#
# 三臂(只换输入矩阵,参数全锁死):
#   raw       : raw_count_cellranger → normalizeData(log1p CP10K@1e4)
#   cellbender: raw_count_cellbender(已 log1p CP10K@1e4 of CellBender)→ 原样进@data,跳过normalizeData
#   decontx   : DecontX 小数去污染 counts → normalizeData
# 锁死:truncatedMean, trim=0.1, population.size=TRUE, nboot=100, seed.use=1,
#       raw.use=TRUE, min.cells=10, 同9细胞类型, sequential, 不跑projectData
###############################################################################

suppressPackageStartupMessages({
  library(zellkonverter); library(SingleCellExperiment); library(Matrix)
  library(CellChat); library(patchwork)
})

options(stringsAsFactors = FALSE)
future::plan("sequential")
set.seed(1)

out_dir <- "/Users/franxy/Documents/博士课题 6.15/严谨/results"
h5      <- "/Users/franxy/Documents/博士课题 6.15/严谨/data/Hahn 人 HFpEF 单核 snRNA/HFpEF_snRNAseq.h5ad"
SEED <- 1; NBOOT <- 100

KEEP <- c(
  "CL:0000746" = "Cardiomyocyte", "CL:0002548" = "Fibroblast",
  "CL:0010008" = "Endothelial1",  "CL:4033076" = "Endothelial2",
  "CL:0000669" = "Pericyte",      "CL:0000359" = "VSMC",
  "CL:0000235" = "Macrophage",    "CL:0000542" = "Lymphocyte",
  "CL:0002350" = "Endocardial"
)

HEADLINE <- c("EPHA","EPHB","ADGRL","UNC5","Netrin","SEMA5","SLIT","NCAM","NRXN","FLRT",
              "LAMININ","COLLAGEN","FN1")

## ---- 建一臂(三臂共用,保证同管线)----------------------------------------
build_cellchat <- function(input_matrix, labels, arm_name, already_normalized) {
  message(sprintf("\n========== 构建臂: %s ==========", arm_name))
  set.seed(SEED)
  if (already_normalized) {
    data_use <- input_matrix
    cat("   已归一化,原样进 @data,跳过 normalizeData\n")
  } else {
    data_use <- normalizeData(input_matrix, scale.factor = 1e4, do.log = TRUE)
    cat("   counts → normalizeData(log1p CP10K@1e4)\n")
  }
  meta <- data.frame(labels = factor(labels), samples = factor("sample1"),
                     row.names = colnames(data_use))
  cc <- createCellChat(object = data_use, meta = meta, group.by = "labels")
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "truncatedMean", trim = 0.1,
                          population.size = TRUE, seed.use = SEED, nboot = NBOOT,
                          raw.use = TRUE)              # 显式 pin
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}

# 两组建对象 + 合并(合并前断言 idents levels 一致)
build_merged <- function(mat, cells_by_group, labels_by_group, arm, already_norm) {
  ccs <- list()
  for (g in names(cells_by_group)) {
    cg <- cells_by_group[[g]]
    ccs[[g]] <- build_cellchat(mat[, cg, drop=FALSE],
                               droplevels(labels_by_group[[g]]),
                               sprintf("%s-%s", arm, g), already_norm)
  }
  ccs <- lapply(ccs, function(x){ x@meta$labels <- droplevels(x@meta$labels); x })
  # ⑧ 合并前确认两组 cell-type levels 一致(否则比较建在错位网络上)
  stopifnot(identical(levels(ccs[[1]]@idents), levels(ccs[[2]]@idents)))
  mergeCellChat(ccs, add.names = names(ccs))
}

## ===========================================================================
## 读入 + 准备三臂输入(同一批细胞)
## ===========================================================================
message("== 读入 h5ad…")
sce <- readH5AD(h5, use_hdf5 = TRUE, reader = "R")
idx <- which(as.character(sce$cell_type) %in% names(KEEP))
sce <- sce[, idx]
celltype <- factor(KEEP[as.character(sce$cell_type)], levels = unname(KEEP))
group    <- ifelse(sce$disease == "MONDO_0005252", "HFpEF", "control")
cell_ids <- colnames(sce)

cr_counts <- as(assay(sce, "raw_count_cellranger"), "CsparseMatrix")  # raw 臂
cb_layer  <- as(assay(sce, "raw_count_cellbender"),  "CsparseMatrix") # (a) 臂(已归一化)
rownames(cr_counts) <- rownames(sce); colnames(cr_counts) <- cell_ids
rownames(cb_layer)  <- rownames(sce); colnames(cb_layer)  <- cell_ids

dx <- readRDS(file.path(out_dir, "decontx_counts.rds"))
# ⑦ decont 行列对齐断言
stopifnot(ncol(dx$counts_decont_raw) == length(dx$cellnames),
          nrow(dx$counts_decont_raw) == nrow(sce))
decont_counts <- as(dx$counts_decont_raw, "CsparseMatrix")  # 小数,不取整
rownames(decont_counts) <- rownames(sce)
colnames(decont_counts) <- dx$cellnames

# ③ 三臂 cell set 一致性
stopifnot(setequal(cell_ids, dx$cellnames))
cat(sprintf(">> 三臂 cell set 一致性 OK:%d 个细胞\n", length(cell_ids)))

## ===========================================================================
## 三臂全部本脚本重建(同一段代码,零管线漂移)
## ===========================================================================
groups <- c("control","HFpEF")
cells_by_group  <- lapply(groups, function(g) cell_ids[group == g]); names(cells_by_group) <- groups
labels_by_group <- lapply(groups, function(g) celltype[group == g]); names(labels_by_group) <- groups

# raw 臂(重建,不复用 merged_final)
merged_raw <- build_merged(cr_counts, cells_by_group, labels_by_group, "raw", FALSE)
saveRDS(merged_raw, file.path(out_dir, "cellchat_merged_raw_v3.rds"))

# (a) CellBender 臂
merged_cb <- build_merged(cb_layer, cells_by_group, labels_by_group, "CellBender", TRUE)
saveRDS(merged_cb, file.path(out_dir, "cellchat_merged_cellbender.rds"))

# (b) DecontX 臂(细胞用 dx 的顺序与标签)
dx_cells_by_group  <- lapply(groups, function(g) dx$cellnames[dx$group == g]); names(dx_cells_by_group) <- groups
dx_labels_by_group <- lapply(groups, function(g) dx$celltype[dx$group == g]); names(dx_labels_by_group) <- groups
merged_dx <- build_merged(decont_counts, dx_cells_by_group, dx_labels_by_group, "DecontX", FALSE)
saveRDS(merged_dx, file.path(out_dir, "cellchat_merged_decontx.rds"))

## ===========================================================================
## 提取每臂:每通路 BH p 值 + 方向(flow_HFpEF − flow_control)
## ===========================================================================
get_rank_dir_bh <- function(merged_obj) {
  gg <- rankNet(merged_obj, mode = "comparison", stacked = TRUE, do.stat = TRUE,
                return.data = TRUE)
  stopifnot("signaling.contribution" %in% names(gg))
  rk <- gg$signaling.contribution
  stopifnot(all(c("name","pvalues","contribution","group") %in% colnames(rk)))
  # p 值(每通路一个)
  pw <- unique(rk[, c("name","pvalues")])
  stopifnot(nrow(pw) == length(unique(rk$name)))            # ⑥ 每通路一个 p
  pw$pvalues_BH <- p.adjust(pw$pvalues, method = "BH")
  # 方向:用原始 contribution(信息流),按通路 pivot 算 HFpEF − control
  flow <- tapply(rk$contribution, list(rk$name, rk$group), sum)
  # group 列里的名字应含 'HFpEF' 和 'control'
  gnames <- colnames(flow)
  hf <- gnames[grepl("HFpEF", gnames)][1]; ct <- gnames[grepl("control", gnames)][1]
  dir_df <- data.frame(name = rownames(flow),
                       flow_HFpEF = flow[, hf], flow_control = flow[, ct])
  dir_df$direction <- ifelse(dir_df$flow_HFpEF > dir_df$flow_control, "HFpEF_up", "HFpEF_down")
  merge(pw, dir_df, by = "name")
}

bh_raw <- get_rank_dir_bh(merged_raw)
bh_cb  <- get_rank_dir_bh(merged_cb)
bh_dx  <- get_rank_dir_bh(merged_dx)

## ④ HEADLINE 通路名 vs CellChatDB 逐字核对
all_names <- sort(unique(c(bh_raw$name, bh_cb$name, bh_dx$name)))
cat("\n========== 所有出现的 pathway 名(核对 HEADLINE 拼写)==========\n")
print(all_names)
missing <- setdiff(HEADLINE, all_names)
if (length(missing)) {
  cat("\n[!] 以下 HEADLINE 名在结果里找不到(可能拼写不符或真缺失,需核对):\n")
  print(missing)
} else {
  cat("\n>> 13 个 HEADLINE 名全部在结果中出现 ✓\n")
}

## ===========================================================================
## 汇总:每头条通路三臂的 BH p + 方向,注册式判定(BH<0.05 且三臂都 HFpEF↑)
## ===========================================================================
summ <- data.frame(pathway = HEADLINE)
look_p   <- function(bh) bh$pvalues_BH[match(HEADLINE, bh$name)]
look_dir <- function(bh) bh$direction[match(HEADLINE, bh$name)]
summ$BH_raw <- look_p(bh_raw); summ$dir_raw <- look_dir(bh_raw)
summ$BH_cb  <- look_p(bh_cb);  summ$dir_cb  <- look_dir(bh_cb)
summ$BH_dx  <- look_p(bh_dx);  summ$dir_dx  <- look_dir(bh_dx)

# 注册式判定:三臂都 BH<0.05 且 三臂方向都 HFpEF_up
sig3 <- with(summ, !is.na(BH_raw)&BH_raw<0.05 & !is.na(BH_cb)&BH_cb<0.05 & !is.na(BH_dx)&BH_dx<0.05)
up3  <- with(summ, dir_raw=="HFpEF_up" & dir_cb=="HFpEF_up" & dir_dx=="HFpEF_up")
up3[is.na(up3)] <- FALSE
summ$robust_3arm <- sig3 & up3

cat("\n========== 三臂稳健性汇总(头条通路)==========\n")
print(summ, row.names = FALSE, digits = 3)
write.csv(summ, file.path(out_dir, "three_arm_robustness.csv"), row.names = FALSE)

cat("\n========== 解读 ==========\n")
cat("robust_3arm=TRUE: 三臂都 BH<0.05 且都 HFpEF↑ → 经得起去污染的稳健头条信号\n")
cat("预期:axon-guidance(EPHA/ADGRL/UNC5/SLIT/SEMA5/NCAM/NRXN)稳健;\n")
cat("     ECM(COLLAGEN/LAMININ/FN1)在去污染臂可能衰减(BH 失显著)或方向变化。\n")
cat("全表如实报告(含衰减者)= 反 cherry-picking 的有效性证据。\n")
cat("Methods 注明:BH 的 family 是各臂检出的全部通路(三臂 family 可能不同)。\n")

cat("\n========== sessionInfo ==========\n")
print(sessionInfo())
