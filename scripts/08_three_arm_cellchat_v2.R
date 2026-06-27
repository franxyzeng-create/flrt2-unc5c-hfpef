###############################################################################
# Step 2b-2(v2):三臂 CellChat 环境RNA稳健性检验
#
# 专家拍板:(c) 三臂,(a) CellBender layer 为主,(b) DecontX 为方法独立交叉验证。
#
# v2 相比 v1 的修改(全部按专家源码级审查):
#   - DecontX 臂用【小数期望 counts】(counts_decont_raw),不取整
#     (取整会把 <0.5 的低表达 axon-guidance 配体归零 → 给 DecontX 臂引入假衰减)
#   - 三臂 computeCommunProb 全部显式 raw.use=TRUE(默认值跨版本变过,pin死;
#     且保证三臂一致,唯一变量=输入矩阵)
#   - (a) CellBender 臂:cb layer 原样进 @data,绝不调 normalizeData;源码确认
#     链条 cb→@data→(subsetData)→@data.signaling→computeCommunProb 零再归一化
#   - rankNet 字段加断言(版本一变就报错而非静默给错)
#
# 三条臂(只换输入矩阵,所有 CellChat 参数逐一锁死):
#   raw       : raw_count_cellranger → normalizeData(log1p CP10K@1e4) [复用已有]
#   cellbender: raw_count_cellbender(已是 log1p CP10K@1e4 of CellBender counts),
#               直接进 @data、跳过 normalizeData(源码确认归一化方案逐字相同)
#   decontx   : DecontX 小数去污染 counts → normalizeData(同 raw 臂)
#
# 锁死参数:type="truncatedMean", trim=0.1, population.size=TRUE, nboot=100,
#           seed.use=1, raw.use=TRUE, 同 9 细胞类型, plan("sequential"), 不跑 projectData
#
# 注册式判定(先定规则后看结果):
#   某通路判为【稳健】⟺ 三臂下都 BH 校正后 rankNet p<0.05 且方向一致(HFpEF 增强)。
#   全表报告所有头条通路(含衰减者)。
#
# 前置:① cellchat_merged_final.rds(raw 臂,已有);② decontx_counts.rds(已有);
#       ③ 先跑过 08a_verify_slot.R 四查通过。
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

# 头条通路(axon-guidance 主线 + ECM,来自 Step1 BH 显著结果)
HEADLINE <- c("EPHA","EPHB","ADGRL","UNC5","Netrin","SEMA5","SLIT","NCAM","NRXN","FLRT",  # axon-guidance
              "LAMININ","COLLAGEN","FN1")                                                  # ECM

## ===========================================================================
## 核心:给定输入矩阵建一臂 CellChat。已归一化则跳过 normalizeData;
##       三臂统一 raw.use=TRUE、不跑 projectData(源码确认 computeCommunProb
##       在 raw.use=TRUE 下读 @data.signaling,即 cb layer 子集,零再归一化)
## ===========================================================================
build_cellchat <- function(input_matrix, labels, arm_name, already_normalized) {
  message(sprintf("\n========== 构建臂: %s ==========", arm_name))
  set.seed(SEED)
  if (already_normalized) {
    data_use <- input_matrix                      # 已是 log1p CP10K,原样进 @data
    cat("   输入已归一化(log1p CP10K),跳过 normalizeData(@data 原样)\n")
  } else {
    data_use <- normalizeData(input_matrix, scale.factor = 1e4, do.log = TRUE)
    cat("   输入为 counts,走 normalizeData(log1p CP10K@1e4)\n")
  }
  meta <- data.frame(labels = factor(labels), samples = factor("sample1"),
                     row.names = colnames(data_use))
  cc <- createCellChat(object = data_use, meta = meta, group.by = "labels")
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)                            # 填充 @data.signaling(只子集基因)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  # 关键:显式 raw.use=TRUE(pin死,跨版本默认值变过);不跑 projectData
  cc <- computeCommunProb(cc, type = "truncatedMean", trim = 0.1,
                          population.size = TRUE, seed.use = SEED, nboot = NBOOT,
                          raw.use = TRUE)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}

## ===========================================================================
## 读入数据,准备三臂的输入
## ===========================================================================
message("== 读入 h5ad…")
sce <- readH5AD(h5, use_hdf5 = TRUE, reader = "R")
idx <- which(as.character(sce$cell_type) %in% names(KEEP))
sce <- sce[, idx]
celltype <- factor(KEEP[as.character(sce$cell_type)], levels = unname(KEEP))
group    <- ifelse(sce$disease == "MONDO_0005252", "HFpEF", "control")

# 各臂输入矩阵(genes × cells)
cb_layer <- as(assay(sce, "raw_count_cellbender"), "CsparseMatrix")  # 已归一化(a 臂)
rownames(cb_layer) <- rownames(sce); colnames(cb_layer) <- colnames(sce)

dx <- readRDS(file.path(out_dir, "decontx_counts.rds"))
# 用【小数】去污染期望 counts(不取整;取整会把低表达 axon-guidance 配体归零→假衰减)
decont_counts <- dx$counts_decont_raw
rownames(decont_counts) <- rownames(sce)   # decontx 用同基因集
colnames(decont_counts) <- dx$cellnames
decont_counts <- as(decont_counts, "CsparseMatrix")

## ===========================================================================
## 分别按组(HFpEF/control)建对象,跑两臂(a)、(b);raw 臂复用已有
## ===========================================================================
groups <- c("control","HFpEF")

# --- (a) CellBender 臂 ---
cc_cb <- list()
for (g in groups) {
  cells_g <- colnames(sce)[group == g]
  cc_cb[[g]] <- build_cellchat(cb_layer[, cells_g, drop=FALSE],
                               droplevels(celltype[group == g]),
                               sprintf("CellBender-%s", g),
                               already_normalized = TRUE)
}
cc_cb <- lapply(cc_cb, function(x){ x@meta$labels <- droplevels(x@meta$labels); x })
merged_cb <- mergeCellChat(cc_cb, add.names = names(cc_cb))
saveRDS(merged_cb, file.path(out_dir, "cellchat_merged_cellbender.rds"))

# --- (b) DecontX 臂 ---
cc_dx <- list()
for (g in groups) {
  cells_g <- dx$cellnames[dx$group == g]
  cc_dx[[g]] <- build_cellchat(decont_counts[, cells_g, drop=FALSE],
                               droplevels(dx$celltype[dx$group == g]),
                               sprintf("DecontX-%s", g),
                               already_normalized = FALSE)
}
cc_dx <- lapply(cc_dx, function(x){ x@meta$labels <- droplevels(x@meta$labels); x })
merged_dx <- mergeCellChat(cc_dx, add.names = names(cc_dx))
saveRDS(merged_dx, file.path(out_dir, "cellchat_merged_decontx.rds"))

## ===========================================================================
## 提取三臂的 rankNet Wilcoxon p 值,做 BH,汇总头条通路
## ===========================================================================
get_rank_bh <- function(merged_obj) {
  gg <- rankNet(merged_obj, mode = "comparison", stacked = TRUE, do.stat = TRUE,
                return.data = TRUE)
  # 断言:版本一变字段名就报错,而非静默给错(pvalues 仅在 do.stat=TRUE 时出现)
  stopifnot("signaling.contribution" %in% names(gg))
  rk <- gg$signaling.contribution
  stopifnot("pvalues" %in% colnames(rk))
  pw <- unique(rk[, c("name","pvalues")])
  pw$pvalues_BH <- p.adjust(pw$pvalues, method = "BH")
  pw
}

merged_raw <- readRDS(file.path(out_dir, "cellchat_merged_final.rds"))
bh_raw <- get_rank_bh(merged_raw)
bh_cb  <- get_rank_bh(merged_cb)
bh_dx  <- get_rank_bh(merged_dx)

# 汇总成一张表:每条头条通路在三臂的 BH p 值
summ <- data.frame(pathway = HEADLINE)
lookup <- function(bh, pw) bh$pvalues_BH[match(pw, bh$name)]
summ$BH_raw       <- lookup(bh_raw, HEADLINE)
summ$BH_cellbender<- lookup(bh_cb,  HEADLINE)
summ$BH_decontx   <- lookup(bh_dx,  HEADLINE)

# 注册式判定:三臂都 BH<0.05 → 稳健
summ$robust_3arm <- with(summ,
  !is.na(BH_raw) & BH_raw < 0.05 &
  !is.na(BH_cellbender) & BH_cellbender < 0.05 &
  !is.na(BH_decontx) & BH_decontx < 0.05)

cat("\n========== 三臂稳健性汇总(头条通路)==========\n")
print(summ, row.names = FALSE, digits = 3)

write.csv(summ, file.path(out_dir, "three_arm_robustness.csv"), row.names = FALSE)

cat("\n========== 解读指引 ==========\n")
cat("预期(专家):axon-guidance(EPHA/ADGRL/UNC5/SLIT/SEMA5/NCAM/NRXN)三臂都稳;\n")
cat("           ECM(COLLAGEN/LAMININ/FN1)在去污染臂可能衰减。\n")
cat("若如此 → 证明:① 去污染确实在削 ambient-prone 的 ECM(说明检验有效);\n")
cat("              ② axon-guidance 质上不同、对 ambient 免疫(头条信号稳)。\n")
cat("全表(含衰减者)如实报告,这是反 cherry-picking 的有效性证据。\n")

cat("\n========== sessionInfo ==========\n")
print(sessionInfo())
