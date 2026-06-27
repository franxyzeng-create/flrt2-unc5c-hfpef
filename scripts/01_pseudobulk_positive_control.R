###############################################################################
# 第 0 步:pseudobulk 阳性对照复刻(成纤维 + 周细胞)
#
# 目的:确认"从 Hahn 的 h5ad → 按 donor 聚合 pseudobulk → limma-voom 差异分析"
#       这条管线在你本地跑得通,且结果与 Circulation Research 2025 原文一致。
#       这是后续 SCENIC / CellChat / 任何新方向的地基。
#
# 严格对应补充材料 Expanded Methods (Page 6-7) 的配方:
#   - 按 donor_id 聚合 counts(每人每细胞类型一个值)
#   - 某 donor 某细胞类型 ≤ 20 核 → 剔除
#   - 只测 HFpEF 和 control 各 ≥ 3 人的细胞类型
#   - 剔除线粒体 + 核糖体基因
#   - filterByExpr(group = disease)
#   - 模型 ~ disease + sex + pool   (pool = biosample_id,8 个测序池)
#   - CellBender 与 CellRanger 两套 counts 都跑,保留两套都显著且方向一致的基因
#
# 本版有意简化的地方(已知,下一步再补):
#   - 未实现背景污染启发式 gene_bkg × mean_ppv ≤ 0.40
#   - 用 TMM 归一化(原文用 DESeq2 size factors,二者高度等价)
#   => 因此 DE 基因数会与原文(成纤维 5905 / 周细胞 1812)同量级但不会完全相同。
#   => 真正的金标准对照 = 你的 logFC vs 原文 Table S8 的 logFC 相关(下一步接上)。
###############################################################################

suppressPackageStartupMessages({
  library(zellkonverter)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(edgeR)
  library(limma)
})

## ---- 0. 配置 --------------------------------------------------------------
h5      <- "/Users/franxy/Documents/博士课题 6.15/严谨/data/Hahn 人 HFpEF 单核 snRNA/HFpEF_snRNAseq.h5ad"
out_dir <- "/Users/franxy/Documents/博士课题 6.15/严谨/results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 目标细胞类型(用 Cell Ontology 编码,稳妥不依赖字符串):
#   CL:0002548 = 成纤维细胞 fibroblast (n=12562)
#   CL:0000669 = 周细胞 pericyte      (n=5100)
TARGETS <- c(Fibroblast = "CL:0002548",
             Pericyte   = "CL:0000669")

MIN_NUCLEI <- 20   # 某 donor 某细胞类型核数 ≤ 此值则剔除
MIN_DONORS <- 3    # 每组(HFpEF / control)至少需要的 donor 数

## ---- 1. 读入 h5ad(磁盘模式,省内存) -------------------------------------
message("== 读入 h5ad(use_hdf5=TRUE)…")
sce <- readH5AD(h5, use_hdf5 = TRUE, reader = "R")

## ---- 2. 锁定映射 + 完整性检查 ---------------------------------------------
message("\n== cell_type 编码 → ontology label 对照(确认翻译)：")
ct_map <- unique(as.data.frame(colData(sce))[, c("cell_type", "cell_type__ontology_label")])
ct_map <- ct_map[order(ct_map$cell_type), ]
print(ct_map, row.names = FALSE)

message("\n== disease 编码:MONDO_0005252 = HFpEF;PATO_0000461 = control")
print(table(sce$disease))

# 确认 donor → pool(biosample_id)是多对一(每个 donor 只属于一个池)
dp <- unique(as.data.frame(colData(sce))[, c("donor_id", "biosample_id")])
stopifnot(!any(duplicated(dp$donor_id)))   # 若报错说明 donor 跨池,假设不成立
message(sprintf("\n== donor 数 = %d;pool(biosample_id)数 = %d;donor→pool 唯一性 OK",
                length(unique(dp$donor_id)), length(unique(dp$biosample_id))))

## ---- 3. pseudobulk 聚合函数 -----------------------------------------------
# 输入某细胞类型编码,返回两套 counts 矩阵(genes × donors)+ 每个 donor 的元数据
make_pseudobulk <- function(sce, ct_code) {
  idx <- which(as.character(sce$cell_type) == ct_code)
  sub <- sce[, idx]
  donor <- factor(as.character(sub$donor_id))

  # 构造 cells × donors 的 0/1 指示矩阵(稀疏)
  ind <- sparseMatrix(i = seq_along(donor), j = as.integer(donor), x = 1,
                      dims = c(length(donor), nlevels(donor)))
  colnames(ind) <- levels(donor)

  counts <- list()
  for (layer in c("raw_count_cellranger", "raw_count_cellbender")) {
    # 把该细胞类型子集的 counts 实化为稀疏矩阵(genes × cells)
    m <- tryCatch(
      as(assay(sub, layer), "CsparseMatrix"),
      error = function(e) {
        stop(sprintf("实化 assay '%s' 失败:%s\n（若是内存问题,告诉我,我给 Python 导出方案）",
                     layer, conditionMessage(e)))
      })
    pb <- as.matrix(m %*% ind)               # genes × donors
    if (layer == "raw_count_cellbender") pb <- round(pb)  # 小数 → 整数
    counts[[layer]] <- pb
  }

  # 每个 donor 的元数据(disease/sex/pool 在同一 donor 内恒定)
  cd <- as.data.frame(colData(sub))
  meta <- cd[!duplicated(cd$donor_id),
             c("donor_id", "disease", "sex", "biosample_id")]
  rownames(meta) <- meta$donor_id
  meta <- meta[colnames(ind), , drop = FALSE]
  meta$n_nuclei <- as.integer(table(donor)[colnames(ind)])

  list(counts = counts, meta = meta)
}

## ---- 4. 差异分析函数(limma-voom) ----------------------------------------
run_de <- function(pb, layer) {
  meta <- pb$meta
  keep_d <- meta$n_nuclei > MIN_NUCLEI         # ≤20 核剔除
  meta <- meta[keep_d, ]
  counts <- pb$counts[[layer]][, rownames(meta), drop = FALSE]

  disease <- factor(ifelse(meta$disease == "MONDO_0005252", "HFpEF", "control"),
                    levels = c("control", "HFpEF"))
  n_hf <- sum(disease == "HFpEF"); n_ct <- sum(disease == "control")
  if (n_hf < MIN_DONORS || n_ct < MIN_DONORS)
    stop(sprintf("可用 donor 不足:HFpEF=%d, control=%d", n_hf, n_ct))

  sex  <- factor(meta$sex)
  pool <- droplevels(factor(meta$biosample_id))

  # 剔除线粒体 + 核糖体基因
  g <- rownames(counts)
  drop_g <- grepl("^MT-", g) | grepl("^RP[LS]", g)
  counts <- counts[!drop_g, , drop = FALSE]

  design <- model.matrix(~ disease + sex + pool)

  dge  <- DGEList(counts)
  keep <- filterByExpr(dge, group = disease)
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- calcNormFactors(dge)                 # TMM
  v    <- voom(dge, design)
  fit  <- eBayes(lmFit(v, design))
  tt   <- topTable(fit, coef = "diseaseHFpEF", number = Inf, sort.by = "none")

  attr(tt, "n_donors") <- c(HFpEF = n_hf, control = n_ct)
  attr(tt, "n_genes_tested") <- nrow(tt)
  tt
}

## ---- 5. 逐细胞类型跑,并做双层一致性 + 命名基因核对 ----------------------
# 原文点名的 top 基因(用于肉眼核对方向是否一致)
NAMED <- list(
  Fibroblast = list(
    up   = c("TLL2","CYS1","NMD3","LDB2","FREM1","NR3C1","PRTFDC1","SYT17"),
    down = c("FKBP5","AOX1","MGST1","ITGA5","NID1","TGFBR3","WASF3","JADE1")),
  Pericyte = list(
    up   = c("ANTXR1","ITIH5","NREP","PLA2G4C","DGKB","NTF3","PITPNM2","TCF7L1"),
    down = c("SPARCL1","LPP","CDC42EP4","CPM","MAP1B","MDM2","TIMP3","SLC39A11"))
)

# 原文报告的 cell-type DE 基因数(供量级对比)
PAPER_N <- c(Fibroblast = 5905, Pericyte = 1812)

results <- list()
for (nm in names(TARGETS)) {
  message(sprintf("\n========== %s (%s) ==========", nm, TARGETS[[nm]]))
  pb <- make_pseudobulk(sce, TARGETS[[nm]])

  de_cr <- run_de(pb, "raw_count_cellranger")
  de_cb <- run_de(pb, "raw_count_cellbender")   # 原文 logFC 以 CellBender 为准

  message(sprintf("可用 donor:HFpEF=%d, control=%d",
                  attr(de_cb,"n_donors")["HFpEF"], attr(de_cb,"n_donors")["control"]))

  # 双层一致性:两套都 FDR<0.05 且 logFC 方向一致
  common <- intersect(rownames(de_cb), rownames(de_cr))
  cb <- de_cb[common, ]; cr <- de_cr[common, ]
  concordant <- cb$adj.P.Val < 0.05 & cr$adj.P.Val < 0.05 &
                sign(cb$logFC) == sign(cr$logFC)
  n_de <- sum(concordant)

  message(sprintf("CellBender 单层 FDR<0.05: %d 个;CellRanger 单层: %d 个",
                  sum(de_cb$adj.P.Val < 0.05), sum(de_cr$adj.P.Val < 0.05)))
  message(sprintf(">> 双层一致 DE 基因数 = %d   (原文 cell-type DE = %d,本版简化故应为同量级)",
                  n_de, PAPER_N[[nm]]))
  message(sprintf("   logFC 两层相关 r = %.3f", cor(cb$logFC, cr$logFC)))

  # 命名基因方向核对(用 CellBender 结果)
  message("   —— 命名基因核对(CellBender logFC / FDR;原文方向见括号)——")
  chk <- function(genes, dir_label) {
    sub <- de_cb[rownames(de_cb) %in% genes, c("logFC","adj.P.Val")]
    if (nrow(sub)) {
      sub <- sub[order(match(rownames(sub), genes)), ]
      for (g in rownames(sub))
        message(sprintf("     %-10s logFC=%+.2f  FDR=%.1e   (原文 %s)",
                        g, sub[g,"logFC"], sub[g,"adj.P.Val"], dir_label))
    }
  }
  chk(NAMED[[nm]]$up,   "上调")
  chk(NAMED[[nm]]$down, "下调")

  # 合并输出表
  merged <- data.frame(
    gene = common,
    logFC_cellbender = cb$logFC, FDR_cellbender = cb$adj.P.Val,
    logFC_cellranger = cr$logFC, FDR_cellranger = cr$adj.P.Val,
    concordant_DE = concordant, row.names = NULL)
  merged <- merged[order(merged$FDR_cellbender), ]
  write.csv(merged, file.path(out_dir, sprintf("DE_%s_pseudobulk.csv", nm)),
            row.names = FALSE)
  results[[nm]] <- list(merged = merged, pb = pb)
}

## ---- 6. 保存 + 结束 --------------------------------------------------------
saveRDS(results, file.path(out_dir, "pseudobulk_positive_control.rds"))
message(sprintf("\n== 完成。结果已保存到:%s", out_dir))
message("   - DE_Fibroblast_pseudobulk.csv")
message("   - DE_Pericyte_pseudobulk.csv")
message("   - pseudobulk_positive_control.rds")

message("\n== sessionInfo ==")
print(sessionInfo())
