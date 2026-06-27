###############################################################################
# Step1续:ADGRL2 的 bulk 反向矛盾检查(load-bearing,推迟三轮、现必做)
#
# 背景:ADGRL2 bulk RVS = -0.30(padj 2.2e-5,高表达,robust下降),
#       但 snRNA CM 内在 = +0.43。"组分稀释"在此【方向反】——CM 是bulk主导
#       贡献者,CM内ADGRL2↑应推高bulk而非压低。故必须查清:
#   Q-A: ADGRL2 跨细胞类型【绝对表达分布】——它在CM是高基底还是低基底?
#        bulk的ADGRL2信号主要来自哪个区室?(若CM低基底+某非CM大区室HFpEF下降→调和)
#   Q-B: 各细胞类型【比例变化】HFpEF vs control(注:snRNA比例估计弱,只作参考)
#
# 两种结局:
#   ADGRL2 在CM低基底+主表达于HFpEF下降的非CM区室 → bulk调和,但旁分泌臂打小受体池=意义薄,带caveat收
#   ADGRL2 确CM主导而bulk仍降 → 矛盾未解 → ADGRL2不进锁定核心,作不一致观察另述
#
# 附:自分泌臂便宜加固——FLRT2&UNC5C在成纤维的【共表达比例】
#   (Fib→Fib自分泌要求同一成纤维同时有配体+受体,非"一部分有配体另一部分有受体")
#
# 数据:Hahn 2025 h5ad,raw_count_cellranger(原始整数),全13类
###############################################################################

suppressPackageStartupMessages({
  library(zellkonverter); library(SingleCellExperiment); library(SummarizedExperiment); library(Matrix)
})
h5      <- "/Users/franxy/Documents/博士课题 6.15/严谨/data/Hahn 人 HFpEF 单核 snRNA/HFpEF_snRNAseq.h5ad"
out_dir <- "/Users/franxy/Documents/博士课题 6.15/严谨/results"

message("== 读入 h5ad…")
sce <- readH5AD(h5, use_hdf5 = TRUE, reader = "R")
cd  <- as.data.frame(colData(sce))
cd$grp <- ifelse(cd$disease=="MONDO_0005252","HFpEF","control")

# ontology code → label
ct_map <- unique(cd[, c("cell_type","cell_type__ontology_label")])
ct_map <- ct_map[order(ct_map$cell_type), ]
code2label <- setNames(as.character(ct_map$cell_type__ontology_label), as.character(ct_map$cell_type))

GENES <- c("ADGRL2","ADGRL3","TENM3","FLRT2","UNC5C")   # 焦点ADGRL2;带几个对照
gi <- match(GENES, rownames(sce))
names(gi) <- GENES
if (any(is.na(gi))) message("[!] 未找到:", paste(GENES[is.na(gi)], collapse=", "))
gi <- gi[!is.na(gi)]

# 取这几个基因的 raw counts(genes × cells),实化
expr <- as.matrix(assay(sce, "raw_count_cellranger")[gi, , drop=FALSE])
rownames(expr) <- names(gi)

## ===== Q-A: 跨细胞类型绝对表达分布(每类的平均表达 + 表达占比) =====
cat("\n=========== Q-A:ADGRL2 等跨细胞类型【绝对表达分布】===========\n")
cat("看 ADGRL2 在哪类细胞高基底;bulk信号主要来自哪个区室(平均表达×细胞数)\n\n")

cell_ct <- as.character(cd$cell_type)
types <- names(sort(table(cell_ct), decreasing=TRUE))
for (g in rownames(expr)) {
  cat(sprintf("---- %s ----\n", g))
  rows <- lapply(types, function(tp){
    idx <- cell_ct==tp
    n <- sum(idx)
    mean_expr <- mean(expr[g, idx])          # 该类平均表达(CP计前的raw均值,看相对高低)
    pct_pos <- mean(expr[g, idx] > 0) * 100  # 表达阳性细胞%
    total_contrib <- sum(expr[g, idx])       # 该类对总表达的绝对贡献(均值×细胞数)
    data.frame(cell_type=code2label[[tp]], n_cells=n,
               mean_expr=mean_expr, pct_pos=pct_pos, total_contrib=total_contrib)
  })
  tab <- do.call(rbind, rows)
  tab$contrib_pct <- 100 * tab$total_contrib / sum(tab$total_contrib)  # 占总表达份额
  tab <- tab[order(-tab$contrib_pct), ]
  print(tab, row.names=FALSE, digits=3)
  # 关键标注:CM 是不是高基底 + CM 贡献份额
  cm_row <- tab[tab$cell_type=="cardiac muscle cell", ]
  if(nrow(cm_row)>0)
    cat(sprintf("   → 心肌(CM):平均表达=%.3f, 阳性%%=%.1f, 占总表达份额=%.1f%%\n\n",
                cm_row$mean_expr, cm_row$pct_pos, cm_row$contrib_pct))
}

## ===== Q-A 续:ADGRL2 在 CM 内 HFpEF vs control 的绝对均值(看是否小基数倍数变化) =====
cat("=========== ADGRL2 在心肌内 HFpEF vs control 绝对均值 ===========\n")
cm_idx <- cell_ct=="CL:0000746"
for(g in c("ADGRL2","TENM3")){
  hf_m <- mean(expr[g, cm_idx & cd$grp=="HFpEF"])
  ct_m <- mean(expr[g, cm_idx & cd$grp=="control"])
  cat(sprintf("  %s @CM: HFpEF均值=%.4f, control均值=%.4f (绝对差%.4f) — 看是否低基数\n",
              g, hf_m, ct_m, hf_m-ct_m))
}

## ===== Q-B: 各细胞类型比例变化(snRNA比例弱,仅参考) =====
cat("\n=========== Q-B:细胞类型比例 HFpEF vs control(snRNA比例估计弱,仅参考)===========\n")
prop_tab <- as.data.frame.matrix(table(cell_ct, cd$grp))
prop_tab$HFpEF_pct  <- 100*prop_tab$HFpEF/sum(prop_tab$HFpEF)
prop_tab$control_pct<- 100*prop_tab$control/sum(prop_tab$control)
prop_tab$delta_pct  <- prop_tab$HFpEF_pct - prop_tab$control_pct
prop_tab$label <- code2label[rownames(prop_tab)]
prop_tab <- prop_tab[order(prop_tab$delta_pct), c("label","HFpEF_pct","control_pct","delta_pct")]
print(prop_tab, row.names=FALSE, digits=3)
cat("   (负delta=HFpEF中该类比例下降;若某ADGRL2高表达的非CM类大幅下降→可能拉低bulk)\n")

## ===== 自分泌臂加固:FLRT2 & UNC5C 在成纤维的共表达 =====
cat("\n=========== 自分泌臂加固:成纤维内 FLRT2 & UNC5C 共表达 ===========\n")
cat("(Fib→Fib自分泌要求同一成纤维同时表达配体+受体;共表达%是金标准)\n")
fib_idx <- cell_ct=="CL:0002548"
f <- expr["FLRT2", fib_idx] > 0
u <- expr["UNC5C", fib_idx] > 0
n_fib <- sum(fib_idx)
cat(sprintf("  成纤维 n=%d\n", n_fib))
cat(sprintf("  FLRT2+: %.1f%%; UNC5C+: %.1f%%; 共表达(both+): %.1f%%\n",
            100*mean(f), 100*mean(u), 100*mean(f & u)))
cat(sprintf("  共表达/期望(若独立)= %.2f (>1 提示倾向共表达,支持自分泌)\n",
            mean(f & u) / (mean(f)*mean(u) + 1e-9)))
# 分HFpEF/control看共表达是否在HFpEF更高
for(grp in c("HFpEF","control")){
  gi2 <- fib_idx & cd$grp==grp
  f2 <- expr["FLRT2", gi2]>0; u2 <- expr["UNC5C", gi2]>0
  cat(sprintf("    %s: 共表达%%=%.1f (n=%d)\n", grp, 100*mean(f2&u2), sum(gi2)))
}

## ===== 判读 =====
cat("\n=========== 判读 ===========\n")
cat("Q-A 决定 ADGRL2 旁分泌臂命运:\n")
cat("  · 若CM是ADGRL2高基底且占bulk贡献大 → CM内↑应推高bulk,但bulk↓ → 矛盾未解 → ADGRL2不进锁定核心\n")
cat("  · 若CM是低基底、ADGRL2 bulk信号主来自某非CM区室、该区室HFpEF比例↓ → 调和bulk↓,但旁分泌臂打小受体池=意义薄,带caveat收\n")
cat("自分泌臂(FLRT2_UNC5C):无论如何已锁;共表达>期望进一步坐实Fib→Fib自分泌真实\n")

saveRDS(list(expr_summary="see stdout"), file.path(out_dir,"adgrl2_distribution_check.rds"))
cat("\n== sessionInfo =="); print(sessionInfo())
