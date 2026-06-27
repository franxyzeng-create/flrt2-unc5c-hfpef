###############################################################################
# Step1续:FLRT2 & UNC5C 成纤维共表达 — per-donor 统计(去 pseudoreplication)
#
# 背景:18脚本报的池化共表达 HFpEF 40.4% vs control 20.3% 是 cell-level 池化
#   = pseudoreplication(几个双阳性多的供者就能拉动总数)。专家要求改 per-donor:
#   每供者内部算成纤维双阳性比例 → Wilcoxon 比组 → 中位/IQR/n/p/效应量。
#   这是自分泌臂【最强单项证据】(单细胞共表达=自分泌金标准),做正式结果。
#
# 技术 caveat(写明):
#   ① 检测阈值≥1 count;snRNA稀疏+dropout会把真双阳性误判单阳性 → 绝对比例对阈值敏感,
#      但同阈值下的相对比较(HFpEF vs control)更稳 → 框成相对 + 附阈值稳健性(≥1/≥2)
#   ② 翻倍可能部分来自两基因整体表达升高→更多细胞过检测线,而非严格协同;
#      但仍支持"更多成纤维具备自分泌能力"
#
# 数据:Hahn 2025 h5ad,raw_count_cellranger,成纤维 CL:0002548;数据集内donor对照
###############################################################################

suppressPackageStartupMessages({
  library(zellkonverter); library(SingleCellExperiment); library(SummarizedExperiment); library(Matrix)
})
h5      <- "/Users/franxy/Documents/博士课题 6.15/严谨/data/Hahn 人 HFpEF 单核 snRNA/HFpEF_snRNAseq.h5ad"
out_dir <- "/Users/franxy/Documents/博士课题 6.15/严谨/results"
MIN_FIB_NUCLEI <- 20   # 某供者成纤维核数<此值则该供者不纳入(比例不稳)

message("== 读入 h5ad…")
sce <- readH5AD(h5, use_hdf5 = TRUE, reader = "R")
cd  <- as.data.frame(colData(sce))
cd$grp <- ifelse(cd$disease=="MONDO_0005252","HFpEF","control")

# 成纤维子集
fib_idx <- as.character(cd$cell_type) == "CL:0002548"
cat(sprintf("成纤维总核数:%d\n", sum(fib_idx)))

# FLRT2 & UNC5C 的 raw counts(成纤维细胞)
gi <- match(c("FLRT2","UNC5C"), rownames(sce)); names(gi) <- c("FLRT2","UNC5C")
expr_fib <- as.matrix(assay(sce, "raw_count_cellranger")[gi, fib_idx, drop=FALSE])
rownames(expr_fib) <- names(gi)
donor_fib <- as.character(cd$donor_id[fib_idx])
grp_fib   <- cd$grp[fib_idx]

## ---- per-donor 双阳性比例(阈值≥1 与 ≥2 两套,稳健性) ----
per_donor_coexpr <- function(thresh) {
  donors <- unique(donor_fib)
  rows <- lapply(donors, function(d){
    sel <- donor_fib == d
    n_fib <- sum(sel)
    if (n_fib < MIN_FIB_NUCLEI) return(NULL)   # 核太少不纳入
    f <- expr_fib["FLRT2", sel] >= thresh
    u <- expr_fib["UNC5C", sel] >= thresh
    data.frame(donor=d, grp=unique(grp_fib[sel]), n_fib=n_fib,
               pct_FLRT2=100*mean(f), pct_UNC5C=100*mean(u),
               pct_coexpr=100*mean(f & u))
  })
  do.call(rbind, rows)
}

for (th in c(1, 2)) {
  cat(sprintf("\n================= 阈值 ≥%d count =================\n", th))
  pd <- per_donor_coexpr(th)
  cat(sprintf("纳入供者:%d(成纤维核≥%d);HFpEF=%d, control=%d\n",
              nrow(pd), MIN_FIB_NUCLEI, sum(pd$grp=="HFpEF"), sum(pd$grp=="control")))

  hf <- pd$pct_coexpr[pd$grp=="HFpEF"]
  ct <- pd$pct_coexpr[pd$grp=="control"]
  # Wilcoxon
  w <- wilcox.test(hf, ct)
  cat(sprintf("\n成纤维双阳性(FLRT2+ & UNC5C+)比例 per-donor:\n"))
  cat(sprintf("  HFpEF : 中位=%.1f%% IQR=[%.1f, %.1f]  (n=%d)\n",
              median(hf), quantile(hf,.25), quantile(hf,.75), length(hf)))
  cat(sprintf("  control: 中位=%.1f%% IQR=[%.1f, %.1f]  (n=%d)\n",
              median(ct), quantile(ct,.25), quantile(ct,.75), length(ct)))
  cat(sprintf("  Wilcoxon p=%.4g\n", w$p.value))
  # 效应量:Cliff's delta(非参,稳健)
  cliff <- function(a,b){ m<-outer(a,b,">"); n<-outer(a,b,"<"); (sum(m)-sum(n))/(length(a)*length(b)) }
  cat(sprintf("  Cliff's delta=%.3f (>0 表 HFpEF 更高;|0.33|中,|0.47|大)\n", cliff(hf,ct)))
  cat(sprintf("  中位倍数:HFpEF/control = %.2f×\n", median(hf)/(median(ct)+1e-9)))

  # 也报单阳性(看双阳性升高是否只是单基因驱动)
  cat(sprintf("\n  对照(单阳性中位):FLRT2+ HFpEF=%.1f%% vs ctrl=%.1f%%;UNC5C+ HFpEF=%.1f%% vs ctrl=%.1f%%\n",
              median(pd$pct_FLRT2[pd$grp=="HFpEF"]), median(pd$pct_FLRT2[pd$grp=="control"]),
              median(pd$pct_UNC5C[pd$grp=="HFpEF"]), median(pd$pct_UNC5C[pd$grp=="control"])))

  if (th==1) { pd_main <- pd; w_main <- w }   # 主结果用≥1
}

## ---- 保存主结果(阈值≥1) ----
cat("\n=========== per-donor 明细(阈值≥1)===========\n")
pd_main <- pd_main[order(pd_main$grp, -pd_main$pct_coexpr), ]
print(pd_main, row.names=FALSE, digits=3)
write.csv(pd_main, file.path(out_dir, "coexpr_per_donor.csv"), row.names=FALSE)

cat("\n---- 判读 ----\n")
cat("· per-donor 后若仍显著(Wilcoxon p<0.05 + Cliff's delta 中-大)→ 共表达翻倍非pseudoreplication,可作正式结果\n")
cat("· 措辞:报相对比较(同阈值HFpEF vs control),绝对比例标注阈值敏感\n")
cat("· 若双阳性升高 > 单阳性升高 → 提示协同(不只单基因驱动);否则诚实写'更多成纤维具备自分泌能力'\n")
cat("· 阈值≥1 vs ≥2 结论一致 → 稳健\n")

cat("\n== sessionInfo =="); print(sessionInfo())
