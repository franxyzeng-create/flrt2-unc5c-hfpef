###############################################################################
# Step1续(CB一致性):用 CellBender 层(去ambient)重跑两个怕ambient的指标
#
# 专家指出:
#   ① 共表达双阳性比例最怕ambient(几个ambient count就把细胞从阴翻阳),且偏倚方向
#      与结论一致(HFpEF FLRT2高→ambient FLRT2多→HFpEF端双阳性被特异性抬高)。
#      19脚本用了 raw_count_cellranger(含ambient) → 必须CB层重跑。
#   ② UNC5C"成纤维富集"措辞需要和沉掉ADGRL2同一把分布尺子,且18也是CR层算的
#      → CB层重确认 UNC5C 是否真成纤维富集(决定"稀释"vs"uninformative"措辞)。
#
# 判读:
#   共表达:CB层撑住(仍p显著+倍数在)→ambient排除;缩水→ambient灌水。
#   UNC5C分布:CB层若成纤维富集(非CM主导)→"bulk-ns=稀释"措辞OK;
#             若非富集→软化为"受体就位;bulk对UNC5C uninformative"。
#
# 数据:raw_count_cellbender(CellBender去污染,已round成整数 — 见01脚本)
###############################################################################

suppressPackageStartupMessages({
  library(zellkonverter); library(SingleCellExperiment); library(SummarizedExperiment); library(Matrix)
})
h5      <- "/Users/franxy/Documents/博士课题 6.15/严谨/data/Hahn 人 HFpEF 单核 snRNA/HFpEF_snRNAseq.h5ad"
out_dir <- "/Users/franxy/Documents/博士课题 6.15/严谨/results"
LAYER   <- "raw_count_cellbender"      # ★去ambient层
MIN_FIB_NUCLEI <- 20

message("== 读入 h5ad…")
sce <- readH5AD(h5, use_hdf5 = TRUE, reader = "R")
cd  <- as.data.frame(colData(sce))
cd$grp <- ifelse(cd$disease=="MONDO_0005252","HFpEF","control")
ct_map <- unique(cd[, c("cell_type","cell_type__ontology_label")])
code2label <- setNames(as.character(ct_map$cell_type__ontology_label), as.character(ct_map$cell_type))
cell_ct <- as.character(cd$cell_type)

# CellBender 层是小数,需round成整数(与01脚本一致)再做阈值
get_expr <- function(genes){
  gi <- match(genes, rownames(sce)); names(gi) <- genes
  m <- as.matrix(assay(sce, LAYER)[gi, , drop=FALSE]); rownames(m) <- genes
  round(m)   # CB层round成整数(同01脚本)
}

## ===== 1. CB层 共表达 per-donor(双阈值) =====
cat("\n############ 1. 共表达 per-donor — CellBender层(去ambient) ############\n")
fib_idx <- cell_ct == "CL:0002548"
expr_fib <- get_expr(c("FLRT2","UNC5C"))[, fib_idx, drop=FALSE]
donor_fib <- as.character(cd$donor_id[fib_idx]); grp_fib <- cd$grp[fib_idx]

per_donor_coexpr <- function(thresh){
  donors <- unique(donor_fib)
  rows <- lapply(donors, function(d){
    sel <- donor_fib==d; if (sum(sel) < MIN_FIB_NUCLEI) return(NULL)
    f <- expr_fib["FLRT2",sel]>=thresh; u <- expr_fib["UNC5C",sel]>=thresh
    data.frame(donor=d, grp=unique(grp_fib[sel]), n_fib=sum(sel),
               pct_FLRT2=100*mean(f), pct_UNC5C=100*mean(u), pct_coexpr=100*mean(f&u))
  })
  do.call(rbind, rows)
}
cliff <- function(a,b){ (sum(outer(a,b,">"))-sum(outer(a,b,"<")))/(length(a)*length(b)) }

for (th in c(1,2)){
  pd <- per_donor_coexpr(th)
  hf <- pd$pct_coexpr[pd$grp=="HFpEF"]; ct <- pd$pct_coexpr[pd$grp=="control"]
  w <- wilcox.test(hf, ct)
  cat(sprintf("\n--- 阈值≥%d (CB层) ---\n", th))
  cat(sprintf("  双阳性: HFpEF中位=%.1f%% [%.1f,%.1f] vs control=%.1f%% [%.1f,%.1f]\n",
              median(hf),quantile(hf,.25),quantile(hf,.75),
              median(ct),quantile(ct,.25),quantile(ct,.75)))
  cat(sprintf("  Wilcoxon p=%.4g; Cliff's delta=%.3f; 倍数=%.2f×\n",
              w$p.value, cliff(hf,ct), median(hf)/(median(ct)+1e-9)))
  cat(sprintf("  单阳性中位: FLRT2+ HFpEF=%.1f%%/ctrl=%.1f%%; UNC5C+ HFpEF=%.1f%%/ctrl=%.1f%%\n",
              median(pd$pct_FLRT2[pd$grp=="HFpEF"]),median(pd$pct_FLRT2[pd$grp=="control"]),
              median(pd$pct_UNC5C[pd$grp=="HFpEF"]),median(pd$pct_UNC5C[pd$grp=="control"])))
  # 验证乘积独立性(专家:共表达≈FLRT2+×UNC5C+,非协同)
  pf<-median(pd$pct_FLRT2[pd$grp=="HFpEF"])/100; pu<-median(pd$pct_UNC5C[pd$grp=="HFpEF"])/100
  cat(sprintf("  [独立性检验]HFpEF: FLRT2+×UNC5C+=%.1f%% vs 实测双阳%.1f%% → %s\n",
              100*pf*pu, median(hf), if(abs(100*pf*pu-median(hf))<5) "≈乘积(无协同,受体就位+配体诱导)" else "偏离乘积"))
  if (th==1) pd_cb <- pd
}
write.csv(pd_cb, file.path(out_dir,"coexpr_per_donor_cellbender.csv"), row.names=FALSE)

## ===== 2. CB层 UNC5C/ADGRL2 跨细胞类型分布(同一把尺子) =====
cat("\n############ 2. UNC5C/ADGRL2 跨细胞类型分布 — CellBender层 ############\n")
cat("(同沉掉ADGRL2的尺子;确认UNC5C是否真成纤维富集→定'稀释'vs'uninformative'措辞)\n")
expr_all <- get_expr(c("UNC5C","ADGRL2"))
types <- names(sort(table(cell_ct), decreasing=TRUE))
for (g in c("UNC5C","ADGRL2")){
  cat(sprintf("\n---- %s (CB层) ----\n", g))
  rows <- lapply(types, function(tp){
    idx <- cell_ct==tp
    data.frame(cell_type=code2label[[tp]], n_cells=sum(idx),
               mean_expr=mean(expr_all[g,idx]), pct_pos=100*mean(expr_all[g,idx]>0),
               total_contrib=sum(expr_all[g,idx]))
  })
  tab <- do.call(rbind, rows)
  tab$contrib_pct <- 100*tab$total_contrib/sum(tab$total_contrib)
  tab <- tab[order(-tab$contrib_pct),]
  print(head(tab,6), row.names=FALSE, digits=3)
  fib <- tab[tab$cell_type=="fibroblast of cardiac tissue",]
  cm  <- tab[tab$cell_type=="cardiac muscle cell",]
  cat(sprintf("  → 成纤维占总表达 %.1f%% (均值%.2f); 心肌占 %.1f%% (均值%.2f)\n",
              fib$contrib_pct, fib$mean_expr, cm$contrib_pct, cm$mean_expr))
  if (g=="UNC5C")
    cat(sprintf("  → UNC5C %s → bulk-ns措辞:%s\n",
        if(fib$contrib_pct>cm$contrib_pct) "成纤维富集" else "非成纤维主导",
        if(fib$contrib_pct>cm$contrib_pct) "'稀释'OK" else "软化为'受体就位;bulk uninformative'"))
}

cat("\n---- 判读 ----\n")
cat("· 共表达CB层撑住(p显著+倍数在)→ambient排除;缩水→ambient灌水(原CR结果不可用)\n")
cat("· 独立性:共表达≈FLRT2+×UNC5C+ → 写'受体就位+配体诱导使自分泌成纤维群扩大',不写协同/第二颗子弹\n")
cat("· UNC5C若CB层仍成纤维富集→'稀释'措辞有据;否则软化\n")
cat("· 注:本步只闭'真实/去ambient',HFpEF特异性仍待BMI检验\n")

cat("\n== sessionInfo =="); print(sessionInfo())
