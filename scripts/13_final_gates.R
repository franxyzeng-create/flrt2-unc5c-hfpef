###############################################################################
# 收口前两道闸:ngenes 检出广度 + sex × disease 平衡
#
# 专家指出的两个 reviewer 必问、手头数据一行就能查的核查,补完再收口 Step 1-3。
#
# ① ngenes(每细胞检出基因数):闭合"检出广度→更多L-R对点亮→总流量虚高"这条机制。
#    比 ncount/UMI 更直接衡量检出广度。HFpEF 深度已略低,ngenes 大概率也相当/更低→强化非技术结论。
# ② sex × disease 平衡:HFpEF 临床女性占多,若队列 HFpEF 偏女、control 性别不同,
#    sex 与 disease 共线,性别二态通路会显出假疾病效应。这是 HFpEF 论文必查项。
###############################################################################

suppressPackageStartupMessages({
  library(zellkonverter); library(SingleCellExperiment); library(Matrix)
})
options(stringsAsFactors = FALSE)

out_dir <- "/Users/franxy/Documents/博士课题 6.15/严谨/results"
h5      <- "/Users/franxy/Documents/博士课题 6.15/严谨/data/Hahn 人 HFpEF 单核 snRNA/HFpEF_snRNAseq.h5ad"

FIXED_SET <- c("Cardiomyocyte","Fibroblast","Endothelial1","Pericyte")
KEEP <- c(
  "CL:0000746"="Cardiomyocyte","CL:0002548"="Fibroblast","CL:0010008"="Endothelial1",
  "CL:4033076"="Endothelial2","CL:0000669"="Pericyte","CL:0000359"="VSMC",
  "CL:0000235"="Macrophage","CL:0000542"="Lymphocyte","CL:0002350"="Endocardial"
)

message("== 读入 h5ad…")
sce <- readH5AD(h5, use_hdf5 = TRUE, reader = "R")
idx <- which(as.character(sce$cell_type) %in% names(KEEP))
sce <- sce[, idx]
celltype <- factor(KEEP[as.character(sce$cell_type)], levels = unname(KEEP))
donor <- as.character(sce$donor_id)
group <- ifelse(sce$disease == "MONDO_0005252", "HFpEF", "control")

# Step 3 保留的 42 病人
pd <- readRDS(file.path(out_dir, "perdonor_final.rds"))
kept <- pd$kept_donors
dg   <- pd$donor_group[kept]
cat(sprintf("核查病人: %d(HFpEF %d, control %d)\n", length(kept),
            sum(dg=="HFpEF"), sum(dg=="control")))

## ===========================================================================
## ① ngenes:检出广度(donor 层面中位,4类集内)
## ===========================================================================
cat("\n========== ① ngenes 检出广度(每细胞检出基因数)==========\n")
ngenes_cell <- as.numeric(colData(sce)[["ngenes"]])
names(ngenes_cell) <- colnames(sce)

ng_by_donor <- sapply(kept, function(d) {
  sel <- donor == d & celltype %in% FIXED_SET
  median(ngenes_cell[sel], na.rm = TRUE)
})
ng_hf <- ng_by_donor[dg=="HFpEF"]; ng_ct <- ng_by_donor[dg=="control"]
cat("HFpEF 中位 ngenes:\n");  print(summary(ng_hf))
cat("control 中位 ngenes:\n");print(summary(ng_ct))
p_ng <- suppressWarnings(wilcox.test(ng_hf, ng_ct)$p.value)
cat(sprintf(">> ngenes 组间 Wilcoxon p = %.3f\n", p_ng))

# 方向判读(专家:深度已control更高,关键看方向不只看p)
dir_ng <- if (median(ng_hf) < median(ng_ct)) "HFpEF更低" else if (median(ng_hf) > median(ng_ct)) "HFpEF更高" else "持平"
cat(sprintf(">> 方向:HFpEF 中位 %.0f vs control 中位 %.0f(%s)\n",
            median(ng_hf), median(ng_ct), dir_ng))
if (median(ng_hf) <= median(ng_ct)) {
  cat("   → HFpEF ngenes ≤ control(与深度一致,bias against)→ 干净闭合检出广度:\n")
  cat("      总流量升高 despite 等量/更少检出基因 → 加固非技术结论\n")
} else {
  cat("   → [注意] HFpEF ngenes 反而更高 despite 深度更低——这是唯一会 complicate 的情形,需细看\n")
}

# 总流量 ~ ngenes 相关
total_flow <- colSums(pd$S)[kept]
rho_ng <- cor(total_flow, ng_by_donor[kept], method="spearman")
p_rho_ng <- suppressWarnings(cor.test(total_flow, ng_by_donor[kept], method="spearman")$p.value)
cat(sprintf(">> 总流量 ~ ngenes:Spearman rho = %.3f, p = %.3g\n", rho_ng, p_rho_ng))

## ===========================================================================
## ② sex × disease 平衡(donor 层面)
## ===========================================================================
cat("\n========== ② sex × disease 平衡 ==========\n")
donor_sex <- tapply(as.character(sce$sex), donor, function(x) x[1])
sex_kept <- donor_sex[kept]
sex_tab <- table(sex = sex_kept, disease = dg)
cat("性别 × 疾病 列联表(42 病人):\n")
print(sex_tab)
p_sex <- suppressWarnings(fisher.test(sex_tab)$p.value)
cat(sprintf(">> Fisher 精确检验 p = %.3f\n", p_sex))

# 各组性别比例 + 偏斜判断(专家:n=18v24 低功效,Fisher p 不充分,必须看实际比例)
cat("\n各组性别构成:\n")
prop_by_grp <- list()
for (g in c("HFpEF","control")) {
  s <- sex_kept[dg==g]; tb <- table(s)
  prop_by_grp[[g]] <- tb / sum(tb)
  cat(sprintf("  %s: %s\n", g, paste(sprintf("%s=%d(%.0f%%)", names(tb), tb, 100*tb/sum(tb)), collapse=", ")))
}

# 偏斜规则:任一组某性别>70%,或两组同一性别比例差>20个百分点 → 需分层敏感性
sexes <- union(names(prop_by_grp[["HFpEF"]]), names(prop_by_grp[["control"]]))
getp <- function(g, s) { v <- prop_by_grp[[g]][s]; ifelse(is.na(v), 0, v) }
max_skew <- max(sapply(c("HFpEF","control"), function(g) max(prop_by_grp[[g]])))
max_diff <- max(sapply(sexes, function(s) abs(getp("HFpEF", s) - getp("control", s))))
cat(sprintf("\n最大单组偏斜 = %.0f%%;两组最大比例差 = %.0f 个百分点\n",
            100*max_skew, 100*max_diff))

need_strat <- (max_skew > 0.70) || (max_diff > 0.20) || (p_sex < 0.05)
if (need_strat) {
  cat(">> 【需要性别分层敏感性】判据触发(单组>70% 或 两组差>20pp 或 Fisher<0.05):\n")
  cat("   对硬核心 UNC5/FLRT/ADGRL 补:按 sex 分层的 van Elteren,或秩回归带 sex 协变量\n")
  cat("   (即便某核心通路分层后减弱,也按 NRXN 同等诚实处理:'overall显著、性别分层后减弱',核心进一步收窄)\n")
} else {
  cat(">> 【性别基本均衡】无单组>70%、两组差≤20pp、Fisher≥0.05\n")
  cat("   → Methods 写'两组性别均衡(Fisher p=X,各组比例见表)',sex 非主混杂,无需分层\n")
}

## ---- 汇总 ----------------------------------------------------------------
cat("\n========== 两道闸汇总 ==========\n")
cat(sprintf("① ngenes: 组间 p=%.3f(%s), 与总流量 rho=%.3f(p=%.3g)\n",
            p_ng, dir_ng, rho_ng, p_rho_ng))
cat(sprintf("② sex×disease: Fisher p=%.3f, 最大单组偏斜%.0f%%, 两组差%.0fpp → %s\n",
            p_sex, 100*max_skew, 100*max_diff,
            if (need_strat) "需性别分层敏感性" else "性别均衡、无需分层"))
if (!need_strat && median(ng_hf) <= median(ng_ct)) {
  cat("两闸都过(ngenes≤control闭合检出广度 + 性别均衡)→ Step 1-3 可收口,无残留混杂洞\n")
} else {
  cat("有闸需跟进(见上),处理后再收口\n")
}

cat("\n========== sessionInfo ==========\n")
print(sessionInfo())
