###############################################################################
# Step 3(最终检验):per-donor 病人层面推断 —— 根除伪重复,冲2区的统计主干
#
# 把分析单元从"细胞对"(rankNet Wilcoxon 的伪重复)提升到"病人"(18 HFpEF vs 24 control)。
#
# 设计(专家逐条拍定):
#   - 固定 4 类集:Cardiomyocyte, Fibroblast, Endothelial1, Pericyte
#     (Endocardial 退出:其达标率组间不对称→差异化丢失会把疾病相关混杂重新引入;
#      且 population.size 加权下它是 UNC5/Netrin/FLRT 承载者→纳入近乎循环论证)
#   - 每病人 4 类都需 ≥10 核;不足者剔除(预期 HFpEF 18, control 24)
#   - summary = apply(netP$prob,3,sum)(= rankNet contribution 同口径);nboot=1 提速
#   - 输入 raw_count_cellranger(与主线 raw 臂一致)
#   - 【绝对流量版(主)】population.size=TRUE,匹配 Step 2
#   - 【相对流量版(co-primary)】每病人每通路 ÷ 该病人全通路总流量(对丰度/类型数稳健)
#   - 报每病人总流量按组分布(混杂核查)
#   - 主检验:每通路 Wilcoxon-on-full(含0,不丢) + BH
#   - 两段式分解:(a) 检出率 Fisher; (b) 检出者中强度 Wilcoxon
#   - BH family:任一组检出 ≥THRESH(非合并),避免误剔 HFpEF 特异通路;13 头条无条件测
#   - pool×disease 列联表核查(共线则提示 van Elteren 敏感性)
#
# 心内膜承载部分(UNC5/Netrin/FLRT)落笔:不在病人层面评估,由 Step 2 支撑(见讨论)。
###############################################################################

suppressPackageStartupMessages({
  library(zellkonverter); library(SingleCellExperiment); library(Matrix); library(CellChat)
})
options(stringsAsFactors = FALSE)
future::plan("sequential"); set.seed(1)

out_dir <- "/Users/franxy/Documents/博士课题 6.15/严谨/results"
h5      <- "/Users/franxy/Documents/博士课题 6.15/严谨/data/Hahn 人 HFpEF 单核 snRNA/HFpEF_snRNAseq.h5ad"

# 固定 4 类集(专家拍定)
FIXED_SET <- c("Cardiomyocyte","Fibroblast","Endothelial1","Pericyte")
MIN_NUC   <- 10          # 每类每病人最少核数
DETECT_THRESH <- 0.25    # BH family:任一组检出比例≥此值(起点,可据直方图调)

KEEP <- c(
  "CL:0000746" = "Cardiomyocyte", "CL:0002548" = "Fibroblast",
  "CL:0010008" = "Endothelial1",  "CL:4033076" = "Endothelial2",
  "CL:0000669" = "Pericyte",      "CL:0000359" = "VSMC",
  "CL:0000235" = "Macrophage",    "CL:0000542" = "Lymphocyte",
  "CL:0002350" = "Endocardial"
)
HEADLINE <- c("EPHA","EPHB","ADGRL","UNC5","Netrin","SEMA5","SLIT","NCAM","NRXN","FLRT",
              "LAMININ","COLLAGEN","FN1")

## ---- 读入 -----------------------------------------------------------------
message("== 读入 h5ad…")
sce <- readH5AD(h5, use_hdf5 = TRUE, reader = "R")
idx <- which(as.character(sce$cell_type) %in% names(KEEP))
sce <- sce[, idx]
celltype <- factor(KEEP[as.character(sce$cell_type)], levels = unname(KEEP))
donor    <- as.character(sce$donor_id)
pool     <- as.character(sce$biosample_id)
group    <- ifelse(sce$disease == "MONDO_0005252", "HFpEF", "control")
counts   <- as(assay(sce, "raw_count_cellranger"), "CsparseMatrix")
rownames(counts) <- rownames(sce); colnames(counts) <- colnames(sce)

donor_group <- tapply(group, donor, function(x) x[1])
donor_pool  <- tapply(pool,  donor, function(x) x[1])

## ---- pool × disease 列联表核查 -------------------------------------------
cat("\n========== pool × disease 列联表(混杂核查)==========\n")
pd_tab <- table(pool = donor_pool, disease = donor_group)
print(pd_tab)
collinear <- any(rowSums(pd_tab > 0) == 1)  # 某池只含单一组
cat(sprintf(">> 是否有池只含单一组(共线): %s\n", collinear))
if (collinear) cat("   → 建议补 van Elteren(按池分层 Wilcoxon)敏感性分析\n") else
  cat("   → 各池组别混合,pool 非主混杂,主检验可不分层\n")

## ===========================================================================
## 逐病人建对象(固定4类集),提取每通路信息流强度
## ===========================================================================
cat(sprintf("\n========== 逐病人 CellChat(固定%d类集, 每类≥%d核, nboot=1)==========\n",
            length(FIXED_SET), MIN_NUC))
donor_list <- names(donor_group)
strength_list <- list()
kept_donors <- character(0)

for (d in donor_list) {
  sel_d  <- donor == d
  lab_d  <- celltype[sel_d]
  ntab   <- table(lab_d)
  # 固定集要求:4 类都 ≥MIN_NUC,否则剔除该病人
  if (!all(FIXED_SET %in% names(ntab)) || any(ntab[FIXED_SET] < MIN_NUC)) {
    cat(sprintf("  病人 %s(%s):4类集未全达标,剔除\n", d, donor_group[d])); next
  }
  cells_d <- colnames(counts)[sel_d & celltype %in% FIXED_SET]
  lab_keep<- droplevels(celltype[sel_d & celltype %in% FIXED_SET])
  lab_keep<- factor(as.character(lab_keep), levels = FIXED_SET)  # 锁定顺序一致
  mat_d   <- counts[, cells_d, drop = FALSE]

  data_n <- normalizeData(mat_d, scale.factor = 1e4, do.log = TRUE)
  meta_d <- data.frame(labels = lab_keep, samples = factor("s1"), row.names = colnames(mat_d))
  cc <- tryCatch({
    cc <- createCellChat(object = data_n, meta = meta_d, group.by = "labels")
    cc@DB <- CellChatDB.human
    cc <- subsetData(cc); cc <- identifyOverExpressedGenes(cc)
    cc <- identifyOverExpressedInteractions(cc)
    cc <- computeCommunProb(cc, type="truncatedMean", trim=0.1, population.size=TRUE,
                            seed.use=1, nboot=1, raw.use=TRUE)
    cc <- filterCommunication(cc, min.cells = MIN_NUC)
    cc <- computeCommunProbPathway(cc); cc
  }, error = function(e){ cat(sprintf("  病人 %s 报错:%s\n",d,conditionMessage(e))); NULL })
  if (is.null(cc)) next

  prob <- cc@netP$prob
  if (is.null(prob) || length(dim(prob)) != 3) { next }
  strength <- apply(prob, 3, sum); names(strength) <- dimnames(prob)[[3]]
  strength_list[[d]] <- strength
  kept_donors <- c(kept_donors, d)
  cat(sprintf("  病人 %s(%s):保留,检出 %d 通路\n", d, donor_group[d], length(strength)))
}

cat(sprintf("\n保留病人: %d(HFpEF %d, control %d)\n",
            length(kept_donors),
            sum(donor_group[kept_donors]=="HFpEF"),
            sum(donor_group[kept_donors]=="control")))

## ---- 构造 通路 × 病人 强度矩阵(缺失=0)-----------------------------------
all_pw <- sort(unique(unlist(lapply(strength_list, names))))
S <- matrix(0, length(all_pw), length(kept_donors), dimnames=list(all_pw, kept_donors))
for (d in kept_donors) S[names(strength_list[[d]]), d] <- strength_list[[d]]

grp <- donor_group[kept_donors]
hf <- kept_donors[grp=="HFpEF"]; ct <- kept_donors[grp=="control"]

## ---- 相对流量版:每病人每通路 ÷ 该病人总流量 ------------------------------
donor_total <- colSums(S)
R <- sweep(S, 2, donor_total, "/"); R[is.na(R)] <- 0

# 混杂核查:每病人总流量按组分布
cat("\n========== 每病人总流量按组分布(混杂核查)==========\n")
cat("HFpEF 总流量:\n");  print(summary(donor_total[hf]))
cat("control 总流量:\n");print(summary(donor_total[ct]))
cat(sprintf("总流量组间 Wilcoxon p = %.3f(若显著,绝对流量增强可能部分是总量差异)\n",
            suppressWarnings(wilcox.test(donor_total[hf], donor_total[ct])$p.value)))

## ===========================================================================
## 检验函数:对一个 通路×病人 矩阵,出每通路 Wilcoxon-on-full + 两段式
##   do_fisher=TRUE 才算检出率Fisher(仅绝对版需要;相对版非零结构与绝对版相同,
##   Fisher会逐字重复,故相对版关掉避免冒充独立证据)
## ===========================================================================
# 安全检验包装:退化输入(全0/样本太少)返回 NA 而非中断整批
safe_wilcox <- function(x, y) tryCatch(
  suppressWarnings(wilcox.test(x, y)$p.value), error = function(e) NA_real_)
safe_HL <- function(x, y) tryCatch(
  suppressWarnings(wilcox.test(x, y, conf.int = TRUE)$estimate), error = function(e) NA_real_)
safe_fisher <- function(tab) tryCatch(
  suppressWarnings(fisher.test(tab)$p.value), error = function(e) NA_real_)

test_matrix <- function(M, hf, ct, label, do_fisher = TRUE) {
  pw <- rownames(M)
  # 先算好命名向量(按 pathway 名索引正确;修复 res[g,] 行名默认"1,2,.."的bug)
  det_h  <- rowSums(M[, hf, drop=FALSE] > 0)        # 命名向量(名=pathway)
  det_c  <- rowSums(M[, ct, drop=FALSE] > 0)
  mean_h <- rowMeans(M[, hf, drop=FALSE])
  mean_c <- rowMeans(M[, ct, drop=FALSE])
  med_h  <- apply(M[, hf, drop=FALSE], 1, median)
  med_c  <- apply(M[, ct, drop=FALSE], 1, median)

  res <- data.frame(pathway = pw, stringsAsFactors = FALSE)
  res$mean_HFpEF     <- mean_h[pw]
  res$mean_control   <- mean_c[pw]
  res$median_HFpEF   <- med_h[pw]      # 展示用,不定方向
  res$median_control <- med_c[pw]
  # 方向用【均值差】定(对零膨胀位移敏感、不会塌成0);中位数仅展示
  res$direction <- ifelse(res$mean_HFpEF > res$mean_control, "HFpEF_up", "HFpEF_down")
  # Hodges-Lehmann 位移估计(与 Wilcoxon p 同源,辅助佐证方向)
  res$HL_shift <- sapply(pw, function(g) safe_HL(M[g, hf], M[g, ct]))
  # 主检验:Wilcoxon-on-full(含0)
  res$p_wilcox <- sapply(pw, function(g) safe_wilcox(M[g, hf], M[g, ct]))
  # 检出率
  res$det_HFpEF    <- det_h[pw]
  res$det_control  <- det_c[pw]
  res$detrate_HFpEF   <- det_h[pw] / length(hf)
  res$detrate_control <- det_c[pw] / length(ct)
  # 两段式 (a):检出率 Fisher(仅 do_fisher 时;用命名向量按名索引)
  if (do_fisher) {
    res$p_fisher_detect <- sapply(pw, function(g) {
      tab <- matrix(c(det_h[g], length(hf)-det_h[g],
                      det_c[g], length(ct)-det_c[g]), nrow = 2)
      safe_fisher(tab)
    })
  } else {
    res$p_fisher_detect <- NA_real_
  }
  # 两段式 (b):检出者中强度 Wilcoxon
  res$p_wilcox_detectors <- sapply(pw, function(g) {
    hv <- M[g, hf][M[g, hf] > 0]; cv <- M[g, ct][M[g, ct] > 0]
    if (length(hv) < 3 || length(cv) < 3) return(NA_real_)
    safe_wilcox(hv, cv)
  })
  rownames(res) <- NULL
  attr(res, "label") <- label
  res
}

## ---- BH family:任一组检出 ≥THRESH(非合并)--------------------------------
detrate_hf <- rowSums(S[, hf, drop=FALSE] > 0) / length(hf)
detrate_ct <- rowSums(S[, ct, drop=FALSE] > 0) / length(ct)
in_family <- (detrate_hf >= DETECT_THRESH) | (detrate_ct >= DETECT_THRESH)
family_pw <- rownames(S)[in_family]
# 13 头条无条件并入 family
family_pw <- union(family_pw, intersect(HEADLINE, rownames(S)))
cat(sprintf("\n========== BH family ==========\n"))
cat(sprintf("任一组检出≥%.0f%% → %d 条;并入13头条后 family = %d 条\n",
            DETECT_THRESH*100, sum(in_family), length(family_pw)))

## ---- 全 169 通路检出直方图(定/核验门槛)----------------------------------
cat("\n========== 全通路检出分布(找自然断点)==========\n")
det_all <- rowSums(S > 0)
cat("各检出病人数的通路计数(检出数:通路数):\n")
print(table(det_all))

## ===========================================================================
## 出结果:绝对流量(主) + 相对流量(co-primary),都在 family 上做 BH
## ===========================================================================
run_and_report <- function(M, label, do_fisher = TRUE) {
  res <- test_matrix(M, hf, ct, label, do_fisher = do_fisher)
  res_fam <- res[res$pathway %in% family_pw, ]
  res_fam$p_wilcox_BH <- p.adjust(res_fam$p_wilcox, "BH")
  if (do_fisher) res_fam$p_fisher_BH <- p.adjust(res_fam$p_fisher_detect, "BH")
  # 合回完整表(family 外的 BH 记 NA)
  res$p_wilcox_BH <- res_fam$p_wilcox_BH[match(res$pathway, res_fam$pathway)]
  res$p_fisher_BH <- if (do_fisher) res_fam$p_fisher_BH[match(res$pathway, res_fam$pathway)] else NA_real_

  cat(sprintf("\n========== 结果:%s ==========\n", label))
  cat("--- 13 头条通路(无条件测;方向由均值差定,中位数仅展示)---\n")
  hl <- res[match(HEADLINE, res$pathway), ]
  hl <- hl[!is.na(hl$pathway), ]
  show_cols <- c("pathway","direction","mean_HFpEF","mean_control","HL_shift",
                 "p_wilcox","p_wilcox_BH","detrate_HFpEF","detrate_control")
  if (do_fisher) show_cols <- c(show_cols, "p_fisher_BH")
  print(hl[, show_cols], row.names = FALSE, digits = 3)
  res
}

cat("\n############ 绝对流量版(主检验,population.size=TRUE)############")
res_abs <- run_and_report(S, "绝对流量")
write.csv(res_abs, file.path(out_dir, "perdonor_test_absolute.csv"), row.names = FALSE)

cat("\n############ 相对流量版(co-primary,占比;Fisher与绝对版相同故不重算)############")
res_rel <- run_and_report(R, "相对流量", do_fisher = FALSE)
write.csv(res_rel, file.path(out_dir, "perdonor_test_relative.csv"), row.names = FALSE)

## ---- 保存全部 ------------------------------------------------------------
saveRDS(list(S=S, R=R, donor_group=donor_group[kept_donors], donor_pool=donor_pool[kept_donors],
             kept_donors=kept_donors, family_pw=family_pw,
             res_abs=res_abs, res_rel=res_rel, pd_tab=pd_tab),
        file.path(out_dir, "perdonor_final.rds"))

cat("\n========== 解读指引 ==========\n")
cat("- 主结论看【绝对流量版】头条:p_wilcox_BH<0.05 且 direction=HFpEF_up = 病人层面显著增强\n")
cat("- 与【相对流量版】对照:两版都成立 → 增强是通路特异、非总量/丰度假象\n")
cat("- 两段式:p_fisher_BH(检出率)说明增强是否靠'更多病人检出'(如 NRXN);\n")
cat("  检出者强度 p_wilcox_detectors 说明检出者中是否更强\n")
cat("- UNC5/Netrin/FLRT 的心内膜承载部分不在此评估(4类集),由 Step 2 支撑\n")
cat("- 总流量组间若不显著 → 绝对流量增强非总量驱动\n")

cat("\n========== sessionInfo ==========\n")
print(sessionInfo())
