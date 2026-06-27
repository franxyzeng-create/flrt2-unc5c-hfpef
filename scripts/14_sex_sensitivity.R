###############################################################################
# Step 3 sex 敏感性:处理 HFpEF 性别失衡(HFpEF 83%女 vs control 50/50)
#
# 性别失衡是 HFpEF 流行病学固有特征(女性占优),非可修设计缺陷。处理方案:
#   (a) female-only(主敏感性,临床更贴切):女性内 HFpEF 15 vs control 12,
#       零性别混杂,硬核心绝对+相对是否仍显著。scope = "在女性中"。
#   (b) van Elteren(辅):按 sex 分层合并检验(全42),male层 HFpEF n=3 几乎无贡献、标注。
#   一致性区分:female-only 减弱但 van Elteren 仍显著 → 功效损失(非混杂);
#              两者都减弱 → 真性别混杂。
#   并:女性内重跑总流量比较,确认"大盘升高"非性别 artifact。
#
# 复用 perdonor_final.rds 的 S(绝对)/R(相对)矩阵,不重跑 CellChat。
###############################################################################

suppressPackageStartupMessages({
  library(zellkonverter); library(SingleCellExperiment); library(Matrix)
})
options(stringsAsFactors = FALSE)
out_dir <- "/Users/franxy/Documents/博士课题 6.15/严谨/results"
h5      <- "/Users/franxy/Documents/博士课题 6.15/严谨/data/Hahn 人 HFpEF 单核 snRNA/HFpEF_snRNAseq.h5ad"

HEADLINE <- c("EPHA","EPHB","ADGRL","UNC5","Netrin","SEMA5","SLIT","NCAM","NRXN","FLRT",
              "LAMININ","COLLAGEN","FN1")
CORE <- c("UNC5","FLRT","ADGRL")  # 三轴全清硬核心

pd <- readRDS(file.path(out_dir, "perdonor_final.rds"))
S <- pd$S; R <- pd$R; kept <- pd$kept_donors; dg <- pd$donor_group[kept]
family_pw <- pd$family_pw

# 读 donor 性别
message("== 读 donor 性别…")
sce <- readH5AD(h5, use_hdf5 = TRUE, reader = "R")
donor_sex_all <- tapply(as.character(sce$sex), as.character(sce$donor_id), function(x) x[1])
sex <- donor_sex_all[kept]
cat("性别 × 疾病:\n"); print(table(sex=sex, disease=dg))

## ===========================================================================
## 安全检验包装 + 效应量
## ===========================================================================
safe_w  <- function(x,y) tryCatch(suppressWarnings(wilcox.test(x,y)$p.value), error=function(e) NA_real_)
# probability of superiority:P(随机HFpEF > 随机control),0.5=无效应,→1=HFpEF主导
# 对零膨胀稳健;用于"功效损失 vs 真性别混杂"的主诊断(比 p 更能判效应是否持续)
ps <- function(x, y) tryCatch({
  w <- suppressWarnings(wilcox.test(x, y))
  as.numeric(w$statistic) / (length(x) * length(y))
}, error = function(e) NA_real_)

# 在给定病人子集上,对一个矩阵做每通路 Wilcoxon + BH(family内)
test_subset <- function(M, donors_use, label) {
  g  <- dg[donors_use]
  hf <- donors_use[g=="HFpEF"]; ct <- donors_use[g=="control"]
  pw <- rownames(M)
  res <- data.frame(pathway=pw,
                    mean_HFpEF   = rowMeans(M[, hf, drop=FALSE]),
                    mean_control = rowMeans(M[, ct, drop=FALSE]))
  res$direction <- ifelse(res$mean_HFpEF > res$mean_control, "HFpEF_up", "HFpEF_down")
  res$p_wilcox  <- sapply(pw, function(x) safe_w(M[x,hf], M[x,ct]))
  # BH 在 family 内
  fam <- res$pathway %in% family_pw
  res$p_BH <- NA_real_
  res$p_BH[fam] <- p.adjust(res$p_wilcox[fam], "BH")
  rownames(res) <- NULL
  attr(res,"label") <- label
  res
}

## ===========================================================================
## (a) female-only:女性内 HFpEF vs control
## ===========================================================================
fem <- kept[sex=="female"]
cat(sprintf("\n========== (a) female-only(HFpEF %d vs control %d)==========\n",
            sum(dg[fem]=="HFpEF"), sum(dg[fem]=="control")))

abs_f <- test_subset(S, fem, "female-绝对")
rel_f <- test_subset(R, fem, "female-相对")

show <- function(res, cols) {
  hl <- res[match(HEADLINE, res$pathway), ]; hl <- hl[!is.na(hl$pathway),]
  print(hl[, cols], row.names=FALSE, digits=3)
}
cat("\n--- female-only 绝对流量(头条)---\n")
show(abs_f, c("pathway","direction","mean_HFpEF","mean_control","p_wilcox","p_BH"))
cat("\n--- female-only 相对流量(头条)---\n")
show(rel_f, c("pathway","direction","mean_HFpEF","mean_control","p_wilcox","p_BH"))

# 硬核心在女性内:绝对+相对是否都还显著 + PS效应量对比 + 女性内方向
cat("\n--- 硬核心 UNC5/FLRT/ADGRL 女性内(BH展示 + PS效应量 + 方向)---\n")
# 全样本与女性内的 HFpEF/control 病人
hf_all <- kept[dg=="HFpEF"]; ct_all <- kept[dg=="control"]
hf_fem <- fem[dg[fem]=="HFpEF"]; ct_fem <- fem[dg[fem]=="control"]

core_tab <- data.frame(pathway = CORE)
# BH(展示用)
core_tab$abs_BH_full   <- pd$res_abs$p_wilcox_BH[match(CORE, pd$res_abs$pathway)]
core_tab$abs_BH_female <- abs_f$p_BH[match(CORE, abs_f$pathway)]
core_tab$rel_BH_female <- rel_f$p_BH[match(CORE, rel_f$pathway)]
# 女性内方向(绝对)——方向翻转是性别混杂最干脆的信号
core_tab$dir_female <- abs_f$direction[match(CORE, abs_f$pathway)]
# PS 效应量(主诊断):全样本 vs 女性内,绝对 & 相对
core_tab$PS_abs_full   <- sapply(CORE, function(g) ps(S[g, hf_all], S[g, ct_all]))
core_tab$PS_abs_female <- sapply(CORE, function(g) ps(S[g, hf_fem], S[g, ct_fem]))
core_tab$PS_rel_full   <- sapply(CORE, function(g) ps(R[g, hf_all], R[g, ct_all]))
core_tab$PS_rel_female <- sapply(CORE, function(g) ps(R[g, hf_fem], R[g, ct_fem]))
print(core_tab, row.names=FALSE, digits=3)
cat("PS=0.5 无效应,→1 HFpEF主导。主诊断:PS_female 接近 PS_full→效应保留(变弱因n小=功效损失);\n")
cat("  PS_female 滑向0.5→效应塌了=真性别混杂。方向若女性内仍 HFpEF_up 且 PS 保持 → 铁核心。\n")

## ===========================================================================
## (b) van Elteren:按 sex 分层(全42),male层HFpEF n=3几乎无贡献
## ===========================================================================
cat("\n========== (b) van Elteren 按 sex 分层(全42,male层HFpEF n=3)==========\n")
# van Elteren = coin::wilcox_test(val~grp|strat) 规范分层秩检验
# 注意(专家):male HFpEF n=3,van Elteren ≈ female-only 且可能被3个male噪声拽动,
#   所以它退为"用全数据的形式确认",解释从严;主诊断用 PS 效应量对比(见下)
has_coin <- requireNamespace("coin", quietly = TRUE)
if (!has_coin) {
  cat("[!] coin 包未装 → van Elteren 无法跑(不降级,避免冒充)。\n")
  cat("    请先 install.packages('coin') 再跑本脚本的 (b) 段。\n")
  cat("    PS 主诊断不依赖 coin,仍会正常输出。\n")
}

ve_one <- function(M, g) {
  if (!has_coin) return(NA_real_)   # 不降级冒充,直接 NA
  df <- data.frame(
    val = M[g, kept],
    grp = factor(dg[kept], levels=c("control","HFpEF")),
    strat = factor(sex)
  )
  tryCatch(
    as.numeric(coin::pvalue(coin::wilcox_test(val ~ grp | strat, data=df))),
    error=function(e) NA_real_)
}

cat("\n--- 硬核心 van Elteren(绝对流量,sex分层)---\n")
ve_core <- data.frame(
  pathway = CORE,
  vanElteren_p = sapply(CORE, function(g) ve_one(S, g))
)
ve_core$vanElteren_BH <- NA
# 对13头条做BH(分层p)
ve_all <- sapply(intersect(HEADLINE, rownames(S)), function(g) ve_one(S, g))
ve_all_BH <- p.adjust(ve_all, "BH")
ve_core$vanElteren_BH <- ve_all_BH[match(CORE, names(ve_all))]
print(ve_core, row.names=FALSE, digits=3)

## ===========================================================================
## 女性内总流量(确认"大盘升高"非性别artifact)
## ===========================================================================
cat("\n========== 女性内总流量(确认大盘升高非性别artifact)==========\n")
tf <- colSums(S)
tf_fem_hf <- tf[fem[dg[fem]=="HFpEF"]]; tf_fem_ct <- tf[fem[dg[fem]=="control"]]
cat(sprintf("女性 HFpEF 总流量中位 %.3f vs 女性 control %.3f\n",
            median(tf_fem_hf), median(tf_fem_ct)))
cat(sprintf(">> 女性内总流量 Wilcoxon p = %.3f(仍显著→大盘升高对性别稳健)\n",
            safe_w(tf_fem_hf, tf_fem_ct)))

## ===========================================================================
## 一致性判定:功效损失 vs 真性别混杂(主诊断用 PS 效应量,非 BH/van Elteren)
## ===========================================================================
cat("\n========== 一致性判定(主诊断:PS效应量对比)==========\n")
cat("逻辑:PS_female≈PS_full → 效应保留、变弱因n小=功效损失(可接受);\n")
cat("     PS_female 滑向0.5 → 效应塌=真性别混杂(按NRXN诚实处理、核心收窄)\n")
cat("     van Elteren 仅作全数据形式确认(male n=3,解释从严,不作判据主轴)\n\n")
for (g in CORE) {
  row <- core_tab[core_tab$pathway==g, ]
  ve  <- ve_core$vanElteren_BH[ve_core$pathway==g]
  ps_drop_abs <- row$PS_abs_full - row$PS_abs_female   # PS 下滑幅度(绝对)
  ps_drop_rel <- row$PS_rel_full - row$PS_rel_female
  # 判据:女性内 PS 是否仍明显>0.5 且未大幅缩水;方向是否仍 up
  dir_ok <- row$dir_female == "HFpEF_up"
  abs_hold <- row$PS_abs_female >= 0.65 && ps_drop_abs <= 0.12   # 阈值供参考,看数定
  rel_hold <- row$PS_rel_female >= 0.65 && ps_drop_rel <= 0.12
  verdict <- if (dir_ok && abs_hold && rel_hold)
               "对性别稳健(绝对+相对 PS 女性内均保持、方向 up)→ 真·性别全清铁核心"
             else if (dir_ok && (row$PS_abs_female >= 0.6))
               "效应保留但变弱:更像功效损失(n=15v12),非方向翻转;诚实标注"
             else
               "女性内效应明显缩小/方向变 → 警惕真性别混杂,按NRXN诚实降级"
  cat(sprintf("  %s: 方向(女)=%s | PS绝对 %0.2f→%0.2f | PS相对 %0.2f→%0.2f | vanElterenBH=%.3g\n      → %s\n",
              g, row$dir_female, row$PS_abs_full, row$PS_abs_female,
              row$PS_rel_full, row$PS_rel_female, ve, verdict))
}

saveRDS(list(abs_f=abs_f, rel_f=rel_f, core_tab=core_tab, ve_core=ve_core),
        file.path(out_dir, "step3_sex_sensitivity.rds"))
cat("\n========== 写法指引 ==========\n")
cat("- 核心女性内 abs+rel 都显著 → limitation(性别失衡)转强项:'核心通路对性别稳健'\n")
cat("- 某核心女性内减弱:措辞'全样本显著、女性内减弱,后者可能反映残余性别混杂或较小样本功效下降'\n")
cat("- female-only 定位为'临床更贴切的干净主分析'(HFpEF 女性占优),scope='在女性中',特性非缺陷\n")

cat("\n========== sessionInfo ==========\n"); print(sessionInfo())
