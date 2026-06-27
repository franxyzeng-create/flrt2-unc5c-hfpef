###############################################################################
# 21_bulk_BMI_sensitivity.R —— FLRT2 的 BMI 去混杂(obesity sensitivity,投稿 Supplementary)
#
# 背景(专家内容审 Q3):队列为肥胖型 HFpEF(BMI 40.8 vs Normal/HFrEF ~26),与肥胖近共线;
#   验证 bulk(Hahn 2021)原文报"BMI 在很大程度上解释 HFpEF 上调基因",FLRT2 属上调类。
#   本脚本把"FLRT2 升高是否单纯肥胖"做成正式去混杂(crude 版结论:非纯肥胖,见下)。
#
# 三个分析:
#   (A) DESeq2 ~ BMI + condition:BMI 校正后 HFpEF-vs-Normal 的 FLRT2 logFC/padj(主)
#   (B) 非 HFpEF 内(Normal+HFrEF)FLRT2 ~ BMI 相关:绕开 condition×BMI 共线
#   (C) HFrEF(非肥胖)vs Normal 的 FLRT2:非肥胖心衰是否也升
#
# crude 预期(本机 log2CPM 实算,作对照):HFpEF vs Normal +0.89(p~7e-7);
#   非HFpEF FLRT2~BMI ρ≈-0.05(p≈0.72);HFrEF 也升(med 2.26 vs Normal 1.45,p~9e-8);
#   全样本含HFpEF才正相关 ρ≈+0.21 = 共线假象。
#
# 数据:严谨/data/Hahn2020_HFpEF_bulk/  (RVS 原始 count + 临床表;FLRT2=ENSG00000185070)
###############################################################################

suppressPackageStartupMessages({ library(readxl); library(DESeq2) })

BASE <- "/Users/franxy/Documents/博士课题 6.15/严谨/data/Hahn2020_HFpEF_bulk"
CNT  <- file.path(BASE, "baderzone-HFpEF_2020-71fc3a4/Data/RawReads_RVS_Ensembl_geneName.xlsx")
CLIN <- file.path(BASE, "ClinicalCharacteristicsRNAseqSamples.csv")
FLRT2 <- "ENSG00000185070"

## ---- 读入 + 对齐 ----------------------------------------------------------
cnt <- as.data.frame(read_excel(CNT)); rownames(cnt) <- cnt[[1]]; cnt[[1]] <- NULL
cnt <- as.matrix(round(cnt)); mode(cnt) <- "integer"          # count 矩阵(基因×样本)
colnames(cnt) <- as.character(colnames(cnt))

cl <- read.csv(CLIN, check.names = FALSE)
cl <- cl[cl[["Tissue Abbr"]] == "RVS", ]
cl$sid <- as.character(cl[["SequenceRunId"]])
cl$condition <- factor(cl[["Disease Abbr"]], levels = c("Normal","HFpEF","HFrEF"))
cl$BMI <- suppressWarnings(as.numeric(cl[["BMI, kg/m2"]]))

common <- intersect(colnames(cnt), cl$sid)
cl <- cl[match(common, cl$sid), ]; cnt <- cnt[, common]
cat(sprintf("matched RVS 样本: %d | 组别: %s\n", length(common),
            paste(names(table(cl$condition)), table(cl$condition), sep="=", collapse=", ")))
cat(sprintf("BMI median: Normal %.1f / HFpEF %.1f / HFrEF %.1f\n",
    median(cl$BMI[cl$condition=="Normal"],na.rm=T),
    median(cl$BMI[cl$condition=="HFpEF"],na.rm=T),
    median(cl$BMI[cl$condition=="HFrEF"],na.rm=T)))

## ---- (A) DESeq2 ~ BMI + condition(BMI 校正)-------------------------------
keep_bmi <- !is.na(cl$BMI)
dds <- DESeqDataSetFromMatrix(cnt[, keep_bmi],
        colData = data.frame(condition = droplevels(cl$condition[keep_bmi]),
                             BMIc = as.numeric(scale(cl$BMI[keep_bmi]))),  # 中心化+标准化:消 DESeq2 collinearity 警告(不改 condition 对比)
        design = ~ BMIc + condition)            # condition 放最后 → results 取其对比
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
res_adj <- results(dds, contrast = c("condition","HFpEF","Normal"))   # BMI 校正后
res_un  <- {
  dds2 <- DESeqDataSetFromMatrix(cnt[, keep_bmi],
            colData = data.frame(condition = droplevels(cl$condition[keep_bmi])),
            design = ~ condition); dds2 <- dds2[rowSums(counts(dds2))>=10,]; dds2 <- DESeq(dds2)
  results(dds2, contrast = c("condition","HFpEF","Normal"))           # 未校正
}
cat("\n=== (A) FLRT2 HFpEF vs Normal:未校正 vs BMI 校正(看 log2FC 是否收缩、padj)===\n")
for (tag in c("未校正","BMI校正")) {
  r <- if (tag=="未校正") res_un else res_adj
  if (FLRT2 %in% rownames(r)) {
    x <- r[FLRT2, ]
    cat(sprintf("  %s: log2FC=%+.3f, lfcSE=%.3f, padj=%.2g\n", tag, x$log2FoldChange, x$lfcSE, x$padj))
  } else cat("  [!] FLRT2 不在该 dds(可能被 filterByExpr 滤掉);改看 (B)/(C)\n")
}
cat("  判读:log2FC 校正后若只收缩、方向/显著性保持 → 非纯肥胖;若塌向 0 → BMI 解释力强(共线致功效有限,需 BMI-matched 子集)\n")

## ---- (B) 非 HFpEF 内 FLRT2 ~ BMI(绕开共线)-------------------------------
vsd <- assay(varianceStabilizingTransformation(dds, blind = TRUE))
f2  <- vsd[FLRT2, ]
sub <- droplevels(cl$condition[keep_bmi]) != "HFpEF"
ct  <- suppressWarnings(cor.test(f2[sub], cl$BMI[keep_bmi][sub], method = "spearman"))
cat(sprintf("\n=== (B) 非 HFpEF(Normal+HFrEF,n=%d)FLRT2(VST)~BMI:Spearman ρ=%+.2f, p=%.2g ===\n",
            sum(sub), ct$estimate, ct$p.value))
cat("  判读:ρ≈0/NS → 剔除 HFpEF 身份后 FLRT2 不随 BMI(反纯肥胖基因)\n")

## ---- (C) HFrEF(非肥胖)vs Normal:非肥胖心衰是否也升 -----------------------
cat("\n=== (C) FLRT2 各组 VST 中位 ===\n")
for (g in c("Normal","HFpEF","HFrEF")) {
  idx <- droplevels(cl$condition[keep_bmi]) == g
  cat(sprintf("  %s: %.2f\n", g, median(f2[idx])))
}
cat("  判读:HFrEF(BMI~26)若≈HFpEF、远高于 Normal → 非肥胖心衰也升,进一步反纯肥胖\n")

## ---- 备选:BMI-matched 子集(最干净,若样本量够)---------------------------
# 取 HFpEF 中 BMI 较低者 + Normal/HFrEF 中 BMI 较高者,做 BMI 分布重叠的子集再比 FLRT2。
# 本队列 HFpEF(40.8)与对照(26)几乎不重叠 → matched 子集可能样本太少;先看 (A)-(C)。

cat("\n== sessionInfo ==\n"); print(sessionInfo())
