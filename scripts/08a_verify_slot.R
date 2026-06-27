###############################################################################
# Step 2b-2 前置验证:CellBender 臂(a)的 slot 行为四查
#
# 专家源码钉死:createCellChat(object=已归一化矩阵) 不会再归一化,
#   链条 cb_layer → @data →(subsetData)→ @data.signaling → computeCommunProb 零变换。
# 本脚本用四查确认这条链在你的 2.2.0.9001 build 上确实成立,再跑三臂。
#
# 四查:
#   ①  createCellChat 后 @data 非零值 == 传入 cb_layer(没被再处理)
#   ②  @data.raw 为空(num[0,0])
#   ③  subsetData 后 @data.signaling == cb_layer 取信号基因子集(只子集不变换)
#   ④  computeCommunProb 必须显式 raw.use=TRUE(默认值跨版本变过,pin死)
###############################################################################

suppressPackageStartupMessages({
  library(CellChat); library(zellkonverter); library(SingleCellExperiment); library(Matrix)
})

h5 <- "/Users/franxy/Documents/博士课题 6.15/严谨/data/Hahn 人 HFpEF 单核 snRNA/HFpEF_snRNAseq.h5ad"
message("== 读入 h5ad,取 cb layer 前 500 细胞做验证…")
sce <- readH5AD(h5, use_hdf5 = TRUE, reader = "R")
cb <- as(assay(sce, "raw_count_cellbender"), "CsparseMatrix")[, 1:500]
rownames(cb) <- rownames(sce); colnames(cb) <- colnames(sce)[1:500]

meta <- data.frame(labels = factor(rep(c("A","B"), 250)),
                   samples = factor("s1"), row.names = colnames(cb))

cat("\n========== 建对象 ==========\n")
cc <- createCellChat(object = cb, meta = meta, group.by = "labels")

## 查①:@data 是否 == 传入 cb(没被再归一化)
cat("\n--- 查①:@data == 传入 cb_layer? ---\n")
d_data <- cc@data
cat("传入 cb 前8非零值:  ", round(cb@x[1:8], 4), "\n")
cat("@data    前8非零值:  ", round(d_data@x[1:8], 4), "\n")
q1 <- isTRUE(all.equal(as.numeric(cb@x), as.numeric(d_data@x)))
cat(">> @data 与 cb 完全相同:", q1, "(TRUE=没被再归一化 ✓)\n")

## 查②:@data.raw 是否为空
cat("\n--- 查②:@data.raw 是否为空 ---\n")
draw_dim <- dim(cc@data.raw)
cat("@data.raw 维度:", paste(draw_dim, collapse=" × "), "\n")
q2 <- (is.null(cc@data.raw) || prod(draw_dim) == 0)
cat(">> @data.raw 为空:", q2, "(预期 TRUE)\n")

## 查③:subsetData 后 @data.signaling == cb 取信号基因子集
cat("\n--- 查③:subsetData 后 @data.signaling 是否只子集不变换 ---\n")
cc@DB <- CellChatDB.human
cc <- subsetData(cc)
ds <- cc@data.signaling
cat("@data.signaling 维度:", paste(dim(ds), collapse=" × "), "\n")
# 取信号基因子集后,这些基因在 cb 里的值应与 @data.signaling 完全一致
common_genes <- intersect(rownames(ds), rownames(cb))
cb_sub <- cb[common_genes, , drop=FALSE]
ds_sub <- ds[common_genes, , drop=FALSE]
q3 <- isTRUE(all.equal(as.numeric(cb_sub@x), as.numeric(ds_sub@x)))
cat(">> @data.signaling == cb 对应基因子集:", q3, "(TRUE=只子集、数值没变 ✓)\n")

## 查④:computeCommunProb 显式 raw.use=TRUE 能正常跑(不报错)
cat("\n--- 查④:computeCommunProb(raw.use=TRUE) 试跑 ---\n")
cc <- identifyOverExpressedGenes(cc)
cc <- identifyOverExpressedInteractions(cc)
cc_test <- tryCatch({
  computeCommunProb(cc, type = "truncatedMean", trim = 0.1,
                    population.size = TRUE, seed.use = 1, nboot = 100,
                    raw.use = TRUE)
}, error = function(e) { cat("computeCommunProb 报错:", conditionMessage(e), "\n"); NULL })
q4 <- !is.null(cc_test)
cat(">> computeCommunProb(raw.use=TRUE) 正常完成:", q4, "\n")

## 总结
cat("\n========== 四查总结 ==========\n")
cat(sprintf("查① @data==cb(没再归一化):     %s\n", q1))
cat(sprintf("查② @data.raw 为空:             %s\n", q2))
cat(sprintf("查③ @data.signaling 只子集不变: %s\n", q3))
cat(sprintf("查④ raw.use=TRUE 正常跑:        %s\n", q4))
if (q1 && q3 && q4) {
  cat("\n>> 四查通过:(a)臂实现正确,cb layer 全程零再归一化,可以跑完整三臂。\n")
} else {
  cat("\n>> 有查未通过,把上面输出贴回,先修 (a) 臂再跑三臂。\n")
}
