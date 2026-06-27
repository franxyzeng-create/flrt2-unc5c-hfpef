###############################################################################
# Step 2b-2 前置验证(v2):CellBender 臂 slot 行为四查
#
# v2 修正:查①查③不再裸比 @x 向量(稀疏矩阵重排/丢全零基因会让 @x 顺序变→假失败),
#         改成按基因名+细胞名对齐后比【整个矩阵】,免疫重排。
#
# 四查:
#   ① createCellChat 后 @data[g,c] == cb_layer[g,c](按名对齐,没被再处理)
#   ② @data.raw 为空(num[0,0])
#   ③ subsetData 后 @data.signaling[g,c] == cb_layer[g,c](只子集不变换)
#   ④ computeCommunProb 显式 raw.use=TRUE 正常跑
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

## 查①:@data[g,c] == cb[g,c](按名对齐比矩阵,免疫重排)
cat("\n--- 查①:@data == 传入 cb_layer(按名对齐)? ---\n")
g <- rownames(cb); cc_data <- cc@data
ok_names <- all(g %in% rownames(cc_data)) && all(colnames(cb) %in% colnames(cc_data))
cat("所有基因/细胞名都在 @data 里:", ok_names, "\n")
q1 <- ok_names && isTRUE(all.equal(cb[g, colnames(cb)], cc_data[g, colnames(cb)]))
cat(">> @data[g,c] 与 cb[g,c] 完全相同:", q1, "(TRUE=没被再归一化 ✓)\n")

## 查②:@data.raw 是否为空
cat("\n--- 查②:@data.raw 是否为空 ---\n")
draw_dim <- dim(cc@data.raw)
cat("@data.raw 维度:", paste(draw_dim, collapse=" × "), "\n")
q2 <- (is.null(cc@data.raw) || prod(draw_dim) == 0)
cat(">> @data.raw 为空:", q2, "(预期 TRUE)\n")

## 查③:subsetData 后 @data.signaling[g,c] == cb[g,c](按名对齐)
cat("\n--- 查③:subsetData 后 @data.signaling 只子集不变换 ---\n")
cc@DB <- CellChatDB.human
cc <- subsetData(cc)
ds <- cc@data.signaling
cat("@data.signaling 维度:", paste(dim(ds), collapse=" × "), "\n")
sig_genes <- rownames(ds)
q3 <- all(sig_genes %in% rownames(cb)) &&
      isTRUE(all.equal(cb[sig_genes, colnames(ds)], ds[sig_genes, colnames(ds)]))
cat(">> @data.signaling[g,c] == cb[g,c] 子集:", q3, "(TRUE=只子集、数值没变 ✓)\n")

## 查④:computeCommunProb 显式 raw.use=TRUE 能跑
cat("\n--- 查④:computeCommunProb(raw.use=TRUE) 试跑 ---\n")
cc <- identifyOverExpressedGenes(cc)
cc <- identifyOverExpressedInteractions(cc)
cc_test <- tryCatch({
  computeCommunProb(cc, type = "truncatedMean", trim = 0.1,
                    population.size = TRUE, seed.use = 1, nboot = 100, raw.use = TRUE)
}, error = function(e) { cat("报错:", conditionMessage(e), "\n"); NULL })
q4 <- !is.null(cc_test)
cat(">> computeCommunProb(raw.use=TRUE) 正常完成:", q4, "\n")

## 总结
cat("\n========== 四查总结 ==========\n")
cat(sprintf("查① @data==cb(按名对齐,没再归一化): %s\n", q1))
cat(sprintf("查② @data.raw 为空:                  %s\n", q2))
cat(sprintf("查③ @data.signaling 只子集不变:      %s\n", q3))
cat(sprintf("查④ raw.use=TRUE 正常跑:             %s\n", q4))
if (q1 && q3 && q4) {
  cat("\n>> 四查通过:(a)臂实现正确,cb layer 全程零再归一化,可以跑完整三臂(08 v3)。\n")
} else {
  cat("\n>> 有查未通过,把上面输出贴回再判。\n")
}
