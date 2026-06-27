###############################################################################
# 第 1 步(出版前必做):可复现性存档 + p 值正确归因 + BH 校正
#
# 依据两位专家对 jinworks/CellChat 源码的审查:
#   1. computeCommunProb 在主进程 set.seed(seed.use) 一次、一次性生成置换矩阵,
#      之后并行 worker 只做确定性算术 → net$pval 跨后端/worker 数逐位可复现。
#      本脚本原样重跑一次并 identical() 断言,把"架构可复现"变成可引用证据。
#   2. computeCommunProb 的 permutation p 值 = nReject/nboot,nboot=100 时
#      分辨率 0.01 → 只能报 "p<0.05",绝不可能是 1e-8。
#   3. 那些 1e-8/1e-10 的极小 p 值来自 rankNet(do.stat=TRUE) 的【配对 Wilcoxon
#      信息流检验】,与 permutation 检验是两回事,稿件必须分开表述。
#   4. CellChat 两个检验都【不做 FDR】(源码无 p.adjust),rankNet 的 p 值
#      由我们自己做 BH 校正后报告。
#
# 前置:已有 03_cellchat_final.R 跑出的:
#   cellchat_control_final.rds / cellchat_HFpEF_final.rds / cellchat_merged_final.rds
###############################################################################

suppressPackageStartupMessages({
  library(zellkonverter)
  library(SingleCellExperiment)
  library(Matrix)
  library(CellChat)
})

options(stringsAsFactors = FALSE)
future::plan("sequential")
set.seed(1)

out_dir <- "/Users/franxy/Documents/博士课题 6.15/严谨/results"
SEED  <- 1
NBOOT <- 100

## ===========================================================================
## A. 可复现性存档:原样重跑 control 组的 computeCommunProb,断言 pval 逐位相同
## ===========================================================================
message("========== A. 可复现性存档 ==========")

# 读回出版版 control 对象作为"基准"
cc_ctrl_saved <- readRDS(file.path(out_dir, "cellchat_control_final.rds"))
pval_saved <- cc_ctrl_saved@net$pval

# 从同一个对象重跑 computeCommunProb(同 seed、同 nboot、单线程)
# 注:重跑前对象需保留 data.signaling 与 over-expressed 结果;_final.rds 里都在
message("== 原样重跑 control 组 computeCommunProb(seed.use=1, nboot=100)…")
set.seed(SEED)
cc_ctrl_rerun <- computeCommunProb(cc_ctrl_saved,
                                   type = "truncatedMean", trim = 0.1,
                                   population.size = TRUE,
                                   seed.use = SEED, nboot = NBOOT)
pval_rerun <- cc_ctrl_rerun@net$pval

# 断言:两次 pval 是否逐位相同
is_identical <- identical(pval_saved, pval_rerun)
max_abs_diff <- max(abs(pval_saved - pval_rerun))
cat(sprintf("\n>> identical(saved pval, rerun pval) = %s\n", is_identical))
cat(sprintf(">> 两次 pval 最大绝对差 = %.3e(应为 0)\n", max_abs_diff))

# 记录可复现性环境信息(写进 reproducibility statement 用)
cat("\n== RNGkind(应为 Mersenne-Twister / Inversion / Rejection)==\n")
print(RNGkind())
repro_log <- file.path(out_dir, "reproducibility_log.txt")
sink(repro_log)
cat("CellChat reproducibility check\n")
cat("Date:", as.character(Sys.time()), "\n")
cat("identical(saved, rerun) pval:", is_identical, "\n")
cat("max abs diff:", max_abs_diff, "\n")
cat("RNGkind:", paste(RNGkind(), collapse=" / "), "\n\n")
print(sessionInfo())
sink()
cat(sprintf("\n可复现性日志已存:%s\n", repro_log))

## ===========================================================================
## B. 两类 p 值彻底分开:permutation(net$pval)vs rankNet Wilcoxon
## ===========================================================================
message("\n========== B. 两类 p 值正确归因 ==========")

cellchat <- readRDS(file.path(out_dir, "cellchat_merged_final.rds"))

# --- B1. permutation p 值的真实分辨率(证明 1e-8 不可能来自这里)---
cat("\n== B1. computeCommunProb permutation p 值的可取值(nboot=100)==\n")
pv <- cc_ctrl_saved@net$pval
pv_nonzero <- pv[pv > 0]
cat("permutation p 值的唯一非零取值(前 10 个):\n")
print(sort(unique(as.vector(pv_nonzero)))[1:10])
cat(sprintf("permutation p 值最小非零值 = %.4f(理论下限 = 1/nboot = %.2f)\n",
            if (length(pv_nonzero)) min(pv_nonzero) else NA, 1/NBOOT))
cat(">> 结论:permutation 检验只能报 p<0.05(分辨率 0.01),不可能是 1e-8。\n")

# --- B2. rankNet 的 Wilcoxon p 值(那些极小 p 值的真正来源)+ BH 校正 ---
cat("\n== B2. rankNet(do.stat=TRUE) Wilcoxon 信息流检验 + BH 校正 ==\n")
gg <- rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE,
              return.data = TRUE)
# rankNet 返回的数据里含每条通路的 Wilcoxon p 值
rk <- gg$signaling.contribution
# 每条通路一个 p 值(两组各一行,p 值相同),去重到通路层面做 BH
pw <- unique(rk[, c("name", "pvalues")])
pw <- pw[order(pw$pvalues), ]
pw$pvalues_BH <- p.adjust(pw$pvalues, method = "BH")
pw$signif_raw <- pw$pvalues < 0.05
pw$signif_BH  <- pw$pvalues_BH < 0.05

cat("\n各通路 Wilcoxon p 值 + BH 校正(按原始 p 排序):\n")
print(pw, row.names = FALSE)

write.csv(pw, file.path(out_dir, "rankNet_wilcox_BH.csv"), row.names = FALSE)
cat(sprintf("\n>> 通路级 Wilcoxon + BH 结果已存:%s\n", file.path(out_dir, "rankNet_wilcox_BH.csv")))

# --- B3. 处理"单组独有通路 p=0"的伪结果(专家提醒)---
cat("\n== B3. 检查 p=0 的通路(可能是单组独有/向量塌缩的构造性 0,非真实检验)==\n")
zero_p <- pw[pw$pvalues == 0, "name"]
if (length(zero_p)) {
  cat("以下通路 p=0,需人工核查是否单组独有(若是,不能当真实显著):\n")
  print(zero_p)
} else {
  cat("无 p=0 通路。\n")
}

## ===========================================================================
## C. 小结
## ===========================================================================
message("\n========== C. 第 1 步完成 ==========")
cat("产物:\n")
cat("  reproducibility_log.txt   : 可复现性断言 + sessionInfo + RNGkind\n")
cat("  rankNet_wilcox_BH.csv     : 通路级 Wilcoxon p 值 + BH 校正后的值\n")
cat("\n稿件表述要点:\n")
cat("  - permutation 检验(net$pval)报 p<0.05(分辨率 0.01),不报具体极小值\n")
cat("  - 1e-8/1e-10 归因于 rankNet 配对 Wilcoxon 信息流检验\n")
cat("  - rankNet p 值已 BH 校正,报告 signif_BH 列\n")
cat("  - CellChat 结果定位为 hypothesis-generating,需 DE/正交方法佐证\n")
