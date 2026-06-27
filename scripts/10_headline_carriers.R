###############################################################################
# Step 3 前置:8 条 robust 头条通路的 sender→receiver 承载分析
#
# 专家指出:定"固定公共类型集"前,必须先确定"头条必需细胞类型"的【精确底线】——
#   即 8 条三臂 robust 通路各自主要由哪些 sender→receiver 细胞对承载,不能凭印象。
#   这份清单是固定集的硬底线(集必须包含这些类型)。
#
# 8 条 robust 通路(三臂都 BH<0.05 且 HFpEF↑):
#   EPHA, ADGRL, UNC5, Netrin, NRXN, FLRT, LAMININ, COLLAGEN
#
# 用 HFpEF 组对象的 netP$prob[,,pathway] 切片(关心信号在 HFpEF 里由谁承载)。
# 不重跑 CellChat,直接读已存对象。
###############################################################################

suppressPackageStartupMessages({ library(CellChat); library(Matrix) })

out_dir <- "/Users/franxy/Documents/博士课题 6.15/严谨/results"

# 8 条 robust 头条通路
ROBUST <- c("EPHA","ADGRL","UNC5","Netrin","NRXN","FLRT","LAMININ","COLLAGEN")

## ---- 读已存的 HFpEF 组对象 -----------------------------------------------
# 优先用 Step2 v3 重建的 merged_raw(和三臂判定同源);取其中 HFpEF 那个 cc 对象
# merged 对象的单组对象在 @ 内不易直接取,改读出版版单组对象 cellchat_HFpEF_final.rds
# (参数与三臂 raw 臂一致:truncatedMean/trim0.1/population.size TRUE/raw.use TRUE)
cc_path <- file.path(out_dir, "cellchat_HFpEF_final.rds")
if (!file.exists(cc_path)) stop("找不到 cellchat_HFpEF_final.rds,请确认 03 脚本已存该对象")
cc <- readRDS(cc_path)

cat("== HFpEF 对象细胞类型:\n"); print(levels(cc@idents))

prob <- cc@netP$prob          # [sender × receiver × pathway]
stopifnot(length(dim(prob)) == 3)
avail_pathways <- dimnames(prob)[[3]]
cat(sprintf("\n对象中通路数: %d\n", length(avail_pathways)))

## ---- 逐条 robust 通路:拆解 sender→receiver 承载 -------------------------
cat("\n========== 8 条 robust 通路的 sender→receiver 承载 ==========\n")
celltypes <- dimnames(prob)[[1]]
all_involved <- character(0)

for (pw in ROBUST) {
  if (!pw %in% avail_pathways) {
    cat(sprintf("\n--- %s:在 HFpEF 对象中未检出(跳过)---\n", pw)); next
  }
  M <- prob[, , pw]                      # sender × receiver 概率矩阵
  total <- sum(M)
  if (total == 0) { cat(sprintf("\n--- %s:总强度 0 ---\n", pw)); next }

  # 找出贡献最大的 sender→receiver 对(累计到 ~80%)
  idx <- which(M > 0, arr.ind = TRUE)
  pairs <- data.frame(
    sender   = celltypes[idx[,1]],
    receiver = celltypes[idx[,2]],
    prob     = M[idx],
    pct      = 100 * M[idx] / total
  )
  pairs <- pairs[order(-pairs$prob), ]
  pairs$cum_pct <- cumsum(pairs$pct)

  cat(sprintf("\n--- %s(总强度 %.4f)主要细胞对(累计到80%%)---\n", pw, total))
  top <- pairs[pairs$cum_pct <= 80 | seq_len(nrow(pairs)) == 1, ]
  print(top, row.names = FALSE, digits = 3)

  # 记录该通路涉及的主要细胞类型(贡献前80%的对里出现的)
  involved <- unique(c(top$sender, top$receiver))
  all_involved <- union(all_involved, involved)
  cat(sprintf("   → 主要涉及细胞类型: %s\n", paste(involved, collapse=", ")))
}

## ---- 汇总:头条必需细胞类型(固定集的硬底线)----------------------------
cat("\n========== 头条必需细胞类型(固定公共类型集的硬底线)==========\n")
cat("以下细胞类型承载了 8 条 robust 头条通路的主要信号(累计前80%),\n")
cat("固定公共类型集【必须】包含它们:\n")
print(sort(all_involved))

# 每个细胞类型作为 sender / receiver 在多少条 robust 通路里是主要承载者
cat("\n--- 各细胞类型在多少条 robust 通路中是主要承载者 ---\n")
involve_count <- setNames(integer(length(celltypes)), celltypes)
for (pw in ROBUST) {
  if (!pw %in% avail_pathways) next
  M <- prob[, , pw]; total <- sum(M)
  if (total == 0) next
  idx <- which(M > 0, arr.ind = TRUE)
  pairs <- data.frame(s=celltypes[idx[,1]], r=celltypes[idx[,2]], p=M[idx])
  pairs <- pairs[order(-pairs$p), ]; pairs$cum <- cumsum(pairs$p)/total
  top <- pairs[pairs$cum <= 0.8 | seq_len(nrow(pairs))==1, ]
  for (ct in unique(c(top$s, top$r))) involve_count[ct] <- involve_count[ct] + 1
}
print(sort(involve_count, decreasing = TRUE))

cat("\n========== 解读指引 ==========\n")
cat("- 上面'头条必需细胞类型'= 固定公共类型集的下限(必含)。\n")
cat("- 结合 09 探查的 Q1(各类型达标病人数)/Q2(保留病人数),据专家 n 优先+底线 gate 规则定集:\n")
cat("  按病人覆盖度从高到低加类型,加到'下个非必需类型丢>10-15%病人'就停;任一组 n<~12 警惕。\n")
cat("- 情形判断:必需类型若都是高丰度(CM/Fib/Endo1)→K小n高两全;\n")
cat("  若含稀有类型(如 Endothelial2)→纳入会砍 n,需权衡或退回 Step2 细胞对证据。\n")

cat("\n========== sessionInfo ==========\n")
print(sessionInfo())
