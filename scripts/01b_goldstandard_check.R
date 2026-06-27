###############################################################################
# 第 0b 步:金标准对照 —— 本地 pseudobulk logFC vs 原文 Table S8 logFC
#
# 目的:把你第 0 步算出的成纤维/周细胞 logFC,与原文补充表 ST8 中
#       同细胞类型的 logFC 按基因对齐,计算 Pearson 相关并画散点。
#       若 r > 0.95,则你的整条 pseudobulk 管线与原文实锤一致。
#
# 前置:先跑过 01_pseudobulk_positive_control.R,已生成
#       results/pseudobulk_positive_control.rds(或两个 DE_*.csv)
#
# ST8 结构(已确认):长表,标题在第 4 行,数据从第 5 行起。
#   第 1 列 Gene、第 5 列 Cell type、
#   第 13 列 logFC CellBender、第 16 列 adj.P.Val CellBender、
#   第 17 列 logFC CellRanger、第 20 列 adj.P.Val CellRanger
###############################################################################

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
})

xl      <- "/Users/franxy/Documents/博士课题 6.15/严谨/data/Hahn 人 HFpEF 单核 snRNA/circres-2025-327433-s02.xlsx"
out_dir <- "/Users/franxy/Documents/博士课题 6.15/严谨/results"

## ---- 1. 读 ST8,用第 4 行做列名,数据从第 5 行起 -------------------------
st8 <- read_excel(xl, sheet = "ST8.HFpEFvsCtrl DEG", skip = 3)  # 跳过前 3 行,第 4 行成为表头

# 清理列名里的换行符,方便引用
names(st8) <- gsub("[\r\n]+", " ", names(st8))
names(st8) <- trimws(names(st8))
cat("=== ST8 实际列名 ===\n"); print(names(st8))

# 定位关键列(按名字模糊匹配,稳妥)
col_gene <- names(st8)[1]                                  # 第 1 列 = Gene
col_ct   <- grep("^Cell type", names(st8), value = TRUE)[1]
col_lfc_cb <- grep("logFC.*CellBender", names(st8), value = TRUE)[1]
col_fdr_cb <- grep("adj.P.Val.*CellBender", names(st8), value = TRUE)[1]
cat(sprintf("\n用作对齐的列:\n  基因: %s\n  细胞类型: %s\n  logFC: %s\n  FDR: %s\n",
            col_gene, col_ct, col_lfc_cb, col_fdr_cb))

## ---- 2. 打印 Cell type 的所有取值(确认筛选字符串) ----------------------
cat("\n=== ST8 中 Cell type 的唯一取值 ===\n")
print(table(st8[[col_ct]], useNA = "ifany"))

## ---- 3. 读取你第 0 步的结果 -----------------------------------------------
res <- readRDS(file.path(out_dir, "pseudobulk_positive_control.rds"))
mine <- list(
  Fibroblast = res$Fibroblast$merged,   # 含 gene, logFC_cellbender, ...
  Pericyte   = res$Pericyte$merged
)

# ST8 里成纤维/周细胞对应的 Cell type 字符串(自动匹配,避免猜错)
ct_pattern <- list(
  Fibroblast = "fibroblast",   # 匹配 "fibroblast of cardiac tissue"
  Pericyte   = "pericyte"
)

## ---- 4. 逐细胞类型对齐 + 相关 + 散点 --------------------------------------
for (nm in names(mine)) {
  # 从 ST8 筛该细胞类型
  ct_vals <- unique(st8[[col_ct]])
  hit <- ct_vals[grepl(ct_pattern[[nm]], ct_vals, ignore.case = TRUE)]
  if (length(hit) == 0) {
    cat(sprintf("\n[!] %s:在 ST8 Cell type 里没匹配到 '%s',跳过。请看上面取值表手动指定。\n",
                nm, ct_pattern[[nm]])); next
  }
  paper <- st8[st8[[col_ct]] %in% hit, c(col_gene, col_lfc_cb, col_fdr_cb)]
  names(paper) <- c("gene", "logFC_paper", "FDR_paper")
  paper$logFC_paper <- as.numeric(paper$logFC_paper)

  # 与你的结果按基因对齐
  m <- merge(mine[[nm]][, c("gene","logFC_cellbender","FDR_cellbender")],
             paper, by = "gene")
  m <- m[is.finite(m$logFC_cellbender) & is.finite(m$logFC_paper), ]

  r_all <- cor(m$logFC_cellbender, m$logFC_paper)
  # 只看两边都显著的基因(更能反映真实信号一致性)
  sig <- m$FDR_cellbender < 0.05 & m$FDR_paper < 0.05
  r_sig <- if (sum(sig) > 2) cor(m$logFC_cellbender[sig], m$logFC_paper[sig]) else NA

  cat(sprintf("\n========== %s ==========\n", nm))
  cat(sprintf("ST8 匹配到的 Cell type: %s\n", paste(hit, collapse = " | ")))
  cat(sprintf("可对齐基因数: %d(其中两边都显著: %d)\n", nrow(m), sum(sig)))
  cat(sprintf(">> logFC Pearson r(全部基因) = %.4f\n", r_all))
  cat(sprintf(">> logFC Pearson r(都显著)   = %.4f\n", r_sig))

  # 散点图
  p <- ggplot(m, aes(logFC_paper, logFC_cellbender)) +
    geom_point(alpha = 0.25, size = 0.6) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
    labs(title = sprintf("%s: 本地 vs 原文 ST8 (CellBender logFC)", nm),
         subtitle = sprintf("r(all)=%.3f, r(sig)=%.3f, n=%d", r_all, r_sig, nrow(m)),
         x = "原文 ST8 logFC (HFpEF vs Ctrl)",
         y = "本地 pseudobulk logFC") +
    theme_bw()
  ggsave(file.path(out_dir, sprintf("goldcheck_%s.png", nm)),
         p, width = 5, height = 5, dpi = 150)
}

cat(sprintf("\n== 完成。散点图已存到 %s（goldcheck_*.png）\n", out_dir))
cat("判读标准:r(all) > 0.95 即认为管线与原文一致;通常会 >0.98。\n")
