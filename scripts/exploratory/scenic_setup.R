###############################################################################
# SCENIC 环境准备(后台挂着装,与 CellChat 探查并行)
#
# 分两部分:
#   A. 安装 R 包(AUCell / RcisTarget / GENIE3 / SCENIC)
#   B. 下载人类 cisTarget motif 数据库(约 1GB,一次性)
#
# 说明:R 版 SCENIC 的 GRN 推断(GENIE3)在 1.4 万核上较吃内存。
#       本脚本先把环境备好;真正跑 SCENIC 时我们会用"按细胞类型抽样"
#       或改用 pySCENIC(更快),到时再决定本地还是上云。
###############################################################################

## ===== A. 安装 R 包 ========================================================
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Bioconductor 三件套
BiocManager::install(c("AUCell", "RcisTarget", "GENIE3"), update = FALSE, ask = FALSE)

# SCENIC 主包(来自 GitHub;若无 remotes 先装)
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("aertslab/SCENIC")

# 验证
for (p in c("AUCell","RcisTarget","GENIE3","SCENIC")) {
  cat(sprintf("%-12s : %s\n", p,
      if (requireNamespace(p, quietly=TRUE)) as.character(packageVersion(p)) else "安装失败"))
}

## ===== B. 下载 cisTarget motif 数据库 ======================================
# 人类 hg38,基于 refseq-r80,500bp-up & 100bp-down + 10kb 两套 feather 文件。
# 这是 SCENIC 打分(RcisTarget)必需的。总共约 1GB,下载一次永久可用。
#
# 注意:文件大、服务器在国外,可能较慢。建议在终端用 wget/curl 挂着下,
#       比在 R 里下更稳、可断点续传。下面给终端命令(推荐),R 内下载备选。

db_dir <- "/Users/franxy/Documents/博士课题 6.15/严谨/cisTarget_db"
dir.create(db_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n========================================================\n")
cat("请在【终端 Terminal】里跑下面的命令下载 motif 数据库(推荐,可续传):\n\n")
cat(sprintf('cd "%s"\n\n', db_dir))
cat('# 两个 feather 排序数据库(基因 × motif 排名)\n')
cat('curl -O -C - https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather\n\n')
cat('curl -O -C - https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather\n\n')
cat('# motif 到 TF 的注释表\n')
cat('curl -O -C - https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl\n')
cat("========================================================\n\n")
cat(sprintf("下载完成后,这些文件应在:%s\n", db_dir))
cat("下一步跑 SCENIC 时,我会用这个目录初始化 RcisTarget。\n")
