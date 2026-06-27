###############################################################################
# fig_common.R  —  主图共享配置 / 主题 / 配色 / 辅助函数
#
# 设计原则(写图前定,遵守专家 meta-原则):
#   - 图内文字一律【英文】(投稿用 + 规避 Mac 上 ggplot 中文字体坑);
#     中文解读放 README_figures.md。
#   - 表达层(measured) vs 交互层(CellChat 推断)在标题/副标题里分开标注;
#     交互层图(Fig2)副标题写明 "CellChat-inferred, contact not verified"。
#   - 不 overclaim:不写 "synergy / second bullet / proof";autocrine 写 "single-compartment"。
#   - 数字尽量在脚本内从真实 CSV 现算(Wilcoxon / Cliff's delta),不硬编码记忆值。
#
# R 易错点(本项目踩过):标量条件用 if/else 不用 ifelse();%||% 先定义后用。
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
})

## ---- 路径(用户本地绝对路径;改这里即可整体迁移) -------------------------
ROOT        <- "/Users/franxy/Documents/博士课题 6.15"
RESULTS_DIR <- file.path(ROOT, "严谨", "results")                 # 读:真实 CSV / rds
DATA_DIR    <- file.path(ROOT, "严谨", "data")                    # 读:h5ad / bulk(若需)
FIG_SRC_DIR <- file.path(ROOT, "FLRT2&UNC5C")                     # 本脚本所在
FIG_DIR     <- file.path(ROOT, "FLRT2&UNC5C", "figures")          # 写:图输出
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

## ---- 配色(Okabe-Ito,色盲友好) ------------------------------------------
GRP_COLORS <- c(control = "#0072B2", HFpEF = "#D55E00")           # 蓝=control, 朱红=HFpEF
GRP_LEVELS <- c("control", "HFpEF")
# 通路角色配色:锁定核心 / 出局 / 其他
ROLE_COLORS <- c(
  "locked core (FLRT2->UNC5C)" = "#D55E00",
  "out of core (ADGRL)"        = "#999999",
  "other guidance cue"         = "#56B4E9",
  "fibrosis context"           = "#009E73"
)
# 9 类细胞顺序(与 script 03 KEEP 一致;绘图固定顺序用)
CELLTYPE_LEVELS <- c("Cardiomyocyte","Fibroblast","Endothelial1","Endothelial2",
                     "Pericyte","VSMC","Macrophage","Lymphocyte","Endocardial")

## ---- 出版主题(投稿统一规格) --------------------------------------------
SUB_FONT <- "Arial"   # 投稿统一字体;系统无 Arial 时 ggplot 自动回退 sans(仅告警不报错)
theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size, base_family = SUB_FONT) +
    theme(
      text             = element_text(family = SUB_FONT, color = "black"),
      axis.text        = element_text(color = "black", size = rel(0.85)),
      axis.title       = element_text(face = "bold"),
      plot.title       = element_text(face = "bold", hjust = 0, size = rel(1.1)),
      plot.subtitle    = element_text(color = "grey30", size = rel(0.8)),
      plot.caption     = element_text(color = "grey40", size = rel(0.7), hjust = 0),
      plot.tag         = element_text(face = "bold", size = rel(1.5)),   # panel A/B 标签
      legend.title     = element_text(face = "bold", size = rel(0.85)),
      legend.text      = element_text(size = rel(0.8)),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold")
    )
}

## ---- 辅助函数 -------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a   # 标量条件 if/else

# 读 CSV(保留原始列名;文件不存在则明确报错)
read_result_csv <- function(fname) {
  fp <- file.path(RESULTS_DIR, fname)
  if (!file.exists(fp)) stop("找不到文件: ", fp, "\n请确认 RESULTS_DIR 是否正确。")
  read.csv(fp, check.names = FALSE, stringsAsFactors = FALSE)
}

# 显著性星号(向量化;基于校正后 p)
p_stars <- function(p) {
  as.character(cut(p,
    breaks = c(-Inf, 1e-4, 1e-3, 1e-2, 5e-2, Inf),
    labels = c("****", "***", "**", "*", "ns")))
}

# Cliff's delta(x vs y;正值=x 倾向更大)。HFpEF 传 x、control 传 y → 正值=HFpEF 更高
cliffs_delta <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  (sum(outer(x, y, ">")) - sum(outer(x, y, "<"))) / (length(x) * length(y))
}

# 两组比较的一行汇总:median、Wilcoxon p、Cliff's delta、倍数
# ★ p:用 wilcox.test 默认(n<50 且无 ties 时走 exact,复现 locked 的精确 p≈2e-7;
#    有 ties 则自动回退渐近近似并 suppressWarnings)。勿写死 exact=FALSE(那会得 ~2.7e-6)。
two_group_summary <- function(value, group, hf = "HFpEF", ct = "control") {
  x <- value[group == hf]; y <- value[group == ct]
  wt <- suppressWarnings(wilcox.test(x, y))
  data.frame(
    n_HFpEF   = length(x), n_control = length(y),
    median_HFpEF = median(x, na.rm = TRUE), median_control = median(y, na.rm = TRUE),
    p_wilcox  = wt$p.value,
    cliffs_d  = cliffs_delta(x, y),
    fold      = median(x, na.rm = TRUE) / median(y, na.rm = TRUE)
  )
}

# 保存图(投稿规格)。本机无 XQuartz → cairo_pdf 不可用;改用 svglite 出 Arial 忠实矢量(SVG,可转 PDF/EPS)+ PNG 600dpi。
# 若已装 XQuartz,可启用下面 cairo_pdf 行直出 Arial PDF。
save_fig <- function(plot, name, w = 7, h = 5) {
  ggsave(file.path(FIG_DIR, paste0(name, ".svg")), plot, width = w, height = h, device = svglite::svglite)
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), plot, width = w, height = h, dpi = 600)
  # ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), plot, width = w, height = h, device = cairo_pdf)  # 需 XQuartz
  message("  saved: ", file.path(FIG_DIR, paste0(name, ".{svg,png}")), "  (PNG 600dpi / SVG Arial 矢量)")
}

message("fig_common.R loaded. FIG_DIR = ", FIG_DIR)
