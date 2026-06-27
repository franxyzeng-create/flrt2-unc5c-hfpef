###############################################################################
# Fig 1 — Per-donor pathway shift heatmap (Step 3, within-snRNA / 表达层↔交互层)
#
# 讲什么:guidance-cue 通路在 per-donor 两个 co-primary(Absolute population.size=TRUE
#         / Relative)下的 HFpEF vs control 方向 + BH 显著性。
#         展示 UNC5/FLRT/ADGRL 在 snRNA 内"双轴全清"(within-modality 硬核心),
#         为后续 Fig3 跨模态正交验证(FLRT2 corroborate / ADGRL2 disconfirm)铺垫。
#
# ★ 定性边界(写 caption 时遵守):本图是 within-modality(snRNA, CellChat 推断)证据;
#   ADGRL 在此通过,不代表跨模态成立——ADGRL 出局的判据在 Fig3(正交 bulk)。
#   "四轴"命名若与你定义不同请指正:此处按 {UNC5,FLRT,ADGRL,Netrin}(+SLIT)理解。
#
# 数据源:perdonor_test_absolute.csv / perdonor_test_relative.csv
#   关键列:pathway, direction(HFpEF_up/down), HL_shift(Hodges-Lehmann), p_wilcox_BH
###############################################################################

source("/Users/franxy/Documents/博士课题 6.15/FLRT2&UNC5C/fig_common.R")
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

## ---- 配置:哪些通路、顺序、角色标注 --------------------------------------
GUIDANCE <- c("UNC5", "FLRT", "ADGRL", "Netrin", "SLIT")   # 主展示(导向线索)
CONTEXT  <- c("COLLAGEN", "LAMININ", "FN1")                 # 纤维化背景(可选,见 SHOW_CONTEXT)
SHOW_CONTEXT <- TRUE

pw_role <- function(pw) {
  if (pw %in% c("UNC5", "FLRT"))       return("locked core (FLRT2->UNC5C)")
  if (pw == "ADGRL")                    return("out of core (ADGRL)")
  if (pw %in% CONTEXT)                  return("fibrosis context")
  "other guidance cue"
}

## ---- 读取 + 整形 ---------------------------------------------------------
load_axis <- function(fname, axis_label) {
  d <- read_result_csv(fname)
  keep <- c(GUIDANCE, if (SHOW_CONTEXT) CONTEXT else character(0))
  miss <- setdiff(keep, d$pathway)
  if (length(miss)) message("[!] ", axis_label, " 缺通路: ", paste(miss, collapse=", "))
  d <- d[d$pathway %in% keep, c("pathway", "direction", "HL_shift", "p_wilcox_BH")]
  d$axis <- axis_label
  d
}

abs_d <- load_axis("perdonor_test_absolute.csv", "Absolute\n(population.size=TRUE)")
rel_d <- load_axis("perdonor_test_relative.csv", "Relative")
dat   <- bind_rows(abs_d, rel_d)

# 编码:sign = +1(HFpEF_up)/ -1(down);fill = sign * -log10(BH)(方向化显著性)
dat <- dat %>%
  mutate(
    sign      = if_else(direction == "HFpEF_up", 1, -1),
    neglog    = -log10(pmax(p_wilcox_BH, 1e-12)),
    signed    = sign * neglog,
    stars     = p_stars(p_wilcox_BH),
    role      = vapply(pathway, pw_role, character(1))
  )

# 行顺序:核心在上 → ADGRL → 其他导向 → 背景
row_order <- c("UNC5", "FLRT", "ADGRL", "Netrin", "SLIT",
               if (SHOW_CONTEXT) CONTEXT else character(0))
row_order <- intersect(row_order, unique(dat$pathway))
dat$pathway <- factor(dat$pathway, levels = rev(row_order))   # rev: 让第一个在顶部
dat$axis    <- factor(dat$axis, levels = c("Absolute\n(population.size=TRUE)", "Relative"))

# 色阶上限(对称),稳健避免极端值吃掉色带
lim <- max(4, ceiling(max(abs(dat$signed), na.rm = TRUE)))

## ---- 画图 -----------------------------------------------------------------
p <- ggplot(dat, aes(axis, pathway, fill = signed)) +
  # ADGRL 行整体淡化(alpha=0.4):within-modality 通过但跨模态出局,降视觉权重(专家 §7q1)
  geom_tile(aes(alpha = pathway == "ADGRL"), color = "white", linewidth = 1.1) +
  # 再给 ADGRL 行加 longdash 灰框,明确"已被 Fig3 裁定"
  geom_tile(data = subset(dat, pathway == "ADGRL"),
            fill = NA, color = "grey25", linewidth = 0.9, linetype = "longdash") +
  geom_text(aes(label = sprintf("%s\nHL=%.3g", stars, HL_shift)),
            size = 3.1, lineheight = 0.9, color = "black") +
  scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.4), guide = "none") +
  scale_fill_gradient2(
    name = expression(bold(paste("signed  ", -log[10], "(BH)"))),
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-lim, lim)) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0),
                   labels = function(x) ifelse(x == "ADGRL", "ADGRL (out)", x)) +
  labs(
    title    = "Figure 1. Per-donor pathway enrichment in HFpEF (Step 3)",
    subtitle = "snRNA / CellChat; two co-primary axes. Cell = BH-Wilcoxon stars + Hodges-Lehmann shift. Red = up in HFpEF.",
    x = NULL, y = NULL,
    caption  = paste0(
      "Within-modality (snRNA) evidence only; color = significance, HL text = effect size (they can diverge).\n",
      "UNC5/FLRT/ADGRL pass both axes here, but ADGRL (faded, dashed, 'out') is adjudicated OUT of the locked core by orthogonal bulk (Fig 3).\n",
      "Stars: ****<1e-4 ***<1e-3 **<1e-2 *<0.05. n = 18 HFpEF vs 24 control donors.")) +
  theme_pub(11) +
  theme(axis.text.y = element_text(face = "bold"),
        axis.ticks  = element_blank(),
        panel.border = element_blank())

# 在核心两行左侧加一个角色色条(用第二个 geom 标注),简洁起见用 y 轴标签颜色
role_map <- dat %>% distinct(pathway, role)
ycols <- ROLE_COLORS[role_map$role[match(levels(dat$pathway), role_map$pathway)]]
p <- p + theme(axis.text.y = element_text(color = ycols, face = "bold"))

save_fig(p, "fig1_fouraxis_heatmap", w = 8, h = 5.2)

## ---- 控制台自检(便于贴回核对) ------------------------------------------
cat("\n[Fig1] 通路 × 轴 × 方向 × BH:\n")
print(dat %>% select(pathway, axis, direction, HL_shift, p_wilcox_BH, stars) %>%
        arrange(pathway, axis), row.names = FALSE)
cat("\n说明:UNC5/FLRT 双轴 ****/***;ADGRL 绝对轴**** 相对轴**(within-modality 通过);",
    "\nNetrin 相对轴 ns(rel BH≈0.10);SLIT 双轴显著但去污染敏感(次级)。\n")
