###############################################################################
# Fig 4 — Per-donor FLRT2/UNC5C co-expression in fibroblasts (表达层, measured)
#
# 讲什么:per-donor 成纤维细胞中 FLRT2+UNC5C 双阳性比例,HFpEF vs control。
#   - 双阳(headline,CellBender de-ambient 层):HFpEF ~41% vs control ~17%(~2.5×)
#   - 拆解:FLRT2+ 翻倍(驱动);UNC5C+ 组成性 ~85% 近饱和(受体预就位)
#   - 措辞:"受体就位 + 配体诱导使自分泌成纤维群扩大 ~2.5×",非协同、非第二颗子弹
#
# ★ 专家原则:per-donor 聚合(非 cell-level pooled,避免 pseudoreplication);
#   报 donor n + effect size(Cliff's delta)+ p。p/delta 脚本内现算,不硬编码。
#
# 数据源(headline):coexpr_per_donor_cellbender.csv(CB 层)
#         (敏感性):coexpr_per_donor.csv(CR 层)— 改 LAYER 切换
###############################################################################

source("/Users/franxy/Documents/博士课题 6.15/FLRT2&UNC5C/fig_common.R")
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(patchwork) })

LAYER <- "CB"   # "CB"(headline,de-ambient) 或 "CR"(敏感性)
fname <- if (LAYER == "CB") "coexpr_per_donor_cellbender.csv" else "coexpr_per_donor.csv"

d <- read_result_csv(fname)
d$grp <- factor(d$grp, levels = GRP_LEVELS)

## ---- 现算统计(双阳 + 两个单阳) -----------------------------------------
stat_co  <- two_group_summary(d$pct_coexpr, d$grp)
stat_fl  <- two_group_summary(d$pct_FLRT2,  d$grp)
stat_un  <- two_group_summary(d$pct_UNC5C,  d$grp)
cat(sprintf("[Fig4/%s] n: HFpEF=%d, control=%d\n", LAYER, stat_co$n_HFpEF, stat_co$n_control))
cat(sprintf("  双阳 coexpr : HFpEF %.1f%% vs control %.1f%%  | p=%.2g, Cliff's d=%.3f, %.2f×\n",
            stat_co$median_HFpEF, stat_co$median_control, stat_co$p_wilcox, stat_co$cliffs_d, stat_co$fold))
cat(sprintf("  FLRT2+      : HFpEF %.1f%% vs control %.1f%%  | p=%.2g, d=%.3f, %.2f×\n",
            stat_fl$median_HFpEF, stat_fl$median_control, stat_fl$p_wilcox, stat_fl$cliffs_d, stat_fl$fold))
cat(sprintf("  UNC5C+      : HFpEF %.1f%% vs control %.1f%%  | p=%.2g, d=%.3f, %.2f×\n",
            stat_un$median_HFpEF, stat_un$median_control, stat_un$p_wilcox, stat_un$cliffs_d, stat_un$fold))
cat("  ★ 注:UNC5C+ δ=0.86(与双阳相同),效应量并不小 → 措辞用 'near-ceiling/high baseline + modest 1.1× fold',\n")
cat("     不用 'constitutive/不变'(否则违反 meta④ 以 effect size 论)。'FLRT2 驱动' 由 fold 论证(2.2× vs 1.1×)。\n")

## ---- ★ 独立性检验(专家要求做实:per-donor obs/expected,不只 median 比 median) -----
# 每 donor:观测双阳% vs 边际乘积期望%(独立假设下)。检验 log2(obs/expected) 是否偏离 0。
di <- d %>%
  mutate(expected  = pmax(pct_FLRT2/100 * pct_UNC5C/100 * 100, 1e-9),
         log2ratio = log2(pmax(pct_coexpr, 1e-9) / expected))
wt_ind  <- suppressWarnings(wilcox.test(di$log2ratio, mu = 0, conf.int = TRUE))
ind_p   <- wt_ind$p.value
ind_med <- 2^median(di$log2ratio)                         # 观测/期望 中位倍数(=1 即完全独立)
ind_ciL <- 2^wt_ind$conf.int[1]; ind_ciR <- 2^wt_ind$conf.int[2]
# per-group 倍数(检验同口径):证 ~2% 偏离【非 HFpEF 特异】(round-2 专家要求)
ind_grp <- di %>% group_by(grp) %>% summarise(ratio = 2^median(log2ratio), .groups = "drop")
ind_hf  <- ind_grp$ratio[ind_grp$grp == "HFpEF"]
ind_ct  <- ind_grp$ratio[ind_grp$grp == "control"]
cat(sprintf("\n[独立性] per-donor obs/expected 倍数 median=%.2f [%.2f,%.2f]; ", ind_med, ind_ciL, ind_ciR))
cat(sprintf("one-sample Wilcoxon(log2ratio vs 0) p=%.2g, n=%d donors\n", ind_p, nrow(di)))
cat(sprintf("  解读★:obs/expected≈1.0(~%.0f%% 高于独立期望),n=42 下形式显著(p≈%.3f)但幅度可忽略 → near-independent。\n",
            (ind_med - 1) * 100, ind_p))
cat(sprintf("       且偏离【非 HFpEF 特异】:HFpEF ratio=%.3f vs control=%.3f(对照不低于病例)→ 强化结论。\n", ind_hf, ind_ct))
cat("       措辞去机制词 'synergy',只作统计陈述。\n")
print(as.data.frame(di %>% group_by(grp) %>%
        summarise(obs = median(pct_coexpr), expected = median(expected),
                  ratio = 2^median(log2ratio), .groups = "drop")),
      row.names = FALSE, digits = 3)

## ---- 通用 boxplot 函数 ---------------------------------------------------
box_one <- function(df, yvar, ylab, st, title) {
  ymax <- max(df[[yvar]], na.rm = TRUE)
  ggplot(df, aes(grp, .data[[yvar]], fill = grp)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.45, color = "grey30") +
    geom_jitter(aes(color = grp), width = 0.12, height = 0, size = 1.8, alpha = 0.9) +
    annotate("segment", x = 1, xend = 2, y = ymax*1.08, yend = ymax*1.08, color = "black") +
    annotate("text", x = 1.5, y = ymax*1.14,
             label = sprintf("p=%.2g | Cliff's d=%.2f | %.1f×",
                             st$p_wilcox, st$cliffs_d, st$fold), size = 3.1) +
    scale_fill_manual(values = GRP_COLORS, guide = "none") +
    scale_color_manual(values = GRP_COLORS, guide = "none") +
    scale_y_continuous(limits = c(0, ymax*1.2)) +
    labs(title = title, x = NULL, y = ylab) +
    theme_pub(11)
}

pA <- box_one(d, "pct_coexpr", "FLRT2+UNC5C+ fibroblasts (% per donor)", stat_co,
              "Double-positive (autocrine-competent)")

# 单阳:长表,facet
dl <- d %>% select(donor, grp, pct_FLRT2, pct_UNC5C) %>%
  pivot_longer(c(pct_FLRT2, pct_UNC5C), names_to = "marker", values_to = "pct") %>%
  mutate(marker = recode(marker, pct_FLRT2 = "FLRT2+ (ligand)", pct_UNC5C = "UNC5C+ (receptor)"))

ann <- bind_rows(
  transform(stat_fl, marker = "FLRT2+ (ligand)"),
  transform(stat_un, marker = "UNC5C+ (receptor)"))

pB <- ggplot(dl, aes(grp, pct, fill = grp)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.45, color = "grey30") +
  geom_jitter(aes(color = grp), width = 0.12, height = 0, size = 1.5, alpha = 0.85) +
  geom_text(data = ann, inherit.aes = FALSE,
            aes(x = 1.5, y = 102, label = sprintf("p=%.2g | d=%.2f | %.1f×", p_wilcox, cliffs_d, fold)),
            size = 2.9) +
  facet_wrap(~ marker) +
  scale_fill_manual(values = GRP_COLORS, guide = "none") +
  scale_color_manual(values = GRP_COLORS, guide = "none") +
  scale_y_continuous(limits = c(0, 110)) +
  labs(title = "Single-positive decomposition", x = NULL,
       y = "% of fibroblasts per donor") +
  theme_pub(11)

p <- (pA | pB) +
  plot_annotation(
    tag_levels = "A",
    title = sprintf("Figure 4. Per-donor co-expression of FLRT2 & UNC5C in fibroblasts (%s layer)", LAYER),
    subtitle = "Receptor UNC5C near-ceiling / high baseline in both groups (higher in HFpEF but modest ~1.1x fold, ceiling); ligand FLRT2 (~2.2x) is the driver -> autocrine-competent pool expands ~2.5x.",
    caption = sprintf(paste0(
      "Per-donor proportions (not cell-level pooled); exact Wilcoxon rank-sum, no ties. Each point = one donor. CellBender layer excludes ambient.\n",
      "Near-independent co-detection: per-donor obs/expected = %.2f [%.2f, %.2f] (~%.0f%% above independence; formally significant at n=%d, p=%.2g, but negligible) and NOT HFpEF-specific (HFpEF %.2f vs control %.2f)."),
      ind_med, ind_ciL, ind_ciR, (ind_med - 1) * 100, nrow(di), ind_p, ind_hf, ind_ct),
    theme = theme(plot.title = element_text(face = "bold", family = SUB_FONT),
                  plot.subtitle = element_text(color = "grey30", size = 9, family = SUB_FONT),
                  plot.tag = element_text(face = "bold", size = 14, family = SUB_FONT))) +
  plot_layout(widths = c(1, 1.5))

save_fig(p, sprintf("fig4_coexpr_perdonor_%s", LAYER), w = 11, h = 5.2)
# 独立性检验已上移到统计段(per-donor obs/expected + 单样本 Wilcoxon),结果进 caption。
