###############################################################################
# Fig 5 — Technical / compositional balance check (Step 3 confound check)
#
# 讲什么:per-donor 技术与构成变量在 HFpEF vs control 是否平衡,
#         以排除"通讯增强是测序深度/核数/细胞构成假象"。
#   - 平衡检查项:total_nuclei, median_depth, fib_frac, cm_frac(应 ns → 平衡)
#   - total_flow = 每 donor 总通讯流量(outcome,非混杂;预期 HFpEF↑,单列标注)
#
# ★ 诚实:若某平衡项显著(如 fib_frac),不掩盖——如实显示 p,留作 limitation。
#         图只呈现数据,结论交给 caption + 专家。
#
# 数据源:step3_confound_check.csv
###############################################################################

source("/Users/franxy/Documents/博士课题 6.15/FLRT2&UNC5C/fig_common.R")
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

d <- read_result_csv("step3_confound_check.csv")
d$group <- factor(d$group, levels = GRP_LEVELS)

## ---- 指标定义(平衡检查 vs outcome) -------------------------------------
metrics <- tibble::tribble(
  ~var,             ~label,                          ~kind,
  "total_nuclei",   "Total nuclei / donor",          "balance",
  "median_depth",   "Median depth (UMI) / donor",    "balance",
  "detected_genes", "Detected genes / donor",        "balance",
  "fib_frac",       "Fibroblast fraction",           "balance",
  "cm_frac",        "Cardiomyocyte fraction",        "balance",
  "total_flow",     "Total communication (outcome)", "outcome"
)
metrics <- metrics %>% filter(var %in% names(d))

## ---- 现算每指标 Wilcoxon p + 方向 ----------------------------------------
stats <- lapply(metrics$var, function(v) {
  s <- two_group_summary(d[[v]], d$group)
  data.frame(var = v, p_wilcox = s$p_wilcox, cliffs_d = s$cliffs_d,
             median_HFpEF = s$median_HFpEF, median_control = s$median_control)
}) %>% bind_rows() %>% left_join(metrics, by = "var")
cat("[Fig5] per-metric HFpEF vs control:\n")
print(stats %>% select(label, kind, median_control, median_HFpEF, cliffs_d, p_wilcox),
      row.names = FALSE, digits = 3)

## ---- 长表 + facet 标签(带 p) -------------------------------------------
lab_map <- setNames(sprintf("%s\n(p=%.2g%s)",
                            stats$label, stats$p_wilcox,
                            ifelse(stats$kind == "outcome", ", outcome", "")),
                    stats$var)
dl <- d %>% select(group, all_of(metrics$var)) %>%
  pivot_longer(-group, names_to = "var", values_to = "value") %>%
  mutate(facet = factor(lab_map[var], levels = lab_map[metrics$var]),
         is_outcome = var == "total_flow")

p <- ggplot(dl, aes(group, value, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.45, color = "grey30") +
  geom_jitter(aes(color = group), width = 0.12, height = 0, size = 1.4, alpha = 0.85) +
  facet_wrap(~ facet, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = GRP_COLORS, guide = "none") +
  scale_color_manual(values = GRP_COLORS, guide = "none") +
  labs(title = "Figure 5. Per-donor technical & compositional balance (Step 3 confound check)",
       subtitle = "First four = balance checks (expect ns). Last = total communication (outcome, expected to differ).",
       x = NULL, y = NULL,
       caption = paste0(
         "Each point = one donor. Median depth trends LOWER in HFpEF (p=0.069, Cliff's d=-0.33); ",
         "since detection / communication are HIGHER despite lower depth, this bias is conservative (against the finding).\n",
         "Any significant balance metric is shown honestly and carried as a limitation.")) +
  theme_pub(11) +
  theme(strip.text = element_text(size = 8.5))

save_fig(p, "fig5_confounder_boxplot", w = 14, h = 4)
