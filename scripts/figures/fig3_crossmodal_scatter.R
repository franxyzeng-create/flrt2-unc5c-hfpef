###############################################################################
# Fig 3 — Cross-modal orthogonal validation (表达层, measured)
#
# 讲什么:区室特异 snRNA pseudobulk log2FC(配体查 sender 区室、受体查 receiver 区室)
#         vs 正交 bulk(Hahn 2021)log2FC。
#         - FLRT2(Fib)两模态一致↑ = IRON corroborate(锁定核心的配体)
#         - ADGRL2(CM)snRNA↑ / bulk 显著↓ = DISCONFIRM(出局判据)
#         - UNC5C(Fib)snRNA↑ / bulk ns = uninformative(稀释,justified,非双标)
#
# ★ 跨模态联合判据(round-2 收敛为 3 类,同一把尺子):
#   corroborate   = snRNA sig↑ & bulk sig↑              (FLRT2, TENM3)
#   disconfirm    = snRNA sig↑ & bulk sig↓ (sign-flip)  (仅 ADGRL2)
#   uninformative = 任一端 ns / 低表达                  (UNC5C, UNC5B, TENM2, TENM4, ADGRL1, ADGRL3)
#   注:ADGRL1/3 虽 bulk 显著↓,但 snRNA ns(CI 跨 0)→ 不称 concordant-down(对 snRNA 端 overclaim,
#       且与 UNC5C[snRNA↑/bulk ns]判法不对称);按 meta② 归 uninformative,形状仍标其 bulk 显著性。
#   disconfirm 是方向性的(只抓 snRNA↑/bulk↓,对应"↑ in HFpEF"假设)。
#   核心只建在两模态一致的分子(FLRT2)上;点形状 = bulk padj<0.05(把"同一把尺子"标到图上)。
#
# 数据源:leg1b_snRNA_compartment_vs_bulk.csv(snRNA 区室特异 + CI + FDR)
#         + leg1_hahn2021_bulk_targets.csv(只取 bulk padj;★ 不再用 call,它只编码 bulk 方向)
###############################################################################

source("/Users/franxy/Documents/博士课题 6.15/FLRT2&UNC5C/fig_common.R")
suppressPackageStartupMessages({ library(dplyr) })

## ---- 读取 + 合并 ---------------------------------------------------------
b <- read_result_csv("leg1b_snRNA_compartment_vs_bulk.csv")
# leg1b 中 FLRT2 出现两次(ligand_FLRT / ligand_UNC5),数值逐位相同、均 Fibroblast(已核),
# 去重安全;保留每基因一行
b <- b %>% distinct(gene, .keep_all = TRUE)

# ★ 专家修正:只取 bulk 的 padj —— 不再用 leg1$call。
#   call 只编码 bulk 方向+显著性,不含 snRNA;旧 classify 据此把"两模态一致↓"
#   (ADGRL1/ADGRL3)误标为 disconfirm(显著反向)= 真实双标 bug。
leg1 <- read_result_csv("leg1_hahn2021_bulk_targets.csv") %>%
  select(gene, bulk_padj = padj_HFpEF)

dat <- b %>%
  left_join(leg1, by = "gene") %>%
  filter(!is.na(snRNA_log2FC_CB), !is.na(bulk_log2FC))   # 去掉 snRNA 低表达 NA(如 TENM4,无 x 坐标)

## ---- ★ 跨模态联合判据(snRNA × bulk;pre-locked,同一把尺子) -------------
#  corroborate     = snRNA 显著↑ 且 bulk 显著↑              → FLRT2, TENM3
#  disconfirm      = snRNA 显著↑ 且 bulk 显著↓ (sign-flip)  → 仅 ADGRL2
#  concordant-down = bulk 显著↓ 且 snRNA 同向↓(降级,非反向)→ ADGRL1, ADGRL3
#  uninformative   = 任一端 ns / 低表达                      → UNC5C, UNC5B, TENM2 …
dat <- dat %>%
  mutate(
    snRNA_sig_up = snRNA_FDR_CB < 0.05 & snRNA_log2FC_CB > 0,
    bulk_sig_up  = !is.na(bulk_padj) & bulk_padj < 0.05 & bulk_log2FC > 0,
    bulk_sig_dn  = !is.na(bulk_padj) & bulk_padj < 0.05 & bulk_log2FC < 0,
    bulk_sig     = !is.na(bulk_padj) & bulk_padj < 0.05,
    cls = dplyr::case_when(                                # case_when 向量化,非标量 ifelse 坑
      snRNA_sig_up & bulk_sig_up ~ "corroborate (both up)",
      snRNA_sig_up & bulk_sig_dn ~ "disconfirm (snRNA up / bulk down)",
      TRUE                        ~ "uninformative (ns / low-expr)"   # 含 ADGRL1/3:bulk↓但 snRNA ns
    )
  )
cls_levels <- c("corroborate (both up)",
                "disconfirm (snRNA up / bulk down)",
                "uninformative (ns / low-expr)")
dat$cls <- factor(dat$cls, levels = cls_levels)
CLS_COL <- c("corroborate (both up)"             = "#D55E00",   # 橙:两模态一致↑
             "disconfirm (snRNA up / bulk down)" = "#2166AC",   # 深蓝:唯一 sign-flip(ADGRL2 独占)
             "uninformative (ns / low-expr)"     = "#999999")   # 灰:任一端 ns(含 ADGRL1/3)

# 标签(带区室);重点基因加粗显示
dat$lab <- sprintf("%s (%s)", dat$gene, substr(dat$compartment, 1, 3))
key <- c("FLRT2", "UNC5C", "ADGRL2", "TENM3", "UNC5B", "ADGRL1", "ADGRL3")
dat$is_key <- dat$gene %in% key

rng <- range(c(dat$snRNA_log2FC_CB, dat$bulk_log2FC, dat$CI_L, dat$CI_R), na.rm = TRUE)
pad <- 0.25 * diff(rng); lims <- c(rng[1]-pad, rng[2]+pad)

## ---- 画图 -----------------------------------------------------------------
p <- ggplot(dat, aes(snRNA_log2FC_CB, bulk_log2FC)) +
  # 象限参考线 + y=x 对角(两模态完全一致)
  annotate("rect", xmin = 0, xmax = lims[2], ymin = 0, ymax = lims[2],
           fill = "#D55E0011") +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey70") +
  # snRNA 95% CI(横向);用 geom_segment 避免 ggplot2 4.0 geom_errorbarh 弃用警告
  geom_segment(aes(x = CI_L, xend = CI_R, y = bulk_log2FC, yend = bulk_log2FC, color = cls),
               linewidth = 0.6, alpha = 0.8) +
  geom_point(aes(color = cls, shape = bulk_sig, size = is_key)) +
  scale_color_manual(name = "cross-modal call", values = CLS_COL, drop = FALSE) +
  # 点形状显式编码 bulk 显著性(实心=padj<0.05,空心=ns)→ 把"bulk padj"标到图上,自证同一把尺子
  scale_shape_manual(name = "bulk padj < 0.05", values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = "significant", `FALSE` = "ns")) +
  scale_size_manual(values = c(`FALSE` = 2.0, `TRUE` = 3.4), guide = "none") +
  coord_equal(xlim = lims, ylim = lims) +
  labs(
    title    = "Figure 3. Cross-modal validation: snRNA vs orthogonal bulk",
    subtitle = "x = compartment-specific snRNA log2FC (CellBender, 95% CI); y = Hahn-2021 bulk log2FC.",
    x = "snRNA log2FC (compartment, CellBender)",
    y = "bulk log2FC (Hahn 2021, HFpEF vs control)",
    caption = paste0(
      "Pre-locked joint rule, same ruler for all genes: class = snRNA(sig,dir) x bulk(sig,dir); shape = bulk padj<0.05; ",
      "disconfirm defined directionally vs the up-in-HFpEF hypothesis. Top-right shaded = concordant up.\n",
      "ADGRL2(CM) is the ONLY cross-modal sign-flip (snRNA up / bulk sig down) -> out of core. ",
      "ADGRL1/3: bulk down but snRNA ns -> uninformative (not reversals). UNC5C(Fib): snRNA up / bulk ns -> uninformative.")) +
  theme_pub(11) +
  # ggplot2 4.0:inside 图例用 legend.position="inside" + legend.position.inside
  theme(legend.position = "inside", legend.position.inside = c(0.99, 0.02),
        legend.justification = c(1, 0),
        legend.background = element_rect(fill = "#FFFFFFAA", color = NA))

# 关键基因标签:有 ggrepel 用防重叠版,否则回退普通 geom_text
if (requireNamespace("ggrepel", quietly = TRUE)) {
  p <- p + ggrepel::geom_text_repel(data = subset(dat, is_key),
             aes(label = lab, color = cls), fontface = "bold", size = 3.3,
             min.segment.length = 0, box.padding = 0.6, seed = 1, show.legend = FALSE)
} else {
  message("[!] 未装 ggrepel,标签用 geom_text(可能重叠)")
  p <- p + geom_text(data = subset(dat, is_key),
             aes(label = lab, color = cls), fontface = "bold", size = 3.1,
             vjust = -0.7, show.legend = FALSE)
}

save_fig(p, "fig3_crossmodal_scatter", w = 7.6, h = 7)

## ---- 控制台自检 -----------------------------------------------------------
cat("\n[Fig3] 跨模态联合判据(snRNA × bulk):\n")
print(dat %>% select(gene, compartment, snRNA_log2FC_CB, snRNA_FDR_CB,
                     bulk_log2FC, bulk_padj, cls),
      row.names = FALSE, digits = 3)
cat("\n关键(同一把尺子):ADGRL2 是唯一 cross-modal sign-flip → disconfirm(独占);",
    "\nADGRL1/3(bulk↓但 snRNA ns)、UNC5C/UNC5B/TENM2(任一端 ns)→ uninformative;FLRT2/TENM3 两模态一致↑ → corroborate。\n")
