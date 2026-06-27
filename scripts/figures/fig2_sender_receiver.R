###############################################################################
# Fig 2 — UNC5 pathway sender→receiver structure (交互层, CellChat-inferred)
#
# 讲什么:UNC5 通路在 HFpEF 的细胞对通讯增强(delta = HFpEF - control)集中在哪;
#         Fibroblast→Fibroblast 单区室自分泌为主导(支持锁定核心的自分泌定位)。
#
# ★ 定性边界(必须遵守,写进 caption):
#   - 这是 CellChat【推断】的通讯,非测量;FLRT-UNC5 属接触依赖 trans-adhesion,
#     CellChat 不验证物理接触 → 副标题写 "inferred; contact not verified"。
#   - 本层 pooled、全 9 类、描述性,不带病人层面 p 值(病人层面在 Step 3 / Fig1)。
#   - 措辞用 "single-compartment autocrine",勿写 synergy / program。
#
# 数据源(主,自包含):mech_sender_receiver.csv(已含 81 个 sender×receiver / 通路)
# 备选(native,见文末):cellchat_*_final.rds → netVisual_aggregate(signaling="UNC5")
###############################################################################

source("/Users/franxy/Documents/博士课题 6.15/FLRT2&UNC5C/fig_common.R")
suppressPackageStartupMessages({ library(dplyr) })

PW <- "UNC5"

sr <- read_result_csv("mech_sender_receiver.csv")
sr <- sr[sr$pathway == PW, ]
if (nrow(sr) == 0) stop("mech_sender_receiver.csv 中无 ", PW, " 通路")

# 固定 9 类顺序(缺失补全为 0,保证 9x9 完整网格)
present <- union(unique(sr$sender), unique(sr$receiver))
lev <- intersect(CELLTYPE_LEVELS, present)
sr$sender   <- factor(sr$sender,   levels = lev)
sr$receiver <- factor(sr$receiver, levels = lev)

## ============================================================================
## (A) sender × receiver delta 热图(最稳健,bulletproof)
## ============================================================================
hot <- sr %>% filter(sender == "Fibroblast", receiver == "Fibroblast")
# 断言:高亮的 Fib→Fib 确为本通路 delta 的 argmax(防换通路/数据后硬编码静默错标)
stopifnot(nrow(hot) == 1, isTRUE(all.equal(hot$delta, max(sr$delta, na.rm = TRUE))))
lim <- max(abs(sr$delta), na.rm = TRUE)

pA <- ggplot(sr, aes(receiver, sender, fill = delta)) +
  geom_tile(color = "white", linewidth = 0.6) +
  # 给 Fib→Fib 加黑框突出
  geom_tile(data = hot, fill = NA, color = "black", linewidth = 1.4) +
  scale_fill_gradient2(name = expression(bold(Delta~prob~(HFpEF-control))),
                       low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-lim, lim)) +
  scale_x_discrete(position = "top", expand = c(0,0)) +
  scale_y_discrete(limits = rev(lev), expand = c(0,0)) +
  labs(title = paste0("Figure 2. ", PW, " pathway: cell-pair communication change (HFpEF - control)"),
       subtitle = "CellChat-inferred; contact not verified. Box = Fibroblast->Fibroblast (autocrine).",
       x = "Receiver", y = "Sender",
       caption = paste0("Descriptive, pooled object, all 9 cell types; no patient-level p (that is Step 3 / Fig 1).\n",
                        "Effect magnitudes are CellChat-INFERRED and small in absolute terms (max delta ~0.003 prob units).")) +
  theme_pub(11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0),
        axis.ticks = element_blank(), panel.border = element_blank())

save_fig(pA, "fig2a_sender_receiver_heatmap", w = 7.4, h = 6)

## ============================================================================
## (B) igraph 网络:仅画增强(delta>0)的边,Fib→Fib 自环突出
## ============================================================================
have_igraph <- requireNamespace("igraph", quietly = TRUE)
if (have_igraph) {
  library(igraph)
  ud <- sr %>% filter(delta > 0) %>% mutate(sender = as.character(sender),
                                            receiver = as.character(receiver))
  v  <- data.frame(name = lev)
  # 节点大小 ∝ 该类作为 sender 的总增量(outgoing delta)
  outdelta <- tapply(ud$delta, ud$sender, sum)
  v$out <- as.numeric(outdelta[v$name]); v$out[is.na(v$out)] <- 0
  v$size <- 8 + 30 * (v$out / max(v$out, 1e-9))   # max(,1e-9) 防该通路无正边时除零
  v$color <- "#CCCCCC"; v$color[v$name == "Fibroblast"] <- ROLE_COLORS[["locked core (FLRT2->UNC5C)"]]

  g <- graph_from_data_frame(ud[, c("sender","receiver","delta")], directed = TRUE, vertices = v)
  E(g)$width <- 0.6 + 9 * (E(g)$delta / max(E(g)$delta, 1e-9))   # 防除零
  # 只高亮 Fibroblast→Fibroblast 自分泌边(橙);其余边含他类自环一律灰,避免过度强调
  el <- ends(g, E(g))
  is_fibloop <- el[, 1] == "Fibroblast" & el[, 2] == "Fibroblast"
  E(g)$color <- ifelse(is_fibloop, ROLE_COLORS[["locked core (FLRT2->UNC5C)"]], "#9E9E9E66")
  set.seed(1)
  # 圆形布局,Fibroblast 排在首位(order = 按期望次序给出的顶点 id 向量)
  desired <- c("Fibroblast", setdiff(lev, "Fibroblast"))
  lay <- layout_in_circle(g, order = match(desired, V(g)$name))

  svglite::svglite(file.path(FIG_DIR, "fig2b_unc5_network.svg"), width = 7, height = 7)  # Arial 矢量(无需 XQuartz)
  par(mar = c(1,1,3,1), family = "Arial")
  plot(g, layout = lay,
       vertex.size = V(g)$size, vertex.color = V(g)$color,
       vertex.frame.color = "white",
       vertex.label = V(g)$name, vertex.label.color = "black",
       vertex.label.cex = 0.8, vertex.label.dist = 1.6,
       edge.arrow.size = 0.5, edge.curved = 0.18,
       edge.loop.angle = 0.6,
       main = "Supplementary Fig. S1. UNC5 pathway enhanced communication (HFpEF > control)")
  mtext("Edge width = delta prob; loop = Fibroblast autocrine. CellChat-inferred, contact not verified.",
        side = 1, cex = 0.7, col = "grey40")
  dev.off()

  png(file.path(FIG_DIR, "fig2b_unc5_network.png"), width = 4200, height = 4200, res = 600)  # 投稿 600dpi
  par(mar = c(1,1,3,1), family = "Arial")
  plot(g, layout = lay,
       vertex.size = V(g)$size, vertex.color = V(g)$color,
       vertex.frame.color = "white",
       vertex.label = V(g)$name, vertex.label.color = "black",
       vertex.label.cex = 0.8, vertex.label.dist = 1.6,
       edge.arrow.size = 0.5, edge.curved = 0.18, edge.loop.angle = 0.6,
       main = "Supplementary Fig. S1. UNC5 pathway enhanced communication (HFpEF > control)")
  dev.off()
  message("  saved: fig2b_unc5_network.{svg,png}")
} else {
  message("[!] 未装 igraph,跳过网络图(B)。热图(A)已生成。")
}

## ---- 控制台自检:UNC5 通路 top 细胞对(按 delta) -------------------------
cat("\n[Fig2] UNC5 pathway sender->receiver, top by delta:\n")
print(sr %>% arrange(desc(delta)) %>%
        select(sender, receiver, prob_HFpEF, prob_control, delta, share_HFpEF) %>%
        head(8), row.names = FALSE, digits = 3)

## ============================================================================
## 备选:CellChat 官方画法(若要 native circle plot;需 ~16GB RAM 逐个加载)
## ----------------------------------------------------------------------------
## library(CellChat)   # 重启后必须先 library 才能用 netVisual_*
## ccH <- readRDS(file.path(RESULTS_DIR, "cellchat_HFpEF_final.rds"))   # ~180MB
## ccC <- readRDS(file.path(RESULTS_DIR, "cellchat_control_final.rds")) # ~260MB
## pdf(file.path(FIG_DIR, "fig2c_native_unc5_circle.pdf"), width = 10, height = 5)
## par(mfrow = c(1,2))
## netVisual_aggregate(ccC, signaling = "UNC5", layout = "circle", main.gradient = TRUE)
## netVisual_aggregate(ccH, signaling = "UNC5", layout = "circle")
## dev.off()
## 注:UNC5 在 control 可能通讯很弱,native 图两侧不可直接比"绝对粗细";
##     差异比较仍以(A)的 delta 热图为准。
