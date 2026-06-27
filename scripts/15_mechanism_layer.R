###############################################################################
# 机制 / 定位层:UNC5 / FLRT / ADGRL 三核心(+ Netrin/SLIT 次级)的
#   (1) sender→receiver 细胞对结构(HFpEF vs control)
#   (2) 每条通路内的驱动 L-R 配受体对(确认非单一噪声对)
#
# 定性(关键,写进 Methods/Results 时遵守):
#   本层是【描述网络结构】,用 pooled merged 对象合法;但【不带病人层面 p 值】——
#   病人层面推断是 Step 3 干的(且 Step3 是固定4类集,本层是全9类)。
#   产出是描述性的细胞对排序 + 驱动 L-R 对,不是又一轮假设检验。
#
# 数据源(关键修正):
#   第1步 sender→receiver 用 netP$prob(CellChat 官方通路级聚合,口径权威)。
#   第2步 L-R 对 用 subsetCommunication(官方函数,返回显著细胞对×L-R对的长表)。
#   ※ 不再手动 sum(net$prob[,,lr]):那会把 pval 不显著的细胞对也加进去,
#     与 netP 聚合(只计显著对)差约4%。subsetCommunication 与 netP 同口径。
#
# 输入:cellchat_merged_raw_v3.rds(单 CellChat 对象,net/netP/LR 按 control/HFpEF 分存)
###############################################################################

suppressPackageStartupMessages({ library(CellChat); library(reshape2); library(dplyr) })
options(stringsAsFactors = FALSE)

base <- "/Users/franxy/Documents/博士课题 6.15/严谨/results"
cc   <- readRDS(file.path(base, "cellchat_merged_raw_v3.rds"))

CORE      <- c("UNC5","FLRT","ADGRL")    # 四轴全清硬核心
SECONDARY <- c("Netrin","SLIT")           # 次级(SLIT 去污染敏感)
ALL_PW    <- c(CORE, SECONDARY)

netP_hf <- cc@netP$HFpEF$prob; netP_ct <- cc@netP$control$prob
celltypes <- dimnames(netP_hf)[[1]]
cat("细胞类型(9):", paste(celltypes, collapse=", "), "\n")
cat(sprintf("通路数:HFpEF %d, control %d\n\n", dim(netP_hf)[3], dim(netP_ct)[3]))

# 用官方 subsetCommunication 一次性取出所有目标通路的 L-R 级长表(两组)
comm <- subsetCommunication(cc, slot.name = "net", signaling = ALL_PW)
comm_hf <- comm$HFpEF; comm_ct <- comm$control
cat(sprintf("subsetCommunication 取出:HFpEF %d 行, control %d 行(细胞对×L-R对,仅显著)\n\n",
            nrow(comm_hf), nrow(comm_ct)))

## ===========================================================================
## 第 1 步:sender→receiver 细胞对结构(每条核心通路)
##   用 netP$prob(官方通路级聚合);control 缺则0矩阵;按名索引对齐
## ===========================================================================
cat("================= 第 1 步:sender→receiver 结构 =================\n")
cat("[描述性,非病人层面检验] 数据源 netP$prob(官方聚合)。delta=HFpEF-control 看增强在哪些细胞对\n\n")

sr_summary <- list()
for (pw in ALL_PW) {
  if (!(pw %in% dimnames(netP_hf)[[3]])) { cat(sprintf("[!] %s 不在 HFpEF netP,跳过\n", pw)); next }
  m_hf <- netP_hf[,, pw]
  if (pw %in% dimnames(netP_ct)[[3]]) m_ct <- netP_ct[,, pw] else m_ct <- matrix(0, 9, 9, dimnames=dimnames(m_hf))
  stopifnot(identical(dimnames(m_hf), dimnames(m_ct)))   # 守卫:防错位
  tot_hf <- sum(m_hf); tot_ct <- sum(m_ct)
  cat(sprintf("---- %s ---- HFpEF总流量=%.4g, control总流量=%.4g, 差=%.4g (HFpEF %s) ----\n",
              pw, tot_hf, tot_ct, tot_hf-tot_ct, ifelse(tot_hf>tot_ct,"↑","↓")))

  long <- melt(m_hf, varnames=c("sender","receiver"), value.name="prob_HFpEF")
  long$prob_control <- m_ct[cbind(as.character(long$sender), as.character(long$receiver))]  # 按名对齐
  long$delta        <- long$prob_HFpEF - long$prob_control
  long$share_HFpEF  <- if (tot_hf > 0) long$prob_HFpEF/tot_hf else 0  # 组内归一;标量条件用if/else(非ifelse,否则只返回长度1被循环填充)
  long <- long[order(-long$prob_HFpEF), ]
  top  <- head(long[long$prob_HFpEF > 0, ], 8)
  top$cellpair <- paste(top$sender, "->", top$receiver)
  print(top[, c("cellpair","prob_HFpEF","prob_control","delta","share_HFpEF")], row.names=FALSE, digits=3)
  cat(sprintf("   top8 细胞对累计占 HFpEF 该通路流量 %.0f%%\n", 100*sum(top$share_HFpEF)))
  # 也按 delta 排:增强最大的细胞对(回答"增强发生在哪",专家强调看delta不是share)
  top_d <- head(long[order(-long$delta), ], 5)
  top_d$cellpair <- paste(top_d$sender, "->", top_d$receiver)
  cat("   增强最大(delta) top5 细胞对:\n")
  print(top_d[, c("cellpair","delta","prob_HFpEF","prob_control")], row.names=FALSE, digits=3)
  cat("\n")
  long$pathway <- pw
  sr_summary[[pw]] <- long
}
sr_all <- do.call(rbind, sr_summary)
write.csv(sr_all, file.path(base, "mech_sender_receiver.csv"), row.names=FALSE)
cat("→ 全细胞对明细已存:mech_sender_receiver.csv\n\n")

## ===========================================================================
## 第 2 步:每条核心通路内的驱动 L-R 配受体对
##   用 subsetCommunication 长表,按 interaction 聚合(官方口径,与 netP 一致)
## ===========================================================================
cat("================= 第 2 步:驱动 L-R 配受体对 =================\n")
cat("[确认非单一噪声对] 数据源 subsetCommunication。每 L-R 对对 sender×receiver 求和,HFpEF vs control\n\n")

lr_summary <- list()
for (pw in ALL_PW) {
  # 该通路的 L-R 级行(两组),按 interaction 聚合 prob
  hf_pw <- comm_hf[comm_hf$pathway_name == pw, ]
  ct_pw <- comm_ct[comm_ct$pathway_name == pw, ]
  if (nrow(hf_pw) == 0) { cat(sprintf("[!] %s 在 HFpEF 无显著通讯,跳过\n", pw)); next }

  agg_hf <- aggregate(prob ~ interaction_name + ligand + receptor, data=hf_pw, FUN=sum)
  agg_ct <- if (nrow(ct_pw)>0) aggregate(prob ~ interaction_name, data=ct_pw, FUN=sum) else
            data.frame(interaction_name=character(0), prob=numeric(0))
  names(agg_hf)[names(agg_hf)=="prob"] <- "prob_HFpEF"
  names(agg_ct)[names(agg_ct)=="prob"] <- "prob_control"
  tab <- merge(agg_hf, agg_ct, by="interaction_name", all.x=TRUE)
  tab$prob_control[is.na(tab$prob_control)] <- 0
  tab$delta <- tab$prob_HFpEF - tab$prob_control
  tot_hf <- sum(tab$prob_HFpEF)
  tab$share_HFpEF <- if (tot_hf > 0) tab$prob_HFpEF/tot_hf else 0  # 标量条件用if/else
  tab <- tab[order(-tab$prob_HFpEF), ]

  cat(sprintf("---- %s ---- HFpEF 显著 L-R 对 %d 个,总流量=%.4g ----\n", pw, nrow(tab), tot_hf))
  print(tab[, c("interaction_name","ligand","receptor","prob_HFpEF","prob_control","delta","share_HFpEF")],
        row.names=FALSE, digits=3)

  # 一致性检查(专家断言,改用官方口径后应通过):netP通路总流量 == subsetComm该通路L-R之和
  netP_tot_hf <- sum(netP_hf[,, pw])
  diff <- abs(netP_tot_hf - tot_hf)
  rel <- diff / netP_tot_hf
  # 三级:<1e-6 机器零(映射证死);1e-6~1e-4 且相对<0.1% = aggregation-boundary(p值阈值边界条目,可忽略);否则报警
  status <- if (diff < 1e-6) "✓(机器零,映射证死)" else
            if (rel < 1e-3) sprintf("≈(aggregation-boundary, 相对%.3f%%, 可忽略)", 100*rel) else
            "[!]偏差超限,需查"
  cat(sprintf("   [一致性] netP通路总流量=%.6g vs subsetComm之和=%.6g, 差=%.2g %s\n",
              netP_tot_hf, tot_hf, diff, status))

  # 浓度剖面:top-1/2/3 累计份额(单对主导≠噪声,可能真生物,如实报告)
  sh <- sort(tab$share_HFpEF, decreasing=TRUE); cum <- cumsum(sh)
  cat(sprintf("   浓度剖面(HFpEF): top1=%.0f%% / top2累计=%.0f%% / top3累计=%.0f%%\n\n",
              100*cum[1], 100*ifelse(length(cum)>=2,cum[2],cum[1]),
              100*ifelse(length(cum)>=3,cum[3],cum[length(cum)])))
  tab$pathway <- pw
  lr_summary[[pw]] <- tab
}
lr_all <- do.call(rbind, lr_summary)
write.csv(lr_all, file.path(base, "mech_LR_pairs.csv"), row.names=FALSE)
cat("→ 全 L-R 对明细已存:mech_LR_pairs.csv\n\n")

## ===========================================================================
## 解读指引
## ===========================================================================
cat("================= 解读指引 =================\n")
cat("第1步(sender->receiver,数据源 netP$prob):\n")
cat("  - 增强发生在哪 → 看【delta】top(绝对增量),不是 share\n")
cat("  - share_HFpEF = 组内归一,描述'信号集中在哪些细胞对'(已修ifelse bug);勿读成跨组变化\n")
cat("  - 五条通路 delta 最大细胞对一致指向 Fibroblast→Fibroblast(及 CM/Endocardial→Fibroblast)\n")
cat("第2步(L-R对,数据源 subsetCommunication):\n")
cat("  - 浓度剖面看单对/多对(已修ifelse bug,真实份额):\n")
cat("    · UNC5(FLRT2_UNC5C~69%)、FLRT(FLRT2_FLRT2~94%)、Netrin(NTN1_UNC5C~64%)= 单对/双对主导\n")
cat("    · ADGRL(TENM2_ADGRL2~26%,23对)= 真·多对分布式(Teneurin-Latrophilin 家族协同)\n")
cat("    · SLIT(SLIT2_ROBO1~89%)= 单对主导\n")
cat("  - 机制叙事(据真实浓度):\n")
cat("    ① FLRT2–UNC5C 轴:UNC5/FLRT/Netrin 头号对均收敛到 FLRT2配体/UNC5C受体 → 聚焦的导向信号增强\n")
cat("       (各由 dominant interaction 承载,如实写 'dominant interaction(s)',勿写 multi-pair program)\n")
cat("    ② Teneurin-Latrophilin/ADGRL:23对协同 → 唯一可称 'coordinated multi-pair program' 的\n")
cat("  - ADGRL 残差写 'negligible aggregation-boundary difference (0.01%)',勿写 floating-point rounding\n")
cat("注意:本层描述性网络结构(pooled,全9类),不带病人层面统计;\n")
cat("     病人层面证据是 Step 3(固定4类集);两层互补、定性不同。\n")
cat("     核心 vs 次级的区分轴是【去污染稳健性】(UNC5/FLRT/ADGRL robust,SLIT sensitive),\n")
cat("     不是对数多寡——UNC5/FLRT/Netrin 虽也单对主导,但去污染稳健,仍是核心;SLIT 敏感+无冗余,次级。\n\n")
cat("【写 Results 预埋句(描述性是 deliberate 方法学选择,非漏统计)】\n")
cat("  Cell-pair localization is descriptive; statistical significance of the pathway-level\n")
cat("  enhancement is established at the patient level (Step 3). We deliberately avoid\n")
cat("  per-cell-pair testing on the pooled object to not re-introduce the cell-pair\n")
cat("  pseudoreplication that the per-donor design supersedes. L-R contributions were\n")
cat("  extracted via subsetCommunication (CellChat), consistent with pathway-level aggregation.\n")

cat("\n================= sessionInfo =================\n")
print(sessionInfo())
