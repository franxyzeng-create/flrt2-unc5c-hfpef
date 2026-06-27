###############################################################################
# 腿 1:Hahn 2021 bulk RV-septal 验证(三角验证第一腿)— corroboration-only 测试
#
# ★定性(专家定调,写论文必遵守):
#   本腿的意义【不在样本独立】(与snRNA队列同组、部分供者重叠,本就不独立),
#   而在【测量模态正交】:单核通讯推断 vs bulk组织DE,即使患者重叠,问的问题/测的层次不同。
#   bulk = 细胞比例 × 每细胞表达 的乘积 → 即使阳性也只说"组织内总量多",
#   不等于"通讯程序激活";且无法区分"每细胞表达↑"vs"该类细胞变多(组分漂移/纤维化)"。
#
# ★可证伪规则(事前锁死,防"怎么都赢"):稀释只把信号衰减向0、不翻号,故——
#   显著上调          = corroboration(强化)
#   ns / 低表达 / NA  = uninformative(与稀释一致,权重转向其他腿,非阴性)
#   显著【下调】      = DISCONFIRMING(稀释解释不了的sign-flip,硬核心须回炉)
#   阈值主要施加在 hub(FLRT2/UNC5C)和核心三联(UNC5/FLRT/ADGRL);次级SLIT/Netrin低权重。
#
# ★BMI 是核心混杂(非polish):此队列中位BMI~41(病态肥胖),Hahn2021原文证明
#   HFpEF上调基因几乎全被BMI单独解释、主上调是代谢/线粒体非ECM。故本腿(未校正)
#   仅看方向;"是否HFpEF特异 vs 肥胖伴随"由 leg2(脚本17,BMI校正)正面回答。
#
# 数据:作者官方DESeq2结果 DEG_allresultsRVS.csv(权威口径,不重算)
#   .pef=HFpEF vs ctrl; l2fc.pef方向; padj.pef BH校正p; external_gene_name=symbol
###############################################################################

options(stringsAsFactors = FALSE)
base <- "/Users/franxy/Documents/博士课题 6.15/严谨/data/Hahn2020_HFpEF_bulk/baderzone-HFpEF_2020-71fc3a4/Data"
rvs  <- read.csv(file.path(base, "DEG_allresultsRVS.csv"), check.names=FALSE)

## ---- 补丁2:列名先 dump + 断言,缺列立刻报(不跑一半才炸) ----
cat("列名:\n"); print(names(rvs))
need <- c("external_gene_name","l2fc.pef","padj.pef","pvalue.pef",
          "control_median","hfpef_median","l2fc.ref","padj.ref")
miss_col <- setdiff(need, names(rvs))
if (length(miss_col) > 0) cat("\n[!] 缺列:", paste(miss_col, collapse=", "), "\n")
stopifnot(length(miss_col) == 0)

# ---- 目标基因:机制层核心+次级全部 L-R 分子,并标注 hub / 核心三联 / 次级 ----
ligands  <- c("FLRT2","FLRT1","FLRT3","NTN1","NTN4","TENM1","TENM2","TENM3","TENM4","SLIT2")
receptors<- c("UNC5A","UNC5B","UNC5C","UNC5D","ADGRL1","ADGRL2","ADGRL3","ROBO1","ROBO2","DCC","NTRK2","DSCAM")
targets  <- c(ligands, receptors)
target_role <- setNames(c(rep("ligand",length(ligands)), rep("receptor",length(receptors))), targets)
HUB  <- c("FLRT2","UNC5C")                         # 机制层收敛核心,主要判读对象
CORE <- c("FLRT2","FLRT1","FLRT3","UNC5A","UNC5B","UNC5C","UNC5D","ADGRL1","ADGRL2","ADGRL3","TENM1","TENM2","TENM3","TENM4")  # UNC5/FLRT/ADGRL
tier <- function(g) if (g %in% HUB) "★hub" else if (g %in% CORE) "core(UNC5/FLRT/ADGRL)" else "secondary(SLIT/Netrin)"
pw_map <- list(
  FLRT2="UNC5/FLRT/ADGRL(hub配体)", FLRT1="UNC5/FLRT/ADGRL", FLRT3="UNC5/FLRT/ADGRL",
  NTN1="Netrin", NTN4="Netrin", SLIT2="SLIT",
  TENM1="ADGRL", TENM2="ADGRL", TENM3="ADGRL", TENM4="ADGRL",
  UNC5A="UNC5", UNC5B="UNC5/Netrin", UNC5C="UNC5/Netrin(hub受体)", UNC5D="UNC5/Netrin",
  ADGRL1="ADGRL", ADGRL2="ADGRL", ADGRL3="ADGRL", ROBO1="SLIT", ROBO2="SLIT",
  DCC="Netrin", NTRK2="Netrin", DSCAM="Netrin")

cat("\n=========== 腿 1:Hahn 2021 bulk RV-septal(corroboration-only)===========\n")
cat(sprintf("目标:%d 基因(配体%d+受体%d);hub=FLRT2,UNC5C;核心三联=UNC5/FLRT/ADGRL\n\n",
            length(targets), length(ligands), length(receptors)))

# ---- 提取 ----
sub <- rvs[rvs$external_gene_name %in% targets,
           c("external_gene_name","l2fc.pef","padj.pef","pvalue.pef",
             "control_median","hfpef_median","l2fc.ref","padj.ref")]
names(sub) <- c("gene","log2FC_HFpEF","padj_HFpEF","pval_HFpEF",
                "ctrl_median","hfpef_median","log2FC_HFrEF","padj_HFrEF")

## ---- 补丁4:重复 symbol 断言(防 hub 循环 recycling) ----
if (any(duplicated(sub$gene))) {
  cat("[!] 重复 symbol:", paste(sub$gene[duplicated(sub$gene)], collapse=", "),
      "→ 按 padj 最小去重\n")
  sub <- sub[order(sub$gene, sub$padj_HFpEF), ]
  sub <- sub[!duplicated(sub$gene), ]
}

sub$role    <- target_role[sub$gene]
sub$tier    <- sapply(sub$gene, tier)
sub$pathway <- sapply(sub$gene, function(g) pw_map[[g]])
sub$direction <- ifelse(sub$log2FC_HFpEF > 0, "up", "down")

## ---- 补丁1+3:三态判读(* / ns / NA低表达),区分"功效不足"与"测了不显著" ----
# padj=NA = DESeq2独立过滤/Cook剔除 = 低表达功效不足(≠ns);低median也降权
sub$call <- ifelse(is.na(sub$padj_HFpEF), "NA(低表达/被过滤)",
             ifelse(sub$padj_HFpEF < 0.05,
                    ifelse(sub$log2FC_HFpEF > 0, "*up(corroborate)", "*DOWN(DISCONFIRM)"),
                    "ns(uninformative)"))
## ---- 补丁5:median 接近0 标记(方向是噪声,降权) ----
sub$lowexpr <- ifelse(pmax(sub$ctrl_median, sub$hfpef_median) < 10, "lowExpr", "")

# 排序:hub优先,再显著上调,再log2FC
sub$tier_rank <- ifelse(sub$tier=="★hub",0, ifelse(grepl("core",sub$tier),1,2))
sub <- sub[order(sub$tier_rank, -(grepl("corroborate",sub$call)), -sub$log2FC_HFpEF), ]

cat("---- 目标基因 HFpEF vs control(RVS),三态判读 ----\n")
print(sub[, c("gene","tier","role","pathway","direction","log2FC_HFpEF","padj_HFpEF","call","ctrl_median","hfpef_median","lowexpr")],
      row.names=FALSE, digits=3)

# ---- 未匹配基因:被低计数过滤(原文做过 low count + blood filtering),≠阴性 ----
missing <- setdiff(targets, rvs$external_gene_name)
if (length(missing) > 0)
  cat(sprintf("\n[未在表中] %d 个(被低计数/血污染过滤,非阴性,uninformative):%s\n",
              length(missing), paste(missing, collapse=", ")))

## ---- 汇总(补丁1:全部 na.rm,且三态分开计数) ----
cat("\n=========== 汇总(三态)===========\n")
up_sig   <- sub$gene[!is.na(sub$padj_HFpEF) & sub$padj_HFpEF<0.05 & sub$log2FC_HFpEF>0]
down_sig <- sub$gene[!is.na(sub$padj_HFpEF) & sub$padj_HFpEF<0.05 & sub$log2FC_HFpEF<0]
ns_g     <- sub$gene[!is.na(sub$padj_HFpEF) & sub$padj_HFpEF>=0.05]
na_g     <- sub$gene[is.na(sub$padj_HFpEF)]
cat(sprintf("检出%d / 共%d 目标基因:\n", nrow(sub), length(targets)))
cat(sprintf("  显著上调(corroborate): %d  %s\n", length(up_sig), paste(up_sig, collapse=", ")))
cat(sprintf("  显著下调(DISCONFIRM!): %d  %s\n", length(down_sig), paste(down_sig, collapse=", ")))
cat(sprintf("  ns(uninformative):     %d  %s\n", length(ns_g), paste(ns_g, collapse=", ")))
cat(sprintf("  NA低表达(被过滤):      %d  %s\n", length(na_g), paste(na_g, collapse=", ")))
cat(sprintf("  未在表中(过滤掉):      %d  %s\n", length(missing), paste(missing, collapse=", ")))

# ---- hub + 核心三联 专项判读(按 assay-aware 三档) ----
cat("\n=========== hub / 核心三联 判读(三角验证权重)===========\n")
for (g in c(HUB, setdiff(CORE, HUB))) {
  if (g %in% sub$gene) {
    r <- sub[sub$gene==g, ]
    verdict <- if (is.na(r$padj_HFpEF)) "uninformative(低表达,权重转其他腿)" else
               if (r$padj_HFpEF<0.05 & r$log2FC_HFpEF>0) "✓corroborate(BMI校正后须仍在才算数→leg2)" else
               if (r$padj_HFpEF<0.05 & r$log2FC_HFpEF<0) "✗DISCONFIRM(硬核心须回炉!)" else
               "uninformative(ns,与稀释一致)"
    star <- if (g %in% HUB) "★" else " "
    cat(sprintf("%s%-7s %-22s log2FC=%6.3f padj=%9.2e %s\n",
                star, g, r$pathway, r$log2FC_HFpEF, r$padj_HFpEF, verdict))
  } else {
    star <- if (g %in% HUB) "★" else " "
    cat(sprintf("%s%-7s 未在表中(被过滤,uninformative)\n", star, g))
  }
}

cat("\n---- 解读规则(事前锁死,照表判读)----\n")
cat("· 显著上调 → corroborate(但须 leg2 BMI校正后仍在,否则只是肥胖签名)\n")
cat("· ns/低表达/NA/未在表 → uninformative(与细胞稀释一致,非阴性,权重转 leg2/leg3)\n")
cat("· 显著下调 → DISCONFIRMING(稀释解释不了,硬核心须回炉重审)\n")
cat("· 主看 hub(FLRT2/UNC5C)+核心三联;次级SLIT/Netrin低权重\n")
cat("· bulk阳性也有上限:无法区分'每细胞表达↑'vs'成纤维变多(组分漂移)';不等于通讯激活\n")

## ---- 存表 ----
dir.create("/Users/franxy/Documents/博士课题 6.15/严谨/results", showWarnings=FALSE, recursive=TRUE)
write.csv(sub, "/Users/franxy/Documents/博士课题 6.15/严谨/results/leg1_hahn2021_bulk_targets.csv", row.names=FALSE)
cat("\n→ 已存:leg1_hahn2021_bulk_targets.csv\n")

cat("\n================= sessionInfo =================\n")
print(sessionInfo())
