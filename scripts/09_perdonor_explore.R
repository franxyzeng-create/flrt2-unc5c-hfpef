###############################################################################
# Step 3(探查):per-donor CellChat —— 摸清病人层面分布,为最终检验定阈值
#
# 目的:Step 3 是推断主干(把分析单元从"细胞对"提升到"病人",根除伪重复)。
#       但纳入哪些细胞类型、哪些通路、最少多少核/多少病人,这些阈值必须靠
#       真实数据定,不能拍脑袋。本脚本只建对象 + 统计分布,【不做最终检验】。
#
# 回答三个问题:
#   Q1 每病人每细胞类型的核数分布 → 哪些细胞类型够做 per-donor
#   Q2 每病人能建出几个细胞类型(min.cells 过滤后)→ 病人间一致性
#   Q3 每通路在多少病人里能检测到(非零信息流)→ 哪些通路够做病人层面检验
#
# 关键设计(专家):
#   - 输入 raw_count_cellranger(与主线 raw 臂一致;去污染稳健性已单独做完)
#   - 每病人单独建 CellChat、computeCommunProb,只提取通路强度,然后丢弃对象省内存
#   - nboot=1 提速:探查只用 prob(信息流强度),prob 与 nboot 无关(nboot 只影响 pval)
#   - 单病人内【不做组间比较】(就一个病人),只提取连续的 pathway 信息流强度
#   - 分析单元 = 每病人每通路的信息流强度(连续值),绝不是 per-donor 的 p 值
#
# 注意:本脚本【有意】让各病人各保留自己≥10核的类型,目的是诊断 Q2(每病人类型数
#       的组间分布)——这正是最终检验必须"固定公共类型集"的依据(类型数若组间不同,
#       sum(prob) 会因可用 sender-receiver 对数不同而系统性偏差)。最终检验会改成
#       所有病人用同一固定类型集 + 绝对/相对双版本,本探查不做检验、只摸分布。
###############################################################################

suppressPackageStartupMessages({
  library(zellkonverter); library(SingleCellExperiment); library(Matrix)
  library(CellChat)
})

options(stringsAsFactors = FALSE)
future::plan("sequential")
set.seed(1)

out_dir <- "/Users/franxy/Documents/博士课题 6.15/严谨/results"
h5      <- "/Users/franxy/Documents/博士课题 6.15/严谨/data/Hahn 人 HFpEF 单核 snRNA/HFpEF_snRNAseq.h5ad"

KEEP <- c(
  "CL:0000746" = "Cardiomyocyte", "CL:0002548" = "Fibroblast",
  "CL:0010008" = "Endothelial1",  "CL:4033076" = "Endothelial2",
  "CL:0000669" = "Pericyte",      "CL:0000359" = "VSMC",
  "CL:0000235" = "Macrophage",    "CL:0000542" = "Lymphocyte",
  "CL:0002350" = "Endocardial"
)

## ---- 读入 -----------------------------------------------------------------
message("== 读入 h5ad…")
sce <- readH5AD(h5, use_hdf5 = TRUE, reader = "R")
idx <- which(as.character(sce$cell_type) %in% names(KEEP))
sce <- sce[, idx]
celltype <- factor(KEEP[as.character(sce$cell_type)], levels = unname(KEEP))
donor    <- as.character(sce$donor_id)
group    <- ifelse(sce$disease == "MONDO_0005252", "HFpEF", "control")

# donor → group 映射(每病人一个组)
donor_group <- tapply(group, donor, function(x) x[1])
cat(sprintf("病人数: %d(HFpEF %d, control %d)\n",
            length(donor_group), sum(donor_group=="HFpEF"), sum(donor_group=="control")))

counts <- as(assay(sce, "raw_count_cellranger"), "CsparseMatrix")
rownames(counts) <- rownames(sce); colnames(counts) <- colnames(sce)

## ===========================================================================
## Q1:每病人每细胞类型的核数矩阵
## ===========================================================================
cat("\n========== Q1:每病人每细胞类型核数 ==========\n")
ct_by_donor <- table(donor, celltype)
cat("(行=病人, 列=细胞类型;只展示汇总统计)\n\n")

# 每个细胞类型:有多少病人拥有 ≥10 / ≥20 / ≥30 核
cat("各细胞类型在多少病人里达到核数阈值(ge10/ge20/ge30=≥10/≥20/≥30核):\n")
thr_tab <- sapply(colnames(ct_by_donor), function(ct) {
  v <- ct_by_donor[, ct]
  c(ge10_donors = sum(v >= 10),
    ge20_donors = sum(v >= 20),
    ge30_donors = sum(v >= 30),
    median_n    = median(v))
})
print(t(thr_tab))

# 按组分别看(HFpEF / control 各有多少病人达标)
cat("\n各细胞类型分组达标病人数(ge20=≥20核):\n")
donor_names <- rownames(ct_by_donor)
dg <- donor_group[donor_names]
grp_tab <- sapply(colnames(ct_by_donor), function(ct) {
  v <- ct_by_donor[, ct]
  c(HFpEF_ge20   = sum(v >= 20 & dg == "HFpEF"),
    control_ge20 = sum(v >= 20 & dg == "control"),
    HFpEF_ge10   = sum(v >= 10 & dg == "HFpEF"),
    control_ge10 = sum(v >= 10 & dg == "control"))
})
print(t(grp_tab))

# 专项:Endocardial 是固定集的唯一开关(5类必需里最小),单独看它每病人核数分布
cat("\n========== Endocardial 专项(固定集的唯一开关)==========\n")
if ("Endocardial" %in% colnames(ct_by_donor)) {
  endo <- ct_by_donor[, "Endocardial"]
  cat("Endocardial 每病人核数分布(全体):\n"); print(summary(as.numeric(endo)))
  cat(sprintf("HFpEF 病人: 达标 ge10=%d/%d, ge20=%d/%d\n",
              sum(endo>=10 & dg=="HFpEF"), sum(dg=="HFpEF"),
              sum(endo>=20 & dg=="HFpEF"), sum(dg=="HFpEF")))
  cat(sprintf("control 病人: 达标 ge10=%d/%d, ge20=%d/%d\n",
              sum(endo>=10 & dg=="control"), sum(dg=="control"),
              sum(endo>=20 & dg=="control"), sum(dg=="control")))
  cat("决策:若多数病人ge20→用ge20(给UNC5/Netrin/FLRT更稳trimmed mean);\n")
  cat("      若ge20砍太多→ge10保n;若ge10都把任一组压到<~12→退4类集(去Endocardial)。\n")
}

## ===========================================================================
## Q2 + Q3:逐病人建 CellChat,提取每通路信息流强度
## ===========================================================================
cat("\n========== Q2+Q3:逐病人建对象、提取通路强度 ==========\n")
cat("(每病人单独建 CellChat,min.cells=10;提取后丢弃对象省内存)\n")

donor_list <- names(donor_group)
# 存:每病人的 (通路 → 信息流强度)
strength_list <- list()
ncells_per_type_list <- list()
ntypes_per_donor <- integer(0)

for (d in donor_list) {
  cells_d <- colnames(counts)[donor == d]
  lab_d   <- droplevels(celltype[donor == d])
  ntab    <- table(lab_d)
  # 只保留该病人内 ≥10 核的细胞类型(CellChat min.cells 会滤,这里先记录)
  keep_types <- names(ntab)[ntab >= 10]
  ntypes_per_donor[d] <- length(keep_types)
  ncells_per_type_list[[d]] <- ntab

  # 少于 2 个细胞类型无法算通讯(没有细胞对),跳过提取
  if (length(keep_types) < 2) {
    strength_list[[d]] <- setNames(numeric(0), character(0))
    cat(sprintf("  病人 %s(%s):仅 %d 类≥10核,跳过\n", d, donor_group[d], length(keep_types)))
    next
  }

  # 子集到保留的细胞类型
  cells_keep <- cells_d[lab_d %in% keep_types]
  lab_keep   <- droplevels(lab_d[lab_d %in% keep_types])
  mat_d      <- counts[, cells_keep, drop = FALSE]

  # 建 CellChat(单病人,只为提取强度)
  data_n <- normalizeData(mat_d, scale.factor = 1e4, do.log = TRUE)
  meta_d <- data.frame(labels = lab_keep, samples = factor("s1"),
                       row.names = colnames(mat_d))
  cc <- tryCatch({
    cc <- createCellChat(object = data_n, meta = meta_d, group.by = "labels")
    cc@DB <- CellChatDB.human
    cc <- subsetData(cc)
    cc <- identifyOverExpressedGenes(cc)
    cc <- identifyOverExpressedInteractions(cc)
    cc <- computeCommunProb(cc, type = "truncatedMean", trim = 0.1,
                            population.size = TRUE, seed.use = 1, nboot = 1,
                            raw.use = TRUE)   # nboot=1:探查只用 prob(与nboot无关),提速
    cc <- filterCommunication(cc, min.cells = 10)
    cc <- computeCommunProbPathway(cc)
    cc
  }, error = function(e) { cat(sprintf("  病人 %s 建对象报错: %s\n", d, conditionMessage(e))); NULL })

  if (is.null(cc)) { strength_list[[d]] <- setNames(numeric(0), character(0)); next }

  # 提取每通路的总信息流强度:netP$prob 是 [sender × receiver × pathway],按通路加总
  prob <- cc@netP$prob
  if (is.null(prob) || length(dim(prob)) != 3) {
    strength_list[[d]] <- setNames(numeric(0), character(0))
  } else {
    pathways <- dimnames(prob)[[3]]
    strength <- apply(prob, 3, sum)             # 每通路:所有 sender-receiver 概率之和
    names(strength) <- pathways
    strength_list[[d]] <- strength
  }
  cat(sprintf("  病人 %s(%s):%d 类, 检出 %d 条通路\n",
              d, donor_group[d], length(keep_types), length(strength_list[[d]])))
}

## ---- Q2 汇总:每病人能建几个细胞类型 -------------------------------------
cat("\n--- Q2:每病人保留的细胞类型数(≥10核)分布 ---\n")
print(summary(ntypes_per_donor))
cat("按组:\n")
print(tapply(ntypes_per_donor[donor_list], donor_group[donor_list], summary))

# 专家①混杂核查(注意:这是【叙事】用,不是【决策】用):
# 固定公共类型集是【结构性必须】的,无论下面这个 p 多少——因为 sum(prob) 在不同
# 数目/身份的 sender-receiver 对上加总,本就不是跨病人可比的同一个量。
#   p<0.05 → 额外拿到"类型数组间确实不同、混杂确实在偏倚"的实锤,写进 manuscript;
#   p≥0.05 → 照样固定集(为可比性+功效),只是另加一句"类型数组间无差异、未偏倚比较"。
# 别让不显著的 p 把我们劝退回"各病人各建各的"。
cat("\n>> 混杂核查(叙事用,非决策):每病人类型数 HFpEF vs control\n")
nt_hf <- ntypes_per_donor[donor_list][donor_group[donor_list]=="HFpEF"]
nt_ct <- ntypes_per_donor[donor_list][donor_group[donor_list]=="control"]
cat(sprintf("   HFpEF 类型数: 中位 %.1f, 均值 %.2f\n", median(nt_hf), mean(nt_hf)))
cat(sprintf("   control 类型数: 中位 %.1f, 均值 %.2f\n", median(nt_ct), mean(nt_ct)))
cat(sprintf("   Wilcoxon p = %.3f\n", suppressWarnings(wilcox.test(nt_hf, nt_ct)$p.value)))
cat("   解读:无论此 p 多少,最终检验都固定公共类型集(结构性必须);\n")
cat("        p<0.05 → 混杂实锤(写进论文);p≥0.05 → 仍固定集 + 一句'未偏倚'宽心。\n")

## ---- Q3 汇总:每通路在多少病人里检出 -------------------------------------
cat("\n--- Q3:每通路在多少病人里检出(非零信息流)---\n")
all_pathways <- unique(unlist(lapply(strength_list, names)))
# 构造 通路 × 病人 的强度矩阵(缺失=0)
strength_mat <- matrix(0, nrow = length(all_pathways), ncol = length(donor_list),
                       dimnames = list(all_pathways, donor_list))
for (d in donor_list) {
  s <- strength_list[[d]]
  if (length(s)) strength_mat[names(s), d] <- s
}

# 每通路:在多少病人检出(非零),分组统计
hf_donors <- donor_list[donor_group[donor_list] == "HFpEF"]
ct_donors <- donor_list[donor_group[donor_list] == "control"]
detect_summary <- data.frame(
  pathway       = all_pathways,
  n_donors_detected = rowSums(strength_mat > 0),
  n_HFpEF_detected  = rowSums(strength_mat[, hf_donors, drop=FALSE] > 0),
  n_control_detected= rowSums(strength_mat[, ct_donors, drop=FALSE] > 0)
)
detect_summary <- detect_summary[order(-detect_summary$n_donors_detected), ]
cat(sprintf("共检出 %d 条通路;按检出病人数排序(前40):\n", nrow(detect_summary)))
print(head(detect_summary, 40), row.names = FALSE)

# 头条通路检出情况
HEADLINE <- c("EPHA","EPHB","ADGRL","UNC5","Netrin","SEMA5","SLIT","NCAM","NRXN","FLRT",
              "LAMININ","COLLAGEN","FN1")
cat("\n--- 头条通路在病人层面的检出情况 ---\n")
hl <- detect_summary[detect_summary$pathway %in% HEADLINE, ]
hl <- hl[match(HEADLINE, hl$pathway), ]
hl <- hl[!is.na(hl$pathway), ]
print(hl, row.names = FALSE)

## ---- 保存(供下一步最终检验用)-------------------------------------------
saveRDS(list(
  strength_mat   = strength_mat,       # 通路 × 病人 信息流强度
  donor_group    = donor_group,
  ntypes_per_donor = ntypes_per_donor,
  ct_by_donor    = ct_by_donor,
  detect_summary = detect_summary
), file.path(out_dir, "perdonor_explore.rds"))

cat(sprintf("\n== 完成。已存 %s\n", file.path(out_dir, "perdonor_explore.rds")))
cat("\n下一步:据 Q1/Q2/Q3 分布定阈值(纳入哪些细胞类型/通路、最少检出病人数),\n")
cat("       再做病人层面 Wilcoxon(19 HFpEF vs 24 control)+ BH。把本输出贴回讨论定阈值。\n")

cat("\n========== sessionInfo ==========\n")
print(sessionInfo())
