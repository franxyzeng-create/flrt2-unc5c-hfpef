###############################################################################
# Step1(验证):核心三联分子在【各自正确区室】的 snRNA pseudobulk 方向
#   配体查 sender 区室、受体查 receiver 区室(专家要求),与 bulk(leg1)配对
#
# 数据源:Hahn 2025 h5ad(数据集【内部】donor 对照,非外部整合LV — 满足专家前提e)
# 口径:完全复用 01_pseudobulk_positive_control.R(per-donor聚合, limma-voom,
#       ~disease+sex+pool, filterByExpr, TMM, MIN_NUCLEI=20, 双层CR/CB一致性)
# 扩展:细胞类型 2类→全9类(因 ADGRL 主导对 TENM2_ADGRL2 的区室是【心肌】非成纤维)
#
# 区室归属(来自 subsetCommunication 的 source→target,已确认):
#   FLRT2_UNC5C(UNC5头号对): sender=Fibroblast, receiver=Fibroblast
#   FLRT2_FLRT2(FLRT头号对): sender=receiver=Fibroblast
#   TENM2_ADGRL2(ADGRL头号对): sender=Cardiomyocyte, receiver=Cardiomyocyte(主)/Endothelial1/Fibroblast
#
# 核心交付:每个核心分子在【其正确区室】的 pseudobulk log2FC vs bulk log2FC 配对表
#   → 两向量一致性 = 最干净的跨模态验证指标(专家c)
#   判读三态对接 bulk(专家f):pseudobulk方向 × bulk方向 → corroborate/组分/降级
###############################################################################

suppressPackageStartupMessages({
  library(zellkonverter); library(SingleCellExperiment); library(SummarizedExperiment)
  library(Matrix); library(edgeR); library(limma)
})

`%||%` <- function(a,b) if (is.null(a)||length(a)==0) b else a   # null合并(前置定义)

## ---- 0. 配置 ----
h5      <- "/Users/franxy/Documents/博士课题 6.15/严谨/data/Hahn 人 HFpEF 单核 snRNA/HFpEF_snRNAseq.h5ad"
out_dir <- "/Users/franxy/Documents/博士课题 6.15/严谨/results"
MIN_NUCLEI <- 20; MIN_DONORS <- 3

# 核心三联 + 次级目标基因,及【各自要查的区室】(基于 source→target)
# role: ligand→查sender区室; receptor→查receiver区室
# 区室直接用 ontology code(避开 CellChat简称 vs h5ad全称 的映射歧义):
#   Fibroblast = CL:0002548 (fibroblast of cardiac tissue)
#   Cardiomyocyte = CL:0000746 (cardiac muscle cell)
CODE_FIB <- "CL:0002548"; CODE_CM <- "CL:0000746"
TARGET_GENES <- data.frame(
  gene = c("FLRT2","FLRT2",  "UNC5C","UNC5B", "TENM2","TENM3","TENM4", "ADGRL1","ADGRL2","ADGRL3"),
  role = c("ligand_FLRT","ligand_UNC5", "receptor_UNC5","receptor_UNC5",
           "ligand_ADGRL","ligand_ADGRL","ligand_ADGRL", "receptor_ADGRL","receptor_ADGRL","receptor_ADGRL"),
  # 该分子的主导区室code(查方向的地方)
  comp_code = c(CODE_FIB, CODE_FIB, CODE_FIB, CODE_FIB,
                CODE_CM, CODE_CM, CODE_CM, CODE_CM, CODE_CM, CODE_CM),
  comp_name = c("Fibroblast","Fibroblast","Fibroblast","Fibroblast",
                "Cardiomyocyte","Cardiomyocyte","Cardiomyocyte","Cardiomyocyte","Cardiomyocyte","Cardiomyocyte"),
  stringsAsFactors = FALSE)
genes_to_query <- unique(TARGET_GENES$gene)

# bulk(leg1)的 log2FC(HFpEF vs control,RVS),用于配对 — 来自 leg1 输出
BULK <- c(FLRT2=1.1362, UNC5C=-0.1587, UNC5B=0.2694, TENM2=-0.3060,
          TENM3=0.2829, TENM4=1.1840, ADGRL1=-0.2421, ADGRL2=-0.2963, ADGRL3=-0.5843)

## ---- 1. 读 h5ad ----
message("== 读入 h5ad…")
sce <- readH5AD(h5, use_hdf5 = TRUE, reader = "R")

# 自动提取全部 cell_type 编码 → label(不手写死)
ct_map <- unique(as.data.frame(colData(sce))[, c("cell_type","cell_type__ontology_label")])
ct_map <- ct_map[order(ct_map$cell_type), ]
message("== 全部细胞类型:"); print(ct_map, row.names=FALSE)

# 确认 donor→pool 唯一
dp <- unique(as.data.frame(colData(sce))[, c("donor_id","biosample_id")])
stopifnot(!any(duplicated(dp$donor_id)))

## ---- 2. pseudobulk 聚合(复用你的函数) ----
make_pseudobulk <- function(sce, ct_code) {
  idx <- which(as.character(sce$cell_type) == ct_code)
  if (length(idx) == 0) return(NULL)
  sub <- sce[, idx]; donor <- factor(as.character(sub$donor_id))
  ind <- sparseMatrix(i=seq_along(donor), j=as.integer(donor), x=1,
                      dims=c(length(donor), nlevels(donor)))
  colnames(ind) <- levels(donor)
  counts <- list()
  for (layer in c("raw_count_cellranger","raw_count_cellbender")) {
    m <- as(assay(sub, layer), "CsparseMatrix")
    pb <- as.matrix(m %*% ind)
    if (layer=="raw_count_cellbender") pb <- round(pb)
    counts[[layer]] <- pb
  }
  cd <- as.data.frame(colData(sub))
  meta <- cd[!duplicated(cd$donor_id), c("donor_id","disease","sex","biosample_id")]
  rownames(meta) <- meta$donor_id; meta <- meta[colnames(ind),,drop=FALSE]
  meta$n_nuclei <- as.integer(table(donor)[colnames(ind)])
  list(counts=counts, meta=meta)
}

## ---- 3. DE(复用你的 limma-voom,但返回目标基因的 logFC+CI) ----
run_de <- function(pb, layer) {
  meta <- pb$meta; keep_d <- meta$n_nuclei > MIN_NUCLEI; meta <- meta[keep_d,]
  counts <- pb$counts[[layer]][, rownames(meta), drop=FALSE]
  disease <- factor(ifelse(meta$disease=="MONDO_0005252","HFpEF","control"),
                    levels=c("control","HFpEF"))
  n_hf <- sum(disease=="HFpEF"); n_ct <- sum(disease=="control")
  if (n_hf < MIN_DONORS || n_ct < MIN_DONORS)
    return(list(tt=NULL, n_hf=n_hf, n_ct=n_ct))   # 功效不足,不强跑
  sex <- factor(meta$sex); pool <- droplevels(factor(meta$biosample_id))
  g <- rownames(counts); drop_g <- grepl("^MT-",g)|grepl("^RP[LS]",g)
  counts <- counts[!drop_g,,drop=FALSE]
  # design:若某区室sex/pool退化(单水平)则去掉该项,避免model.matrix报错
  terms <- "~ disease"
  if (nlevels(sex) > 1) terms <- paste(terms, "+ sex")
  if (nlevels(pool) > 1) terms <- paste(terms, "+ pool")
  design <- model.matrix(as.formula(terms))
  dge <- DGEList(counts); keep <- filterByExpr(dge, group=disease)
  dge <- dge[keep,,keep.lib.sizes=FALSE]; dge <- calcNormFactors(dge)
  v <- voom(dge, design); fit <- eBayes(lmFit(v, design))
  tt <- topTable(fit, coef="diseaseHFpEF", number=Inf, sort.by="none", confint=TRUE)
  list(tt=tt, n_hf=n_hf, n_ct=n_ct, design=terms)
}

## ---- 4. 对目标区室跑 pseudobulk DE,提取目标基因 ----
# 用 ontology code 直接跑(核心:Fibroblast + Cardiomyocyte)
comp_codes <- unique(TARGET_GENES[, c("comp_code","comp_name")])
message(sprintf("\n== 要查的区室:%s", paste(comp_codes$comp_name, collapse=", ")))

de_by_comp <- list()   # 以 comp_name 为键
for (k in seq_len(nrow(comp_codes))) {
  code <- comp_codes$comp_code[k]; cname <- comp_codes$comp_name[k]
  message(sprintf("\n---------- 区室 %s (%s) ----------", cname, code))
  if (!(code %in% as.character(sce$cell_type))) {
    message(sprintf("  [!] code %s 不在 h5ad cell_type 里,跳过", code)); next
  }
  pb <- make_pseudobulk(sce, code)
  if (is.null(pb)) { message("  无细胞,跳过"); next }
  de_cb <- run_de(pb, "raw_count_cellbender")
  de_cr <- run_de(pb, "raw_count_cellranger")
  message(sprintf("  可用 donor:HFpEF=%d, control=%d (design: %s)",
                  de_cb$n_hf, de_cb$n_ct, if(!is.null(de_cb$design)) de_cb$design else "NA"))
  de_by_comp[[cname]] <- list(cb=de_cb, cr=de_cr)
}

## ---- 5. 核心交付:目标基因在【正确区室】的 pseudobulk log2FC vs bulk ----
message("\n=========== 核心交付:区室特异 pseudobulk vs bulk 配对 ===========")
rows <- list()
for (i in seq_len(nrow(TARGET_GENES))) {
  g <- TARGET_GENES$gene[i]; comp <- TARGET_GENES$comp_name[i]; role <- TARGET_GENES$role[i]
  dec <- de_by_comp[[comp]]
  if (is.null(dec) || is.null(dec$cb$tt)) {
    rows[[length(rows)+1]] <- data.frame(gene=g, compartment=comp, role=role,
      snRNA_log2FC_CB=NA, CI_L=NA, CI_R=NA, snRNA_FDR_CB=NA, snRNA_log2FC_CR=NA,
      n_hf=if(!is.null(dec)) dec$cb$n_hf else NA, n_ct=if(!is.null(dec)) dec$cb$n_ct else NA,
      bulk_log2FC=BULK[[g]] %||% NA); next
  }
  tt_cb <- dec$cb$tt; tt_cr <- dec$cr$tt
  if (g %in% rownames(tt_cb)) {
    r_cb <- tt_cb[g, ]; r_cr <- if (!is.null(tt_cr) && g %in% rownames(tt_cr)) tt_cr[g,"logFC"] else NA
    rows[[length(rows)+1]] <- data.frame(
      gene=g, compartment=comp, role=role,
      snRNA_log2FC_CB=r_cb$logFC, CI_L=r_cb$CI.L, CI_R=r_cb$CI.R,
      snRNA_FDR_CB=r_cb$adj.P.Val, snRNA_log2FC_CR=r_cr,
      n_hf=dec$cb$n_hf, n_ct=dec$cb$n_ct, bulk_log2FC=BULK[[g]] %||% NA)
  } else {
    rows[[length(rows)+1]] <- data.frame(gene=g, compartment=comp, role=role,
      snRNA_log2FC_CB=NA, CI_L=NA, CI_R=NA, snRNA_FDR_CB=NA, snRNA_log2FC_CR=NA,
      n_hf=dec$cb$n_hf, n_ct=dec$cb$n_ct, bulk_log2FC=BULK[[g]] %||% NA)  # 被filterByExpr滤掉=低表达
  }
}
tab <- do.call(rbind, rows)
# 跨模态一致性判读(专家f三态)
tab$concord <- with(tab, ifelse(is.na(snRNA_log2FC_CB), "snRNA低表达/功效不足",
  ifelse(sign(snRNA_log2FC_CB)==sign(bulk_log2FC),
         ifelse(snRNA_log2FC_CB>0,"两模态一致↑(corroborate)","两模态一致↓(降级信号)"),
         ifelse(snRNA_log2FC_CB>0 & bulk_log2FC<0,"snRNA↑/bulk↓(查组分:他细胞拉低bulk?)",
                "snRNA↓/bulk↑(罕见,需查)"))))
print(tab, row.names=FALSE, digits=3)

## ---- 6. 关键基因专项解读 ----
message("\n=========== 关键判读 ===========")
say <- function(g) {
  r <- tab[tab$gene==g, ][1,]
  if (is.na(r$snRNA_log2FC_CB)) {
    message(sprintf("%s @ %s: snRNA 低表达/功效不足(n_hf=%s,n_ct=%s) — uninformative",
                    g, r$compartment, r$n_hf, r$n_ct)); return(invisible())
  }
  message(sprintf("%s @ %s(%s): snRNA log2FC=%+.3f [%.2f,%.2f] FDR=%.2e | bulk=%+.3f → %s",
                  g, r$compartment, r$role, r$snRNA_log2FC_CB, r$CI_L, r$CI_R,
                  r$snRNA_FDR_CB, r$bulk_log2FC, r$concord))
}
message("【FLRT2 — hub配体,表达层claim命门】")
say("FLRT2")
message("\n【UNC5C — hub受体(头号对FLRT2_UNC5C的receiver),交互层claim命门】")
say("UNC5C")
message("\n【ADGRL 受体 — bulk全下调,查心肌区室内在方向定降级与否】")
for (g in c("ADGRL2","ADGRL1","ADGRL3","TENM2")) say(g)
message("\n【其余核心】")
for (g in c("UNC5B","TENM3","TENM4")) say(g)

message("\n---- 判读规则(对接bulk,专家f)----")
message("· 配体查sender区室、受体查receiver区室(已按source→target归属)")
message("· snRNA↓+bulk↓ 一致 → 该分子非成纤维/区室内在上调,原CellChat信号是概率空间产物 → 降级")
message("· snRNA↑+bulk↓ → 查组分(他细胞高表达且HFpEF变少拉低bulk);但比例估计弱,组分解释只能写'可能'")
message("· snRNA↑+bulk↑ → 跨模态corroborate(强)")
message("· snRNA低表达/功效不足 → uninformative,看donor n,别把欠功效读成'内在不升'")
message("· 主交付=snRNA_log2FC vs bulk_log2FC 两向量一致性(配对表/散点)")

write.csv(tab, file.path(out_dir, "leg1b_snRNA_compartment_vs_bulk.csv"), row.names=FALSE)
message(sprintf("\n→ 已存:leg1b_snRNA_compartment_vs_bulk.csv"))
message("\n== sessionInfo =="); print(sessionInfo())
