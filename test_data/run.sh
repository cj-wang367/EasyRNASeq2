#!/bin/bash
# ============================================================
# EasyRNASeq2 一键运行脚本
#
# 使用前请修改下面 "必填参数" 部分的路径
# 然后执行: bash run.sh
#
# 缓存说明:
#   .snakemake/  → 工作目录 (本地, 每个项目独立)
#   Conda 环境   → 流程目录 .snakemake/conda/ (多用户共享)
#   重置缓存:    rm -rf .snakemake/
#   完全重置:    rm -rf .snakemake/ results/ tmp/
# ============================================================
set -e
export PATH=/home/mycovirus/miniconda3/envs/snakemake/bin:$PATH

# ============================================================
# 1. 激活 Snakemake 环境
# ============================================================
conda activate snakemake

# 确认 snakemake 在 PATH 中
SNK=$(which snakemake 2>/dev/null || echo "/home/mycovirus/miniconda3/envs/snakemake/bin/snakemake")
echo "[INFO] Snakemake: $SNK"

# ============================================================
# 2. 必填参数 — 请根据你的数据修改
# ============================================================
REF_FA="GCF_000146945.2_ASM14694v2_genomic.fna"   # 参考基因组 FASTA
REF_GTF="genome.gtf"                               # 参考基因组注释 GTF
SPECIES="Sclerotinia_sclerotiorum"                 # 物种名
SAMPLE_INFO="sample_info.tsv"                      # 样本信息表
COMPARE_FILE="compare.tsv"                         # 比较组文件

# ============================================================
# 3. 可选参数 — 按需调整
# ============================================================

# --- 分析工具 ---
ALIGN_TOOL="hisat2"              # hisat2 / star / salmon
EXP_TOOL="featureCounts"         # featureCounts / stringtie / rsem
STRAND="unstranded"              # unstranded / firststrand / secondstrand
COUNT_NORM="fpkm"                # fpkm / tpm
BASES_DATA=6                     # 合同数据量 (G)
HISAT2_INDEX_ONLYFA=""           # 大基因组 (>3G) 加: --hisat2_index_onlyfa

# --- 差异分析 ---
DIFF_PVALUE=0.05                 # 显著性阈值
LOG_FC=1.0                       # log2 差异倍数阈值
DIFF_PTYPE="padj"               # padj (校正) / pvalue (原始)
BATCH_EFFECT=""                  # 启用批次校正: --batch_effect

# --- 富集分析 ---
KEGG_ORG=""                      # KEGG 物种代码 (如 ath/hsa/osa, 留空自动推断)
RUN_PPI=""                       # 启用 PPI: --run_ppi
RUN_WGCNA=""                     # 启用 WGCNA: --run_wgcna

# --- 运行控制 ---
JOBS=8                           # 最大并发任务数
OUTDIR="results"                 # 输出目录
TMPDIR="tmp"                     # 临时目录
FORCEALL=""                      # 强制重跑: --forceall
VERBOSE=""                       # 详细日志: --verbose
DRYRUN=""                        # 仅预览: --dryrun
VALIDATE=""                      # 仅校验: --validate

# --- 数据库路径 (留空则自动扫描) ---
KEGG_DB=""                       # KEGG 本地数据库
GO_DB=""                         # GO 数据库 (go-basic.obo)
NR_DB=""                         # NR 蛋白数据库
UNIPROT_DB=""                    # UniProt/Swiss-Prot
PFAM_DB=""                       # Pfam-A.hmm
STRING_DB=""                     # STRING PPI 数据库
EGGNOG_DB=""                     # eggNOG-mapper 数据库目录

# ============================================================
# 4. 拼接命令并运行
# ============================================================

CMD="python EasyRNASeq2.py \
    --ref_fa \"$REF_FA\" \
    --ref_gtf \"$REF_GTF\" \
    --species \"$SPECIES\" \
    --sample_info \"$SAMPLE_INFO\" \
    --compare_file \"$COMPARE_FILE\" \
    --align_tool \"$ALIGN_TOOL\" \
    --exp_tool \"$EXP_TOOL\" \
    --strand \"$STRAND\" \
    --count_norm \"$COUNT_NORM\" \
    --bases_data $BASES_DATA \
    --diff_pvalue $DIFF_PVALUE \
    --logFC $LOG_FC \
    --diff_ptype \"$DIFF_PTYPE\" \
    --jobs $JOBS \
    --outdir \"$OUTDIR\" \
    --tmpdir \"$TMPDIR\""

# 可选参数
[ -n "$HISAT2_INDEX_ONLYFA" ]  && CMD="$CMD $HISAT2_INDEX_ONLYFA"
[ -n "$KEGG_ORG" ]             && CMD="$CMD --kegg_org \"$KEGG_ORG\""
[ -n "$KEGG_DB" ]              && CMD="$CMD --kegg_db \"$KEGG_DB\""
[ -n "$GO_DB" ]                && CMD="$CMD --go_db \"$GO_DB\""
[ -n "$NR_DB" ]                && CMD="$CMD --nr_db \"$NR_DB\""
[ -n "$UNIPROT_DB" ]           && CMD="$CMD --uniprot_db \"$UNIPROT_DB\""
[ -n "$PFAM_DB" ]              && CMD="$CMD --pfam_db \"$PFAM_DB\""
[ -n "$STRING_DB" ]            && CMD="$CMD --string_db \"$STRING_DB\""
[ -n "$EGGNOG_DB" ]            && CMD="$CMD --eggnog_db \"$EGGNOG_DB\""
[ -n "$BATCH_EFFECT" ]         && CMD="$CMD $BATCH_EFFECT"
[ -n "$RUN_PPI" ]              && CMD="$CMD $RUN_PPI"
[ -n "$RUN_WGCNA" ]            && CMD="$CMD $RUN_WGCNA"
[ -n "$FORCEALL" ]             && CMD="$CMD $FORCEALL"
[ -n "$VERBOSE" ]              && CMD="$CMD $VERBOSE"
[ -n "$DRYRUN" ]               && CMD="$CMD $DRYRUN"
[ -n "$VALIDATE" ]             && CMD="$CMD $VALIDATE"

echo ""
echo "=========================================="
echo "  EasyRNASeq2 运行命令:"
echo "  $CMD"
echo "=========================================="
echo ""

eval $CMD
