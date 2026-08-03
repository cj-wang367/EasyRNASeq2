# EasyRNASeq2 — 二代有参转录组分析流程

基于 **Snakemake 9 + Conda** 的一站式 RNA-seq 分析流程，涵盖从原始数据到差异表达、功能富集、可变剪切、LncRNA 鉴定的全部分析。

---

## 快速开始

### 最简用法（用户只需修改 5 个参数）

```bash
conda activate snakemake

python EasyRNASeq2.py \
    --ref_fa genome.fna \
    --ref_gtf genome.gtf \
    --species Sclerotinia_sclerotiorum \
    --sample_info sample_info.tsv \
    --compare_file compare.tsv \
    --jobs 8
```

流程会自动：
- 生成配置文件
- 校验所有输入文件格式
- 运行完整分析
- 输出日志（含进度条，不再打印完整 shell 命令）
- 生成 HTML 运行报告 (`results/run_report.html`)
- 如有失败任务，生成错误摘要 (`results/run_error_summary.txt`)

### 一键脚本（可保存为 `run.sh` 反复使用）

```bash
# 复制模板
cp run.sh my_project.sh

# 编辑必填参数
vim my_project.sh    # 修改 REF_FA / REF_GTF / SPECIES / SAMPLE_INFO / COMPARE_FILE

# 运行
bash my_project.sh
```

`run.sh` 中所有可调参数都有中文注释，按需取消注释即可。

### 使用参考基因组目录（自动识别）

```bash
python EasyRNASeq2.py \
    --ref_anno_dir /data/genomes/rice \
    --species Oryza_sativa \
    --sample_info sample_info.tsv \
    --compare_file compare.tsv \
    --jobs 8
```

---

## 必选参数

| 参数 | 说明 | 示例 |
|------|------|------|
| `--ref_anno_dir` | **推荐**。参考基因组目录，自动查找 FASTA/GTF/HISAT2索引 | `/data/genomes/rice` |
| 或分别指定: | | |
| `--ref_fa` | 参考基因组 FASTA | `genome.fa` |
| `--ref_gtf` | 参考基因组注释 GTF（支持 `.gff/.gff3` 自动转换） | `genes.gtf` |
| `--species` | 物种名（用于自动匹配 KEGG/GO 数据库） | `Oryza_sativa` |
| `--sample_info` | 样本信息表（TSV 格式，默认 `sample_info.tsv`） | |
| `--compare_file` | 比较组文件（TSV 格式，默认 `compare.tsv`） | |

---

## 输入文件格式

### sample_info.tsv（制表符分隔）

```
group	sample	ngs
WT	WT_1	/path/to/WT_1_data
WT	WT_2	/path/to/WT_2_data
WT	WT_3	/path/to/WT_3_data
Treat	T_1	/path/to/T_1_data
Treat	T_2	/path/to/T_2_data
Treat	T_3	/path/to/T_3_data
```

| 列 | 说明 |
|----|------|
| `group` | 分组名（同一处理的生物学重复用相同的组名） |
| `sample` | 样本名（**唯一**，不重复） |
| `ngs` | 测序数据目录路径（目录下需包含 `*R1.fq.gz` 和 `*R2.fq.gz` 文件） |
| `batch` | （可选）批次信息，添加后可用 `--batch_effect` 启用批次校正 |

### compare.tsv

每行一对比较组，支持两种格式：

```
# 格式1（推荐）：TSV
Control	Treat

# 格式2（兼容）：下划线
Control_vs_Treat
```

`_vs_` 前为对照组，后为实验组。FC = 实验组 / 对照组。

---

## 可选参数

### 分析工具

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--align_tool` | `hisat2` | 比对工具：`hisat2` / `star` / `salmon` |
| `--exp_tool` | `featureCounts` | 定量工具：`featureCounts` / `stringtie` / `rsem` |
| `--strand` | `unstranded` | 链特异性：`unstranded` / `firststrand` / `secondstrand` |
| `--count_norm` | `fpkm` | 标准化方式：`fpkm` / `tpm` |
| `--bases_data` | `6` | 合同数据量（G），影响数据截取量 |
| `--hisat2_index_onlyfa` | 否 | 大基因组(>3G)建议开启 |

> **Salmon 模式** (`--align_tool salmon`)：自动禁用比对质控、新转录本发现、可变剪切、LncRNA 模块。适合仅需快速定量+差异分析。

### 差异分析

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--diff_pvalue` | `0.05` | 显著性阈值 |
| `--logFC` | `1.0` | 差异倍数阈值（log2） |
| `--diff_ptype` | `padj` | 显著性类型：`padj` / `pvalue` |
| `--batch_effect` | 关闭 | 启用批次校正（需 `batch` 列） |

> 自动选择方法: 每组重复≥2 → DESeq2；某组仅1个样本 → edgeR

### 富集分析 / PPI / WGCNA

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--kegg_org` | 自动推断 | KEGG 物种代码（如 `ath` / `hsa` / `osa`） |
| `--run_ppi` | 关闭 | 启用 STRING PPI 蛋白互作网络 |
| `--run_wgcna` | 关闭 | 启用 WGCNA 共表达网络分析 |

### 数据库路径（留空自动扫描）

| 参数 | 用途 |
|------|------|
| `--go_db` | GO 数据库 (`go-basic.obo`) |
| `--kegg_db` | KEGG 本地数据库目录 |
| `--string_db` | STRING PPI 数据库 |
| `--rrna_db` | rRNA 参考序列 (Silva) |
| `--nr_db` | NR 蛋白数据库 |
| `--uniprot_db` | UniProt/Swiss-Prot |
| `--pfam_db` | Pfam-A.hmm |
| `--eggnog_db` | eggNOG-mapper 数据库目录 |

自动扫描路径：流程目录 `data/database/`、当前目录 `database/`、`~/easynaseq/databases/`

### 运行控制

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--jobs` | `8` | 最大并发任务数 |
| `--outdir` | `results` | 结果输出目录 |
| `--tmpdir` | `tmp` | 临时文件目录 |
| `--conda_prefix` | — | Conda 环境共享目录（复用加速） |
| `--dryrun` | — | 仅预览任务不执行 |
| `--validate` | — | 仅校验输入参数和文件 |
| `--unlock` | — | 解锁工作目录 |
| `--forceall` | — | 强制重跑所有任务 |
| `--verbose` | — | 详细日志（打印完整 shell 命令，调试用） |

> **默认日志已简化**：只显示进度条 + 错误信息，不再打印每个任务的 Shell 命令。需要调试时加 `--verbose`。

---

## 运行示例

```bash
# 基础运行（HISAT2 + featureCounts）
python EasyRNASeq2.py \
    --ref_fa genome.fa --ref_gtf genes.gtf \
    --species Oryza_sativa --jobs 8

# 大基因组（>3G）
python EasyRNASeq2.py \
    --ref_anno_dir /data/genomes/wheat \
    --species Triticum_aestivum \
    --hisat2_index_onlyfa --jobs 4

# Salmon 快速模式
python EasyRNASeq2.py \
    --ref_anno_dir /data/genomes/human \
    --species Homo_sapiens \
    --align_tool salmon --jobs 8

# 启用全部高级分析 + 批次校正
python EasyRNASeq2.py \
    --ref_fa genome.fa --ref_gtf genes.gtf \
    --species Oryza_sativa \
    --run_wgcna --run_ppi --batch_effect --jobs 8

# 修改差异阈值
python EasyRNASeq2.py \
    --ref_fa genome.fa --ref_gtf genes.gtf \
    --species Oryza_sativa \
    --diff_pvalue 0.01 --logFC 2.0

# 仅校验输入文件
python EasyRNASeq2.py \
    --ref_fa genome.fa --ref_gtf genes.gtf \
    --species Oryza_sativa --validate

# 异常退出后解锁
python EasyRNASeq2.py --unlock
```

---

## 输出结果目录

```
results/
├── 01.raw_data/                     # 原始数据统计
│   └── data_stat_all.xls
├── 02.clean_data/                   # 质控后数据 (fastp + FastQC)
│   ├── {sample}/
│   └── multiqc_report.html
├── 03.align/                        # 比对结果
│   ├── {sample}/
│   ├── align_stat_all.txt
│   ├── geneBody_coverage/
│   ├── read_distribution/
│   └── junction_saturation/
├── 04.transcripts/                  # 新转录本
│   ├── all_transcripts.fa / .gtf
│   ├── all.gene_trans_map.txt
│   ├── new_transcripts/
│   │   ├── new_transcripts.fa / .pep.fa / .cds.fa
│   │   ├── new_code_stat.pdf / .xls
│   │   └── gffcompare_stats.txt
│   ├── Number_stats.pdf
│   └── anno/ (go / kegg / kog / pfam)
├── 05.expression/                   # 表达定量
│   ├── gene/
│   │   ├── gene_count.xls / gene_fpkm.xls / gene_tpm.xls
│   │   ├── box_density/
│   │   ├── pca/
│   │   ├── corr/
│   │   └── heatmap/
│   └── transcript/
│       └── transcript_count.xls / transcript_cpm.xls
├── 06.de_enrich/                    # 差异分析 + 富集
│   ├── PCA.pdf
│   ├── correlation_heatmap.pdf
│   ├── DE_feature_counts.pdf
│   ├── venn_upset/
│   └── {Ctrl_vs_Treat}/
│       ├── {comp}.DE_results.All.xls
│       ├── {comp}.DE_results.Significant.xls
│       ├── {comp}.volcano.pdf / {comp}_MA.pdf
│       └── enrich/
│           ├── go/
│           ├── kegg/
│           └── gsea/
├── 07.as/                           # 可变剪切 (rMATS)
├── 08.lnc/                          # LncRNA 分析
└── 09.wgcna/                        # WGCNA 共表达网络
```

---

## 常见问题

### 1. 改代码后仍报同样的错误

首先确认 `.snakemake/` 缓存在**工作目录**（运行脚本的目录），不在流程目录。

```bash
rm -rf .snakemake/     # 清除工作目录中的 Snakemake 缓存
```

> **注意**：仅删除 `.snakemake/` 不会触发重跑（输出文件仍在 `results/` 中）。
> 如需完全重置：`rm -rf .snakemake/ results/ tmp/`

### 2. 为什么修改代码后全部重算？

Snakemake 跟踪代码哈希值。修改以下文件会导致相关规则全部重跑：
- `envs/*.yaml` — Conda 环境哈希变化 → 使用该环境的所有规则失效
- `rules/*.smk` — 规则代码变化 → 对应规则失效
- `Snakefile` — 工作流定义变化

这是 Snakemake 的正常行为，确保结果可复现。

### 3. 缓存目录说明

| 目录 | 位置 | 用途 | 共享? |
|------|------|------|-------|
| `.snakemake/` | **工作目录** | 元数据、锁、代码哈希 | 每个项目独立 |
| Conda 环境 | **流程目录** `.snakemake/conda/` | 复用预构建的 conda 环境 | 多用户共享 |

- 删除工作目录的 `.snakemake/` = 清除运行时元数据（安全）
- 删除流程目录的 `.snakemake/conda/` = 清除所有 conda 环境（所有用户需重新构建）
- 流程目录本身保持干净，每次运行的缓存都在各自的工作目录中

### 4. 运行速度慢

使用 `--conda_prefix` 复用预构建 conda 环境：
```bash
python EasyRNASeq2.py ... --conda_prefix /path/to/shared_envs
```

### 5. 日志太长

默认已使用 `--quiet progress` 模式，只显示进度条和错误。如需要完整 shell 命令用于调试，加上 `--verbose`。

### 6. StringTie 报错 "no valid ID found for GFF record"

GTF/GFF 文件缺少 `gene_id` 或 `transcript_id` 属性：
```bash
gffread -T input.gff -o output.gtf
```

### 7. Snakemake 锁定

```bash
python EasyRNASeq2.py --unlock
```

### 8. 查看运行后错误摘要

```bash
cat results/run_error_summary.txt      # 文本摘要
# 或
python scripts/parse_run_log.py EasyRNASeq_run_XXXXXX.log   # 详细解析
# 或
open results/run_report.html           # HTML 可视化报告
```

---

## 依赖

- **Python**: ≥ 3.9
- **Snakemake**: ≥ 8.0（推荐 9.x）
- **Conda / Mamba**: 用于环境管理
- **Snakemake 环境路径**: `/home/mycovirus/miniconda3/envs/snakemake/bin/snakemake`
- **所有生物信息学工具**: 通过 Conda 自动安装，无需手动配置
