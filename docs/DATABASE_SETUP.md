# EasyRNASeq — 数据库准备指南

本文档列出 EasyRNASeq 流程运行所需/可选的数据库文件，供提前下载准备。

---

## 目录结构建议

将所有数据库放在统一目录下：

```
/data/databases/
├── genomes/              # 参考基因组
│   └── athaliana/
│       ├── genome.fa
│       ├── genes.gtf
│       ├── gene.pep.fa
│       ├── hisat2_index/  # HISAT2 索引
│       └── star_index/    # STAR 索引
├── rrna/                  # rRNA 数据库
│   └── silva-138-99-seqs.fa
├── kegg/                  # KEGG 数据库 (可选)
├── go/                    # GO 数据库
│   └── go-basic.obo
├── nr/                    # NR 蛋白数据库 (可选)
│   └── nr.fasta
├── uniprot/               # UniProt (可选)
│   └── uniprot_sprot.fasta
├── pfam/                  # Pfam (可选)
│   └── Pfam-A.hmm
└── string/                # STRING PPI (可选)
    └── 9606.protein.links.v12.0.txt
```

---

## 1. 参考基因组 FASTA + GTF

| 物种 | 来源 | 下载链接 |
|------|------|---------|
| 拟南芥 | Ensembl Plants | `https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/arabidopsis_thaliana/dna/` |
| 人类 | Ensembl | `https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/` |
| 小鼠 | Ensembl | `https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/` |
| 水稻 | Ensembl Plants | `https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/oryza_sativa/dna/` |
| 玉米 | Ensembl Plants | `https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/zea_mays/dna/` |
| 其他 | Ensembl / NCBI | 根据物种查找 |

```bash
# 拟南芥示例
mkdir -p /data/databases/genomes/athaliana
cd /data/databases/genomes/athaliana

wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.59.gtf.gz

gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gunzip Arabidopsis_thaliana.TAIR10.59.gtf.gz

ln -s Arabidopsis_thaliana.TAIR10.dna.toplevel.fa genome.fa
ln -s Arabidopsis_thaliana.TAIR10.59.gtf genes.gtf
```

---

## 2. HISAT2 / STAR 索引

**何时需要**: 使用 `--align_tool hisat2` 或 `star` 时

索引可**预先构建** (推荐) 或**流程自动构建** (耗时)。

### 2.1 预先构建 HISAT2 索引

```bash
mkdir -p /data/databases/genomes/athaliana/hisat2_index

# 含剪接位点 (小基因组 < 3G)
hisat2_extract_splice_sites.py genes.gtf > splicesites.txt
hisat2_extract_exons.py genes.gtf > exons.txt
hisat2-build -p 8 --ss splicesites.txt --exon exons.txt genome.fa hisat2_index/genome.fa

# 仅 FASTA (大基因组 > 3G)
hisat2-build -p 4 genome.fa hisat2_index/genome.fa
```

### 2.2 预先构建 STAR 索引

```bash
mkdir -p /data/databases/genomes/athaliana/star_index

STAR --runMode genomeGenerate \
     --genomeDir star_index/ \
     --genomeFastaFiles genome.fa \
     --sjdbGTFfile genes.gtf \
     --sjdbOverhang 149 \
     --runThreadN 8
```

### 2.3 使用 pre-built ref_anno_dir

索引构建好后，通过 `--ref_anno_dir` 指向即可，流程自动识别：

```bash
EasyRNASeq2.py --ref_anno_dir /data/databases/genomes/athaliana \
               --species Arabidopsis_thaliana
```

`ref_anno_dir` 目录下自动识别的文件：

| 文件/目录 | 自动映射为 |
|-----------|----------|
| `*.fa` 或 `*.fasta` | `ref_fa` |
| `*.gtf` | `ref_gtf` |
| `*.pep.fa` | `ref_pep` |
| `hisat2_index/` | HISAT2 索引 |
| `star_index/` | STAR 索引 |

---

## 3. rRNA 数据库 (用于 rRNA 残留评估)

```bash
mkdir -p /data/databases/rrna
cd /data/databases/rrna

# Silva LSU + SSU rRNA (推荐)
wget https://www.arb-silva.de/fileadmin/silva_databases/release_138/Exports/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/release_138/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz

# 合并
zcat SILVA_138.1_*_NR99_tax_silva.fasta.gz > silva-138-rrna.fa

# 构建 Bowtie2 索引
bowtie2-build silva-138-rrna.fa silva_rrna
```

---

## 4. KEGG 数据库

### 4.1 KEGG REST API (免费，无需许可证)

流程默认使用 KEGG API，只需在配置中指定 `kegg_org`（物种代码）：

```bash
# 常见物种代码
# hsa = 人类, mmu = 小鼠, ath = 拟南芥, dme = 果蝇
# sce = 酵母, cel = 线虫, dre = 斑马鱼, rno = 大鼠
# osa = 水稻, zma = 玉米, gmx = 大豆, sly = 番茄

EasyRNASeq2.py ... --kegg_org ath
```

### 4.2 本地 KEGG 数据库 (可选，需许可证)

如已获取 KEGG FTP 许可证：

```bash
mkdir -p /data/databases/kegg
# 联系 KEGG 获取 ftp 账号后下载
# wget -r -np -nH --cut-dirs=2 ftp://ftp.genome.jp/pub/kegg/
```

---

## 5. GO 数据库

```bash
mkdir -p /data/databases/go
cd /data/databases/go

# GO 基础文件 (obo 格式)
wget https://purl.obolibrary.org/obo/go/go-basic.obo

# GO 注释文件 (以拟南芥为例)
wget https://ftp.ebi.ac.uk/pub/databases/GO/goa/ARATH/goa_arabidopsis.gaf.gz
gunzip goa_arabidopsis.gaf.gz
```

---

## 6. NR 蛋白数据库 (可选，用于功能注释)

```bash
mkdir -p /data/databases/nr
cd /data/databases/nr

# NCBI NR 数据库 (很大 ~100GB+)
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
gunzip nr.gz

# 构建 Diamond 索引 (比 BLAST 快 100 倍)
diamond makedb --in nr.fasta --db nr
```

---

## 7. UniProt / Swiss-Prot (可选)

```bash
mkdir -p /data/databases/uniprot
cd /data/databases/uniprot

wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
```

---

## 8. Pfam 蛋白结构域 (可选)

```bash
mkdir -p /data/databases/pfam
cd /data/databases/pfam

wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz

# 构建 hmmpress 索引
hmmpress Pfam-A.hmm
```

---

## 9. STRING PPI 数据库 (可选)

```bash
mkdir -p /data/databases/string
cd /data/databases/string

# 以人类 (9606) 为例
wget https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz
wget https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz
gunzip *.gz
```

---

## 10. 快速准备清单

按分析需求选择：

| 分析需求 | 必需数据库 | 大小 |
|---------|-----------|------|
| **基础 RNA-seq** (比对+定量+差异) | 基因组 FASTA + GTF | ~500M-3G |
| + rRNA 评估 | Bowtie2 索引 + Silva rRNA | ~1G |
| + GO 富集 | go-basic.obo (可选: GO 注释文件) | ~50M |
| + KEGG 富集 | KEGG API (免费) + kegg_org 参数 | 0 |
| + PPI 网络 | STRING 数据库 | ~500M |
| + NR 注释 | NR 数据库 + Diamond 索引 | ~150G |
| + Pfam 注释 | Pfam-A.hmm | ~2G |
| + UniProt 注释 | uniprot_sprot.fasta | ~100M |

---

## 11. 配置示例

```bash
# 完整配置示例
EasyRNASeq2.py \
    --ref_anno_dir /data/databases/genomes/athaliana \
    --species Arabidopsis_thaliana \
    --kingdom Plants \
    --kegg_org ath \
    --rrna_db /data/databases/rrna/silva-138-rrna.fa \
    --go_db /data/databases/go/go-basic.obo \
    --string_db /data/databases/string/3702.protein.links.v12.0.txt \
    --run_ppi \
    --sample_info ./samples.tsv \
    --compare_file ./compare.tsv \
    --jobs 8
```
