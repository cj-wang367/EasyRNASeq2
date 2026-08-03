# ============================================================
# EasyRNASeq v2 - 二代有参转录组分析流程
# 主 Snakefile
#
# 运行方式:
#   python EasyRNASeq2.py --ref_fa genome.fa --ref_gtf genes.gtf ...
#   或
#   snakemake --snakefile Snakefile --configfile config/config.yaml
# ============================================================

import os
import sys
import shutil
import subprocess as _sp
import yaml
import pandas as pd

# ---- 流程根目录 (用于 conda 环境路径解析) ----
# 注意: Snakemake 中 __file__ 不指向 Snakefile，必须用 workflow.basedir
PIPELINE_DIR = workflow.basedir


def run_r_script(r_code):
    """通过 subprocess 运行 R 脚本, 替代 shell heredoc (解决 Linux 兼容性).

    用法: 在 rule 的 run: 块中调用 run_r_script(r_code)
    """
    import subprocess as _rsp
    result = _rsp.run(
        ["Rscript", "-"],
        input=r_code,
        text=True,
        capture_output=True,
    )
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        raise Exception(f"Rscript failed (exit {result.returncode})")
    # R 的 message/warning 输出到 stdout, 打印以便调试
    if result.stdout:
        print(result.stdout, end="")

# ---- 加载配置 ----
# 优先使用 EasyRNASeq2.py 生成的配置文件
_generated_config = "config/config.generated.yaml"
if os.path.exists(_generated_config):
    configfile: _generated_config
    print(f"[INFO] 使用自动生成的配置文件: {_generated_config}")
else:
    configfile: "config/config.yaml"
    print(f"[INFO] 使用默认配置文件: config/config.yaml")

# ---- 预检: ref_anno_dir 自动识别 ----
_ref_anno_dir = config["paths"].get("ref_anno_dir")
if _ref_anno_dir and os.path.isdir(str(_ref_anno_dir)):
    import glob as _g
    _ad = str(_ref_anno_dir)
    print(f"[INFO] ref_anno_dir 指定: {_ad}")
    # 自动找 FASTA
    _fa_files = sorted(_g.glob(os.path.join(_ad, "*.fa")) + _g.glob(os.path.join(_ad, "*.fasta")))
    for _f in _fa_files:
        _bn = os.path.basename(_f).lower()
        if "transcript" in _bn: continue
        if "pep" in _bn: continue
        if not config["paths"].get("ref_fa"):
            config["paths"]["ref_fa"] = _f
            print(f"  [auto] ref_fa ← {_f}")
            break
    # 自动找 GTF
    _gtf_files = sorted(_g.glob(os.path.join(_ad, "*.gtf")) + _g.glob(os.path.join(_ad, "*.gtf.gz")))
    if _gtf_files and not config["paths"].get("ref_gtf"):
        config["paths"]["ref_gtf"] = _gtf_files[0]
        print(f"  [auto] ref_gtf ← {_gtf_files[0]}")
    # 自动找 pep
    for _f in _fa_files:
        if "pep" in os.path.basename(_f).lower() and not config["paths"].get("ref_pep"):
            config["paths"]["ref_pep"] = _f
            print(f"  [auto] ref_pep ← {_f}")
    # 自动找 HISAT2 索引
    _hisat2_idx_dir = os.path.join(_ad, "hisat2_index")
    if os.path.isdir(_hisat2_idx_dir) and not config["alignment"].get("index"):
        _ht2 = _g.glob(os.path.join(_hisat2_idx_dir, "*.1.ht2"))
        if _ht2:
            config["alignment"]["index"] = _hisat2_idx_dir
            print(f"  [auto] hisat2 索引 ← {_hisat2_idx_dir}")
    # 自动找 STAR 索引
    _star_idx_dir = os.path.join(_ad, "star_index")
    if os.path.isdir(_star_idx_dir) and not config["alignment"].get("index"):
        if os.path.isfile(os.path.join(_star_idx_dir, "Genome")):
            config["alignment"]["index"] = _star_idx_dir
            print(f"  [auto] star 索引 ← {_star_idx_dir}")

# ---- 预检: 验证必填项 ----
_errors = []

# 路径解析辅助函数: 相对路径基于 CWD 解析
def _resolve_path(p):
    if p and not os.path.isabs(str(p)):
        return os.path.abspath(str(p))
    return str(p) if p else None

# 必填路径
_ref_fa = config["paths"].get("ref_fa")
_ref_gtf = config["paths"].get("ref_gtf")
_sample_info = config["paths"].get("sample_info")
_compare_file = config["paths"].get("compare_file")

if not _ref_fa:
    _errors.append("paths.ref_fa 未设置! 请通过 EasyRNASeq2.py --ref_fa 指定或在 config.yaml 中填写")
else:
    _ref_fa_resolved = _resolve_path(_ref_fa)
    if not os.path.exists(_ref_fa_resolved):
        _errors.append(f"参考基因组 FASTA 不存在: {_ref_fa_resolved}  (原始值: {_ref_fa})")

if not _ref_gtf:
    _errors.append("paths.ref_gtf 未设置! 请通过 EasyRNASeq2.py --ref_gtf 指定或在 config.yaml 中填写")
else:
    _ref_gtf_resolved = _resolve_path(_ref_gtf)
    if not os.path.exists(_ref_gtf_resolved):
        _errors.append(f"参考基因组 GTF 不存在: {_ref_gtf_resolved}  (原始值: {_ref_gtf})")

if not _sample_info:
    _errors.append("paths.sample_info 未设置!")
else:
    _si_resolved = _resolve_path(_sample_info)
    if not os.path.exists(_si_resolved):
        _errors.append(f"样本信息表不存在: {_si_resolved}  (原始值: {_sample_info})")

if not _compare_file:
    _errors.append("paths.compare_file 未设置!")
else:
    _cf_resolved = _resolve_path(_compare_file)
    if not os.path.exists(_cf_resolved):
        _errors.append(f"比较组文件不存在: {_cf_resolved}  (原始值: {_compare_file})")

# 验证样本信息表格式 (使用解析后的路径)
_si_path = config["paths"]["sample_info"]  # 已被上面更新为绝对路径
if _si_path and os.path.exists(str(_si_path)):
    try:
        _test_df = pd.read_csv(_si_path, sep="\t", comment="#")
        _required_cols = ["group", "sample", "ngs"]
        _missing_cols = [c for c in _required_cols if c not in _test_df.columns]
        if _missing_cols:
            _errors.append(f"样本信息表缺少必需列: {_missing_cols} (需要: group, sample, ngs)")
    except Exception as _e:
        _errors.append(f"无法解析样本信息表: {_e}")

# 预检 GTF 格式: 检查 gene_id / transcript_id 属性是否存在
# StringTie 要求 GTF 记录必须有这些属性，否则报 "no valid ID found for GFF record"
_gtf_path = config["paths"].get("ref_gtf")
if _gtf_path and os.path.exists(str(_resolve_path(_gtf_path))):
    _gtf_resolved = str(_resolve_path(_gtf_path))
    _gid_ok, _tid_ok = False, False
    try:
        import re
        _checked = 0
        with open(_gtf_resolved) as _gf:
            for _line in _gf:
                if _line.startswith('#'): continue
                _parts = _line.strip().split('\t')
                if len(_parts) < 9: continue
                _checked += 1
                if not _gid_ok and re.search(r'gene_id\s+"?([^\s;"]+)"?', _parts[8]):
                    _gid_ok = True
                if not _tid_ok and re.search(r'transcript_id\s+"?([^\s;"]+)"?', _parts[8]):
                    _tid_ok = True
                if _checked >= 500 or (_gid_ok and _tid_ok):
                    break
        if not _gid_ok:
            _errors.append(f"GTF 文件中未检测到 gene_id 属性! StringTie/featureCounts 将无法正常工作。请检查 GTF 格式: {_gtf_resolved}")
        if not _tid_ok:
            _errors.append(f"GTF 文件中未检测到 transcript_id 属性! StringTie 组装将失败。请检查 GTF 格式: {_gtf_resolved}")
    except Exception as _e:
        _errors.append(f"读取 GTF 文件失败: {_gtf_resolved} ({_e})")

# GFF → GTF 自动转换: 如果用户提供了 .gff/.gff3 文件, 自动用 gffread 转换为 GTF
# 大部分生物信息学工具 (StringTie, featureCounts, rMATS) 严格要求 GTF 格式
_gtf_val = str(config["paths"].get("ref_gtf", ""))
_gff_exts = (".gff", ".gff3", ".gff.gz", ".gff3.gz")
if _gtf_val and _gtf_val.lower().endswith(_gff_exts):
    _tmp_dir = _resolve_path(config["paths"].get("tmp_dir", "tmp"))
    _gff_resolved = str(_resolve_path(_gtf_val))
    _gtf_converted = os.path.join(str(_tmp_dir), "ref", os.path.basename(_gff_resolved).replace(".gz","").replace(".gff3",".gtf").replace(".gff",".gtf"))
    os.makedirs(os.path.dirname(_gtf_converted), exist_ok=True)
    # gffread 可在 Snakefile 解析阶段直接调用 (它通常已通过 conda 安装)
    _gffread = shutil.which("gffread")
    if _gffread:
        print(f"\n  [INFO] 检测到 GFF 格式输入, 自动转换为 GTF...")
        print(f"  [INFO]   gffread -T {_gff_resolved} -o {_gtf_converted}")
        try:
            _sp.run([_gffread, "-T", _gff_resolved, "-o", _gtf_converted],
                    check=True, capture_output=True, text=True, timeout=300)
            config["paths"]["ref_gtf"] = _gtf_converted
            print(f"  [OK]   GFF → GTF 转换成功: {_gtf_converted}")
            # 重新检查转换后的 GTF
            _gid_ok2, _tid_ok2 = False, False
            with open(_gtf_converted) as _gf2:
                for _line in _gf2:
                    if _line.startswith('#'): continue
                    _parts = _line.strip().split('\t')
                    if len(_parts) < 9: continue
                    if not _gid_ok2 and re.search(r'gene_id\s', _parts[8]): _gid_ok2 = True
                    if not _tid_ok2 and re.search(r'transcript_id\s', _parts[8]): _tid_ok2 = True
                    if _gid_ok2 and _tid_ok2: break
            if not _gid_ok2:
                _errors.append(f"转换后的 GTF 仍缺少 gene_id 属性! 原始 GFF 格式可能不兼容: {_gff_resolved}")
            if not _tid_ok2:
                _errors.append(f"转换后的 GTF 仍缺少 transcript_id 属性! 原始 GFF 格式可能不兼容: {_gff_resolved}")
        except Exception as _e:
            _errors.append(f"GFF → GTF 自动转换失败! 请手动用 gffread 转换: gffread -T {_gff_resolved} -o output.gtf \n错误: {_e}")
    else:
        _errors.append(f"检测到 GFF 格式输入但未找到 gffread! 请:\n"
                        f"  1. 手动转换: gffread -T {_gff_resolved} -o output.gtf\n"
                        f"  2. 或将 ref_gtf 指向已转换的 .gtf 文件")

# 将关键路径更新为解析后的绝对路径 (方便后续 rule 使用)
if _ref_fa: config["paths"]["ref_fa"] = _resolve_path(_ref_fa)
if _ref_gtf: config["paths"]["ref_gtf"] = _resolve_path(_ref_gtf)
if _sample_info: config["paths"]["sample_info"] = _resolve_path(_sample_info)
if _compare_file: config["paths"]["compare_file"] = _resolve_path(_compare_file)

# 打印错误并退出
if _errors:
    print("\n" + "=" * 60)
    print("  配置校验失败! 请修正以下问题:")
    print("=" * 60)
    for i, err in enumerate(_errors, 1):
        print(f"    [{i}] {err}")
    print("=" * 60)
    print("\n  使用 EasyRNASeq2.py 可自动校验并生成配置:")
    print("    python EasyRNASeq2.py --ref_fa genome.fa --ref_gtf genes.gtf --species Athaliana --validate")
    print()
    sys.exit(1)

# ---- 全局变量 ----
RESULT_DIR = _resolve_path(config["paths"]["result_dir"])
TMP_DIR = _resolve_path(config["paths"]["tmp_dir"])
SPECIES = config["project"]["species"]

# ---- 模式检测 ----
_align_tool = config["alignment"]["tool"]
IS_SALMON = (_align_tool == "salmon")
print(f"[INFO] 比对/定量工具: {_align_tool}")

if IS_SALMON:
    print("=" * 60)
    print("  Salmon 模式")
    print("  注意: Salmon 不产生 BAM 文件，以下模块将自动禁用:")
    print("    - 比对质控 (geneBody覆盖度/Read分布/Junction饱和度)")
    print("    - 新转录本发现")
    print("    - 可变剪切 (rMATS)")
    print("    - LncRNA 鉴定")
    print("=" * 60)

    # 自动禁用不兼容模块
    _disabled_modules = []
    if config.get("new_transcript", {}).get("run", False):
        _disabled_modules.append("新转录本发现")
    if config.get("alternative_splicing", {}).get("run", False):
        _disabled_modules.append("可变剪切 (rMATS)")
    if config.get("lncrna", {}).get("run", False):
        _disabled_modules.append("LncRNA 鉴定")
    if config.get("alignment_qc", {}).get("run_coverage", False) or \
       config.get("alignment_qc", {}).get("run_distribution", False):
        _disabled_modules.append("比对质控")

    if _disabled_modules:
        print(f"  [INFO] 已自动禁用: {', '.join(_disabled_modules)}")
        print("=" * 60)

# ---- 通配符定义 ----
samples_df = pd.read_csv(config["paths"]["sample_info"], sep="\t", comment="#")
SAMPLES = samples_df["sample"].tolist()
GROUPS = samples_df["group"].unique().tolist()
HAS_BATCH = "batch" in samples_df.columns
if HAS_BATCH:
    print(f"[INFO] 检测到批次(batch)列: {samples_df['batch'].unique().tolist()}")

# SAMPLES_DATA: sample → {r1_files, r2_files, group}
# 鲁棒的 FASTQ 发现: 支持 _R1/.R1/_1 分隔符, .fq/.fastq/.fq.gz/.fastq.gz 扩展名,
# 以及多 lane 命名 (sample_L001_R1.fq.gz), 同时避免短样本名误匹配 (A6 ≠ A61)
SAMPLES_DATA = {}
import glob as _glob

# 支持的扩展名 (压缩 + 非压缩)
_FASTQ_EXTS = ['.fq.gz', '.fastq.gz', '.fq', '.fastq']

for _, row in samples_df.iterrows():
    sample = row["sample"]
    ngs_dir = _resolve_path(str(row["ngs"]))
    ngs_dir_str = str(ngs_dir)

    r1_files, r2_files = [], []

    for ext in _FASTQ_EXTS:
        # ---- 精确匹配 (无通配符, 零碰撞风险) ----
        for r1_sfx, r2_sfx in [
            (f"{sample}_R1{ext}", f"{sample}_R2{ext}"),   # sample_R1.fq.gz
            (f"{sample}.R1{ext}", f"{sample}.R2{ext}"),   # sample.R1.fq.gz
            (f"{sample}_1{ext}",  f"{sample}_2{ext}"),    # sample_1.fq.gz
        ]:
            r1_path = os.path.join(ngs_dir_str, r1_sfx)
            r2_path = os.path.join(ngs_dir_str, r2_sfx)
            if os.path.isfile(r1_path):
                r1_files.append(r1_path)
            if os.path.isfile(r2_path):
                r2_files.append(r2_path)

        # ---- 多 lane 匹配 (双下划线确保样本边界, A6 不会误匹配 A61) ----
        # E.g.: sample_L001_R1.fq.gz, sample_S1_R1.fq.gz
        for r1_pat, r2_pat in [
            (f"{sample}_*_R1{ext}", f"{sample}_*_R2{ext}"),
            (f"{sample}_*_1{ext}",  f"{sample}_*_2{ext}"),
        ]:
            r1_files.extend(sorted(_glob.glob(os.path.join(ngs_dir_str, r1_pat))))
            r2_files.extend(sorted(_glob.glob(os.path.join(ngs_dir_str, r2_pat))))

    # 去重 + 排序
    r1_files = sorted(set(r1_files))
    r2_files = sorted(set(r2_files))

    if not r1_files:
        print(f"[WARNING] 样本 {sample}: 在 {ngs_dir_str} 中未找到 FASTQ 文件")
        print(f"         支持的命名格式: {{sample}}_R1.fq.gz / {{sample}}.R1.fq.gz / {{sample}}_1.fq.gz")
        print(f"         支持的扩展名:   .fq.gz / .fastq.gz / .fq / .fastq")
        # 保留占位符让流程能继续 DAG 构建 (实际运行时文件校验会报告具体缺失)
        r1_files = [os.path.join(ngs_dir_str, f"{sample}_R1.fq.gz")]
        r2_files = [os.path.join(ngs_dir_str, f"{sample}_R2.fq.gz")]

    SAMPLES_DATA[sample] = {
        "group": row["group"],
        "r1_files": r1_files,
        "r2_files": r2_files,
    }
print(f"[INFO] 共识别 {len(SAMPLES_DATA)} 个样本的 FASTQ 文件")

_diff_config = config.get("diff", {})
if _diff_config.get("batch_effect") and not HAS_BATCH:
    print("[WARNING] diff.batch_effect=true 但样本信息表中没有 batch 列，将忽略批次效应")

# ---- 比较组解析 ----
comparisons_raw = []
with open(config["paths"]["compare_file"]) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        comparisons_raw.append(line)

_comp_list = []
for line in comparisons_raw:
    if '\t' in line:
        parts = line.split('\t')
        ctrl, treat = parts[0].strip(), parts[1].strip()
    elif '_vs_' in line:
        parts = line.split('_vs_')
        ctrl, treat = parts[0].strip(), parts[1].strip()
        print(f"[WARNING] 比较组使用旧格式 '{line}'，建议改用 TSV 格式: '{ctrl}\\t{treat}'")
    else:
        raise ValueError(
            f"比较组文件格式错误，无法解析: '{line}'\n"
            f"支持的格式:\n"
            f"  格式1 (推荐): ctrl_group<TAB>treat_group\n"
            f"  格式2 (兼容):  ctrl_vs_treat"
        )
    comp_name = f"{ctrl}_vs_{treat}"
    _comp_list.append({"ctrl": ctrl, "treat": treat, "name": comp_name})

COMPARISONS = [c["name"] for c in _comp_list]
COMPARISON_MAP = {c["name"]: c for c in _comp_list}

print(f"[INFO] 解析到 {len(COMPARISONS)} 对比较组:")
for c in _comp_list:
    print(f"       {c['ctrl']}  vs  {c['treat']}")

# ---- 根据模式加载规则模块 ----
# 通用模块 (两种模式都需要)
include: "rules/file_check.smk"
include: "rules/data_ready.smk"
include: "rules/qc.smk"

if not IS_SALMON:
    # === 基因组比对模式 (HISAT2/STAR) ===
    include: "rules/align.smk"
    include: "rules/align_qc.smk"
    include: "rules/expression.smk"
    include: "rules/transcript.smk"
    include: "rules/annotation.smk"
    include: "rules/as.smk"
    include: "rules/lncrna.smk"
else:
    # === Salmon 伪比对模式 ===
    include: "rules/salmon.smk"

# 通用下游模块
include: "rules/diff.smk"
include: "rules/enrich.smk"
include: "rules/wgcna.smk"
# report.smk 已移除 — 用户未准备报告模板


# ============================================================
# 主规则: all — 动态构建目标列表
# ============================================================
# 注: 若某模块配置为 run=false，其输出不会加入 all，Snakemake 不执行该模块

def _build_all_targets(salmon_mode=False):
    """根据配置动态构建 all 规则的输出目标列表"""
    T = []  # targets

    # ---- 数据质控 (通用) ----
    T += expand(f"{RESULT_DIR}/02.clean_data/{{sample}}/{{sample}}.R1.fastp.fq.gz", sample=SAMPLES)
    T += expand(f"{RESULT_DIR}/02.clean_data/{{sample}}/{{sample}}.R1.fastp_fastqc.html", sample=SAMPLES)
    T += expand(f"{RESULT_DIR}/02.clean_data/{{sample}}/{{sample}}.R2.fastp_fastqc.html", sample=SAMPLES)
    T.append(f"{RESULT_DIR}/02.clean_data/multiqc_report.html")
    T.append(f"{RESULT_DIR}/01.raw_data/data_stat_all.xls")

    if salmon_mode:
        # ---- Salmon 模式 ----
        T.append(f"{RESULT_DIR}/05.expression/gene/gene_count.xls")
        T.append(f"{RESULT_DIR}/05.expression/gene/gene_tpm.xls")
        T.append(f"{RESULT_DIR}/05.expression/gene/gene_fpkm.xls")
        T.append(f"{RESULT_DIR}/05.expression/gene/gene_norm.xls")
    else:
        # ---- 基因组比对模式 ----
        # 比对
        T += expand(f"{RESULT_DIR}/03.align/{{sample}}/{{sample}}.sorted.bam", sample=SAMPLES)
        T += expand(f"{RESULT_DIR}/03.align/{{sample}}/{{sample}}.sorted.bam.bai", sample=SAMPLES)
        T.append(f"{RESULT_DIR}/03.align/align_stat_all.txt")

        # 比对质控
        _aqc = config.get("alignment_qc", {})
        if _aqc.get("run_coverage", True):
            T.append(f"{RESULT_DIR}/03.align/geneBody_coverage/geneBodyCoverage_all.pdf")
            T.append(f"{RESULT_DIR}/03.align/geneBody_coverage/geneBodyCoverage_all.png")
        if _aqc.get("run_distribution", True):
            T.append(f"{RESULT_DIR}/03.align/read_distribution/read_distribution_all.txt")
            # Per-sample read distribution
            T += expand(f"{RESULT_DIR}/03.align/read_distribution/{{sample}}/{{sample}}.read_distribution.pdf",
                        sample=SAMPLES)
            T.append(f"{RESULT_DIR}/03.align/read_distribution/read_distribution_all.pdf")
        if _aqc.get("run_junction_saturation", False):
            T.append(f"{RESULT_DIR}/03.align/junction_saturation/junction_saturation_all.pdf")
            T.append(f"{RESULT_DIR}/03.align/junction_saturation/junction_saturation_all.png")

        # 表达矩阵 (gene + transcript)
        T.append(f"{RESULT_DIR}/05.expression/gene/gene_count.xls")
        T.append(f"{RESULT_DIR}/05.expression/gene/gene_tpm.xls")
        T.append(f"{RESULT_DIR}/05.expression/gene/gene_fpkm.xls")
        T.append(f"{RESULT_DIR}/05.expression/gene/gene_norm.xls")
        T.append(f"{RESULT_DIR}/05.expression/transcript/transcript_count.xls")
        T.append(f"{RESULT_DIR}/05.expression/transcript/transcript_cpm.xls")
        T.append(f"{RESULT_DIR}/05.expression/transcript/transcript_fpkm.xls")
        T.append(f"{RESULT_DIR}/05.expression/transcript/transcript_norm.xls")

        # 表达可视化
        T.append(f"{RESULT_DIR}/05.expression/gene/box_density/boxplot.pdf")
        T.append(f"{RESULT_DIR}/05.expression/gene/box_density/density.pdf")
        T.append(f"{RESULT_DIR}/05.expression/gene/pca/PCA_2D.pdf")
        T.append(f"{RESULT_DIR}/05.expression/gene/pca/PCA_coordinates.xls")
        T.append(f"{RESULT_DIR}/05.expression/gene/corr/pearson_cor.pdf")
        T.append(f"{RESULT_DIR}/05.expression/gene/corr/spearman_cor.pdf")
        T.append(f"{RESULT_DIR}/05.expression/gene/heatmap/heatmap.pdf")

        # 新转录本发现
        if config.get("new_transcript", {}).get("run", True):
            T.append(f"{RESULT_DIR}/04.transcripts/new_transcripts/new_transcripts.fa")
            T.append(f"{RESULT_DIR}/04.transcripts/new_transcripts/new_transcripts.pep.fa")
            T.append(f"{RESULT_DIR}/04.transcripts/new_transcripts/new_transcripts.gtf")
            T.append(f"{RESULT_DIR}/04.transcripts/all_transcripts.fa")
            T.append(f"{RESULT_DIR}/04.transcripts/all_transcripts.gtf")
            T.append(f"{RESULT_DIR}/04.transcripts/all.gene_trans_map.txt")
            # 转录本统计图表
            T.append(f"{RESULT_DIR}/04.transcripts/new_transcripts/new_transcripts_cds.length.pdf")
            T.append(f"{RESULT_DIR}/04.transcripts/new_transcripts/new_code_stat.pdf")
            T.append(f"{RESULT_DIR}/04.transcripts/all_transcripts.length.pdf")
            T.append(f"{RESULT_DIR}/04.transcripts/Number_stats.pdf")
            T.append(f"{RESULT_DIR}/04.transcripts/trans_for_gene_number_stats.pdf")

        # 功能注释 (仅当有新转录本 + eggNOG 数据库可用时)
        _has_eggnog = bool(config["paths"].get("eggnog_db"))
        if config.get("new_transcript", {}).get("run", True) and _has_eggnog:
            T.append(f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/new_transcripts.all_annotation.xls")
            T.append(f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/new_transcripts.stat.xls")
            T.append(f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/go/go_plot.pdf")
            T.append(f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/kog/kog_plot.pdf")
            T.append(f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/pfam/pfam_plot.pdf")

        # 可变剪切
        if config.get("alternative_splicing", {}).get("run", True):
            T += expand(f"{RESULT_DIR}/07.as/{{comp}}/AS_summary.xls", comp=COMPARISONS)
            if len(COMPARISONS) >= 1:
                T.append(f"{RESULT_DIR}/07.as/AS_event_distribution.pdf")

        # LncRNA 鉴定
        if config.get("lncrna", {}).get("run", True):
            T.append(f"{RESULT_DIR}/08.lnc/novel_lncrna_list.txt")
            T.append(f"{RESULT_DIR}/08.lnc/lncrna_classification.xls")
            T.append(f"{RESULT_DIR}/08.lnc/lncrna_pca.pdf")
            T.append(f"{RESULT_DIR}/08.lnc/lncrna_corr.pdf")
            if config.get("lncrna", {}).get("run_de", False):
                T += expand(f"{RESULT_DIR}/08.lnc/de/{{comp}}/{{comp}}.lnc_DE_results.All.xls", comp=COMPARISONS)
            if config.get("lncrna", {}).get("run_target", False):
                T.append(f"{RESULT_DIR}/08.lnc/target/cis_target_genes.xls")
                T.append(f"{RESULT_DIR}/08.lnc/target/trans_target_genes.xls")

    # ---- 差异分析 (通用) ----
    T += expand(f"{RESULT_DIR}/06.de_enrich/{{comp}}/{{comp}}.DE_results.All.xls", comp=COMPARISONS)
    T += expand(f"{RESULT_DIR}/06.de_enrich/{{comp}}/{{comp}}.DE_results.Significant.xls", comp=COMPARISONS)
    T += expand(f"{RESULT_DIR}/06.de_enrich/{{comp}}/{{comp}}.volcano.pdf", comp=COMPARISONS)
    T += expand(f"{RESULT_DIR}/06.de_enrich/{{comp}}/{{comp}}.heatmap.pdf", comp=COMPARISONS)
    T += expand(f"{RESULT_DIR}/06.de_enrich/{{comp}}/{{comp}}_MA.pdf", comp=COMPARISONS)
    T.append(f"{RESULT_DIR}/06.de_enrich/PCA.pdf")
    T.append(f"{RESULT_DIR}/06.de_enrich/correlation_heatmap.pdf")
    T.append(f"{RESULT_DIR}/06.de_enrich/DE_feature_counts.pdf")
    T.append(f"{RESULT_DIR}/06.de_enrich/DE_feature_counts.xls")
    if len(COMPARISONS) >= 2:
        T.append(f"{RESULT_DIR}/06.de_enrich/venn_upset/venn.pdf")
        T.append(f"{RESULT_DIR}/06.de_enrich/venn_upset/upset.pdf")
        T.append(f"{RESULT_DIR}/06.de_enrich/venn_upset/venn.txt")

    # ---- 富集分析 (通用) ----
    _enr = config.get("enrich", {})
    if _enr.get("run_go", True):
        T += expand(f"{RESULT_DIR}/06.de_enrich/{{comp}}/enrich/go/go_enrich.xls", comp=COMPARISONS)
        T += expand(f"{RESULT_DIR}/06.de_enrich/{{comp}}/enrich/go/go_barplot.pdf", comp=COMPARISONS)
        T += expand(f"{RESULT_DIR}/06.de_enrich/{{comp}}/enrich/go/BP.result.xls", comp=COMPARISONS)
        T += expand(f"{RESULT_DIR}/06.de_enrich/{{comp}}/enrich/go/BP.top10.pdf", comp=COMPARISONS)
    if _enr.get("run_kegg", True):
        T += expand(f"{RESULT_DIR}/06.de_enrich/{{comp}}/enrich/kegg/kegg_enrich.xls", comp=COMPARISONS)
    if _enr.get("run_gsea", True):
        T += expand(f"{RESULT_DIR}/06.de_enrich/{{comp}}/enrich/gsea/gsea_results.xls", comp=COMPARISONS)
    if _enr.get("run_ppi", False):
        T += expand(f"{RESULT_DIR}/06.de_enrich/{{comp}}/enrich/ppi/ppi_network.pdf", comp=COMPARISONS)

    # ---- WGCNA (可选) ----
    if config.get("wgcna", {}).get("run", False):
        T.append(f"{RESULT_DIR}/09.wgcna/module_genes.xls")
        T.append(f"{RESULT_DIR}/09.wgcna/module_trait_cor.xls")
        T.append(f"{RESULT_DIR}/09.wgcna/network_heatmap.pdf")

    # ---- 报告: 已移除 (用户未准备报告模板) ----
    return T


# ---- 主规则 ----
rule all:
    input:
        rules.file_check.output.check_done,
        _build_all_targets(salmon_mode=IS_SALMON),


# ---- 运行完成 / 异常 ----
onsuccess:
    print("=" * 60)
    print("  EasyRNASeq v2 分析完成!")
    print(f"  结果目录: {os.path.abspath(RESULT_DIR)}")
    print(f"  分析模式: {'Salmon 伪比对' if IS_SALMON else _align_tool + ' 基因组比对'}")
    print("=" * 60)
    # 删除可能残留的错误摘要 (全成功则不需要)
    _err_summary = os.path.join(os.path.abspath(RESULT_DIR), "run_error_summary.txt")
    if os.path.exists(_err_summary):
        os.remove(_err_summary)


onerror:
    import traceback
    _err_summary = os.path.join(os.path.abspath(RESULT_DIR), "run_error_summary.txt")
    _total = 0
    _failed = 0
    _errors = []
    try:
        for _j in workflow.jobs:
            _total += 1
            # 检查 job 的 log 是否有错误
            _log = getattr(_j, "log", None)
            _log_path = _log[0] if _log and len(_log) > 0 else None
            _has_err = False
            if _log_path and os.path.exists(str(_log_path)):
                try:
                    with open(str(_log_path), errors='replace') as _lf:
                        _content = _lf.read()
                        if 'Error' in _content or 'error' in _content:
                            _has_err = True
                except:
                    pass
            if _j.failed or _has_err:
                _failed += 1
                _rule = getattr(_j, "rule", "?")
                _wildcards = getattr(_j, "wildcards", {})
                _wc_str = ", ".join(f"{k}={v}" for k,v in _wildcards.items()) if _wildcards else "—"
                _errors.append(f"  [{_rule}]  wildcards: {_wc_str}")
                if _log_path:
                    _errors.append(f"         log: {_log_path}")
    except:
        pass

    with open(_err_summary, 'w') as _es:
        _es.write("=" * 60 + "\\n")
        _es.write("  EasyRNASeq v2 — 运行错误摘要\\n")
        _es.write(f"  时间: {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\\n")
        _es.write(f"  失败/总数: {_failed}/{_total}\\n")
        _es.write("=" * 60 + "\\n\\n")
        if _errors:
            _es.write("失败规则:\\n")
            for _e in _errors:
                _es.write(_e + "\\n")
        else:
            _es.write("(无法自动收集详细错误, 请查看上方 stderr 输出)\\n")
        _es.write("\\n" + "=" * 60 + "\\n")
        _es.write("提示: 查看详细日志请去掉 --quiet 后重新运行\\n")

    print("\\n" + "=" * 60)
    print("  ⚠ 部分任务失败!")
    print(f"  错误摘要已写入: {_err_summary}")
    print(f"  失败/总数: {_failed}/{_total}")
    print("=" * 60)
    if _errors:
        print("\\n失败规则:")
        for _e in _errors:
            print(_e)


# ============================================================
# 样本信息解析 rule (通用)
# ============================================================
rule parse_sample_info:
    """解析样本信息表，生成样本->fq的映射"""
    input:
        sample_info = config["paths"]["sample_info"],
    output:
        sample_map = f"{TMP_DIR}/parsed/sample_fq_map.txt",
    params:
        result_dir = RESULT_DIR,
    conda:
        "envs/qc.yaml"
    script:
        "scripts/prep_sample_info.py"
