#!/usr/bin/env python3
"""
EasyRNASeq2 — 二代有参转录组分析流程 (基于 Snakemake + Conda)

使用方法:
    # 方式1: 指定 ref_anno_dir (推荐, 自动识别参考文件+索引)
    EasyRNASeq2.py --ref_anno_dir /data/genomes/athaliana --species Athaliana

    # 方式2: 分别指定
    EasyRNASeq2.py --ref_fa genome.fa --ref_gtf genes.gtf --species Athaliana

    # Salmon 快速模式
    EasyRNASeq2.py --ref_anno_dir /data/genomes/athaliana --species Athaliana --align_tool salmon

    # 其他模式
    EasyRNASeq2.py ... --validate     # 仅校验
    EasyRNASeq2.py ... --dryrun       # 预览
    EasyRNASeq2.py --unlock           # 解锁
"""

import argparse, os, sys, subprocess, datetime, glob as _globmod
from pathlib import Path

# ============================================================
# 默认值
# ============================================================
KINGDOMS = ["All", "Plants", "Animals", "Fungi", "Bacteria", "Archaea", "Virus", "Human"]


def get_pipeline_dir():
    return Path(__file__).resolve().parent


def check_file(label, path, required=True):
    if not path:
        msg = f"[ERROR] {label}: 未设置 (必填)" if required else f"[SKIP] {label}: 未设置"
        return (False, msg) if required else (True, msg)
    if os.path.exists(path):
        return True, f"[OK]   {label}: {path}"
    msg = f"[ERROR] {label}: 文件不存在 — {path}" if required else f"[SKIP] {label}: 文件不存在"
    return (not required, msg)


def check_command(cmd):
    try:
        subprocess.run([cmd, "--version"], capture_output=True, timeout=10)
        return True
    except Exception:
        return False


def auto_detect_ref_anno(ref_anno_dir):
    """从 ref_anno_dir 自动识别参考文件"""
    result = {"ref_fa": None, "ref_gtf": None, "ref_pep": None,
              "hisat2_index": None, "star_index": None, "ref_transcriptome": None}
    if not ref_anno_dir or not os.path.isdir(ref_anno_dir):
        return result

    d = str(ref_anno_dir)

    # 找 FASTA
    fa_files = (list(_globmod.glob(os.path.join(d, "*.fa"))) +
                list(_globmod.glob(os.path.join(d, "*.fasta"))))
    for f in fa_files:
        if "transcript" in os.path.basename(f).lower():
            result["ref_transcriptome"] = f
        elif "pep" in os.path.basename(f).lower():
            result["ref_pep"] = f
        elif not result["ref_fa"]:
            result["ref_fa"] = f

    # 找 GTF
    gtf_files = list(_globmod.glob(os.path.join(d, "*.gtf"))) + \
                list(_globmod.glob(os.path.join(d, "*.gtf.gz")))
    if gtf_files:
        result["ref_gtf"] = gtf_files[0]

    # 找 HISAT2 索引
    hisat2_dir = os.path.join(d, "hisat2_index")
    if os.path.isdir(hisat2_dir):
        ht2_files = _globmod.glob(os.path.join(hisat2_dir, "*.1.ht2"))
        if ht2_files:
            result["hisat2_index"] = hisat2_dir

    # 找 STAR 索引
    star_dir = os.path.join(d, "star_index")
    if os.path.isdir(star_dir) and os.path.isfile(os.path.join(star_dir, "Genome")):
        result["star_index"] = star_dir

    return result


def detect_species_from_ref_anno(ref_anno_dir):
    """尝试从目录名推断物种"""
    if ref_anno_dir:
        return os.path.basename(str(ref_anno_dir).rstrip("/\\"))
    return None


# KEGG organism 自动推断表
SPECIES_TO_KEGG = {
    "arabidopsis_thaliana": "ath", "homo_sapiens": "hsa", "mus_musculus": "mmu",
    "rattus_norvegicus": "rno", "danio_rerio": "dre", "drosophila_melanogaster": "dme",
    "saccharomyces_cerevisiae": "sce", "caenorhabditis_elegans": "cel",
    "oryza_sativa": "osa", "zea_mays": "zma", "glycine_max": "gmx",
    "solanum_lycopersicum": "sly", "triticum_aestivum": "tae",
    "brachypodium_distachyon": "bdi", "sorghum_bicolor": "sbi",
}


def infer_kegg_org(species):
    if not species:
        return None
    return SPECIES_TO_KEGG.get(species.lower().replace(" ", "_"))


# ---- 数据库自动查找 ----
DB_SCAN_PATHS = [
    "./data/database",              # 项目本地
    "./database",                   # 简写
    "~/easynaseq/databases",        # 用户目录
]

def auto_detect_databases(pipeline_dir):
    """扫描常见路径，自动发现数据库文件，返回 dict"""
    import glob as _gm
    result = {}

    # 构建搜索路径列表
    scan_roots = [os.path.abspath(p) for p in DB_SCAN_PATHS]
    # 也搜索流程安装目录下的 data/database
    scan_roots.append(os.path.join(str(pipeline_dir), "data", "database"))

    # 展开 ~
    scan_roots = [os.path.expanduser(p) for p in scan_roots]

    found_roots = [p for p in scan_roots if os.path.isdir(p)]

    if not found_roots:
        return result

    print(f"\n  --- 数据库自动查找 ---")
    for root in found_roots:
        print(f"  扫描: {root}")
        # GO: go-basic.obo
        for pattern in ["go-basic.obo", "go*.obo", "gene_ontology*.obo"]:
            hits = _gm.glob(os.path.join(root, "go", pattern))
            if hits and not result.get("go_db"):
                result["go_db"] = hits[0]
                print(f"    → GO: {hits[0]}")
                break
        # rRNA: *.fa / *.fasta
        if not result.get("rrna_db"):
            for pattern in ["*.fa", "*.fasta"]:
                hits = _gm.glob(os.path.join(root, "rrna", pattern))
                if hits:
                    result["rrna_db"] = hits[0]
                    print(f"    → rRNA: {hits[0]}")
                    break
        # STRING: protein.links*.txt
        if not result.get("string_db"):
            for pattern in ["*protein.links*.txt", "*protein.links*.txt.gz", "*.links*.txt"]:
                hits = _gm.glob(os.path.join(root, "string", pattern))
                if hits:
                    result["string_db"] = hits[0]
                    print(f"    → STRING: {hits[0]}")
                    break
        # KEGG
        kegg_dir = os.path.join(root, "kegg")
        if os.path.isdir(kegg_dir) and not result.get("kegg_db"):
            result["kegg_db"] = kegg_dir
            print(f"    → KEGG: {kegg_dir}")
        # NR
        for pattern in ["nr.fa", "nr.fasta", "nr.gz", "nr"]:
            hits = _gm.glob(os.path.join(root, "nr", pattern))
            if hits and not result.get("nr_db"):
                result["nr_db"] = hits[0]
                print(f"    → NR: {hits[0]}")
                break
        # UniProt
        for pattern in ["uniprot_sprot.fasta", "uniprot*.fasta", "uniprot*.fa"]:
            hits = _gm.glob(os.path.join(root, "uniprot", pattern))
            if hits and not result.get("uniprot_db"):
                result["uniprot_db"] = hits[0]
                print(f"    → UniProt: {hits[0]}")
                break
        # Pfam
        for pattern in ["Pfam-A.hmm", "Pfam-A.hmm.gz"]:
            hits = _gm.glob(os.path.join(root, "pfam", pattern))
            if hits and not result.get("pfam_db"):
                result["pfam_db"] = hits[0]
                print(f"    → Pfam: {hits[0]}")
                break
        # eggNOG: 查找含 eggnog.db 的目录
        for pattern in ["eggnog.db"]:
            if not result.get("eggnog_db"):
                # 直接在 root 下找，或在 eggnog/ 子目录找
                for sub in ["", "eggnog"]:
                    hits = _gm.glob(os.path.join(root, sub, pattern))
                    if hits:
                        result["eggnog_db"] = os.path.dirname(hits[0])
                        print(f"    → eggNOG: {os.path.dirname(hits[0])}")
                        break
        # genomes: 查找子目录作为 ref_anno_dir 候选
        genomes_dir = os.path.join(root, "genomes")
        if os.path.isdir(genomes_dir) and not result.get("genomes_dir"):
            subdirs = [os.path.join(genomes_dir, d) for d in os.listdir(genomes_dir)
                       if os.path.isdir(os.path.join(genomes_dir, d))]
            if subdirs:
                result["genomes_dir"] = genomes_dir
                result["ref_anno_dir_candidates"] = subdirs
                print(f"    → Genomes: {len(subdirs)} 个物种目录")

    if not result:
        print(f"    (未找到数据库目录, 将跳过相关分析)")

    print()
    return result


def generate_config(args, pipeline_dir, ref_info):
    """生成 config.generated.yaml"""
    config_path = pipeline_dir / "config" / "config.generated.yaml"

    is_salmon = (args.align_tool == "salmon")

    # 自动推断缺失值
    species = args.species or ref_info.get("species_guess") or "Unknown"
    ref_fa = args.ref_fa or ref_info.get("ref_fa", "null")
    ref_gtf = args.ref_gtf or ref_info.get("ref_gtf", "null")
    ref_pep = args.ref_pep or ref_info.get("ref_pep", "null")
    ref_trans = args.ref_transcriptome or ref_info.get("ref_transcriptome", "null")
    hisat2_idx = args.hisat2_index or ref_info.get("hisat2_index", "null")
    star_idx = args.star_index or ref_info.get("star_index", "null")
    kegg_org = args.kegg_org or infer_kegg_org(species) or "null"

    def _fmt(v):
        if v is None or v == "" or v == "null":
            return "null"
        return f'"{v}"'

    # 模块开关
    run_new_transcript = "false" if is_salmon else "true"
    run_as = "false" if is_salmon else "true"
    run_lncrna = "false" if is_salmon else "true"
    run_align_qc = "false" if is_salmon else "true"
    exp_tool = "salmon" if is_salmon else args.exp_tool

    # 确定 index 设置
    if is_salmon:
        index_val = "null"
    elif args.align_tool == "hisat2" and hisat2_idx != "null":
        index_val = hisat2_idx
    elif args.align_tool == "star" and star_idx != "null":
        index_val = star_idx
    else:
        index_val = "null"

    content = f"""# 由 EasyRNASeq2.py 自动生成
# 时间: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

project:
  id: "{species}"
  title: "{species} RNA-seq"
  species: "{species}"
  kingdom: "{args.kingdom}"
  data_type: "Illumina"

paths:
  sample_info: "{args.sample_info}"
  compare_file: "{args.compare_file}"
  ref_anno_dir: {_fmt(args.ref_anno_dir)}
  ref_fa: {_fmt(ref_fa)}
  ref_gtf: {_fmt(ref_gtf)}
  ref_gff: null
  ref_pep: {_fmt(ref_pep)}
  ref_transcriptome: {_fmt(ref_trans)}
  result_dir: "{args.outdir}"
  tmp_dir: "{args.tmpdir}"
  kegg_db: {_fmt(args.kegg_db)}
  go_db: {_fmt(args.go_db)}
  nr_db: {_fmt(args.nr_db)}
  uniprot_db: {_fmt(args.uniprot_db)}
  pfam_db: {_fmt(args.pfam_db)}
  string_db: {_fmt(args.string_db)}
  rrna_db: {_fmt(args.rrna_db)}
  eggnog_db: {_fmt(args.eggnog_db)}
  hisat2_path: ""
  star_path: ""

data_preprocess:
  run: true
  bases_data: {args.bases_data}
  cut_sample_reads: 50000
  run_rmrRNA: {"true" if args.rrna_db else "false"}
  max_rrna_rate: 5
  run_extract: true
  target_mapped_ratio: 90
  min_mapped_ratio: 85
  minmin_mapped_ratio: 50
  float_percent: 1
  lower_bound_ratio: 0.95
  upper_bound_ratio: 1.05

qc:
  phred: 33
  fastp_min_len: 30
  adapter_r1: null
  adapter_r2: null
  trim_poly_g: true
  poly_g_min_len: 10

alignment:
  tool: "{args.align_tool}"
  index: {index_val}
  hisat2_index_onlyFa: {str(args.hisat2_index_onlyfa).lower()}
  strand: "{args.strand}"
  bam_index_tool: "sambamba"
  hisat2_cpu: 6
  hisat2_mem_gb: 16
  hisat2_index_mem_gb: 50
  star_cpu: 6
  star_mem_gb: 32

alignment_qc:
  run_coverage: {run_align_qc}
  run_distribution: {run_align_qc}
  run_junction_saturation: {run_align_qc}
  run_infer_experiment: false

new_transcript:
  run: {run_new_transcript}
  assemble_tool: "stringtie"
  merge_tool: "stringtie"
  orf_tool: "TransDecoder"
  circl_window: 500000

expression:
  tool: "{exp_tool}"
  count_norm: "{args.count_norm}"
  run_exon: false

diff:
  fc_calculate_type: 2
  diff_pvalue: {args.diff_pvalue}
  logFC: {args.logFC}
  pOrFDR: "{args.diff_ptype}"
  method: "auto"
  dispersion: 0.05
  batch_effect: {str(args.batch_effect).lower()}

enrich:
  run_go: true
  run_kegg: true
  run_gsea: true
  run_gsva: false
  run_reactome: false
  run_ppi: {str(args.run_ppi).lower()}
  enrich_plot_x: "rich_factor"
  enrich_plot_p: "pvalue"
  kegg_org: {_fmt(kegg_org)}
  ppi_score_threshold: 400

alternative_splicing:
  run: {run_as}
  read_length: 150

lncrna:
  run: {run_lncrna}
  run_de: {run_lncrna}
  run_target: {run_lncrna}
  run_enrich: false
  cis_window: 100000

wgcna:
  run: {str(args.run_wgcna).lower()}
  power: 9
  min_module_size: 100
  merge_cut_height: 0.25
  line: 0.85

salmon:
  kmer_len: 31
  libtype: "A"
  gc_bias: true
  seq_bias: true
  validate_mappings: true

variation:
  run: false
  maf: 0.05
  max_missing: 0.75

report:
  run: false
  language: "cn"
  company: "cn"

runtime:
  max_jobs: {args.jobs}
  latancy_wait: 60
  use_conda: true
  cluster: null
"""
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(content, encoding="utf-8")
    return config_path


def validate_inputs(args, ref_info):
    errors, warnings = [], []

    # ref_fa / ref_gtf 可从 ref_anno_dir 自动获取
    ref_fa = args.ref_fa or ref_info.get("ref_fa")
    ref_gtf = args.ref_gtf or ref_info.get("ref_gtf")

    if not ref_fa:
        errors.append("参考基因组 FASTA 未设置 (--ref_fa 或 --ref_anno_dir)")
    elif not os.path.exists(str(ref_fa)):
        errors.append(f"参考基因组 FASTA 不存在: {ref_fa}")

    if not ref_gtf:
        errors.append("参考基因组 GTF 未设置 (--ref_gtf 或 --ref_anno_dir)")
    elif not os.path.exists(str(ref_gtf)):
        errors.append(f"参考基因组 GTF 不存在: {ref_gtf}")

    if not args.species:
        species_guess = ref_info.get("species_guess")
        if species_guess:
            args.species = species_guess
            warnings.append(f"未指定 --species, 已从目录名自动推断: {species_guess}")
        else:
            errors.append("物种名未设置 (--species)")

    for label, val, required in [
        ("样本信息表", args.sample_info, True),
        ("比较组文件", args.compare_file, True),
        ("转录本序列", args.ref_transcriptome, False),
        ("rRNA 数据库", args.rrna_db, False),
        ("KEGG 数据库", args.kegg_db, False),
        ("GO 数据库", args.go_db, False),
        ("STRING PPI 数据库", args.string_db, False),
    ]:
        ok, msg = check_file(label, val, required)
        if not ok:
            errors.append(msg)
        else:
            print(f"  {msg}")

    # ref_anno_dir 信息
    if args.ref_anno_dir:
        print(f"  [REF_ANNO_DIR] {args.ref_anno_dir}")
        if ref_info.get("ref_fa"):
            print(f"     → ref_fa: {ref_info['ref_fa']}")
        if ref_info.get("ref_gtf"):
            print(f"     → ref_gtf: {ref_info['ref_gtf']}")
        if ref_info.get("ref_pep"):
            print(f"     → ref_pep: {ref_info['ref_pep']}")
        if ref_info.get("hisat2_index"):
            print(f"     → hisat2 索引: {ref_info['hisat2_index']}")
        if ref_info.get("star_index"):
            print(f"     → star 索引: {ref_info['star_index']}")

    # 样本信息表
    if args.sample_info and os.path.exists(args.sample_info):
        try:
            import pandas as pd
            df = pd.read_csv(args.sample_info, sep="\t", comment="#")
            missing = [c for c in ["group", "sample", "ngs"] if c not in df.columns]
            if missing:
                errors.append(f"样本信息表缺少列: {missing}")
            else:
                print(f"  [OK]   样本信息表: {len(df)} 个样本, {df['group'].nunique()} 个分组")
                if "batch" in df.columns:
                    print(f"  [INFO] 检测到 batch 列 — 可启用 --batch_effect")
        except Exception as e:
            errors.append(f"无法解析样本信息表: {e}")

    # 比较组文件
    if args.compare_file and os.path.exists(args.compare_file):
        n_comp = 0
        with open(args.compare_file) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    n_comp += 1
        print(f"  [OK]   比较组文件: {n_comp} 对比较 (格式: treat_vs_ctrl)")

    # conda / snakemake
    if check_command("snakemake"):
        print(f"  [OK]   Snakemake 可用")
    else:
        warnings.append("Snakemake 未安装 (请先 conda activate snakemake)")
    if check_command("conda"):
        print(f"  [OK]   Conda 可用")

    # 线程检查
    if args.jobs > 16:
        warnings.append(f"并发任务数 {args.jobs} 偏高，64线程机器建议 ≤8，避免资源耗尽")

    # Salmon 模式
    if args.align_tool == "salmon":
        print(f"\n  --- Salmon 模式 ---")
        print(f"  [INFO] 自动禁用: 比对质控 / 新转录本 / 可变剪切 / LncRNA")
        print(f"  [INFO] 保留: 差异分析 / GO-KEGG富集 / GSEA / WGCNA")

    # PPI 模式
    if args.run_ppi and not args.string_db:
        warnings.append("--run_ppi 已启用但未指定 --string_db, PPI 将跳过")

    if warnings:
        print(f"\n  --- 警告 ({len(warnings)}) ---")
        for w in warnings:
            print(f"  [WARN] {w}")

    return len(errors) == 0, errors


def run_snakemake(args, pipeline_dir, config_path):
    snakefile = pipeline_dir / "Snakefile"
    # 确保输出目录存在
    cwd = Path.cwd()
    outdir_abs = cwd / args.outdir
    outdir_abs.mkdir(parents=True, exist_ok=True)

    # 日志文件 (每次运行生成新的, 放在工作目录方便查看)
    import time
    _ts = time.strftime("%Y%m%d_%H%M%S")
    log_file = cwd / f"EasyRNASeq_run_{_ts}.log"

    cmd = [
        "snakemake", "--snakefile", str(snakefile),
        "--configfile", str(config_path),
        "--cores", str(args.jobs), "--use-conda",
        "--latency-wait", "60", "--keep-going",
    ]
    # 默认安静模式 (不打印 shell 命令), --verbose 可恢复详细输出
    if not args.verbose:
        cmd += ["--quiet", "progress"]
    if args.forceall:
        cmd.append("--forceall")
    # Conda 环境始终复用流程目录中的缓存 (多用户共享, 避免重复下载)
    _conda_prefix = str(args.conda_prefix) if args.conda_prefix else str(pipeline_dir / ".snakemake" / "conda")
    cmd += ["--conda-prefix", _conda_prefix]
    if args.dryrun:
        cmd.append("--dry-run")
    elif args.unlock:
        cmd.append("--unlock")
    if args.verbose:
        cmd.append("--verbose")

    print(f"\n{'=' * 60}")
    print(f"  Snakemake 命令:")
    print(f"    {' '.join(cmd)}")
    print(f"  运行日志: {log_file}")
    print(f"  缓存位置: .snakemake/ → {cwd / '.snakemake'}")
    print(f"  Conda 环境: {_conda_prefix}")
    print(f"{'=' * 60}\n")

    # 同时输出到终端和日志文件
    import subprocess as _sp
    with open(str(log_file), "w", encoding="utf-8") as _log_fh:
        _log_fh.write(f"# EasyRNASeq2 Run Log\n")
        _log_fh.write(f"# Time: {_ts}\n")
        _log_fh.write(f"# CMD: {' '.join(cmd)}\n")
        _log_fh.write(f"# CWD: {cwd}\n")
        _log_fh.write(f"# {'=' * 60}\n\n")
        _log_fh.flush()

        # 注意: 不设置 cwd, Snakemake 在用户工作目录运行
        # .snakemake/ 缓存放在用户工作目录, 流程目录保持干净
        _proc = _sp.Popen(
            cmd,
            stdout=_sp.PIPE,
            stderr=_sp.STDOUT,
            text=True,
            bufsize=1,
        )
        for _line in _proc.stdout:
            print(_line, end="")         # 终端输出
            _log_fh.write(_line)         # 日志文件
            _log_fh.flush()
        _proc.wait()

    rc = _proc.returncode

    # 生成运行摘要
    print(f"\n{'=' * 60}")
    if rc == 0:
        print(f"  ✅ 全部任务成功完成!")
    else:
        print(f"  ⚠️  部分任务失败 (退出码: {rc})")
        print(f"  查看错误摘要: {args.outdir}/run_error_summary.txt")
        print(f"  查看 HTML 报告: {args.outdir}/run_report.html")
    print(f"  完整日志: {log_file}")
    print(f"{'=' * 60}")

    # 复制 Snakemake 内部日志到工作目录
    # .snakemake/ 现在位于用户工作目录
    _snakemake_log = cwd / ".snakemake" / "log"
    if _snakemake_log.exists():
        _smk_log_dest = cwd / f"EasyRNASeq_snakemake_{_ts}.log"
        try:
            _latest_log = sorted(_snakemake_log.glob("*.log"),
                                key=lambda p: p.stat().st_mtime, reverse=True)
            if _latest_log:
                import shutil
                shutil.copy(_latest_log[0], _smk_log_dest)
                print(f"\n[INFO] Snakemake 内部日志: {_smk_log_dest}")
        except Exception:
            pass

    print(f"\n[INFO] 完整运行日志: {log_file}")
    return rc


def main():
    parser = argparse.ArgumentParser(
        description="EasyRNASeq2 — 二代有参转录组分析流程",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  EasyRNASeq2.py --ref_anno_dir /data/genomes/athaliana --species Athaliana
  EasyRNASeq2.py --ref_fa genome.fa --ref_gtf genes.gtf --species Athaliana --align_tool salmon
  EasyRNASeq2.py --ref_anno_dir /data/genomes/athaliana --species Athaliana --validate

自定义 run.sh:
  #!/bin/bash
  conda activate snakemake
  EasyRNASeq2.py --ref_anno_dir /data/genomes/athaliana --species Athaliana --jobs 8
        """)

    # ---- 必填参数 ----
    required = parser.add_argument_group("参考基因组 (任选其一)")
    grp = required.add_mutually_exclusive_group()
    grp.add_argument("--ref_anno_dir", help="参考基因组注释目录 (含 FASTA/GTF/索引)")
    grp_ref = parser.add_argument_group("或分别指定")
    grp_ref.add_argument("--ref_fa", help="参考基因组 FASTA")
    grp_ref.add_argument("--ref_gtf", help="参考基因组 GTF")
    required.add_argument("--species", help="物种名 (如 Arabidopsis_thaliana/Homo_sapiens)")

    # ---- 输入文件 ----
    io = parser.add_argument_group("输入/输出")
    io.add_argument("--sample_info", default="sample_info.tsv", help="样本信息表 (默认: sample_info.tsv)")
    io.add_argument("--compare_file", default="compare.tsv", help="比较组文件 (默认: compare.tsv, 格式: treat_vs_ctrl)")
    io.add_argument("--ref_transcriptome", help="转录本 FASTA (Salmon 用, 留空自动提取)")
    io.add_argument("--ref_pep", help="蛋白序列 FASTA (注释用)")
    io.add_argument("--hisat2_index", help="已有 HISAT2 索引路径")
    io.add_argument("--star_index", help="已有 STAR 索引路径")
    io.add_argument("--kingdom", default="Plants", choices=KINGDOMS, help="界分类 (默认: Plants)")
    io.add_argument("--outdir", default="results", help="输出目录 (默认: results)")
    io.add_argument("--tmpdir", default="tmp", help="临时目录 (默认: tmp)")
    io.add_argument("--jobs", type=int, default=8, help="最大并发数 (默认: 8, 64线程机器建议≤10)")
    io.add_argument("--conda_prefix", help="Conda 环境共享目录 (复用预构建环境, 默认: .snakemake/conda)")

    # ---- 分析工具 ----
    tool = parser.add_argument_group("分析工具")
    tool.add_argument("--align_tool", default="hisat2", choices=["hisat2","star","salmon"],
                      help="比对/定量工具 (默认: hisat2)")
    tool.add_argument("--exp_tool", default="featureCounts", choices=["featureCounts","stringtie","rsem"],
                      help="定量工具 (hisat2/star 模式, 默认: featureCounts)")
    tool.add_argument("--strand", default="unstranded", choices=["unstranded","firststrand","secondstrand"])
    tool.add_argument("--count_norm", default="fpkm", choices=["fpkm","tpm"])
    tool.add_argument("--bases_data", type=int, default=6, help="合同数据量 G (默认: 6)")
    tool.add_argument("--hisat2_index_onlyfa", action="store_true", help="大基因组仅FASTA建索引")

    # ---- 差异分析 ----
    diff = parser.add_argument_group("差异分析")
    diff.add_argument("--diff_pvalue", type=float, default=0.05)
    diff.add_argument("--logFC", type=float, default=1.0)
    diff.add_argument("--diff_ptype", default="padj", choices=["padj","pvalue"])
    diff.add_argument("--batch_effect", action="store_true", help="批次效应校正")

    # ---- 富集 + PPI + WGCNA ----
    enrich_grp = parser.add_argument_group("富集分析/PPI/WGCNA")
    enrich_grp.add_argument("--kegg_org", help="KEGG 物种代码 (如 ath/hsa/mmu, 留空自动推断)")
    enrich_grp.add_argument("--run_ppi", action="store_true", help="启用 STRING PPI 蛋白互作网络")
    enrich_grp.add_argument("--run_wgcna", action="store_true", help="启用 WGCNA 共表达网络")

    # ---- 数据库 ----
    db = parser.add_argument_group("数据库路径")
    db.add_argument("--kegg_db", help="KEGG 本地数据库")
    db.add_argument("--go_db", help="GO 数据库 (go-basic.obo)")
    db.add_argument("--nr_db", help="NR 蛋白数据库")
    db.add_argument("--uniprot_db", help="UniProt/Swiss-Prot")
    db.add_argument("--pfam_db", help="Pfam-A.hmm")
    db.add_argument("--string_db", help="STRING PPI 数据库 (如 3702.protein.links.v12.0.txt)")
    db.add_argument("--rrna_db", help="rRNA 参考序列 (Silva)")
    db.add_argument("--eggnog_db", help="eggNOG-mapper 数据库目录 (含 eggnog.db + eggnog_proteins.dmnd)")

    # ---- 运行模式 ----
    mode = parser.add_argument_group("运行模式")
    mode.add_argument("--dryrun", action="store_true", help="预览任务")
    mode.add_argument("--unlock", action="store_true", help="解锁工作目录")
    mode.add_argument("--validate", action="store_true", help="仅校验")
    mode.add_argument("--forceall", action="store_true",
                      help="强制重跑所有任务 (忽略已有结果)")
    mode.add_argument("--verbose", action="store_true",
                      help="详细输出 (打印每个任务的完整 shell 命令, 调试用)")

    args = parser.parse_args()
    pipeline_dir = get_pipeline_dir()

    # ---- 数据库自动查找 (用户未手动指定时自动扫描) ----
    db_auto = auto_detect_databases(pipeline_dir)
    if not args.go_db: args.go_db = db_auto.get("go_db")
    if not args.rrna_db: args.rrna_db = db_auto.get("rrna_db")
    if not args.string_db: args.string_db = db_auto.get("string_db")
    if not args.kegg_db: args.kegg_db = db_auto.get("kegg_db")
    if not args.nr_db: args.nr_db = db_auto.get("nr_db")
    if not args.uniprot_db: args.uniprot_db = db_auto.get("uniprot_db")
    if not args.pfam_db: args.pfam_db = db_auto.get("pfam_db")
    if not args.eggnog_db: args.eggnog_db = db_auto.get("eggnog_db")

    cwd = Path.cwd()

    # 相对→绝对路径 (所有路径参数统一基于 CWD 解析)
    _path_attrs = [
        "ref_fa", "ref_gtf", "ref_pep", "ref_transcriptome", "ref_anno_dir",
        "sample_info", "compare_file", "hisat2_index", "star_index",
        "kegg_db", "go_db", "nr_db", "uniprot_db", "pfam_db", "string_db", "rrna_db", "eggnog_db",
        "outdir", "tmpdir", "conda_prefix",
    ]
    for attr in _path_attrs:
        val = getattr(args, attr, None)
        if val and not os.path.isabs(str(val)):
            setattr(args, attr, str(cwd / val))

    # ---- ref_anno_dir 自动识别 ----
    ref_info = {}
    if args.ref_anno_dir:
        ref_info = auto_detect_ref_anno(args.ref_anno_dir)
        if not ref_info.get("species_guess"):
            ref_info["species_guess"] = detect_species_from_ref_anno(args.ref_anno_dir)

    # 自动推断物种
    if not args.species and ref_info.get("species_guess"):
        args.species = ref_info["species_guess"]

    print(f"\n{'=' * 60}")
    print(f"  EasyRNASeq2 — 二代有参转录组分析流程")
    print(f"{'=' * 60}")
    print(f"  流程目录: {pipeline_dir}")
    print(f"  工作目录: {cwd}")
    print(f"  物种: {args.species or '(待指定)'}")
    print(f"  工具: {args.align_tool}")
    print(f"  并发数: {args.jobs}")
    print()

    # ---- unlock 模式: 跳过所有校验, 直接解锁 ----
    if args.unlock:
        print(f"\n  --- 解锁工作目录 ---")
        snakefile = pipeline_dir / "Snakefile"
        cmd = ["snakemake", "--snakefile", str(snakefile), "--unlock"]
        import subprocess as _sp2
        _sp2.run(cmd)  # 不设 cwd, 使用用户工作目录中的 .snakemake/
        print(f"  [OK] 解锁完成")
        return

    # ---- 校验 ----
    print(f"  --- 输入校验 ---")
    all_ok, errors = validate_inputs(args, ref_info)

    if not all_ok:
        print(f"\n{'=' * 60}")
        print(f"  校验失败: {len(errors)} 个错误")
        print(f"{'=' * 60}")
        for i, err in enumerate(errors, 1):
            print(f"  [{i}] {err}")
        print(f"\n  提示: python EasyRNASeq2.py --help 查看完整参数")
        print(f"{'=' * 60}\n")
        sys.exit(1)

    print(f"\n  [OK] 校验通过!")

    if args.validate:
        print(f"\n{'=' * 60}")
        print(f"  所有输入就绪!")
        print(f"{'=' * 60}\n")
        return

    # ---- 生成配置 + 运行 ----
    print(f"\n  --- 生成配置 ---")
    config_path = generate_config(args, pipeline_dir, ref_info)
    print(f"  [OK] 配置: {config_path}")

    rc = run_snakemake(args, pipeline_dir, config_path)
    print()
    if rc == 0:
        print(f"{'=' * 60}")
        print(f"  分析完成! 结果: {cwd / args.outdir}")
        print(f"{'=' * 60}")
    else:
        print(f"{'=' * 60}")
        print(f"  异常退出 (code: {rc})")
        print(f"  排查: 1) conda info --envs  2) EasyRNASeq2.py --unlock")
        print(f"{'=' * 60}")
    sys.exit(rc)


if __name__ == "__main__":
    main()
