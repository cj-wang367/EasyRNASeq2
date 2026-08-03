#!/usr/bin/env python3
"""
EasyRNASeq2 — 输入文件预检模块
================================
对参考基因组 FASTA、GTF 注释、样本信息表、比较组文件进行完整性校验。
全部通过后生成 check.done 哨兵文件，否则 exit(1) 并输出详细错误信息。

检查项目:
  1. 基因组 FASTA  — 格式 / 碱基组成 / 染色体名
  2. GTF 注释     — gene_id/transcript_id / 染色体匹配 / 9 列
  3. GFF3 伪装检测 — 扩展名 .gtf 但内部是 GFF3 格式
  4. 样本信息表   — TSV 格式 / 必需列 / FASTQ 存在
  5. 比较组文件   — 组名匹配 / 格式正确
  6. 跨文件互检   — 染色体交叉匹配 / 样本-分组对应
"""

import os
import sys
import re
import gzip
import logging
import argparse
import time
from pathlib import Path

__version__ = "2.0"
AUTHOR = "EasyRNASeq2"


def fmt_time(seconds):
    """格式化耗时"""
    m, s = divmod(int(seconds), 60)
    if m > 0:
        return f"{m}m{s}s"
    return f"{s}s"


# ============================================================
#  基因组 FASTA 检查
# ============================================================
def check_genome_fasta(fa_path):
    """
    检查项:
      - 文件存在且非空
      - 第一行以 '>' 开头
      - 最后一行不以 '>' 开头
      - 碱基组成: ATCGN >= 90%
      - 无非法碱基字符
      - 染色体名不重复
      - 每条染色体都有序列
    返回: (error_count, error_list, chr_lengths)
    """
    errors = []
    chr_lengths = {}

    if not os.path.exists(fa_path):
        return 1, [f"基因组文件不存在: {fa_path}"], chr_lengths
    if os.path.getsize(fa_path) == 0:
        return 1, [f"基因组文件为空: {fa_path}"], chr_lengths

    # 检测 gzip
    opener = gzip.open if fa_path.endswith('.gz') else open

    try:
        with opener(fa_path, 'rt' if fa_path.endswith('.gz') else 'r') as fh:
            lines = fh.readlines()
    except Exception as e:
        return 1, [f"无法读取基因组文件: {fa_path} ({e})"], chr_lengths

    if not lines:
        return 1, [f"基因组文件内容为空: {fa_path}"], chr_lengths

    # 第一行检查
    first_line = lines[0].strip()
    if not first_line.startswith('>'):
        errors.append(f"基因组文件第一行不是 '>' 开头: {first_line[:80]}...")

    # 最后一行检查
    last_line = lines[-1].strip()
    if last_line.startswith('>'):
        errors.append("基因组文件最后一行不能以 '>' 开头（染色体名后无序列）")

    # 碱基统计
    valid_bases = set('ACDEFGHIKLMNPQRSTVWYXUN*-')
    chr_names = []
    current_chr = None
    nom_base = 0  # ATCGN bases
    all_base = 0
    base_errors = 0
    chr_seqs = {}

    for i, line in enumerate(lines, 1):
        l = line.strip()
        if l.startswith('>'):
            chr_name = l[1:].split()[0]
            if chr_name in chr_names:
                errors.append(f"重复的染色体名: {chr_name} (第 {i} 行)")
            chr_names.append(chr_name)
            current_chr = chr_name
            if current_chr not in chr_seqs:
                chr_seqs[current_chr] = 0
        elif current_chr:
            upper_l = l.upper()
            chr_seqs[current_chr] += len(l)
            # 前 500 行做碱基统计
            if i <= 500:
                for ch in upper_l:
                    if ch in 'ATCGN':
                        nom_base += 1
                    all_base += 1
                    if ch not in valid_bases and base_errors < 5:
                        errors.append(f"非法碱基字符 '{ch}' 在第 {i} 行")
                        base_errors += 1

    # 染色体名检查
    if not chr_names:
        errors.append("基因组文件中未检测到任何染色体（无 '>' 行）")

    # 空染色体检查
    for cn, clen in chr_seqs.items():
        chr_lengths[cn] = clen
        if clen == 0:
            errors.append(f"染色体 {cn} 没有序列")

    # 碱基比例检查
    if all_base > 0 and nom_base / all_base < 0.9:
        pct = round(100 * nom_base / all_base, 1)
        errors.append(f"ATCGN 碱基占比仅 {pct}%（< 90%），文件可能不是基因组 FASTA")

    return len(errors), errors, chr_lengths


# ============================================================
#  GTF 注释检查
# ============================================================
def check_gtf(gtf_path, chr_lengths=None):
    """
    检查项:
      - 文件存在且非空
      - 每行 9 列
      - gene_id 和 transcript_id 属性
      - 染色体名与基因组一致
      - 坐标不超过染色体长度
      - exon / transcript 特征存在
      - 检测 GFF3 伪装（ID= 模式）
    返回: (error_count, error_list, gtf_info_dict)
    """
    errors = []
    info = {
        'gene_ids': set(),
        'transcript_ids': set(),
        'has_exon': False,
        'has_transcript': False,
        'is_gff3': False,
        'gff3_lines': 0,
        'total_lines': 0,
    }

    if not os.path.exists(gtf_path):
        return 1, [f"GTF 文件不存在: {gtf_path}"], info
    if os.path.getsize(gtf_path) == 0:
        return 1, [f"GTF 文件为空: {gtf_path}"], info

    # 检测 GFF3 文件头
    with open(gtf_path, 'r') as fh:
        for line in fh:
            if line.startswith('##gff-version'):
                info['is_gff3'] = True
                errors.append(
                    "GTF 文件检测到 '##gff-version' 头，这是 GFF3 格式！\n"
                    f"  请用 gffread 转换: gffread -T {gtf_path} -o output.gtf\n"
                    "  然后重新指定 --ref_gtf output.gtf"
                )
                return len(errors), errors, info
            if not line.startswith('#'):
                break

    # 逐行解析
    with open(gtf_path, 'r') as fh:
        gene_lines_with_transcript_id = 0
        chr_errors = 0
        coord_errors = 0

        for i, line in enumerate(fh, 1):
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.strip().split('\t')
            if len(parts) != 9:
                errors.append(f"第 {i} 行不是 9 列（实际 {len(parts)} 列）")
                if len(errors) >= 20:
                    break
                continue

            info['total_lines'] += 1
            attr = parts[8]
            feature = parts[2]

            # 检测 GFF3 模式 (ID=xxx)
            if re.search(r'\bID=', attr) and not re.search(r'gene_id|transcript_id', attr):
                info['gff3_lines'] += 1
            if info['total_lines'] <= 500 and info['gff3_lines'] / max(info['total_lines'], 1) > 0.8:
                info['is_gff3'] = True
                errors.append(
                    "GTF 文件实际是 GFF3 格式（>80% 记录使用 'ID=' 属性）！\n"
                    f"  请用 gffread -T {gtf_path} -o output.gtf 转换后重新指定"
                )
                return len(errors), errors, info

            # 提取 gene_id / transcript_id
            gid_m = re.search(r'gene_id\s+"?([^\s;"]+)"?', attr)
            tid_m = re.search(r'transcript_id\s+"?([^\s;"]+)"?', attr)

            if gid_m:
                info['gene_ids'].add(gid_m.group(1))
            if tid_m:
                info['transcript_ids'].add(tid_m.group(1))

            # gene 特征不应有 transcript_id
            if feature == 'gene':
                if tid_m:
                    gene_lines_with_transcript_id += 1
            else:
                # 非 gene 行必须有 gene_id 和 transcript_id
                if not gid_m and not tid_m:
                    if info['total_lines'] <= 5:
                        errors.append(
                            f"第 {i} 行缺少 gene_id 和 transcript_id ({feature})"
                        )
                elif not gid_m:
                    if info['total_lines'] <= 5:
                        errors.append(f"第 {i} 行缺少 gene_id ({feature})")
                elif not tid_m:
                    if info['total_lines'] <= 5:
                        errors.append(f"第 {i} 行缺少 transcript_id ({feature})")

            # 特征类型统计
            if feature == 'exon':
                info['has_exon'] = True
            if feature == 'transcript' or feature == 'mRNA':
                info['has_transcript'] = True

            # 染色体名检查
            if chr_lengths and chr_errors < 5:
                if parts[0] not in chr_lengths:
                    errors.append(f"GTF 染色体 {parts[0]} (第 {i} 行) 不在基因组中")
                    chr_errors += 1

            # 坐标检查
            if chr_lengths and parts[0] in chr_lengths and coord_errors < 5:
                try:
                    if int(parts[4]) > chr_lengths[parts[0]]:
                        errors.append(
                            f"GTF 第 {i} 行坐标 ({parts[4]}) 超过染色体"
                            f" {parts[0]} 长度 ({chr_lengths[parts[0]]})"
                        )
                        coord_errors += 1
                except ValueError:
                    pass

            if len(errors) >= 20:
                break

    # 汇总检查
    if not info['has_exon']:
        errors.append("GTF 文件中未检测到 'exon' 特征")
    if not info['has_transcript']:
        errors.append("GTF 文件中未检测到 'transcript' 特征（StringTie 需要）")

    if gene_lines_with_transcript_id > 0:
        logging.warning(
            f"GTF 中有 {gene_lines_with_transcript_id} 行 gene 特征包含"
            f" transcript_id（建议去除）"
        )

    return len(errors), errors, info


# ============================================================
#  样本信息表检查
# ============================================================
def check_sample_info(tsv_path):
    """
    检查项:
      - 文件存在且非空
      - TSV 格式
      - 必需列: sample, group
      - 样本名唯一
      - 每个 group 至少有样本
      - 可选 FASTQ 检查
    返回: (error_count, error_list, samples_dict)
    """
    errors = []
    samples_dict = {'samples': [], 'groups': set(), 'sample_groups': {}}

    if not os.path.exists(tsv_path):
        return 1, [f"样本信息表不存在: {tsv_path}"], samples_dict
    if os.path.getsize(tsv_path) == 0:
        return 1, [f"样本信息表为空: {tsv_path}"], samples_dict

    try:
        with open(tsv_path, 'r') as fh:
            lines = []
            for line in fh:
                line = line.strip()
                if line and not line.startswith('#'):
                    lines.append(line)
    except Exception as e:
        return 1, [f"无法读取样本信息表: {tsv_path} ({e})"], samples_dict

    if not lines:
        return 1, [f"样本信息表无有效数据行: {tsv_path}"], samples_dict

    # 解析 header
    header = lines[0].split('\t')
    if 'sample' not in header:
        errors.append(f"样本信息表缺少 'sample' 列。现有列: {header}")
    if 'group' not in header:
        errors.append(f"样本信息表缺少 'group' 列。现有列: {header}")

    if errors:
        return len(errors), errors, samples_dict

    si = header.index('sample')
    gi = header.index('group')
    has_ngs = 'ngs' in header
    has_fq1 = 'fq1' in header
    ngs_idx = header.index('ngs') if has_ngs else -1
    fq1_idx = header.index('fq1') if has_fq1 else -1
    fq2_idx = header.index('fq2') if 'fq2' in header else -1

    seen_samples = set()
    for i, line in enumerate(lines[1:], 2):
        cols = line.split('\t')
        if len(cols) < max(si, gi) + 1:
            errors.append(f"样本信息表第 {i} 行列数不足 ({len(cols)} 列)")
            continue

        sample = cols[si].strip()
        group = cols[gi].strip()

        if not sample:
            errors.append(f"样本信息表第 {i} 行 sample 为空")
            continue
        if not group:
            errors.append(f"样本信息表第 {i} 行 group 为空 (sample={sample})")
            continue

        if sample in seen_samples:
            errors.append(f"样本名重复: {sample} (第 {i} 行)")
        seen_samples.add(sample)

        samples_dict['samples'].append(sample)
        samples_dict['groups'].add(group)
        samples_dict['sample_groups'][sample] = group

        # FASTQ 检查 (如果提供了 ngs 或 fq1/fq2 列)
        if has_ngs and ngs_idx < len(cols):
            ngs_dir = cols[ngs_idx].strip()
            if ngs_dir and os.path.isdir(ngs_dir):
                r1_files = list(Path(ngs_dir).glob('*R1*.fq.gz')) + \
                           list(Path(ngs_dir).glob('*R1*.fastq.gz')) + \
                           list(Path(ngs_dir).glob('*_1.fq.gz')) + \
                           list(Path(ngs_dir).glob('*_1.fastq.gz'))
                if not r1_files:
                    errors.append(f"样本 {sample}: 目录 {ngs_dir} 中未找到 R1 FASTQ 文件")
        elif has_fq1:
            if fq1_idx < len(cols):
                fq1_path = cols[fq1_idx].strip()
                if fq1_path and not os.path.exists(fq1_path):
                    errors.append(f"样本 {sample}: fq1 文件不存在: {fq1_path}")
            if fq2_idx > 0 and fq2_idx < len(cols):
                fq2_path = cols[fq2_idx].strip()
                if fq2_path and not os.path.exists(fq2_path):
                    errors.append(f"样本 {sample}: fq2 文件不存在: {fq2_path}")

    if not samples_dict['samples']:
        errors.append("样本信息表中无有效样本")

    return len(errors), errors, samples_dict


# ============================================================
#  比较组文件检查
# ============================================================
def check_compare_file(compare_path, valid_groups=None):
    """
    支持两种格式:
      1. TSV 格式: ctrl\ttreat
      2. 旧格式: treat_vs_ctrl (每行一个)
    返回: (error_count, error_list, comparisons)
    """
    errors = []
    comparisons = []

    if not os.path.exists(compare_path):
        return 1, [f"比较组文件不存在: {compare_path}"], comparisons
    if os.path.getsize(compare_path) == 0:
        return 1, [f"比较组文件为空: {compare_path}"], comparisons

    with open(compare_path, 'r') as fh:
        for i, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            if '\t' in line:
                parts = line.split('\t')
                if len(parts) >= 2:
                    comparisons.append((parts[0].strip(), parts[1].strip()))
                else:
                    errors.append(f"比较组文件第 {i} 行格式错误 (TSV 需两列): {line}")
            elif '_vs_' in line:
                parts = line.split('_vs_')
                if len(parts) == 2:
                    comparisons.append((parts[0].strip(), parts[1].strip()))
                else:
                    errors.append(f"比较组文件第 {i} 行格式错误: {line}")
            else:
                errors.append(f"比较组文件第 {i} 行无法解析: {line}")

    if not comparisons:
        errors.append("比较组文件中无有效比较组")
        return len(errors), errors, comparisons

    # 验证组名
    if valid_groups:
        for ctrl, treat in comparisons:
            if ctrl not in valid_groups:
                errors.append(f"比较组 '{ctrl}_vs_{treat}': 对照组 '{ctrl}' 不在样本分组中 ({sorted(valid_groups)})")
            if treat not in valid_groups:
                errors.append(f"比较组 '{ctrl}_vs_{treat}': 实验组 '{treat}' 不在样本分组中 ({sorted(valid_groups)})")
            if ctrl == treat:
                errors.append(f"比较组对照组和实验组相同: '{ctrl}'")

    return len(errors), errors, comparisons


# ============================================================
#  主入口
# ============================================================
def main():
    parser = argparse.ArgumentParser(
        description=f"EasyRNASeq2 输入文件预检 v{__version__}",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--ref_fa', required=True, help='参考基因组 FASTA')
    parser.add_argument('--ref_gtf', required=True, help='参考基因组 GTF 注释')
    parser.add_argument('--sample_info', required=True, help='样本信息表 (TSV)')
    parser.add_argument('--compare_file', default=None, help='比较组文件')
    parser.add_argument('--species', default='unknown', help='物种名')
    parser.add_argument('--outdir', default='.', help='输出目录')
    parser.add_argument('--log_level', default='INFO',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'])
    args = parser.parse_args()

    # 初始化日志
    os.makedirs(args.outdir, exist_ok=True)
    log_file = os.path.join(args.outdir, 'file_check.log')
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format='%(asctime)s [%(levelname)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=[
            logging.FileHandler(log_file, mode='w'),
            logging.StreamHandler(sys.stderr),
        ]
    )

    total_errors = 0
    all_errors = []

    logging.info("=" * 60)
    logging.info(f"  EasyRNASeq2 文件预检 — 物种: {args.species}")
    logging.info("=" * 60)

    # ---- 1. 基因组 FASTA ----
    t0 = time.time()
    logging.info("[1/4] 检查基因组 FASTA ...")
    n_err, errs, chr_lengths = check_genome_fasta(args.ref_fa)
    total_errors += n_err
    all_errors.extend(errs)
    if n_err == 0:
        logging.info(f"  [OK] 基因组 FASTA: {len(chr_lengths)} 条染色体, "
                     f"总长 {sum(chr_lengths.values()):,} bp ({fmt_time(time.time() - t0)})")
    else:
        for e in errs:
            logging.error(f"  {e}")

    # ---- 2. GTF 注释 ----
    t0 = time.time()
    logging.info("[2/4] 检查 GTF 注释 ...")
    n_err, errs, gtf_info = check_gtf(args.ref_gtf, chr_lengths)
    total_errors += n_err
    all_errors.extend(errs)
    if n_err == 0:
        logging.info(f"  [OK] GTF: {len(gtf_info['gene_ids'])} 个基因, "
                     f"{len(gtf_info['transcript_ids'])} 个转录本 "
                     f"({fmt_time(time.time() - t0)})")
    else:
        for e in errs:
            logging.error(f"  {e}")
        if gtf_info.get('is_gff3'):
            logging.error(
                "  ============================================\n"
                "  解决方案: 在服务器上运行以下命令转换 GFF3 → GTF:\n"
                f"  gffread -T {args.ref_gtf} -o output.gtf\n"
                "  然后重新运行 EasyRNASeq2.py --ref_gtf output.gtf ...\n"
                "  ============================================"
            )

    # ---- 3. 样本信息表 ----
    t0 = time.time()
    logging.info("[3/4] 检查样本信息表 ...")
    n_err, errs, samples_info = check_sample_info(args.sample_info)
    total_errors += n_err
    all_errors.extend(errs)
    if n_err == 0:
        logging.info(f"  [OK] 样本信息: {len(samples_info['samples'])} 个样本, "
                     f"{len(samples_info['groups'])} 个分组 "
                     f"({fmt_time(time.time() - t0)})")
    else:
        for e in errs:
            logging.error(f"  {e}")

    # ---- 4. 比较组文件 ----
    if args.compare_file:
        t0 = time.time()
        logging.info("[4/4] 检查比较组文件 ...")
        valid_groups = samples_info.get('groups', set())
        n_err, errs, comparisons = check_compare_file(args.compare_file, valid_groups)
        total_errors += n_err
        all_errors.extend(errs)
        if n_err == 0:
            logging.info(f"  [OK] 比较组: {len(comparisons)} 对 "
                         f"({fmt_time(time.time() - t0)})")
            for ctrl, treat in comparisons:
                logging.info(f"    {ctrl} vs {treat}")
        else:
            for e in errs:
                logging.error(f"  {e}")
    else:
        logging.info("[4/4] 跳过 (未指定 --compare_file)")

    # ---- 输出结果 ----
    sentinel = os.path.join(args.outdir, 'check.done')
    if total_errors == 0:
        logging.info("=" * 60)
        logging.info("  [PASS] 所有输入文件检查通过!")
        logging.info("=" * 60)
        # 写哨兵文件
        with open(sentinel, 'w') as fh:
            fh.write(f"# EasyRNASeq2 文件预检通过\n")
            fh.write(f"# 物种: {args.species}\n")
            fh.write(f"# 基因组: {args.ref_fa}\n")
            fh.write(f"# GTF: {args.ref_gtf}\n")
            fh.write(f"# 染色体数: {len(chr_lengths)}\n")
            fh.write(f"# 基因数: {len(gtf_info.get('gene_ids', set()))}\n")
            fh.write(f"# 转录本数: {len(gtf_info.get('transcript_ids', set()))}\n")
            fh.write(f"# 样本数: {len(samples_info.get('samples', []))}\n")
            fh.write(f"# 分组数: {len(samples_info.get('groups', set()))}\n")
        print(f"\n[OK] 哨兵文件: {sentinel}")
        return 0
    else:
        logging.error("=" * 60)
        logging.error(f"  [FAIL] 共 {total_errors} 个错误!")
        logging.error("=" * 60)
        for i, e in enumerate(all_errors[:50], 1):
            logging.error(f"  [{i}] {e}")
        if len(all_errors) > 50:
            logging.error(f"  ... 还有 {len(all_errors) - 50} 个错误未显示")
        # 确保哨兵文件不存在
        if os.path.exists(sentinel):
            os.remove(sentinel)
        return 1


if __name__ == '__main__':
    sys.exit(main())
