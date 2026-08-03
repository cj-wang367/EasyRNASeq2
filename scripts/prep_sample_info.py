#!/usr/bin/env python3
"""
解析样本信息表，匹配 FASTQ 文件，生成样本->fq 映射表
"""
import os
import sys
import glob
import pandas as pd

def main():
    sample_info = snakemake.input.sample_info
    output_file = snakemake.output.sample_map

    df = pd.read_csv(sample_info, sep="\t", comment="#")

    records = []
    for _, row in df.iterrows():
        sample = row["sample"]
        group = row["group"]
        ngs_dir = row["ngs"]

        # 匹配 R1/R2 FASTQ 文件
        r1_pattern = os.path.join(ngs_dir, "*R1*.f*q.gz")
        r2_pattern = os.path.join(ngs_dir, "*R2*.f*q.gz")

        r1_files = sorted(glob.glob(r1_pattern))
        r2_files = sorted(glob.glob(r2_pattern))

        if not r1_files or not r2_files:
            # 尝试其他模式
            r1_pattern = os.path.join(ngs_dir, "*_1.f*q.gz")
            r2_pattern = os.path.join(ngs_dir, "*_2.f*q.gz")
            r1_files = sorted(glob.glob(r1_pattern))
            r2_files = sorted(glob.glob(r2_pattern))

        if not r1_files or not r2_files:
            print(f"警告: 样本 {sample} 未找到配对 FASTQ 文件, 跳过")
            continue

        # 多个文件用逗号分隔（多 lane 合并）
        records.append({
            "sample": sample,
            "group": group,
            "r1_files": ",".join(r1_files),
            "r2_files": ",".join(r2_files),
            "n_files": len(r1_files),
        })

    result_df = pd.DataFrame(records)
    result_df.to_csv(output_file, sep="\t", index=False)
    print(f"样本信息解析完成: {len(result_df)} 个样本")
    print(f"分组: {result_df['group'].unique().tolist()}")

if __name__ == "__main__":
    main()
