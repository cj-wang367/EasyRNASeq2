#!/usr/bin/env python3
"""
合并各样本的表达定量结果为统一矩阵
"""
import os
import sys
import pandas as pd

def merge_counts(input_dir, output_file, data_type="gene"):
    """合并 featureCounts 或 StringTie 的各样本 count 文件

    支持两种目录结构:
      - 扁平: input_dir/*.gene_count.txt
      - 嵌套: input_dir/{sample}/{sample}.gene_count.txt
    """
    # 递归搜索所有子目录
    files = []
    for root, dirs, filenames in os.walk(input_dir):
        for f in filenames:
            if f.endswith(f".{data_type}_count.txt"):
                files.append(os.path.join(root, f))

    if not files:
        print(f"警告: 在 {input_dir} (含子目录) 中未找到 {data_type} count 文件")
        # 列出目录内容帮助排查
        for item in os.listdir(input_dir):
            item_path = os.path.join(input_dir, item)
            if os.path.isdir(item_path):
                sub_items = os.listdir(item_path)
                print(f"  {item}/: {sub_items[:5]}")
        return

    dfs = []
    for f in sorted(files):
        sample = os.path.basename(f).replace(f".{data_type}_count.txt", "")
        df = pd.read_csv(f, sep="\t", comment="#", header=None)

        # 根据 data_type 和列数判断列名
        ncols = df.shape[1]
        if data_type == "gene":
            df.columns = ["gene_id", "count"]
        elif ncols == 2:
            df.columns = ["transcript_id", "count"]
        elif ncols >= 3:
            # featureCounts 转录本级输出可能有 gene_id + transcript_id + count
            # 取倒数第二列为 ID, 最后一列为 count
            df = df.iloc[:, [-2, -1]]
            df.columns = ["transcript_id", "count"]
        else:
            print(f"WARNING: {f} has {ncols} columns, expected >=2. Skipping.")
            continue

        df = df.set_index(df.columns[0])
        df.columns = [sample]
        dfs.append(df)

    # 合并所有样本
    merged = pd.concat(dfs, axis=1).fillna(0).astype(int)
    merged.index.name = "id"
    merged.to_csv(output_file, sep="\t")
    print(f"合并完成: {output_file} ({len(merged)} features x {len(dfs)} samples)")

if __name__ == "__main__":
    input_dir = sys.argv[1]
    output_file = sys.argv[2]
    data_type = sys.argv[3] if len(sys.argv) > 3 else "gene"
    merge_counts(input_dir, output_file, data_type)
