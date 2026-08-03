# ============================================================
# Salmon 模块: 转录本伪比对 + 定量
# 当 ALIGN_TOOL=salmon 时, 替代 align.smk + expression.smk
#
# Salmon 跳过基因组比对, 直接将 reads 比对到转录本序列:
#   原始 FASTQ → Salmon 伪比对 → count/TPM 矩阵
#
# 注意: Salmon 模式不产生 BAM 文件, 以下模块将自动禁用:
#   - 比对质控 (align_qc.smk)
#   - 新转录本发现 (transcript.smk, 需 BAM + StringTie)
#   - 可变剪切 (as.smk, 需 BAM)
#   - LncRNA 鉴定 (lncrna.smk, 需新转录本)
# ============================================================

import os

# ---- 转录本 FASTA 来源 ----
def get_transcriptome_fa():
    """获取转录本 FASTA: 优先用户提供, 否则从基因组+GTF自动提取"""
    ref_transcriptome = config["paths"].get("ref_transcriptome")
    if ref_transcriptome and str(ref_transcriptome) and os.path.exists(str(ref_transcriptome)):
        return str(ref_transcriptome)
    # 由 prepare_transcriptome 规则自动生成
    return f"{TMP_DIR}/ref/transcriptome.fa"


# ---- 提取转录本序列 ----
rule prepare_transcriptome:
    """从基因组 FASTA + GTF 提取转录本序列"""
    input:
        fa  = config["paths"]["ref_fa"],
        gtf = config["paths"]["ref_gtf"],
    output:
        transcriptome_fa = f"{TMP_DIR}/ref/transcriptome.fa",
    params:
        outdir = f"{TMP_DIR}/ref",
    conda:
        f"{PIPELINE_DIR}/envs/transcript.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        echo "[INFO] 从基因组+GTF提取转录本序列..."
        gffread -w {output.transcriptome_fa} -g {input.fa} {input.gtf}
        n_trans=$(grep -c '^>' {output.transcriptome_fa} || echo 0)
        echo "[INFO] 提取到 $n_trans 条转录本序列"
        if [ "$n_trans" -eq 0 ]; then
            echo "ERROR: 未能从GTF提取转录本序列, 请检查ref_fa和ref_gtf是否匹配"
            exit 1
        fi
        """


# ---- Salmon 索引构建 ----
rule salmon_index:
    """Salmon 转录本索引"""
    input:
        fa = get_transcriptome_fa(),
    output:
        idx = touch(f"{TMP_DIR}/index/salmon/complete_ref_lens.bin"),
    params:
        outdir = f"{TMP_DIR}/index/salmon",
        kmer   = config.get("salmon", {}).get("kmer_len", 31),
    threads: 8
    conda:
        f"{PIPELINE_DIR}/envs/salmon.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        echo "[INFO] 构建 Salmon 索引 (k={params.kmer})..."
        salmon index \
            -t {input.fa} \
            -i {params.outdir} \
            --kmerLen {params.kmer} \
            -p {threads}
        echo "[INFO] 索引构建完成"
        """


# ---- Salmon 定量 ----
rule salmon_quant:
    """Salmon 转录本水平定量"""
    input:
        r1  = rules.fastp_trim.output.r1_clean,
        r2  = rules.fastp_trim.output.r2_clean,
        idx = rules.salmon_index.output.idx,
    output:
        quant_sf = f"{TMP_DIR}/05.expression/{{sample}}/quant.sf",
    params:
        outdir   = f"{TMP_DIR}/05.expression/{{sample}}",
        libtype  = config.get("salmon", {}).get("libtype", "A"),
        gc_bias  = "--gcBias" if config.get("salmon", {}).get("gc_bias", True) else "",
        seq_bias = "--seqBias" if config.get("salmon", {}).get("seq_bias", True) else "",
        validate = "--validateMappings" if config.get("salmon", {}).get("validate_mappings", True) else "",
    threads: 8
    conda:
        f"{PIPELINE_DIR}/envs/salmon.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        echo "[INFO] Salmon 定量: {wildcards.sample}"
        salmon quant \
            -i {TMP_DIR}/index/salmon \
            -l {params.libtype} \
            -1 {input.r1} \
            -2 {input.r2} \
            -o {params.outdir} \
            -p {threads} \
            {params.gc_bias} \
            {params.seq_bias} \
            {params.validate} \
            2>&1 | tee {params.outdir}/salmon.log

        if [ ! -s {output.quant_sf} ]; then
            echo "ERROR: Salmon 定量失败, quant.sf 为空"
            exit 1
        fi
        echo "[INFO] Salmon 定量完成: $(wc -l < {output.quant_sf}) 条转录本"
        """


# ---- 合并 Salmon 定量为基因级矩阵 ----
# 规则名保持 merge_gene_counts 以兼容下游 diff.smk
rule merge_gene_counts:
    """合并各样本 Salmon 定量结果 → 基因级 count/TPM/FPKM 矩阵"""
    input:
        gtf       = config["paths"]["ref_gtf"],
        quant_sfs = expand(f"{TMP_DIR}/05.expression/{{sample}}/quant.sf", sample=SAMPLES),
    output:
        count_mat = f"{RESULT_DIR}/05.expression/gene/gene_count.xls",
        tpm_mat   = f"{RESULT_DIR}/05.expression/gene/gene_tpm.xls",
        fpkm_mat  = f"{RESULT_DIR}/05.expression/gene/gene_fpkm.xls",
    params:
        outdir      = f"{RESULT_DIR}/05.expression/gene",
        samples_str = lambda w, samples=SAMPLES: ",".join(samples),
    conda:
        f"{PIPELINE_DIR}/envs/salmon.yaml"
    run:
        import os, subprocess, sys
        os.makedirs(f"{params.outdir}", exist_ok=True)
        r_code = """
        import pandas as pd
        import re
        import os
        
        # ---- 1. 从 GTF 构建 tx2gene 映射 ----
        gtf_path = '__SNK_5__'
        tx2gene = {}
        gene_names = {}
        with open(gtf_path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\\t')
                if len(parts) < 9 or parts[2] != 'transcript':
                    continue
                tid_m = re.search(r'transcript_id\\s+"([^"]+)"', parts[8])
                gid_m = re.search(r'gene_id\\s+"([^"]+)"', parts[8])
                gn_m  = re.search(r'gene_name\\s+"([^"]+)"', parts[8])
                if tid_m and gid_m:
                    tx2gene[tid_m.group(1)] = gid_m.group(1)
                    if gn_m:
                        gene_names[gid_m.group(1)] = gn_m.group(1)
        
        print(f'[INFO] tx2gene 映射: {len(tx2gene)} transcripts → {len(set(tx2gene.values()))} genes')
        
        if len(tx2gene) == 0:
            print('ERROR: 无法从GTF解析transcript_id/gene_id, 请检查GTF格式')
            exit(1)
        
        # ---- 2. 读取每个样本的 quant.sf ----
        samples = '__SNK_0__'.split(',')
        gene_counts = {}
        gene_tpms = {}
        
        for sample in samples:
            sf_path = os.path.join('__SNK_6__/05.expression', sample, 'quant.sf')
            if not os.path.exists(sf_path):
                print(f'ERROR: quant.sf not found for {sample}: {sf_path}')
                continue
        
            df = pd.read_csv(sf_path, sep='\\t')
            df['gene_id'] = df['Name'].map(tx2gene)
            unmapped = df['gene_id'].isna().sum()
            if unmapped > 0:
                print(f'  [WARNING] {sample}: {unmapped}/{len(df)} transcripts 未映射到基因')
            df = df.dropna(subset=['gene_id'])
        
            # 基因级汇总: NumReads 求和, TPM 求和
            gene_count = df.groupby('gene_id')['NumReads'].sum()
            gene_tpm   = df.groupby('gene_id')['TPM'].sum()
            gene_counts[sample] = gene_count
            gene_tpms[sample] = gene_tpm
            print(f'  {sample}: {len(df)} transcripts → {len(gene_count)} genes')
        
        # ---- 3. 构建矩阵 ----
        count_df = pd.DataFrame(gene_counts).fillna(0).astype(int)
        tpm_df   = pd.DataFrame(gene_tpms).fillna(0.0)
        
        # FPKM 近似: FPKM ≈ TPM * (total_reads / 1e6) → 简化处理: 从 Count 计算
        total_reads = count_df.sum(axis=0)
        # 获取基因长度
        gene_len = {}
        with open(gtf_path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\\t')
                if len(parts) < 9 or parts[2] != 'exon':
                    continue
                gid_m = re.search(r'gene_id\\s+"([^"]+)"', parts[8])
                if not gid_m:
                    continue
                gid = gid_m.group(1)
                if gid not in gene_len:
                    gene_len[gid] = 0
                start, end = int(parts[3]), int(parts[4])
                gene_len[gid] += end - start + 1
        
        # FPKM = count / (gene_length/1000) / (total_reads/1e6)
        fpkm_data = {}
        for gid in count_df.index:
            length = gene_len.get(gid, 1000)
            fpkm_data[gid] = count_df.loc[gid] / (length / 1000.0) / (total_reads / 1e6)
        fpkm_df = pd.DataFrame(fpkm_data).T.fillna(0.0)
        
        # ---- 4. 保存 ----
        count_df.sort_index().to_csv('__SNK_1__', sep='\\t', float_format='%.0f')
        tpm_df.sort_index().to_csv('__SNK_3__', sep='\\t', float_format='%.2f')
        fpkm_df.sort_index().to_csv('__SNK_2__', sep='\\t', float_format='%.2f')
        
        print(f'\\n=== Salmon 定量汇总完成 ===')
        print(f'  基因数: {count_df.shape[0]}')
        print(f'  样本数: {count_df.shape[1]}')
        print(f'  Count 矩阵: __SNK_1__')
        print(f'  TPM 矩阵:   __SNK_3__')
        print(f'  FPKM 矩阵:  __SNK_2__')
        
        """
        r_code = r_code.replace("__SNK_0__", str(params.samples_str))
        r_code = r_code.replace("__SNK_1__", str(output.count_mat))
        r_code = r_code.replace("__SNK_2__", str(output.fpkm_mat))
        r_code = r_code.replace("__SNK_3__", str(output.tpm_mat))
        r_code = r_code.replace("__SNK_4__", str(params.outdir))
        r_code = r_code.replace("__SNK_5__", str(input.gtf))
        r_code = r_code.replace("__SNK_6__", str(TMP_DIR))
        result = subprocess.run(["python3", "-c", r_code], text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- 标准化表达矩阵 (兼容下游 diff.smk 的依赖) ----
rule normalize_expression:
    """提供标准化表达矩阵 (Salmon 模式直接使用 TPM/FPKM)"""
    input:
        tpm_mat   = rules.merge_gene_counts.output.tpm_mat,
        fpkm_mat  = rules.merge_gene_counts.output.fpkm_mat,
    output:
        norm_mat  = f"{RESULT_DIR}/05.expression/gene/gene_norm.xls",
    params:
        norm = config["expression"]["count_norm"],
    conda:
        f"{PIPELINE_DIR}/envs/salmon.yaml"
    shell:
        """
        if [ "{params.norm}" = "tpm" ]; then
            cp {input.tpm_mat} {output.norm_mat}
        else
            cp {input.fpkm_mat} {output.norm_mat}
        fi
        echo "[INFO] 标准化矩阵 ({params.norm}) → {output.norm_mat}"
        """
