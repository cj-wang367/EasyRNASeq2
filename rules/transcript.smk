# ============================================================
# 新转录本发现模块: StringTie 组装 + GffCompare 合并 + TransDecoder ORF
# ============================================================

# ---- StringTie 转录本组装 ----
rule stringtie_assemble:
    """StringTie 各样本独立组装（非 -e 模式）"""
    input:
        bam = f"{RESULT_DIR}/03.align/{{sample}}/{{sample}}.sorted.bam",
        gtf = config["paths"]["ref_gtf"],
    output:
        gtf_out = f"{TMP_DIR}/04.transcripts/{{sample}}/{{sample}}.stringtie.assemble.gtf",
    params:
        outdir  = f"{TMP_DIR}/04.transcripts/{{sample}}",
        strand = {
            "unstranded": "",
            "firststrand": "--rf",
            "secondstrand": "--fr",
        }[config["alignment"]["strand"]],
    threads: 4
    conda:
        f"{PIPELINE_DIR}/envs/transcript.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        stringtie {input.bam} \\
            -G {input.gtf} \\
            -p {threads} \\
            {params.strand} \\
            -o {output.gtf_out}
        """


# ---- StringTie merge ----
rule stringtie_merge:
    """合并所有样本的 StringTie GTF"""
    input:
        gtfs = expand(f"{TMP_DIR}/04.transcripts/{{sample}}/{{sample}}.stringtie.assemble.gtf", sample=SAMPLES),
        ref_gtf = config["paths"]["ref_gtf"],
    output:
        merged_gtf = f"{TMP_DIR}/04.transcripts/merged.stringtie.gtf",
    params:
        outdir = f"{TMP_DIR}/04.transcripts",
    threads: 8
    conda:
        f"{PIPELINE_DIR}/envs/transcript.yaml"
    shell:
        """
        # 创建 merge list
        printf '%s\\n' {input.gtfs} > {params.outdir}/mergelist.txt
        stringtie --merge \\
            -G {input.ref_gtf} \\
            -p {threads} \\
            -o {output.merged_gtf} \\
            {params.outdir}/mergelist.txt
        """


# ---- GffCompare 与参考比较 ----
rule gffcompare:
    """GffCompare 比较合并转录本与参考注释"""
    input:
        merged_gtf = rules.stringtie_merge.output.merged_gtf,
        ref_gtf    = config["paths"]["ref_gtf"],
    output:
        compared   = f"{TMP_DIR}/04.transcripts/gffcmp.annotated.gtf",
        tracking   = f"{TMP_DIR}/04.transcripts/gffcmp.tracking",
        stats      = f"{TMP_DIR}/04.transcripts/gffcmp.stats",
    params:
        prefix = f"{TMP_DIR}/04.transcripts/gffcmp",
    conda:
        f"{PIPELINE_DIR}/envs/transcript.yaml"
    shell:
        """
        gffcompare -r {input.ref_gtf} \\
            -o {params.prefix} \\
            {input.merged_gtf}
        """


# ---- 提取新转录本 FASTA ----
rule extract_transcript_fasta:
    """从新转录本 GTF 提取转录本序列"""
    input:
        gtf = rules.gffcompare.output.compared,
        fa  = config["paths"]["ref_fa"],
    output:
        new_fa   = f"{TMP_DIR}/04.transcripts/new_transcripts.fa",
        all_fa   = f"{TMP_DIR}/04.transcripts/all_transcripts.fa",
    conda:
        f"{PIPELINE_DIR}/envs/transcript.yaml"
    run:
        import os, subprocess, sys
        # [TODO] 原 shell 命令: gffread -w {output.all_fa} -g {input.fa} {input.gtf} 2>&1 || {{
        # [TODO] 原 shell 命令: echo "[WARN] gffread failed — creating empty FASTA placeholders"
        # [TODO] 原 shell 命令: echo "" > {output.all_fa}
        # [TODO] 原 shell 命令: echo "" > {output.new_fa}
        # [TODO] 原 shell 命令: exit 0
        # [TODO] 原 shell 命令: }}
        # [TODO] 原 shell 命令: if [ ! -s {output.all_fa} ] || [ $(grep -c '^>' {output.all_fa} 2>/dev/null || echo 0) -eq 0 ]; then
        # [TODO] 原 shell 命令: echo "[WARN] No transcripts in FASTA — creating empty placeholder"
        # [TODO] 原 shell 命令: echo "" > {output.all_fa}
        # [TODO] 原 shell 命令: echo "" > {output.new_fa}
        # [TODO] 原 shell 命令: exit 0
        # [TODO] 原 shell 命令: fi
        r_code = """
        import re
        new_ids = set()
        with open('__SNK_2__') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                m = re.search(r'class_code\\s+"([uijx])"', line)
                if m:
                    m2 = re.search(r'transcript_id\\s+"([^"]+)"', line)
                    if m2:
                        new_ids.add(m2.group(1))
        print(f'Found {len(new_ids)} novel transcripts')
        
        # 提取新转录本
        written = 0
        with open('__SNK_1__') as fa_in, open('__SNK_0__', 'w') as fa_out:
            write = False
            for line in fa_in:
                if line.startswith('>'):
                    tid = line[1:].split()[0]
                    write = tid in new_ids
                if write:
                    fa_out.write(line)
                    written += 1
        print(f'Extracted {written} novel transcript sequences')
        
        """
        r_code = r_code.replace("__SNK_0__", str(output.new_fa))
        r_code = r_code.replace("__SNK_1__", str(output.all_fa))
        r_code = r_code.replace("__SNK_2__", str(input.gtf))
        r_code = r_code.replace("__SNK_3__", str(input.fa))
        result = subprocess.run(["python3", "-c", r_code], text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- TransDecoder ORF 预测 ----
rule transdecoder_orf:
    """TransDecoder 预测开放阅读框"""
    input:
        fa = rules.extract_transcript_fasta.output.new_fa,
    output:
        pep    = f"{TMP_DIR}/04.transcripts/new_transcripts.fa.transdecoder.pep",
        cds    = f"{TMP_DIR}/04.transcripts/new_transcripts.fa.transdecoder.cds",
        gff3   = f"{TMP_DIR}/04.transcripts/new_transcripts.fa.transdecoder.gff3",
    threads: 8
    conda:
        f"{PIPELINE_DIR}/envs/transcript.yaml"
    shell:
        """
        # TransDecoder 会在输入文件所在目录写入临时文件,
        # 为避免污染 Snakemake 规划的路径, 先将输入复制到当前工作目录操作
        workdir={TMP_DIR}/04.transcripts
        mkdir -p $workdir
        fa_base=$(basename {input.fa})
        cp {input.fa} $workdir/$fa_base
        cd $workdir
        TransDecoder.LongOrfs -t $fa_base 2>&1 | tail -5
        TransDecoder.Predict \\
            -t $fa_base \\
            --cpu {threads} 2>&1 | tail -5
        # 将输出重命名为 Snakemake 期望的文件名
        mv ${{fa_base}}.transdecoder.pep {output.pep}
        mv ${{fa_base}}.transdecoder.cds {output.cds}
        mv ${{fa_base}}.transdecoder.gff3 {output.gff3}
        # 清理 TransDecoder 产生的中间文件
        rm -f $fa_base ${{fa_base}}.transdecoder_dir/ ${{fa_base}}.transdecoder.*.__checkpoints ${{fa_base}}.transdecoder.*.gff3 ${{fa_base}}.transdecoder.*.fasta ${{fa_base}}.transdecoder.*.scores ${{fa_base}}.transdecoder.*.train 2>/dev/null
        rm -rf ${{fa_base}}.transdecoder_dir 2>/dev/null || true
        """


# ---- 汇总新转录本到结果目录 ----
# 注: 必须放在 transdecoder_orf 之后，Snakemake 不支持前向 rules.X.output 引用
rule collect_transcript_results:
    """将新转录本和 ORF 结果复制到结果目录，供 all 规则和用户使用"""
    input:
        new_fa     = rules.extract_transcript_fasta.output.new_fa,
        all_fa     = rules.extract_transcript_fasta.output.all_fa,
        pep        = rules.transdecoder_orf.output.pep,
        cds        = rules.transdecoder_orf.output.cds,
        gff        = rules.transdecoder_orf.output.gff3,
        stats      = rules.gffcompare.output.stats,
        tracking   = rules.gffcompare.output.tracking,
        merged_gtf = rules.stringtie_merge.output.merged_gtf,
    output:
        novel_fa      = f"{RESULT_DIR}/04.transcripts/new_transcripts/new_transcripts.fa",
        novel_pep     = f"{RESULT_DIR}/04.transcripts/new_transcripts/new_transcripts.pep.fa",
        novel_cds     = f"{RESULT_DIR}/04.transcripts/new_transcripts/new_transcripts.cds.fa",
        novel_gff3    = f"{RESULT_DIR}/04.transcripts/new_transcripts/new_transcripts.gff3",
        novel_gtf     = f"{RESULT_DIR}/04.transcripts/new_transcripts/new_transcripts.gtf",
        cmp_stats     = f"{RESULT_DIR}/04.transcripts/new_transcripts/gffcompare_stats.txt",
        all_fa_out    = f"{RESULT_DIR}/04.transcripts/all_transcripts.fa",
        all_gtf_out   = f"{RESULT_DIR}/04.transcripts/all_transcripts.gtf",
        gene_trans_map = f"{RESULT_DIR}/04.transcripts/all.gene_trans_map.txt",
    conda:
        f"{PIPELINE_DIR}/envs/qc.yaml"
    run:
        import os, subprocess, sys
        os.makedirs(f"{RESULT_DIR}/04.transcripts/new_transcripts", exist_ok=True)
        # [TODO] 原 shell 命令: cp {input.new_fa}     {output.novel_fa}
        # [TODO] 原 shell 命令: cp {input.pep}        {output.novel_pep}
        # [TODO] 原 shell 命令: cp {input.cds}        {output.novel_cds}
        # [TODO] 原 shell 命令: cp {input.gff}        {output.novel_gff3}
        # [TODO] 原 shell 命令: cp {input.merged_gtf} {output.novel_gtf}
        # [TODO] 原 shell 命令: cp {input.stats}      {output.cmp_stats}
        # [TODO] 原 shell 命令: cp {input.all_fa}     {output.all_fa_out}
        # [TODO] 原 shell 命令: cp {input.merged_gtf} {output.all_gtf_out}
        r_code = """
        import re
        gene_trans = {}
        with open('__SNK_5__') as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.strip().split('\\t')
                if len(parts) < 9 or parts[2] not in ('transcript','mRNA'): continue
                gid_m = re.search(r'gene_id\\s+\"([^\"]+)\"', parts[8])
                tid_m = re.search(r'transcript_id\\s+\"([^\"]+)\"', parts[8])
                if gid_m and tid_m:
                    gid = gid_m.group(1)
                    tid = tid_m.group(1)
                    if gid not in gene_trans: gene_trans[gid] = set()
                    gene_trans[gid].add(tid)
        with open('__SNK_0__', 'w') as out:
            out.write('gene_id\\ttranscript_ids\\n')
            for gid, tids in sorted(gene_trans.items()):
                out.write(f'{gid}\\t{\",\".join(sorted(tids))}\\n')
        
        """
        r_code = r_code.replace("__SNK_0__", str(output.gene_trans_map))
        r_code = r_code.replace("__SNK_1__", str(output.all_gtf_out))
        r_code = r_code.replace("__SNK_2__", str(output.novel_gff3))
        r_code = r_code.replace("__SNK_3__", str(output.all_fa_out))
        r_code = r_code.replace("__SNK_4__", str(output.novel_cds))
        r_code = r_code.replace("__SNK_5__", str(input.merged_gtf))
        r_code = r_code.replace("__SNK_6__", str(output.novel_gtf))
        r_code = r_code.replace("__SNK_7__", str(output.novel_pep))
        r_code = r_code.replace("__SNK_8__", str(output.cmp_stats))
        r_code = r_code.replace("__SNK_9__", str(output.novel_fa))
        r_code = r_code.replace("__SNK_10__", str(input.all_fa))
        r_code = r_code.replace("__SNK_11__", str(input.new_fa))
        r_code = r_code.replace("__SNK_12__", str(input.stats))
        r_code = r_code.replace("__SNK_13__", str(RESULT_DIR))
        r_code = r_code.replace("__SNK_14__", str(input.cds))
        r_code = r_code.replace("__SNK_15__", str(input.pep))
        r_code = r_code.replace("__SNK_16__", str(input.gff))
        result = subprocess.run(["python3", "-c", r_code], text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- 新转录本 CDS 长度分布 ----
rule transcript_cds_length:
    """新转录本 CDS 长度分布图"""
    input:
        cds = rules.collect_transcript_results.output.novel_cds,
    output:
        pdf = f"{RESULT_DIR}/04.transcripts/new_transcripts/new_transcripts_cds.length.pdf",
        png = f"{RESULT_DIR}/04.transcripts/new_transcripts/new_transcripts_cds.length.png",
    conda:
        f"{PIPELINE_DIR}/envs/qc.yaml"
    run:
        import os, subprocess, sys
        r_code = """
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np
        
        lengths = []
        with open('__SNK_2__') as f:
            cur_len = 0
            for line in f:
                if line.startswith('>'):
                    if cur_len > 0:
                        lengths.append(cur_len)
                    cur_len = 0
                else:
                    cur_len += len(line.strip())
            if cur_len > 0:
                lengths.append(cur_len)
        
        if not lengths:
            print('No CDS sequences found')
            lengths = [0]
        
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # 直方图
        axes[0].hist(lengths, bins=50, color='#3498db', edgecolor='white', alpha=0.8)
        axes[0].set_xlabel('CDS Length (nt)')
        axes[0].set_ylabel('Frequency')
        axes[0].set_title(f'CDS Length Distribution\\n({len(lengths)} transcripts)')
        if len(lengths) > 1:
            axes[0].axvline(np.median(lengths), color='red', linestyle='--',
                            label=f'Median={np.median(lengths):.0f} nt')
            axes[0].legend()
        
        # 累积分布
        sorted_lens = np.sort(lengths)
        cumulative = np.arange(1, len(sorted_lens)+1) / len(sorted_lens)
        axes[1].plot(sorted_lens, cumulative, color='#e74c3c', lw=2)
        axes[1].set_xlabel('CDS Length (nt)')
        axes[1].set_ylabel('Cumulative Fraction')
        axes[1].set_title('CDS Length Cumulative Distribution')
        axes[1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('__SNK_0__', dpi=150, bbox_inches='tight')
        plt.savefig('__SNK_1__', dpi=150, bbox_inches='tight')
        print(f'CDS length plot: {len(lengths)} transcripts, median={np.median(lengths) if len(lengths)>1 else 0:.0f} nt')
        
        """
        r_code = r_code.replace("__SNK_0__", str(output.pdf))
        r_code = r_code.replace("__SNK_1__", str(output.png))
        r_code = r_code.replace("__SNK_2__", str(input.cds))
        result = subprocess.run(["python3", "-c", r_code], text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- 新转录本分类统计 ----
rule transcript_class_stats:
    """新转录本 class_code 分类统计图+表 (GffCompare)"""
    input:
        cmp_gtf = rules.gffcompare.output.compared,
    output:
        pdf = f"{RESULT_DIR}/04.transcripts/new_transcripts/new_code_stat.pdf",
        png = f"{RESULT_DIR}/04.transcripts/new_transcripts/new_code_stat.png",
        xls = f"{RESULT_DIR}/04.transcripts/new_transcripts/new_code_stat.xls",
    conda:
        f"{PIPELINE_DIR}/envs/qc.yaml"
    run:
        import os, subprocess, sys
        r_code = """
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import re
        from collections import Counter
        
        class_labels = {
            '=': 'Complete match',
            'c': 'Contained',
            'j': 'Novel isoform',
            'e': 'Single exon overlap',
            'i': 'Within intron',
            'o': 'Generic overlap',
            'p': 'Polymerase run-on',
            'u': 'Intergenic',
            'x': 'Opposite strand',
            's': 'Exon overlap opp strand',
        }
        
        # Parse GffCompare annotated GTF for class_code
        class_counts = Counter()
        with open('__SNK_0__') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                m = re.search(r'class_code\\s+"([^"]+)"', line)
                if m:
                    cc = m.group(1)
                    class_counts[cc] += 1
        
        if not class_counts:
            print('No class_code annotations found, creating placeholder')
            class_counts = {'=': 0, 'u': 0}
        
        # Sort by count desc
        sorted_cc = sorted(class_counts.items(), key=lambda x: -x[1])
        codes = [x[0] for x in sorted_cc]
        counts = [x[1] for x in sorted_cc]
        labels = [f'{c}\\n({class_labels.get(c,"?")})' for c in codes]
        
        # Colors
        cmap = plt.cm.Set3
        colors = [cmap(i/max(1,len(codes))) for i in range(len(codes))]
        
        fig, ax = plt.subplots(figsize=(10, 6))
        bars = ax.bar(range(len(codes)), counts, color=colors, edgecolor='white')
        ax.set_xticks(range(len(codes)))
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_ylabel('Number of Transcripts')
        ax.set_title('Novel Transcript Classification (GffCompare class_code)')
        
        # Add value labels on bars
        for bar, count in zip(bars, counts):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(counts)*0.01,
                    str(count), ha='center', va='bottom', fontsize=9)
        
        plt.tight_layout()
        plt.savefig('__SNK_2__', dpi=150, bbox_inches='tight')
        plt.savefig('__SNK_3__', dpi=150, bbox_inches='tight')
        
        # Write XLS
        with open('__SNK_1__', 'w') as out:
            out.write('Class\\tLabel\\tCount\\n')
            for cc, count in sorted_cc:
                out.write(f'{cc}\\t{class_labels.get(cc,"?")}\\t{count}\\n')
        
        print(f'Classification stats: {sum(counts)} transcripts, {len(sorted_cc)} classes')
        
        """
        r_code = r_code.replace("__SNK_0__", str(input.cmp_gtf))
        r_code = r_code.replace("__SNK_1__", str(output.xls))
        r_code = r_code.replace("__SNK_2__", str(output.pdf))
        r_code = r_code.replace("__SNK_3__", str(output.png))
        result = subprocess.run(["python3", "-c", r_code], text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- 全部转录本长度分布 ----
rule transcript_length_stats:
    """全部转录本长度分布图+表"""
    input:
        all_fa = rules.collect_transcript_results.output.all_fa_out,
    output:
        pdf = f"{RESULT_DIR}/04.transcripts/all_transcripts.length.pdf",
        png = f"{RESULT_DIR}/04.transcripts/all_transcripts.length.png",
        txt = f"{RESULT_DIR}/04.transcripts/all_transcripts.length.txt",
    conda:
        f"{PIPELINE_DIR}/envs/qc.yaml"
    run:
        import os, subprocess, sys
        r_code = """
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np
        
        lengths = []
        with open('__SNK_0__') as f:
            cur_len = 0
            for line in f:
                if line.startswith('>'):
                    if cur_len > 0:
                        lengths.append(cur_len)
                    cur_len = 0
                else:
                    cur_len += len(line.strip())
            if cur_len > 0:
                lengths.append(cur_len)
        
        if not lengths:
            lengths = [0]
        
        lengths_arr = np.array(lengths)
        
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        axes[0].hist(lengths, bins=80, color='#2ecc71', edgecolor='white', alpha=0.8)
        axes[0].set_xlabel('Transcript Length (nt)')
        axes[0].set_ylabel('Frequency')
        axes[0].set_title(f'All Transcript Length Distribution\\n({len(lengths)} transcripts)')
        if len(lengths) > 1:
            axes[0].axvline(np.median(lengths), color='red', linestyle='--',
                            label=f'Median={np.median(lengths):.0f} nt')
            axes[0].legend()
        
        sorted_lens = np.sort(lengths)
        cumulative = np.arange(1, len(sorted_lens)+1) / len(sorted_lens)
        axes[1].plot(sorted_lens, cumulative, color='#e74c3c', lw=2)
        axes[1].set_xlabel('Transcript Length (nt)')
        axes[1].set_ylabel('Cumulative Fraction')
        axes[1].set_title('Cumulative Length Distribution')
        axes[1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('__SNK_2__', dpi=150, bbox_inches='tight')
        plt.savefig('__SNK_3__', dpi=150, bbox_inches='tight')
        
        # Write TXT
        with open('__SNK_1__', 'w') as out:
            out.write('total_transcripts\\tmin_len\\tmax_len\\tmean_len\\tmedian_len\\n')
            out.write(f'{len(lengths)}\\t{lengths_arr.min():.0f}\\t{lengths_arr.max():.0f}\\t{lengths_arr.mean():.1f}\\t{np.median(lengths_arr):.0f}\\n')
        print(f'Length stats: {len(lengths)} transcripts, mean={lengths_arr.mean():.0f} nt')
        
        """
        r_code = r_code.replace("__SNK_0__", str(input.all_fa))
        r_code = r_code.replace("__SNK_1__", str(output.txt))
        r_code = r_code.replace("__SNK_2__", str(output.pdf))
        r_code = r_code.replace("__SNK_3__", str(output.png))
        result = subprocess.run(["python3", "-c", r_code], text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- 基因/转录本数量统计 ----
rule gene_transcript_number_stats:
    """基因和转录本数量统计图+表"""
    input:
        gene_trans_map = rules.collect_transcript_results.output.gene_trans_map,
    output:
        number_pdf  = f"{RESULT_DIR}/04.transcripts/Number_stats.pdf",
        number_png  = f"{RESULT_DIR}/04.transcripts/Number_stats.png",
        number_txt  = f"{RESULT_DIR}/04.transcripts/Number_stats.txt",
        trans_per_gene_pdf = f"{RESULT_DIR}/04.transcripts/trans_for_gene_number_stats.pdf",
        trans_per_gene_png = f"{RESULT_DIR}/04.transcripts/trans_for_gene_number_stats.png",
        trans_per_gene_txt = f"{RESULT_DIR}/04.transcripts/trans_for_gene_number_stats.txt",
    conda:
        f"{PIPELINE_DIR}/envs/qc.yaml"
    run:
        import os, subprocess, sys
        r_code = """
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        
        # Parse gene-transcript mapping
        gene_trans = {}
        with open('__SNK_3__') as f:
            next(f)  # skip header
            for line in f:
                parts = line.strip().split('\\t')
                if len(parts) >= 2:
                    gene_trans[parts[0]] = parts[1].split(',') if parts[1] else []
        
        num_genes = len(gene_trans)
        all_trans = set()
        for tids in gene_trans.values():
            all_trans.update(tids)
        num_trans = len(all_trans)
        
        # Transcripts per gene distribution
        trans_per_gene = [len(tids) for tids in gene_trans.values()]
        from collections import Counter
        tpg_counts = Counter(trans_per_gene)
        
        # --- Number stats ---
        fig1, ax1 = plt.subplots(figsize=(6, 5))
        bars1 = ax1.bar(['Genes', 'Transcripts'], [num_genes, num_trans],
                        color=['#3498db', '#e74c3c'], edgecolor='white', width=0.5)
        for bar, val in zip(bars1, [num_genes, num_trans]):
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(num_genes, num_trans)*0.01,
                     str(val), ha='center', fontsize=14, fontweight='bold')
        ax1.set_title('Gene and Transcript Counts')
        ax1.set_ylabel('Count')
        plt.tight_layout()
        plt.savefig('__SNK_6__', dpi=150, bbox_inches='tight')
        plt.savefig('__SNK_4__', dpi=150, bbox_inches='tight')
        
        with open('__SNK_5__', 'w') as out:
            out.write(f'num_genes\\t{num_genes}\\n')
            out.write(f'num_transcripts\\t{num_trans}\\n')
        
        # --- Transcripts per gene ---
        fig2, ax2 = plt.subplots(figsize=(10, 5))
        tpg_sorted = sorted(tpg_counts.items())
        x_vals = [x[0] for x in tpg_sorted]
        y_vals = [x[1] for x in tpg_sorted]
        ax2.bar(x_vals, y_vals, color='#9b59b6', edgecolor='white')
        ax2.set_xlabel('Number of Transcripts per Gene')
        ax2.set_ylabel('Number of Genes')
        ax2.set_title(f'Transcripts per Gene Distribution\\n({num_genes} genes, {num_trans} transcripts)')
        ax2.set_xticks(x_vals)
        plt.tight_layout()
        plt.savefig('__SNK_0__', dpi=150, bbox_inches='tight')
        plt.savefig('__SNK_1__', dpi=150, bbox_inches='tight')
        
        with open('__SNK_2__', 'w') as out:
            out.write('transcripts_per_gene\\tnum_genes\\n')
            for k, v in tpg_sorted:
                out.write(f'{k}\\t{v}\\n')
        
        print(f'Number stats: {num_genes} genes, {num_trans} transcripts')
        
        """
        r_code = r_code.replace("__SNK_0__", str(output.trans_per_gene_pdf))
        r_code = r_code.replace("__SNK_1__", str(output.trans_per_gene_png))
        r_code = r_code.replace("__SNK_2__", str(output.trans_per_gene_txt))
        r_code = r_code.replace("__SNK_3__", str(input.gene_trans_map))
        r_code = r_code.replace("__SNK_4__", str(output.number_png))
        r_code = r_code.replace("__SNK_5__", str(output.number_txt))
        r_code = r_code.replace("__SNK_6__", str(output.number_pdf))
        result = subprocess.run(["python3", "-c", r_code], text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")
