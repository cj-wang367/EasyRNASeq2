# ============================================================
# data_ready 模块: 数据预处理
# 功能:
#   1. 多 lane FASTQ 合并
#   2. 数据统计 (reads 数, 碱基数)
#   3. 抽小数据 (50000 reads) 用于 rRNA 和基因组比对快速评估
#   4. rRNA 评估 (Bowtie2 vs Silva/rRNA DB)
#   5. 比对率评估 (HISAT2 vs 参考基因组)
#   6. 基于比对率和数据量决定 cut/extract/raw 策略
#   7. 数据截取 (按合同量) 或 数据提取 (低比对率)
# ============================================================

# ---- 多 lane FASTQ 合并 ----
rule merge_lanes:
    """合并同一样本的多 lane FASTQ 文件"""
    input:
        r1_files = lambda wildcards: SAMPLES_DATA[wildcards.sample]["r1_files"],
        r2_files = lambda wildcards: SAMPLES_DATA[wildcards.sample]["r2_files"],
    output:
        r1_merged = f"{TMP_DIR}/data_ready/merged/{{sample}}/{{sample}}.R1.fq.gz",
        r2_merged = f"{TMP_DIR}/data_ready/merged/{{sample}}/{{sample}}.R2.fq.gz",
    params:
        outdir = f"{TMP_DIR}/data_ready/merged/{{sample}}",
    conda:
        f"{PIPELINE_DIR}/envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        # R1 文件: 使用 set -- 正确计数空格分隔的 Snakemake 列表
        set -- {input.r1_files}
        if [ $# -eq 1 ]; then
            ln -sf "$(realpath "$1")" {output.r1_merged}
        else
            cat "$@" > {output.r1_merged}
        fi
        # R2 文件
        set -- {input.r2_files}
        if [ $# -eq 1 ]; then
            ln -sf "$(realpath "$1")" {output.r2_merged}
        else
            cat "$@" > {output.r2_merged}
        fi
        """


# ---- 抽取少量 reads (用于快速评估) ----
rule sample_reads:
    """抽取少量 reads 用于 rRNA 和比对率评估"""
    input:
        r1 = rules.merge_lanes.output.r1_merged,
        r2 = rules.merge_lanes.output.r2_merged,
    output:
        r1_sample = f"{TMP_DIR}/data_ready/sample/{{sample}}/{{sample}}.sample.R1.fq.gz",
        r2_sample = f"{TMP_DIR}/data_ready/sample/{{sample}}/{{sample}}.sample.R2.fq.gz",
    params:
        n_reads   = config["data_preprocess"]["cut_sample_reads"],
        outdir    = f"{TMP_DIR}/data_ready/sample/{{sample}}",
    threads: 4
    conda:
        f"{PIPELINE_DIR}/envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        seqkit sample -n {params.n_reads} -s 42 -2 \\
            -j {threads} \\
            -o {output.r1_sample} {input.r1}
        seqkit sample -n {params.n_reads} -s 42 -2 \\
            -j {threads} \\
            -o {output.r2_sample} {input.r2}
        """


# ---- 完整数据统计 ----
rule data_stat:
    """统计每个样本的总 reads 和总碱基数"""
    input:
        r1 = rules.merge_lanes.output.r1_merged,
        r2 = rules.merge_lanes.output.r2_merged,
    output:
        stat = f"{TMP_DIR}/data_ready/stat/{{sample}}.stat.txt",
    conda:
        f"{PIPELINE_DIR}/envs/qc.yaml"
    shell:
        """
        mkdir -p {TMP_DIR}/data_ready/stat
        r1_reads=$(echo $(seqkit stats -T {input.r1} | tail -1 | awk -F'\\t' '{{print $4}}'))
        r1_bases=$(echo $(seqkit stats -T {input.r1} | tail -1 | awk -F'\\t' '{{print $5}}'))
        r2_bases=$(echo $(seqkit stats -T {input.r2} | tail -1 | awk -F'\\t' '{{print $5}}'))
        total_bases=$((r1_bases + r2_bases))
        total_bases_g=$(echo "scale=2; $total_bases / 1000000000" | bc)
        echo -e "sample\\treads\\tbases\\tbases_G" > {output.stat}
        echo -e "{wildcards.sample}\\t$r1_reads\\t$total_bases\\t$total_bases_g" >> {output.stat}
        """


# ---- 合并数据统计 ----
rule merge_data_stats:
    """汇总所有样本的数据统计"""
    input:
        stats = expand(f"{TMP_DIR}/data_ready/stat/{{sample}}.stat.txt", sample=SAMPLES),
    output:
        all_stat = f"{RESULT_DIR}/01.raw_data/data_stat_all.xls",
    conda:
        f"{PIPELINE_DIR}/envs/qc.yaml"
    shell:
        """
        head -1 $(echo {input.stats} | awk '{{print $1}}') > {output.all_stat}
        for f in {input.stats}; do
            tail -n +2 $f >> {output.all_stat}
        done
        """


# ---- rRNA 评估 (Bowtie2 vs rRNA database) ----
rule rrna_evaluation:
    """Bowtie2 比对 rRNA 数据库评估 rRNA 残留"""
    input:
        r1 = rules.sample_reads.output.r1_sample,
        r2 = rules.sample_reads.output.r2_sample,
    output:
        rrna_rate = f"{TMP_DIR}/data_ready/rrna/{{sample}}.rrna_rate.txt",
    params:
        rrna_db    = config["paths"].get("rrna_db", ""),
        outdir     = f"{TMP_DIR}/data_ready/rrna/{{sample}}",
        max_rrna   = config["data_preprocess"]["max_rrna_rate"],
    threads: 4
    conda:
        f"{PIPELINE_DIR}/envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir}

        if [ -n "{params.rrna_db}" ] && [ -f "{params.rrna_db}" ]; then
            # 如果 rRNA DB 未建索引，先建
            if [ ! -f "{params.rrna_db}.1.bt2" ]; then
                bowtie2-build --threads {threads} {params.rrna_db} {params.rrna_db}
            fi

            bowtie2 -x {params.rrna_db} \\
                -1 {input.r1} -2 {input.r2} \\
                -p {threads} \\
                --very-sensitive-local \\
                --no-unal \\
                -S {params.outdir}/rrna.sam 2> {params.outdir}/bowtie2.log

            total=$(grep "reads; of these:" {params.outdir}/bowtie2.log | awk '{{print $1}}')
            aligned=$(grep "aligned concordantly exactly 1 time" {params.outdir}/bowtie2.log | awk '{{print $1}}')
            if [ -z "$aligned" ]; then aligned=0; fi
            if [ -z "$total" ] || [ "$total" -eq 0 ]; then total=1; fi
            rate=$(echo "scale=4; $aligned / $total * 100" | bc)
            echo -e "sample\\trrna_reads\\ttotal_reads\\trrna_rate" > {output.rrna_rate}
            echo -e "{wildcards.sample}\\t$aligned\\t$total\\t$rate" >> {output.rrna_rate}

            if (( $(echo "$rate > {params.max_rrna}" | bc -l) )); then
                echo "WARNING: High rRNA rate ($rate%) for {wildcards.sample}" >> {output.rrna_rate}
            fi
        else
            echo "rRNA database not configured, skipping" > {output.rrna_rate}
            echo -e "sample\\trrna_reads\\ttotal_reads\\trrna_rate" > {output.rrna_rate}
            echo -e "{wildcards.sample}\\t0\\t0\\t0" >> {output.rrna_rate}
        fi
        """


# ---- HISAT2 索引前缀 (预构建 > 自动构建) ----
def _get_hisat2_index_prefix():
    """获取有效的 HISAT2 索引前缀"""
    import glob as _g2
    idx = config["alignment"].get("index")
    if idx and str(idx) not in ("", "null", "None"):
        idx = str(idx)
        # 目录: 扫描 *.1.ht2 推断前缀
        if os.path.isdir(idx):
            ht2_files = sorted(_g2.glob(os.path.join(idx, "*.1.ht2")))
            if ht2_files:
                return ht2_files[0].rsplit(".1.ht2", 1)[0]
        return idx
    return f"{TMP_DIR}/index/hisat2/genome.fa"


# ---- 比对率评估 (用小数据快速比对) ----
rule mapping_rate_evaluation:
    """HISAT2 比对小数据评估基因组比对率"""
    input:
        r1    = rules.sample_reads.output.r1_sample,
        r2    = rules.sample_reads.output.r2_sample,
        index = _get_hisat2_index_prefix() + ".1.ht2",
    output:
        map_rate = f"{TMP_DIR}/data_ready/mapping/{{sample}}.map_rate.txt",
    params:
        idx_prefix   = _get_hisat2_index_prefix(),
        outdir       = f"{TMP_DIR}/data_ready/mapping/{{sample}}",
        target_map   = config["data_preprocess"]["target_mapped_ratio"],
        min_map      = config["data_preprocess"]["min_mapped_ratio"],
    threads: 4
    conda:
        f"{PIPELINE_DIR}/envs/align.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        hisat2 -x {params.idx_prefix} \\
            -1 {input.r1} -2 {input.r2} \\
            -p {threads} \\
            --no-unal \\
            -S {params.outdir}/test.sam 2> {params.outdir}/hisat2.log

        rate=$(grep "overall alignment rate" {params.outdir}/hisat2.log | awk '{{print $1}}' | tr -d '%')
        echo -e "sample\\tmapping_rate" > {output.map_rate}
        echo -e "{wildcards.sample}\\t$rate" >> {output.map_rate}

        target={params.target_map}
        min_map={params.min_map}
        if (( $(echo "$rate >= $target" | bc -l) )); then
            echo "OK: mapping rate $rate% >= target $target%" >> {output.map_rate}
        elif (( $(echo "$rate >= $min_map" | bc -l) )); then
            echo "WARNING: mapping rate $rate% < target $target% but >= min $min_map%" >> {output.map_rate}
        else
            echo "LOW: mapping rate $rate% < min $min_map% -- consider data extraction" >> {output.map_rate}
        fi
        """


# ---- 数据截取策略判定 ----
rule determine_strategy:
    """根据数据量和比对率决定 cut/extract/raw 策略"""
    input:
        data_stat  = f"{TMP_DIR}/data_ready/stat/{{sample}}.stat.txt",
        map_rate   = f"{TMP_DIR}/data_ready/mapping/{{sample}}.map_rate.txt",
    output:
        strategy   = f"{TMP_DIR}/data_ready/strategy/{{sample}}.strategy.txt",
    params:
        cut_base   = config["data_preprocess"]["bases_data"],
        lower      = config["data_preprocess"]["lower_bound_ratio"],
        upper      = config["data_preprocess"]["upper_bound_ratio"],
        target_map = config["data_preprocess"]["target_mapped_ratio"],
        min_map    = config["data_preprocess"]["min_mapped_ratio"],
    conda:
        f"{PIPELINE_DIR}/envs/qc.yaml"
    run:
        import subprocess

        # 解析数据量
        data_stat_file = str(input.data_stat)
        map_rate_file = str(input.map_rate)

        with open(data_stat_file) as f:
            lines = f.readlines()
            if len(lines) > 1:
                parts = lines[1].strip().split("\t")
                bases_g = float(parts[3])
            else:
                bases_g = 0

        with open(map_rate_file) as f:
            lines = f.readlines()
            if len(lines) > 1:
                parts = lines[1].strip().split("\t")
                map_rate_val = float(parts[1])
            else:
                map_rate_val = 0

        cut_base = params.cut_base
        lower = params.lower
        upper = params.upper
        target = params.target_map
        min_map = params.min_map

        # 策略逻辑
        strategy = "raw"
        data_status = "within_range"

        if bases_g > cut_base * upper:
            data_status = "too_much"
        elif bases_g < cut_base * lower:
            data_status = "insufficient"

        if map_rate_val < min_map:
            strategy = "extract"
        elif data_status == "too_much":
            # run_cut 默认 True (config.data_preprocess 无此字段时直接裁剪)
            strategy = "cut"

        with open(str(output.strategy), "w") as f:
            f.write(f"sample\t{wildcards.sample}\n")
            f.write(f"bases_G\t{bases_g:.2f}\n")
            f.write(f"mapping_rate\t{map_rate_val:.2f}\n")
            f.write(f"data_status\t{data_status}\n")
            f.write(f"strategy\t{strategy}\n")

        print(f"[{wildcards.sample}] bases={bases_g:.2f}G, map_rate={map_rate_val:.2f}%, strategy={strategy}")


# ---- 数据截取 (按合同量) ----
rule cut_data:
    """按合同数据量截取 reads"""
    input:
        r1       = rules.merge_lanes.output.r1_merged,
        r2       = rules.merge_lanes.output.r2_merged,
        strategy = rules.determine_strategy.output.strategy,
    output:
        r1_cut   = f"{RESULT_DIR}/02.clean_data/{{sample}}/{{sample}}.R1.final.fq.gz",
        r2_cut   = f"{RESULT_DIR}/02.clean_data/{{sample}}/{{sample}}.R2.final.fq.gz",
    params:
        cut_base   = config["data_preprocess"]["bases_data"],
        read_len   = 150,  # 默认 reads 长度
        outdir     = f"{RESULT_DIR}/02.clean_data/{{sample}}",
    threads: 4
    conda:
        f"{PIPELINE_DIR}/envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir}

        # 读取策略
        strategy=$(grep "^strategy" {input.strategy} | awk -F'\\t' '{{print $2}}')

        if [ "$strategy" == "cut" ]; then
            # 计算需要的 reads 数
            reads_per_file=$(echo "{params.cut_base} * 1000000000 / {params.read_len} / 2" | bc)
            seqkit sample -n $reads_per_file -s 666 -2 \\
                -j {threads} \\
                -o {output.r1_cut} {input.r1}
            seqkit sample -n $reads_per_file -s 666 -2 \\
                -j {threads} \\
                -o {output.r2_cut} {input.r2}
        else
            # raw/extract: 直接使用原始数据
            ln -sf "$(realpath {input.r1})" {output.r1_cut}
            ln -sf "$(realpath {input.r2})" {output.r2_cut}
        fi
        """
