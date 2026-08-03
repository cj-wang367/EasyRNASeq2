# ============================================================
# QC 模块: 数据质控 (fastp + fastqc + seqkit)
# 输入: data_ready 模块的 cut_data 输出 (或原始 FASTQ)
# ============================================================

# ---- 获取 fastp 输入文件 ----
# 如果 data_ready 模块运行了，使用其输出；否则使用原始 FASTQ
def get_fastp_input_r1(wildcards):
    """获取 fastp 的 R1 输入 (多文件需先合并)"""
    if config.get("data_preprocess", {}).get("run", False):
        return rules.cut_data.output.r1_cut
    elif len(SAMPLES_DATA[wildcards.sample]["r1_files"]) > 1:
        return rules.merge_lanes.output.r1_merged
    else:
        return SAMPLES_DATA[wildcards.sample]["r1_files"][0]

def get_fastp_input_r2(wildcards):
    """获取 fastp 的 R2 输入 (多文件需先合并)"""
    if config.get("data_preprocess", {}).get("run", False):
        return rules.cut_data.output.r2_cut
    elif len(SAMPLES_DATA[wildcards.sample]["r2_files"]) > 1:
        return rules.merge_lanes.output.r2_merged
    else:
        return SAMPLES_DATA[wildcards.sample]["r2_files"][0]


# ---- 数据统计 ----
rule seqkit_stat:
    """统计原始数据 reads 数和碱基数"""
    input:
        r1 = lambda w: SAMPLES_DATA[w.sample]["r1_files"][0],
        r2 = lambda w: SAMPLES_DATA[w.sample]["r2_files"][0],
    output:
        stat = f"{RESULT_DIR}/01.raw_data/{{sample}}/{{sample}}.seqkit_stat.txt",
    params:
        outdir = f"{RESULT_DIR}/01.raw_data/{{sample}}",
    conda:
        f"{PIPELINE_DIR}/envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        seqkit stat -j 4 {input.r1} {input.r2} > {output.stat}
        """


# ---- 数据清洗: fastp ----
rule fastp_trim:
    """fastp 去接头 + 质控过滤 (输入来自 data_ready 或原始 FASTQ)"""
    input:
        r1 = get_fastp_input_r1,
        r2 = get_fastp_input_r2,
    output:
        r1_clean = f"{RESULT_DIR}/02.clean_data/{{sample}}/{{sample}}.R1.fastp.fq.gz",
        r2_clean = f"{RESULT_DIR}/02.clean_data/{{sample}}/{{sample}}.R2.fastp.fq.gz",
        json      = f"{RESULT_DIR}/02.clean_data/{{sample}}/{{sample}}.fastp.json",
        html      = f"{RESULT_DIR}/02.clean_data/{{sample}}/{{sample}}.fastp.html",
    params:
        min_len      = config["qc"]["fastp_min_len"],
        adapter_r1   = config["qc"]["adapter_r1"],
        adapter_r2   = config["qc"]["adapter_r2"],
        phred        = config["qc"]["phred"],
        trim_poly_g  = "--trim_poly_g --poly_g_min_len " + str(config["qc"]["poly_g_min_len"]) if config["qc"]["trim_poly_g"] else "",
        outdir       = f"{RESULT_DIR}/02.clean_data/{{sample}}",
    threads: 4
    conda:
        f"{PIPELINE_DIR}/envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        echo "[INFO] fastp input R1: {input.r1}"
        echo "[INFO] fastp input R2: {input.r2}"
        fastp \\
            -i {input.r1} -I {input.r2} \\
            -o {output.r1_clean} -O {output.r2_clean} \\
            -j {output.json} -h {output.html} \\
            -l {params.min_len} \\
            {params.trim_poly_g} \\
            -w {threads} \\
            --detect_adapter_for_pe \\
            2>&1 | tee {params.outdir}/{wildcards.sample}.fastp.log

        # 检查输出
        if [ ! -s {output.r1_clean} ] || [ ! -s {output.r2_clean} ]; then
            echo "ERROR: fastp 输出文件为空!"
            exit 1
        fi
        """


# ---- FastQC 质控 ----
rule fastqc:
    """FastQC 质量评估"""
    input:
        r1 = rules.fastp_trim.output.r1_clean,
        r2 = rules.fastp_trim.output.r2_clean,
    output:
        html_r1 = f"{RESULT_DIR}/02.clean_data/{{sample}}/{{sample}}.R1.fastp_fastqc.html",
        html_r2 = f"{RESULT_DIR}/02.clean_data/{{sample}}/{{sample}}.R2.fastp_fastqc.html",
        zip_r1  = f"{RESULT_DIR}/02.clean_data/{{sample}}/{{sample}}.R1.fastp_fastqc.zip",
        zip_r2  = f"{RESULT_DIR}/02.clean_data/{{sample}}/{{sample}}.R2.fastp_fastqc.zip",
    params:
        outdir = f"{RESULT_DIR}/02.clean_data/{{sample}}",
    threads: 4
    conda:
        f"{PIPELINE_DIR}/envs/qc.yaml"
    shell:
        """
        fastqc -t {threads} -o {params.outdir} {input.r1} {input.r2}
        """


# ---- MultiQC 汇总 ----
rule multiqc:
    """MultiQC 汇总所有 QC 结果"""
    input:
        fastqc_r1 = expand(f"{RESULT_DIR}/02.clean_data/{{sample}}/{{sample}}.R1.fastp_fastqc.html", sample=SAMPLES),
    output:
        report = f"{RESULT_DIR}/02.clean_data/multiqc_report.html",
    params:
        outdir = f"{RESULT_DIR}/02.clean_data",
    conda:
        f"{PIPELINE_DIR}/envs/qc.yaml"
    shell:
        """
        multiqc -f -o {params.outdir} {params.outdir}
        """
