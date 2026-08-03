# ============================================================
# 输入文件预检模块
#
# 所有正式分析前的前置检查, 确保输入文件完整可靠。
# 如果检查失败, 流程立即停止, 避免运行到一半才发现输入问题。
#
# 检查项:
#   1. 基因组 FASTA — 格式 / 碱基组成 / 染色体
#   2. GTF 注释     — gene_id / transcript_id / GFF3 伪装检测
#   3. 样本信息表   — TSV 格式 / 必需列 / FASTQ 文件存在
#   4. 比较组文件   — 格式正确 / 组名匹配样本
#   5. 跨文件互检   — 染色体交叉匹配
# ============================================================

rule file_check:
    """输入文件完整性预检 — 所有规则的前置依赖"""
    input:
        ref_fa      = config["paths"]["ref_fa"],
        ref_gtf     = config["paths"]["ref_gtf"],
        sample_info = config["paths"]["sample_info"],
        compare_raw = config["paths"].get("compare_file", ""),
    output:
        check_done = f"{TMP_DIR}/file_check/check.done",
    params:
        outdir   = f"{TMP_DIR}/file_check",
        species  = config.get("species", "unknown"),
        compare_file = config["paths"].get("compare_file", ""),
    conda:
        f"{PIPELINE_DIR}/envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        python3 {PIPELINE_DIR}/scripts/file_check.py \
            --ref_fa {input.ref_fa} \
            --ref_gtf {input.ref_gtf} \
            --sample_info {input.sample_info} \
            --compare_file {params.compare_file} \
            --species "{params.species}" \
            --outdir {params.outdir}
        """
