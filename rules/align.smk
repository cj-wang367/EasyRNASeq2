# ============================================================
# 比对模块: HISAT2/STAR 索引构建 + 比对 + 排序
#
# 索引策略 (参考 Nextflow data_ready/map_hisat2 设计):
#   1. 优先使用预构建索引 (config alignment.index 或 ref_anno_dir)
#   2. 不存在时自动构建: hisat2_index_onlyFa=false → 含剪接位点
#                        hisat2_index_onlyFa=true  → 仅 FASTA
#   3. 索引输出到 tmp/index/ 下, 可手动复制到 ref_anno_dir 持久化
# ============================================================

import glob as _gm

# ---- 预检: 是否使用已有索引 ----
_HAS_PREBUILT_INDEX = False
_IDX_DIR = ""
_index_cfg = config["alignment"].get("index")
if _index_cfg and str(_index_cfg) not in ("", "null", "None"):
    _IDX_DIR = str(_index_cfg)
    _tool = config["alignment"]["tool"]
    if _tool == "hisat2":
        if os.path.isdir(_IDX_DIR):
            _ht2 = _gm.glob(os.path.join(_IDX_DIR, "*.1.ht2"))
            if _ht2:
                _HAS_PREBUILT_INDEX = True
        elif os.path.exists(_IDX_DIR):
            _HAS_PREBUILT_INDEX = True
        if _HAS_PREBUILT_INDEX:
            print(f"[INFO] 使用已有 HISAT2 索引: {_IDX_DIR}")
    elif _tool == "star":
        if os.path.isfile(os.path.join(_IDX_DIR, "Genome")):
            _HAS_PREBUILT_INDEX = True
            print(f"[INFO] 使用已有 STAR 索引: {_IDX_DIR}")

if not _HAS_PREBUILT_INDEX:
    print(f"[INFO] 未找到预构建索引, 将自动构建 {config['alignment']['tool']} 索引")


# ---- HISAT2 索引构建 (合并小基因组/大基因组为单规则) ----
if not _HAS_PREBUILT_INDEX and config["alignment"]["tool"] == "hisat2":

    rule hisat2_build_index:
        """构建 HISAT2 索引 — 根据 hisat2_index_onlyFa 决定是否含剪接信息"""
        input:
            fa  = config["paths"]["ref_fa"],
            gtf = config["paths"]["ref_gtf"],
        output:
            idx = touch(f"{TMP_DIR}/index/hisat2/genome.fa.1.ht2"),
        params:
            outdir  = f"{TMP_DIR}/index/hisat2",
            only_fa = config["alignment"].get("hisat2_index_onlyFa", False),
        threads: 6
        conda:
            f"{PIPELINE_DIR}/envs/align.yaml"
        shell:
            """
            mkdir -p {params.outdir}
            if [ "{params.only_fa}" = "True" ] || [ "{params.only_fa}" = "true" ]; then
                echo "[INFO] 构建 HISAT2 索引 (仅 FASTA 模式 — 大基因组)"
                hisat2-build -p {threads} {input.fa} {params.outdir}/genome.fa
            else
                echo "[INFO] 构建 HISAT2 索引 (含剪接位点+外显子 — 小基因组)"
                hisat2_extract_splice_sites.py {input.gtf} > {params.outdir}/splicesites.txt
                hisat2_extract_exons.py {input.gtf} > {params.outdir}/exons.txt
                hisat2-build -p {threads} \\
                    --ss {params.outdir}/splicesites.txt \\
                    --exon {params.outdir}/exons.txt \\
                    {input.fa} {params.outdir}/genome.fa
            fi
            echo "[INFO] HISAT2 索引构建完成: {params.outdir}/genome.fa"
            """


# ---- STAR 索引构建 ----
if not _HAS_PREBUILT_INDEX and config["alignment"]["tool"] == "star":

    rule star_build_index:
        """构建 STAR 索引"""
        input:
            fa  = config["paths"]["ref_fa"],
            gtf = config["paths"]["ref_gtf"],
        output:
            idx = touch(f"{TMP_DIR}/index/star/Genome"),
        params:
            outdir = f"{TMP_DIR}/index/star",
            mem_gb = config["alignment"]["star_mem_gb"],
        threads: 6
        conda:
            f"{PIPELINE_DIR}/envs/align.yaml"
        shell:
            """
            mkdir -p {params.outdir}
            STAR \\
                --runMode genomeGenerate \\
                --runThreadN {threads} \\
                --genomeDir {params.outdir} \\
                --genomeFastaFiles {input.fa} \\
                --sjdbGTFfile {input.gtf} \\
                --sjdbOverhang 149 \\
                --limitGenomeGenerateRAM {params.mem_gb}000000000
            echo "[INFO] STAR 索引构建完成: {params.outdir}"
            """


# ---- 索引路径解析 ----
def _get_hisat2_index_input(wildcards):
    """返回 HISAT2 索引 .1.ht2 文件路径 (自动判断预构建/自构建)"""
    if _HAS_PREBUILT_INDEX:
        return sorted(_gm.glob(os.path.join(_IDX_DIR, "*.1.ht2")))
    return [f"{TMP_DIR}/index/hisat2/genome.fa.1.ht2"]

def _get_star_index_input(wildcards):
    """返回 STAR 索引 Genome 文件路径"""
    if _HAS_PREBUILT_INDEX:
        return [os.path.join(_IDX_DIR, "Genome")]
    return [f"{TMP_DIR}/index/star/Genome"]

def _get_align_index_dir():
    """返回比对索引目录"""
    if _HAS_PREBUILT_INDEX:
        return _IDX_DIR
    _tool = config["alignment"]["tool"]
    return f"{TMP_DIR}/index/{_tool}"


# ---- HISAT2 比对 ----
if config["alignment"]["tool"] == "hisat2":

    rule hisat2_align:
        """HISAT2 双端比对"""
        input:
            r1    = rules.fastp_trim.output.r1_clean,
            r2    = rules.fastp_trim.output.r2_clean,
            index = _get_hisat2_index_input,
        output:
            sam  = f"{TMP_DIR}/03.align/{{sample}}/{{sample}}.sam",
            stat = f"{RESULT_DIR}/03.align/{{sample}}/{{sample}}.align_stat.txt",
        params:
            idx_dir = _get_align_index_dir(),
            outdir  = f"{RESULT_DIR}/03.align/{{sample}}",
            phred   = config["qc"]["phred"],
            strand  = {
                "unstranded": "",
                "firststrand": "--rna-strandness RF",
                "secondstrand": "--rna-strandness FR",
            }[config["alignment"]["strand"]],
        threads: 6
        conda:
            f"{PIPELINE_DIR}/envs/align.yaml"
        shell:
            """
            mkdir -p {params.outdir} {TMP_DIR}/03.align/{wildcards.sample}
            IDX_DIR="{params.idx_dir}"
            # Find index basename (strip .1.ht2 suffix)
            IDX_BASE=$(ls "$IDX_DIR"/*.1.ht2 2>/dev/null | head -1 | sed 's/\\.1\\.ht2$//')
            if [ -z "$IDX_BASE" ]; then
                echo "ERROR: HISAT2 index not found in $IDX_DIR"
                echo "Files in directory:"; ls -la "$IDX_DIR"/
                exit 1
            fi
            echo "[INFO] Using HISAT2 index: $IDX_BASE"
            hisat2 -x "$IDX_BASE" \
                -1 {input.r1} -2 {input.r2} \
                -S {output.sam} \
                -p {threads} \
                --phred{params.phred} \
                {params.strand} \
                --dta \
                --new-summary \
                2> {output.stat}
            """

    # SAM → sorted BAM
    rule samtools_sort_index:
        """SAM → sorted BAM + index"""
        input:
            sam = f"{TMP_DIR}/03.align/{{sample}}/{{sample}}.sam",
        output:
            bam      = f"{RESULT_DIR}/03.align/{{sample}}/{{sample}}.sorted.bam",
            bai      = f"{RESULT_DIR}/03.align/{{sample}}/{{sample}}.sorted.bam.bai",
            flagstat = f"{RESULT_DIR}/03.align/{{sample}}/{{sample}}.flagstat.txt",
        params:
            outdir = f"{RESULT_DIR}/03.align/{{sample}}",
        threads: 4
        conda:
            f"{PIPELINE_DIR}/envs/align.yaml"
        shell:
            """
            mkdir -p {params.outdir}
            samtools sort -@ {threads} -o {output.bam} {input.sam}
            samtools index {output.bam}
            samtools flagstat {output.bam} > {output.flagstat}
            """


# ---- STAR 比对 ----
if config["alignment"]["tool"] == "star":

    rule star_align:
        """STAR 双端比对"""
        input:
            r1    = rules.fastp_trim.output.r1_clean,
            r2    = rules.fastp_trim.output.r2_clean,
            index = _get_star_index_input,
        output:
            bam  = f"{TMP_DIR}/03.align/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam",
            stat = f"{RESULT_DIR}/03.align/{{sample}}/{{sample}}.align_stat.txt",
        params:
            idx_dir = _get_align_index_dir(),
            outdir  = f"{RESULT_DIR}/03.align/{{sample}}",
        threads: 6
        conda:
            f"{PIPELINE_DIR}/envs/align.yaml"
        shell:
            """
            mkdir -p {params.outdir} {TMP_DIR}/03.align/{wildcards.sample}
            echo "[INFO] Using STAR index: {params.idx_dir}"
            STAR \\
                --genomeDir {params.idx_dir} \\
                --readFilesIn {input.r1} {input.r2} \\
                --readFilesCommand zcat \\
                --runThreadN {threads} \\
                --outFileNamePrefix {TMP_DIR}/03.align/{wildcards.sample}/{wildcards.sample}. \\
                --outSAMtype BAM SortedByCoordinate \\
                --outSAMunmapped Within \\
                --outSAMattributes Standard \\
                --quantMode TranscriptomeSAM GeneCounts \\
                2> {output.stat}
            """


# ---- 比对率汇总 (通用) ----
rule merge_align_stats:
    input:
        stats = expand(f"{RESULT_DIR}/03.align/{{sample}}/{{sample}}.align_stat.txt", sample=SAMPLES),
    output:
        merged = f"{RESULT_DIR}/03.align/align_stat_all.txt",
    conda:
        f"{PIPELINE_DIR}/envs/qc.yaml"
    run:
        import os
        T = chr(9)   # tab
        N = chr(10)  # newline
        with open(output.merged, "w") as out_fh:
            header = T.join(["sample", "total_reads", "mapped_reads",
                             "uniqu_mapped", "multi_mapped", "mapping_rate"])
            out_fh.write(header + N)
            for stat_file in input.stats:
                sample = os.path.basename(stat_file).replace(".align_stat.txt", "")
                total = mapped = multi = rate = "0"
                with open(stat_file) as in_fh:
                    for line in in_fh:
                        if "Total pairs" in line:
                            total = line.strip().split()[-1]
                        elif "Aligned concordantly exactly 1 time" in line:
                            mapped = line.strip().split()[-1].strip("()%")
                        elif "Aligned concordantly >1 times" in line:
                            multi = line.strip().split()[-1].strip("()%")
                        elif "Overall alignment rate" in line:
                            rate = line.strip().split()[0].rstrip("%")
                row = T.join([sample, total, mapped, mapped, multi, rate])
                out_fh.write(row + N)
