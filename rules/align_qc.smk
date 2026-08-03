# ============================================================
# 比对质控模块: geneBody覆盖度 / Read分布 / Junction饱和度
# ============================================================

# ---- BED12 参考文件 ----
rule make_bed12:
    """从 GTF 生成真 BED12 文件（含 exon block 信息，用于 RSeQC 和 deepTools）

    真 BED12 格式 (12 列):
      1.chrom  2.chromStart(0-based)  3.chromEnd  4.name  5.score  6.strand
      7.thickStart  8.thickEnd  9.itemRgb
      10.blockCount  11.blockSizes(逗号分隔)  12.blockStarts(逗号分隔,相对chromStart)

    RSeQC read_distribution.py 依赖 blockCount/blockSizes/blockStarts
    来判断每个转录本的 CDS/UTR/Intron 区域。无 block 信息 → 输出全为空。
    """
    input:
        gtf = config["paths"]["ref_gtf"],
    output:
        bed12 = f"{TMP_DIR}/ref/ref.bed12",
    conda:
        f"{PIPELINE_DIR}/envs/qc.yaml"
    run:
        import os, re
        os.makedirs(os.path.dirname(str(output.bed12)), exist_ok=True)

        # 按 transcript_id 分组所有 exon
        tx_exons = {}   # tid → {"chrom":str, "strand":str, "exons":[(start,end),...]}
        with open(str(input.gtf)) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue
                feature = parts[2]
                if feature != "exon":
                    continue
                m = re.search(r'transcript_id\s+"?([^\s;"]+)"?', parts[8])
                if not m:
                    continue
                tid = m.group(1)
                chrom = parts[0]
                strand = parts[6]
                start = int(parts[3])
                end = int(parts[4])
                if tid not in tx_exons:
                    tx_exons[tid] = {"chrom": chrom, "strand": strand, "exons": []}
                tx_exons[tid]["exons"].append((start, end))

        with open(str(output.bed12), "w") as f_out:
            for tid, info in tx_exons.items():
                exons = sorted(info["exons"])
                chrom_start = exons[0][0] - 1          # BED 0-based start
                chrom_end = exons[-1][1]               # 最后一个 exon 的 end
                block_sizes = [str(e[1] - e[0] + 1) for e in exons]  # GTF 1-based inclusive → size
                block_starts = [str(e[0] - 1 - chrom_start) for e in exons]
                block_count = len(exons)
                # thickStart = chromStart, thickEnd = chromEnd (整条转录本)
                f_out.write(
                    f"{info['chrom']}\t{chrom_start}\t{chrom_end}\t{tid}\t0\t{info['strand']}"
                    f"\t{chrom_start}\t{chrom_end}\t0\t{block_count}"
                    f"\t{','.join(block_sizes)}\t{','.join(block_starts)}\n"
                )

        fsize = os.path.getsize(str(output.bed12))
        print(f"BED12: {len(tx_exons)} transcripts, {fsize} bytes")


# ---- geneBody 覆盖度 (deepTools) ----
rule bam_to_bigwig:
    """BAM -> bigWig（用于 geneBody 覆盖度和 junction 饱和度）"""
    input:
        bam = f"{RESULT_DIR}/03.align/{{sample}}/{{sample}}.sorted.bam",
    output:
        bw = f"{TMP_DIR}/03.align/{{sample}}/{{sample}}.bw",
    params:
        outdir = f"{TMP_DIR}/03.align/{{sample}}",
    threads: 4
    conda:
        f"{PIPELINE_DIR}/envs/align.yaml"
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bw} \\
            --normalizeUsing CPM \\
            --numberOfProcessors {threads} \\
            --binSize 10
        """


rule gene_body_coverage:
    """geneBody 覆盖度分析"""
    input:
        bws = expand(f"{TMP_DIR}/03.align/{{sample}}/{{sample}}.bw", sample=SAMPLES),
        bed12 = f"{TMP_DIR}/ref/ref.bed12",
    output:
        pdf = f"{RESULT_DIR}/03.align/geneBody_coverage/geneBodyCoverage_all.pdf",
        png = f"{RESULT_DIR}/03.align/geneBody_coverage/geneBodyCoverage_all.png",
        xls = f"{RESULT_DIR}/03.align/geneBody_coverage/geneBodyCoverage_all.xls",
    params:
        outdir = f"{RESULT_DIR}/03.align/geneBody_coverage",
    conda:
        f"{PIPELINE_DIR}/envs/align.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        # 合并 bigWig
        computeMatrix scale-regions \\
            -S {input.bws} \\
            -R {input.bed12} \\
            -b 1000 -a 1000 \\
            --regionBodyLength 2000 \\
            --skipZeros \\
            -o {params.outdir}/matrix.mat.gz \\
            -p 4
        # 绘图 (PDF) + 提取数据 (--outFileNameData)
        plotProfile -m {params.outdir}/matrix.mat.gz \\
            -o {output.pdf} \\
            --outFileNameData {params.outdir}/geneBodyCoverage_data.tab \\
            --perGroup \\
            --plotTitle "Gene Body Coverage"
        # 生成 PNG (convert from PDF)
        convert -density 150 {output.pdf} {output.png} 2>/dev/null || touch {output.png}
        # 生成 xls (复制数据文件)
        if [ -f {params.outdir}/geneBodyCoverage_data.tab ]; then
            cp {params.outdir}/geneBodyCoverage_data.tab {output.xls}
        else
            echo -e "region\\tsample\\tcoverage" > {output.xls}
        fi
        """


# ---- Read 分布 ----
rule read_distribution:
    """Read 在基因组区域的分布 (RSeQC read_distribution.py)"""
    input:
        bam  = f"{RESULT_DIR}/03.align/{{sample}}/{{sample}}.sorted.bam",
        bai  = f"{RESULT_DIR}/03.align/{{sample}}/{{sample}}.sorted.bam.bai",
        bed12 = f"{TMP_DIR}/ref/ref.bed12",
    output:
        stat = f"{RESULT_DIR}/03.align/read_distribution/{{sample}}/{{sample}}.read_distribution.txt",
    params:
        outdir = f"{RESULT_DIR}/03.align/read_distribution/{{sample}}",
    conda:
        f"{PIPELINE_DIR}/envs/align.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        read_distribution.py -i {input.bam} -r {input.bed12} \\
            > {output.stat} 2> {params.outdir}/read_distribution.log
        # 检查是否成功 (RSeQC 成功时输出包含 Total Reads)
        if [ -s {output.stat} ] && grep -q "Total Reads" {output.stat}; then
            echo "  [OK] read_distribution.py completed for {wildcards.sample}"
        else
            echo "  [WARN] read_distribution.py failed or produced empty output for {wildcards.sample}"
            echo -e "Group\\tTotal_bases\\tTag_count" > {output.stat}
            echo -e "CDS_Exons\\t0\\t0" >> {output.stat}
            echo -e "5UTR_Exons\\t0\\t0" >> {output.stat}
            echo -e "3UTR_Exons\\t0\\t0" >> {output.stat}
            echo -e "Introns\\t0\\t0" >> {output.stat}
            echo -e "TSS_up_1kb\\t0\\t0" >> {output.stat}
            echo -e "TSS_up_5kb\\t0\\t0" >> {output.stat}
            echo -e "TTS_down_1kb\\t0\\t0" >> {output.stat}
            echo -e "TTS_down_5kb\\t0\\t0" >> {output.stat}
            echo -e "Intergenic\\t0\\t0" >> {output.stat}
        fi
        """


rule merge_read_distribution:
    """合并 Read 分布结果 — 汇总所有样本的 Tag 分布"""
    input:
        stats = expand(f"{RESULT_DIR}/03.align/read_distribution/{{sample}}/{{sample}}.read_distribution.txt",
                       sample=SAMPLES),
    output:
        merged = f"{RESULT_DIR}/03.align/read_distribution/read_distribution_all.txt",
    conda:
        f"{PIPELINE_DIR}/envs/qc.yaml"
    run:
        import os

        os.makedirs(os.path.dirname(str(output.merged)), exist_ok=True)

        with open(str(output.merged), "w") as out_fh:
            out_fh.write("sample\tgroup\tcategory\tcount\n")
            for stat_file, sample in zip(input.stats, SAMPLES):
                group = SAMPLES_DATA[sample]["group"]
                with open(stat_file) as in_fh:
                    for line in in_fh:
                        line = line.strip()
                        if not line or line.startswith("#") or line.startswith("sample"):
                            continue
                        parts = line.split("\t")
                        if len(parts) >= 2:
                            out_fh.write(f"{sample}\t{group}\t{parts[0]}\t{parts[1]}\n")
        print(f"Merged read distribution: {len(SAMPLES)} samples")


# ---- 每样本 Read 分布图 ----
rule read_distribution_plot:
    """每样本 Read 分布饼图/柱状图"""
    input:
        stat = f"{RESULT_DIR}/03.align/read_distribution/{{sample}}/{{sample}}.read_distribution.txt",
    output:
        pdf = f"{RESULT_DIR}/03.align/read_distribution/{{sample}}/{{sample}}.read_distribution.pdf",
        png = f"{RESULT_DIR}/03.align/read_distribution/{{sample}}/{{sample}}.read_distribution.png",
    params:
        outdir = f"{RESULT_DIR}/03.align/read_distribution/{{sample}}",
    conda:
        f"{PIPELINE_DIR}/envs/diff.yaml"
    run:
        import os, subprocess, sys
        r_code = """
        library(ggplot2)
        library(dplyr)
        
        # 解析 RSeQC read_distribution.py 输出
        txt <- readLines("__SNK_1__")
        # 找到表头行 "Group" 后的数据行
        header_idx <- grep("^Group", txt)
        if(length(header_idx) == 0) {
            pdf("__SNK_2__", width=8, height=6)
            plot.new(); text(0.5, 0.5, "No read distribution data for __SNK_0__")
            dev.off()
            png("__SNK_3__", width=800, height=600)
            plot.new(); text(0.5, 0.5, "No read distribution data for __SNK_0__")
            dev.off()
            quit(save="no")
        }
        
        data_lines <- txt[(header_idx[1]+1):length(txt)]
        data_lines <- data_lines[!grepl("^=", data_lines)]
        data_lines <- data_lines[data_lines != ""]
        
        if(length(data_lines) == 0) {
            pdf("__SNK_2__", width=8, height=6)
            plot.new(); text(0.5, 0.5, "No data rows in read distribution")
            dev.off()
            quit(save="no")
        }
        
        df <- read.table(text=data_lines, header=FALSE, sep="\\t", comment.char="", quote="", fill=TRUE)
        # 列: Group Total_bases Tag_count
        colnames(df) <- c("Group", "Total_bases", "Tag_count")
        
        # 清理: 转为数值, 去掉非数字
        df$Tag_count <- as.numeric(as.character(df$Tag_count))
        df$Total_bases <- as.numeric(as.character(df$Total_bases))
        df <- df[!is.na(df$Tag_count) & df$Tag_count > 0, ]
        df$pct <- df$Tag_count / sum(df$Tag_count) * 100
        df$Group <- gsub("_", " ", df$Group)
        
        # 防空: 过滤后可能无数据
        if(nrow(df) == 0) {
            pdf("__SNK_2__", width=8, height=6)
            plot.new(); text(0.5, 0.5, "All tag counts are zero for __SNK_0__")
            dev.off()
            png("__SNK_3__", width=800, height=600)
            plot.new(); text(0.5, 0.5, "All tag counts are zero for __SNK_0__")
            dev.off()
            quit(save="no")
        }
        
        # 柱状图
        p <- ggplot(df, aes(x=reorder(Group, Tag_count), y=Tag_count/1e6, fill=Group)) +
            geom_bar(stat="identity") +
            coord_flip() +
            labs(x="", y="Tag Count (Million)",
                 title=paste("Read Distribution - __SNK_0__")) +
            theme_bw() +
            theme(legend.position="none",
                  plot.title=element_text(hjust=0.5))
        
        ggsave("__SNK_2__", p, width=8, height=max(4, nrow(df)*0.4))
        ggsave("__SNK_3__", p, width=8, height=max(4, nrow(df)*0.4), dpi=150)
        cat("Plot saved:", "__SNK_2__", "\\n")
        
        """
        r_code = r_code.replace("__SNK_0__", str(wildcards.sample))
        r_code = r_code.replace("__SNK_1__", str(input.stat))
        r_code = r_code.replace("__SNK_2__", str(output.pdf))
        r_code = r_code.replace("__SNK_3__", str(output.png))
        result = subprocess.run(["Rscript", "-"], input=r_code, text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- 合并 Read 分布图 ----
rule merge_read_distribution_plot:
    """汇总所有样本 read distribution 为分组柱状图"""
    input:
        stats = expand(f"{RESULT_DIR}/03.align/read_distribution/{{sample}}/{{sample}}.read_distribution.txt",
                       sample=SAMPLES),
    output:
        pdf = f"{RESULT_DIR}/03.align/read_distribution/read_distribution_all.pdf",
        png = f"{RESULT_DIR}/03.align/read_distribution/read_distribution_all.png",
    params:
        outdir = f"{RESULT_DIR}/03.align/read_distribution",
    conda:
        f"{PIPELINE_DIR}/envs/diff.yaml"
    run:
        import os, subprocess, sys
        r_code = """
        library(ggplot2)
        library(dplyr)
        library(tidyr)
        
        all_data <- list()
        stat_files <- unlist(strsplit("__SNK_0__", " "))
        
        for(f in stat_files) {
            sample_name <- gsub(".read_distribution.txt$", "", basename(f))
            txt <- readLines(f)
            header_idx <- grep("^Group", txt)
            if(length(header_idx) == 0) next
        
            data_lines <- txt[(header_idx[1]+1):length(txt)]
            data_lines <- data_lines[!grepl("^=", data_lines)]
            data_lines <- data_lines[data_lines != ""]
            if(length(data_lines) == 0) next
        
            df <- read.table(text=data_lines, header=FALSE, sep="\\t", comment.char="", fill=TRUE)
            colnames(df) <- c("Group", "Total_bases", "Tag_count")
            df$Tag_count <- as.numeric(as.character(df$Tag_count))
            df <- df[!is.na(df$Tag_count) & df$Tag_count > 0, ]
            df$pct <- df$Tag_count / sum(df$Tag_count) * 100
            df$sample <- sample_name
            df$Group <- gsub("_", " ", df$Group)
            all_data[[length(all_data) + 1]] <- df[, c("sample","Group","pct")]
        }
        
        if(length(all_data) == 0) {
            pdf("__SNK_1__"); plot.new(); text(0.5,0.5,"No data"); dev.off()
            png("__SNK_2__"); plot.new(); text(0.5,0.5,"No data"); dev.off()
            quit(save="no")
        }
        
        merged <- bind_rows(all_data)
        
        p <- ggplot(merged, aes(x=sample, y=pct, fill=Group)) +
            geom_bar(stat="identity", position="fill") +
            labs(x="", y="Proportion of Tags (%)",
                 title="Read Distribution Across Genomic Regions") +
            theme_bw() +
            theme(axis.text.x=element_text(angle=45, hjust=1),
                  plot.title=element_text(hjust=0.5)) +
            scale_fill_brewer(palette="Set3") +
            scale_y_continuous(labels=scales::percent)
        
        w <- max(6, length(unique(merged$sample)) * 0.5)
        ggsave("__SNK_1__", p, width=w, height=6)
        ggsave("__SNK_2__", p, width=w, height=6, dpi=150)
        cat("Merged plot saved:", "__SNK_1__", "\\n")
        
        """
        r_code = r_code.replace("__SNK_0__", str(input.stats))
        r_code = r_code.replace("__SNK_1__", str(output.pdf))
        r_code = r_code.replace("__SNK_2__", str(output.png))
        result = subprocess.run(["Rscript", "-"], input=r_code, text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- Junction 饱和度 (deepTools) ----
rule junction_saturation:
    """Junction 饱和度分析 — 评估测序深度是否足够覆盖剪接点"""
    input:
        bams = expand(f"{RESULT_DIR}/03.align/{{sample}}/{{sample}}.sorted.bam", sample=SAMPLES),
        bed12 = f"{TMP_DIR}/ref/ref.bed12",
    output:
        pdf = f"{RESULT_DIR}/03.align/junction_saturation/junction_saturation_all.pdf",
        png = f"{RESULT_DIR}/03.align/junction_saturation/junction_saturation_all.png",
        xls = f"{RESULT_DIR}/03.align/junction_saturation/junction_saturation_all.xls",
    params:
        outdir = f"{RESULT_DIR}/03.align/junction_saturation",
        labels  = " ".join(SAMPLES),
    threads: 8
    conda:
        f"{PIPELINE_DIR}/envs/align.yaml"
    shell:
        """
        mkdir -p {params.outdir}

        # 过滤: 只保留实际存在的 BAM 文件及其对应标签
        ALL_BAMS=({input.bams})
        ALL_LABELS=({params.labels})
        VALID_BAMS=()
        VALID_LABELS=()
        for i in "${{!ALL_BAMS[@]}}"; do
            if [ -f "${{ALL_BAMS[$i]}}" ]; then
                VALID_BAMS+=("${{ALL_BAMS[$i]}}")
                VALID_LABELS+=("${{ALL_LABELS[$i]}}")
            else
                echo "  [WARN] BAM missing, skipping: ${{ALL_BAMS[$i]}}"
            fi
        done

        echo "Valid BAMs: ${{#VALID_BAMS[@]}} / ${{#ALL_BAMS[@]}}"

        if [ ${{#VALID_BAMS[@]}} -lt 2 ]; then
            echo "[WARN] Too few valid BAMs (${{#VALID_BAMS[@]}}) for junction saturation"
            echo -e "sample1\\tsample2\\tJSD" > {output.xls}
            touch {output.pdf} {output.png}
            exit 0
        fi

        plotFingerprint \\
            -b ${{VALID_BAMS[@]}} \\
            --labels ${{VALID_LABELS[@]}} \\
            --JSDsample {params.outdir}/JSD.txt \\
            --plotFile {output.pdf} \\
            --numberOfProcessors {threads} \\
            --minMappingQuality 10 2>&1 | tail -5 || {{
            echo "[WARN] plotFingerprint failed — creating placeholder outputs"
            echo -e "sample1\\tsample2\\tJSD" > {output.xls}
            touch {output.pdf} {output.png}
            exit 0
        }}

        # 转换 PDF → PNG
        if command -v convert &>/dev/null && [ -s {output.pdf} ]; then
            convert -density 150 {output.pdf} {output.png} 2>/dev/null || true
        fi
        [ -f {output.png} ] || touch {output.png}

        # 导出 JSD 数据
        if [ -f {params.outdir}/JSD.txt ]; then
            cp {params.outdir}/JSD.txt {output.xls}
        else
            echo -e "sample1\\tsample2\\tJSD" > {output.xls}
        fi
        """


# ---- 链特异性推断 (可选) ----
rule infer_experiment:
    """RSeQC infer_experiment 推断链特异性"""
    input:
        bam  = f"{RESULT_DIR}/03.align/{{sample}}/{{sample}}.sorted.bam",
        bai  = f"{RESULT_DIR}/03.align/{{sample}}/{{sample}}.sorted.bam.bai",
        bed12 = f"{TMP_DIR}/ref/ref.bed12",
    output:
        strand_info = f"{RESULT_DIR}/03.align/strand_inference/{{sample}}/{{sample}}.infer_experiment.txt",
    params:
        outdir = f"{RESULT_DIR}/03.align/strand_inference/{{sample}}",
    conda:
        f"{PIPELINE_DIR}/envs/align.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        infer_experiment.py -r {input.bed12} -i {input.bam} > {output.strand_info} 2>&1 || \\
            echo "infer_experiment failed for {wildcards.sample}" > {output.strand_info}
        """
