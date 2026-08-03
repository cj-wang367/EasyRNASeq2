# ============================================================
# 表达定量模块: featureCounts / StringTie / RSEM
# ============================================================

# ---- gene 长度提取 ----
rule get_gene_length:
    """从 GTF 提取 gene 长度 (合并重叠外显子)"""
    input:
        gtf = config["paths"]["ref_gtf"],
    output:
        gene_len = f"{TMP_DIR}/ref/gene_length.txt",
    conda:
        f"{PIPELINE_DIR}/envs/expression.yaml"
    run:
        import os, re

        os.makedirs(os.path.dirname(str(output.gene_len)), exist_ok=True)

        genes = {}
        with open(str(input.gtf)) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                if parts[2] != 'exon':
                    continue
                start = int(parts[3])
                end = int(parts[4])
                # 提取 gene_id
                m = re.search(r'gene_id\s+"?([^\s;"]+)"?', parts[8])
                if m:
                    gid = m.group(1)
                    if gid not in genes:
                        genes[gid] = set()
                    genes[gid].add((start, end))

        with open(str(output.gene_len), 'w') as o:
            o.write('gene_id\tlength\n')
            for gid, exons in genes.items():
                sorted_exons = sorted(exons)
                merged_len = 0
                c_start, c_end = sorted_exons[0]
                for s, e in sorted_exons[1:]:
                    if s <= c_end:
                        c_end = max(c_end, e)
                    else:
                        merged_len += c_end - c_start + 1
                        c_start, c_end = s, e
                merged_len += c_end - c_start + 1
                o.write(f'{gid}\t{merged_len}\n')


# ---- featureCounts 定量 ----
rule featurecounts:
    """featureCounts 基因和转录本定量"""
    input:
        bam      = f"{RESULT_DIR}/03.align/{{sample}}/{{sample}}.sorted.bam",
        gtf      = config["paths"]["ref_gtf"],
    output:
        gene_count   = f"{TMP_DIR}/05.expression/{{sample}}/{{sample}}.gene_count.txt",
        tran_count   = f"{TMP_DIR}/05.expression/{{sample}}/{{sample}}.transcript_count.txt",
    params:
        outdir = f"{TMP_DIR}/05.expression/{{sample}}",
        strand = {
            "unstranded": 0,
            "firststrand": 2,
            "secondstrand": 1,
        }[config["alignment"]["strand"]],
    threads: 4
    conda:
        f"{PIPELINE_DIR}/envs/expression.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        # 基因水平
        featureCounts -a {input.gtf} \\
            -o {params.outdir}/{wildcards.sample}.gene.xls \\
            -t exon -g gene_id \\
            -s {params.strand} \\
            -p --primary \\
            -T {threads} \\
            {input.bam} 2>&1 | tail -5

        # 提取 count 列
        tail -n +3 {params.outdir}/{wildcards.sample}.gene.xls | \\
            awk -F'\\t' '{{print $1"\\t"$NF}}' > {output.gene_count}

        # 转录本水平
        featureCounts -a {input.gtf} \\
            -o {params.outdir}/{wildcards.sample}.transcript.xls \\
            -t exon -g transcript_id \\
            -s {params.strand} \\
            -p --primary \\
            -T {threads} \\
            {input.bam} 2>&1 | tail -5

        tail -n +3 {params.outdir}/{wildcards.sample}.transcript.xls | \\
            awk -F'\\t' '{{print $1"\\t"$NF}}' > {output.tran_count}
        """


# ---- StringTie 定量 ----
rule stringtie_quant:
    """StringTie 表达定量"""
    input:
        bam = f"{RESULT_DIR}/03.align/{{sample}}/{{sample}}.sorted.bam",
        gtf = config["paths"]["ref_gtf"],
    output:
        gtf_out     = f"{TMP_DIR}/05.expression/{{sample}}/{{sample}}.stringtie.gtf",
        gene_abund  = f"{TMP_DIR}/05.expression/{{sample}}/{{sample}}.gene_abund.tab",
    params:
        outdir = f"{TMP_DIR}/05.expression/{{sample}}",
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
            -e -B \\
            -p {threads} \\
            {params.strand} \\
            -o {output.gtf_out} \\
            -A {output.gene_abund}
        """


# ---- 合并表达矩阵 ----
rule merge_gene_counts:
    """合并各样本 gene count 为矩阵"""
    input:
        counts = expand(f"{TMP_DIR}/05.expression/{{sample}}/{{sample}}.gene_count.txt", sample=SAMPLES),
    output:
        count_mat = f"{RESULT_DIR}/05.expression/gene/gene_count.xls",
    conda:
        f"{PIPELINE_DIR}/envs/expression.yaml"
    shell:
        "python3 {PIPELINE_DIR}/scripts/merge_counts.py {TMP_DIR}/05.expression {output.count_mat} gene"


# ---- 合并转录本表达矩阵 ----
rule merge_transcript_counts:
    """合并各样本 transcript count 为矩阵"""
    input:
        counts = expand(f"{TMP_DIR}/05.expression/{{sample}}/{{sample}}.transcript_count.txt", sample=SAMPLES),
    output:
        count_mat = f"{RESULT_DIR}/05.expression/transcript/transcript_count.xls",
    conda:
        f"{PIPELINE_DIR}/envs/expression.yaml"
    shell:
        "python3 {PIPELINE_DIR}/scripts/merge_counts.py {TMP_DIR}/05.expression {output.count_mat} transcript"


# ---- FPKM/TPM 标准化 (gene 级) ----
rule normalize_expression:
    """gene count -> FPKM + TPM 标准化"""
    input:
        count_mat  = rules.merge_gene_counts.output.count_mat,
        gene_len   = rules.get_gene_length.output.gene_len,
    output:
        norm_mat   = f"{RESULT_DIR}/05.expression/gene/gene_norm.xls",
        tpm_mat    = f"{RESULT_DIR}/05.expression/gene/gene_tpm.xls",
        fpkm_mat   = f"{RESULT_DIR}/05.expression/gene/gene_fpkm.xls",
    params:
        norm_type  = config['expression']['count_norm'],
    conda:
        f"{PIPELINE_DIR}/envs/expression.yaml"
    run:
        import os, subprocess, sys
        r_code = """
        import pandas as pd
        import numpy as np
        import shutil
        
        # 读取 count 矩阵
        counts = pd.read_csv('__SNK_2__', sep='\\t', index_col=0)
        lengths = pd.read_csv('__SNK_4__', sep='\\t', index_col=0)
        
        # 找到共同的基因
        common_genes = counts.index.intersection(lengths.index)
        counts = counts.loc[common_genes]
        lengths = lengths.loc[common_genes]
        
        # FPKM
        total_reads = counts.sum(axis=0)
        fpkm = counts.div(lengths['length'], axis=0) * 1e3
        fpkm = fpkm.div(total_reads, axis=1) * 1e6
        
        # TPM (从 FPKM)
        tpm = fpkm.div(fpkm.sum(axis=0), axis=1) * 1e6
        
        fpkm.to_csv('__SNK_3__', sep='\\t')
        tpm.to_csv('__SNK_5__', sep='\\t')
        
        # 根据配置选择标准化矩阵作为 norm_mat
        norm_type = '__SNK_0__'
        if norm_type == 'fpkm':
            shutil.copy('__SNK_3__', '__SNK_1__')
        else:
            shutil.copy('__SNK_5__', '__SNK_1__')
        print(f'Normalization: {norm_type}, files saved.')
        
        """
        r_code = r_code.replace("__SNK_0__", str(params.norm_type))
        r_code = r_code.replace("__SNK_1__", str(output.norm_mat))
        r_code = r_code.replace("__SNK_2__", str(input.count_mat))
        r_code = r_code.replace("__SNK_3__", str(output.fpkm_mat))
        r_code = r_code.replace("__SNK_4__", str(input.gene_len))
        r_code = r_code.replace("__SNK_5__", str(output.tpm_mat))
        result = subprocess.run(["python3", "-c", r_code], text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- FPKM/TPM 标准化 (transcript 级) ----
rule normalize_transcript_expression:
    """transcript count -> CPM 标准化 (无长度信息时的降级方案)"""
    input:
        count_mat = rules.merge_transcript_counts.output.count_mat,
    output:
        norm_mat  = f"{RESULT_DIR}/05.expression/transcript/transcript_norm.xls",
        cpm_mat   = f"{RESULT_DIR}/05.expression/transcript/transcript_cpm.xls",
        fpkm_mat  = f"{RESULT_DIR}/05.expression/transcript/transcript_fpkm.xls",
    params:
        norm_type = config['expression']['count_norm'],
    conda:
        f"{PIPELINE_DIR}/envs/expression.yaml"
    run:
        import os, subprocess, sys
        r_code = """
        import pandas as pd
        import numpy as np
        import shutil
        
        # 读取 transcript count 矩阵
        counts = pd.read_csv('__SNK_2__', sep='\\t', index_col=0)
        
        # CPM (Counts Per Million): 无转录本长度信息时的降级标准化
        # 标准 CPM = (count / total_count) * 1e6
        cpm = counts.div(counts.sum(axis=0), axis=1) * 1e6
        cpm = cpm.round(2)
        cpm.to_csv('__SNK_4__', sep='\\t')
        
        # FPKM 对于 transcript 级: 需要 transcript 长度, 如无则用 CPM 代替
        fpkm = cpm.copy()
        fpkm.to_csv('__SNK_3__', sep='\\t')
        
        # norm_mat: 根据配置选择
        norm_type = '__SNK_0__'
        if norm_type == 'fpkm':
            shutil.copy('__SNK_3__', '__SNK_1__')
        else:
            shutil.copy('__SNK_4__', '__SNK_1__')
        print(f'Transcript normalization (CPM fallback): {norm_type}, {counts.shape[0]} transcripts x {counts.shape[1]} samples')
        print(f'Note: TPM=CPM here — transcript lengths unavailable; for true TPM, enable transcript-level length extraction')
        
        """
        r_code = r_code.replace("__SNK_0__", str(params.norm_type))
        r_code = r_code.replace("__SNK_1__", str(output.norm_mat))
        r_code = r_code.replace("__SNK_2__", str(input.count_mat))
        r_code = r_code.replace("__SNK_3__", str(output.fpkm_mat))
        r_code = r_code.replace("__SNK_4__", str(output.cpm_mat))
        result = subprocess.run(["python3", "-c", r_code], text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ============================================================
# 表达量可视化模块
# ============================================================

rule expression_box_density:
    """基因和转录本水平表达量箱线图+密度图"""
    input:
        gene_norm   = rules.normalize_expression.output.norm_mat,
        tran_norm   = rules.normalize_transcript_expression.output.norm_mat,
        sample_desc = f"{TMP_DIR}/diff/sample_desc.txt",
    output:
        gene_box     = f"{RESULT_DIR}/05.expression/gene/box_density/boxplot.pdf",
        gene_density = f"{RESULT_DIR}/05.expression/gene/box_density/density.pdf",
        tran_box     = f"{RESULT_DIR}/05.expression/transcript/box_density/boxplot.pdf",
        tran_density = f"{RESULT_DIR}/05.expression/transcript/box_density/density.pdf",
    conda:
        f"{PIPELINE_DIR}/envs/diff.yaml"
    run:
        import os, subprocess, sys
        os.makedirs(f"{RESULT_DIR}/05.expression/gene/box_density {RESULT_DIR}/05.expression/transcript/box_density", exist_ok=True)
        r_code = """
        library(ggplot2)
        library(tidyr)
        library(dplyr)
        samples <- read.table("__SNK_2__", header=TRUE, sep="\\t", comment.char="", check.names=FALSE)
        plot_bd <- function(f, box_f, dens_f, label) {
            expr <- read.table(f, header=TRUE, row.names=1, sep="\\t", check.names=FALSE)
            expr_log <- as.data.frame(log2(expr + 1))
            expr_long <- expr_log %>% tibble::rownames_to_column("gene") %>% pivot_longer(-gene, names_to="sample", values_to="expr")
            expr_long <- merge(expr_long, samples, by="sample", all.x=TRUE)
            p1 <- ggplot(expr_long, aes(x=sample, y=expr, fill=group)) + geom_boxplot(outlier.size=0.3, outlier.alpha=0.3) +
                labs(x="", y="log2(Expression + 1)", title=paste0(label, " - Expression Boxplot")) +
                theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="bottom", plot.title=element_text(hjust=0.5))
            ggsave(box_f, p1, width=max(8, length(unique(expr_long$sample))*0.8), height=6)
            p2 <- ggplot(expr_long, aes(x=expr, color=sample, group=sample)) + geom_density(size=0.8) +
                labs(x="log2(Expression + 1)", y="Density", title=paste0(label, " - Expression Density")) +
                theme_bw() + theme(legend.position="bottom", plot.title=element_text(hjust=0.5))
            ggsave(dens_f, p2, width=10, height=6)
        }
        plot_bd("__SNK_4__", "__SNK_5__", "__SNK_1__", "Gene")
        plot_bd("__SNK_3__", "__SNK_6__", "__SNK_0__", "Transcript")
        
        """
        r_code = r_code.replace("__SNK_0__", str(output.tran_density))
        r_code = r_code.replace("__SNK_1__", str(output.gene_density))
        r_code = r_code.replace("__SNK_2__", str(input.sample_desc))
        r_code = r_code.replace("__SNK_3__", str(input.tran_norm))
        r_code = r_code.replace("__SNK_4__", str(input.gene_norm))
        r_code = r_code.replace("__SNK_5__", str(output.gene_box))
        r_code = r_code.replace("__SNK_6__", str(output.tran_box))
        r_code = r_code.replace("__SNK_7__", str(RESULT_DIR))
        result = subprocess.run(["Rscript", "-"], input=r_code, text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


rule expression_pca:
    """基因水平 PCA 分析 (含坐标+方差)"""
    input:
        gene_norm   = rules.normalize_expression.output.norm_mat,
        sample_desc = f"{TMP_DIR}/diff/sample_desc.txt",
    output:
        pca_pdf   = f"{RESULT_DIR}/05.expression/gene/pca/PCA_2D.pdf",
        pca_png   = f"{RESULT_DIR}/05.expression/gene/pca/PCA_2D.png",
        pca_coord = f"{RESULT_DIR}/05.expression/gene/pca/PCA_coordinates.xls",
        pca_prop  = f"{RESULT_DIR}/05.expression/gene/pca/PCA_proportion.xls",
        pca_var   = f"{RESULT_DIR}/05.expression/gene/pca/PCA_variance.pdf",
    conda:
        f"{PIPELINE_DIR}/envs/diff.yaml"
    run:
        import os, subprocess, sys
        os.makedirs(f"{RESULT_DIR}/05.expression/gene/pca", exist_ok=True)
        r_code = """
        library(ggplot2); library(ggrepel)
        samples <- read.table("__SNK_0__", header=TRUE, sep="\\t", comment.char="", check.names=FALSE)
        expr <- read.table("__SNK_3__", header=TRUE, row.names=1, sep="\\t", check.names=FALSE)
        expr_log <- log2(expr + 1)
        row_vars <- apply(expr_log, 1, var); expr_log <- expr_log[row_vars > 0, ]
        pca <- prcomp(t(expr_log), scale.=TRUE)
        pca_df <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], sample=colnames(expr_log))
        pca_df <- merge(pca_df, samples, by="sample")
        ve <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
        p <- ggplot(pca_df, aes(x=PC1, y=PC2, color=group, label=sample)) + geom_point(size=3) +
            geom_text_repel(size=3, show.legend=FALSE) + stat_ellipse(level=0.95) +
            labs(x=paste0("PC1 (",ve[1],"%)"), y=paste0("PC2 (",ve[2],"%)"), title="Gene Expression PCA") +
            theme_bw() + theme(plot.title=element_text(hjust=0.5))
        ggsave("__SNK_5__", p, width=8, height=6)
        ggsave("__SNK_6__", p, width=8, height=6, dpi=150)
        write.table(pca_df, "__SNK_1__", sep="\\t", row.names=FALSE, quote=FALSE)
        write.table(data.frame(PC=1:length(pca$sdev), Variance=ve), "__SNK_2__", sep="\\t", row.names=FALSE, quote=FALSE)
        scree <- data.frame(PC=1:min(20,length(pca$sdev)), Variance=ve[1:min(20,length(pca$sdev))])
        p_var <- ggplot(scree, aes(x=factor(PC), y=Variance)) + geom_bar(stat="identity", fill="#3498db") +
            labs(x="PC", y="Variance (%)", title="PCA Variance Explained") + theme_bw()
        ggsave("__SNK_4__", p_var, width=10, height=5)
        
        """
        r_code = r_code.replace("__SNK_0__", str(input.sample_desc))
        r_code = r_code.replace("__SNK_1__", str(output.pca_coord))
        r_code = r_code.replace("__SNK_2__", str(output.pca_prop))
        r_code = r_code.replace("__SNK_3__", str(input.gene_norm))
        r_code = r_code.replace("__SNK_4__", str(output.pca_var))
        r_code = r_code.replace("__SNK_5__", str(output.pca_pdf))
        r_code = r_code.replace("__SNK_6__", str(output.pca_png))
        r_code = r_code.replace("__SNK_7__", str(RESULT_DIR))
        result = subprocess.run(["Rscript", "-"], input=r_code, text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


rule expression_correlation:
    """基因水平样本相关性 (Pearson + Spearman)"""
    input:
        gene_norm   = rules.normalize_expression.output.norm_mat,
        sample_desc = f"{TMP_DIR}/diff/sample_desc.txt",
    output:
        pearson_pdf  = f"{RESULT_DIR}/05.expression/gene/corr/pearson_cor.pdf",
        pearson_xls  = f"{RESULT_DIR}/05.expression/gene/corr/pearson_cor.xls",
        spearman_pdf = f"{RESULT_DIR}/05.expression/gene/corr/spearman_cor.pdf",
        spearman_xls = f"{RESULT_DIR}/05.expression/gene/corr/spearman_cor.xls",
    conda:
        f"{PIPELINE_DIR}/envs/diff.yaml"
    run:
        import os, subprocess, sys
        os.makedirs(f"{RESULT_DIR}/05.expression/gene/corr", exist_ok=True)
        r_code = """
        library(pheatmap)
        samples <- read.table("__SNK_4__", header=TRUE, sep="\\t", comment.char="", check.names=FALSE)
        rownames(samples) <- samples$sample
        expr <- read.table("__SNK_5__", header=TRUE, row.names=1, sep="\\t", check.names=FALSE)
        expr_log <- log2(expr + 1)
        common <- intersect(colnames(expr_log), samples$sample)
        ann_col <- data.frame(Group=samples[common, "group"]); rownames(ann_col) <- common
        for(m in c("pearson","spearman")) {
            cm <- cor(expr_log[, common], method=m)
            pdf(if(m=="pearson") "__SNK_3__" else "__SNK_1__", width=10, height=8)
            pheatmap(cm, annotation_col=ann_col, display_numbers=TRUE, number_format="%.2f", fontsize_number=8,
                main=paste("Sample Correlation (", m, ")"))
            dev.off()
            write.table(cm, if(m=="pearson") "__SNK_2__" else "__SNK_0__", sep="\\t", quote=FALSE)
        }
        
        """
        r_code = r_code.replace("__SNK_0__", str(output.spearman_xls))
        r_code = r_code.replace("__SNK_1__", str(output.spearman_pdf))
        r_code = r_code.replace("__SNK_2__", str(output.pearson_xls))
        r_code = r_code.replace("__SNK_3__", str(output.pearson_pdf))
        r_code = r_code.replace("__SNK_4__", str(input.sample_desc))
        r_code = r_code.replace("__SNK_5__", str(input.gene_norm))
        r_code = r_code.replace("__SNK_6__", str(RESULT_DIR))
        result = subprocess.run(["Rscript", "-"], input=r_code, text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


rule expression_heatmap:
    """Top 5000 高变基因样本聚类热图"""
    input:
        gene_norm   = rules.normalize_expression.output.norm_mat,
        sample_desc = f"{TMP_DIR}/diff/sample_desc.txt",
    output:
        heatmap_pdf  = f"{RESULT_DIR}/05.expression/gene/heatmap/heatmap.pdf",
        heatmap_png  = f"{RESULT_DIR}/05.expression/gene/heatmap/heatmap.png",
        heatmap_data = f"{RESULT_DIR}/05.expression/gene/heatmap/heatmap_plot.xls",
    conda:
        f"{PIPELINE_DIR}/envs/diff.yaml"
    run:
        import os, subprocess, sys
        os.makedirs(f"{RESULT_DIR}/05.expression/gene/heatmap", exist_ok=True)
        r_code = """
        library(pheatmap)
        samples <- read.table("__SNK_3__", header=TRUE, sep="\\t", comment.char="", check.names=FALSE)
        rownames(samples) <- samples$sample
        expr <- read.table("__SNK_4__", header=TRUE, row.names=1, sep="\\t", check.names=FALSE)
        expr_log <- log2(expr + 1)
        rv <- apply(expr_log, 1, var)
        sel <- order(rv, decreasing=TRUE)[1:min(5000, nrow(expr_log))]
        esub <- expr_log[sel, , drop=FALSE]
        ez <- t(scale(t(esub))); ez <- ez[complete.cases(ez), , drop=FALSE]
        common <- intersect(colnames(ez), samples$sample)
        ann <- data.frame(Group=samples[common, "group"]); rownames(ann) <- common
        if(nrow(ez) < 2) {
            pdf("__SNK_1__", width=8, height=6)
            plot.new(); text(0.5, 0.5, paste0("Too few variable genes (n=", nrow(ez), ") for heatmap"))
            dev.off()
            png("__SNK_2__", width=800, height=600)
            plot.new(); text(0.5, 0.5, paste0("Too few variable genes (n=", nrow(ez), ") for heatmap"))
            dev.off()
            write.table(data.frame(Note="Too few variable genes after Z-score"), "__SNK_0__", sep="\\t")
            quit(status=0)
        }
        pdf("__SNK_1__", width=10, height=8)
        pheatmap(ez[,common,drop=FALSE], annotation_col=ann, show_rownames=FALSE, cluster_rows=TRUE,
            cluster_cols=TRUE, scale="none", main=paste0("Sample Clustering (Top ", nrow(ez)," Genes)"))
        dev.off()
        png("__SNK_2__", width=1200, height=960)
        pheatmap(ez[,common,drop=FALSE], annotation_col=ann, show_rownames=FALSE, cluster_rows=TRUE,
            cluster_cols=TRUE, scale="none", main=paste0("Sample Clustering (Top ", nrow(ez)," Genes)"))
        dev.off()
        write.table(ez[,common,drop=FALSE], "__SNK_0__", sep="\\t", quote=FALSE)
        
        """
        r_code = r_code.replace("__SNK_0__", str(output.heatmap_data))
        r_code = r_code.replace("__SNK_1__", str(output.heatmap_pdf))
        r_code = r_code.replace("__SNK_2__", str(output.heatmap_png))
        r_code = r_code.replace("__SNK_3__", str(input.sample_desc))
        r_code = r_code.replace("__SNK_4__", str(input.gene_norm))
        r_code = r_code.replace("__SNK_5__", str(RESULT_DIR))
        result = subprocess.run(["Rscript", "-"], input=r_code, text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")
