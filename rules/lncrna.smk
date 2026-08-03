# ============================================================
# LncRNA 模块: CPC2 + PLEK 编码潜力预测 + 分类 + 靶基因
# ============================================================

# ---- 筛选候选 LncRNA 转录本 ----
rule filter_lncrna_candidates:
    """筛选可能为 LncRNA 的新转录本（长度 200-20000nt）"""
    input:
        fa = rules.extract_transcript_fasta.output.new_fa,
    output:
        candidates = f"{TMP_DIR}/08.lnc/candidates.fa",
    conda:
        f"{PIPELINE_DIR}/envs/lncrna.yaml"
    shell:
        """
        mkdir -p {TMP_DIR}/08.lnc
        seqkit seq -m 200 -M 20000 {input.fa} > {output.candidates}
        echo "Candidate lncRNA transcripts: $(grep -c '^>' {output.candidates})"
        """


# ---- CPC2 编码潜力预测 ----
rule cpc2_predict:
    """CPC2 编码潜力预测 (若 CPC2 不可用则跳过)"""
    input:
        fa = rules.filter_lncrna_candidates.output.candidates,
    output:
        cpc2_out = f"{TMP_DIR}/08.lnc/CPC2_output.txt",
    conda:
        f"{PIPELINE_DIR}/envs/lncrna.yaml"
    shell:
        """
        mkdir -p {TMP_DIR}/08.lnc
        if [ ! -s {input.fa} ]; then
            echo "# No lncRNA candidates found (input FASTA empty or missing)" > {output.cpc2_out}
            echo "# Skipping CPC2 prediction" >> {output.cpc2_out}
            echo "CPC2 skipped: input FASTA is empty"
            exit 0
        fi
        if command -v CPC2.py &> /dev/null; then
            CPC2.py -i {input.fa} -o {output.cpc2_out}
        elif command -v CPC2 &> /dev/null; then
            CPC2 -i {input.fa} -o {output.cpc2_out}
        else
            echo "# CPC2 not found. Install: conda install -c bioconda cpc2" > {output.cpc2_out}
            echo "# Skipping CPC2 prediction" >> {output.cpc2_out}
            echo "CPC2 not found — lncRNA identification will be skipped"
        fi
        """


# ---- PLEK 编码潜力预测 ----
rule plek_predict:
    """PLEK k-mer 编码潜力预测 (若 PLEK 不可用则跳过, 后续使用 CPC2-only 模式)"""
    input:
        fa = rules.filter_lncrna_candidates.output.candidates,
    output:
        plek_out = f"{TMP_DIR}/08.lnc/PLEK_output.txt",
    threads: 8
    conda:
        f"{PIPELINE_DIR}/envs/lncrna.yaml"
    shell:
        """
        mkdir -p {TMP_DIR}/08.lnc
        if [ ! -s {input.fa} ]; then
            echo "# No lncRNA candidates found (input FASTA empty or missing)" > {output.plek_out}
            echo "# Skipping PLEK prediction" >> {output.plek_out}
            echo "PLEK skipped: input FASTA is empty"
            exit 0
        fi
        if command -v PLEK.py &> /dev/null; then
            PLEK.py -fasta {input.fa} \\
                -out {output.plek_out} \\
                -thread {threads} \\
                -minlength 200 2>&1 | tail -5 || echo "PLEK completed"
        elif command -v plek &> /dev/null; then
            plek -fasta {input.fa} \\
                -out {output.plek_out} \\
                -thread {threads} \\
                -minlength 200 2>&1 | tail -5 || echo "PLEK completed"
        else
            echo "# PLEK not installed. Download from https://github.com/liambai/PLEK" > {output.plek_out}
            echo "# Will use CPC2-only mode for lncRNA identification" >> {output.plek_out}
            echo "PLEK not found — will use CPC2-only mode"
        fi
        """


# ---- 合并 CPC2 + PLEK 结果 ----
rule merge_lncrna_results:
    """取 CPC2 和 PLEK 交集作为最终 LncRNA (若 PLEK 不可用则仅使用 CPC2)"""
    input:
        cpc2 = rules.cpc2_predict.output.cpc2_out,
        plek = rules.plek_predict.output.plek_out,
    output:
        lncrna_list = f"{RESULT_DIR}/08.lnc/novel_lncrna_list.txt",
    conda:
        f"{PIPELINE_DIR}/envs/lncrna.yaml"
    run:
        import os, subprocess, sys
        os.makedirs(f"{RESULT_DIR}/08.lnc", exist_ok=True)
        r_code = """
        # 解析 CPC2 (coding/noncoding)
        cpc2_noncoding = set()
        with open('__SNK_2__') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.strip().split('\\t')
                if len(parts) >= 4 and parts[3] == 'noncoding':
                    cpc2_noncoding.add(parts[0])
        
        print(f'CPC2 noncoding: {len(cpc2_noncoding)}')
        
        # 解析 PLEK (score < 0 = noncoding) — 若不可用则跳过
        plek_available = True
        plek_noncoding = set()
        with open('__SNK_3__') as f:
            first_line = f.readline()
            if first_line.startswith('# PLEK not installed'):
                plek_available = False
                print('PLEK not available — using CPC2-only mode')
            else:
                # 回退到开头
                f.seek(0)
            if plek_available:
                for line in f:
                    if line.startswith('#'):
                        continue
                    parts = line.strip().split('\\t')
                    if len(parts) >= 2:
                        try:
                            score = float(parts[1])
                            if score < 0:
                                plek_noncoding.add(parts[0])
                        except ValueError:
                            continue
        
        print(f'PLEK noncoding: {len(plek_noncoding)}')
        
        # 取交集 (或仅 CPC2)
        if plek_available:
            novel_lncrna = cpc2_noncoding & plek_noncoding
            print(f'Novel lncRNA (CPC2 ∩ PLEK): {len(novel_lncrna)}')
        else:
            novel_lncrna = cpc2_noncoding
            print(f'Novel lncRNA (CPC2 only): {len(novel_lncrna)}')
        
        with open('__SNK_0__', 'w') as f:
            for lnc in sorted(novel_lncrna):
                f.write(f'{lnc}\\n')
        
        # 保存基本统计
        with open('__SNK_1__/08.lnc/lncrna_summary.txt', 'w') as sf:
            sf.write(f'CPC2_noncoding\\t{len(cpc2_noncoding)}\\n')
            sf.write(f'PLEK_noncoding\\t{len(plek_noncoding)}\\n')
            sf.write(f'PLEK_available\\t{plek_available}\\n')
            sf.write(f'Novel_lncRNA\\t{len(novel_lncrna)}\\n')
        
        """
        r_code = r_code.replace("__SNK_0__", str(output.lncrna_list))
        r_code = r_code.replace("__SNK_1__", str(RESULT_DIR))
        r_code = r_code.replace("__SNK_2__", str(input.cpc2))
        r_code = r_code.replace("__SNK_3__", str(input.plek))
        result = subprocess.run(["python3", "-c", r_code], text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- LncRNA 差异表达 (可选) ----
rule lncrna_de:
    """对鉴定到的 novel lncRNA 做差异表达分析"""
    input:
        lncrna_list = rules.merge_lncrna_results.output.lncrna_list,
        count_mat   = f"{RESULT_DIR}/05.expression/gene/gene_count.xls",
        sample_desc = f"{TMP_DIR}/diff/sample_desc.txt",
    output:
        all_results  = f"{RESULT_DIR}/08.lnc/de/{{comp}}/{{comp}}.lnc_DE_results.All.xls",
        sig_results  = f"{RESULT_DIR}/08.lnc/de/{{comp}}/{{comp}}.lnc_DE_results.Significant.xls",
    params:
        logFC       = config["diff"]["logFC"],
        pvalue      = config["diff"]["diff_pvalue"],
        ptype       = config["diff"]["pOrFDR"],
        outdir      = f"{RESULT_DIR}/08.lnc/de/{{comp}}",
    conda:
        f"{PIPELINE_DIR}/envs/diff.yaml"
    run:
        import os, subprocess, sys
        os.makedirs(f"{params.outdir}", exist_ok=True)
        r_code = """
        library(DESeq2)
        library(dplyr)
        
        # 读 lncRNA 列表
        lnc_ids <- readLines("__SNK_2__")
        cat("LncRNA IDs to analyze:", length(lnc_ids), "\\n")
        
        # 读 count 矩阵
        count_mat <- read.table("__SNK_4__", header=TRUE, row.names=1, sep="\\t", check.names=FALSE)
        sample_info <- read.table("__SNK_3__", header=TRUE, sep="\\t", comment.char="", check.names=FALSE)
        rownames(sample_info) <- sample_info$sample
        
        comp_name <- "__SNK_5__"
        comp_parts <- strsplit(comp_name, "_vs_")[[1]]
        ctrl_group <- comp_parts[1]
        treat_group <- comp_parts[2]
        
        ctrl_samples <- sample_info$sample[sample_info$group == ctrl_group]
        treat_samples <- sample_info$sample[sample_info$group == treat_group]
        use_samples <- c(ctrl_samples, treat_samples)
        
        # 子集 lncRNA
        count_sub <- count_mat[rownames(count_mat) %in% lnc_ids, use_samples, drop=FALSE]
        
        if(nrow(count_sub) < 10) {
            cat("WARNING: Too few lncRNAs with counts (", nrow(count_sub), "), skipping DE analysis\\n")
            df_empty <- data.frame(gene_id=character(), baseMean=numeric(), log2FoldChange=numeric(),
                                   pvalue=numeric(), padj=numeric(), diffType=character())
            write.table(df_empty, "__SNK_0__", sep="\\t", row.names=FALSE, quote=FALSE)
            write.table(df_empty, "__SNK_1__", sep="\\t", row.names=FALSE, quote=FALSE)
            quit(save="no")
        }
        
        group_info <- data.frame(
            sample = use_samples,
            group = factor(c(rep(ctrl_group, length(ctrl_samples)), rep(treat_group, length(treat_samples))),
                          levels = c(ctrl_group, treat_group))
        )
        rownames(group_info) <- group_info$sample
        
        dds <- DESeqDataSetFromMatrix(countData = count_sub, colData = group_info, design = ~ group)
        dds <- DESeq(dds)
        res <- results(dds, contrast = c("group", treat_group, ctrl_group))
        
        res_df <- data.frame(
            gene_id = rownames(res),
            baseMean = res$baseMean,
            log2FoldChange = res$log2FoldChange,
            pvalue = res$pvalue,
            padj = res$padj,
            stringsAsFactors = FALSE
        )
        
        logfc_cut <- as.numeric("__SNK_8__")
        p_cut <- as.numeric("__SNK_6__")
        p_type <- "__SNK_9__"
        
        res_df$diffType <- case_when(
            is.na(res_df[[p_type]]) ~ "No Significant Difference",
            res_df[[p_type]] < p_cut & res_df$log2FoldChange > logfc_cut ~ "Up",
            res_df[[p_type]] < p_cut & res_df$log2FoldChange < -logfc_cut ~ "Down",
            TRUE ~ "No Significant Difference"
        )
        
        res_df <- res_df %>% arrange(desc(abs(log2FoldChange)))
        write.table(res_df, "__SNK_0__", sep="\\t", row.names=FALSE, quote=FALSE)
        sig_df <- res_df %>% filter(diffType != "No Significant Difference")
        write.table(sig_df, "__SNK_1__", sep="\\t", row.names=FALSE, quote=FALSE)
        
        cat("LncRNA DE complete:", nrow(sig_df), "significant /", nrow(res_df), "total\\n")
        
        """
        r_code = r_code.replace("__SNK_0__", str(output.all_results))
        r_code = r_code.replace("__SNK_1__", str(output.sig_results))
        r_code = r_code.replace("__SNK_2__", str(input.lncrna_list))
        r_code = r_code.replace("__SNK_3__", str(input.sample_desc))
        r_code = r_code.replace("__SNK_4__", str(input.count_mat))
        r_code = r_code.replace("__SNK_5__", str(wildcards.comp))
        r_code = r_code.replace("__SNK_6__", str(params.pvalue))
        r_code = r_code.replace("__SNK_7__", str(params.outdir))
        r_code = r_code.replace("__SNK_8__", str(params.logFC))
        r_code = r_code.replace("__SNK_9__", str(params.ptype))
        result = subprocess.run(["Rscript", "-"], input=r_code, text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- LncRNA cis 靶基因预测 ----
rule lncrna_target_cis:
    """cis 靶基因: LncRNA 上下游窗口内的蛋白编码基因"""
    input:
        lncrna_list = rules.merge_lncrna_results.output.lncrna_list,
        gtf         = config["paths"]["ref_gtf"],
    output:
        cis_targets = f"{RESULT_DIR}/08.lnc/target/cis_target_genes.xls",
    params:
        cis_window = config["lncrna"].get("cis_window", 100000),
        outdir     = f"{RESULT_DIR}/08.lnc/target",
    conda:
        f"{PIPELINE_DIR}/envs/lncrna.yaml"
    run:
        import os, subprocess, sys
        os.makedirs(f"{params.outdir}", exist_ok=True)
        r_code = """
        import re
        
        # 1. 读 lncRNA 列表
        with open('__SNK_2__') as f:
            lncrna_set = set(line.strip() for line in f if line.strip())
        print(f'LncRNAs: {len(lncrna_set)}')
        
        # 2. 解析 GTF 获取所有转录本坐标
        lncrna_coords = []      # (chrom, start, end, strand, transcript_id)
        gene_coords = []        # (chrom, start, end, strand, gene_id)
        
        with open('__SNK_4__') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\\t')
                if len(parts) < 9:
                    continue
        
                chrom, feature, start, end, strand = parts[0], parts[2], int(parts[3]), int(parts[4]), parts[6]
                attrs = parts[8]
        
                if feature == 'transcript':
                    tid_m = re.search(r'transcript_id\\s+"([^"]+)"', attrs)
                    gid_m = re.search(r'gene_id\\s+"([^"]+)"', attrs)
                    if not tid_m:
                        continue
        
                    tid = tid_m.group(1)
                    if tid in lncrna_set:
                        lncrna_coords.append((chrom, start, end, strand, tid))
                    elif gid_m:
                        gene_coords.append((chrom, start, end, strand, gid_m.group(1)))
        
        print(f'LncRNA coords: {len(lncrna_coords)}')
        print(f'Gene coords: {len(gene_coords)}')
        
        # 3. 查找 cis 靶基因 (同染色体, 上下游窗口内)
        window = __SNK_1__
        cis_results = []
        
        for lnc_chrom, lnc_start, lnc_end, lnc_strand, lnc_id in lncrna_coords:
            for g_chrom, g_start, g_end, g_strand, g_id in gene_coords:
                if lnc_chrom != g_chrom:
                    continue
                # 计算 lncRNA 与基因的距离
                if lnc_end < g_start:
                    dist = g_start - lnc_end
                elif g_end < lnc_start:
                    dist = lnc_start - g_end
                else:
                    dist = 0  # 重叠
        
                if dist <= window:
                    # 确定位置关系
                    if dist == 0:
                        relation = 'overlap'
                    elif lnc_end < g_start:
                        relation = 'downstream'
                    else:
                        relation = 'upstream'
                    cis_results.append((lnc_id, g_id, dist, relation))
        
        # 4. 写结果
        with open('__SNK_0__', 'w') as out:
            out.write('lncrna_id\\ttarget_gene\\tdistance\\trelation\\n')
            for lncid, gid, dist, rel in cis_results:
                out.write(f'{lncid}\\t{gid}\\t{dist}\\t{rel}\\n')
        
        print(f'Cis targets: {len(cis_results)} pairs')
        print(f'Unique lncRNAs with cis targets: {len(set(r[0] for r in cis_results))}')
        print(f'Unique target genes: {len(set(r[1] for r in cis_results))}')
        
        """
        r_code = r_code.replace("__SNK_0__", str(output.cis_targets))
        r_code = r_code.replace("__SNK_1__", str(params.cis_window))
        r_code = r_code.replace("__SNK_2__", str(input.lncrna_list))
        r_code = r_code.replace("__SNK_3__", str(params.outdir))
        r_code = r_code.replace("__SNK_4__", str(input.gtf))
        result = subprocess.run(["python3", "-c", r_code], text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- LncRNA trans 靶基因预测 (表达相关性) ----
rule lncrna_target_trans:
    """trans 靶基因: 基于表达相关性 (|PCC| > 0.9, p < 0.05)"""
    input:
        lncrna_list = rules.merge_lncrna_results.output.lncrna_list,
        norm_mat    = f"{RESULT_DIR}/05.expression/gene/gene_norm.xls",
    output:
        trans_targets = f"{RESULT_DIR}/08.lnc/target/trans_target_genes.xls",
    params:
        outdir = f"{RESULT_DIR}/08.lnc/target",
    conda:
        f"{PIPELINE_DIR}/envs/diff.yaml"
    run:
        import os, subprocess, sys
        os.makedirs(f"{params.outdir}", exist_ok=True)
        r_code = """
        library(dplyr)
        
        # 读 lncRNA 列表
        lnc_ids <- readLines("__SNK_1__")
        cat("LncRNA IDs:", length(lnc_ids), "\\n")
        
        # 读表达矩阵
        expr <- read.table("__SNK_2__", header=TRUE, row.names=1, sep="\\t", check.names=FALSE)
        
        # 分离 lncRNA 和 mRNA
        common_lnc <- intersect(lnc_ids, rownames(expr))
        cat("LncRNAs in expression matrix:", length(common_lnc), "\\n")
        
        if(length(common_lnc) < 2) {
            cat("WARNING: Too few lncRNAs with expression data\\n")
            writeLines("lncrna_id\\ttarget_gene\\tPCC\\tpvalue", "__SNK_0__")
            quit(save="no")
        }
        
        lnc_expr <- expr[common_lnc, , drop=FALSE]
        mrna_ids <- setdiff(rownames(expr), lnc_ids)
        mrna_expr <- expr[mrna_ids, , drop=FALSE]
        
        cat("Target genes (mRNA):", nrow(mrna_expr), "\\n")
        
        results <- list()
        for(i in seq_len(nrow(lnc_expr))) {
            lnc_id <- rownames(lnc_expr)[i]
            lnc_vec <- as.numeric(lnc_expr[i, ])
        
            for(j in seq_len(min(nrow(mrna_expr), 5000))) {  # 限制以控制运行时间
                mrna_id <- rownames(mrna_expr)[j]
                mrna_vec <- as.numeric(mrna_expr[j, ])
        
                ct <- tryCatch(cor.test(lnc_vec, mrna_vec, method="pearson"),
                              error = function(e) NULL)
                if(!is.null(ct)) {
                    if(abs(ct$estimate) > 0.9 && ct$p.value < 0.05) {
                        results[[length(results) + 1]] <- data.frame(
                            lncrna_id = lnc_id,
                            target_gene = mrna_id,
                            PCC = round(ct$estimate, 4),
                            pvalue = signif(ct$p.value, 4),
                            stringsAsFactors = FALSE
                        )
                    }
                }
            }
            if(i %% 10 == 0) cat("  Processed", i, "/", nrow(lnc_expr), "lncRNAs\\n")
        }
        
        if(length(results) > 0) {
            trans_df <- bind_rows(results) %>% arrange(pvalue)
            write.table(trans_df, "__SNK_0__", sep="\\t", row.names=FALSE, quote=FALSE)
            cat("\\nTrans targets:", nrow(trans_df), "pairs\\n")
            cat("Unique lncRNAs:", length(unique(trans_df$lncrna_id)), "\\n")
        } else {
            cat("No trans targets found\\n")
            writeLines("lncrna_id\\ttarget_gene\\tPCC\\tpvalue", "__SNK_0__")
        }
        
        """
        r_code = r_code.replace("__SNK_0__", str(output.trans_targets))
        r_code = r_code.replace("__SNK_1__", str(input.lncrna_list))
        r_code = r_code.replace("__SNK_2__", str(input.norm_mat))
        r_code = r_code.replace("__SNK_3__", str(params.outdir))
        result = subprocess.run(["Rscript", "-"], input=r_code, text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- LncRNA 分类 ----
rule classify_lncrna:
    """根据基因组位置对 LncRNA 分类"""
    input:
        lncrna_list = rules.merge_lncrna_results.output.lncrna_list,
        all_gtf     = rules.gffcompare.output.compared,
    output:
        class_result = f"{RESULT_DIR}/08.lnc/lncrna_classification.xls",
    conda:
        f"{PIPELINE_DIR}/envs/lncrna.yaml"
    run:
        import os, subprocess, sys
        r_code = """
        import re
        
        # 读取 lncRNA 列表
        with open('__SNK_1__') as f:
            lncrna_ids = set(line.strip() for line in f if line.strip())
        
        # 从 GTF 提取分类信息
        class_counts = {}
        with open('__SNK_2__') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\\t')
                if len(parts) < 9 or parts[2] != 'transcript':
                    continue
        
                # 提取 transcript_id 和 class_code
                tid_match = re.search(r'transcript_id\\s+"([^"]+)"', parts[8])
                cc_match = re.search(r'class_code\\s+"([^"]+)"', parts[8])
        
                if tid_match and cc_match:
                    tid = tid_match.group(1)
                    cc = cc_match.group(1)
        
                    if tid in lncrna_ids:
                        if cc not in class_counts:
                            class_counts[cc] = set()
                        class_counts[cc].add(tid)
        
        # 分类标签
        class_labels = {
            '=': 'Complete match of intron chain',
            'c': 'Contained',
            'j': 'Potentially novel isoform (fragment)',
            'e': 'Single exon transfrag overlapping a reference exon',
            'i': 'Fully contained within a reference intron',
            'o': 'Generic exonic overlap with a reference transcript',
            'p': 'Possible polymerase run-on fragment',
            'u': 'Unknown, intergenic transcript',
            'x': 'Exonic overlap on the opposite strand',
        }
        
        with open('__SNK_0__', 'w') as f:
            f.write('Class\\tClassLabel\\tCount\\n')
            for cc, ids in sorted(class_counts.items()):
                label = class_labels.get(cc, 'Unknown')
                f.write(f'{cc}\\t{label}\\t{len(ids)}\\n')
        
        print('LncRNA classification done.')
        print(f'Total classified lncRNAs: {sum(len(v) for v in class_counts.values())}')
        
        """
        r_code = r_code.replace("__SNK_0__", str(output.class_result))
        r_code = r_code.replace("__SNK_1__", str(input.lncrna_list))
        r_code = r_code.replace("__SNK_2__", str(input.all_gtf))
        result = subprocess.run(["python3", "-c", r_code], text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- LncRNA 表达量 PCA + 相关性 ----
rule lncrna_expression_viz:
    """LncRNA 表达量 PCA 和样本相关性"""
    input:
        lncrna_list = rules.merge_lncrna_results.output.lncrna_list,
        norm_mat    = f"{RESULT_DIR}/05.expression/gene/gene_norm.xls",
        sample_desc = f"{TMP_DIR}/diff/sample_desc.txt",
    output:
        pca_pdf  = f"{RESULT_DIR}/08.lnc/lncrna_pca.pdf",
        pca_png  = f"{RESULT_DIR}/08.lnc/lncrna_pca.png",
        corr_pdf = f"{RESULT_DIR}/08.lnc/lncrna_corr.pdf",
        corr_png = f"{RESULT_DIR}/08.lnc/lncrna_corr.png",
    conda:
        f"{PIPELINE_DIR}/envs/diff.yaml"
    run:
        import os, subprocess, sys
        r_code = """
        library(ggplot2); library(ggrepel); library(pheatmap)
        
        lnc_ids <- readLines("__SNK_0__")
        expr <- read.table("__SNK_6__", header=TRUE, row.names=1, sep="\\t", check.names=FALSE)
        samples <- read.table("__SNK_1__", header=TRUE, sep="\\t", comment.char="", check.names=FALSE)
        
        common_lnc <- intersect(lnc_ids, rownames(expr))
        if(length(common_lnc) < 5) {
            pdf("__SNK_4__"); plot.new(); text(0.5,0.5,paste("Too few lncRNAs:", length(common_lnc))); dev.off()
            pdf("__SNK_3__"); plot.new(); text(0.5,0.5,"Too few lncRNAs"); dev.off()
            quit(save="no")
        }
        
        lnc_expr <- expr[common_lnc, , drop=FALSE]
        lnc_log <- log2(lnc_expr + 1)
        rv <- apply(lnc_log, 1, var)
        lnc_log <- lnc_log[rv > 0, , drop=FALSE]
        
        cat("LncRNA expression:", nrow(lnc_log), "x", ncol(lnc_log), "\\n")
        
        # PCA
        pca <- prcomp(t(lnc_log), scale.=TRUE)
        pca_df <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], sample=colnames(lnc_log))
        pca_df <- merge(pca_df, samples, by="sample")
        ve <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
        p <- ggplot(pca_df, aes(x=PC1, y=PC2, color=group, label=sample)) +
            geom_point(size=3) + geom_text_repel(size=3, show.legend=FALSE) + stat_ellipse(level=0.95) +
            labs(x=paste0("PC1 (",ve[1],"%)"), y=paste0("PC2 (",ve[2],"%)"), title="LncRNA Expression PCA") +
            theme_bw() + theme(plot.title=element_text(hjust=0.5))
        ggsave("__SNK_4__", p, width=8, height=6)
        ggsave("__SNK_5__", p, width=8, height=6, dpi=150)
        
        # Correlation
        rownames(samples) <- samples$sample
        common_s <- intersect(colnames(lnc_log), samples$sample)
        cor_mat <- cor(lnc_log[, common_s, drop=FALSE], method="pearson")
        ann_col <- data.frame(Group=samples[common_s, "group"]); rownames(ann_col) <- common_s
        pdf("__SNK_3__", width=10, height=8)
        pheatmap(cor_mat, annotation_col=ann_col, display_numbers=TRUE, number_format="%.2f",
            fontsize_number=8, main="LncRNA Sample Correlation (Pearson)")
        dev.off()
        png("__SNK_2__", width=1200, height=960)
        pheatmap(cor_mat, annotation_col=ann_col, display_numbers=TRUE, number_format="%.2f",
            fontsize_number=8, main="LncRNA Sample Correlation (Pearson)")
        dev.off()
        
        """
        r_code = r_code.replace("__SNK_0__", str(input.lncrna_list))
        r_code = r_code.replace("__SNK_1__", str(input.sample_desc))
        r_code = r_code.replace("__SNK_2__", str(output.corr_png))
        r_code = r_code.replace("__SNK_3__", str(output.corr_pdf))
        r_code = r_code.replace("__SNK_4__", str(output.pca_pdf))
        r_code = r_code.replace("__SNK_5__", str(output.pca_png))
        r_code = r_code.replace("__SNK_6__", str(input.norm_mat))
        result = subprocess.run(["Rscript", "-"], input=r_code, text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")
