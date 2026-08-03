# ============================================================
# 差异分析模块: DESeq2 / edgeR + 可视化
# 支持: 多对比较、批次效应校正(可选)
# ============================================================

# ---- 样本分组文件 ----
rule make_sample_desc:
    """生成差异分析用的样本分组文件 (含可选的 batch 列)"""
    input:
        sample_info = config["paths"]["sample_info"],
    output:
        sample_desc = f"{TMP_DIR}/diff/sample_desc.txt",
    conda:
        f"{PIPELINE_DIR}/envs/qc.yaml"
    run:
        import os
        import pandas as pd

        os.makedirs(os.path.dirname(str(output.sample_desc)), exist_ok=True)

        df = pd.read_csv(str(input.sample_info), sep='\t', comment='#')
        cols = ['sample', 'group']
        if 'batch' in df.columns:
            cols.append('batch')
        df[cols].to_csv(str(output.sample_desc), sep='\t', index=False)
        print(f"Sample desc: {len(df)} samples, {df['group'].nunique()} groups")
        if 'batch' in df.columns:
            print(f"Batch info: {df['batch'].nunique()} batches")


# ---- PCA 分析 ----
rule pca_analysis:
    """PCA 主成分分析"""
    input:
        expr_mat    = rules.normalize_expression.output.norm_mat,
        sample_desc = rules.make_sample_desc.output.sample_desc,
    output:
        pca_pdf   = f"{RESULT_DIR}/06.de_enrich/PCA.pdf",
        pca_png   = f"{RESULT_DIR}/06.de_enrich/PCA.png",
    conda:
        f"{PIPELINE_DIR}/envs/diff.yaml"
    run:
        import os, subprocess, sys
        os.makedirs(f"{RESULT_DIR}/06.de_enrich", exist_ok=True)
        r_code = """
        library(ggplot2)
        library(ggrepel)
        
        expr <- read.table("__SNK_1__", header=TRUE, row.names=1, sep="\\t", check.names=FALSE)
        samples <- read.table("__SNK_0__", header=TRUE, sep="\\t", comment.char="", check.names=FALSE)
        
        expr_log <- log2(expr + 1)
        
        # 过滤掉方差为 0 的基因（所有样本表达值相同），避免 prcomp(scale.=TRUE) 报错
        row_vars <- apply(expr_log, 1, var)
        expr_log <- expr_log[row_vars > 0, ]
        cat("Genes retained after variance filter:", nrow(expr_log), "/", length(row_vars), "\\n")
        
        pca <- prcomp(t(expr_log), scale.=TRUE)
        pca_df <- data.frame(
            PC1 = pca$x[,1],
            PC2 = pca$x[,2],
            sample = colnames(expr_log)
        )
        pca_df <- merge(pca_df, samples, by="sample")
        
        var_explained <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
        
        # 批次信息用于形状区分
        if("batch" %in% colnames(pca_df)) {
            pca_df$batch <- factor(pca_df$batch)
            p <- ggplot(pca_df, aes(x=PC1, y=PC2, color=group, shape=batch, label=sample)) +
                geom_point(size=3) +
                geom_text_repel(size=3, show.legend=FALSE) +
                stat_ellipse(aes(group=group), level=0.95) +
                labs(
                    x = paste0("PC1 (", var_explained[1], "%)"),
                    y = paste0("PC2 (", var_explained[2], "%)"),
                    title = "PCA Analysis (shape = batch)"
                ) +
                theme_bw() + theme(plot.title = element_text(hjust=0.5))
        } else {
            p <- ggplot(pca_df, aes(x=PC1, y=PC2, color=group, label=sample)) +
                geom_point(size=3) +
                geom_text_repel(size=3, show.legend=FALSE) +
                stat_ellipse(level=0.95) +
                labs(
                    x = paste0("PC1 (", var_explained[1], "%)"),
                    y = paste0("PC2 (", var_explained[2], "%)"),
                    title = "PCA Analysis"
                ) +
                theme_bw() + theme(plot.title = element_text(hjust=0.5))
        }
        
        ggsave("__SNK_3__", p, width=8, height=6)
        ggsave("__SNK_2__", p, width=8, height=6, dpi=150)
        
        # 如果有多批次，额外生成按批次着色的 PCA
        if("batch" %in% colnames(pca_df) && length(unique(pca_df$batch)) > 1) {
            p2 <- ggplot(pca_df, aes(x=PC1, y=PC2, color=batch, label=sample)) +
                geom_point(size=3) +
                geom_text_repel(size=3, show.legend=FALSE) +
                stat_ellipse(level=0.95) +
                labs(
                    x = paste0("PC1 (", var_explained[1], "%)"),
                    y = paste0("PC2 (", var_explained[2], "%)"),
                    title = "PCA Analysis (colored by batch)"
                ) +
                theme_bw() + theme(plot.title = element_text(hjust=0.5))
            ggsave(sub("[.]pdf$", "_by_batch.pdf", "__SNK_3__"), p2, width=8, height=6)
        }
        
        """
        r_code = r_code.replace("__SNK_0__", str(input.sample_desc))
        r_code = r_code.replace("__SNK_1__", str(input.expr_mat))
        r_code = r_code.replace("__SNK_2__", str(output.pca_png))
        r_code = r_code.replace("__SNK_3__", str(output.pca_pdf))
        r_code = r_code.replace("__SNK_4__", str(RESULT_DIR))
        result = subprocess.run(["Rscript", "-"], input=r_code, text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- 相关性热图 ----
rule correlation_heatmap:
    """样本相关性热图"""
    input:
        expr_mat    = rules.normalize_expression.output.norm_mat,
        sample_desc = rules.make_sample_desc.output.sample_desc,
    output:
        corr_pdf = f"{RESULT_DIR}/06.de_enrich/correlation_heatmap.pdf",
        corr_png = f"{RESULT_DIR}/06.de_enrich/correlation_heatmap.png",
    conda:
        f"{PIPELINE_DIR}/envs/diff.yaml"
    run:
        import os, subprocess, sys
        os.makedirs(f"{RESULT_DIR}/06.de_enrich", exist_ok=True)
        r_code = """
        library(pheatmap)
        library(tidyverse)
        
        expr <- read.table("__SNK_3__", header=TRUE, row.names=1, sep="\\t", check.names=FALSE)
        samples <- read.table("__SNK_0__", header=TRUE, sep="\\t", comment.char="", check.names=FALSE)
        
        expr_log <- log2(expr + 1)
        cor_mat <- cor(expr_log, method="pearson")
        
        rownames(samples) <- samples$sample
        ann_col <- data.frame(Group = samples[colnames(cor_mat), "group"])
        rownames(ann_col) <- colnames(cor_mat)
        if("batch" %in% colnames(samples)) {
            ann_col$Batch <- factor(samples[colnames(cor_mat), "batch"])
        }
        
        pdf("__SNK_2__", width=10, height=8)
        pheatmap(cor_mat,
            annotation_col = ann_col,
            display_numbers = TRUE,
            number_format = "%.2f",
            fontsize_number = 8,
            main = "Sample Correlation (Pearson)")
        dev.off()
        
        png("__SNK_1__", width=1200, height=960)
        pheatmap(cor_mat,
            annotation_col = ann_col,
            display_numbers = TRUE,
            number_format = "%.2f",
            fontsize_number = 8,
            main = "Sample Correlation (Pearson)")
        dev.off()
        
        """
        r_code = r_code.replace("__SNK_0__", str(input.sample_desc))
        r_code = r_code.replace("__SNK_1__", str(output.corr_png))
        r_code = r_code.replace("__SNK_2__", str(output.corr_pdf))
        r_code = r_code.replace("__SNK_3__", str(input.expr_mat))
        r_code = r_code.replace("__SNK_4__", str(RESULT_DIR))
        result = subprocess.run(["Rscript", "-"], input=r_code, text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- 差异分析 (DESeq2/edgeR) ----
# 使用 COMPARISON_MAP 获取正确的 ctrl/treat 组名
rule diff_analysis:
    """差异表达分析 — 支持多对比较和批次效应校正"""
    input:
        count_mat   = rules.merge_gene_counts.output.count_mat,
        sample_desc = rules.make_sample_desc.output.sample_desc,
    output:
        all_results  = f"{RESULT_DIR}/06.de_enrich/{{comp}}/{{comp}}.DE_results.All.xls",
        sig_results  = f"{RESULT_DIR}/06.de_enrich/{{comp}}/{{comp}}.DE_results.Significant.xls",
        up_id_list   = f"{RESULT_DIR}/06.de_enrich/{{comp}}/{{comp}}.Significant.Up.id_list",
        down_id_list = f"{RESULT_DIR}/06.de_enrich/{{comp}}/{{comp}}.Significant.Down.id_list",
    params:
        logFC       = config["diff"]["logFC"],
        pvalue      = config["diff"]["diff_pvalue"],
        ptype       = config["diff"]["pOrFDR"],
        method      = config["diff"]["method"],
        dispersion  = config["diff"]["dispersion"],
        fc_type     = config["diff"]["fc_calculate_type"],
        batch_effect = config["diff"].get("batch_effect", False),
        outdir      = f"{RESULT_DIR}/06.de_enrich/{{comp}}",
    conda:
        f"{PIPELINE_DIR}/envs/diff.yaml"
    run:
        import os, subprocess, sys
        os.makedirs(params.outdir, exist_ok=True)
        r_code = """
        library(DESeq2)
        library(edgeR)
        library(tidyverse)
        
        # ---- 检查输入 ----
        if(!file.exists("__SNK_7__")) {
            stop("ERROR: Count matrix not found: __SNK_7__")
        }
        if(!file.exists("__SNK_4__")) {
            stop("ERROR: Sample description not found: __SNK_4__")
        }
        
        # ---- 读数据 ----
        count_mat <- read.table("__SNK_7__", header=TRUE, row.names=1, sep="\\t", check.names=FALSE)
        sample_info <- read.table("__SNK_4__", header=TRUE, sep="\\t", comment.char="", check.names=FALSE)
        rownames(sample_info) <- sample_info$sample
        
        comp_name <- "__SNK_9__"
        comp_parts <- strsplit(comp_name, "_vs_")[[1]]
        ctrl_group <- comp_parts[1]
        treat_group <- comp_parts[2]
        use_batch <- as.logical("__SNK_0__") && ("batch" %in% colnames(sample_info))
        
        cat("Comparison:", ctrl_group, "vs", treat_group, "\\n")
        cat("Batch correction:", use_batch, "\\n")
        
        # ---- 提取比较组样本 ----
        ctrl_samples <- sample_info$sample[sample_info$group == ctrl_group]
        treat_samples <- sample_info$sample[sample_info$group == treat_group]
        use_samples <- c(ctrl_samples, treat_samples)
        
        if(length(use_samples) < 3) {
            cat("ERROR: 比较组的样本总数不足 (需要 >= 3) \\n")
            cat("  对照组", ctrl_group, ":", length(ctrl_samples), "个样本\\n")
            cat("  实验组", treat_group, ":", length(treat_samples), "个样本\\n")
            quit(status=1)
        }
        
        # ---- 子集数据 ----
        count_sub <- count_mat[, use_samples, drop=FALSE]
        group_info <- data.frame(
            sample = use_samples,
            group = c(rep(ctrl_group, length(ctrl_samples)),
                      rep(treat_group, length(treat_samples)))
        )
        if(use_batch) {
            group_info$batch <- factor(sample_info[use_samples, "batch"])
            cat("Batch levels:", levels(group_info$batch), "\\n")
        }
        group_info$group <- factor(group_info$group,
            levels = c(ctrl_group, treat_group))
        rownames(group_info) <- group_info$sample
        
        # ---- 选择方法 ----
        method <- "__SNK_12__"
        if(method == "auto") {
            if(length(ctrl_samples) >= 2 && length(treat_samples) >= 2) {
                method <- "DESeq2"
            } else {
                method <- "edgeR"
            }
        }
        cat("Using method:", method, "\\n")
        
        fc_type <- as.integer("__SNK_8__")
        logfc_cut <- as.numeric("__SNK_13__")
        p_cut <- as.numeric("__SNK_11__")
        p_type <- "__SNK_14__"
        
        # ---- 分析（带错误捕获，失败时输出诊断信息）----
        sink("__SNK_10__/diff_debug.log", split=TRUE)
        run_success <- tryCatch({
            if(method == "DESeq2") {
                if(use_batch) {
                    cat("DESeq2 design: ~ batch + group\\n")
                    dds <- DESeqDataSetFromMatrix(
                        countData = count_sub,
                        colData = group_info,
                        design = ~ batch + group
                    )
                } else {
                    cat("DESeq2 design: ~ group\\n")
                    dds <- DESeqDataSetFromMatrix(
                        countData = count_sub,
                        colData = group_info,
                        design = ~ group
                    )
                }
                cat("Running DESeq2...\\n")
                dds <- DESeq(dds)
                cat("DESeq2 done. Extracting results...\\n")
                res <- results(dds, contrast = c("group", treat_group, ctrl_group))
                if(fc_type == 1) {
                    res$log2FoldChange <- -res$log2FoldChange
                }
                norm_counts <- counts(dds, normalized=TRUE)
                baseMean_ctrl <- rowMeans(norm_counts[, ctrl_samples, drop=FALSE])
                baseMean_treat <- rowMeans(norm_counts[, treat_samples, drop=FALSE])
            } else {
                if(use_batch) {
                    cat("edgeR with batch correction via removeBatchEffect\\n")
                }
                y <- DGEList(counts = count_sub, group = group_info$group)
                keep <- filterByExpr(y)
                y <- y[keep, , keep.lib.sizes=FALSE]
                y <- calcNormFactors(y)
        
                if(length(ctrl_samples) >= 2 && length(treat_samples) >= 2) {
                    if(use_batch) {
                        design_mat <- model.matrix(~ batch + group, data=group_info)
                        y <- estimateDisp(y, design=design_mat)
                        fit <- glmQLFit(y, design=design_mat)
                        qlf <- glmQLFTest(fit, coef=ncol(design_mat))
                    } else {
                        y <- estimateDisp(y)
                        fit <- glmQLFit(y)
                        qlf <- glmQLFTest(fit, contrast = c(-1, 1))
                    }
                } else {
                    disp <- as.numeric("__SNK_5__")
                    et <- exactTest(y, dispersion = disp^2)
                    qlf <- et
                }
        
                res <- data.frame(
                    log2FoldChange = if(fc_type == 1) -qlf$table$logFC else qlf$table$logFC,
                    pvalue = qlf$table$PValue,
                    row.names = rownames(qlf$table)
                )
                res$padj <- p.adjust(res$pvalue, method="BH")
                res$baseMean <- rowMeans(y$counts)
        
                norm_counts <- cpm(y, log=FALSE)
                baseMean_ctrl <- rowMeans(norm_counts[, ctrl_samples, drop=FALSE])
                baseMean_treat <- rowMeans(norm_counts[, treat_samples, drop=FALSE])
            }
            TRUE  # 返回成功标记
        }, error = function(e) {
            cat("ERROR in differential analysis:", conditionMessage(e), "\\n")
            cat("Call:", paste(deparse(conditionCall(e)), collapse="\\n"), "\\n")
            return(FALSE)
        })
        sink()
        
        if(!run_success) {
            quit(status=1)
        }
        
        # ---- 构建结果表 ----
        res_df <- data.frame(
            gene_id = rownames(res),
            baseMean = res$baseMean,
            baseMean_ctrl = baseMean_ctrl[rownames(res)],
            baseMean_treat = baseMean_treat[rownames(res)],
            log2FoldChange = res$log2FoldChange,
            pvalue = res$pvalue,
            padj = res$padj,
            stringsAsFactors = FALSE
        )
        
        # 判断差异类型
        res_df$diffType <- case_when(
            is.na(res_df[[p_type]]) ~ "No Significant Difference",
            res_df[[p_type]] < p_cut & res_df$log2FoldChange > logfc_cut ~ "Up",
            res_df[[p_type]] < p_cut & res_df$log2FoldChange < -logfc_cut ~ "Down",
            TRUE ~ "No Significant Difference"
        )
        
        # 排序
        res_df <- res_df %>% arrange(desc(abs(log2FoldChange)))
        
        # 保存
        write.table(res_df, "__SNK_3__", sep="\\t", row.names=FALSE, quote=FALSE)
        cat("Written:", "__SNK_3__", "(", nrow(res_df), "rows )\\n")
        
        sig_df <- res_df %>% filter(diffType != "No Significant Difference")
        write.table(sig_df, "__SNK_2__", sep="\\t", row.names=FALSE, quote=FALSE)
        cat("Written:", "__SNK_2__", "(", nrow(sig_df), "rows )\\n")
        
        up_ids <- sig_df %>% filter(diffType == "Up") %>% pull(gene_id)
        down_ids <- sig_df %>% filter(diffType == "Down") %>% pull(gene_id)
        writeLines(up_ids, "__SNK_6__")
        writeLines(down_ids, "__SNK_1__")
        
        cat("\\n=== 差异分析完成 ===\\n")
        cat("  Comparison:", ctrl_group, "vs", treat_group, "\\n")
        cat("  Method:", method, "\\n")
        cat("  Batch corrected:", use_batch, "\\n")
        cat("  Total genes:", nrow(res_df), "\\n")
        cat("  Up-regulated:", length(up_ids), "\\n")
        cat("  Down-regulated:", length(down_ids), "\\n")
        cat("  Threshold: |log2FC| >", logfc_cut, ",", p_type, "<", p_cut, "\\n")
        
        """
        r_code = r_code.replace("__SNK_0__", str(params.batch_effect))
        r_code = r_code.replace("__SNK_1__", str(output.down_id_list))
        r_code = r_code.replace("__SNK_2__", str(output.sig_results))
        r_code = r_code.replace("__SNK_3__", str(output.all_results))
        r_code = r_code.replace("__SNK_4__", str(input.sample_desc))
        r_code = r_code.replace("__SNK_5__", str(params.dispersion))
        r_code = r_code.replace("__SNK_6__", str(output.up_id_list))
        r_code = r_code.replace("__SNK_7__", str(input.count_mat))
        r_code = r_code.replace("__SNK_8__", str(params.fc_type))
        r_code = r_code.replace("__SNK_9__", str(wildcards.comp))
        r_code = r_code.replace("__SNK_10__", str(params.outdir))
        r_code = r_code.replace("__SNK_11__", str(params.pvalue))
        r_code = r_code.replace("__SNK_12__", str(params.method))
        r_code = r_code.replace("__SNK_13__", str(params.logFC))
        r_code = r_code.replace("__SNK_14__", str(params.ptype))
        result = subprocess.run(["Rscript", "-"], input=r_code, text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- 火山图 ----
rule volcano_plot:
    """差异表达火山图"""
    input:
        de_results = f"{RESULT_DIR}/06.de_enrich/{{comp}}/{{comp}}.DE_results.All.xls",
    output:
        volcano_pdf = f"{RESULT_DIR}/06.de_enrich/{{comp}}/{{comp}}.volcano.pdf",
        volcano_png = f"{RESULT_DIR}/06.de_enrich/{{comp}}/{{comp}}.volcano.png",
    params:
        logFC  = config["diff"]["logFC"],
        pvalue = config["diff"]["diff_pvalue"],
        ptype  = config["diff"]["pOrFDR"],
    conda:
        f"{PIPELINE_DIR}/envs/diff.yaml"
    run:
        import os, subprocess, sys
        r_code = """
        library(ggplot2)
        library(ggrepel)
        library(dplyr)
        
        df <- read.table("__SNK_2__", header=TRUE, sep="\\t", check.names=FALSE)
        logfc_cut <- as.numeric("__SNK_4__")
        p_cut <- as.numeric("__SNK_3__")
        p_type <- "__SNK_5__"
        
        if(!"pvalue" %in% colnames(df)) {
            cat("ERROR: Missing pvalue/padj column in DE results\\n")
            quit(status=1)
        }
        
        df$log10p <- -log10(df[[p_type]])
        df$log10p[is.infinite(df$log10p)] <- max(df$log10p[is.finite(df$log10p)], na.rm=TRUE)
        
        top_up <- df %>% filter(diffType == "Up") %>% arrange(pvalue) %>% head(10)
        top_down <- df %>% filter(diffType == "Down") %>% arrange(pvalue) %>% head(10)
        top_genes <- bind_rows(top_up, top_down)
        
        p <- ggplot(df, aes(x=log2FoldChange, y=log10p, color=diffType)) +
            geom_point(size=1, alpha=0.6) +
            scale_color_manual(values = c(
                "Up" = "#E74C3C",
                "Down" = "#2980B9",
                "No Significant Difference" = "grey70"
            )) +
            geom_vline(xintercept = c(-logfc_cut, logfc_cut), linetype = "dashed", alpha=0.5) +
            geom_hline(yintercept = -log10(p_cut), linetype = "dashed", alpha=0.5) +
            geom_text_repel(
                data = top_genes,
                aes(label = gene_id),
                size = 3,
                max.overlaps = 20
            ) +
            labs(
                x = "log2(Fold Change)",
                y = paste0("-log10(", p_type, ")"),
                title = "Volcano Plot"
            ) +
            theme_bw() +
            theme(
                legend.position = "bottom",
                plot.title = element_text(hjust=0.5)
            )
        
        ggsave("__SNK_0__", p, width=8, height=7)
        ggsave("__SNK_1__", p, width=8, height=7, dpi=150)
        
        """
        r_code = r_code.replace("__SNK_0__", str(output.volcano_pdf))
        r_code = r_code.replace("__SNK_1__", str(output.volcano_png))
        r_code = r_code.replace("__SNK_2__", str(input.de_results))
        r_code = r_code.replace("__SNK_3__", str(params.pvalue))
        r_code = r_code.replace("__SNK_4__", str(params.logFC))
        r_code = r_code.replace("__SNK_5__", str(params.ptype))
        result = subprocess.run(["Rscript", "-"], input=r_code, text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- 差异基因聚类热图 ----
rule deg_heatmap:
    """差异基因表达聚类热图 (Top 50)"""
    input:
        de_results = f"{RESULT_DIR}/06.de_enrich/{{comp}}/{{comp}}.DE_results.All.xls",
        norm_mat   = rules.normalize_expression.output.norm_mat,
        sample_desc = rules.make_sample_desc.output.sample_desc,
    output:
        heatmap_pdf = f"{RESULT_DIR}/06.de_enrich/{{comp}}/{{comp}}.heatmap.pdf",
        heatmap_png = f"{RESULT_DIR}/06.de_enrich/{{comp}}/{{comp}}.heatmap.png",
    params:
        logFC  = config["diff"]["logFC"],
        ptype  = config["diff"]["pOrFDR"],
        pvalue = config["diff"]["diff_pvalue"],
        top_n  = 50,
    conda:
        f"{PIPELINE_DIR}/envs/diff.yaml"
    run:
        import os, subprocess, sys
        os.makedirs(f"{RESULT_DIR}/06.de_enrich/{wildcards.comp}", exist_ok=True)
        r_code = """
        library(pheatmap)
        library(dplyr)
        
        df <- read.table("__SNK_3__", header=TRUE, sep="\\t", check.names=FALSE)
        expr <- read.table("__SNK_4__", header=TRUE, row.names=1, sep="\\t", check.names=FALSE)
        samples <- read.table("__SNK_2__", header=TRUE, sep="\\t", comment.char="", check.names=FALSE)
        
        logfc_cut <- as.numeric("__SNK_8__")
        p_cut <- as.numeric("__SNK_6__")
        p_type <- "__SNK_9__"
        top_n <- as.integer("__SNK_7__")
        
        # 筛选显著差异基因，按显著性排序取 top N
        sig_df <- df %>%
            filter(diffType != "No Significant Difference") %>%
            filter(!is.na(!!sym(p_type))) %>%
            arrange(!!sym(p_type))
        
        if(nrow(sig_df) > top_n) {
            sig_df <- head(sig_df, top_n)
        }
        if(nrow(sig_df) < 2) {
            cat("WARNING: Too few significant DEGs for heatmap (", nrow(sig_df), ")\\n")
            pdf("__SNK_1__", width=10, height=6)
            plot.new()
            text(0.5, 0.5, paste0("Not enough DEGs for heatmap\\n(", nrow(sig_df), " significant genes)"))
            dev.off()
            png("__SNK_0__", width=1200, height=720)
            plot.new()
            text(0.5, 0.5, paste0("Not enough DEGs for heatmap\\n(", nrow(sig_df), " significant genes)"))
            dev.off()
            quit(status=0)
        }
        
        # 提取差异基因的表达矩阵
        deg_ids <- sig_df$gene_id
        heat_data <- expr[rownames(expr) %in% deg_ids, , drop=FALSE]
        
        # Z-score 标准化 (按行)
        heat_data <- t(scale(t(heat_data)))
        heat_data <- heat_data[complete.cases(heat_data), , drop=FALSE]
        
        # 样本注释
        rownames(samples) <- samples$sample
        common_samples <- intersect(colnames(heat_data), samples$sample)
        heat_data <- heat_data[, common_samples, drop=FALSE]
        ann_col <- data.frame(Group = samples[common_samples, "group"])
        rownames(ann_col) <- common_samples
        if("batch" %in% colnames(samples)) {
            ann_col$Batch <- factor(samples[common_samples, "batch"])
        }
        
        cat("Heatmap:", nrow(heat_data), "genes x", ncol(heat_data), "samples\\n")
        
        pdf("__SNK_1__", width=10, height=max(6, nrow(heat_data)*0.2))
        pheatmap(heat_data,
            annotation_col = ann_col,
            show_rownames = (nrow(heat_data) <= 50),
            fontsize_row = 7,
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            scale = "none",
            color = colorRampPalette(c("#2980B9", "white", "#E74C3C"))(100),
            main = paste0("Top ", nrow(heat_data), " DEGs - Z-score"))
        dev.off()
        
        png("__SNK_0__", width=1200, height=max(720, nrow(heat_data)*20))
        pheatmap(heat_data,
            annotation_col = ann_col,
            show_rownames = (nrow(heat_data) <= 50),
            fontsize_row = 7,
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            scale = "none",
            color = colorRampPalette(c("#2980B9", "white", "#E74C3C"))(100),
            main = paste0("Top ", nrow(heat_data), " DEGs - Z-score"))
        dev.off()
        
        """
        r_code = r_code.replace("__SNK_0__", str(output.heatmap_png))
        r_code = r_code.replace("__SNK_1__", str(output.heatmap_pdf))
        r_code = r_code.replace("__SNK_2__", str(input.sample_desc))
        r_code = r_code.replace("__SNK_3__", str(input.de_results))
        r_code = r_code.replace("__SNK_4__", str(input.norm_mat))
        r_code = r_code.replace("__SNK_5__", str(wildcards.comp))
        r_code = r_code.replace("__SNK_6__", str(params.pvalue))
        r_code = r_code.replace("__SNK_7__", str(params.top_n))
        r_code = r_code.replace("__SNK_8__", str(params.logFC))
        r_code = r_code.replace("__SNK_9__", str(params.ptype))
        r_code = r_code.replace("__SNK_10__", str(RESULT_DIR))
        result = subprocess.run(["Rscript", "-"], input=r_code, text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- DE 特征统计 (所有比较组汇总) ----
rule de_feature_counts:
    """汇总所有比较组的差异基因数量统计"""
    input:
        sig_results = expand(f"{RESULT_DIR}/06.de_enrich/{{comp}}/{{comp}}.DE_results.Significant.xls", comp=COMPARISONS),
    output:
        pdf = f"{RESULT_DIR}/06.de_enrich/DE_feature_counts.pdf",
        png = f"{RESULT_DIR}/06.de_enrich/DE_feature_counts.png",
        xls = f"{RESULT_DIR}/06.de_enrich/DE_feature_counts.xls",
    conda:
        f"{PIPELINE_DIR}/envs/diff.yaml"
    run:
        import os, subprocess, sys
        r_code = """
        library(ggplot2)
        library(dplyr)
        all_data <- list()
        sig_files <- unlist(strsplit("__SNK_0__", " "))
        for(f in sig_files) {
            comp_name <- basename(dirname(f))
            df <- tryCatch(read.table(f, header=TRUE, sep="\\t", check.names=FALSE), error=function(e) NULL)
            if(is.null(df) || nrow(df)==0) {
                all_data[[length(all_data)+1]] <- data.frame(Comparison=comp_name, Direction=c("Up","Down"), Count=c(0,0))
            } else {
                up_n <- sum(df$diffType == "Up", na.rm=TRUE)
                down_n <- sum(df$diffType == "Down", na.rm=TRUE)
                all_data[[length(all_data)+1]] <- data.frame(Comparison=comp_name, Direction=c("Up","Down"), Count=c(up_n, down_n))
            }
        }
        if(length(all_data)==0) {
            pdf("__SNK_3__"); plot.new(); text(0.5,0.5,"No DEG data"); dev.off()
            png("__SNK_2__"); plot.new(); text(0.5,0.5,"No DEG data"); dev.off()
            quit(save="no")
        }
        plot_df <- bind_rows(all_data)
        write.table(plot_df, "__SNK_1__", sep="\\t", row.names=FALSE, quote=FALSE)
        p <- ggplot(plot_df, aes(x=Comparison, y=Count, fill=Direction)) +
            geom_bar(stat="identity", position="dodge") +
            geom_text(aes(label=Count), position=position_dodge(0.9), vjust=-0.2, size=3) +
            scale_fill_manual(values=c(Up="#E74C3C", Down="#2980B9")) +
            labs(x="", y="Number of DEGs", title="Differential Expression Gene Counts") +
            theme_bw() + theme(axis.text.x=element_text(angle=30, hjust=1), plot.title=element_text(hjust=0.5))
        w <- max(6, length(unique(plot_df$Comparison))*1.5)
        ggsave("__SNK_3__", p, width=w, height=6)
        ggsave("__SNK_2__", p, width=w, height=6, dpi=150)
        
        """
        r_code = r_code.replace("__SNK_0__", str(input.sig_results))
        r_code = r_code.replace("__SNK_1__", str(output.xls))
        r_code = r_code.replace("__SNK_2__", str(output.png))
        r_code = r_code.replace("__SNK_3__", str(output.pdf))
        result = subprocess.run(["Rscript", "-"], input=r_code, text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- Venn/Upset 图 (多比较组) ----
rule venn_upset:
    """多比较组间差异基因 Venn/Upset 图"""
    input:
        sig_results = expand(f"{RESULT_DIR}/06.de_enrich/{{comp}}/{{comp}}.DE_results.Significant.xls", comp=COMPARISONS),
    output:
        venn_pdf  = f"{RESULT_DIR}/06.de_enrich/venn_upset/venn.pdf",
        venn_png  = f"{RESULT_DIR}/06.de_enrich/venn_upset/venn.png",
        upset_pdf = f"{RESULT_DIR}/06.de_enrich/venn_upset/upset.pdf",
        upset_png = f"{RESULT_DIR}/06.de_enrich/venn_upset/upset.png",
        venn_txt  = f"{RESULT_DIR}/06.de_enrich/venn_upset/venn.txt",
    conda:
        f"{PIPELINE_DIR}/envs/diff.yaml"
    run:
        import os, subprocess, sys
        os.makedirs(f"{RESULT_DIR}/06.de_enrich/venn_upset", exist_ok=True)
        r_code = """
        library(dplyr)
        sig_files <- unlist(strsplit("__SNK_0__", " "))
        deg_sets <- list()
        for(f in sig_files) {
            comp_name <- basename(dirname(f))
            df <- tryCatch(read.table(f, header=TRUE, sep="\\t", check.names=FALSE), error=function(e) NULL)
            if(!is.null(df) && nrow(df)>0) deg_sets[[comp_name]] <- df$gene_id[df$diffType != "No Significant Difference"]
            else deg_sets[[comp_name]] <- character(0)
        }
        n_sets <- length(deg_sets)
        if(n_sets < 2) {
            pdf("__SNK_4__"); plot.new(); text(0.5,0.5,"Need >= 2 comparisons"); dev.off()
            pdf("__SNK_2__"); plot.new(); text(0.5,0.5,"Need >= 2 comparisons"); dev.off()
            quit(save="no")
        }
        # Upset plot
        library(UpSetR)
        all_genes <- unique(unlist(deg_sets))
        upset_mat <- as.data.frame(sapply(deg_sets, function(s) as.integer(all_genes %in% s)))
        rownames(upset_mat) <- all_genes
        if(nrow(upset_mat)>0 && ncol(upset_mat)>=2) {
            pdf("__SNK_2__", width=10, height=6)
            upset(upset_mat, nsets=ncol(upset_mat), order.by="freq", main.bar.color="#3498db", sets.bar.color="#e74c3c")
            dev.off()
            png("__SNK_1__", width=1200, height=720)
            upset(upset_mat, nsets=ncol(upset_mat), order.by="freq", main.bar.color="#3498db", sets.bar.color="#e74c3c")
            dev.off()
        }
        # Venn (max 5)
        if(n_sets <= 5) {
            library(VennDiagram)
            vc <- c("#3498db","#e74c3c","#2ecc71","#f39c12","#9b59b6")[1:n_sets]
            venn_list <- lapply(deg_sets, unique)
            pdf("__SNK_4__", width=8, height=8)
            futile.logger::flog.threshold(futile.logger::ERROR)
            grid.newpage()
            vp <- venn.diagram(venn_list, filename=NULL, fill=vc, category.names=names(venn_list), cat.cex=1.2, main="DEG Overlap")
            grid.draw(vp)
            dev.off()
            png("__SNK_5__", width=800, height=800)
            grid.newpage(); grid.draw(vp); dev.off()
        }
        sink("__SNK_3__")
        for(nm in names(deg_sets)) cat(nm, ":", length(deg_sets[[nm]]), "genes\\n")
        sink()
        
        """
        r_code = r_code.replace("__SNK_0__", str(input.sig_results))
        r_code = r_code.replace("__SNK_1__", str(output.upset_png))
        r_code = r_code.replace("__SNK_2__", str(output.upset_pdf))
        r_code = r_code.replace("__SNK_3__", str(output.venn_txt))
        r_code = r_code.replace("__SNK_4__", str(output.venn_pdf))
        r_code = r_code.replace("__SNK_5__", str(output.venn_png))
        r_code = r_code.replace("__SNK_6__", str(RESULT_DIR))
        result = subprocess.run(["Rscript", "-"], input=r_code, text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- MA 图 ----
rule ma_plot:
    """差异表达 MA 图"""
    input:
        de_results = f"{RESULT_DIR}/06.de_enrich/{{comp}}/{{comp}}.DE_results.All.xls",
    output:
        ma_pdf = f"{RESULT_DIR}/06.de_enrich/{{comp}}/{{comp}}_MA.pdf",
        ma_png = f"{RESULT_DIR}/06.de_enrich/{{comp}}/{{comp}}_MA.png",
    params:
        logFC  = config["diff"]["logFC"],
        pvalue = config["diff"]["diff_pvalue"],
        ptype  = config["diff"]["pOrFDR"],
    conda:
        f"{PIPELINE_DIR}/envs/diff.yaml"
    run:
        import os, subprocess, sys
        r_code = """
        library(ggplot2)
        df <- read.table("__SNK_0__", header=TRUE, sep="\\t", check.names=FALSE)
        logfc_cut <- as.numeric("__SNK_5__")
        p_cut <- as.numeric("__SNK_3__")
        df$log10baseMean <- log10(df$baseMean + 1)
        p <- ggplot(df, aes(x=log10baseMean, y=log2FoldChange, color=diffType)) +
            geom_point(size=1, alpha=0.5) +
            scale_color_manual(values=c(Up="#E74C3C", Down="#2980B9", "No Significant Difference"="grey70")) +
            geom_hline(yintercept=c(-logfc_cut, logfc_cut), linetype="dashed", alpha=0.5) +
            labs(x="log10(baseMean + 1)", y="log2(Fold Change)", title=paste("MA Plot - __SNK_1__")) +
            theme_bw() + theme(legend.position="bottom", plot.title=element_text(hjust=0.5))
        ggsave("__SNK_2__", p, width=8, height=7)
        ggsave("__SNK_4__", p, width=8, height=7, dpi=150)
        
        """
        r_code = r_code.replace("__SNK_0__", str(input.de_results))
        r_code = r_code.replace("__SNK_1__", str(wildcards.comp))
        r_code = r_code.replace("__SNK_2__", str(output.ma_pdf))
        r_code = r_code.replace("__SNK_3__", str(params.pvalue))
        r_code = r_code.replace("__SNK_4__", str(output.ma_png))
        r_code = r_code.replace("__SNK_5__", str(params.logFC))
        result = subprocess.run(["Rscript", "-"], input=r_code, text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")
