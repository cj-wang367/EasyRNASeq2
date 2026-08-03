# ============================================================
# WGCNA 模块: 加权基因共表达网络分析
# ============================================================

# ---- 数据预处理 ----
rule wgcna_prepare:
    """准备 WGCNA 输入数据"""
    input:
        expr_mat = rules.normalize_expression.output.norm_mat,
        sample_info = config["paths"]["sample_info"],
    output:
        expr_csv   = f"{TMP_DIR}/09.wgcna/expression.csv",
        trait_csv  = f"{TMP_DIR}/09.wgcna/traits.csv",
    conda:
        f"{PIPELINE_DIR}/envs/wgcna.yaml"
    run:
        import os, subprocess, sys
        os.makedirs(f"{TMP_DIR}/09.wgcna", exist_ok=True)
        r_code = """
        import pandas as pd
        
        # 读表达矩阵
        expr = pd.read_csv('__SNK_3__', sep='\\t', index_col=0)
        
        # 过滤低表达基因 (至少在一半样本中表达 > 1)
        keep = (expr > 1).sum(axis=1) >= (expr.shape[1] / 2)
        expr_filtered = expr[keep]
        print(f'Genes after filter: {expr_filtered.shape[0]} / {expr.shape[0]}')
        
        # 转置 (WGCNA 需要 gene x sample)
        expr_filtered.T.to_csv('__SNK_2__')
        
        # 生成 trait 矩阵
        samples = pd.read_csv('__SNK_0__', sep='\\t', comment='#')
        samples = samples[['sample', 'group']]
        
        # One-hot 编码分组
        trait_df = pd.get_dummies(samples.set_index('sample')['group'])
        trait_df.to_csv('__SNK_1__')
        print(f'Traits: {list(trait_df.columns)}')
        
        """
        r_code = r_code.replace("__SNK_0__", str(input.sample_info))
        r_code = r_code.replace("__SNK_1__", str(output.trait_csv))
        r_code = r_code.replace("__SNK_2__", str(output.expr_csv))
        r_code = r_code.replace("__SNK_3__", str(input.expr_mat))
        r_code = r_code.replace("__SNK_4__", str(TMP_DIR))
        result = subprocess.run(["python3", "-c", r_code], text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- WGCNA 主分析 ----
rule wgcna_analysis:
    """WGCNA 共表达网络构建"""
    input:
        expr      = rules.wgcna_prepare.output.expr_csv,
        trait     = rules.wgcna_prepare.output.trait_csv,
    output:
        module_genes    = f"{RESULT_DIR}/09.wgcna/module_genes.xls",
        module_trait    = f"{RESULT_DIR}/09.wgcna/module_trait_cor.xls",
        module_trait_p  = f"{RESULT_DIR}/09.wgcna/module_trait_pvalue.xls",
        dendrogram      = f"{RESULT_DIR}/09.wgcna/sample_dendrogram.pdf",
        power_plot      = f"{RESULT_DIR}/09.wgcna/soft_threshold_power.pdf",
        network_heatmap = f"{RESULT_DIR}/09.wgcna/network_heatmap.pdf",
    params:
        power      = config["wgcna"]["power"],
        min_size   = config["wgcna"]["min_module_size"],
        merge_h    = config["wgcna"]["merge_cut_height"],
        outdir     = f"{RESULT_DIR}/09.wgcna",
    conda:
        f"{PIPELINE_DIR}/envs/wgcna.yaml"
    run:
        import os, subprocess, sys
        os.makedirs(f"{params.outdir}", exist_ok=True)
        r_code = """
        library(WGCNA)
        library(ggplot2)
        
        options(stringsAsFactors = FALSE)
        enableWGCNAThreads(nThreads = 4)
        
        # 读数据
        expr_data <- read.csv("__SNK_14__", row.names=1, check.names=FALSE)
        trait_data <- read.csv("__SNK_13__", row.names=1, check.names=FALSE)
        
        cat("Expression data:", ncol(expr_data), "samples x", nrow(expr_data), "genes\\n")
        cat("Traits:", ncol(trait_data), "\\n")
        
        # 防御: WGCNA 至少需要 4 个样本 (建议 >= 15)
        if(ncol(expr_data) < 4) {
            cat("WARNING: Too few samples (", ncol(expr_data), ") for WGCNA (need >= 4)\\n")
            # 生成占位输出满足 Snakemake
            pdf("__SNK_7__", width=8, height=6)
            plot.new(); text(0.5, 0.5, paste0("Insufficient samples (n=", ncol(expr_data), ") for WGCNA (need >=4)"))
            dev.off()
            write.table(data.frame(Note="Insufficient samples"), "__SNK_8__", sep="\\t")
            write.table(data.frame(Note="Insufficient samples"), "__SNK_4__", sep="\\t")
            write.table(data.frame(Note="Insufficient samples"), "__SNK_2__", sep="\\t")
            write.table(data.frame(Note="Insufficient samples"), "__SNK_5__", sep="\\t")
            pdf("__SNK_3__", width=8, height=6); plot.new(); text(0.5, 0.5, "Insufficient samples"); dev.off()
            pdf("__SNK_0__", width=8, height=6); plot.new(); text(0.5, 0.5, "Insufficient samples"); dev.off()
            quit(status=0)
        }
        
        # 样本聚类
        sample_tree <- hclust(dist(t(expr_data)), method="average")
        
        pdf("__SNK_7__", width=12, height=6)
        plot(sample_tree, main="Sample Clustering Dendrogram", xlab="", sub="")
        dev.off()
        
        # 软阈值选择
        powers <- c(1:30)
        sft <- pickSoftThreshold(t(expr_data),
            powerVector = powers,
            networkType = "signed",
            verbose = 5)
        
        pdf("__SNK_6__", width=10, height=5)
        par(mfrow=c(1,2))
        plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
            xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit",
            main="Scale independence", type="n")
        text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
            labels=powers, col="red")
        abline(h=0.85, col="red", lty=2)
        
        plot(sft$fitIndices[,1], sft$fitIndices[,5],
            xlab="Soft Threshold (power)", ylab="Mean Connectivity",
            main="Mean connectivity", type="n")
        text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, col="red")
        dev.off()
        
        # 确定 power
        power <- as.numeric("__SNK_12__")
        if(is.na(sft$powerEstimate)) {
            cat("No power estimate found, using configured power:", power, "\\n")
        } else {
            cat("Estimated power:", sft$powerEstimate, "\\n")
        }
        
        # 一步法网络构建
        net <- blockwiseModules(
            t(expr_data),
            power = power,
            minModuleSize = as.numeric("__SNK_9__"),
            mergeCutHeight = as.numeric("__SNK_10__"),
            networkType = "signed",
            numericLabels = TRUE,
            pamRespectsDendro = FALSE,
            saveTOMs = TRUE,
            saveTOMFileBase = "__SNK_11__/TOM",
            verbose = 3,
            nThreads = 4
        )
        
        # 模块-性状关联
        ME_list <- net$MEs
        colnames(ME_list) <- gsub("^ME", "", colnames(ME_list))
        MEs <- ME_list
        
        # 确保 trait 顺序与 MEs 一致
        common_samples <- intersect(rownames(MEs), rownames(trait_data))
        MEs <- MEs[common_samples, ]
        trait_data_sub <- trait_data[common_samples, , drop=FALSE]
        
        module_trait_cor <- cor(MEs, trait_data_sub, use="p")
        module_trait_pvalue <- corPvalueStudent(module_trait_cor, nrow(MEs))
        
        write.table(module_trait_cor, "__SNK_4__", sep="\\t", quote=FALSE)
        write.table(module_trait_pvalue, "__SNK_2__", sep="\\t", quote=FALSE)
        
        # 模块-基因对应表
        module_labels <- net$colors
        gene_modules <- data.frame(
            gene = colnames(expr_data),
            module = paste0("M", module_labels),
            stringsAsFactors = FALSE
        )
        write.table(gene_modules, "__SNK_5__", sep="\\t", row.names=FALSE, quote=FALSE)
        
        cat("\\nModule sizes:\\n")
        print(table(gene_modules$module))
        
        # 网络热图
        pdf("__SNK_1__", width=10, height=10)
        labeledHeatmap(
            Matrix = module_trait_cor,
            xLabels = colnames(trait_data_sub),
            yLabels = colnames(MEs),
            ySymbols = colnames(MEs),
            colorLabels = FALSE,
            colors = blueWhiteRed(50),
            textMatrix = round(module_trait_cor, 2),
            setStdMargins = FALSE,
            cex.text = 0.7,
            main = "Module-Trait Relationships"
        )
        dev.off()
        
        cat("\\nWGCNA analysis completed.\\n")
        
        """
        r_code = r_code.replace("__SNK_0__", str(output.module_genes_heatmap))
        r_code = r_code.replace("__SNK_1__", str(output.network_heatmap))
        r_code = r_code.replace("__SNK_2__", str(output.module_trait_p))
        r_code = r_code.replace("__SNK_3__", str(output.module_heatmap))
        r_code = r_code.replace("__SNK_4__", str(output.module_trait))
        r_code = r_code.replace("__SNK_5__", str(output.module_genes))
        r_code = r_code.replace("__SNK_6__", str(output.power_plot))
        r_code = r_code.replace("__SNK_7__", str(output.dendrogram))
        r_code = r_code.replace("__SNK_8__", str(output.power_est))
        r_code = r_code.replace("__SNK_9__", str(params.min_size))
        r_code = r_code.replace("__SNK_10__", str(params.merge_h))
        r_code = r_code.replace("__SNK_11__", str(params.outdir))
        r_code = r_code.replace("__SNK_12__", str(params.power))
        r_code = r_code.replace("__SNK_13__", str(input.trait))
        r_code = r_code.replace("__SNK_14__", str(input.expr))
        result = subprocess.run(["Rscript", "-"], input=r_code, text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")
