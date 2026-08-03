# ============================================================
# 富集分析模块: GO / KEGG 富集 + GSEA
# 支持多物种 OrgDb 动态选择
# ============================================================

# ---- 物种 → OrgDb / KEGG code 映射 ----
# 可根据需要扩展
ORGDB_MAP = {
    "arabidopsis_thaliana": {"orgdb": "org.At.tair.db", "keytype": "TAIR", "kegg": "ath"},
    "homo_sapiens":         {"orgdb": "org.Hs.eg.db",   "keytype": "ENSEMBL", "kegg": "hsa"},
    "mus_musculus":         {"orgdb": "org.Mm.eg.db",   "keytype": "ENSEMBL", "kegg": "mmu"},
    "rattus_norvegicus":    {"orgdb": "org.Rn.eg.db",   "keytype": "ENSEMBL", "kegg": "rno"},
    "danio_rerio":          {"orgdb": "org.Dr.eg.db",   "keytype": "ENSEMBL", "kegg": "dre"},
    "drosophila_melanogaster": {"orgdb": "org.Dm.eg.db", "keytype": "ENSEMBL", "kegg": "dme"},
    "saccharomyces_cerevisiae": {"orgdb": "org.Sc.sgd.db", "keytype": "ENSEMBL", "kegg": "sce"},
    "caenorhabditis_elegans":  {"orgdb": "org.Ce.eg.db", "keytype": "ENSEMBL", "kegg": "cel"},
    "oryza_sativa":         {"orgdb": "org.At.tair.db", "keytype": "TAIR", "kegg": "osa"},  # 水稻用 osa
    "zea_mays":             {"orgdb": "org.At.tair.db", "keytype": "TAIR", "kegg": "zma"},
    "glycine_max":          {"orgdb": "org.At.tair.db", "keytype": "TAIR", "kegg": "gmx"},
    "solanum_lycopersicum": {"orgdb": "org.At.tair.db", "keytype": "TAIR", "kegg": "sly"},
}


# ---- GO 富集 ----
def _get_go_t2g(wildcards):
    """如果 eggNOG 已配置则返回注释文件路径"""
    if config["paths"].get("eggnog_db"):
        return f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/go/new_transcripts.go.xls"
    return []

def _get_kegg_g2k(wildcards):
    """如果 eggNOG 已配置则返回 KO 注释文件路径"""
    if config["paths"].get("eggnog_db"):
        return f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/kegg/new_transcripts.id2K_ko.xls"
    return []

rule go_enrich:
    """GO 富集分析 — 优先用 OrgDb, 无 OrgDb 则用 eggNOG 注释"""
    input:
        up_ids   = f"{RESULT_DIR}/06.de_enrich/{{comp}}/{{comp}}.Significant.Up.id_list",
        down_ids = f"{RESULT_DIR}/06.de_enrich/{{comp}}/{{comp}}.Significant.Down.id_list",
        all_ids  = f"{RESULT_DIR}/06.de_enrich/{{comp}}/{{comp}}.DE_results.All.xls",
        go_t2g   = _get_go_t2g,
    output:
        enrich_go    = f"{RESULT_DIR}/06.de_enrich/{{comp}}/enrich/go/go_enrich.xls",
        go_bar       = f"{RESULT_DIR}/06.de_enrich/{{comp}}/enrich/go/go_barplot.pdf",
        go_dot       = f"{RESULT_DIR}/06.de_enrich/{{comp}}/enrich/go/go_dotplot.pdf",
        go_bp_result = f"{RESULT_DIR}/06.de_enrich/{{comp}}/enrich/go/BP.result.xls",
        go_cc_result = f"{RESULT_DIR}/06.de_enrich/{{comp}}/enrich/go/CC.result.xls",
        go_mf_result = f"{RESULT_DIR}/06.de_enrich/{{comp}}/enrich/go/MF.result.xls",
        go_bp_top10  = f"{RESULT_DIR}/06.de_enrich/{{comp}}/enrich/go/BP.top10.pdf",
        go_cc_top10  = f"{RESULT_DIR}/06.de_enrich/{{comp}}/enrich/go/CC.top10.pdf",
        go_mf_top10  = f"{RESULT_DIR}/06.de_enrich/{{comp}}/enrich/go/MF.top10.pdf",
    params:
        species  = config["project"]["species"],
        kingdom  = config["project"]["kingdom"],
        outdir   = f"{RESULT_DIR}/06.de_enrich/{{comp}}/enrich/go",
        orgdb    = ORGDB_MAP.get(
            config["project"]["species"].lower(), {}
        ).get("orgdb"),
        keytype  = ORGDB_MAP.get(
            config["project"]["species"].lower(), {}
        ).get("keytype", "SYMBOL"),
        kegg_org = ORGDB_MAP.get(
            config["project"]["species"].lower(), {}
        ).get("kegg"),
    conda:
        f"{PIPELINE_DIR}/envs/diff.yaml"
    run:
        import os, subprocess, sys
        os.makedirs(f"{params.outdir}", exist_ok=True)
        r_code = """
        library(clusterProfiler)
        library(ggplot2)
        library(dplyr)
        
        # ---- 读数据 ----
        up_genes <- readLines("__SNK_14__")
        down_genes <- readLines("__SNK_7__")
        all_df <- read.table("__SNK_10__", header=TRUE, sep="\\t", check.names=FALSE)
        all_genes <- unique(all_df$gene_id)
        diff_genes <- c(up_genes, down_genes)
        
        cat("Up genes:", length(up_genes), "\\n")
        cat("Down genes:", length(down_genes), "\\n")
        cat("Total diff genes:", length(diff_genes), "\\n")
        
        species <- "__SNK_9__"
        cat("Species:", species, "\\n")
        
        # ---- 自动选择 OrgDb ----
        orgdb_name <- "__SNK_15__"
        keytype <- "__SNK_8__"
        
        has_orgdb <- FALSE
        if(orgdb_name != "None" && orgdb_name != "") {
            # 尝试加载 OrgDb (失败则自动安装)
            tryCatch({
                library(orgdb_name, character.only=TRUE)
                orgdb <- get(orgdb_name)
                has_orgdb <- TRUE
                cat("Using OrgDb:", orgdb_name, "with keyType:", keytype, "\\n")
            }, error = function(e) {
                cat("OrgDb not found, auto-installing", orgdb_name, "via BiocManager...\\n")
                tryCatch({
                    BiocManager::install(orgdb_name, update=FALSE, ask=FALSE)
                    library(orgdb_name, character.only=TRUE)
                    orgdb <- get(orgdb_name)
                    has_orgdb <<- TRUE
                    cat("Installed and loaded:", orgdb_name, "\\n")
                }, error = function(e2) {
                    cat("Auto-install failed:", e2$message, "\\n")
                    cat("Please install", orgdb_name, "manually or provide GO annotation file.\\n")
                })
            })
        }
        
        # ---- 检查 eggNOG GO 注释 ----
        go_t2g_file <- "__SNK_16__"
        has_go_anno <- file.exists(go_t2g_file) &&
                       file.info(go_t2g_file)$size > 80 &&
                       !grepl("^# (No GO|Placeholder)", readLines(go_t2g_file, n=1))
        
        if(has_go_anno) {
            cat("Using eggNOG-derived GO annotation:", go_t2g_file, "\\n")
            go_term2gene <- read.table(go_t2g_file, header=TRUE, sep="\\t", stringsAsFactors=FALSE)
            cat("  GO TERM2GENE entries:", nrow(go_term2gene), "\\n")
        }
        
        # ---- GO 富集 ----
        run_go_enrich <- function(genes, type_label, ont="BP") {
            if(length(genes) == 0) return(NULL)
        
            # 方法1: enrichGO + OrgDb (模式物种)
            if(has_orgdb) {
                tryCatch({
                    ego <- enrichGO(
                        gene = genes,
                        OrgDb = orgdb,
                        ont = ont,
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2,
                        keyType = keytype
                    )
                    if(!is.null(ego) && nrow(ego) > 0) {
                        ego_df <- as.data.frame(ego)
                        ego_df$type <- type_label
                        ego_df$ontology <- ont
                        return(ego_df)
                    }
                }, error = function(e) {
                    cat("  enrichGO error for", type_label, ont, ":", e$message, "\\n")
                })
            }
        
            # 方法2: enricher + eggNOG TERM2GENE (非模式物种)
            if(!has_orgdb && has_go_anno) {
                tryCatch({
                    ego <- enricher(
                        gene = genes,
                        TERM2GENE = go_term2gene,
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2
                    )
                    if(!is.null(ego) && nrow(ego) > 0) {
                        ego_df <- as.data.frame(ego)
                        ego_df$type <- type_label
                        ego_df$ontology <- ont
                        cat("  GO (enricher) for", type_label, ont, ":", nrow(ego_df), "terms\\n")
                        return(ego_df)
                    }
                }, error = function(e) {
                    cat("  GO enricher error for", type_label, ont, ":", e$message, "\\n")
                })
            }
        
            if(!has_orgdb && !has_go_anno) {
                cat("  No GO annotation available — need OrgDb or eggNOG DB\\n")
            }
            return(NULL)
        }
        
        # Run BP for main result (backward compatible)
        ego_up_df <- run_go_enrich(up_genes, "Up", "BP")
        ego_down_df <- run_go_enrich(down_genes, "Down", "BP")
        ego_all <- rbind(ego_up_df, ego_down_df)
        
        # Additional CC + MF enrichments (all DEGs combined)
        ego_bp <- rbind(run_go_enrich(up_genes,"Up","BP"), run_go_enrich(down_genes,"Down","BP"))
        ego_cc <- rbind(run_go_enrich(up_genes,"Up","CC"), run_go_enrich(down_genes,"Down","CC"))
        ego_mf <- rbind(run_go_enrich(up_genes,"Up","MF"), run_go_enrich(down_genes,"Down","MF"))
        
        # ---- 输出 ----
        if(!is.null(ego_all) && nrow(ego_all) > 0) {
            write.table(ego_all, "__SNK_6__", sep="\\t", row.names=FALSE, quote=FALSE)
            cat("Total enriched GO terms:", nrow(ego_all), "\\n")
        
            # 柱状图
            ego_top <- ego_all %>%
                group_by(type) %>%
                slice_min(p.adjust, n=10) %>%
                ungroup()
            p <- ggplot(ego_top, aes(x=reorder(Description, -log10(p.adjust)),
                                     y=-log10(p.adjust), fill=type)) +
                geom_bar(stat="identity") +
                coord_flip() +
                labs(x="", y="-log10(p.adjust)", title="GO Enrichment (BP)") +
                theme_bw()
            ggsave("__SNK_13__", p, width=10, height=max(6, nrow(ego_top)*0.25))
        
            # 气泡图
            p2 <- ggplot(ego_all, aes(x=type, y=Description, size=Count, color=p.adjust)) +
                geom_point() +
                scale_color_gradient(low="red", high="blue") +
                labs(title="GO Enrichment Dotplot") +
                theme_bw()
            ggsave("__SNK_12__", p2, width=12, height=max(6, nrow(ego_all)*0.2))
        
            # ---- BP/CC/MF 单独输出 ----
            write_ontology <- function(ego_df, result_file, top10_file, ont_label) {
                if(!is.null(ego_df) && nrow(ego_df) > 0) {
                    write.table(ego_df, result_file, sep="\\t", row.names=FALSE, quote=FALSE)
                    ego_top10 <- ego_df %>% slice_min(p.adjust, n=10)
                    p_top <- ggplot(ego_top10, aes(x=reorder(Description, -log10(p.adjust)),
                                                   y=-log10(p.adjust), fill=type)) +
                        geom_bar(stat="identity") + coord_flip() +
                        labs(x="", y="-log10(p.adjust)", title=paste0("GO ", ont_label, " Top10")) + theme_bw()
                    ggsave(top10_file, p_top, width=10, height=max(4, nrow(ego_top10)*0.3))
                    cat("  GO", ont_label, ": wrote", nrow(ego_df), "terms\\n")
                } else {
                    write.table(data.frame(Note=paste("No GO", ont_label, "terms")), result_file, sep="\\t", row.names=FALSE)
                    pdf(top10_file); plot.new(); text(0.5,0.5,paste("No GO", ont_label, "terms")); dev.off()
                }
            }
            write_ontology(ego_bp, "__SNK_1__", "__SNK_4__", "BP")
            write_ontology(ego_cc, "__SNK_2__", "__SNK_3__", "CC")
            write_ontology(ego_mf, "__SNK_0__", "__SNK_5__", "MF")
        } else {
            msg <- data.frame(Note = c(
                "No enriched GO terms found",
                paste("Species:", species),
                paste("OrgDb:", orgdb_name),
                paste("Up genes:", length(up_genes)),
                paste("Down genes:", length(down_genes)),
                "Tip: For non-model species, provide GO annotation file (go_db in config)"
            ))
            write.table(msg, "__SNK_6__", sep="\\t", row.names=FALSE, quote=FALSE)
            pdf("__SNK_13__"); plot.new(); text(0.5, 0.5, "No enriched GO terms", cex=1.5); dev.off()
            pdf("__SNK_12__"); plot.new(); text(0.5, 0.5, "No enriched GO terms", cex=1.5); dev.off()
            cat("WARNING: No enriched GO terms. For non-model species, use a GO annotation file.\\n")
        }
        
        """
        r_code = r_code.replace("__SNK_0__", str(output.go_mf_result))
        r_code = r_code.replace("__SNK_1__", str(output.go_bp_result))
        r_code = r_code.replace("__SNK_2__", str(output.go_cc_result))
        r_code = r_code.replace("__SNK_3__", str(output.go_cc_top10))
        r_code = r_code.replace("__SNK_4__", str(output.go_bp_top10))
        r_code = r_code.replace("__SNK_5__", str(output.go_mf_top10))
        r_code = r_code.replace("__SNK_6__", str(output.enrich_go))
        r_code = r_code.replace("__SNK_7__", str(input.down_ids))
        r_code = r_code.replace("__SNK_8__", str(params.keytype))
        r_code = r_code.replace("__SNK_9__", str(params.species))
        r_code = r_code.replace("__SNK_10__", str(input.all_ids))
        r_code = r_code.replace("__SNK_11__", str(params.outdir))
        r_code = r_code.replace("__SNK_12__", str(output.go_dot))
        r_code = r_code.replace("__SNK_13__", str(output.go_bar))
        r_code = r_code.replace("__SNK_14__", str(input.up_ids))
        r_code = r_code.replace("__SNK_15__", str(params.orgdb))
        r_code = r_code.replace("__SNK_16__", str(input.go_t2g))
        result = subprocess.run(["Rscript", "-"], input=r_code, text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- KEGG 富集 ----
rule kegg_enrich:
    """KEGG 通路富集分析 — 优先用 enrichKEGG, 无物种代码则用 eggNOG KO 注释"""
    input:
        up_ids    = f"{RESULT_DIR}/06.de_enrich/{{comp}}/{{comp}}.Significant.Up.id_list",
        down_ids  = f"{RESULT_DIR}/06.de_enrich/{{comp}}/{{comp}}.Significant.Down.id_list",
        kegg_g2k  = _get_kegg_g2k,
    output:
        enrich_kegg  = f"{RESULT_DIR}/06.de_enrich/{{comp}}/enrich/kegg/kegg_enrich.xls",
        kegg_bar     = f"{RESULT_DIR}/06.de_enrich/{{comp}}/enrich/kegg/kegg_barplot.pdf",
    params:
        kegg_org = config["enrich"].get("kegg_org") or \
            ORGDB_MAP.get(config["project"]["species"].lower(), {}).get("kegg", "ath"),
        outdir   = f"{RESULT_DIR}/06.de_enrich/{{comp}}/enrich/kegg",
    conda:
        f"{PIPELINE_DIR}/envs/diff.yaml"
    run:
        import os, subprocess, sys
        os.makedirs(f"{params.outdir}", exist_ok=True)
        r_code = """
        library(clusterProfiler)
        library(ggplot2)
        library(dplyr)
        
        up_genes <- readLines("__SNK_6__")
        down_genes <- readLines("__SNK_3__")
        
        # ---- 检查 eggNOG KEGG KO 注释 ----
        kegg_g2k_file <- "__SNK_4__"
        has_kegg_anno <- file.exists(kegg_g2k_file) &&
                         file.info(kegg_g2k_file)$size > 80 &&
                         !grepl("^# (No KEGG|Placeholder)", readLines(kegg_g2k_file, n=1))
        
        if(has_kegg_anno) {
            cat("Using eggNOG-derived KEGG KO annotation:", kegg_g2k_file, "\\n")
            kegg_gene2ko <- read.table(kegg_g2k_file, header=TRUE, sep="\\t", stringsAsFactors=FALSE)
            cat("  KEGG KO entries:", nrow(kegg_gene2ko), "\\n")
        }
        
        cat("KEGG organism code:", "__SNK_1__", "\\n")
        cat("Up genes:", length(up_genes), "\\n")
        cat("Down genes:", length(down_genes), "\\n")
        
        run_kegg_enrich <- function(genes, type_label) {
            if(length(genes) == 0) return(NULL)
        
            # 方法1: enrichKEGG (有物种代码)
            has_kegg_code <- "__SNK_1__" != "null" && "__SNK_1__" != ""
            if(has_kegg_code) {
                tryCatch({
                    ekegg <- enrichKEGG(
                        gene = genes,
                        organism = "__SNK_1__",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2
                    )
                    if(!is.null(ekegg) && nrow(ekegg) > 0) {
                        ekegg_df <- as.data.frame(ekegg)
                        ekegg_df$type <- type_label
                        return(ekegg_df)
                    }
                }, error = function(e) {
                    cat("  enrichKEGG error for", type_label, ":", e$message, "\\n")
                })
            }
        
            # 方法2: enricher + eggNOG KO→gene (无物种代码但有 eggNOG 注释)
            if(!has_kegg_code && has_kegg_anno) {
                tryCatch({
                    # 用 KO → gene 做富集
                    ekegg <- enricher(
                        gene = genes,
                        TERM2GENE = kegg_gene2ko,
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2
                    )
                    if(!is.null(ekegg) && nrow(ekegg) > 0) {
                        ekegg_df <- as.data.frame(ekegg)
                        ekegg_df$type <- type_label
                        cat("  KEGG KO (enricher) for", type_label, ":", nrow(ekegg_df), "terms\\n")
                        return(ekegg_df)
                    }
                }, error = function(e) {
                    cat("  KEGG enricher error for", type_label, ":", e$message, "\\n")
                })
            }
        
            if(!has_kegg_code && !has_kegg_anno) {
                cat("  No KEGG annotation available\\n")
            }
            return(NULL)
        }
        
        ekegg_up_df <- run_kegg_enrich(up_genes, "Up")
        ekegg_down_df <- run_kegg_enrich(down_genes, "Down")
        
        ekegg_all <- rbind(ekegg_up_df, ekegg_down_df)
        
        if(!is.null(ekegg_all) && nrow(ekegg_all) > 0) {
            write.table(ekegg_all, "__SNK_0__", sep="\\t", row.names=FALSE, quote=FALSE)
            cat("Total enriched KEGG pathways:", nrow(ekegg_all), "\\n")
        
            ekegg_top <- ekegg_all %>%
                group_by(type) %>%
                slice_min(p.adjust, n=10) %>%
                ungroup()
            p <- ggplot(ekegg_top, aes(x=reorder(Description, -log10(p.adjust)),
                                       y=-log10(p.adjust), fill=type)) +
                geom_bar(stat="identity") +
                coord_flip() +
                labs(x="", y="-log10(p.adjust)", title="KEGG Pathway Enrichment") +
                theme_bw()
            ggsave("__SNK_2__", p, width=10, height=max(6, nrow(ekegg_top)*0.25))
        } else {
            cat("WARNING: No enriched KEGG pathways found\\n")
            msg <- data.frame(Note=c(
                paste("No enriched KEGG pathways"),
                paste("Organism code:", "__SNK_1__"),
                paste("Up genes:", length(up_genes)),
                paste("Down genes:", length(down_genes))
            ))
            write.table(msg, "__SNK_0__", sep="\\t", row.names=FALSE, quote=FALSE)
            pdf("__SNK_2__"); plot.new(); text(0.5, 0.5, "No enriched pathways", cex=1.5); dev.off()
        }
        
        """
        r_code = r_code.replace("__SNK_0__", str(output.enrich_kegg))
        r_code = r_code.replace("__SNK_1__", str(params.kegg_org))
        r_code = r_code.replace("__SNK_2__", str(output.kegg_bar))
        r_code = r_code.replace("__SNK_3__", str(input.down_ids))
        r_code = r_code.replace("__SNK_4__", str(input.kegg_g2k))
        r_code = r_code.replace("__SNK_5__", str(params.outdir))
        r_code = r_code.replace("__SNK_6__", str(input.up_ids))
        result = subprocess.run(["Rscript", "-"], input=r_code, text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- PPI 蛋白互作网络 ----
rule ppi_network:
    """STRING 蛋白互作网络 — 对显著差异基因构建 PPI 网络"""
    input:
        sig_results = f"{RESULT_DIR}/06.de_enrich/{{comp}}/{{comp}}.DE_results.Significant.xls",
    output:
        ppi_edges   = f"{RESULT_DIR}/06.de_enrich/{{comp}}/enrich/ppi/ppi_edges.tsv",
        ppi_nodes   = f"{RESULT_DIR}/06.de_enrich/{{comp}}/enrich/ppi/ppi_nodes.tsv",
        ppi_pdf     = f"{RESULT_DIR}/06.de_enrich/{{comp}}/enrich/ppi/ppi_network.pdf",
    params:
        string_db     = config["paths"].get("string_db", ""),
        score_threshold = config["enrich"].get("ppi_score_threshold", 400),
        outdir        = f"{RESULT_DIR}/06.de_enrich/{{comp}}/enrich/ppi",
    conda:
        f"{PIPELINE_DIR}/envs/diff.yaml"
    run:
        import os, subprocess, sys
        os.makedirs(f"{params.outdir}", exist_ok=True)
        r_code = """
        library(dplyr)
        library(ggplot2)
        library(igraph)
        
        # 读取显著差异基因
        sig_df <- read.table("__SNK_1__", header=TRUE, sep="\\t", check.names=FALSE)
        
        if(nrow(sig_df) == 0) {
            cat("No significant DEGs for PPI analysis\\n")
            write.table(data.frame(Note="No significant DEGs"), "__SNK_4__", sep="\\t", row.names=FALSE)
            write.table(data.frame(Note="No significant DEGs"), "__SNK_3__", sep="\\t", row.names=FALSE)
            pdf("__SNK_6__"); plot.new(); text(0.5,0.5,"No significant DEGs for PPI"); dev.off()
            quit(save="no")
        }
        
        cat("Total significant DEGs:", nrow(sig_df), "\\n")
        cat("  Up:", sum(sig_df$diffType == "Up"), "\\n")
        cat("  Down:", sum(sig_df$diffType == "Down"), "\\n")
        
        string_db <- "__SNK_2__"
        has_string_db <- (string_db != "" && string_db != "null" && file.exists(string_db))
        
        # ---- 如果 STRING 数据库可用 ----
        if(has_string_db) {
            cat("Using STRING database:", string_db, "\\n")
        
            tryCatch({
                string_links <- read.table(string_db, header=TRUE, sep=" ", check.names=FALSE)
        
                # 过滤得分
                score_col <- if("combined_score" %in% colnames(string_links)) "combined_score"
                            else if("score" %in% colnames(string_links)) "score"
                            else NULL
        
                if(!is.null(score_col)) {
                    string_links <- string_links[string_links[[score_col]] >= __SNK_0__, ]
                    cat("STRING interactions after score filter:", nrow(string_links), "\\n")
                }
        
                # 提取 DEG 之间的互作
                deg_ids <- sig_df$gene_id
                deg_links <- string_links[
                    string_links$protein1 %in% deg_ids &
                    string_links$protein2 %in% deg_ids, ]
        
                cat("PPI edges among DEGs:", nrow(deg_links), "\\n")
        
                if(nrow(deg_links) > 0) {
                    write.table(deg_links, "__SNK_4__", sep="\\t", row.names=FALSE, quote=FALSE)
        
                    # 构建网络
                    g <- graph_from_data_frame(deg_links[, c("protein1","protein2")], directed=FALSE)
                    V(g)$degree <- degree(g)
                    top_nodes <- sort(V(g)$degree, decreasing=TRUE)
                    cat("Top hub genes:", paste(names(head(top_nodes, 10)), collapse=", "), "\\n")
        
                    # 保存节点信息
                    node_df <- data.frame(
                        gene = V(g)$name,
                        degree = V(g)$degree,
                        stringsAsFactors = FALSE
                    )
                    # 合并差异表达信息
                    deg_info <- sig_df[, c("gene_id","log2FoldChange","diffType")]
                    colnames(deg_info)[1] <- "gene"
                    node_df <- merge(node_df, deg_info, by="gene", all.x=TRUE)
                    write.table(node_df, "__SNK_3__", sep="\\t", row.names=FALSE, quote=FALSE)
        
                    # 简化绘图 (仅画 top 基因)
                    if(vcount(g) > 100) {
                        keep_nodes <- names(sort(V(g)$degree, decreasing=TRUE)[1:100])
                        g <- induced_subgraph(g, keep_nodes)
                    }
        
                    pdf("__SNK_6__", width=12, height=12)
                    set.seed(42)
                    plot(g,
                        vertex.size = log(V(g)$degree + 1) * 2,
                        vertex.label = ifelse(V(g)$degree > quantile(V(g)$degree, 0.9),
                                              V(g)$name, ""),
                        vertex.label.cex = 0.6,
                        vertex.label.color = "black",
                        vertex.color = "lightblue",
                        vertex.frame.color = "grey50",
                        edge.width = 0.3,
                        edge.color = "grey80",
                        layout = layout_with_fr(g),
                        main = paste0("PPI Network\\n(__SNK_5__)"))
                    dev.off()
                } else {
                    cat("No PPI edges found among DEGs\\n")
                    write.table(data.frame(Note="No PPI edges found"), "__SNK_4__", sep="\\t", row.names=FALSE)
                    write.table(data.frame(Note="No PPI edges found"), "__SNK_3__", sep="\\t", row.names=FALSE)
                    pdf("__SNK_6__"); plot.new(); text(0.5,0.5,"No PPI edges among DEGs"); dev.off()
                }
            }, error = function(e) {
                cat("PPI Error:", e$message, "\\n")
                write.table(data.frame(Error=e$message), "__SNK_4__", sep="\\t", row.names=FALSE)
                write.table(data.frame(Error=e$message), "__SNK_3__", sep="\\t", row.names=FALSE)
                pdf("__SNK_6__"); plot.new(); text(0.5,0.5,paste("Error:", substr(e$message,1,80)), cex=0.8); dev.off()
            })
        } else {
            cat("STRING database not available. Use --string_db to specify path.\\n")
            cat("Download from: https://stringdb-downloads.org/\\n")
            write.table(data.frame(Note="STRING database not provided (use --string_db)"), "__SNK_4__", sep="\\t", row.names=FALSE)
            write.table(data.frame(Note="STRING database not provided"), "__SNK_3__", sep="\\t", row.names=FALSE)
            pdf("__SNK_6__"); plot.new(); text(0.5,0.5,"STRING database not available", cex=1); dev.off()
        }
        
        """
        r_code = r_code.replace("__SNK_0__", str(params.score_threshold))
        r_code = r_code.replace("__SNK_1__", str(input.sig_results))
        r_code = r_code.replace("__SNK_2__", str(params.string_db))
        r_code = r_code.replace("__SNK_3__", str(output.ppi_nodes))
        r_code = r_code.replace("__SNK_4__", str(output.ppi_edges))
        r_code = r_code.replace("__SNK_5__", str(wildcards.comp))
        r_code = r_code.replace("__SNK_6__", str(output.ppi_pdf))
        r_code = r_code.replace("__SNK_7__", str(params.outdir))
        result = subprocess.run(["Rscript", "-"], input=r_code, text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- GSEA (基因集富集分析) ----
rule gsea_analysis:
    """GSEA 基因集富集分析"""
    input:
        de_all = f"{RESULT_DIR}/06.de_enrich/{{comp}}/{{comp}}.DE_results.All.xls",
    output:
        gsea_results = f"{RESULT_DIR}/06.de_enrich/{{comp}}/enrich/gsea/gsea_results.xls",
        gsea_plot    = f"{RESULT_DIR}/06.de_enrich/{{comp}}/enrich/gsea/gsea_plot.pdf",
    params:
        kegg_org = config["enrich"].get("kegg_org") or \
            ORGDB_MAP.get(config["project"]["species"].lower(), {}).get("kegg", "ath"),
        outdir   = f"{RESULT_DIR}/06.de_enrich/{{comp}}/enrich/gsea",
    conda:
        f"{PIPELINE_DIR}/envs/diff.yaml"
    run:
        import os, subprocess, sys
        os.makedirs(f"{params.outdir}", exist_ok=True)
        r_code = """
        library(clusterProfiler)
        library(enrichplot)
        library(ggplot2)
        
        de_df <- read.table("__SNK_4__", header=TRUE, sep="\\t", check.names=FALSE)
        
        if(nrow(de_df) == 0) {
            cat("ERROR: Empty DE results\\n")
            write.table(data.frame(Error="Empty DE results"), "__SNK_0__", sep="\\t")
            pdf("__SNK_1__"); plot.new(); text(0.5,0.5,"Empty DE results"); dev.off()
            quit(save="no")
        }
        
        # 创建排序基因列表 (按 log2FC 降序)
        gene_list <- de_df$log2FoldChange
        names(gene_list) <- de_df$gene_id
        gene_list <- gene_list[!is.na(gene_list)]
        gene_list <- sort(gene_list, decreasing = TRUE)
        
        cat("Gene list length:", length(gene_list), "\\n")
        cat("Range:", range(gene_list), "\\n")
        
        tryCatch({
            gsea_res <- gseKEGG(
                geneList = gene_list,
                organism = "__SNK_2__",
                pvalueCutoff = 1.0,
                minGSSize = 10,
                maxGSSize = 500,
                eps = 0
            )
        
            if(!is.null(gsea_res) && nrow(gsea_res) > 0) {
                gsea_df <- as.data.frame(gsea_res)
                write.table(gsea_df, "__SNK_0__", sep="\\t", row.names=FALSE, quote=FALSE)
                cat("GSEA enriched gene sets:", nrow(gsea_df), "\\n")
        
                pdf("__SNK_1__", width=14, height=10)
                if(nrow(gsea_df) >= 5) {
                    top_paths <- head(gsea_df[order(gsea_df$pvalue), "ID"], 5)
                    print(gseaplot2(gsea_res, geneSetID = top_paths, pvalue_table = TRUE))
                } else if(nrow(gsea_df) >= 1) {
                    print(gseaplot2(gsea_res, geneSetID = 1:nrow(gsea_df), pvalue_table = TRUE))
                }
                dev.off()
            } else {
                cat("No enriched gene sets found\\n")
                write.table(data.frame(Note="No enriched gene sets"), "__SNK_0__", sep="\\t")
                pdf("__SNK_1__"); plot.new(); text(0.5,0.5,"No enriched gene sets"); dev.off()
            }
        }, error = function(e) {
            cat("GSEA Error:", e$message, "\\n")
            write.table(data.frame(Error=e$message), "__SNK_0__", sep="\\t")
            pdf("__SNK_1__"); plot.new(); text(0.5,0.5,paste("GSEA Error:", substr(e$message,1,80)), cex=0.8); dev.off()
        })
        
        """
        r_code = r_code.replace("__SNK_0__", str(output.gsea_results))
        r_code = r_code.replace("__SNK_1__", str(output.gsea_plot))
        r_code = r_code.replace("__SNK_2__", str(params.kegg_org))
        r_code = r_code.replace("__SNK_3__", str(params.outdir))
        r_code = r_code.replace("__SNK_4__", str(input.de_all))
        result = subprocess.run(["Rscript", "-"], input=r_code, text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")
