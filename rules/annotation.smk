# ============================================================
# 功能注释模块: eggNOG-mapper 一站式注释
#
# eggNOG-mapper 可同时提供:
#   GO / KEGG KO / KEGG Pathway / COG / Pfam / CAZy / EC
#   比 Diamond BLAST + InterProScan 方案更简单高效
#
# 数据库要求:
#   --data_dir 指向 eggNOG 数据库目录, 至少包含:
#     eggnog.db               # 注释 SQLite 数据库
#     eggnog_proteins.dmnd     # Diamond 格式蛋白序列库
#     eggnog.taxa.db           # 分类学数据库
# ============================================================

# ---- eggNOG-mapper 功能注释 ----
rule eggnog_mapper:
    """eggNOG-mapper 一站式功能注释
    输出列:
      query / seed_eggNOG_ortholog / GO_terms / KEGG_ko /
      COG_category / Preferred_name / PFAMs / EC / CAZy ...
    """
    input:
        pep = rules.transdecoder_orf.output.pep,
    output:
        annotations = f"{TMP_DIR}/05.annotation/eggnog.emapper.annotations",
    params:
        outdir   = f"{TMP_DIR}/05.annotation",
        data_dir = config["paths"].get("eggnog_db", ""),
    threads: 8
    conda:
        f"{PIPELINE_DIR}/envs/annotation.yaml"
    shell:
        """
        mkdir -p {params.outdir}

        # 检查输入蛋白序列是否为空 (TransDecoder 可能未预测出任何蛋白)
        if [ ! -s {input.pep} ]; then
            echo "[WARN] No protein sequences found (input .pep is empty or missing)"
            echo "# No protein sequences — eggNOG annotation skipped" > {output.annotations}
            echo "query" > {TMP_DIR}/05.annotation/empty_columns.txt
            exit 0
        fi
        PEPTIDE_COUNT=$(grep -c '^>' {input.pep} 2>/dev/null || echo 0)
        if [ "$PEPTIDE_COUNT" -eq 0 ]; then
            echo "[WARN] Zero protein sequences in .pep file"
            echo "# Zero protein sequences — eggNOG annotation skipped" > {output.annotations}
            exit 0
        fi

        # 检查 emapper.py 可用性
        if ! command -v emapper.py &> /dev/null; then
            echo "[ERROR] emapper.py not found! Please install eggnog-mapper: conda install -c bioconda eggnog-mapper"
            echo "# eggnog-mapper not installed" > {output.annotations}
            exit 0
        fi

        # 检查数据库
        DATA_DIR="{params.data_dir}"
        if [ -z "$DATA_DIR" ] || [ ! -d "$DATA_DIR" ]; then
            echo "[ERROR] eggNOG database not found: $DATA_DIR"
            echo "  Please specify --eggnog_db /path/to/eggnog_data/"
            echo "# eggNOG database not found at $DATA_DIR" > {output.annotations}
            exit 0
        fi
        if [ ! -f "$DATA_DIR/eggnog.db" ]; then
            echo "[ERROR] eggnog.db not found in $DATA_DIR"
            echo "# eggnog.db missing" > {output.annotations}
            exit 0
        fi
        if [ ! -f "$DATA_DIR/eggnog_proteins.dmnd" ]; then
            echo "[ERROR] eggnog_proteins.dmnd not found in $DATA_DIR"
            echo "# Diamond DB missing" > {output.annotations}
            exit 0
        fi

        echo "[INFO] eggNOG database: $DATA_DIR"
        echo "[INFO] Input proteins: $(grep -c '^>' {input.pep}) sequences"
        echo "[INFO] Running eggNOG-mapper (may take 30-120 min)..."

        emapper.py \
            -i {input.pep} \
            --output {params.outdir}/eggnog \
            --data_dir {params.data_dir} \
            --cpu {threads} \
            -m diamond \
            --itype proteins \
            --go_evidence non-electronic \
            --target_orthologs all \
            --tax_scope auto \
            2>&1 | tail -20

        # 检查输出
        if [ -f {params.outdir}/eggnog.emapper.annotations ]; then
            cp {params.outdir}/eggnog.emapper.annotations {output.annotations}
            n_annot=$(tail -n +2 {output.annotations} | wc -l)
            echo "[OK] eggNOG-mapper done: $n_annot sequences annotated"
        else
            echo "[ERROR] eggNOG-mapper produced no output"
            echo "# No annotations produced" > {output.annotations}
            exit 1
        fi
        """


# ---- 提取 GO 注释 (gene → GO mapping) ----
rule extract_go_annotation:
    """从 eggNOG 结果提取基因-GO映射 (TERM2GENE 格式, 用于 clusterProfiler::enricher)"""
    input:
        eggnog = rules.eggnog_mapper.output.annotations,
    output:
        go_term2gene = f"{TMP_DIR}/05.annotation/go_term2gene.txt",
        go_gene2term = f"{TMP_DIR}/05.annotation/go_gene2term.txt",
    conda:
        f"{PIPELINE_DIR}/envs/annotation.yaml"
    run:
        import os, subprocess, sys
        r_code = """
        import pandas as pd
        import sys
        
        try:
            df = pd.read_csv('__SNK_2__', sep='\\t', comment='#')
        
            # 自动检测 GO 列 (列名可能是 GO_terms / GOs / go_terms)
            go_col = None
            for col in df.columns:
                if 'go' in col.lower() and ('term' in col.lower() or col.lower() in ['gos', 'go_terms']):
                    go_col = col
                    break
        
            if not go_col:
                print('# No GO column found in eggNOG output', file=sys.stderr)
                print('# Columns:', list(df.columns), file=sys.stderr)
                open('__SNK_1__', 'w').write('# No GO terms found\\n')
                open('__SNK_0__', 'w').write('# No GO terms found\\n')
                sys.exit(0)
        
            print(f'[INFO] Using GO column: {go_col}', file=sys.stderr)
        
            # 构建 TERM2GENE (GO term → gene_id) 和 GENE2TERM (gene_id → GO term)
            term2gene = []  # GO, gene
            gene2term = []  # gene, GO_terms (comma-separated)
        
            gene_col = df.columns[0]  # query name = first column
        
            for _, row in df.iterrows():
                gene_id = str(row[gene_col])
                go_terms = str(row[go_col]) if pd.notna(row[go_col]) else ''
        
                if not go_terms or go_terms == 'nan' or go_terms == '-':
                    continue
        
                # GO terms 格式: "GO:0005524,GO:0004672"
                terms = [t.strip() for t in go_terms.split(',') if t.strip() and t.strip() != '-']
                if terms:
                    gene2term.append({'gene': gene_id, 'GO': ','.join(terms)})
                    for term in terms:
                        term2gene.append({'GO': term, 'gene': gene_id})
        
            # 写 TERM2GENE — clusterProfiler::enricher() 的标准输入
            if term2gene:
                pd.DataFrame(term2gene).to_csv('__SNK_1__', sep='\\t', index=False)
                print(f'[INFO] GO TERM2GENE: {len(term2gene)} entries', file=sys.stderr)
            else:
                open('__SNK_1__', 'w').write('# No GO terms found\\n')
        
            # 写 GENE2TERM — 备用格式
            if gene2term:
                pd.DataFrame(gene2term).to_csv('__SNK_0__', sep='\\t', index=False)
                print(f'[INFO] GO GENE2TERM: {len(gene2term)} entries', file=sys.stderr)
            else:
                open('__SNK_0__', 'w').write('# No GO terms found\\n')
        
        except Exception as e:
            print(f'# Error extracting GO: {e}', file=sys.stderr)
            open('__SNK_1__', 'w').write(f'# Error: {e}\\n')
            open('__SNK_0__', 'w').write(f'# Error: {e}\\n')
        
        """
        r_code = r_code.replace("__SNK_0__", str(output.go_gene2term))
        r_code = r_code.replace("__SNK_1__", str(output.go_term2gene))
        r_code = r_code.replace("__SNK_2__", str(input.eggnog))
        result = subprocess.run(["python3", "-c", r_code], text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- 提取 KEGG 注释 (gene → KO mapping) ----
rule extract_kegg_annotation:
    """从 eggNOG 结果提取基因-KO映射 (用于 KEGG 富集)"""
    input:
        eggnog = rules.eggnog_mapper.output.annotations,
    output:
        kegg_gene2ko   = f"{TMP_DIR}/05.annotation/kegg_gene2ko.txt",
        kegg_gene2pathway = f"{TMP_DIR}/05.annotation/kegg_gene2pathway.txt",
    conda:
        f"{PIPELINE_DIR}/envs/annotation.yaml"
    run:
        import os, subprocess, sys
        r_code = """
        import pandas as pd
        import sys
        
        try:
            df = pd.read_csv('__SNK_2__', sep='\\t', comment='#')
        
            # 自动检测 KEGG KO 列
            ko_col = None
            for col in df.columns:
                if 'kegg_ko' in col.lower() or col.lower() == 'ko':
                    ko_col = col
                    break
        
            # 自动检测 KEGG Pathway 列
            path_col = None
            for col in df.columns:
                if 'kegg_pathway' in col.lower() or col.lower() == 'pathway':
                    path_col = col
                    break
        
            if not ko_col and not path_col:
                print('# No KEGG KO/Pathway column found', file=sys.stderr)
                open('__SNK_1__', 'w').write('# No KEGG KO found\\n')
                open('__SNK_0__', 'w').write('# No KEGG Pathway found\\n')
                sys.exit(0)
        
            gene_col = df.columns[0]
        
            # Gene → KO
            if ko_col:
                gene2ko = []
                for _, row in df.iterrows():
                    gene_id = str(row[gene_col])
                    ko_val = str(row[ko_col]) if pd.notna(row[ko_col]) else ''
                    if ko_val and ko_val != 'nan' and ko_val != '-':
                        # KO 格式: "ko:K00001,ko:K00002" 或 "K00001"
                        kos = [k.strip().replace('ko:', '') for k in ko_val.split(',') if k.strip() and k.strip() != '-']
                        for ko in kos:
                            if ko:
                                gene2ko.append({'gene': gene_id, 'KO': ko})
        
                if gene2ko:
                    pd.DataFrame(gene2ko).to_csv('__SNK_1__', sep='\\t', index=False)
                    print(f'[INFO] KEGG gene→KO: {len(gene2ko)} entries', file=sys.stderr)
                else:
                    open('__SNK_1__', 'w').write('# No KO assignments\\n')
            else:
                open('__SNK_1__', 'w').write('# No KO column\\n')
        
            # Gene → Pathway
            if path_col:
                gene2path = []
                for _, row in df.iterrows():
                    gene_id = str(row[gene_col])
                    path_val = str(row[path_col]) if pd.notna(row[path_col]) else ''
                    if path_val and path_val != 'nan' and path_val != '-':
                        # Pathway 格式: "ko00010,ko00020" 或 "00010,00020"
                        paths = [p.strip().replace('ko', '') for p in path_val.split(',') if p.strip() and p.strip() != '-']
                        for p in paths:
                            if p:
                                gene2path.append({'gene': gene_id, 'pathway': f'map{p}'})
        
                if gene2path:
                    pd.DataFrame(gene2path).to_csv('__SNK_0__', sep='\\t', index=False)
                    print(f'[INFO] KEGG gene→pathway: {len(gene2path)} entries', file=sys.stderr)
                else:
                    open('__SNK_0__', 'w').write('# No pathway assignments\\n')
            else:
                open('__SNK_0__', 'w').write('# No Pathway column\\n')
        
        except Exception as e:
            print(f'# Error extracting KEGG: {e}', file=sys.stderr)
            open('__SNK_1__', 'w').write(f'# Error: {e}\\n')
            open('__SNK_0__', 'w').write(f'# Error: {e}\\n')
        
        """
        r_code = r_code.replace("__SNK_0__", str(output.kegg_gene2pathway))
        r_code = r_code.replace("__SNK_1__", str(output.kegg_gene2ko))
        r_code = r_code.replace("__SNK_2__", str(input.eggnog))
        result = subprocess.run(["python3", "-c", r_code], text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- 提取 KOG/COG 注释 ----
rule extract_kog_annotation:
    """从 eggNOG 结果提取 COG/KOG 分类 (COG_category 列)"""
    input:
        eggnog = rules.eggnog_mapper.output.annotations,
    output:
        kog_gene2cog = f"{TMP_DIR}/05.annotation/kog_gene2cog.txt",
        kog_cog_stats = f"{TMP_DIR}/05.annotation/kog_cog_stats.xls",
    conda:
        f"{PIPELINE_DIR}/envs/annotation.yaml"
    run:
        import os, subprocess, sys
        r_code = """
        import pandas as pd
        import sys
        
        # COG 分类定义
        COG_DESC = {
            'A': 'RNA processing and modification',
            'B': 'Chromatin structure and dynamics',
            'C': 'Energy production and conversion',
            'D': 'Cell cycle control, cell division',
            'E': 'Amino acid transport and metabolism',
            'F': 'Nucleotide transport and metabolism',
            'G': 'Carbohydrate transport and metabolism',
            'H': 'Coenzyme transport and metabolism',
            'I': 'Lipid transport and metabolism',
            'J': 'Translation, ribosomal structure and biogenesis',
            'K': 'Transcription',
            'L': 'Replication, recombination and repair',
            'M': 'Cell wall/membrane/envelope biogenesis',
            'N': 'Cell motility',
            'O': 'Posttranslational modification, protein turnover',
            'P': 'Inorganic ion transport and metabolism',
            'Q': 'Secondary metabolites biosynthesis',
            'R': 'General function prediction only',
            'S': 'Function unknown',
            'T': 'Signal transduction mechanisms',
            'U': 'Intracellular trafficking, secretion',
            'V': 'Defense mechanisms',
            'W': 'Extracellular structures',
            'X': 'Mobilome: prophages, transposons',
            'Y': 'Nuclear structure',
            'Z': 'Cytoskeleton',
        }
        
        try:
            df = pd.read_csv('__SNK_2__', sep='\\t', comment='#')
            gene_col = df.columns[0]
        
            # 找 COG_category 列
            cog_col = None
            for col in df.columns:
                if 'cog_category' in col.lower() or col.lower() == 'cog':
                    cog_col = col
                    break
        
            if not cog_col:
                print('# No COG_category column found')
                open('__SNK_1__', 'w').write('# No COG info\\n')
                open('__SNK_0__', 'w').write('category\\tdescription\\tcount\\n')
                sys.exit(0)
        
            from collections import Counter
            gene2cog = []
            cat_counter = Counter()
            for _, row in df.iterrows():
                gene_id = str(row[gene_col])
                cog_val = str(row[cog_col]) if pd.notna(row[cog_col]) else ''
                if cog_val and cog_val != 'nan' and cog_val != '-':
                    cats = [c.strip() for c in cog_val if c.strip() and c.strip() != '-']
                    if cats:
                        gene2cog.append({'gene': gene_id, 'COG_category': cog_val})
                        for c in cats:
                            cat_counter[c] += 1
        
            if gene2cog:
                pd.DataFrame(gene2cog).to_csv('__SNK_1__', sep='\\t', index=False)
                print(f'[INFO] KOG gene→COG: {len(gene2cog)} entries')
            else:
                open('__SNK_1__', 'w').write('# No COG assignments\\n')
        
            with open('__SNK_0__', 'w') as f:
                f.write('category\\tdescription\\tcount\\n')
                for cat, count in cat_counter.most_common():
                    desc = COG_DESC.get(cat, 'Unknown')
                    f.write(f'{cat}\\t{desc}\\t{count}\\n')
        
        except Exception as e:
            print(f'# Error extracting COG: {e}')
        
        """
        r_code = r_code.replace("__SNK_0__", str(output.kog_cog_stats))
        r_code = r_code.replace("__SNK_1__", str(output.kog_gene2cog))
        r_code = r_code.replace("__SNK_2__", str(input.eggnog))
        result = subprocess.run(["python3", "-c", r_code], text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- 提取 Pfam 注释 ----
rule extract_pfam_annotation:
    """从 eggNOG 结果提取 Pfam 结构域信息 (PFAMs 列)"""
    input:
        eggnog = rules.eggnog_mapper.output.annotations,
    output:
        pfam_gene2pfam = f"{TMP_DIR}/05.annotation/pfam_gene2pfam.txt",
        pfam_stats     = f"{TMP_DIR}/05.annotation/pfam_stats.xls",
    conda:
        f"{PIPELINE_DIR}/envs/annotation.yaml"
    run:
        import os, subprocess, sys
        r_code = """
        import pandas as pd
        import sys
        
        try:
            df = pd.read_csv('__SNK_2__', sep='\\t', comment='#')
            gene_col = df.columns[0]
        
            # 找 PFAMs 列
            pfam_col = None
            for col in df.columns:
                if 'pfam' in col.lower():
                    pfam_col = col
                    break
        
            if not pfam_col:
                print('# No PFAMs column found')
                open('__SNK_0__', 'w').write('# No Pfam info\\n')
                open('__SNK_1__', 'w').write('domain\\tcount\\n')
                sys.exit(0)
        
            from collections import Counter
            gene2pfam = []
            domain_counter = Counter()
            for _, row in df.iterrows():
                gene_id = str(row[gene_col])
                pfam_val = str(row[pfam_col]) if pd.notna(row[pfam_col]) else ''
                if pfam_val and pfam_val != 'nan' and pfam_val != '-':
                    domains = [d.strip().split('.')[0] for d in pfam_val.split(',') if d.strip() and d.strip() != '-']
                    if domains:
                        gene2pfam.append({'gene': gene_id, 'PFAMs': pfam_val})
                        for d in domains:
                            domain_counter[d] += 1
        
            if gene2pfam:
                pd.DataFrame(gene2pfam).to_csv('__SNK_0__', sep='\\t', index=False)
                print(f'[INFO] Pfam gene→domain: {len(gene2pfam)} entries')
            else:
                open('__SNK_0__', 'w').write('# No Pfam assignments\\n')
        
            with open('__SNK_1__', 'w') as f:
                f.write('domain\\tcount\\n')
                for dom, count in domain_counter.most_common():
                    f.write(f'{dom}\\t{count}\\n')
        
        except Exception as e:
            print(f'# Error extracting Pfam: {e}')
        
        """
        r_code = r_code.replace("__SNK_0__", str(output.pfam_gene2pfam))
        r_code = r_code.replace("__SNK_1__", str(output.pfam_stats))
        r_code = r_code.replace("__SNK_2__", str(input.eggnog))
        result = subprocess.run(["python3", "-c", r_code], text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- 注释可视化 ----
rule annotation_plots:
    """生成 GO/KEGG/KOG/Pfam 注释分类图"""
    input:
        eggnog       = rules.eggnog_mapper.output.annotations,
        go_gene2term = rules.extract_go_annotation.output.go_gene2term,
        kog_stats    = rules.extract_kog_annotation.output.kog_cog_stats,
        pfam_stats   = rules.extract_pfam_annotation.output.pfam_stats,
    output:
        go_plot     = f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/go/go_plot.pdf",
        go_plot_png = f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/go/go_plot.png",
        kog_plot    = f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/kog/kog_plot.pdf",
        kog_plot_png= f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/kog/kog_plot.png",
        pfam_plot   = f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/pfam/pfam_plot.pdf",
        pfam_plot_png=f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/pfam/pfam_plot.png",
    conda:
        f"{PIPELINE_DIR}/envs/diff.yaml"
    run:
        import os, subprocess, sys
        os.makedirs(f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/go", exist_ok=True)
        os.makedirs(f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/kog", exist_ok=True)
        os.makedirs(f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/pfam", exist_ok=True)
        r_code = """
        library(ggplot2)
        library(dplyr)
        
        # --- GO Classification ---
        go_file <- "__SNK_2__"
        if(file.exists(go_file) && file.info(go_file)$size > 80) {
            go_df <- tryCatch(read.table(go_file, header=TRUE, sep="\\t", stringsAsFactors=FALSE), error=function(e) NULL)
            if(!is.null(go_df) && nrow(go_df) > 0) {
                p <- ggplot(go_df, aes(x=GO)) +
                    geom_bar(fill="#3498db") +
                    labs(title="GO Annotation", x="GO Term Count", y="Gene Count") +
                    theme_bw() + coord_flip()
                ggsave("__SNK_8__", p, width=10, height=8)
                ggsave("__SNK_3__", p, width=10, height=8, dpi=150)
            } else {
                pdf("__SNK_8__"); plot.new(); text(0.5,0.5,"No GO data"); dev.off()
                png("__SNK_3__"); plot.new(); text(0.5,0.5,"No GO data"); dev.off()
            }
        } else {
            pdf("__SNK_8__"); plot.new(); text(0.5,0.5,"No GO annotation"); dev.off()
            png("__SNK_3__"); plot.new(); text(0.5,0.5,"No GO annotation"); dev.off()
        }
        
        # --- KOG Classification ---
        kog_file <- "__SNK_7__"
        if(file.exists(kog_file) && file.info(kog_file)$size > 80) {
            kog_df <- tryCatch(read.table(kog_file, header=TRUE, sep="\\t", stringsAsFactors=FALSE), error=function(e) NULL)
            if(!is.null(kog_df) && nrow(kog_df) > 0) {
                # Truncate long descriptions
                kog_df$label <- paste0(kog_df$category, ": ", substr(kog_df$description, 1, 40))
                p <- ggplot(kog_df, aes(x=reorder(label, count), y=count)) +
                    geom_bar(stat="identity", fill="#e67e22") +
                    coord_flip() +
                    labs(x="", y="Gene Count", title="KOG/COG Classification") +
                    theme_bw() + theme(axis.text.y=element_text(size=7))
                ggsave("__SNK_6__", p, width=12, height=max(6, nrow(kog_df)*0.25))
                ggsave("__SNK_1__", p, width=12, height=max(6, nrow(kog_df)*0.25), dpi=150)
            } else {
                pdf("__SNK_6__"); plot.new(); text(0.5,0.5,"No KOG data"); dev.off()
                png("__SNK_1__"); plot.new(); text(0.5,0.5,"No KOG data"); dev.off()
            }
        } else {
            pdf("__SNK_6__"); plot.new(); text(0.5,0.5,"No KOG annotation"); dev.off()
            png("__SNK_1__"); plot.new(); text(0.5,0.5,"No KOG annotation"); dev.off()
        }
        
        # --- Pfam Top 20 ---
        pfam_file <- "__SNK_4__"
        if(file.exists(pfam_file) && file.info(pfam_file)$size > 80) {
            pfam_df <- tryCatch(read.table(pfam_file, header=TRUE, sep="\\t", stringsAsFactors=FALSE), error=function(e) NULL)
            if(!is.null(pfam_df) && nrow(pfam_df) > 0) {
                pfam_top <- head(pfam_df, 20)
                p <- ggplot(pfam_top, aes(x=reorder(domain, count), y=count)) +
                    geom_bar(stat="identity", fill="#9b59b6") +
                    coord_flip() +
                    labs(x="", y="Gene Count", title="Top 20 Pfam Domains") +
                    theme_bw()
                ggsave("__SNK_5__", p, width=10, height=6)
                ggsave("__SNK_0__", p, width=10, height=6, dpi=150)
            } else {
                pdf("__SNK_5__"); plot.new(); text(0.5,0.5,"No Pfam data"); dev.off()
                png("__SNK_0__"); plot.new(); text(0.5,0.5,"No Pfam data"); dev.off()
            }
        } else {
            pdf("__SNK_5__"); plot.new(); text(0.5,0.5,"No Pfam annotation"); dev.off()
            png("__SNK_0__"); plot.new(); text(0.5,0.5,"No Pfam annotation"); dev.off()
        }
        
        cat("Annotation plots generated\\n")
        
        """
        r_code = r_code.replace("__SNK_0__", str(output.pfam_plot_png))
        r_code = r_code.replace("__SNK_1__", str(output.kog_plot_png))
        r_code = r_code.replace("__SNK_2__", str(input.go_gene2term))
        r_code = r_code.replace("__SNK_3__", str(output.go_plot_png))
        r_code = r_code.replace("__SNK_4__", str(input.pfam_stats))
        r_code = r_code.replace("__SNK_5__", str(output.pfam_plot))
        r_code = r_code.replace("__SNK_6__", str(output.kog_plot))
        r_code = r_code.replace("__SNK_7__", str(input.kog_stats))
        r_code = r_code.replace("__SNK_8__", str(output.go_plot))
        r_code = r_code.replace("__SNK_9__", str(RESULT_DIR))
        result = subprocess.run(["Rscript", "-"], input=r_code, text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- 结果汇总到输出目录 ----
rule collect_annotation_results:
    """将注释结果复制到结果目录的结构化子目录"""
    input:
        eggnog      = rules.eggnog_mapper.output.annotations,
        go_t2g      = rules.extract_go_annotation.output.go_term2gene,
        go_g2t      = rules.extract_go_annotation.output.go_gene2term,
        kegg_g2k    = rules.extract_kegg_annotation.output.kegg_gene2ko,
        kegg_g2p    = rules.extract_kegg_annotation.output.kegg_gene2pathway,
        kog_g2c     = rules.extract_kog_annotation.output.kog_gene2cog,
        kog_stats   = rules.extract_kog_annotation.output.kog_cog_stats,
        pfam_g2p    = rules.extract_pfam_annotation.output.pfam_gene2pfam,
        pfam_stats  = rules.extract_pfam_annotation.output.pfam_stats,
    output:
        # Structured annotation outputs
        anno_all        = f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/new_transcripts.all_annotation.xls",
        anno_stat       = f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/new_transcripts.stat.xls",
        anno_go_xls     = f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/go/new_transcripts.go.xls",
        anno_go_stat    = f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/go/new_transcripts.go_stat.xls",
        anno_go_wego    = f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/go/new_transcripts.weGO.xls",
        anno_kegg_xls   = f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/kegg/new_transcripts.kegg.xls",
        anno_kegg_sum   = f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/kegg/kegg_summary.xls",
        anno_kegg_path  = f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/kegg/new_transcripts.pathway.xls",
        anno_kegg_id2ko = f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/kegg/new_transcripts.id2K_ko.xls",
        anno_kog_xls    = f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/kog/new_transcripts.kog.xls",
        anno_pfam_xls   = f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/pfam/new_transcripts.pfam.xls",
    conda:
        f"{PIPELINE_DIR}/envs/qc.yaml"
    run:
        import os, subprocess, sys
        os.makedirs(f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/go", exist_ok=True)
        os.makedirs(f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/kegg", exist_ok=True)
        os.makedirs(f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/kog", exist_ok=True)
        os.makedirs(f"{RESULT_DIR}/04.transcripts/new_transcripts/anno/pfam", exist_ok=True)
        # [TODO] 原 shell 命令: cp {input.eggnog}      {output.anno_all}
        # [TODO] 原 shell 命令: cp {input.go_t2g}      {output.anno_go_wego}
        # [TODO] 原 shell 命令: cp {input.go_g2t}      {output.anno_go_xls}
        # [TODO] 原 shell 命令: python3 -c "
        # [TODO] 原 shell 命令: import pandas as pd
        # [TODO] 原 shell 命令: try:
        # [TODO] 原 shell 命令: df = pd.read_csv('{input.go_g2t}', sep='\\t', comment='#')
        # [TODO] 原 shell 命令: annotated = len(df) if 'gene' in df.columns else 0
        # [TODO] 原 shell 命令: except:
        # [TODO] 原 shell 命令: annotated = 0
        # [TODO] 原 shell 命令: with open('{output.anno_go_stat}', 'w') as f:
        # [TODO] 原 shell 命令: f.write(f'category\\tcount\\n')
        # [TODO] 原 shell 命令: f.write(f'Annotated\\t{{annotated}}\\n')
        # [TODO] 原 shell 命令: "
        # [TODO] 原 shell 命令: cp {input.kegg_g2k}    {output.anno_kegg_id2ko}
        # [TODO] 原 shell 命令: cp {input.kegg_g2p}    {output.anno_kegg_path}
        # [TODO] 原 shell 命令: python3 -c "
        # [TODO] 原 shell 命令: import pandas as pd
        # [TODO] 原 shell 命令: try:
        # [TODO] 原 shell 命令: df = pd.read_csv('{input.kegg_g2p}', sep='\\t', comment='#')
        # [TODO] 原 shell 命令: summary = df['pathway'].value_counts().reset_index()
        # [TODO] 原 shell 命令: summary.columns = ['pathway','count']
        # [TODO] 原 shell 命令: summary.to_csv('{output.anno_kegg_sum}', sep='\\t', index=False)
        # [TODO] 原 shell 命令: except Exception as e:
        # [TODO] 原 shell 命令: with open('{output.anno_kegg_sum}', 'w') as f:
        # [TODO] 原 shell 命令: f.write(f'# Error: {{e}}\\n')
        # [TODO] 原 shell 命令: "
        # [TODO] 原 shell 命令: python3 -c "
        # [TODO] 原 shell 命令: import pandas as pd
        # [TODO] 原 shell 命令: try:
        # [TODO] 原 shell 命令: df = pd.read_csv('{input.kegg_g2k}', sep='\\t', comment='#')
        # [TODO] 原 shell 命令: df.to_csv('{output.anno_kegg_xls}', sep='\\t', index=False)
        # [TODO] 原 shell 命令: except Exception as e:
        # [TODO] 原 shell 命令: with open('{output.anno_kegg_xls}', 'w') as f:
        # [TODO] 原 shell 命令: f.write(f'# Error: {{e}}\\n')
        # [TODO] 原 shell 命令: "
        # [TODO] 原 shell 命令: cp {input.kog_g2c}     {output.anno_kog_xls}
        # [TODO] 原 shell 命令: cp {input.pfam_g2p}    {output.anno_pfam_xls}
        r_code = """
        import os
        stats = []
        for fpath, label in [
            ('__SNK_8__', 'GO'),
            ('__SNK_0__', 'KEGG_KO'),
            ('__SNK_1__', 'KEGG_Pathway'),
            ('__SNK_6__', 'KOG'),
            ('__SNK_2__', 'Pfam'),
        ]:
            try:
                with open(fpath) as f:
                    lines = [l for l in f if not l.startswith('#') and l.strip()]
                stats.append(f'{label}\\t{len(lines)-1 if len(lines)>1 else 0}')
            except:
                stats.append(f'{label}\\t0')
        
        with open('__SNK_9__', 'w') as f:
            f.write('database\\tannotated_count\\n')
            for s in stats:
                f.write(s + '\\n')
        print('Annotation stat:', stats)
        
        """
        r_code = r_code.replace("__SNK_0__", str(output.anno_kegg_id2ko))
        r_code = r_code.replace("__SNK_1__", str(output.anno_kegg_path))
        r_code = r_code.replace("__SNK_2__", str(output.anno_pfam_xls))
        r_code = r_code.replace("__SNK_3__", str(output.anno_kegg_xls))
        r_code = r_code.replace("__SNK_4__", str(output.anno_kegg_sum))
        r_code = r_code.replace("__SNK_5__", str(output.anno_go_stat))
        r_code = r_code.replace("__SNK_6__", str(output.anno_kog_xls))
        r_code = r_code.replace("__SNK_7__", str(output.anno_go_wego))
        r_code = r_code.replace("__SNK_8__", str(output.anno_go_xls))
        r_code = r_code.replace("__SNK_9__", str(output.anno_stat))
        r_code = r_code.replace("__SNK_10__", str(output.anno_all))
        r_code = r_code.replace("__SNK_11__", str(input.pfam_g2p))
        r_code = r_code.replace("__SNK_12__", str(input.kegg_g2p))
        r_code = r_code.replace("__SNK_13__", str(input.kegg_g2k))
        r_code = r_code.replace("__SNK_14__", str(input.kog_g2c))
        r_code = r_code.replace("__SNK_15__", str(input.eggnog))
        r_code = r_code.replace("__SNK_16__", str(input.go_t2g))
        r_code = r_code.replace("__SNK_17__", str(input.go_g2t))
        r_code = r_code.replace("__SNK_18__", str(RESULT_DIR))
        result = subprocess.run(["python3", "-c", r_code], text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")
