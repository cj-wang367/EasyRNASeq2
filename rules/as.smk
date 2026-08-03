# ============================================================
# 可变剪切模块: rMATS (SE/MXE/A3SS/A5SS/RI)
# 结果结构:
#   07.as/{comp}/
#     ├── AS_summary.xls          # 各事件类型统计
#     ├── rmats.log               # rMATS 运行日志
#     ├── SE.MATS.JC.txt          # 各类型差异可变剪切结果
#     ├── MXE.MATS.JC.txt
#     ├── A3SS.MATS.JC.txt
#     ├── A5SS.MATS.JC.txt
#     └── RI.MATS.JC.txt
# ============================================================

# ---- rMATS 差异可变剪切 ----
rule rmats_diff:
    """rMATS 差异可变剪切分析"""
    input:
        gtf     = config["paths"]["ref_gtf"],
        bams    = expand(f"{RESULT_DIR}/03.align/{{sample}}/{{sample}}.sorted.bam", sample=SAMPLES),
    output:
        se_jc    = f"{RESULT_DIR}/07.as/{{comp}}/SE.MATS.JC.txt",
        mxe_jc   = f"{RESULT_DIR}/07.as/{{comp}}/MXE.MATS.JC.txt",
        a3ss_jc  = f"{RESULT_DIR}/07.as/{{comp}}/A3SS.MATS.JC.txt",
        a5ss_jc  = f"{RESULT_DIR}/07.as/{{comp}}/A5SS.MATS.JC.txt",
        ri_jc    = f"{RESULT_DIR}/07.as/{{comp}}/RI.MATS.JC.txt",
    params:
        read_len    = config["alternative_splicing"]["read_length"],
        sample_info = config["paths"]["sample_info"],
        strand_lib  = {
            "unstranded": "fr-unstranded",
            "firststrand": "fr-firststrand",
            "secondstrand": "fr-secondstrand",
        }[config["alignment"]["strand"]],
        outdir      = f"{RESULT_DIR}/07.as/{{comp}}",
        tmpdir      = f"{TMP_DIR}/07.as/{{comp}}",
    threads: 4
    conda:
        f"{PIPELINE_DIR}/envs/as.yaml"
    shell:
        """
        mkdir -p {params.outdir} {params.tmpdir}

        # 生成 b1.txt (对照组 BAM 列表) 和 b2.txt (实验组 BAM 列表)
        python3 -c "
import pandas as pd
import os

# 从比较组名解析 ctrl/treat 组名（格式: ctrl_vs_treat）
comp = '{wildcards.comp}'
parts = comp.split('_vs_')
if len(parts) != 2:
    raise ValueError(f'Invalid comparison name: {{comp}} (expected ctrl_vs_treat)')

ctrl_group = parts[0]
treat_group = parts[1]
result_dir = '{RESULT_DIR}'
outdir = '{params.outdir}'
tmpdir = '{params.tmpdir}'

df = pd.read_csv('{params.sample_info}', sep='\\t', comment='#')

ctrl_samples = df[df['group'] == ctrl_group]['sample'].tolist()
treat_samples = df[df['group'] == treat_group]['sample'].tolist()

if len(ctrl_samples) == 0:
    raise ValueError(f'No samples found for control group: {{ctrl_group}}')
if len(treat_samples) == 0:
    raise ValueError(f'No samples found for treatment group: {{treat_group}}')

ctrl_bams = [f'{{result_dir}}/03.align/{{s}}/{{s}}.sorted.bam' for s in ctrl_samples]
treat_bams = [f'{{result_dir}}/03.align/{{s}}/{{s}}.sorted.bam' for s in treat_samples]

for bam in ctrl_bams + treat_bams:
    if not os.path.exists(bam):
        print(f'WARNING: BAM file not found: {{bam}}')

with open(f'{{tmpdir}}/b1.txt', 'w') as f:
    f.write(','.join(ctrl_bams))
with open(f'{{tmpdir}}/b2.txt', 'w') as f:
    f.write(','.join(treat_bams))

print(f'b1 ({{ctrl_group}}): {{len(ctrl_samples)}} samples, {{len(ctrl_bams)}} BAMs')
print(f'b2 ({{treat_group}}): {{len(treat_samples)}} samples, {{len(treat_bams)}} BAMs')
"

        # 检查 b1/b2 是否生成成功
        if [ ! -s {params.tmpdir}/b1.txt ] || [ ! -s {params.tmpdir}/b2.txt ]; then
            echo "ERROR: Failed to generate b1/b2 files for rMATS"
            exit 1
        fi

        # 在临时目录运行 rMATS（会产生大量中间文件）
        rmats.py \
            --b1 {params.tmpdir}/b1.txt \
            --b2 {params.tmpdir}/b2.txt \
            --gtf {input.gtf} \
            --od {params.tmpdir} \
            --tmp {params.tmpdir}/tmp \
            -t paired \
            --readLength {params.read_len} \
            --cstat 0.1 \
            --libType {params.strand_lib} \
            --nthread {threads} \
            --task both 2>&1 | tee {params.tmpdir}/rmats.log

        # 只复制最终结果文件到输出目录，中间文件留在临时目录
        echo "=== Extracting final results ==="
        for event in SE MXE A3SS A5SS RI; do
            for suffix in MATS.JC.txt MATS.JCEC.txt; do
                src="{params.tmpdir}/${{event}}.${{suffix}}"
                dst="{params.outdir}/${{event}}.${{suffix}}"
                if [ -f "$src" ] && [ -s "$src" ]; then
                    cp "$src" "$dst"
                    echo "  OK: ${{event}}.${{suffix}} ($(wc -l < "$dst") lines)"
                else
                    echo "  SKIP: ${{event}}.${{suffix}} (empty or missing)"
                fi
            done
        done

        # 确保 JC 文件至少是一个空占位，避免下游 rmats_stats 报错
        for f in {output.se_jc} {output.mxe_jc} {output.a3ss_jc} {output.a5ss_jc} {output.ri_jc}; do
            if [ ! -f "$f" ]; then
                echo -e "ID\\tGeneID\\tgeneSymbol\\tchr\\tstrand\\texonStart_0base\\texonEnd\\tupstreamES\\tupstreamEE\\tdownstreamES\\tdownstreamEE\\tIJC_SAMPLE_1\\tSJC_SAMPLE_1\\tIncFormLen\\tSkipFormLen\\tPValue\\tFDR\\tIncLevel1\\tIncLevel2\\tIncLevelDifference" > "$f"
                echo "  EMPTY: $(basename $f) (placeholder created)"
            fi
        done

        # 复制日志
        cp {params.tmpdir}/rmats.log {params.outdir}/rmats.log
        """


# ---- rMATS 结果统计 ----
rule rmats_stats:
    """汇总 rMATS 各事件类型统计"""
    input:
        se   = f"{RESULT_DIR}/07.as/{{comp}}/SE.MATS.JC.txt",
        mxe  = f"{RESULT_DIR}/07.as/{{comp}}/MXE.MATS.JC.txt",
        a3ss = f"{RESULT_DIR}/07.as/{{comp}}/A3SS.MATS.JC.txt",
        a5ss = f"{RESULT_DIR}/07.as/{{comp}}/A5SS.MATS.JC.txt",
        ri   = f"{RESULT_DIR}/07.as/{{comp}}/RI.MATS.JC.txt",
    output:
        summary = f"{RESULT_DIR}/07.as/{{comp}}/AS_summary.xls",
    conda:
        f"{PIPELINE_DIR}/envs/as.yaml"
    run:
        import os, subprocess, sys
        r_code = """
        import pandas as pd
        import os
        import sys
        
        event_desc = {
            'SE':   'Skipped Exon',
            'MXE':  'Mutually Exclusive Exons',
            'A3SS': "Alternative 3-prime Splice Site",
            'A5SS': "Alternative 5-prime Splice Site",
            'RI':   'Retained Intron',
        }
        
        event_files = {
            'SE':   '__SNK_4__',
            'MXE':  '__SNK_3__',
            'A3SS': '__SNK_1__',
            'A5SS': '__SNK_2__',
            'RI':   '__SNK_5__',
        }
        
        records = []
        for event_type, path in event_files.items():
            if os.path.exists(path) and os.path.getsize(path) > 0:
                try:
                    df = pd.read_csv(path, sep='\\t')
                    total = len(df)
                    # FDR 列名可能是 FDR 或 fdr
                    fdr_col = 'FDR' if 'FDR' in df.columns else 'fdr'
                    if fdr_col in df.columns:
                        sig = int((df[fdr_col] < 0.05).sum())
                    else:
                        sig = 0
                    records.append({
                        'EventType': event_type,
                        'Description': event_desc[event_type],
                        'TotalEvents': total,
                        'Significant_FDR0.05': sig,
                    })
                except Exception as e:
                    print(f"Parse error for {event_type}: {e}", file=sys.stderr)
                    records.append({
                        'EventType': event_type,
                        'Description': event_desc[event_type],
                        'TotalEvents': -1,
                        'Significant_FDR0.05': -1,
                        'Note': f'Parse error: {str(e)[:80]}'
                    })
            else:
                records.append({
                    'EventType': event_type,
                    'Description': event_desc[event_type],
                    'TotalEvents': 0,
                    'Significant_FDR0.05': 0,
                    'Note': 'File empty or not found'
                })
        
        summary = pd.DataFrame(records)
        summary.to_csv('__SNK_0__', sep='\\t', index=False)
        print(summary.to_string(index=False))
        print()
        print(f'AS summary saved to: __SNK_0__')
        
        """
        r_code = r_code.replace("__SNK_0__", str(output.summary))
        r_code = r_code.replace("__SNK_1__", str(input.a3ss))
        r_code = r_code.replace("__SNK_2__", str(input.a5ss))
        r_code = r_code.replace("__SNK_3__", str(input.mxe))
        r_code = r_code.replace("__SNK_4__", str(input.se))
        r_code = r_code.replace("__SNK_5__", str(input.ri))
        result = subprocess.run(["python3", "-c", r_code], text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Script failed (exit {result.returncode})")
        if result.stdout:
            print(result.stdout, end="")


# ---- AS 事件类型分布图 ----
rule as_event_plot:
    """汇总所有比较组的 AS 事件类型分布图"""
    input:
        summaries = expand(f"{RESULT_DIR}/07.as/{{comp}}/AS_summary.xls", comp=COMPARISONS),
    output:
        pdf = f"{RESULT_DIR}/07.as/AS_event_distribution.pdf",
        png = f"{RESULT_DIR}/07.as/AS_event_distribution.png",
    conda:
        f"{PIPELINE_DIR}/envs/diff.yaml"
    run:
        import subprocess, sys
        r_code = """
        library(ggplot2)
        library(dplyr)
        library(tidyr)
        all_data <- list()
        sum_files <- unlist(strsplit("__SUMMARIES__", " "))
        for(f in sum_files) {
            comp_name <- basename(dirname(f))
            df <- tryCatch(read.table(f, header=TRUE, sep="\\t", check.names=FALSE), error=function(e) NULL)
            if(!is.null(df) && nrow(df) > 0 && all(c("EventType","TotalEvents","Significant_FDR0.05") %in% colnames(df))) {
                df$Comparison <- comp_name
                all_data[[length(all_data)+1]] <- df[,c("EventType","TotalEvents","Significant_FDR0.05","Comparison")]
            }
        }
        if(length(all_data) == 0) {
            pdf("__PDF__"); plot.new(); text(0.5,0.5,"No AS data"); dev.off()
            png("__PNG__"); plot.new(); text(0.5,0.5,"No AS data"); dev.off()
            quit(save="no")
        }
        plot_df <- bind_rows(all_data)
        if(nrow(plot_df) == 0) {
            pdf("__PDF__"); plot.new(); text(0.5,0.5,"Empty AS data"); dev.off()
            png("__PNG__"); plot.new(); text(0.5,0.5,"Empty AS data"); dev.off()
            quit(save="no")
        }
        events_long <- plot_df %>%
            pivot_longer(c(TotalEvents, Significant_FDR0.05),
                         names_to="Category", values_to="Count")
        events_long$Category <- dplyr::recode(events_long$Category,
            TotalEvents="All Events", Significant_FDR0.05="Significant (FDR<0.05)")
        w <- max(8, length(unique(plot_df$Comparison))*2)
        p <- ggplot(events_long, aes(x=EventType, y=Count, fill=Category)) +
            geom_bar(stat="identity", position="dodge") +
            facet_wrap(~Comparison, scales="free_y") +
            scale_fill_manual(values=c("All Events"="#3498db","Significant (FDR<0.05)"="#e74c3c")) +
            labs(x="", y="Event Count", title="Alternative Splicing Event Distribution") +
            theme_bw() + theme(axis.text.x=element_text(angle=30, hjust=1), plot.title=element_text(hjust=0.5))
        ggsave("__PDF__", p, width=w, height=6)
        ggsave("__PNG__", p, width=w, height=6, dpi=150)
        """
        r_code = r_code.replace('__SUMMARIES__', ' '.join(input.summaries))
        r_code = r_code.replace('__PDF__', output.pdf)
        r_code = r_code.replace('__PNG__', output.png)
        result = subprocess.run(["Rscript", "-"], input=r_code, text=True,
                                capture_output=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise Exception(f"Rscript failed (exit {result.returncode})")
