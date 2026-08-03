#!/usr/bin/env python3
"""
将 .smk 文件中的 shell heredoc (Rscript << 'REOF' / python3 << 'PYEOF') 转换为 run 块。
用法: python scripts/convert_heredoc_to_run.py [--dry-run] [rules/align_qc.smk ...]
"""
import re, sys, os

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def find_heredoc_blocks(content):
    """找到所有 shell: + heredoc 块, 返回 (start, end, info) 列表."""
    blocks = []
    pos = 0

    # 匹配 shell: 指令后跟 triple-quoted 字符串
    # shell_indent = shell: 本身的缩进; indent = """ 的缩进 (深一级)
    shell_pat = re.compile(
        r'(?P<shell_indent>[ \t]*)shell:\s*\n'
        r'(?P<indent>\s+)"""\s*\n)',
        re.MULTILINE
    )

    while pos < len(content):
        m = shell_pat.search(content, pos)
        if not m:
            break

        shell_start = m.start()
        shell_indent = m.group('shell_indent')
        indent = m.group('indent')  # """ 的缩进 = 内容缩进
        body_start = m.end()

        # 找到对应的 """ (同缩进)
        close_marker = f'\n{indent}"""'
        close_pos = content.find(close_marker, body_start)
        if close_pos < 0:
            # 尝试只匹配 """ (某些文件可能缩进不一致)
            close_pos2 = content.find(f'{indent}"""', body_start)
            if close_pos2 >= 0:
                close_pos = close_pos2 - 1  # 回退换行符
            else:
                pos = body_start
                continue

        shell_body = content[body_start:close_pos]
        shell_end = close_pos + len(close_marker)

        # 检查是否包含 heredoc
        heredoc_match = re.search(
            r"(?:Rscript|python3?)\s*<<\s*'(?P<delim>[A-Z]+)'\s*\n",
            shell_body
        )

        if heredoc_match:
            delim = heredoc_match.group('delim')
            lang = 'R' if 'Rscript' in heredoc_match.group() else 'Python'

            # 找到 heredoc body (在 Rscript/python3 << 'XXX' 之后到 XXX 之前)
            hd_start_in_shell = heredoc_match.end()
            # 找 REOF/PYEOF 结束标记 (单独一行, 有相同或更少缩进)
            hd_end_pat = re.compile(
                r'^(\s*)' + re.escape(delim) + r'\s*$',
                re.MULTILINE
            )
            hd_end_match = hd_end_pat.search(shell_body, hd_start_in_shell)

            if hd_end_match:
                hd_body = shell_body[hd_start_in_shell:hd_end_match.start()]
                hd_end_in_shell = hd_end_match.end()

                # heredoc 之前的 bash 命令
                pre_heredoc = shell_body[:heredoc_match.start()].strip()
                # heredoc 之后的内容 (通常为空或只有空白)
                post_heredoc = shell_body[hd_end_in_shell:].strip()

                blocks.append({
                    'shell_start': shell_start,
                    'shell_end': shell_end,
                    'shell_indent': shell_indent,  # shell: 关键字缩进
                    'indent': indent,              # """ 内容缩进 (深一级)
                    'language': lang,
                    'delimiter': delim,
                    'pre_heredoc': pre_heredoc,
                    'hd_body': hd_body,
                    'post_heredoc': post_heredoc,
                })

        pos = shell_end

    return blocks


def extract_snakemake_vars(text):
    """提取 Snakemake 变量引用, 排除双花括号转义."""
    # 先移走所有 {{ }} 转义
    cleaned = re.sub(r'\{\{', '  ', text)
    cleaned = re.sub(r'\}\}', '  ', cleaned)

    vars_found = set()
    # 匹配 {input.xxx}, {output.xxx}, {params.xxx}, {wildcards.xxx},
    # {RESULT_DIR}, {TMP_DIR}, {PIPELINE_DIR}, {config...}, {threads}
    for m in re.finditer(
        r'\{('
        r'input\.[a-zA-Z_][\w.]*'
        r'|output\.[a-zA-Z_][\w.]*'
        r'|params\.[a-zA-Z_][\w.]*'
        r'|wildcards\.[a-zA-Z_][\w.]*'
        r'|RESULT_DIR'
        r'|TMP_DIR'
        r'|PIPELINE_DIR'
        r'|config\[[^\]]+\](?:\[[^\]]+\])*'
        r'|threads'
        r')\}',
        cleaned
    ):
        vars_found.add(m.group(0))

    # 按长度降序排序 (长变量先替换, 避免子串冲突)
    return sorted(vars_found, key=len, reverse=True)


def convert_block(content, block_info):
    """转换单个 heredoc 块为 run 块."""
    indent = block_info['indent']              # """ 内容缩进 (给 Python 代码用)
    shell_indent = block_info['shell_indent']  # shell: 关键字缩进 (给 run: 用)
    hd_body = block_info['hd_body']
    pre_heredoc = block_info['pre_heredoc']

    # 1. 提取 Snakemake 变量
    snk_vars = extract_snakemake_vars(hd_body)
    # 也检查 pre_heredoc 中的变量
    if pre_heredoc:
        snk_vars = extract_snakemake_vars(pre_heredoc + '\n' + hd_body)

    # 2. 为每个 Snakemake 变量创建占位符
    var_to_placeholder = {}
    for i, var in enumerate(snk_vars):
        placeholder = f'__SNK_{i}__'
        var_to_placeholder[var] = placeholder

    # 3. 替换 R/Python 代码中的 Snakemake 变量为占位符
    r_code = hd_body
    for var, placeholder in var_to_placeholder.items():
        r_code = r_code.replace(var, placeholder)

    # 4. 还原双花括号转义 {{ → {, }} → }
    r_code = r_code.replace('{{', '{').replace('}}', '}')

    # 5. 处理缩进: 去掉公共缩进
    lines = r_code.split('\n')
    # 找到非空行的最小缩进
    min_indent = None
    for line in lines:
        if line.strip():
            stripped = len(line) - len(line.lstrip())
            if min_indent is None or stripped < min_indent:
                min_indent = stripped
    if min_indent and min_indent > 0:
        lines = [l[min_indent:] if l.startswith(' ' * min_indent) else l for l in lines]
    r_code_clean = '\n'.join(lines)

    # 6. 生成 run 块
    run_lines = []
    run_lines.append(f'{shell_indent}run:')     # run: 与 shell: 同级缩进

    # 导入
    run_lines.append(f'{indent}import os, subprocess, sys')

    # pre_heredoc 中的 bash 命令转换为 Python
    # 使用 f-string 保留原始 {var} 语法 (f-string 和 Snakemake 共用 {xx} 语法)
    if pre_heredoc:
        for cmd_line in pre_heredoc.split('\n'):
            cmd_line = cmd_line.strip()
            if not cmd_line or cmd_line.startswith('#'):
                continue
            # mkdir -p → os.makedirs, 用 f-string 保留变量
            if cmd_line.startswith('mkdir -p '):
                dir_path = cmd_line[len('mkdir -p '):].strip()
                run_lines.append(f'{indent}os.makedirs(f"{dir_path}", exist_ok=True)')
            else:
                # 其他 bash 命令保留为注释
                run_lines.append(f'{indent}# [TODO] 原 shell 命令: {cmd_line}')

    # R/Python 代码字符串
    run_lines.append(f'{indent}r_code = """')
    for line in r_code_clean.split('\n'):
        run_lines.append(f'{indent}{line}')
    run_lines.append(f'{indent}"""')

    # 占位符替换语句
    for var, placeholder in var_to_placeholder.items():
        # {input.xxx} → input.xxx
        py_expr = var[1:-1]  # 去掉花括号
        run_lines.append(f'{indent}r_code = r_code.replace("{placeholder}", str({py_expr}))')

    # 执行语句
    if block_info['language'] == 'R':
        run_lines.append(f'{indent}result = subprocess.run(["Rscript", "-"], input=r_code, text=True, capture_output=True)')
    else:
        run_lines.append(f'{indent}result = subprocess.run(["python3", "-c", r_code], text=True, capture_output=True)')

    run_lines.append(f'{indent}if result.returncode != 0:')
    run_lines.append(f'{indent}    print(result.stderr, file=sys.stderr)')
    run_lines.append(f'{indent}    raise Exception(f"Script failed (exit {{result.returncode}})")')
    run_lines.append(f'{indent}if result.stdout:')
    run_lines.append(f'{indent}    print(result.stdout, end="")')

    return '\n'.join(run_lines)


def process_file(filepath, dry_run=False):
    """处理单个 .smk 文件."""
    with open(filepath, encoding='utf-8') as f:
        content = f.read()

    blocks = find_heredoc_blocks(content)

    if not blocks:
        print(f"  {os.path.basename(filepath)}: 未找到 heredoc shell 块")
        return 0

    # 从后往前替换 (保持位置偏移正确)
    conversions = 0
    for block_info in reversed(blocks):
        new_block = convert_block(content, block_info)
        content = (content[:block_info['shell_start']] +
                   new_block +
                   content[block_info['shell_end']:])
        conversions += 1

    if dry_run:
        print(f"  {os.path.basename(filepath)}: 将转换 {conversions} 个 heredoc → run 块")
        return 0

    with open(filepath, 'w', encoding='utf-8', newline='\n') as f:
        f.write(content)

    print(f"  {os.path.basename(filepath)}: 已转换 {conversions} 个 heredoc → run 块")
    return conversions


if __name__ == '__main__':
    dry_run = '--dry-run' in sys.argv
    files = [a for a in sys.argv[1:] if a != '--dry-run']

    if not files:
        # 默认处理所有 .smk 文件
        files = sorted(
            os.path.join(SRC, 'rules', f) for f in os.listdir(os.path.join(SRC, 'rules'))
            if f.endswith('.smk')
        )

    print("=" * 60)
    print("  转换 shell heredoc → run 块")
    if dry_run:
        print("  (预览模式, 不实际修改)")
    print("=" * 60)

    total = 0
    for fp in files:
        if os.path.exists(fp):
            total += process_file(fp, dry_run)
        else:
            print(f"  [SKIP] 文件不存在: {fp}")

    print(f"\n总计转换: {total} 个 heredoc 块\n")
    if not dry_run:
        print("下一步:")
        print("  1. 检查改动: git diff rules/")
        print("  2. 验证 DAG: snakemake --dryrun")
        print("  3. 同步到服务器后重跑")
