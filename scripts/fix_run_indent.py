#!/usr/bin/env python3
"""修复 run: 块缩进 — 从 conda 内容级提升到 rule 指令级."""
import re, os

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def fix_file(filepath):
    with open(filepath, encoding='utf-8') as f:
        text = f.read()

    original = text

    # 模式: conda:\n<8sp>env_path\n<12sp>run:\n<12sp>code...
    # 需要把 run: 提升到 conda: 同级 (4sp), 代码降到 8sp
    # Step 1: 把 12-space run: → 4-space run:
    text = re.sub(r'\n            run:', r'\n    run:', text)

    # Step 2: 把紧接在 run: 后面的 12-space 代码行 → 8-space
    # 匹配 run:\n 后面所有 12-space 开头的行 (直到遇到非12-space行或空行)
    def fix_run_content(m):
        """修复 run: 块内的代码缩进 (12sp→8sp, 16sp→12sp)"""
        block = m.group(0)
        lines = block.split('\n')
        fixed_lines = []
        for line in lines:
            if line.startswith('                '):  # 16 spaces → 12
                fixed_lines.append(line[4:])
            elif line.startswith('            '):     # 12 spaces → 8
                fixed_lines.append(line[4:])
            else:
                fixed_lines.append(line)
        return '\n'.join(fixed_lines)

    # 匹配从 run: 开始到空行或 rule/conda/input/output/params 结束的块
    text = re.sub(
        r'(    run:\n(?:(?:            .*\n)|(?:                .*\n))*)',
        fix_run_content,
        text
    )

    if text != original:
        with open(filepath, 'w', encoding='utf-8', newline='\n') as f:
            f.write(text)
        return True
    return False


if __name__ == '__main__':
    import sys
    files = sorted(
        os.path.join(SRC, 'rules', f) for f in os.listdir(os.path.join(SRC, 'rules'))
        if f.endswith('.smk')
    )

    total = 0
    for fp in files:
        if os.path.exists(fp):
            if fix_file(fp):
                print(f"  {os.path.basename(fp)}: 已修复")
                total += 1
            else:
                print(f"  {os.path.basename(fp)}: 无需修复")

    print(f"\n修复: {total} 个文件")
