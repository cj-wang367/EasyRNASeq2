#!/usr/bin/env python3
"""
解析 Snakemake 运行日志，生成简洁的失败摘要。

用法:
  python scripts/parse_run_log.py test_data/test2/EasyRNASeq_run_XXXXXX.log

输出:
  - 管道进度 (完成/失败/总计)
  - 每个失败规则的名称、原因、wildcards
  - 跳过规则原因 (如 "输入文件缺失")
"""
import re, sys, os
from collections import defaultdict

def parse_log(log_path):
    with open(log_path, errors='replace') as f:
        content = f.read()

    lines = content.split('\n')

    # 提取 Job stats
    total_jobs = 0
    stats_section = False
    for line in lines:
        if 'Job stats:' in line:
            stats_section = True
            continue
        if stats_section:
            if line.strip().startswith('total'):
                m = re.search(r'total\s+(\d+)', line)
                if m:
                    total_jobs = int(m.group(1))
                break

    # 提取所有错误
    errors = []
    in_error = False
    current_error = None
    error_lines = []

    for i, line in enumerate(lines):
        # 检测错误开始
        if 'Error in rule ' in line:
            if current_error:
                errors.append(current_error)
            m = re.search(r'Error in rule (\S+):', line)
            rule_name = m.group(1) if m else '?'
            current_error = {
                'rule': rule_name,
                'jobid': '',
                'reason': '',
                'output': '',
            }
            in_error = True
            error_lines = []
            continue

        if in_error and current_error:
            error_lines.append(line)

            m = re.search(r'jobid:\s*(\S+)', line)
            if m:
                current_error['jobid'] = m.group(1)

            m = re.search(r'output:\s*(.+)', line)
            if m:
                current_error['output'] = m.group(1).strip()[:100]

            m = re.search(r'reason:\s*(.+)', line)
            if m:
                current_error['reason'] = m.group(1).strip()[:150]

            # Check for specific error messages in surrounding lines
            if 'ModuleNotFoundError' in line:
                current_error['reason'] = 'ModuleNotFoundError: ' + line.strip()[:100]
            elif 'returned non-zero exit status' in line:
                current_error['reason'] = 'Shell command failed (non-zero exit)'
            elif 'NameError' in line:
                m = re.search(r"NameError:.*", line)
                if m:
                    current_error['reason'] = m.group(0)[:150]

            if 'message:' in line and 'None' not in line:
                m = re.search(r'message:\s*(.+)', line)
                if m and m.group(1).strip() and m.group(1).strip() != 'None':
                    current_error['reason'] = m.group(1).strip()[:150]

    if current_error:
        errors.append(current_error)

    # 提取中间性 error 信息 (深层次原因)
    deep_errors = []
    for i, line in enumerate(lines):
        if 'RuleException:' in line:
            deep_errors.append(line.strip()[:200])

    # 输出报告
    print("=" * 65)
    print("  EasyRNASeq v2 — 运行摘要")
    print(f"  日志: {os.path.basename(log_path)}")
    print(f"  计划任务数: {total_jobs}")
    print(f"  失败任务数: {len(errors)}")
    print("=" * 65)

    if not errors:
        print("\n✓ 所有任务成功完成!")
        return

    print(f"\n--- 失败规则 ({len(errors)}) ---\n")
    for e in errors:
        print(f"  [{e['rule']}]")
        print(f"    jobid: {e['jobid']}")
        if e['reason']:
            print(f"    reason: {e['reason']}")
        if e['output']:
            print(f"    output: {e['output']}")
        print()

    if deep_errors:
        print(f"--- 深层错误 ({len(deep_errors)}) ---\n")
        for de in deep_errors[:15]:
            print(f"  {de}")
        if len(deep_errors) > 15:
            print(f"  ... 还有 {len(deep_errors)-15} 条")
        print()

    # 错误类型汇总
    error_types = defaultdict(list)
    for e in errors:
        reason = e['reason'] or '未知原因'
        # 简化归类
        if 'ModuleNotFoundError' in reason:
            cat = '缺少Python模块'
        elif 'NameError' in reason:
            cat = 'Snakemake变量未定义'
        elif 'non-zero exit' in reason:
            cat = 'Shell命令执行失败'
        elif 'Missing output' in reason:
            cat = '输出文件缺失'
        elif 'Missing input' in reason:
            cat = '输入文件缺失'
        elif 'Insufficient' in reason or 'Too few' in reason:
            cat = '数据不足(防御性跳过)'
        elif 'Code has changed' in reason:
            cat = '代码变更触发重新运行'
        else:
            cat = '其他'
        error_types[cat].append(e['rule'])

    print("--- 错误类型汇总 ---\n")
    for cat, rules in sorted(error_types.items()):
        print(f"  {cat}: {len(rules)} 个")
        for r in rules:
            print(f"    - {r}")

    print(f"\n详细错误请查看原始日志: {log_path}")
    print("使用 --report report.html 可生成 HTML 可视化报告")


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("用法: python parse_run_log.py <log_file>")
        sys.exit(1)
    parse_log(sys.argv[1])
