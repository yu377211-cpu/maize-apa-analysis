#!/usr/bin/env python3
import subprocess
import argparse
from multiprocessing import Pool, cpu_count
import os
import re

def calculate_ld_pair(args):
    """计算单个SNP对的LD"""
    snp_pair, vcf_file, plink_path, output_dir, results_file = args
    snp1, snp2 = snp_pair
    temp_prefix = os.path.join(output_dir, f"{snp1}.{snp2}.ld")
    
    try:
        cmd = [
            plink_path,
            '--vcf', vcf_file,
            '--ld', snp1, snp2,
            '--double-id',
            '--out', temp_prefix
        ]

        result = subprocess.run(cmd, check=True, capture_output=True, text=True)

        # 解析R2值
        r2_value = parse_ld_output(temp_prefix)

        # 实时写入结果
        with open(results_file, 'a') as f:
            f.write(f"{snp1}\t{snp2}\t{r2_value}\n")

        # 只删除.nosex文件，保留.log文件
        try:
            nosex_file = f"{temp_prefix}.nosex"
            if os.path.exists(nosex_file):
                os.remove(nosex_file)
        except:
            pass

        return (snp1, snp2, r2_value, True)

    except Exception as e:
        # 错误时也写入NA
        with open(results_file, 'a') as f:
            f.write(f"{snp1}\t{snp2}\tNA\n")
        
        # 清理.nosex文件（如果存在）
        try:
            nosex_file = f"{temp_prefix}.nosex"
            if os.path.exists(nosex_file):
                os.remove(nosex_file)
        except:
            pass
        return (snp1, snp2, 'NA', False)

def parse_ld_output(temp_prefix):
    """从PLINK日志文件中解析R2值，支持科学计数法"""
    log_file = f"{temp_prefix}.log"
    
    if not os.path.exists(log_file):
        return 'NA'
    
    try:
        with open(log_file, 'r') as f:
            log_content = f.read()
            
            # 改进的R2提取模式，支持科学计数法
            # 匹配 "R-sq = 0.207689" 或 "R-sq = 2.85923e-05"
            r2_patterns = [
                r"R-sq\s*=\s*([0-9]+\.[0-9]+(?:[eE][-+]?[0-9]+)?)",  # 科学计数法
                r"R-sq\s*=\s*([0-9]\.[0-9]+(?:[eE][-+]?[0-9]+)?)",   # 0.xxx格式
                r"R-sq\s*=\s*(\d*\.?\d+(?:[eE][-+]?\d+)?)"           # 通用格式
            ]
            
            for pattern in r2_patterns:
                match = re.search(pattern, log_content)
                if match:
                    r2_str = match.group(1)
                    # 将科学计数法转换为浮点数，再转回字符串避免精度问题
                    try:
                        r2_float = float(r2_str)
                        # 如果R2值太小，直接返回科学计数法表示
                        if r2_float < 0.0001:
                            return f"{r2_float:.2e}"
                        else:
                            return str(r2_float)
                    except ValueError:
                        return r2_str
            
            # 如果没找到R-sq，尝试其他可能的格式
            alternative_patterns = [
                r"R-squared\s*=\s*([0-9]+\.[0-9]+(?:[eE][-+]?[0-9]+)?)",
                r"R2\s*=\s*([0-9]+\.[0-9]+(?:[eE][-+]?[0-9]+)?)"
            ]
            
            for pattern in alternative_patterns:
                match = re.search(pattern, log_content, re.IGNORECASE)
                if match:
                    r2_str = match.group(1)
                    try:
                        r2_float = float(r2_str)
                        if r2_float < 0.0001:
                            return f"{r2_float:.2e}"
                        else:
                            return str(r2_float)
                    except ValueError:
                        return r2_str
            
            return 'NA'
                
    except Exception as e:
        print(f"解析日志文件错误: {e}")
        return 'NA'

def read_snp_pairs(pairs_file):
    """读取SNP对文件"""
    pairs = []
    with open(pairs_file, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split()
                if len(parts) >= 2:
                    snp1, snp2 = parts[0], parts[1]
                    pairs.append((snp1, snp2))
    return pairs

def test_r2_parsing():
    """测试R2解析功能"""
    test_cases = [
        "R-sq = 0.207689       D' = 0.947291",
        "R-sq = 2.85923e-05    D' = 0.0728061",
        "R-sq = 1.5e-10        D' = 0.001234",
        "R-sq = 0.000123       D' = 0.045678",
        "R-squared = 0.456789  D-prime = 0.987654"
    ]
    
    for i, test_case in enumerate(test_cases):
        # 创建临时日志内容测试
        temp_log = f"PLINK log content\n{test_case}\nEnd of log"
        print(f"测试用例 {i+1}: {test_case}")
        
        # 模拟解析过程
        r2_pattern = r"R-sq\s*=\s*([0-9]+\.[0-9]+(?:[eE][-+]?[0-9]+)?)"
        match = re.search(r2_pattern, test_case)
        if match:
            r2_value = match.group(1)
            print(f"  提取的R2: {r2_value}")
            try:
                r2_float = float(r2_value)
                if r2_float < 0.0001:
                    formatted = f"{r2_float:.2e}"
                else:
                    formatted = str(r2_float)
                print(f"  格式化后: {formatted}")
            except ValueError as e:
                print(f"  转换错误: {e}")
        print()

def main():
    parser = argparse.ArgumentParser(description='并行计算SNP对的LD值')
    parser.add_argument('--vcf', required=True, help='VCF文件路径')
    parser.add_argument('--pairs', required=True, help='SNP对文件路径 (output_500kb.cis.pair格式)')
    parser.add_argument('--output', required=True, help='结果文件路径')
    parser.add_argument('--plink', required=True, help='PLINK可执行文件路径')
    parser.add_argument('--tmpdir', default='./ld_temp', help='临时文件目录')
    parser.add_argument('--processes', type=int, default=None, help='进程数 (默认: CPU核心数)')
    parser.add_argument('--test', action='store_true', help='测试R2解析功能')

    args = parser.parse_args()

    if args.test:
        print("测试R2解析功能...")
        test_r2_parsing()
        return

    # 创建临时目录
    os.makedirs(args.tmpdir, exist_ok=True)

    # 设置进程数
    if args.processes is None:
        num_processes = min(cpu_count(), 8)
    else:
        num_processes = args.processes

    # 读取SNP对
    pairs = read_snp_pairs(args.pairs)
    print(f"读取到 {len(pairs)} 对SNP")

    # 初始化输出文件
    with open(args.output, 'w') as f:
        f.write("SNP_A\tSNP_B\tR2\n")

    print(f"开始计算 {len(pairs)} 对SNP的LD，使用 {num_processes} 个进程")
    print(f"临时文件目录: {args.tmpdir}")
    print(f"结果文件: {args.output}")

    # 准备任务参数
    tasks = [(pair, args.vcf, args.plink, args.tmpdir, args.output) for pair in pairs]

    # 并行计算
    success_count = 0
    with Pool(processes=num_processes) as pool:
        for i, result in enumerate(pool.imap(calculate_ld_pair, tasks)):
            snp1, snp2, r2_value, success = result
            if success:
                success_count += 1

            # 进度显示
            if (i + 1) % 100 == 0:
                print(f"已完成 {i + 1}/{len(pairs)}，成功: {success_count}")

    print(f"\n计算完成!")
    print(f"总SNP对: {len(pairs)}")
    print(f"成功计算: {success_count}")
    print(f"失败: {len(pairs) - success_count}")
    print(f"结果文件: {args.output}")
    print(f"Log文件保存在: {args.tmpdir}")

    # 显示前几个结果作为样例
    print("\n前5个结果样例:")
    with open(args.output, 'r') as f:
        for i, line in enumerate(f):
            if i < 6:  # 头部+5个结果
                print(line.strip())

if __name__ == "__main__":
    main()
