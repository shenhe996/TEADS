# -*- coding: utf-8 -*-
import sys
import math
from collections import defaultdict

def parse_gff(gff_file):
    """解析GFF文件，提取重复序列信息和染色体长度"""
    chrom_data = defaultdict(list)
    chrom_lengths = {}
    
    with open(gff_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '###' or not line.strip():
                continue
                
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
                
            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            attrs = parts[8]
            
            # 从attributes中提取classification，不区分classification的大小写
            classification = None
            for attr in attrs.split(';'):
                if attr.lower().startswith('classification='):
                    classification = attr.split('=')[1]
                    break
            if not classification:
                continue
           
                
            # 提取一级分类
            te_type = classification.split('/')[0]
            
            # 只保留需要的类型：LTR, LINE, SINE, DNA, unknown
            if te_type in ['LTR', 'LINE', 'SINE', 'DNA', 'unknown']:
                chrom_data[chrom].append((start, end, te_type))
                # 更新染色体长度
                if chrom not in chrom_lengths or end > chrom_lengths[chrom]:
                    chrom_lengths[chrom] = end
    
    return chrom_data, chrom_lengths

def merge_intervals(intervals):
    """合并重叠的区间并计算总长度"""
    if not intervals:
        return 0
        
    # 按起始位置排序
    intervals.sort(key=lambda x: x[0])
    
    merged = []
    current_start, current_end = intervals[0]
    
    for start, end in intervals[1:]:
        if start <= current_end + 1:  # 允许相邻区间合并
            current_end = max(current_end, end)
        else:
            merged.append((current_start, current_end))
            current_start, current_end = start, end
    
    merged.append((current_start, current_end))
    
    # 计算总长度
    total_length = 0
    for s, e in merged:
        total_length += e - s + 1
    
    return total_length

def calculate_overlap(start1, end1, start2, end2):
    """计算两个区间的重叠长度"""
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    return max(0, overlap_end - overlap_start + 1) if overlap_start <= overlap_end else 0

def generate_windows(chrom_len):
    """生成100个3Mb窗口，自适应步长"""
    WINDOW_SIZE = 3000000  # 3Mb窗口
    
    if chrom_len <= WINDOW_SIZE:
        # 对于小于3Mb的染色体，只生成一个窗口
        return [(1, chrom_len)]
    
    # 计算自适应步长
    step = (chrom_len - WINDOW_SIZE) / 99.0
    windows = []
    for i in range(100):
        start_pos = 1 + round(i * step)
        end_pos = start_pos + WINDOW_SIZE - 1
        if end_pos > chrom_len:
            end_pos = chrom_len
        windows.append((int(start_pos), int(end_pos)))
    return windows

def main():
    # 将输出重定向到文件
    sys.stdout = open('###output.txt', 'w', encoding='utf-8')
    # 硬编码输入文件路径
    gff_file = r"E:\基因组装\ltr\###.gff3"
    
    try:
        chrom_data, chrom_lengths = parse_gff(gff_file)
    except FileNotFoundError:
        print(f"错误：找不到文件 {gff_file}")
        print("请检查文件路径是否正确")
        sys.exit(1)
    except Exception as e:
        print(f"处理文件时出错: {e}")
        sys.exit(1)
    
    # 输出表头
    print("Chr\tFrom\tTo\tLTR\tLINE\tSINE\tDNA\tunknown\tTE")
    
    # 1. 全基因组统计
    genome_total = 0
    type_genome = defaultdict(int)
    te_genome_intervals = []  # 用于去重计算总TE
    
    # 计算全基因组总长度
    for chrom, length in chrom_lengths.items():
        genome_total += length
    
    # 收集全基因组重复序列
    for chrom, repeats in chrom_data.items():
        for start, end, te_type in repeats:
            type_genome[te_type] += (end - start + 1)
            te_genome_intervals.append((start, end))


    # 合并全基因组重复区间（去重）
    te_genome_merged_length = merge_intervals(te_genome_intervals)
    
    # 计算全基因组比例
    ltr_ratio = type_genome['LTR'] / genome_total
    line_ratio = type_genome['LINE'] / genome_total
    sine_ratio = type_genome['SINE'] / genome_total
    dna_ratio = type_genome['DNA'] / genome_total
    unk_ratio = type_genome['unknown'] / genome_total
    te_ratio = ltr_ratio + line_ratio + sine_ratio + dna_ratio + unk_ratio
    
    # 输出全基因组结果
    print(f"whole_genome\t1\t{genome_total}\t"
          f"{ltr_ratio:.4f}\t{line_ratio:.4f}\t{sine_ratio:.4f}\t"
          f"{dna_ratio:.4f}\t{unk_ratio:.4f}\t{te_ratio:.4f}")
    
    # 2. 处理每个染色体的滑动窗口
    for chrom, repeats in chrom_data.items():
        chrom_len = chrom_lengths.get(chrom, 0)
        if chrom_len == 0:
            continue  # 跳过长度为0的染色体
            
        windows = generate_windows(chrom_len)
        
        # 处理每个窗口
        for start_win, end_win in windows:
            window_len = end_win - start_win + 1
            
            # 收集当前窗口的重复序列
            type_intervals = defaultdict(list)
            te_window_intervals = []  # 用于去重计算总TE
            
            for r_start, r_end, r_type in repeats:
                overlap = calculate_overlap(start_win, end_win, r_start, r_end)
                if overlap > 0:
                    # 计算实际重叠部分
                    overlap_start = max(r_start, start_win)
                    overlap_end = min(r_end, end_win)
                    
                    # 添加到类型特定列表
                    type_intervals[r_type].append((overlap_start, overlap_end))
                    # 添加到总TE列表（用于去重）
                    te_window_intervals.append((overlap_start, overlap_end))
            
            # 计算每个类型的比例（合并重叠区间）
            type_lengths = {}
            for te_type in ['LTR', 'LINE', 'SINE', 'DNA', 'unknown']:
                if te_type in type_intervals:
                    type_lengths[te_type] = merge_intervals(type_intervals[te_type])
                else:
                    type_lengths[te_type] = 0
            
            # 计算总TE比例（去重）
            te_merged_length = merge_intervals(te_window_intervals)
            
            # 计算比例
            ltr_ratio = type_lengths['LTR'] / window_len if window_len > 0 else 0.0
            line_ratio = type_lengths['LINE'] / window_len if window_len > 0 else 0.0
            sine_ratio = type_lengths['SINE'] / window_len if window_len > 0 else 0.0
            dna_ratio = type_lengths['DNA'] / window_len if window_len > 0 else 0.0
            unk_ratio = type_lengths['unknown'] / window_len if window_len > 0 else 0.0
            te_ratio = te_merged_length / window_len if window_len > 0 else 0.0
            
            # 输出结果（保留4位小数）
            print(f"{chrom}\t{start_win}\t{end_win}\t"
                  f"{ltr_ratio:.4f}\t{line_ratio:.4f}\t{sine_ratio:.4f}\t"
                  f"{dna_ratio:.4f}\t{unk_ratio:.4f}\t{te_ratio:.4f}")
        # 关闭文件
    sys.stdout.close()

if __name__ == "__main__":
    main()
