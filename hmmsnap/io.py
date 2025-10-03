# hmmcall/io.py
import pyBigWig
import pysam
import pandas as pd
import numpy as np
import math
from .utils import STD_CHR_WITH_CHR, STD_CHR_WITHOUT_CHR, STD_CHR_MAPPING

def normalize_chromosomes(df):
    def convert_chrom(chrom):
        chrom_str = str(chrom)
        if chrom_str.startswith('chr'):
            return chrom_str
        else:
            return STD_CHR_MAPPING.get(chrom_str, chrom_str)
    df_copy = df.copy()
    df_copy['Chromosome'] = df_copy['Chromosome'].apply(convert_chrom)
    return df_copy

def read_ratio_bw(ratio_bw, keep=STD_CHR_WITH_CHR):
    bw = pyBigWig.open(ratio_bw)
    file_chroms = list(bw.chroms().keys())
    if any(c.startswith('chr') for c in file_chroms):
        target_chroms = STD_CHR_WITH_CHR
    else:
        target_chroms = STD_CHR_WITHOUT_CHR
    chrs = [c for c in file_chroms if c in target_chroms]
    rows = []
    for chr_ in chrs:
        ivs = bw.intervals(chr_)
        if not ivs: continue
        for (s,e,v) in ivs:
            if v is None or np.isnan(v): continue
            rows.append((chr_, s, e, float(v)))
    bw.close()
    if not rows: raise RuntimeError("No data read from ratio bigWig.")
    df = pd.DataFrame(rows, columns=["Chromosome","Start","End","Ratio"])
    return df.sort_values(["Chromosome","Start"], ignore_index=True)

def calculate_coverage_from_bam(bam_file, bin_size, chroms, effective_genome_size=None):
    """
    从 BAM 文件计算 coverage，使用高效的方法
    """
    import pysam
    
    # 打开 BAM 文件
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    # 获取染色体长度
    chrom_lengths = {}
    for ref, length in zip(bam.references, bam.lengths):
        if ref in chroms:
            chrom_lengths[ref] = length
    
    # 计算总 mapped reads（用于标准化）
    total_mapped = 0
    if effective_genome_size is not None:
        print(f"[info] Counting total mapped reads in {bam_file}...")
        for ref in chroms:
            if ref in chrom_lengths:
                # 只计算标准染色体上的 mapped reads
                total_mapped += bam.count(
                    reference=ref,
                    read_callback=lambda read: not (read.is_unmapped or read.is_duplicate or read.is_secondary)
                )
        scaling_factor = effective_genome_size / total_mapped if total_mapped > 0 else 1.0
        print(f"[info] Total mapped reads: {total_mapped}, scaling_factor: {scaling_factor:.6f}")
    else:
        scaling_factor = 1.0
    
    rows = []
    
    for chrom in chroms:
        if chrom not in chrom_lengths:
            continue
            
        chrom_len = chrom_lengths[chrom]
        # 生成 bins
        starts = list(range(0, chrom_len, bin_size))
        ends = [min(start + bin_size, chrom_len) for start in starts]
        
        print(f"[info] Processing {chrom} ({len(starts)} bins)...")
        
        # 为每个 bin 计算 coverage
        for i, (start, end) in enumerate(zip(starts, ends)):
            # 使用 count() 获取该区域的 read 数量
            read_count = bam.count(
                reference=chrom,
                start=start,
                end=end,
                read_callback=lambda read: not (read.is_unmapped or read.is_duplicate or read.is_secondary)
            )
            
            # 简单的 coverage 估算：reads per bin
            coverage = read_count
            
            # 应用 scaling
            normalized_coverage = coverage * scaling_factor
            rows.append((chrom, start, end, normalized_coverage))
            
            # 进度提示（每 1000 个 bins 显示一次）
            if (i + 1) % 1000 == 0:
                print(f"  ... processed {i + 1}/{len(starts)} bins")
    
    bam.close()
    return rows

def read_treatment_control_bam(treatment_bam, control_bam, bin_size=10000, pseudocnt=1e-6, 
                              effective_genome_size=None):
    """从 Treatment 和 Control BAM 文件计算 log2(Treatment/Control)"""
    
    # 获取染色体列表
    treatment_header = pysam.AlignmentFile(treatment_bam, "rb").header
    control_header = pysam.AlignmentFile(control_bam, "rb").header
    
    treatment_chroms = set(treatment_header.references)
    control_chroms = set(control_header.references)
    
    if any(c.startswith('chr') for c in treatment_chroms):
        target_chroms = [c for c in STD_CHR_WITH_CHR if c in treatment_chroms and c in control_chroms]
    else:
        target_chroms = [c for c in STD_CHR_WITHOUT_CHR if c in treatment_chroms and c in control_chroms]
    
    if not target_chroms:
        raise RuntimeError("No common standard chromosomes found in BAM files.")
    
    print(f"[info] Using chromosomes: {target_chroms}")
    
    # 计算 coverage
    print(f"[info] Calculating Treatment coverage...")
    treatment_rows = calculate_coverage_from_bam(
        treatment_bam, bin_size, target_chroms, effective_genome_size
    )
    
    print(f"[info] Calculating Control coverage...")
    control_rows = calculate_coverage_from_bam(
        control_bam, bin_size, target_chroms, effective_genome_size
    )
    
    # 转为 DataFrame 并合并
    treatment_df = pd.DataFrame(treatment_rows, columns=["Chromosome", "Start", "End", "Coverage"])
    control_df = pd.DataFrame(control_rows, columns=["Chromosome", "Start", "End", "Coverage"])
    
    # 合并两个 DataFrame
    merged_df = pd.merge(treatment_df, control_df, on=["Chromosome", "Start", "End"], suffixes=('_t', '_c'))
    
    if len(merged_df) == 0:
        raise RuntimeError("No overlapping bins between Treatment and Control.")
    
    # 计算 ratio
    merged_df["Ratio"] = np.log2(
        (merged_df["Coverage_t"] + pseudocnt) / (merged_df["Coverage_c"] + pseudocnt)
    )
    
    result_df = merged_df[["Chromosome", "Start", "End", "Ratio"]].sort_values(
        ["Chromosome", "Start"]
    ).reset_index(drop=True)
    
    print(f"[info] Generated {len(result_df)} ratio bins")
    return result_df

def read_treatment_control_bw(treatment_bw, control_bw, bin_size=10000, pseudocnt=1e-6):
    """
    从 Treatment 和 Control bigWig 文件计算 log2(Treatment/Control)
    使用标准基因组位置对齐，解决 bin 不匹配问题
    """
    bt = pyBigWig.open(treatment_bw)
    bc = pyBigWig.open(control_bw)
    
    # 获取染色体信息
    treatment_chroms = set(bt.chroms().keys())
    control_chroms = set(bc.chroms().keys())
    
    if any(c.startswith('chr') for c in treatment_chroms):
        target_chroms = [c for c in STD_CHR_WITH_CHR if c in treatment_chroms and c in control_chroms]
    else:
        target_chroms = [c for c in STD_CHR_WITHOUT_CHR if c in treatment_chroms and c in control_chroms]
    
    if not target_chroms:
        raise RuntimeError("No common standard chromosomes found in bigWig files.")
    
    print(f"[info] Using chromosomes: {target_chroms}")
    
    rows = []
    for chrom in target_chroms:
        chrom_len = bt.chroms(chrom)
        if chrom_len is None or bc.chroms(chrom) is None:
            print(f"[warning] {chrom}: chromosome length not found in one of the bigWig files, skipping.")
            continue
            
        # 使用相同的 bin 策略生成标准位置
        starts = list(range(0, chrom_len, bin_size))
        ends = [min(start + bin_size, chrom_len) for start in starts]
        
        valid_bins = 0
        for start, end in zip(starts, ends):
            try:
                # 从 bigWig 获取该区域的值（自动处理缺失值）
                v1 = bt.stats(chrom, start, end, type="mean")[0]
                v2 = bc.stats(chrom, start, end, type="mean")[0]
                
                # 处理缺失值和 NaN
                if v1 is None or v2 is None or math.isnan(v1) or math.isnan(v2):
                    continue
                    
                ratio = math.log2((v1 + pseudocnt) / (v2 + pseudocnt))
                rows.append((chrom, start, end, ratio))
                valid_bins += 1
                
            except (RuntimeError, IndexError, ValueError) as e:
                # 某些区域可能无法获取值，跳过
                continue
        
        print(f"[info] {chrom}: processed {valid_bins} bins out of {len(starts)}")
    
    bt.close()
    bc.close()
    
    if not rows:
        raise RuntimeError("No valid data extracted from Treatment/Control bigWig files.")
    
    df = pd.DataFrame(rows, columns=["Chromosome", "Start", "End", "Ratio"])
    return df.sort_values(["Chromosome", "Start"], ignore_index=True)