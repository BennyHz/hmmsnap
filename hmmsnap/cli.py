# hmmcall/cli.py
import argparse
import os
import sys
import numpy as np
import pyranges as pr
from .io import (
    read_ratio_bw,
    read_treatment_control_bam,
    read_treatment_control_bw,
    normalize_chromosomes
)
from .smoothing import detect_bin, rebin_data_simple, running_median, smooth_signal_conv
from .hmm import hmm_segment, get_state_names
from .segmentation import segments_from_states_all, bed_from_segments, merge_segments
from .utils import STD_CHR_WITH_CHR, GENOME_SIZES

def main():
    ap = argparse.ArgumentParser(
        description="HMMsnap: Hidden Markov Model with signal smoothing for nuclear domain (NAD/LAD) detection"
    )
    ap.add_argument("--out-dir", "-o", required=True,
                    help="Output directory (required)")

    # 输入选项（四选一）
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--ratio-bw", 
                       help="Precomputed log2(Treatment/Control) bigWig file")
    group.add_argument("--treatment-bam", 
                       help="Treatment BAM file")
    group.add_argument("--treatment-bw", 
                       help="Treatment bigWig file")

    ap.add_argument("--control-bw", 
                    help="Control bigWig file (required when using --treatment-bw)")
    ap.add_argument("--control-bam", 
                    help="Control BAM file (required when using --treatment-bam)")

    # 分辨率参数（区分 BAM 和 bigWig）
    ap.add_argument("--bin-size", "-bs", type=int, default=10000,
                    help="Bin size (bp) for coverage calculation from BAM files (default: 10000)")
    ap.add_argument("--rebin-size", "-rs", type=int, default=None,
                    help="Rebin data to this size (bp) before processing; "
                         "applies to all input types including precomputed ratio (default: None)")

    # 基因组参数（二选一）
    genome_group = ap.add_mutually_exclusive_group()
    genome_group.add_argument("--genome", choices=["hg38", "hg19", "mm10", "mm39"],
                             help="Reference genome for automatic effective genome size selection")
    genome_group.add_argument("--effective-genome-size", type=int, default=None,
                             help="Effective genome size for normalization (e.g., 2652463870 for hg38)")

    # HMM/平滑/后处理
    ap.add_argument("--states", "-s", type=int, nargs="+", default=[2], 
                    help="List of state counts to run (2-5), e.g., --states 2 3 (default: [2])")
    ap.add_argument("--min-len", "-ml", type=int, default=10_000,
                    help="Minimum length of called regions in base pairs (default: 10000)")
    ap.add_argument("--gap-merge", "-gm", type=int, default=None,
                    help="Bridge gaps between regions if ≤ this size (bp) (default: None)")
    ap.add_argument("--stay-strength", "-ss", type=int, default=10,
                    help="Higher value -> more persistence in HMM state (default: 10)")
    ap.add_argument("--seed", type=int, default=1,
                    help="Random seed for HMM initialization (default: 1)")
    ap.add_argument("--pseudocount", "-pc", type=float, default=1e-6,
                    help="Pseudocount added to coverage values to avoid log(0) (default: 1e-6)")

    ap.add_argument("--smooth-method", "-sm", choices=["median", "conv"], default="median",
                    help="Smoothing method: median (default) or conv (convolution)")
    ap.add_argument("--smooth-span", "-sp", type=int, default=5,
                    help="Span (number of bins) for smoothing (default: 5)")
    
    ap.add_argument("--sample", "-n", type=str, default="sample",
                    help="Sample name for output files (default: sample)")

    args = ap.parse_args()

    # 验证 states 参数
    for ns in args.states:
        if ns < 2 or ns > 5:
            sys.exit(f"Error: --states must be between 2 and 5, got {ns}")

    # 验证输入参数
    if args.treatment_bam and not args.control_bam:
        sys.exit("When using --treatment-bam, --control-bam is required.")
    if args.treatment_bw and not args.control_bw:
        sys.exit("When using --treatment-bw, --control-bw is required.")
    
    # 确定 effective genome size
    effective_genome_size = None
    if args.effective_genome_size is not None:
        effective_genome_size = args.effective_genome_size
        print(f"[info] Using provided effective genome size: {effective_genome_size}")
    elif args.genome is not None:
        from .utils import GENOME_SIZES
        effective_genome_size = GENOME_SIZES[args.genome]
        print(f"[info] Using genome {args.genome} with effective genome size: {effective_genome_size}")
    else:
        print("[info] No genome or effective genome size provided. Using raw coverage (no normalization).")
    
    os.makedirs(args.out_dir, exist_ok=True)

    # 读取数据
    if args.ratio_bw:
        df = read_ratio_bw(args.ratio_bw)
    elif args.treatment_bam:
        # 从 BAM 文件生成 ratio
        df = read_treatment_control_bam(
            args.treatment_bam, 
            args.control_bam, 
            bin_size=args.bin_size,
            pseudocnt=args.pseudocount,
            effective_genome_size=effective_genome_size
        )
    else:  # args.treatment_bw is not None
        # 从 bigWig 文件生成 ratio
        df = read_treatment_control_bw(
            args.treatment_bw,
            args.control_bw,
            bin_size=args.bin_size,
            pseudocnt=args.pseudocount
        )

    df = normalize_chromosomes(df)
    df = df[df.Chromosome.isin(STD_CHR_WITH_CHR)].reset_index(drop=True)

    if args.rebin_size is not None:
        print(f"[info] rebining data to {args.rebin_size} bp bins")
        df = rebin_data_simple(df, args.rebin_size)
        bin_bp = args.rebin_size
    else:
        bin_bp = detect_bin(df)
    
    k = args.smooth_span
    print(f"[info] bin={bin_bp} bp, smooth span={k} bins (~{k*bin_bp/1000:.1f} kb)")

    for ns in args.states:
        print(f"[run] HMM with {ns} state(s)")
        
        state_segments = {}

        for chr_ in [c for c in STD_CHR_WITH_CHR if c in df.Chromosome.unique()]:
            cdf = df[df.Chromosome==chr_].reset_index(drop=True)
            if cdf.empty: continue
            y0 = cdf["Ratio"].to_numpy(dtype="float64")
            y0 = y0 - np.median(y0)

            if args.smooth_method == "conv":
                ys = smooth_signal_conv(y0, k)
                print(f"[info] {chr_}: using convolution smoothing with span={k} bins")
            else:
                ys = running_median(y0, k)
                print(f"[info] {chr_}: using running median smoothing with k={k} bins")

            states, means = hmm_segment(
                ys, n_states=ns, stay_strength=args.stay_strength,
                seed=args.seed
            )
            print(f"[debug] {chr_}: {ns}-state means = {means}")
            
            # 获取状态名称
            state_names = get_state_names(ns)
            sorted_states = np.argsort(means)
            state_mapping = {sorted_states[i]: state_names[i] for i in range(ns)}
            print(f"[debug] {chr_} state mapping: {state_mapping}")

            # 提取所有状态的 segments
            for original_state in range(ns):
                segs = segments_from_states_all(cdf, states, original_state, args.min_len)
                state_chr = bed_from_segments(cdf, segs, score_from="Ratio")
                if len(state_chr) > 0:
                    mapped_state_name = state_mapping[original_state]
                    if mapped_state_name not in state_segments:
                        state_segments[mapped_state_name] = []
                    state_segments[mapped_state_name].append(state_chr)

        # 计算 slack（仅来自 --gap-merge）
        slack = args.gap_merge if args.gap_merge is not None else 0

        # 构建 tag
        tag = f"{ns}-state"
        if args.rebin_size is not None:
            tag += f"_{args.rebin_size//1000}kb"
        else:
            tag += f"_{bin_bp//1000}kb"
        smooth_method_name = "median" if args.smooth_method == "median" else "conv"
        tag += f"_{smooth_method_name}_span{k}"
        if args.gap_merge is not None:
            tag += f"_gap{int(round(args.gap_merge/1000))}kb"

        # 输出每个状态
        for state_name, segments_list in state_segments.items():
            if segments_list:
                state_regions = merge_segments(segments_list, slack)
                if len(state_regions) > 0:
                    bed_file = os.path.join(args.out_dir, f"{args.sample}-{state_name}_HMM_{tag}.bed")
                    output_regions = normalize_chromosomes(state_regions[["Chromosome","Start","End"]])
                    pr.PyRanges(output_regions).to_bed(bed_file)
                    print(f"[done] {ns}-state {state_name}: {len(state_regions)} -> {bed_file}")

if __name__ == "__main__":
    main()