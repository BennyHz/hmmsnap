# hmmcall/segmentation.py
import numpy as np
import pandas as pd
import pyranges as pr

def segments_from_states_all(chr_df, states, target_state, min_len):
    starts, ends = chr_df.Start.to_numpy(), chr_df.End.to_numpy()
    chrom = chr_df.Chromosome.iloc[0]
    pos = (states == target_state)
    segs = []
    i, n = 0, len(states)
    while i < n:
        if not pos[i]: i += 1; continue
        j = i+1
        while j < n and pos[j] and starts[j]==ends[j-1]: j += 1
        s, e = starts[i], ends[j-1]
        if e - s >= min_len: segs.append((chrom, int(s), int(e), i, j))
        i = j
    return segs

def bed_from_segments(df, segs, score_from="Ratio"):
    if not segs: return pd.DataFrame(columns=["Chromosome","Start","End","Score"])
    vals = df[score_from].to_numpy()
    out = [(c,s,e, float(np.mean(vals[i:j]))) for (c,s,e,i,j) in segs]
    return pd.DataFrame(out, columns=["Chromosome","Start","End","Score"])

def merge_segments(segments_list, slack=0):
    if not segments_list:
        return pd.DataFrame(columns=["Chromosome", "Start", "End", "Score"])
    all_segments = pd.concat(segments_list, ignore_index=True)
    if len(all_segments) == 0:
        return pd.DataFrame(columns=["Chromosome", "Start", "End", "Score"])
    merged_pr = pr.PyRanges(all_segments[["Chromosome", "Start", "End"]]).merge(slack=slack)
    merged_df = merged_pr.df[["Chromosome", "Start", "End"]]
    if len(merged_df) > 0:
        merged_df["seg_id"] = np.arange(len(merged_df))
        merged_pr_with_id = pr.PyRanges(merged_df)
        ov = pr.PyRanges(all_segments[["Chromosome","Start","End","Score"]]).join(merged_pr_with_id).df
        if len(ov) > 0:
            seg_id_col = next((c for c in ov.columns if c.startswith("seg_id")), "seg_id")
            merged_with_score = (
                ov.groupby(seg_id_col)["Score"]
                  .mean()
                  .reset_index()
                  .merge(merged_df, on=seg_id_col, how="left")
                  [["Chromosome", "Start", "End", "Score"]]
            )
            return merged_with_score
    merged_df["Score"] = 0.0
    return merged_df[["Chromosome", "Start", "End", "Score"]]