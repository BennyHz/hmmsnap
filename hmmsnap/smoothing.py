# hmmcall/smoothing.py
import pandas as pd
import numpy as np

def detect_bin(df):
    """从数据中自动检测 bin 大小（取中位数）"""
    return int(pd.Series((df.End - df.Start)).median())

def rebin_data_simple(df, bin_size):
    """将数据重新分 bin 到指定大小"""
    df_copy = df.copy()
    df_copy['NewStart'] = (df_copy['Start'] // bin_size) * bin_size
    df_copy['NewEnd'] = df_copy['NewStart'] + bin_size
    result = df_copy.groupby(['Chromosome', 'NewStart', 'NewEnd'])['Ratio'].mean().reset_index()
    result.columns = ['Chromosome', 'Start', 'End', 'Ratio']
    return result.sort_values(['Chromosome', 'Start']).reset_index(drop=True)

def running_median(signal, window_size):
    """使用滑动中位数平滑信号"""
    s = pd.Series(signal, dtype="float64")
    return s.rolling(window=window_size, center=True, min_periods=1).median().to_numpy()

def smooth_signal_conv(signal, span):
    """使用均值卷积平滑信号"""
    if span <= 1:
        return signal
    kernel = np.ones(span) / span
    pad = span // 2
    padded = np.pad(signal, (pad, span - pad - 1), mode='edge')
    return np.convolve(padded, kernel, mode='valid')