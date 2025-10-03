# hmmcall/hmm.py
import numpy as np
from hmmlearn.hmm import GaussianHMM

def hmm_segment(y, n_states=2, stay_strength=10, seed=1):
    y = y.reshape(-1,1)
    hmm = GaussianHMM(
        n_components=n_states,
        covariance_type="diag",
        n_iter=200,
        tol=1e-3,
        random_state=seed,
        init_params="mc"
    )
    
    # === 自适应初始均值计算（方案1：基于分位数）===
    y_flat = y.flatten()
    
    if n_states == 2:
        # 2-state: 使用 25% 和 75% 分位数
        q25 = np.percentile(y_flat, 25)
        q75 = np.percentile(y_flat, 75)
        means_init = np.array([[q25], [q75]])
        
    elif n_states == 3:
        # 3-state: 使用 20%, 50%, 80% 分位数
        q20 = np.percentile(y_flat, 20)
        q50 = np.percentile(y_flat, 50)
        q80 = np.percentile(y_flat, 80)
        means_init = np.array([[q20], [q50], [q80]])
        
    elif n_states == 4:
        # 4-state: 使用 15%, 35%, 65%, 85% 分位数
        percentiles = [15, 35, 65, 85]
        means_init = np.percentile(y_flat, percentiles).reshape(-1, 1)
        
    elif n_states == 5:
        # 5-state: 使用 10%, 30%, 50%, 70%, 90% 分位数
        percentiles = [10, 30, 50, 70, 90]
        means_init = np.percentile(y_flat, percentiles).reshape(-1, 1)
        
    else:
        # 通用方案：均匀分布的分位数（适用于6+ states）
        step = 100.0 / (n_states + 1)
        percentiles = [step * (i + 1) for i in range(n_states)]
        means_init = np.percentile(y_flat, percentiles).reshape(-1, 1)
    
    # 添加最小间隔保证状态分离（避免数值问题）
    if n_states >= 2:
        min_gap = 0.05  # 最小状态间隔，可根据需要调整
        for i in range(1, n_states):
            if means_init[i, 0] - means_init[i-1, 0] < min_gap:
                means_init[i, 0] = means_init[i-1, 0] + min_gap
    
    hmm.means_init = means_init
    print(f"[debug] Adaptive means_init: {means_init.flatten()}")
    # === 自适应初始化结束 ===
    
    A = np.full((n_states,n_states), 1.0)
    np.fill_diagonal(A, stay_strength)
    A /= A.sum(axis=1, keepdims=True)
    hmm.transmat_ = A
    hmm.startprob_ = np.full(n_states, 1.0/n_states)
    hmm.fit(y)
    states = hmm.predict(y)
    means = np.array([y[states==s].mean() if np.any(states==s) else -1e9 for s in range(n_states)]).flatten()
    return states, means

def get_state_names(n_states):
    """返回状态名称列表，用于命名输出文件"""
    if n_states == 2:
        return ["Background", "Peak"]
    elif n_states == 3:
        return ["Background", "Intermediate", "Peak"]
    elif n_states == 4:
        return ["Background", "Low", "High", "Peak"]
    elif n_states == 5:
        return ["Background", "Low", "Medium", "High", "Peak"]
    else:
        return ["Background"] + [f"State{i}" for i in range(1, n_states-1)] + ["Peak"]