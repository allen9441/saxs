import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
from data_loader import load_data
import os
import gc
from analysis import SchulzSphereAnalyzer, init_worker, worker_analysis_opt

st.set_page_config(page_title="SAXS Analysis", layout="wide")

st.title("SAXS Data Analysis")

# sidebar
st.sidebar.header("設定")
uploaded_file = st.sidebar.file_uploader("上傳 Excel 檔案", type=["xlsx"])
default_file = "all_data.xlsx"

filepath = uploaded_file if uploaded_file else default_file

@st.cache_data
def get_data(path):
    return load_data(path)

try:
    with st.spinner("載入數據中..."):
        q_nm, intensities = get_data(filepath)
    st.success(f"已載入 {intensities.shape[1]} 幀數據，包含 {len(q_nm)} 個 Q 點。")
except Exception as e:
    st.error(f"數據載入失敗: {e}")
    st.stop()

# parameters
st.sidebar.subheader("Guinier Parameters")
q_min = st.sidebar.number_input("Q Min (nm^-1)", value=0.1075, format="%.4f")
q_max = st.sidebar.number_input("Q Max (nm^-1)", value=0.3103, format="%.4f")

st.sidebar.subheader("運算設定")
cpu_count = os.cpu_count() or 4
default_workers = min(4, cpu_count)
max_workers = st.sidebar.number_input(
    "並行處理核心數 (Max Workers)", 
    min_value=1, 
    max_value=cpu_count, 
    value=default_workers,
    help="減少核心數可降低記憶體使用量，增加核心數可加快分析速度(但可能導致記憶體不足)。"
)

st.sidebar.subheader("後處理")
auto_fix = st.sidebar.checkbox("自動修復離群點", value=True)
outlier_threshold = st.sidebar.slider("離群點判定閾值 (%)", 10, 100, 20, 5) / 100.0
window_size = st.sidebar.number_input("平滑窗口大小", min_value=3, value=15, step=2, help="增加窗口大小可過濾連續的異常點")
min_rg_limit = st.sidebar.number_input("最小合理 Rg (nm)", value=3.0, help="低於此數值將被視為異常")

def remove_outliers(df, column='Rg_Model', window=15, threshold=0.2, min_val=3.0):
    """
    使用移動中位數檢測並修復離群點及填補缺值
    """
    df_clean = df.copy()
    series = df_clean[column]
    
    # 計算移動中位數
    rolling_median = series.rolling(window=window, center=True).median()
    
    # 填充邊緣 NaN
    rolling_median = rolling_median.bfill().ffill()
    
    # 檢測差異
    diff = np.abs(series - rolling_median)
    rel_diff = diff / rolling_median
    
    # 標記離群點 (相對差異超過閾值)
    # 注意：如果 series 是 NaN，比較結果為 False，所以需額外處理 NaN
    outliers = rel_diff > threshold
    
    # 物理範圍過濾
    physical_outliers = (series < min_val) | (series > 100.0)
    
    # 包含 NaN (缺值)
    nan_mask = series.isna()
    
    mask = outliers | physical_outliers | nan_mask
    
    if np.sum(mask) == 0:
        return df_clean, 0
    
    # 替換為插值
    df_clean.loc[mask, column] = np.nan
    df_clean[column] = df_clean[column].interpolate(method='linear', limit_direction='both')
    
    return df_clean, np.sum(mask)

# analysis button
if st.button("開始分析"):
    st.write("開始分析... 請稍候 (可能需要幾分鐘)。")
    
    # init analyzers (for retry logic only)
    schulz = SchulzSphereAnalyzer()
    
    results = []
    
    # progress Bar
    progress_bar = st.progress(0)
    total_frames = intensities.shape[1]
    
    last_good_params = None
    last_good_rg = None

    # data for parallel execution
    intensity_list = [intensities.iloc[:, i].values for i in range(total_frames)]

    with ProcessPoolExecutor(max_workers=max_workers, initializer=init_worker, initargs=(q_nm, (q_min, q_max), (0.98, 1.023))) as executor:
        results_iter = executor.map(worker_analysis_opt, intensity_list)
        
        # iterate through results as they come in (in order)
        for i, (rg_g, r2_g, bg, res_m) in enumerate(results_iter):
            # update progress
            if i % 10 == 0:
                progress_bar.progress(i / total_frames)
            
            # reconstruct local vars for logic
            intensity = intensity_list[i]
            intensity_corr = intensity - bg
            
            # 檢測單一離群值
            is_outlier = False
            if res_m['success'] and last_good_rg is not None:
                 change = abs(res_m['Rg'] - last_good_rg) / last_good_rg
                 if change > 0.3:
                     is_outlier = True
            
            # 如果失敗或疑似離群，且有歷史紀錄，嘗試使用上一幀的結果作為初始猜測重試
            if (not res_m['success'] or is_outlier) and last_good_params is not None:
                 res_retry = schulz.fit_frame(q_nm, intensity_corr, initial_guess=last_good_params)
                 if res_retry['success']:
                     # 如果重試成功，採納重試結果 (避免平滑噪訊導致的單點跳動)
                     res_m = res_retry
    
            if res_m['success']:
                rg_m = res_m['Rg']
                rav_m = res_m['R_avg']
                p_m = res_m['polydispersity']
                
                # 更新歷史紀錄
                last_good_rg = rg_m
                # 確保 params 長度符合 (fit_frame 可能回傳 bg=0)
                last_good_params = [res_m['R_avg'], res_m['z'], res_m['scale'], res_m['background']]
            else:
                rg_m = np.nan
                rav_m = np.nan
                p_m = np.nan
                
            results.append({
                "Time (s)": i,
                "Rg_Guinier": rg_g,
                "Guinier_R2": r2_g,
                "Rg_Model": rg_m,
                "Rav_Model": rav_m,
                "Polydispersity": p_m
            })
        
    progress_bar.progress(1.0)
    
    # convert to DataFrame
    results_df = pd.DataFrame(results)
    
    # auto-fix outliers
    if auto_fix:
        try:
            # Fix Model Rg
            results_df, n_fixed_m = remove_outliers(
                results_df, 
                column='Rg_Model', 
                window=window_size, 
                threshold=outlier_threshold, 
                min_val=min_rg_limit
            )
            
            # Fix Guinier Rg (also filling NaNs)
            results_df, n_fixed_g = remove_outliers(
                results_df, 
                column='Rg_Guinier', 
                window=window_size, 
                threshold=outlier_threshold, 
                min_val=min_rg_limit
            )
            
            if n_fixed_m > 0 or n_fixed_g > 0:
                st.info(f"已自動修復: {n_fixed_m} 個 Model Rg, {n_fixed_g} 個 Guinier Rg 異常點/缺值。")
                
        except Exception as e:
            st.warning(f"修復離群點時發生錯誤: {e}")
            
    st.session_state['results'] = results_df
    st.success("分析完成！")

# display results if available
if 'results' in st.session_state:
    df_res = st.session_state['results']
    
    # Rearrange Logic
    # sort by size (Model Rg)
    df_sorted = df_res.sort_values(by="Rg_Model").reset_index(drop=True)
    df_sorted['Sorted_Index'] = df_sorted.index
    
    # visualization
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("尺寸隨時間變化 (Size vs Time)")
        fig1, ax1 = plt.subplots()
        ax1.plot(df_res["Time (s)"], df_res["Rg_Model"], label="Model Rg", alpha=0.7)
        ax1.plot(df_res["Time (s)"], df_res["Rg_Guinier"], label="Guinier Rg", alpha=0.5)
        ax1.set_xlabel("Time (s)")
        ax1.set_ylabel("Size (nm)")
        ax1.legend()
        st.pyplot(fig1)
        plt.close(fig1)
        
    with col2:
        st.subheader("重新排列後的尺寸分佈")
        st.write("依據尺寸大小排序後的數據")
        fig2, ax2 = plt.subplots()
        # plotting the sorted values vs index (0-1200)
        ax2.plot(df_sorted["Sorted_Index"], df_sorted["Rg_Model"], label="Sorted Model Rg", color='green')
        ax2.set_xlabel("Sorted Index (0-1200)")
        ax2.set_ylabel("Size (nm)")
        ax2.legend()
        st.pyplot(fig2)
        plt.close(fig2)

    # garbage collection
    gc.collect()
        
    # data table
    st.subheader("結果表格")
    st.dataframe(df_res)
    
    # download
    csv = df_res.to_csv(index=False).encode('utf-8')
    st.download_button(
        "下載 CSV",
        csv,
        "saxs_results.csv",
        "text/csv",
        key='download-csv'
    )
    
    # excel download
    import io
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
        df_res.to_excel(writer, index=False, sheet_name='Results')
        df_sorted.to_excel(writer, index=False, sheet_name='Sorted_Results')
        
    st.download_button(
        label="下載 Excel (.xlsx)",
        data=buffer,
        file_name="saxs_results.xlsx",
        mime="application/vnd.ms-excel"
    )
