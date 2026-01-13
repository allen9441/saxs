import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
from data_loader import load_data
import os
from analysis import SchulzSphereAnalyzer, init_worker, worker_analysis_opt

st.set_page_config(page_title="SAXS Analysis", layout="wide")

st.title("SAXS Data Analysis")

# sidebar
st.sidebar.header("設定")
uploaded_file = st.sidebar.file_uploader("上傳 Excel 檔案", type=["xlsx"])
default_file = "lq_data.xlsx"

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

    # use ProcessPoolExecutor to bypass GIL for CPU-bound tasks
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
            
            # 檢測離群值
            is_outlier = False
            if res_m['success'] and last_good_rg is not None:
                 # 變化超過 30% 且 Rg 小於 5 nm (通常離群點是向下掉)，視為sus
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
    st.success("分析完成 (Analysis Complete)!")
    
    # convert to DataFrame
    results_df = pd.DataFrame(results)
    st.session_state['results'] = results_df

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
