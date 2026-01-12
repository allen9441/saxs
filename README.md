# SAXS Data Analysis

分析小角 X 光散射 (SAXS) 的隨時間變化數據。依據文獻方法，實現了兩種分析模式來計算粒徑大小。

## 功能特點

1.  **數據讀取與前處理**：
    - 自動讀取 Excel (`.xlsx`) 格式的原始數據。
    - 提取 Q 向量與 1200 筆隨時間變化的強度 (Intensity) 數據。
    - 自動轉換 Q 單位 (由 $\text{\AA}^{-1}$ 轉為 $\text{nm}^{-1}$)。

2.  **分析方法**：
    - **方法一：與模型無關的分析 (Model-Independent / Guinier Analysis)**
        - 依據文獻，利用高 Q 區間 ($0.98 - 1.02 \text{nm}^{-1}$) 計算並扣除 Cluster/背景訊號。
        - 在低 Q 區間 ($0.1075 - 0.3103 \text{nm}^{-1}$) 進行 Guinier 線性擬合 ($ln(I)$ vs $q^2$)。
        - 計算迴轉半徑 ($R_g$)。
    - **方法二：與模型相關的分析 (Model-Dependent / Schulz Spheres)**
        - 使用多分散球體模型擬合數據。
        - 利用非線性最小平方法擬合每一幀數據。
        - 計算平均半徑 ($R_{av}$)、多分散度 ($p$) 及對應的迴轉半徑 ($R_g$)。

3.  **數據視覺化與重排**：
    - **趨勢圖**：顯示尺寸 ($R_g$) 隨時間 (0-1200s) 的變化。
    - **重排分佈圖**：將計算出的尺寸大小進行排序 (Sorting)，呈現尺寸分佈的趨勢 (橫軸為 0-1200 的排序索引)。

4.  **結果輸出**：
    - 提供 CSV 與 Excel 格式的分析結果下載。

## 安裝說明 (建議在虛擬環境下使用)

建立虛擬環境：

```bash
python -m venv .venv
```

進入虛擬環境：

Windows
```cmd
.venv\Scripts\activate.bat
```
Linux
```bash
source .venv/bin/activate
```

請確保使用 Python 3.8 或以上版本，並安裝所需套件：

```bash
pip install -r requirements.txt
```

## 使用方法

1.  將原始數據檔案 `lq_data.xlsx` 置於專案資料夾中（或在程式介面中上傳）。
2.  執行 Streamlit 應用程式：

```bash
streamlit run app.py
```

3.  開啟瀏覽器（預設為 `http://localhost:8501`）。
4.  點擊左側欄位的參數設定（如有需要），然後點擊主畫面的 **"開始分析"** 按鈕開始。
5.  等待進度條完成後，即可查看圖表並下載結果。

## 檔案結構

- `app.py`: 主程式，包含 UI 介面與流程控制，以及多進程分配。
- `analysis.py`: 包含 Guinier 分析與 Schulz 球體模型擬合的核心演算法。
- `data_loader.py`: 負責讀取 Excel 檔案與數據格式化。
