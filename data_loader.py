import pandas as pd

def load_data(filepath):
    """
    從 Excel 讀取 SAXS 數據。
    將 Q 轉換為 nm^-1。

    Parameters
    ----------
    filepath : str
        input data excel 檔案路徑

    Returns
    -------
    q : numpy.ndarray
        Q 向量 (nm^-1)
    intensities : pandas.DataFrame
        強度數據的 DataFrame (每欄代表一個時間點)
    """
    try:
        df = pd.read_excel(filepath)
        q_raw = df.iloc[:, 0].values
        
        # 將 Q 從 A^-1 轉換為 nm^-1
        q_nm = q_raw * 10.0
        
        # 排除第一欄 Q
        intensities = df.iloc[:, 1:]
        
        return q_nm, intensities
    except Exception as e:
        raise RuntimeError(f"讀取數據失敗: {e}")
