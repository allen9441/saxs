import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import linregress
import saxs_core

class GuinierAnalyzer:
    def __init__(self, q_range=(0.1075, 0.3103), bg_q_range=(0.98, 1.023)):
        self.q_min, self.q_max = q_range
        self.bg_q_min, self.bg_q_max = bg_q_range

    def get_background(self, q, intensity):
        """
        從高 Q 區間的平坦區計算背景值。

        Parameters
        ----------
        q : numpy.ndarray
            散射向量 (nm^-1)
        intensity : numpy.ndarray
            散射強度

        Returns
        -------
        float
            計算出的背景強度值
        """
        mask = (q >= self.bg_q_min) & (q <= self.bg_q_max)
        if not np.any(mask):
            return 0.0
        return np.mean(intensity[mask])

    def analyze(self, q, intensity):
        """
        執行 Guinier 分析。

        Parameters
        ----------
        q : numpy.ndarray
            散射向量 (nm^-1)
        intensity : numpy.ndarray
            散射強度

        Returns
        -------
        Rg : float
            迴轉半徑 (Radius of Gyration, nm)
        I0 : float
            零角度散射強度 (Forward Scattering Intensity)
        quality : float
            擬合品質 (R-squared)
        """
        # 1. 背景扣除
        bg = self.get_background(q, intensity)
        I_corr = intensity - bg
        
        # 2. 選擇 Guinier 區間
        mask = (q >= self.q_min) & (q <= self.q_max)
        q_sel = q[mask]
        I_sel = I_corr[mask]
        
        if len(q_sel) < 3 or np.any(I_sel <= 0):
            return np.nan, np.nan, 0.0
            
        # 3. 線性擬合: ln(I) vs q^2
        x = q_sel**2
        y = np.log(I_sel)
        
        try:
            p = np.polyfit(x, y, 1) # 斜率, 截距
            slope = p[0]
            intercept = p[1]
            
            # slope = -Rg^2 / 3
            # Rg = sqrt(-3 * slope)
            if slope >= 0:
                Rg = 0.0 # 對 Guinier 而言無效
            else:
                Rg = np.sqrt(-3 * slope)
                
            I0 = np.exp(intercept)
            
            # R-squared (判定係數)
            y_pred = np.polyval(p, x)
            ss_res = np.sum((y - y_pred)**2)
            ss_tot = np.sum((y - np.mean(y))**2)
            r2 = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0
            
            return Rg, I0, r2
            
        except Exception:
            return np.nan, np.nan, 0.0

class SchulzSphereAnalyzer:
    def __init__(self):
        pass

    def model_intensity(self, q, r_avg, z, scale, bg):
        """
        完整的模型強度計算。
        
        Parameters
        ----------
        q : numpy.ndarray
            散射向量 (nm^-1)
        r_avg : float
            平均半徑 (nm)
        z : float
            Schulz 分佈的寬度參數 (與多分散度有關，z 越大分佈越窄)
        scale : float
            縮放因子 (與粒子濃度、體積有關)
        bg : float
            背景強度

        Returns
        -------
        numpy.ndarray
            計算出的模型散射強度 I(q)
        """
        return saxs_core.model_intensity(q, r_avg, z, scale, bg)

    def estimate_initial_guess(self, q, intensity):
        """
        使用 Guinier 近似動態估算初始半徑。
        """
        # 簡單過濾
        valid = (intensity > 1e-10) & (~np.isnan(intensity))
        if np.sum(valid) < 5:
            return 15.0 # Fallback

        q_c = q[valid]
        I_c = intensity[valid]
        
        # 取前 15 點 (假設是 Low Q)
        n = min(len(q_c), 15)
        x = q_c[:n]**2
        y = np.log(I_c[:n])
        
        try:
            slope, intercept, r_value, p_value, std_err = linregress(x, y)
            if slope >= 0:
                 return 15.0
            
            Rg = np.sqrt(-3 * slope)
            # Sphere Rg = sqrt(3/5) R. So R = sqrt(5/3) Rg
            R_guess = np.sqrt(5/3) * Rg
            return np.clip(R_guess, 1.0, 1000.0)
        except:
            return 15.0

    def fit_frame(self, q, intensity, initial_guess=None, q_max_fit=0.35, fix_background=True):
        """
        對單一幀數據進行模型擬合。

        Parameters
        ----------
        q : numpy.ndarray
            散射向量 (nm^-1)
        intensity : numpy.ndarray
            散射強度 (已扣除背景建議)
        initial_guess : list, optional
            初始猜測值 [R_avg, z, scale, bg] (預設為 None)
        q_max_fit : float, optional
            擬合使用的最大 Q 值 (nm^-1)。預設 0.35 以接近 Guinier 區間，確保擬合主要粒子訊號。
        fix_background : bool, optional
            是否固定背景為 0 (預設為 True，假設已外部扣除)

        Returns
        -------
        dict
            包含擬合結果的字典:
            - 'R_avg': 平均半徑 (nm)
            - 'z': 分佈參數
            - 'scale': 縮放因子
            - 'background': 背景值
            - 'Rg': 計算出的迴轉半徑 (nm)
            - 'polydispersity': 多分散度 (sigma/R)
            - 'success': 是否成功 (bool)
            - 'error': 錯誤訊息 (若失敗)
        """
        # 過濾 NaN/Inf
        valid_mask = np.isfinite(q) & np.isfinite(intensity)
        q_clean = q[valid_mask]
        i_clean = intensity[valid_mask]

        if len(q_clean) < 5:
            return {'success': False, 'error': 'Not enough valid data points'}

        # 動態計算初始猜測值
        guess_R = self.estimate_initial_guess(q_clean, i_clean)
        
        # 確保 scale 為正
        scale_guess = np.mean(i_clean[:5]) if len(i_clean) > 5 else 1e-5
        if scale_guess <= 0: scale_guess = 1e-5

        default_guess = [guess_R, 50.0, scale_guess, 0.0]
        
        # init p0
        if initial_guess is None:
            p0_full = default_guess
        else:
            p0_full = initial_guess

        # 限制 fitting 範圍
        mask = q_clean <= q_max_fit
        q_fit = q_clean[mask]
        i_fit = i_clean[mask]

        if len(q_fit) < 10:
             return {'success': False, 'error': 'Not enough data points in q range'}

        if fix_background:
            # 只擬合前 3 個參數: R_avg, z, scale
            p0 = p0_full[:3]
            bounds = ([0.1, 0.1, 0], [np.inf, np.inf, np.inf])
            
            def fit_func(q_val, r_av, z_val, sc):
                return self.model_intensity(q_val, r_av, z_val, sc, 0.0)
        else:
            p0 = p0_full
            bounds = ([0.1, 0.1, 0, 0], [np.inf, np.inf, np.inf, np.inf])
            
            def fit_func(q_val, r_av, z_val, sc, b):
                return self.model_intensity(q_val, r_av, z_val, sc, b)

        try:
            popt, pcov = curve_fit(fit_func, q_fit, i_fit, p0=p0, bounds=bounds, maxfev=2000)
            
            # 重試機制
            retry = False
            if initial_guess is not None:
                # check if at the rift
                p_check = initial_guess[:len(popt)]
                if np.allclose(popt, p_check, rtol=1e-4):
                    retry = True
            
            if retry:
                 # 如果卡住，使用動態計算的 default_guess 重試
                 p0_retry = default_guess[:len(popt)]
                 popt, pcov = curve_fit(fit_func, q_fit, i_fit, p0=p0_retry, bounds=bounds, maxfev=2000)

            # 整理結果
            if fix_background:
                r_avg, z, scale = popt
                bg = 0.0
            else:
                r_avg, z, scale, bg = popt
            
            # Rg from eq11
            term = (3 * (z + 8) * (z + 7)) / (5 * (z + 1)**2)
            Rg = r_avg * np.sqrt(term)
            
            # 多分散度 p
            p_index = 1 / np.sqrt(z + 1)
            
            return {
                'R_avg': r_avg,
                'z': z,
                'scale': scale,
                'background': bg,
                'Rg': Rg,
                'polydispersity': p_index,
                'success': True
            }
        except Exception as e:
            # 如果失敗，嘗試使用預設值再試一次（如果尚未嘗試）
            if initial_guess is not None:
                try:
                    p0_retry = default_guess[:len(bounds[0])]
                    popt, pcov = curve_fit(fit_func, q_fit, i_fit, p0=p0_retry, bounds=bounds, maxfev=2000)
                    
                    if fix_background:
                        r_avg, z, scale = popt
                        bg = 0.0
                    else:
                        r_avg, z, scale, bg = popt
                        
                    term = (3 * (z + 8) * (z + 7)) / (5 * (z + 1)**2)
                    Rg = r_avg * np.sqrt(term)
                    p_index = 1 / np.sqrt(z + 1)
                    return {
                        'R_avg': r_avg, 'z': z, 'scale': scale, 'background': bg,
                        'Rg': Rg, 'polydispersity': p_index, 'success': True
                    }
                except:
                    pass
                    
            return {'success': False, 'error': str(e)}

def worker_analysis(q, intensity, guinier_q_range, guinier_bg_q_range):
    """
    單次工作單元，用來多進程加速安排。
    """
    guinier = GuinierAnalyzer(q_range=guinier_q_range, bg_q_range=guinier_bg_q_range)
    schulz = SchulzSphereAnalyzer()
    
    # 1. Guinier
    rg_g, i0_g, r2_g = guinier.analyze(q, intensity)
    
    # 2. Schulz Sphere Background
    bg = guinier.get_background(q, intensity)
    
    # 3. Schulz Sphere Fit
    intensity_corr = intensity - bg
    res_m = schulz.fit_frame(q, intensity_corr)
    
    return rg_g, r2_g, bg, res_m

# reduce memory overhead
_worker_q = None
_worker_q_range = None
_worker_bg_q_range = None

def init_worker(q, q_range, bg_q_range):

    global _worker_q, _worker_q_range, _worker_bg_q_range
    _worker_q = q
    _worker_q_range = q_range
    _worker_bg_q_range = bg_q_range

def worker_analysis_opt(intensity):

    global _worker_q, _worker_q_range, _worker_bg_q_range
    if _worker_q is None:
        raise RuntimeError("Worker not initialized properly")
        
    return worker_analysis(_worker_q, intensity, _worker_q_range, _worker_bg_q_range)
