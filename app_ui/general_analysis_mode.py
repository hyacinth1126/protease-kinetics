import sys
from datetime import datetime
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
import numpy as np
import os
from pathlib import Path
from scipy.optimize import curve_fit

def _debug_log(msg: str) -> None:
    """Cloud 로그용: stderr에 출력 후 flush"""
    try:
        ts = datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")
        sys.stderr.write(f"[{ts}] [general_analysis_mode] {msg}\n")
        sys.stderr.flush()
    except Exception:
        pass

_debug_log("module load: starting")

_debug_log("module load: importing mode_general_analysis.analysis")
from mode_general_analysis.analysis import (
    UnitStandardizer,
    DataNormalizer,
    RegionDivider,
    ModelA_SubstrateDepletion,
    ModelB_EnzymeDeactivation,
    ModelC_MassTransfer,
    ModelD_ConcentrationDependentFmax,
    ModelE_ProductInhibition,
    ModelF_EnzymeSurfaceSequestration
)
_debug_log("module load: importing mode_general_analysis.plot")
from mode_general_analysis.plot import Visualizer
_debug_log("module load: importing mode_prep_raw_data.prep (michaelis_menten_calibration)")
from mode_prep_raw_data.prep import michaelis_menten_calibration
_debug_log("module load: general_analysis_mode imports complete")


def _wide_mm_curves_to_long(df):
    """Convert wide-format MM curves (concentration, time_min, rfu_interpolated repeated) to long format."""
    cols = list(df.columns)
    if len(cols) < 3 or len(cols) % 3 != 0:
        return df
    n_blocks = len(cols) // 3
    long_rows = []
    for i in range(len(df)):
        for b in range(n_blocks):
            conc_val = df.iloc[i, b * 3]
            time_val = df.iloc[i, b * 3 + 1]
            rfu_val = df.iloc[i, b * 3 + 2]
            long_rows.append({
                'Concentration': conc_val,
                'Time_min': time_val,
                'RFU_Interpolated': rfu_val
            })
    out = pd.DataFrame(long_rows)
    try:
        out['Concentration [μM]'] = pd.to_numeric(out['Concentration'].astype(str).str.replace('μM', '').str.replace('μg/mL', '').str.strip(), errors='coerce')
        if out['Concentration [μM]'].notna().any():
            out = out.drop(columns=['Concentration'], errors='ignore')
        else:
            out = out.drop(columns=['Concentration [μM]'], errors='ignore')
    except Exception:
        pass
    return out


def verbose_callback(message: str, level: str = "info"):
    """Callback function for logging from analysis modules"""
    if level == "error":
        st.error(message)
    elif level == "warning":
        st.warning(message)
    elif level == "debug":
        st.code(message)
    else:
        st.info(message)


# Known enzyme molecular weights (kDa) for raw_* filename convention
_KNOWN_ENZYME_MW = {"RgpA": 56.0, "RgpB": 43.3}


def _parse_filename_for_model_simulation(filename: str):
    """
    Parse CSV/XLSX filename to extract Enzyme Name, MW (kDa), Substrate Name.
    - raw_*.ext: * is enzyme name (e.g. raw_RgpA.csv → RgpA); MW from _KNOWN_ENZYME_MW if present.
    - EnzymeName_MW_SubstrateName.ext: (e.g. RgpA_56.60_Dabcyl-HEK-K(FITC).csv)
    Returns (enzyme_name, mw_float or None, substrate_name).
    """
    if not filename or not isinstance(filename, str):
        return None, None, None
    base = os.path.splitext(filename.strip())[0]
    if base.lower().startswith('raw_') and len(base) > 4:
        enzyme_name = base[4:].strip()
        mw = _KNOWN_ENZYME_MW.get(enzyme_name)
        return enzyme_name, mw, None
    parts = base.split('_')
    enzyme_name = None
    mw = None
    substrate_name = None
    if len(parts) >= 3:
        enzyme_name = parts[0].strip() or None
        try:
            mw = float(parts[1].strip())
        except (ValueError, TypeError):
            pass
        substrate_name = '_'.join(parts[2:]).strip() or None
    elif len(parts) == 2:
        enzyme_name = parts[0].strip() or None
        try:
            mw = float(parts[1].strip())
        except (ValueError, TypeError):
            pass
    elif len(parts) == 1 and parts[0].strip():
        enzyme_name = parts[0].strip()
    return enzyme_name, mw, substrate_name


def general_analysis_mode(st):
    """Model Simulation Mode - Standard FRET Analysis"""
    _debug_log("general_analysis_mode(): entered")

    # Default values; will be overridden from session_state if set by filename parsing
    default_mw = 56.6
    default_enzyme_name = "Kgp"
    default_substrate_name = "Dabcyl-HEK-K(FITC)"

    # Sidebar configuration
    enzyme_mw = st.sidebar.number_input(
        "Enzyme Molecular Weight (kDa)",
        min_value=1.0,
        max_value=500.0,
        value=st.session_state.get('parsed_enzyme_mw', default_mw),
        step=0.1,
        key="model_sim_enzyme_mw",
        help="Enter enzyme molecular weight required for concentration conversion. Auto-filled from filename (EnzymeName_MW_SubstrateName)."
    )
    
    enzyme_name = st.sidebar.text_input(
        "Enzyme Name (Optional)",
        value=st.session_state.get('parsed_enzyme_name', default_enzyme_name),
        placeholder="enzyme",
        key="model_sim_enzyme_name",
        help="Enzyme name displayed in graph legend. Auto-filled from filename (EnzymeName_MW_SubstrateName)."
    )
    if enzyme_name.strip() == "":
        enzyme_name = "enzyme"
    
    substrate_name = st.sidebar.text_input(
        "Substrate Name (Optional)",
        value=st.session_state.get('parsed_substrate_name', default_substrate_name),
        placeholder="substrate",
        key="model_sim_substrate_name",
        help="Substrate name displayed in graph legend. Auto-filled from filename (EnzymeName_MW_SubstrateName)."
    )
    if substrate_name.strip() == "":
        substrate_name = "substrate"

    # Separator before data source section
    st.sidebar.markdown("---")
    st.sidebar.subheader("📁 Data Source")
    
    uploaded_file = st.sidebar.file_uploader(
        "Upload CSV/XLSX File (Fitted Curves)",
        type=['csv', 'xlsx'],
        key="model_sim_upload",
        help="Result file from Data Load Mode. Filename format 'EnzymeName_MW_SubstrateName.csv' auto-fills the fields above."
    )

    # When a file is loaded, parse filename and update session_state so inputs show parsed values on next run
    if uploaded_file is not None and getattr(uploaded_file, 'name', None):
        current_name = uploaded_file.name
        last_parsed = st.session_state.get('last_parsed_filename_for_model')
        if last_parsed != current_name:
            ename, mw, sname = _parse_filename_for_model_simulation(current_name)
            if ename is not None:
                st.session_state['parsed_enzyme_name'] = ename
                st.session_state['model_sim_enzyme_name'] = ename
            if mw is not None:
                st.session_state['parsed_enzyme_mw'] = mw
                st.session_state['model_sim_enzyme_mw'] = mw
            if sname is not None:
                st.session_state['parsed_substrate_name'] = sname
                st.session_state['model_sim_substrate_name'] = sname
            st.session_state['last_parsed_filename_for_model'] = current_name
            st.rerun()
    
    # Download sample Fitted Curves (Data Load Mode results)
    try:
        with open("data_interpolation_mode/results/MM_interpolated_curves.csv", "rb") as f:
            sample_bytes = f.read()
        st.sidebar.download_button(
            label="📥 Download Data Load Result CSV",
            data=sample_bytes,
            file_name="MM_interpolated_curves.csv",
            mime="text/csv",
            help="Result CSV file generated from Data Load Mode"
        )
    except Exception:
        pass
    
    # Step 1: Load Fitted Curves data (원본 데이터 플롯용)
    df_fitted = None
    rfu_col = None
    
    # 0순위: Session State 확인 (Data Load 모드에서 방금 실행된 경우)
    if 'interpolation_results' in st.session_state and st.session_state.get('mm_data_ready', False):
        try:
            results = st.session_state['interpolation_results']
            df_fitted = results['interp_df'].copy()
            rfu_col = 'RFU_Interpolated' if 'RFU_Interpolated' in df_fitted.columns else 'RFU_Calculated'
            st.session_state['data_source_type'] = 'Data Load result (in memory)'
            # Check experiment_type to display basis
            experiment_type = results.get('experiment_type', 'Substrate Concentration Variation (Standard Michaelis-Menten)')
            if experiment_type == "Substrate Concentration Variation (Standard Michaelis-Menten)":
                result_type = "Substrate-based"
            else:
                result_type = "Enzyme-based"
            st.success(f"Results Applied ({result_type})")
        except Exception as e:
            # 메모리 로드 실패 시 파일 로드 시도
            pass

    if uploaded_file is not None:
        # 업로드된 파일 처리 (업로드 시 항상 파일 기준으로 덮어씀)
        import tempfile
        file_extension = uploaded_file.name.split('.')[-1].lower()
        
        with tempfile.NamedTemporaryFile(delete=False, suffix=f'.{file_extension}', mode='wb') as tmp_file:
            tmp_file.write(uploaded_file.getbuffer())
            tmp_path = tmp_file.name
        
        st.session_state['data_source_type'] = f"Uploaded file: {uploaded_file.name}"
        try:
            if file_extension == 'xlsx':
                # XLSX 파일: "Time–FLU Interpolated curves" 또는 구 시트명 시트 읽기
                xl = pd.ExcelFile(tmp_path, engine='openpyxl')
                curve_sheet = ('Time–FLU Interpolated curves' if 'Time–FLU Interpolated curves' in xl.sheet_names else
                               ('Michaelis-Menten Curves' if 'Michaelis-Menten Curves' in xl.sheet_names else
                                ('Michaelis-Menten curves' if 'Michaelis-Menten curves' in xl.sheet_names else None)))
                df_fitted = pd.read_excel(tmp_path, sheet_name=curve_sheet, engine='openpyxl') if curve_sheet else pd.read_excel(tmp_path, sheet_name=1, engine='openpyxl')
                c0, c1, c2 = (df_fitted.columns[0], df_fitted.columns[1], df_fitted.columns[2]) if len(df_fitted.columns) >= 3 else (None, None, None)
                if c0 is not None and 'concentration' in str(c0).lower() and 'time_min' in str(c1).lower() and 'rfu_interpolated' in str(c2).lower() and len(df_fitted.columns) % 3 == 0:
                    df_fitted = _wide_mm_curves_to_long(df_fitted)
                rfu_col = 'RFU_Interpolated' if 'RFU_Interpolated' in df_fitted.columns else ('rfu_interpolated' if 'rfu_interpolated' in df_fitted.columns else 'RFU_Calculated')
            else:
                # CSV 파일
                df_fitted = pd.read_csv(tmp_path)
                rfu_col = 'RFU_Interpolated' if 'RFU_Interpolated' in df_fitted.columns else 'RFU_Calculated'
        finally:
            os.unlink(tmp_path)
    else:
        # 업로드 없음: Data Load 결과(df_fitted)가 이미 있으면 유지, 없을 때만 파일에서 로드
        import os
        from pathlib import Path
        
        if df_fitted is None:
            # 1순위: XLSX 파일 (Michaelis-Menten_calibration_results.xlsx)
            xlsx_paths = [
                'Michaelis-Menten_calibration_results.xlsx',
                str(Path(__file__).parent.parent / 'Michaelis-Menten_calibration_results.xlsx'),
            ]
            
            for path in xlsx_paths:
                try:
                    if os.path.exists(path):
                        xl = pd.ExcelFile(path, engine='openpyxl')
                        curve_sheet = ('Time–FLU Interpolated curves' if 'Time–FLU Interpolated curves' in xl.sheet_names else
                                       ('Michaelis-Menten Curves' if 'Michaelis-Menten Curves' in xl.sheet_names else
                                        ('Michaelis-Menten curves' if 'Michaelis-Menten curves' in xl.sheet_names else None)))
                        df_fitted = pd.read_excel(path, sheet_name=curve_sheet or 'Time–FLU Interpolated curves', engine='openpyxl') if curve_sheet else pd.read_excel(path, sheet_name=1, engine='openpyxl')
                        c0, c1, c2 = (df_fitted.columns[0], df_fitted.columns[1], df_fitted.columns[2]) if len(df_fitted.columns) >= 3 else (None, None, None)
                        if c0 is not None and 'concentration' in str(c0).lower() and 'time_min' in str(c1).lower() and 'rfu_interpolated' in str(c2).lower() and len(df_fitted.columns) % 3 == 0:
                            df_fitted = _wide_mm_curves_to_long(df_fitted)
                        rfu_col = 'RFU_Interpolated' if 'RFU_Interpolated' in df_fitted.columns else ('rfu_interpolated' if 'rfu_interpolated' in df_fitted.columns else 'RFU_Calculated')
                        break
                except Exception:
                    continue
        
            # 2순위: CSV 파일
            if df_fitted is None:
                csv_paths = [
                    'data_interpolation_mode/results/MM_interpolated_curves.csv',
                    str(Path(__file__).parent.parent / 'data_interpolation_mode' / 'results' / 'MM_interpolated_curves.csv'),
                ]
                for path in csv_paths:
                    try:
                        if os.path.exists(path):
                            df_fitted = pd.read_csv(path)
                            rfu_col = 'RFU_Interpolated' if 'RFU_Interpolated' in df_fitted.columns else 'RFU_Calculated'
                            break
                    except Exception:
                        continue
        
        if df_fitted is None:
            st.error("Data Load Mode result file not found. Please run 'Data Load Mode' to download results or upload a CSV/XLSX file.")
            st.stop()
        
        # rfu_col이 아직 설정되지 않았으면 설정
        if rfu_col is None and df_fitted is not None:
            if 'RFU_Interpolated' in df_fitted.columns:
                rfu_col = 'RFU_Interpolated'
            elif 'RFU_Calculated' in df_fitted.columns:
                rfu_col = 'RFU_Calculated'
            else:
                rfu_col = 'RFU_Interpolated'  # 기본값
    
    # 엑셀 파일의 보간된 곡선 데이터 사용
    # Detect RFU column name
    rfu_col = None
    if 'RFU_Calculated' in df_fitted.columns:
        rfu_col = 'RFU_Calculated'
    elif 'RFU_Interpolated' in df_fitted.columns:
        rfu_col = 'RFU_Interpolated'
    else:
        st.error("RFU data column not found. (RFU_Calculated or RFU_Interpolated)")
        st.stop()
    
    # 엑셀 파일의 데이터 변환
    df_raw_converted = []
    unique_times = sorted(df_fitted['Time_min'].unique())
    
    # 농도 컬럼 이름 감지 (우선순위: ug/mL -> uM -> Concentration)
    conc_col_name = 'Concentration'
    if 'Concentration [ug/mL]' in df_fitted.columns:
        conc_col_name = 'Concentration [ug/mL]'
    elif 'Concentration [μM]' in df_fitted.columns:
        conc_col_name = 'Concentration [μM]'
    elif 'Concentration' in df_fitted.columns:
        conc_col_name = 'Concentration'
    
    for time in unique_times:
        time_data = df_fitted[df_fitted['Time_min'] == time]
        
        # Create row for each concentration
        for _, row in time_data.iterrows():
            conc_val = row.get(conc_col_name, 0)
            rfu = row[rfu_col]
            
            df_raw_converted.append({
                'time_min': time,
                'enzyme_ugml': conc_val,
                'conc_col_name': conc_col_name, # 원래 컬럼 이름 저장
                'FL_intensity': rfu,
                'SD': 0  # 보간된 곡선 데이터는 SD 없음
            })
    
    df_raw = pd.DataFrame(df_raw_converted)
    
    # 시간 범위 저장
    original_time_max = df_raw['time_min'].max()
    
    # 데이터 정보
    unique_times = sorted(df_raw['time_min'].unique())
    unique_concs = sorted(df_raw['enzyme_ugml'].unique())
    
    # Store data source type for later use
    st.session_state['data_source_type'] = 'Fitted Curves (from Data Load mode)'
    st.session_state['original_time_max'] = original_time_max
    # 원본 fitted 데이터 저장 (Data Load 모드와 동일한 그래프를 그리기 위해)
    # df_fitted는 보간된 곡선 데이터이므로 원본 데이터 플롯에 사용
    if df_fitted is not None:
        st.session_state['df_fitted_original'] = df_fitted
        # rfu_col도 저장 (원본 데이터 플롯용)
        if rfu_col is not None:
            st.session_state['rfu_col'] = rfu_col
        else:
            # rfu_col이 없으면 기본값 사용
            st.session_state['rfu_col'] = 'RFU_Interpolated'
    
    # 우선순위: 1) normalization_results의 exponential 식 값, 2) interpolated 값의 최소/최대값
    fitted_params = None
    xlsx_path_for_mm_results = None
    
    # 0순위: Session State의 normalization_results 확인 (Data Load 모드의 exponential 식에서 나온 F0, Fmax)
    if 'interpolation_results' in st.session_state and st.session_state.get('mm_data_ready', False):
        try:
            results = st.session_state['interpolation_results']
            
            # Normalization results에서 가져오기 (가장 정확함 - exponential 식에서 나온 값)
            if 'normalization_results' in results:
                fitted_params = {}
                norm_results = results['normalization_results']
                for conc_name, data in norm_results.items():
                    conc_val = data['concentration']
                    fitted_params[float(conc_val)] = {
                        'F0': float(data['F0']),
                        'Fmax': float(data['Fmax'])
                    }
            # 없으면 dataframe에서 가져오기
            elif 'mm_results_df' in results:
                df_mm = results['mm_results_df']
                fitted_params = {}
                # 농도 컬럼 찾기
                conc_col = None
                for col in df_mm.columns:
                    if 'Concentration' in col:
                        conc_col = col
                        break
                
                if conc_col:
                    for _, row in df_mm.iterrows():
                        if pd.notna(row['F0']) and pd.notna(row['Fmax']):
                            fitted_params[float(row[conc_col])] = {
                                'F0': float(row['F0']),
                                'Fmax': float(row['Fmax'])
                            }
            
            if fitted_params and len(fitted_params) > 0:
                st.session_state['fitted_params'] = fitted_params
        except Exception as e:
            pass

    # 1순위: Interpolated 값에서 F0, Fmax 계산 (농도별 최소값/최대값) - exponential 식이 없을 때만 사용
    if fitted_params is None or len(fitted_params) == 0:
        fitted_params_from_interp = {}
        if df_fitted is not None and rfu_col in df_fitted.columns:
            for conc_val in unique_concs:
                # 같은 농도의 모든 interpolated 값 가져오기
                conc_data = df_fitted[df_fitted[conc_col_name] == conc_val]
                if len(conc_data) > 0:
                    rfu_values = conc_data[rfu_col].values
                    F0_interp = float(np.min(rfu_values))  # 최소값 = F0
                    Fmax_interp = float(np.max(rfu_values))  # 최대값 = Fmax
                    fitted_params_from_interp[float(conc_val)] = {
                        'F0': F0_interp,
                        'Fmax': Fmax_interp
                    }
        
        if fitted_params_from_interp and len(fitted_params_from_interp) > 0:
            fitted_params = fitted_params_from_interp
            st.session_state['fitted_params'] = fitted_params
            st.sidebar.success(f"✅ F0, Fmax computed from interpolated values ({len(fitted_params)} concentrations)")

    # 2순위: MM Results 시트에서 읽기 (exponential 식이나 interpolated 값이 없을 때만)
    if fitted_params is None or len(fitted_params) == 0:
        # 업로드된 파일 또는 자동 로드된 파일 경로 확인
        if uploaded_file is not None:
            import tempfile
            file_extension = uploaded_file.name.split('.')[-1].lower()
            if file_extension == 'xlsx':
                with tempfile.NamedTemporaryFile(delete=False, suffix='.xlsx', mode='wb') as tmp_file:
                    tmp_file.write(uploaded_file.getbuffer())
                    xlsx_path_for_mm_results = tmp_file.name
        else:
            # 자동 로드된 파일 경로 사용
            xlsx_paths = [
                'Michaelis-Menten_calibration_results.xlsx',
                str(Path(__file__).parent.parent / 'Michaelis-Menten_calibration_results.xlsx'),
            ]
            for path in xlsx_paths:
                if os.path.exists(path):
                    xlsx_path_for_mm_results = path
                    break
    
    # MM Results 시트 읽기
    if xlsx_path_for_mm_results is not None:
        try:
            # 1. Model simulation input (또는 구파일 호환 MM Results) — F0, Fmax, v0
            xl_mm = pd.ExcelFile(xlsx_path_for_mm_results, engine='openpyxl')
            mm_sheet = 'Model simulation input' if 'Model simulation input' in xl_mm.sheet_names else ('Analysis mode input' if 'Analysis mode input' in xl_mm.sheet_names else ('MM Results' if 'MM Results' in xl_mm.sheet_names else None))
            df_mm_results = pd.read_excel(xlsx_path_for_mm_results, sheet_name=mm_sheet, engine='openpyxl') if mm_sheet else None
            
            if df_mm_results is not None and 'F0' in df_mm_results.columns and 'Fmax' in df_mm_results.columns:
                fitted_params = {}
                conc_col_name = 'Concentration [ug/mL]' if 'Concentration [ug/mL]' in df_mm_results.columns else 'Concentration'
                
                # v0 data extraction lists
                v0_concs = []
                v0_vals = []
                
                for _, row in df_mm_results.iterrows():
                    conc_value = row[conc_col_name]
                    if pd.notna(conc_value):
                        try:
                            conc_float = float(conc_value)
                            # F0, Fmax
                            if pd.notna(row['F0']) and pd.notna(row['Fmax']):
                                fitted_params[conc_float] = {
                                    'F0': float(row['F0']),
                                    'Fmax': float(row['Fmax'])
                                }
                            
                            # v0
                            v0_val = row.get('v0', row.get('v0 (RFU/min)', row.get('v₀ (RFU/min)')))
                            if pd.notna(v0_val):
                                v0_concs.append(conc_float)
                                v0_vals.append(float(v0_val))
                        except (ValueError, TypeError):
                            continue
                
                if len(fitted_params) > 0:
                    st.sidebar.success(f"✅ F0, Fmax parameters loaded ({len(fitted_params)} concentrations, MM Results sheet)")
                    # Interpolated 값이 없을 때만 MM Results 사용
                    if 'fitted_params' not in st.session_state or len(st.session_state.get('fitted_params', {})) == 0:
                        st.session_state['fitted_params'] = fitted_params
                else:
                    if 'fitted_params' not in st.session_state or len(st.session_state.get('fitted_params', {})) == 0:
                        fitted_params = None
                        st.session_state['fitted_params'] = None
                
                # Store v0 data from file
                if v0_concs and v0_vals:
                    st.session_state['v0_data_from_file'] = {
                        'concentrations': v0_concs,
                        'v0_values': v0_vals
                    }
            else:
                fitted_params = None
                st.session_state['fitted_params'] = None

            # 2. Michaelis-Menten Fit Results (Vmax, Km)
            try:
                xl = pd.ExcelFile(xlsx_path_for_mm_results)
                fit_sheet = 'Michaelis-Menten Fit Results' if 'Michaelis-Menten Fit Results' in xl.sheet_names else ('Fit results' if 'Fit results' in xl.sheet_names else ('MM Fit Results' if 'MM Fit Results' in xl.sheet_names else None))
                if fit_sheet:
                    df_fit = pd.read_excel(xlsx_path_for_mm_results, sheet_name=fit_sheet, engine='openpyxl')
                    mm_fit_from_file = {}
                    
                    # Determine columns
                    p_col = '파라미터' if '파라미터' in df_fit.columns else 'Parameter'
                    v_col = '값' if '값' in df_fit.columns else 'Value'
                    
                    if p_col in df_fit.columns and v_col in df_fit.columns:
                        params = dict(zip(df_fit[p_col], df_fit[v_col]))
                        
                        def get_param(keys):
                            for k in keys:
                                found = next((p for p in params if k.lower() in str(p).lower()), None)
                                if found: return params[found]
                            return None
                        
                        vmax = get_param(['Vmax'])
                        km = get_param(['Km'])
                        r2 = get_param(['R²', 'R2', 'R_squared'])
                        slope = get_param(['Slope'])
                        intercept = get_param(['Intercept'])
                        
                        def to_float(x):
                            try: return float(x)
                            except: return None
                            
                        if vmax is not None and km is not None:
                            mm_fit_from_file = {
                                'Vmax': to_float(vmax),
                                'Km': to_float(km),
                                'R_squared': to_float(r2),
                                'fit_success': True,
                                'experiment_type': "Substrate Concentration Variation (Standard MM)",
                                'equation': f"v₀ = {to_float(vmax):.2f}[S] / ({to_float(km):.2f} + [S])"
                            }
                        elif slope is not None:
                            mm_fit_from_file = {
                                'slope': to_float(slope),
                                'intercept': to_float(intercept) if intercept else 0,
                                'R_squared': to_float(r2),
                                'fit_success': True,
                                'experiment_type': "Enzyme Concentration Variation",
                                'equation': f"v₀ = {to_float(slope):.4f}[E] + {to_float(intercept) if intercept else 0:.4f}"
                            }
                    
                    if mm_fit_from_file:
                         st.session_state['mm_fit_from_file'] = mm_fit_from_file
            except Exception:
                pass

        except Exception:
            fitted_params = None
            st.session_state['fitted_params'] = None
        finally:
            # 임시 파일 삭제
            if uploaded_file is not None and xlsx_path_for_mm_results and os.path.exists(xlsx_path_for_mm_results):
                try:
                    os.unlink(xlsx_path_for_mm_results)
                except:
                    pass
    else:
        fitted_params = None
        st.session_state['fitted_params'] = None
    
    # Step 2: Standardize units
    standardizer = UnitStandardizer(enzyme_mw=enzyme_mw)
    df_standardized = standardizer.standardize(df_raw)
    
    # Store time unit for later use
    time_unit = 'min' if 'time_min' in df_raw.columns else 's'
    st.session_state['time_unit'] = time_unit
    
    # Step 3-4: Normalization and region division
    normalizer = DataNormalizer()
    region_divider = RegionDivider()
    
    # Step 3-1: Initial temporary normalization (model-free threshold or fitted params)
    df_current = normalizer.normalize_temporary(df_standardized, fitted_params=fitted_params)
    
    # Step 4: Divide regions
    df_current = region_divider.divide_regions(df_current)
    
    # Step 3-2: Final normalization (using region information or fitted params)
    df_current = normalizer.normalize_final(df_current, fitted_params=fitted_params)
    
    df = df_current
    
    # Display data (aligned with Data Load Mode: 📋 Data Preview + 4 metrics + expander)
    st.subheader("📋 Data Preview")
    
    # Data Load에서 가져온 결과가 있으면 동일한 값 사용 (Data Preview 일치)
    interp_results = st.session_state.get('interpolation_results') if st.session_state.get('mm_data_ready') else None

    # Detect original column names for display
    time_unit = st.session_state.get('time_unit', 'min')
    original_time_max = st.session_state.get('original_time_max', df['time_s'].max())
    if interp_results is not None and interp_results.get('reaction_time_max_min') is not None:
        reaction_max = interp_results['reaction_time_max_min']
        time_display = f"{reaction_max:.0f} min"
        time_label = "Time (min)"
    elif time_unit == 'min':
        time_display = f"{original_time_max:.0f} min"
        time_label = "Time (min)"
    else:
        time_display = f"{original_time_max:.0f} s" if original_time_max < 100 else f"{original_time_max/60:.1f} min"
        time_label = "Time (s)"
    conc_col = 'enzyme_ugml'

    experiment_type = None
    if interp_results is not None and 'mm_fit_results' in interp_results:
        experiment_type = interp_results['mm_fit_results'].get('experiment_type')
    if experiment_type is None and 'mm_fit_from_file' in st.session_state:
        experiment_type = st.session_state['mm_fit_from_file'].get('experiment_type')

    original_conc_col = df['conc_col_name'].iloc[0] if 'conc_col_name' in df.columns else 'Concentration [ug/mL]'
    _substrate_std = ("Substrate 농도 변화 (표준 MM)", "Substrate Concentration Variation (Standard MM)")
    if experiment_type in _substrate_std:
        conc_unit = "μM"
    elif 'uM' in original_conc_col or 'μM' in original_conc_col:
        conc_unit = "μM"
    elif 'nM' in original_conc_col:
        conc_unit = "nM"
    else:
        conc_unit = "μg/mL"

    st.session_state['time_label'] = time_label
    st.session_state['conc_unit'] = conc_unit

    unique_concs = sorted(df[conc_col].unique())
    num_concs = len(unique_concs)
    # Data Load 모드와 동일하게: interp 있으면 그값, 없으면 Data Load 샘플 기준(8 points, N=50)
    if interp_results is not None and interp_results.get('data_points_per_concentration') is not None:
        points_per_conc = interp_results['data_points_per_concentration']
    elif num_concs > 0:
        points_per_conc = int(df.groupby(conc_col).size().iloc[0])
    else:
        points_per_conc = 8  # Data Load 샘플 raw_enzyme.csv 기준
    if interp_results is not None and interp_results.get('n_replicates') is not None:
        n_replicates = interp_results['n_replicates']
    else:
        n_replicates = st.session_state.get('n_replicates', 50)  # Data Load 샘플 기준 N=50
        if n_replicates == 'N/A' or n_replicates is None:
            n_replicates = 50
    
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Number of Concentrations", num_concs)
    with col2:
        st.metric("Data Points per Concentration", points_per_conc)
    with col3:
        st.metric("Reaction Time", time_display)
    with col4:
        st.metric("N (Number of Replicates)", n_replicates)
    
    with st.expander("Concentration Data Information", expanded=False):
        if num_concs > 0:
            alpha_col = 'alpha' if 'alpha' in df.columns else 'alpha_temp'
            preview_cols = [conc_col, 'time_s', alpha_col]
            preview_df = df[[c for c in preview_cols if c in df.columns]].copy()
            rename_map = {'time_s': time_label, 'alpha': 'α (Normalized)', 'alpha_temp': 'α (Normalized)'}
            preview_df = preview_df.rename(columns={k: v for k, v in rename_map.items() if k in preview_df.columns})
            st.dataframe(preview_df, use_container_width=True, hide_index=True, height=400)
        else:
            st.info("No data available.")
    
    # Tabs for different views
    tab1, tab_alpha, tab_evsalpha, tab2, tab3, tab4 = st.tabs([
        "📊 v₀ vs [S] Fit",
        "📈 Alpha Calculation",
        "📊 [E] vs α Plot",
        "🔬 Model Fitting",
        "📉 Model Comparison",
        "💡 Diagnostic Analysis"
    ])
    
    with tab1:
        # v0 vs [S] Michaelis-Menten Fit Graph
        st.subheader("v₀ vs [S] Michaelis-Menten Fit")
        
        # Data preparation
        v0_data = None
        mm_fit = None
        norm_results_data = None
        
        # 1. Try from session state (Memory from Data Load mode)
        if 'interpolation_results' in st.session_state and st.session_state.get('mm_data_ready', False):
            results = st.session_state['interpolation_results']
            if 'v0_vs_concentration' in results and 'mm_fit_results' in results:
                v0_data = results['v0_vs_concentration']
                mm_fit = results['mm_fit_results']
            
            # 정규화 결과 데이터 가져오기 (농도 값으로 변환)
            if 'normalization_results' in results:
                norm_results_raw = results['normalization_results']
                norm_results_data = {}
                for conc_name, data in norm_results_raw.items():
                    if 'concentration' in data:
                        conc_val = float(data['concentration'])
                        norm_results_data[conc_val] = {
                            'concentration': conc_val,
                            'F0': data.get('F0', None),
                            'Fmax': data.get('Fmax', None),
                            'k_obs': data.get('k_obs', None),
                            'tau': data.get('tau', None),
                            'R_squared': data.get('R_squared', None),
                            'equation': data.get('equation', None)
                        }
        
        # 2. Try from session state (Loaded from file in this mode)
        if (v0_data is None or mm_fit is None) and 'v0_data_from_file' in st.session_state:
            v0_data = st.session_state['v0_data_from_file']
            if 'mm_fit_from_file' in st.session_state:
                mm_fit = st.session_state['mm_fit_from_file']

        # 3. Fallback: 기본 예시 샘플 — 프로젝트 내 calibration xlsx 또는 내장 샘플
        if (v0_data is None or mm_fit is None) and uploaded_file is None:
            _sample_paths = [
                Path(__file__).parent.parent / 'Michaelis-Menten_calibration_results.xlsx',
                Path('Michaelis-Menten_calibration_results.xlsx'),
            ]
            for _path in _sample_paths:
                if not getattr(_path, 'exists', lambda: os.path.exists(str(_path)))():
                    continue
                try:
                    xl = pd.ExcelFile(str(_path), engine='openpyxl')
                    mm_sheet = next((s for s in ['Model simulation input', 'Analysis mode input', 'MM Results'] if s in xl.sheet_names), None)
                    fit_sheet = next((s for s in ['Michaelis-Menten Fit Results', 'Fit results', 'MM Fit Results'] if s in xl.sheet_names), None)
                    if not mm_sheet or not fit_sheet:
                        continue
                    df_mm = pd.read_excel(str(_path), sheet_name=mm_sheet, engine='openpyxl')
                    conc_col = 'Concentration [ug/mL]' if 'Concentration [ug/mL]' in df_mm.columns else 'Concentration'
                    v0_col = next((c for c in ['v0', 'v0 (RFU/min)', 'v₀ (RFU/min)'] if c in df_mm.columns), None)
                    if conc_col in df_mm.columns and v0_col:
                        v0_concs = df_mm[conc_col].dropna().astype(float).tolist()
                        v0_vals = df_mm[v0_col].dropna().astype(float).tolist()
                        if len(v0_concs) == len(v0_vals) and len(v0_concs) > 0:
                            v0_data = {'concentrations': v0_concs, 'v0_values': v0_vals}
                            st.session_state['v0_data_from_file'] = v0_data
                    df_fit = pd.read_excel(str(_path), sheet_name=fit_sheet, engine='openpyxl')
                    p_col = 'Parameter' if 'Parameter' in df_fit.columns else '파라미터'
                    v_col = 'Value' if 'Value' in df_fit.columns else '값'
                    if p_col in df_fit.columns and v_col in df_fit.columns:
                        params = dict(zip(df_fit[p_col], df_fit[v_col]))
                        def _p(*keys):
                            for k in keys:
                                f = next((x for x in params if k.lower() in str(x).lower()), None)
                                if f is not None:
                                    return params[f]
                            return None
                        vmax, km, r2 = _p('Vmax'), _p('Km'), _p('R²', 'R2', 'R_squared')
                        slope, intercept = _p('Slope'), _p('Intercept')
                        try:
                            vmax = float(vmax) if vmax is not None else None
                            km = float(km) if km is not None else None
                        except (TypeError, ValueError):
                            vmax, km = None, None
                        if vmax is not None and km is not None:
                            mm_fit = {'Vmax': vmax, 'Km': km, 'R_squared': float(r2) if r2 is not None else None, 'fit_success': True,
                                      'experiment_type': 'Substrate Concentration Variation (Standard Michaelis-Menten)',
                                      'equation': f'v₀ = {vmax:.2f}[S] / ({km:.2f} + [S])'}
                            st.session_state['mm_fit_from_file'] = mm_fit
                            break
                        if slope is not None:
                            mm_fit = {'slope': float(slope), 'intercept': float(intercept) if intercept is not None else 0,
                                      'R_squared': float(r2) if r2 is not None else None, 'fit_success': True,
                                      'experiment_type': 'Enzyme Concentration Variation (Fixed substrate)',
                                      'equation': f'v₀ = {float(slope):.4f}[E] + {float(intercept) if intercept else 0:.4f}'}
                            st.session_state['mm_fit_from_file'] = mm_fit
                            break
                except Exception:
                    continue
                break

        if v0_data is None or mm_fit is None:
            if 'sample_v0_data' not in st.session_state:
                st.session_state['sample_v0_data'] = {'concentrations': [0.31, 0.63, 1.25, 2.5, 5.0], 'v0_values': [16000., 25800., 35500., 30000., 35200.]}
                st.session_state['sample_mm_fit'] = {'slope': 2808.79, 'intercept': 23046.9, 'R_squared': 0.4266, 'fit_success': True,
                    'experiment_type': 'Enzyme Concentration Variation (Fixed substrate)', 'equation': 'v₀ = 2808.7897[E] + 23046.8969'}
            v0_data = st.session_state.get('sample_v0_data')
            mm_fit = st.session_state.get('sample_mm_fit')
            if v0_data and mm_fit:
                st.caption("📌 Showing **example sample data** (Enzyme concentration variation). Run Data Load Mode or load a result file for your own data.")

        # 4. 파일에서 정규화 결과 읽기 (Normalization Results 시트 또는 MM Results 시트)
        if norm_results_data is None:
            xlsx_path_for_norm = None
            if uploaded_file is not None:
                import tempfile
                file_extension = uploaded_file.name.split('.')[-1].lower()
                if file_extension == 'xlsx':
                    with tempfile.NamedTemporaryFile(delete=False, suffix='.xlsx', mode='wb') as tmp_file:
                        tmp_file.write(uploaded_file.getbuffer())
                        xlsx_path_for_norm = tmp_file.name
            else:
                xlsx_paths = [
                    'Michaelis-Menten_calibration_results.xlsx',
                    str(Path(__file__).parent.parent / 'Michaelis-Menten_calibration_results.xlsx'),
                ]
                for path in xlsx_paths:
                    if os.path.exists(path):
                        xlsx_path_for_norm = path
                        break
            
            if xlsx_path_for_norm is not None:
                try:
                    # Normalization Results 시트 시도
                    xl = pd.ExcelFile(xlsx_path_for_norm)
                    if 'Normalization Results' in xl.sheet_names or 'Normalization results' in xl.sheet_names:
                        norm_sheet = 'Normalization Results' if 'Normalization Results' in xl.sheet_names else 'Normalization results'
                        df_norm = pd.read_excel(xlsx_path_for_norm, sheet_name=norm_sheet, engine='openpyxl')
                        norm_results_data = {}
                        for _, row in df_norm.iterrows():
                            # 농도 추출 (농도 컬럼 찾기)
                            conc_col = None
                            for col in df_norm.columns:
                                if '농도' in col or 'Concentration' in col:
                                    conc_col = col
                                    break
                            
                            if conc_col and pd.notna(row.get(conc_col)):
                                try:
                                    conc_val = float(str(row[conc_col]).replace('μM', '').replace('μg/mL', '').strip())
                                    norm_results_data[conc_val] = {
                                        'concentration': conc_val,
                                        'F0': row.get('F₀', row.get('F0', None)),
                                        'Fmax': row.get('F_max', row.get('Fmax', None)),
                                        'k_obs': row.get('k_obs', None),
                                        'tau': row.get('τ', row.get('tau', None)),
                                        'R_squared': row.get('R²', row.get('R_squared', None)),
                                        'equation': row.get('방정식', row.get('equation', None))
                                    }
                                except (ValueError, TypeError):
                                    continue
                    # Normalization Results 시트가 없으면 Model simulation input / MM Results 시트에서 가져오기
                    elif 'Model simulation input' in xl.sheet_names or 'Analysis mode input' in xl.sheet_names or 'MM Results' in xl.sheet_names:
                        fallback_sheet = 'Model simulation input' if 'Model simulation input' in xl.sheet_names else ('Analysis mode input' if 'Analysis mode input' in xl.sheet_names else 'MM Results')
                        df_mm = pd.read_excel(xlsx_path_for_norm, sheet_name=fallback_sheet, engine='openpyxl')
                        norm_results_data = {}
                        conc_col_name = 'Concentration [ug/mL]' if 'Concentration [ug/mL]' in df_mm.columns else 'Concentration'
                        for _, row in df_mm.iterrows():
                            if pd.notna(row.get(conc_col_name)):
                                try:
                                    conc_val = float(row[conc_col_name])
                                    norm_results_data[conc_val] = {
                                        'concentration': conc_val,
                                        'F0': row.get('F0', None),
                                        'Fmax': row.get('Fmax', None),
                                        'k_obs': row.get('k_obs', None),
                                        'tau': row.get('tau', row.get('τ', None)),
                                        'R_squared': row.get('R²', row.get('R_squared', None)),
                                        'equation': row.get('방정식', row.get('equation', None))
                                    }
                                except (ValueError, TypeError):
                                    continue
                except Exception:
                    pass
                finally:
                    if uploaded_file is not None and xlsx_path_for_norm and os.path.exists(xlsx_path_for_norm):
                        try:
                            os.unlink(xlsx_path_for_norm)
                        except:
                            pass

        # Plotting
        if v0_data and mm_fit:
            # Determine exp type
            exp_type = mm_fit.get('experiment_type', 'Substrate Concentration Variation (Standard MM)')
            
            fig_v0 = go.Figure()
            
            # Experimental Points
            fig_v0.add_trace(go.Scatter(
                x=v0_data['concentrations'],
                y=v0_data['v0_values'],
                mode='markers',
                name='Experimental v₀',
                marker=dict(size=10, color='red', line=dict(width=2, color='black'))
            ))
            
            # Fit Line
            if mm_fit.get('fit_success'):
                conc_min = min(v0_data['concentrations'])
                conc_max = max(v0_data['concentrations'])
                conc_range = np.linspace(conc_min * 0.5, conc_max * 1.5, 200)
                
                if exp_type in ("Substrate 농도 변화 (표준 MM)", "Substrate Concentration Variation (Standard MM)") and mm_fit.get('Vmax') is not None:
                     v0_fitted = michaelis_menten_calibration(conc_range, mm_fit['Vmax'], mm_fit['Km'])
                     line_name = mm_fit.get('equation', 'MM Fit')
                     
                     # Stats text
                     stats_text = f"Vmax = {mm_fit['Vmax']:.2f}<br>"
                     stats_text += f"Km = {mm_fit['Km']:.4f} μM<br>"
                     if mm_fit.get('R_squared'):
                        stats_text += f"R² = {mm_fit['R_squared']:.4f}"
                     
                     xaxis_title = '[S] (μM)'
                     title = 'Initial Velocity (v₀) vs Substrate Concentration [S]'
                     
                else: # Linear/Enzyme
                     slope = mm_fit.get('slope', 0)
                     intercept = mm_fit.get('intercept', 0)
                     v0_fitted = slope * conc_range + intercept
                     line_name = mm_fit.get('equation', 'Linear Fit')
                     
                     # Stats text
                     stats_text = f"Slope = {slope:.4f}<br>"
                     stats_text += f"Intercept = {intercept:.4f}<br>"
                     if mm_fit.get('R_squared'):
                        stats_text += f"R² = {mm_fit['R_squared']:.4f}<br>"
                     stats_text += "<br><b>⚠️ Cannot calculate Km</b>"
                     
                     xaxis_title = '[E] (μg/mL)'
                     title = 'Initial Velocity (v₀) vs Enzyme Concentration [E] (Constant Substrate)'

                fig_v0.add_trace(go.Scatter(
                    x=conc_range,
                    y=v0_fitted,
                    mode='lines',
                    name=line_name,
                    line=dict(width=2.5, color='blue')
                ))
                
                fig_v0.add_annotation(
                    xref="paper", yref="paper",
                    x=0.05, y=0.95,
                    xanchor='left', yanchor='top',
                    text=stats_text,
                    showarrow=False,
                    bgcolor="rgba(0,0,0,0)",
                    bordercolor="rgba(0,0,0,0)",
                    borderwidth=0,
                    font=dict(size=11)
                )

            fig_v0.update_layout(
                title=title,
                xaxis_title=xaxis_title,
                yaxis_title='Initial Velocity v₀ (Fluorescence Units / Time)',
                template='plotly_white',
                height=500,
                hovermode='x unified'
            )
            st.plotly_chart(fig_v0, use_container_width=True)
            
            # Show table with additional columns
            st.subheader("📋 Experimental Data")
            
            # 테이블 데이터 준비
            table_data = {
                xaxis_title: v0_data['concentrations'],
                'v₀ (RFU/min)': v0_data['v0_values']
            }
            
            # 정규화 결과 데이터 추가
            if norm_results_data:
                fmax_list = []
                r2_list = []
                k_obs_list = []
                tau_list = []
                equation_list = []
                
                for conc in v0_data['concentrations']:
                    # 농도 매칭 (부동소수점 오차 고려)
                    matched_data = None
                    for norm_conc, norm_data in norm_results_data.items():
                        if abs(float(conc) - float(norm_conc)) < 0.001:
                            matched_data = norm_data
                            break
                    
                    if matched_data:
                        fmax_list.append(matched_data.get('Fmax', None))
                        r2_list.append(matched_data.get('R_squared', None))
                        k_obs_list.append(matched_data.get('k_obs', None))
                        tau_list.append(matched_data.get('tau', None))
                        equation_list.append(matched_data.get('equation', None))
                    else:
                        fmax_list.append(None)
                        r2_list.append(None)
                        k_obs_list.append(None)
                        tau_list.append(None)
                        equation_list.append(None)
                
                table_data['Fmax'] = fmax_list
                table_data['R²'] = r2_list
                table_data['k_obs'] = k_obs_list
                table_data['τ'] = tau_list
                table_data['Equation'] = equation_list
            
            df_table = pd.DataFrame(table_data).sort_values(xaxis_title)
            st.dataframe(df_table, use_container_width=True, hide_index=True)
                 
        else:
            st.info("⚠️ No Michaelis-Menten fitting data. Run analysis in Data Load Mode or load a result file that includes the 'Michaelis-Menten Fit Results' sheet.")
    
    with tab_alpha:
        st.subheader("📈 Alpha (α) Calculation")

        with st.expander("📖 What is Alpha (α)? — definition and formula", expanded=False):
            st.markdown("""
            **What is Alpha (α)?**
            Normalized cleavage ratio, ranging from 0 (no cleavage) to 1 (full cleavage).

            **Formula**: α(t) = (F(t) - F₀) / (Fmax - F₀)
            - **F(t)**: Fluorescence at time t
              - With Data Load Mode results: interpolated values from the normalization exponential curve (RFU_Interpolated)
              - Direct calculation: raw fluorescence from the data
            - **F₀**: Initial fluorescence
              - If available from Data Load Mode: F0 from interpolated values (MM Results sheet)
              - Otherwise: minimum fluorescence per concentration (min(F))
            - **Fmax**: Maximum fluorescence
              - If available from Data Load Mode: Fmax from interpolated values (MM Results sheet)
              - Otherwise: Region-based normalization
                1. If a plateau region exists: mean fluorescence of the plateau (mean(F_plateau))
                2. If sufficient exponential rise (≥3 points): F∞ from exponential fit (F(t) = F₀ + A·(1 - e^(-k·t)), Fmax = F₀ + A)
                3. Otherwise: maximum fluorescence (max(F))
            """)

        # Check if alpha column exists
        if 'alpha' not in df.columns:
            st.error("❌ Alpha was not computed. Data normalization is required.")
            st.info("💡 Data is not normalized. Please check the data load and normalization steps.")
        else:
            # Alpha vs Time Plot
            st.subheader("📊 Normalized Data: α(t) vs Time")
            
            fig_alpha = Visualizer.plot_normalized_data(df, conc_unit, time_label, 
                                                       use_lines=True,
                                                       enzyme_name=enzyme_name,
                                                       substrate_name=substrate_name,
                                                       experiment_type=experiment_type)
            # 원본 시간 범위로 xaxis 설정
            original_time_max = st.session_state.get('original_time_max', df['time_s'].max())
            if time_unit == 'min':
                fig_alpha.update_xaxes(range=[0, original_time_max])
            else:
                fig_alpha.update_xaxes(range=[0, original_time_max])
            st.plotly_chart(fig_alpha, use_container_width=True)
            
            # Alpha Statistics
            st.subheader("📋 Alpha Statistics by Concentration")
            
            conc_col = 'enzyme_ugml' if 'enzyme_ugml' in df.columns else df['conc_col_name'].iloc[0] if 'conc_col_name' in df.columns else None
            
            if conc_col:
                alpha_stats = []
                for conc in sorted(df[conc_col].unique()):
                    subset = df[df[conc_col] == conc]
                    alpha_stats.append({
                        f'Concentration ({conc_unit})': conc,
                        'Alpha mean': f"{subset['alpha'].mean():.4f}",
                        'Alpha std': f"{subset['alpha'].std():.4f}"
                    })
                
                st.dataframe(pd.DataFrame(alpha_stats), use_container_width=True, hide_index=True)
            
            # F0, Fmax 정보
            st.subheader("🔬 Normalization Parameters (F₀, Fmax)")
            
            # Check if fitted parameters are being used
            fitted_params_used = st.session_state.get('fitted_params', None)
            using_fitted_params = fitted_params_used is not None and len(fitted_params_used) > 0
            
            if using_fitted_params:
                st.success(f"✅ F0, Fmax parameters loaded ({len(fitted_params_used)} concentrations)")
                st.info("💡 F0, Fmax are constants from the normalization exponential in Data Load Mode.")
                st.info("📊 Formula used: F(t) = F₀ + (Fmax - F₀)·[1 - exp(-k_obs·t)]")
            else:
                st.info("ℹ️ Using default normalization (Region-based calculation)")
            
            # F0, Fmax 테이블
            if conc_col and 'F0' in df.columns and 'Fmax' in df.columns:
                f0_fmax_data = []
                for conc in sorted(df[conc_col].unique()):
                    subset = df[df[conc_col] == conc]
                    fmax_method = subset['Fmax_method'].iloc[0] if 'Fmax_method' in subset.columns else "N/A"
                    
                    # F0, Fmax 값의 출처 확인
                    df_F0 = subset['F0'].iloc[0]
                    df_Fmax = subset['Fmax'].iloc[0]
                    
                    # fitted_params에서 값 확인
                    source_info = "Region-based calculation"
                    if using_fitted_params and fitted_params_used:
                        # 농도 매칭 (부동소수점 오차 고려)
                        matched_conc = None
                        for fitted_conc in fitted_params_used.keys():
                            if abs(float(conc) - float(fitted_conc)) < 0.001:
                                matched_conc = fitted_conc
                                break
                        
                        if matched_conc:
                            fitted_F0 = fitted_params_used[matched_conc]['F0']
                            fitted_Fmax = fitted_params_used[matched_conc]['Fmax']
                            
                            # 값이 일치하는지 확인
                            if abs(df_F0 - fitted_F0) < 0.01 and abs(df_Fmax - fitted_Fmax) < 0.01:
                                if fmax_method == 'fitted_from_data_load':
                                    source_info = "Normalization exponential"
                                else:
                                    source_info = "Normalization exponential (from normalization step)"
                            else:
                                source_info = f"From normalization (F0={fitted_F0:.2f}, Fmax={fitted_Fmax:.2f})"
                    
                    f0_fmax_data.append({
                        f'Concentration ({conc_unit})': conc,
                        'F₀ (initial)': f"{df_F0:.2f}",
                        'Fmax (max)': f"{df_Fmax:.2f}",
                        'Fmax method': fmax_method,
                        'Source': source_info
                    })
                
                st.dataframe(pd.DataFrame(f0_fmax_data), use_container_width=True, hide_index=True)
            
            # Normalization method description
            with st.expander("📖 Normalization method details", expanded=False):
                if using_fitted_params:
                    st.markdown("""
                    **Using F0, Fmax from normalization exponential:**
                    - F0, Fmax: Constants from the normalization exponential in Data Load Mode
                    - **Formula**: F(t) = F₀ + (Fmax - F₀)·[1 - exp(-k_obs·t)]
                      - F₀: Initial fluorescence (from normalization)
                      - Fmax: Maximum fluorescence (from normalization)
                      - k_obs: Observed rate constant
                    - **Alpha**: α(t) = (F(t) − F₀) / (Fmax − F₀)
                    - Values from iterative normalization in Data Load Mode
                    """)
                else:
                    st.markdown("""
                    **Default normalization (Region-based):**
                    
                    1. **Temporary normalization**
                       - F0_temp = min(F), Fmax_temp = max(F)
                       - α_temp = (F - F0_temp) / (Fmax_temp - F0_temp)
                    
                    2. **Region division**
                       - Initial linear region
                       - Exponential growth region
                       - Plateau region
                    
                    3. **Final normalization**
                       - F0 = F0_temp
                       - Fmax: plateau mean if available; else exponential fit for F∞; else max(F)
                       - α = (F - F0) / (Fmax - F0)
                    
                    **Fmax methods:**
                    - `plateau_avg`: Plateau mean
                    - `exponential_fit`: F∞ from exponential fit
                    - `fallback_max`: max(F) (fallback)
                    """)
            
            # Download Alpha data
            st.subheader("💾 Download Alpha Data")
            
            # Alpha 데이터 준비
            alpha_download_df = df[['time_s', conc_col, 'alpha', 'F0', 'Fmax']].copy() if conc_col else df[['time_s', 'alpha', 'F0', 'Fmax']].copy()
            alpha_download_df = alpha_download_df.sort_values(['time_s', conc_col] if conc_col else 'time_s')
            
            csv_alpha = alpha_download_df.to_csv(index=False)
            st.download_button(
                label="📥 Download Alpha Data (CSV)",
                data=csv_alpha,
                file_name="alpha_calculation_results.csv",
                mime="text/csv",
                help="CSV with time, concentration, alpha, F0, Fmax"
            )

    with tab_evsalpha:
        st.subheader("📊 [E] vs α Plot")
        st.markdown("Enzyme concentration [E] vs cleavage fraction α (α mean per concentration with fit).")
        # 데이터 소스 표시 (로드된 데이터에 따라 플롯이 바뀌는지 확인용)
        if st.session_state.get('mm_data_ready') and st.session_state.get('interpolation_results'):
            st.caption("📌 Data source: **Data Load result (in memory)** — change raw data in Data Load and run again to update this plot.")
        elif st.session_state.get('data_source_type'):
            st.caption(f"📌 Data source: **{st.session_state.get('data_source_type', 'File')}**")
        if 'alpha' not in df.columns:
            st.info("ℹ️ Compute alpha in the **Alpha Calculation** tab first.")
        else:
            conc_col_ev = 'enzyme_ugml' if 'enzyme_ugml' in df.columns else (df['conc_col_name'].iloc[0] if 'conc_col_name' in df.columns else None)
            if conc_col_ev is None:
                st.warning("No concentration column found.")
            else:
                agg = df.groupby(conc_col_ev)['alpha'].agg(['max', 'mean', 'min']).reset_index()
                agg = agg.rename(columns={'max': 'α max', 'mean': 'α mean', 'min': 'α min'})
                x_label = f"[E] ({conc_unit})"
                x_ev = agg[conc_col_ev].values.astype(float)
                y_mean = agg['α mean'].values.astype(float)
                n_pts = len(x_ev)
                ss_tot = np.sum((y_mean - np.mean(y_mean)) ** 2) if n_pts > 1 else 1.0
                fig_ev = go.Figure()
                # α mean only (no α max)
                fig_ev.add_trace(go.Scatter(
                    x=x_ev,
                    y=y_mean,
                    mode='markers',
                    name='α mean',
                    marker=dict(size=10, color='#ff7f0e', symbol='diamond')
                ))
                x_smooth = np.linspace(0, max(x_ev) * 1.05, 200)
                fit_results = []  # (name, R², AIC, params_text, y_fit, color, dash)
                # 1) Exponential: α = α_max (1 - e^{-k[E]})
                def _exp_sat(E, a_max, k):
                    return a_max * (1 - np.exp(-k * np.maximum(E, 0)))
                try:
                    p0 = (float(np.max(y_mean)), 1.0)
                    bounds = ([0.001, 0.001], [2.0, 1e4])
                    popt, _ = curve_fit(_exp_sat, x_ev, y_mean, p0=p0, bounds=bounds, maxfev=5000)
                    a_max, k = float(popt[0]), float(popt[1])
                    y_exp = _exp_sat(x_smooth, a_max, k)
                    ss_res = np.sum((y_mean - _exp_sat(x_ev, a_max, k)) ** 2)
                    r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
                    aic = n_pts * np.log(ss_res / n_pts + 1e-12) + 2 * 2 if n_pts > 0 else 0
                    fit_results.append((
                        'Exponential',
                        r2, aic,
                        f'α_max = {a_max:.4f}, k = {k:.4f}',
                        y_exp, '#2ca02c', 'dash'
                    ))
                except Exception:
                    pass
                # 2) Hyperbolic (MM-type): α = α_max [E] / (K_half + [E])
                def _hyp_sat(E, a_max, K_half):
                    E_safe = np.maximum(E, 0)
                    return a_max * E_safe / (K_half + E_safe + 1e-12)
                try:
                    p0 = (float(np.max(y_mean)), 0.5 * (np.max(x_ev) or 1))
                    bounds = ([0.001, 0.001], [2.0, 1e4])
                    popt, _ = curve_fit(_hyp_sat, x_ev, y_mean, p0=p0, bounds=bounds, maxfev=5000)
                    a_max, K_half = float(popt[0]), float(popt[1])
                    y_hyp = _hyp_sat(x_smooth, a_max, K_half)
                    ss_res = np.sum((y_mean - _hyp_sat(x_ev, a_max, K_half)) ** 2)
                    r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
                    aic = n_pts * np.log(ss_res / n_pts + 1e-12) + 2 * 2 if n_pts > 0 else 0
                    fit_results.append((
                        'Hyperbolic (MM)',
                        r2, aic,
                        f'α_max = {a_max:.4f}, K_half = {K_half:.4f}',
                        y_hyp, '#9467bd', 'dot'
                    ))
                except Exception:
                    pass
                for name, r2, aic, params_text, y_fit, color, dash in fit_results:
                    label = f'{name}: R²={r2:.4f}'
                    fig_ev.add_trace(go.Scatter(
                        x=x_smooth,
                        y=y_fit,
                        mode='lines',
                        name=label,
                        line=dict(width=2, color=color, dash=dash)
                    ))
                fig_ev.update_layout(
                    title='α vs [E] (Enzyme Concentration)',
                    xaxis_title=x_label,
                    yaxis_title='α (cleavage fraction)',
                    template='plotly_white',
                    height=500,
                    hovermode='x unified',
                    legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='right', x=1)
                )
                st.plotly_chart(fig_ev, use_container_width=True)
                if fit_results:
                    st.markdown("**Fit models**")
                    st.latex(r"\text{Exponential: } \alpha = \alpha_{max}(1 - e^{-k[E]}) \quad \text{Hyperbolic: } \alpha = \frac{\alpha_{max}[E]}{K_{half}+[E]}")
                    for name, r2, aic, params_text, _, _, _ in fit_results:
                        st.caption(f"**{name}**: {params_text} → R² = {r2:.4f}, AIC = {aic:.2f}")
                    if len(fit_results) == 2:
                        aic_exp, aic_hyp = fit_results[0][2], fit_results[1][2]
                        if aic_exp < aic_hyp:
                            st.caption("Lower AIC → **Exponential** fits slightly better.")
                        elif aic_hyp < aic_exp:
                            st.caption("Lower AIC → **Hyperbolic** fits slightly better.")
                        else:
                            st.caption("Similar AIC → choose by interpretation (K_half = half-saturation [E] for Hyperbolic).")
                with st.expander("📋 [E] vs α data", expanded=False):
                    st.dataframe(agg, use_container_width=True, hide_index=True)
    
    with tab2:
        st.subheader("🔬 Global Model Fitting")
        
        with st.expander("📚 Kinetic Model Description", expanded=False):
            st.markdown(r"""
            This simulator provides six kinetic models for analyzing peptide substrate–enzyme reactions.

            #### 1. Basic Models (A–C)
            Based on classical enzyme kinetics; Fmax is assumed independent of enzyme concentration.

            **📌 Model A: Substrate Depletion**
            - **Overview**: First-order reaction; rate decreases as substrate [S] is consumed.
            - **Equations**:
              $$ \alpha(t) = 1 - e^{-k_{obs} \cdot t} $$
              $$ k_{obs} = \frac{k_{cat}}{K_M} \cdot [E] $$
            - **Parameters**:
              - $k_{cat}/K_M$: Catalytic efficiency ($M^{-1}s^{-1}$)
            - **Notes**:
              - At low [E], initial rate v₀ is linear in [E].
              - At long times, all substrate is cleaved and α → 1.

            **📌 Model B: Enzyme Deactivation**
            - **Overview**: Enzyme gradually loses activity during the reaction.
            - **Equations**:
              Enzyme concentration decays exponentially: $[E]_t = [E]_0 \cdot e^{-k_d t}$
              $$ \alpha(t) = 1 - \exp\left[-\frac{k_{cat}/K_M \cdot [E]_0}{k_d} (1 - e^{-k_d t})\right] $$
            - **Parameters**:
              - $k_{cat}/K_M$: Catalytic efficiency ($M^{-1}s^{-1}$)
              - $k_d$: Deactivation rate constant ($s^{-1}$)
            - **Notes**: Curve plateaus earlier than expected; reaction can stop with substrate left ($\alpha_{\infty} < 1$).

            **📌 Model C: Mass Transfer Limitation**
            - **Overview**: Diffusion of enzyme to the substrate surface is slower than the reaction.
            - **Equations**:
              Surface enzyme concentration $[E]_s$ is lower than bulk $[E]_b$
              $$ [E]_s \approx \frac{[E]_b}{1 + Da}, \quad Da = \frac{k_{cat} \Gamma_0}{K_M k_m} $$
              $$ \alpha(t) = 1 - e^{-k_{obs} \cdot t}, \quad k_{obs} = \frac{k_{cat}}{K_M} [E]_s $$
            - **Parameters**:
              - $k_{cat}/K_M$: Catalytic efficiency
              - $k_m$: Mass transfer coefficient ($m/s$)
              - $\Gamma_0$: Initial surface substrate density ($pmol/cm^2$)
            - **Notes**: At high [E], rate saturates due to diffusion limit.

            ---

            #### 2. Extended Models (D–F)
            Fmax (maximum fluorescence) depends on [E]; for more complex surface reactions.

            **📌 Model D: Concentration-Dependent Fmax**
            - **Overview**: Higher [E] allows access to more substrate (e.g. deeper gel penetration).
            - **Equations**:
              Maximum cleavage $\alpha_{max}$ depends on [E]
              $$ \alpha(t) = \alpha_{max}([E]) \cdot (1 - e^{-k_{obs} t}) $$
              $$ \alpha_{max}([E]) = \alpha_{\infty} \cdot (1 - e^{-k_{access} [E]}) $$
            - **Parameters**:
              - $k_{cat}/K_M$: Catalytic efficiency
              - $\alpha_{\infty}$: Theoretical maximum cleavage
              - $k_{access}$: Accessibility coefficient ($M^{-1}$)
            - **Notes**: At low [E] only surface is cleaved; at high [E], Fmax increases.

            **📌 Model E: Product Inhibition**
            - **Overview**: Reaction product binds the active site and inhibits the enzyme.
            - **Equations**:
              Competitive inhibition
              $$ \frac{d\alpha}{dt} = \frac{k_{obs}(1-\alpha)}{1 + K_{i,eff}\alpha} $$
            - **Parameters**:
              - $k_{cat}/K_M$: Catalytic efficiency
              - $K_{i,eff}$: Effective inhibition constant (dimensionless, $[S]_0/K_i$)
            - **Notes**: Rate drops sharply in the later phase.

            **📌 Model F: Enzyme Surface Sequestration**
            - **Overview**: Enzyme is irreversibly adsorbed to the surface or gel and unavailable for reaction.
            - **Equations**:
              Enzyme depleted by surface adsorption ($k_{ads}$)
              $$ \alpha(t) \approx \frac{(k_{cat}/K_M)[E]}{k_{ads}(1+K_{ads}[E])} (1 - e^{-k_{ads} t}) $$
            - **Parameters**:
              - $k_{cat}/K_M$: Catalytic efficiency
              - $k_{ads}$: Adsorption rate constant ($s^{-1}$)
              - $K_{ads}$: Adsorption equilibrium constant ($M^{-1}$)
            - **Notes**: Reactivity can be lower than expected even at high [E].

            ### 📊 Model fit: AIC

            **Akaike Information Criterion (AIC)**
            Balances goodness of fit and model complexity; lower is better.

            **Formula**:
            $$ AIC = 2k - 2\ln(\hat{L}) $$
            where:
            - $k$: number of parameters
            - $\hat{L}$: maximum likelihood

            Here we use RSS-based:
            $$ AIC = n \ln\left(\frac{RSS}{n}\right) + 2k + C $$
            - $n$: number of data points
            - $RSS$: residual sum of squares ($\sum (y_{obs} - y_{pred})^2$)
            - $C$: constant

            **Interpretation**:
            - **ΔAIC < 2**: No significant difference between models
            - **ΔAIC > 10**: Lower-AIC model is strongly preferred
            """)
        
        st.markdown("**Basic models (A–C)**: Classical enzyme kinetics")
        
        # Model selection — 한동안 Model A만 활성화, 나머지 B–F 숨김
        fit_model_a = st.checkbox("Model A: Substrate Depletion", value=True)
        st.caption("✓ First-order reaction & substrate depletion")
        fit_model_b = False
        fit_model_c = False
        fit_model_d = False
        fit_model_e = False
        fit_model_f = False
        
        if st.button("🚀 Run Global Fitting", type="primary"):
            # 데이터 상태 확인 및 검증
            with st.expander("🔍 Data status", expanded=False):
                st.write("**Required columns:**")
                required_cols = ['alpha', 'time_s', 'FL_intensity']
                missing_cols = [col for col in required_cols if col not in df.columns]
                if missing_cols:
                    st.error(f"❌ Missing columns: {missing_cols}")
                    st.stop()
                else:
                    st.success(f"✅ Required columns present: {required_cols}")
                
                st.write("**Data statistics:**")
                st.write(f"- Total data points: {len(df)}")
                st.write(f"- Alpha range: {df['alpha'].min():.4f} ~ {df['alpha'].max():.4f}")
                st.write(f"- Alpha mean: {df['alpha'].mean():.4f}")
                st.write(f"- Alpha std: {df['alpha'].std():.4f}")
                st.write(f"- Time range: {df['time_s'].min():.2f} ~ {df['time_s'].max():.2f} s")
                
                conc_col = 'enzyme_ugml' if 'enzyme_ugml' in df.columns else df['conc_col_name'].iloc[0] if 'conc_col_name' in df.columns else None
                if conc_col:
                    st.write("**Alpha statistics by concentration:**")
                    conc_stats = df.groupby(conc_col)['alpha'].agg(['count', 'min', 'max', 'mean', 'std'])
                    st.dataframe(conc_stats, use_container_width=True)
                
                if df['alpha'].max() < 0.1:
                    st.warning("⚠️ All alpha values are below 0.1. Normalization may have failed.")
                if df['alpha'].std() < 0.01:
                    st.warning("⚠️ Alpha variability is very low. Data may not be properly normalized.")
            
            results = []
            
            # Create a status container
            status_container = st.empty()
            result_container = st.container()
            
            # Model A
            if fit_model_a:
                with status_container:
                    with st.spinner("🔄 Fitting Model A: Substrate Depletion..."):
                        model_a = ModelA_SubstrateDepletion(enzyme_mw=enzyme_mw)
                        result_a = model_a.fit_global(df, verbose_callback=verbose_callback)
                        results.append(result_a)
                
                if result_a:
                    with result_container:
                        st.success(f"✅ Model A done: R² = {result_a.r_squared:.4f}, AIC = {result_a.aic:.2f}")
                else:
                    with result_container:
                        st.error("❌ Model A fitting failed")
            
            # Model B
            if fit_model_b:
                with status_container:
                    with st.spinner("🔄 Fitting Model B: Enzyme Deactivation..."):
                        model_b = ModelB_EnzymeDeactivation(enzyme_mw=enzyme_mw)
                        result_b = model_b.fit_global(df, verbose_callback=verbose_callback)
                        results.append(result_b)
                
                if result_b:
                    with result_container:
                        st.success(f"✅ Model B done: R² = {result_b.r_squared:.4f}, AIC = {result_b.aic:.2f}")
                else:
                    with result_container:
                        st.error("❌ Model B fitting failed")
            
            # Model C
            if fit_model_c:
                with status_container:
                    with st.spinner("🔄 Fitting Model C: Mass Transfer Limitation..."):
                        model_c = ModelC_MassTransfer(enzyme_mw=enzyme_mw)
                        result_c = model_c.fit_global(df, verbose_callback=verbose_callback)
                        results.append(result_c)
                
                if result_c:
                    with result_container:
                        st.success(f"✅ Model C done: R² = {result_c.r_squared:.4f}, AIC = {result_c.aic:.2f}")
                else:
                    with result_container:
                        st.error("❌ Model C fitting failed")
            
            # Model D
            if fit_model_d:
                with status_container:
                    with st.spinner("🔄 Fitting Model D: Concentration-Dependent Fmax..."):
                        model_d = ModelD_ConcentrationDependentFmax(enzyme_mw=enzyme_mw)
                        result_d = model_d.fit_global(df, verbose_callback=verbose_callback)
                        results.append(result_d)
                
                if result_d:
                    with result_container:
                        st.success(f"✅ Model D done: R² = {result_d.r_squared:.4f}, AIC = {result_d.aic:.2f}")
                else:
                    with result_container:
                        st.error("❌ Model D fitting failed")
            
            # Model E
            if fit_model_e:
                with status_container:
                    with st.spinner("🔄 Fitting Model E: Product Inhibition..."):
                        model_e = ModelE_ProductInhibition(enzyme_mw=enzyme_mw)
                        result_e = model_e.fit_global(df, verbose_callback=verbose_callback)
                        results.append(result_e)
                
                if result_e:
                    with result_container:
                        st.success(f"✅ Model E done: R² = {result_e.r_squared:.4f}, AIC = {result_e.aic:.2f}")
                else:
                    with result_container:
                        st.error("❌ Model E fitting failed")
            
            # Model F
            if fit_model_f:
                with status_container:
                    with st.spinner("🔄 Fitting Model F: Enzyme Adsorption/Sequestration..."):
                        model_f = ModelF_EnzymeSurfaceSequestration(enzyme_mw=enzyme_mw)
                        result_f = model_f.fit_global(df, verbose_callback=verbose_callback)
                        results.append(result_f)
                
                if result_f:
                    with result_container:
                        st.success(f"✅ Model F done: R² = {result_f.r_squared:.4f}, AIC = {result_f.aic:.2f}")
                else:
                    with result_container:
                        st.error("❌ Model F fitting failed")
            
            # Clear status container after all done
            status_container.empty()
            
            # Store results in session state
            st.session_state['fit_results'] = results
            st.session_state['df'] = df
            
            # Show completion message
            with result_container:
                st.success("🎉 All model fitting complete! Check results in the 'Model Comparison' tab.")
    
    with tab3:
        if 'fit_results' in st.session_state:
            results = st.session_state['fit_results']
            df = st.session_state['df']
            
            st.subheader("📊 Model Comparison")
            
            # Comparison table
            comparison_df = Visualizer.create_comparison_table(results)
            st.dataframe(comparison_df, use_container_width=True)
            
            # Determine best model
            valid_results = [r for r in results if r is not None]
            if valid_results:
                best_aic = min(r.aic for r in valid_results)
                best_model = [r for r in valid_results if r.aic == best_aic][0]
                
                st.success(f"🏆 Best model (lowest AIC): **{best_model.name}** (AIC = {best_model.aic:.2f})")
                
                # Parameter details for best model
                st.subheader(f"Best model parameters: {best_model.name}")
                param_data = []
                for param, value in best_model.params.items():
                    std = best_model.params_std.get(param, 0)
                    param_data.append({
                        'Parameter': param,
                        'Value': f"{value:.4e}",
                        'Std. Error': f"{std:.4e}",
                        'Rel. Error': f"{(std/value*100):.2f}%" if value != 0 else "N/A"
                    })
                st.dataframe(pd.DataFrame(param_data), use_container_width=True)
            
            # Plot all model fits
            st.subheader("📈 Overall Model Fitting Results")
            fig_models = Visualizer.plot_model_fits(df, results, conc_unit, time_label,
                                                    enzyme_name=enzyme_name,
                                                    substrate_name=substrate_name)
            # 원본 시간 범위로 xaxis 설정
            original_time_max = st.session_state.get('original_time_max', df['time_s'].max())
            if time_unit == 'min':
                fig_models.update_xaxes(range=[0, original_time_max])
            else:
                fig_models.update_xaxes(range=[0, original_time_max])
            st.plotly_chart(fig_models, use_container_width=True)
            
            # Individual model plots
            st.subheader("📊 Individual Model Comparison")
            st.markdown("Compare raw data and fitted results for each model.")
            
            # Create tabs for each model
            model_names = [r.name for r in results if r is not None]
            
            if len(model_names) > 0:
                model_tabs_ui = st.tabs(model_names)
                
                for idx, (tab, result) in enumerate(zip(model_tabs_ui, [r for r in results if r is not None])):
                    with tab:
                        # Color scheme for each model
                        model_colors = ['#FF6B6B', '#4ECDC4', '#FFD93D']
                        color = model_colors[idx % len(model_colors)]
                        
                        # Display individual model plot
                        fig_ind = Visualizer.plot_individual_model(df, result, conc_unit, time_label, color)
                        # 원본 시간 범위로 xaxis 설정
                        original_time_max = st.session_state.get('original_time_max', df['time_s'].max())
                        if time_unit == 'min':
                            fig_ind.update_xaxes(range=[0, original_time_max])
                        else:
                            fig_ind.update_xaxes(range=[0, original_time_max])
                        st.plotly_chart(fig_ind, use_container_width=True)
                        
                        # Display parameters
                        st.markdown(f"**{result.name} parameters**")
                        param_cols = st.columns(len(result.params))
                        for col_idx, (param, value) in enumerate(result.params.items()):
                            with param_cols[col_idx]:
                                std = result.params_std.get(param, 0)
                                st.metric(
                                    label=param,
                                    value=f"{value:.4e}",
                                    delta=f"±{std:.4e}" if std > 0 else None
                                )
            
            # Download results
            st.subheader("💾 Download Results")
            csv = comparison_df.to_csv(index=False)
            st.download_button(
                label="Download comparison table (CSV)",
                data=csv,
                file_name="model_comparison.csv",
                mime="text/csv"
            )
        else:
            st.info("👈 Please run fitting in the 'Model Fitting' tab first.")
    
    with tab4:
        st.subheader("💡 Diagnostic Analysis")
        
        # Initial rate analysis
        st.plotly_chart(
            Visualizer.plot_initial_rates(df, conc_unit, time_unit), 
            use_container_width=True
        )
        
        st.markdown("""
        ### 📋 Model selection guidelines
        
        #### Basic models (A–C)
        
        **Model A (Substrate depletion)** when:
        - Initial rate v₀ is linear in [E] at low [E]
        - Saturation fluorescence F∞ ≈ constant (α → 1)
        - No significant enzyme deactivation
        
        **Model B (Enzyme deactivation)** when:
        - F∞ < theoretical maximum (α < 1 at saturation)
        - Fast initial rise then plateau lower than expected
        - kd > 0 with significant contribution
        
        **Model C (Mass transfer)** when:
        - Initial burst (0–5 s) then slower approach
        - Sensitive to stirring/flow
        - v₀ vs [E] saturates at high [E]
        
        #### Extended models (D–F): **Fmax depends on [E]**
        
        **Model D (Concentration-dependent Fmax)** when:
        - α_max increases at higher [E] (more substrate access)
        - Gel penetration depth (thick/dense gel)
        - Secondary cleavage increases product release
        - **Parameters**: α_∞, k_access
        
        **Model E (Product inhibition)** when:
        - Fast initial rise then slowdown (product buildup)
        - Stronger inhibition at low [E]
        - Rate recovers when product is removed
        - **Parameters**: Ki_eff
        
        **Model F (Enzyme adsorption/sequestration)** when:
        - Less affected at high [E] (saturation)
        - Negatively charged surface / PDA coating, dense gel
        - Enzyme activity decreases over time (irreversible)
        - **Parameters**: k_ads, K_ads
        
        ### 📊 Statistics
        - **AIC/BIC**: Lower is better (parameter penalty)
        - **R²**: Higher is better (>0.95 good)
        - **RMSE**: Lower is better
        - **Δ AIC > 10**: Strong evidence against the higher-AIC model
        - **Δ AIC < 2**: No significant difference between models
        """)
        
        # Experimental suggestions
        st.subheader("🧪 Suggested Follow-up Experiments")
        
        st.markdown("""
        ### 🔍 Experiments to check Fmax vs [E]
        
        1. **Long-time measurement at various [E]** (30 min–1 h)
           - Measure saturation fluorescence (Fmax) per concentration
           - Plot [E] vs Fmax → linear vs saturated
           - **Linear increase** → Model D likely
           - **Constant** → Basic models A–C
        
        2. **Gel thickness** (Model D)
           - Thin (50 μm) vs thick (500 μm) gel
           - Thick gel: stronger [E] dependence → diffusion/penetration limit
           - Thin gel: [E]-independent → surface reaction dominant
        
        3. **Product addition** (Model E)
           - Add pre-cleaved peptide
           - Initial rate drops → product inhibition
           - α_max decreases at high [product]
        
        4. **Surface treatment** (Model F)
           - Positive vs negative (PDA) vs neutral (PEG) surface
           - Negative surface: stronger [E] dependence → adsorption
           - PEG: less adsorption → consider D/E
        
        ### 🧬 Classical mechanism tests
        
        5. **Pulse-chase** (Model B)
           - Add fresh enzyme at t=5 min
           - Curve rises again → substrate left (Model A)
           - No change → enzyme deactivation (Model B)
        
        6. **Stirring/flow** (Model C)
           - Static vs rotation (100 rpm) vs flow (1 mL/min)
           - Higher flow → higher α → mass transfer limit
           - No change → reaction-limited (A/B)
        
        7. **Substrate density** (Model A)
           - 0.5×, 1×, 2× peptide immobilization
           - Fmax scales → substrate depletion
           - Fmax unchanged → other mechanism
        
        8. **Solution-phase control**
           - Soluble substrate (same concentration)
           - Full cleavage (α→1) → surface/diffusion issue
           - Incomplete → intrinsic inhibition/deactivation
        
        ### 🎯 Model decision tree
        
        ```
        Does Fmax increase with [E]?
        ├─ YES → Test extended models (D–F)
        │   ├─ Gel thickness sensitive? → Model D (penetration)
        │   ├─ Product addition reduces rate? → Model E (inhibition)
        │   └─ Surface charge sensitive? → Model F (adsorption)
        │
        └─ NO → Test basic models (A–C)
            ├─ Pulse-chase response? → Model A (substrate)
            ├─ α_max decreases with time? → Model B (deactivation)
            └─ Flow sensitive? → Model C (diffusion)
        ```
        """)


