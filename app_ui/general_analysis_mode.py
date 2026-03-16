import sys
from datetime import datetime
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
import numpy as np
import os
import zipfile
from io import BytesIO
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
    fit_progress_curves_pseudo_first_order,
    ProgressCurveFitResult,
    KobsVsELinearResult,
)
_debug_log("module load: importing mode_general_analysis.plot")
from mode_general_analysis.plot import Visualizer
_debug_log("module load: importing mode_prep_raw_data.prep (michaelis_menten_calibration)")
from mode_prep_raw_data.prep import michaelis_menten_calibration
from app_ui.data_load_mode import (
    _build_v0_fig as _dl_build_v0_fig,
    _build_v0_supplementary_fig as _dl_build_v0_supplementary_fig,
    _apply_display_layout_like_export as _dl_apply_display_layout_like_export,
    _export_fig_to_png_bytes as _dl_export_fig_to_png_bytes,
    _render_client_side_png_download as _dl_render_client_side_png_download,
    _render_client_side_download_all as _dl_render_client_side_download_all,
)
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
    
    st.sidebar.markdown("---")
    st.sidebar.subheader("α (Alpha) normalization")
    use_shared_Finf = st.sidebar.checkbox(
        "Use shared F_∞ (same for all [E])",
        value=True,
        help="Use one F_∞ (F_∞,max = highest plateau) for all [E]. When plateau height depends on [E] (incomplete cleavage in time window), this global normalization α(t)=(F_t−F₀)/(F_∞,max−F₀) is the usual paper approach and addresses reviewer concerns."
    )
    
    # Step 1: Load Fitted Curves data (원본 데이터 플롯용)
    df_fitted = None
    rfu_col = None
    
    # 0순위: Session State 확인 (Data Load 모드에서 방금 실행된 경우)
    if 'interpolation_results' in st.session_state and st.session_state.get('mm_data_ready', False):
        results = st.session_state['interpolation_results']
        # Results Applied 메시지는 항상 표시 (예외 시에도 유지)
        experiment_type = results.get('experiment_type', 'Substrate Concentration Variation (Standard Michaelis-Menten)')
        result_type = "Substrate-based" if experiment_type == "Substrate Concentration Variation (Standard Michaelis-Menten)" else "Enzyme-based"
        raw_display = results.get('raw_data_source_display') or results.get('uploaded_filename') or ''
        if raw_display:
            st.success(f"Results Applied ({result_type}) — Raw data: **{raw_display}**")
        else:
            st.success(f"Results Applied ({result_type})")
        try:
            df_fitted = results['interp_df'].copy()
            rfu_col = 'RFU_Interpolated' if 'RFU_Interpolated' in df_fitted.columns else 'RFU_Calculated'
            st.session_state['data_source_type'] = 'Data Load result (in memory)'
            # Sidebar 값이 Data Load 결과를 따르도록 동기화 (파일명 파싱으로 효소/기질 정보 반영)
            uploaded_filename = results.get('uploaded_filename')
            if uploaded_filename:
                ename, mw, sname = _parse_filename_for_model_simulation(uploaded_filename)
                if ename is not None:
                    st.session_state['parsed_enzyme_name'] = ename
                    st.session_state['model_sim_enzyme_name'] = ename
                if mw is not None:
                    st.session_state['parsed_enzyme_mw'] = mw
                    st.session_state['model_sim_enzyme_mw'] = mw
                if sname is not None:
                    st.session_state['parsed_substrate_name'] = sname
                    st.session_state['model_sim_substrate_name'] = sname
                st.session_state['last_parsed_filename_for_model'] = uploaded_filename
        except Exception:
            # 메모리 로드 실패 시 파일 로드 시도 (아래 uploaded_file/파일 경로 분기로 진행)
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
    df_current = normalizer.normalize_final(df_current, fitted_params=fitted_params, use_shared_Finf=use_shared_Finf)
    
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
    if experiment_type is None and 'sample_mm_fit' in st.session_state:
        experiment_type = st.session_state['sample_mm_fit'].get('experiment_type')

    def _is_substrate_experiment(exp_type):
        return exp_type in (
            "Substrate 농도 변화 (표준 MM)",
            "Substrate Concentration Variation (Standard MM)",
            "Substrate Concentration Variation (Standard Michaelis-Menten)",
        )

    original_conc_col = df['conc_col_name'].iloc[0] if 'conc_col_name' in df.columns else 'Concentration [ug/mL]'
    if _is_substrate_experiment(experiment_type):
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
    # Default to enzyme-based label when experiment type is not yet available.
    v0_fit_tab_title = "📊 v₀ vs [S] Fit" if _is_substrate_experiment(experiment_type) else "📊 v₀ vs [E] Linear Fit"
    tab1, tab_v0_window, tab_kobs, tab_alpha, tab_evsalpha, tab_export = st.tabs([
        v0_fit_tab_title,
        "🕐 v₀ vs [E] by window",
        "📉 k_obs vs [E]",
        "📈 Alpha Calculation",
        "📊 [E] vs α Plot",
        "📤 Export Tables & Plots",
    ])

    # Export tab payloads (filled by each tab below)
    fig_v0_export = None
    fig_v0_supp_export = None
    df_v0_table_export = None
    fig_alpha_export = None
    df_alpha_stats_export = None
    df_alpha_download_export = None
    fig_ev_export = None
    df_ev_export = None
    
    with tab1:
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
            exp_type = mm_fit.get('experiment_type', 'Enzyme Concentration Variation (Fixed substrate)')
            xaxis_title = '[S] (μM)' if _is_substrate_experiment(exp_type) else '[E] (μg/mL)'

            if _is_substrate_experiment(exp_type):
                st.subheader("v₀ vs [S] Michaelis-Menten Fit")
            else:
                st.subheader("v₀ vs [E] Linear Fit (Constant Substrate)")
                st.warning("⚠️ This is not the standard Michaelis-Menten model. v is linearly related to [E].")

            dl_like_results = {
                'v0_vs_concentration': v0_data,
                'mm_fit_results': mm_fit,
            }
            fig_v0 = _dl_build_v0_fig(dl_like_results)
            fig_v0 = _dl_apply_display_layout_like_export(fig_v0)
            st.plotly_chart(fig_v0, use_container_width=True)
            fig_v0_export = go.Figure(fig_v0)

            if not _is_substrate_experiment(exp_type):
                fig_v0_supp = _dl_build_v0_supplementary_fig(dl_like_results)
                if fig_v0_supp is not None:
                    st.info("Overall linear fit R² is below 0.8. Showing supplementary low-3 linear fit and (if 5 points exist) a constant fit for 3rd-5th [E].")
                    fig_v0_supp = _dl_apply_display_layout_like_export(fig_v0_supp)
                    st.plotly_chart(fig_v0_supp, use_container_width=True)
                    fig_v0_supp_export = go.Figure(fig_v0_supp)
            
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
            df_v0_table_export = df_table.copy()
                 
        else:
            st.info("⚠️ No Michaelis-Menten fitting data. Run analysis in Data Load Mode or load a result file that includes the 'Michaelis-Menten Fit Results' sheet.")
    
    with tab_v0_window:
        st.subheader("🕐 v₀ vs [E] linearity by initial-velocity window")
        st.caption("Compare v₀ vs [E] linear fit when v₀ is computed over different time windows (0–30 s, 0–1 min, 0–2 min, …). Useful to check substrate depletion: short windows → better linearity; long windows → high [E] points drop below line.")
        with st.expander("📋 Kinetic 해석 가이드 (정규화 속도 · semi-log · k_obs vs [E] · 고농도 local slope)"):
            st.markdown("""
**A. 각 농도별 정규화된 속도 비교**  
$v(t)/v(0)$ 또는 최대값으로 나눠서 겹쳐 보세요.
- **고농도만** 더 빨리 감소하면 → **depletion** 쪽.
- 모든 조건이 비슷한 형태면 → 단순한 1차형 progress curve 가능성.

**B. Semi-log plot**  
$\\ln v(t)$ vs $t$
- 1차형이면 거의 **직선**입니다.
- 고농도에서만 더 빨리 꺾이면 → depletion/artifact를 의심할 수 있습니다.

**C. 가장 좋은 방법**  
v(t) 직접 미분보다 **progress curve fitting**이 더 낫습니다.  
$F(t) = F_\\infty (1 - e^{-k_{\\mathrm{obs}} t})$  
각 농도에서 $k_{\\mathrm{obs}}$를 구하고, **$k_{\\mathrm{obs}}$ vs [E]**를 보세요.  
- 이게 **linear**면, 적어도 현재 범위에서는 **pseudo-first-order kinetics**로 설명이 잘 됩니다.

**고농도만 따로 확인하는 실무적인 방법**  
고농도 2–3개 조건만 골라서, **0–15 sec**, **15–30 sec**, **30–60 sec**처럼 아주 짧은 **local slope**를 계산해 보세요.
- 고농도에서만 0–15 sec slope가 비정상적으로 크고 바로 감소 → **depletion/mixing artifact** 후보  
- 그냥 매끈한 exponential 감소 → 일반적인 fast kinetics

**지금 단계에서 결론**  
현재 그림만 보면, "고농도에서 substrate depletion 때문에 $v_0$가 과소추정됐다"는 주장은 아직 약합니다.  
오히려 먼저 확인할 건:
1. 수치미분 방식이 **초기 spike를 과장**했는지  
2. **$k_{\\mathrm{obs}}$ vs [E]**가 linear인지  
3. 정규화된 **$v(t)$ shape**가 농도별로 다른지  
입니다.
""")
        exp_type_for_window = (mm_fit or {}).get('experiment_type', '')
        is_enzyme_var = exp_type_for_window in ("Enzyme Concentration Variation (Fixed substrate)", "Enzyme Concentration Variation", "Enzyme 농도 변화 (기질 고정)")
        if not is_enzyme_var:
            st.info("ℹ️ This tab is for **Enzyme Concentration Variation** (fixed substrate, v₀ vs [E]). For substrate variation (MM) experiments this comparison is not applicable.")
        else:
            conc_col_raw = 'enzyme_ugml' if 'enzyme_ugml' in df_raw.columns else next((c for c in ['Concentration [ug/mL]', 'Concentration'] if c in df_raw.columns), None)
            if conc_col_raw is None or 'time_min' not in df_raw.columns:
                st.warning("Time–fluorescence curves not found. Upload a file that contains **Time–FLU Interpolated curves** (or run Data Load Mode) to use this comparison.")
            else:
                time_max_avail = float(df_raw['time_min'].max())
                # Windows: 30 s, 1 min, 2 min, 3 min (in minutes)
                window_minutes = [0.5, 1.0, 2.0, 3.0]
                window_minutes = [w for w in window_minutes if w <= time_max_avail + 0.01]
                if not window_minutes:
                    window_minutes = [time_max_avail] if time_max_avail > 0 else [1.0]
                window_labels = [f"0–{int(w*60)} s" if w < 1 else f"0–{w:.0f} min" for w in window_minutes]
                results_by_window = []
                flu_col = 'FL_intensity' if 'FL_intensity' in df_raw.columns else 'RFU_Interpolated'
                if flu_col not in df_raw.columns and 'RFU_Interpolated' in df_raw.columns:
                    flu_col = 'RFU_Interpolated'

                def _central_diff_velocity(t, F):
                    """v(t_i) = (F_{i+1} - F_{i-1}) / (t_{i+1} - t_{i-1}). Interior points only."""
                    n = len(t)
                    if n < 3:
                        return t, np.full_like(t, np.nan)
                    v = np.full(n, np.nan)
                    for i in range(1, n - 1):
                        dt = t[i + 1] - t[i - 1]
                        if abs(dt) > 1e-12:
                            v[i] = (F[i + 1] - F[i - 1]) / dt
                    return t, v

                for t_max_min, wlabel in zip(window_minutes, window_labels):
                    v0_per_conc = []
                    concs_used = []
                    for conc in df_raw[conc_col_raw].unique():
                        sub = df_raw[(df_raw[conc_col_raw] == conc) & (df_raw['time_min'] <= t_max_min)].sort_values('time_min')
                        if len(sub) < 2:
                            continue
                        t = sub['time_min'].values.astype(float)
                        y = sub[flu_col].values.astype(float)
                        if np.ptp(t) < 1e-6:
                            continue
                        coeffs = np.polyfit(t, y, 1)
                        v0_per_conc.append(coeffs[0])
                        concs_used.append(float(conc))
                    if len(concs_used) < 2:
                        results_by_window.append({'window': wlabel, 't_max_min': t_max_min, 'R2': None, 'slope': None, 'intercept': None, 'concs': concs_used, 'v0': v0_per_conc})
                        continue
                    x = np.array(concs_used)
                    y = np.array(v0_per_conc)
                    coeffs = np.polyfit(x, y, 1)
                    fit = np.polyval(coeffs, x)
                    ss_res = np.sum((y - fit) ** 2)
                    ss_tot = np.sum((y - np.mean(y)) ** 2)
                    r2 = float(1 - ss_res / ss_tot) if ss_tot > 0 else 0
                    results_by_window.append({'window': wlabel, 't_max_min': t_max_min, 'R2': r2, 'slope': coeffs[0], 'intercept': coeffs[1], 'concs': concs_used, 'v0': v0_per_conc})
                table_rows = [{'Window': r['window'], 'R²': f"{r['R2']:.4f}" if r['R2'] is not None else "—", 'Slope': f"{r['slope']:.4f}" if r['slope'] is not None else "—", 'Intercept': f"{r['intercept']:.2f}" if r['intercept'] is not None else "—"} for r in results_by_window]
                st.dataframe(pd.DataFrame(table_rows), use_container_width=True, hide_index=True)
                fig_w = go.Figure()
                colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
                for i, r in enumerate(results_by_window):
                    if r['concs'] and r['v0'] and r['slope'] is not None:
                        fig_w.add_trace(go.Scatter(x=r['concs'], y=r['v0'], mode='markers', name=f"{r['window']} (R²={r['R2']:.4f})", marker=dict(size=10, color=colors[i % len(colors)])))
                        x_smooth = np.linspace(min(r['concs']), max(r['concs']), 50)
                        fig_w.add_trace(go.Scatter(x=x_smooth.tolist(), y=(r['slope'] * x_smooth + r['intercept']).tolist(), mode='lines', name=f"{r['window']} fit", line=dict(dash='dash', color=colors[i % len(colors)])))
                fig_w.update_layout(title="v₀ vs [E] by initial-velocity window", xaxis_title="[E] (μg/mL)", yaxis_title="Initial Velocity v₀ (RFU/min)", template="plotly_white", height=500, hovermode="x unified")
                st.plotly_chart(fig_w, use_container_width=True)
                st.markdown("**해석**: 짧은 window에서 선형성(R²)이 좋고, 긴 window로 갈수록 고농도 [E]가 아래로 처지면 → **substrate depletion** 또는 초기 burst 이후 빠른 기질 소모 가능성이 큽니다.")
                st.markdown("---")

                concs_sorted = sorted(df_raw[conc_col_raw].dropna().unique(), key=lambda x: float(x))
                colors_abc = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f']

                # --- A. 정규화된 속도 비교 v(t)/v(0) 또는 v(t)/v_max ---
                st.subheader("A. 각 농도별 정규화된 속도 비교 v(t)/v(0)")
                st.caption("v(t)/v(0) 또는 최대값으로 나눠서 겹쳐 보세요. 고농도만 더 빨리 감소하면 depletion; 모든 조건이 비슷한 형태면 단순 1차형 progress curve 가능성.")
                norm_by = st.radio("정규화 기준", ["v(0) (첫 유효 속도)", "최대값 v_max"], key="norm_velocity_radio", horizontal=True)
                fig_norm = go.Figure()
                for i, conc in enumerate(concs_sorted):
                    sub = df_raw[df_raw[conc_col_raw] == conc].sort_values('time_min')
                    if len(sub) < 3:
                        continue
                    t = sub['time_min'].values.astype(float)
                    F = sub[flu_col].values.astype(float).astype(np.float64)
                    t_plot, v_t = _central_diff_velocity(t, F)
                    valid = ~np.isnan(v_t) & (v_t > 0)
                    if not np.any(valid):
                        continue
                    v_valid = v_t[valid]
                    t_valid = t_plot[valid]
                    if norm_by.startswith("v(0)"):
                        v0_ref = v_valid[0]
                        if v0_ref <= 0:
                            continue
                        y_plot = v_valid / v0_ref
                    else:
                        v_max = np.nanmax(v_valid)
                        if v_max <= 0:
                            continue
                        y_plot = v_valid / v_max
                    fig_norm.add_trace(go.Scatter(
                        x=t_valid, y=y_plot, mode='lines+markers', name=f"[E] = {float(conc):.2f} μg/mL",
                        line=dict(width=2, color=colors_abc[i % len(colors_abc)]), marker=dict(size=4, color=colors_abc[i % len(colors_abc)])
                    ))
                fig_norm.update_layout(
                    title="Normalized velocity v(t)/v(0) or v(t)/v_max vs Time",
                    xaxis_title="Time (min)", yaxis_title="v(t) / v(ref)",
                    template="plotly_white", height=400, hovermode="x unified"
                )
                st.plotly_chart(fig_norm, use_container_width=True)

                # --- B. Semi-log ln v(t) vs t ---
                st.subheader("B. Semi-log: ln v(t) vs t")
                st.caption("1차형이면 거의 직선. 고농도에서만 더 빨리 꺾이면 depletion/artifact 의심.")
                v_min_cut = st.number_input("v(t) 최소값 (이하 제외, ln 안정성)", value=1e-6, format="%.0e", key="v_min_semilog")
                fig_log = go.Figure()
                for i, conc in enumerate(concs_sorted):
                    sub = df_raw[df_raw[conc_col_raw] == conc].sort_values('time_min')
                    if len(sub) < 3:
                        continue
                    t = sub['time_min'].values.astype(float)
                    F = sub[flu_col].values.astype(float).astype(np.float64)
                    t_plot, v_t = _central_diff_velocity(t, F)
                    valid = ~np.isnan(v_t) & (v_t > v_min_cut)
                    if not np.any(valid):
                        continue
                    fig_log.add_trace(go.Scatter(
                        x=t_plot[valid], y=np.log(v_t[valid]), mode='lines+markers', name=f"[E] = {float(conc):.2f} μg/mL",
                        line=dict(width=2, color=colors_abc[i % len(colors_abc)]), marker=dict(size=4, color=colors_abc[i % len(colors_abc)])
                    ))
                fig_log.update_layout(
                    title="ln v(t) vs Time (by [E])",
                    xaxis_title="Time (min)", yaxis_title="ln v(t)",
                    template="plotly_white", height=400, hovermode="x unified"
                )
                st.plotly_chart(fig_log, use_container_width=True)

                # --- C. k_obs vs [E] (progress curve fitting) — reviewer-facing ---
                st.subheader("C. k_obs vs [E] (progress curve fitting)")
                st.caption("F(t) = F₀ + (F∞ − F₀)(1 − e^(−k_obs·t)), time in seconds → k_obs (s⁻¹). Linear k_obs vs [E] supports pseudo-first-order and addresses substrate depletion / diffusion limitation concerns.")
                with st.expander("📌 Plateau가 [E]에 따라 다를 때 (리뷰어 대응)"):
                    st.markdown("""
                    - **현상**: 최종 형광(plateau)이 효소 농도가 높을수록 커짐.
                    - **해석**: 실험 시간 내 완전 절단 미도달이 가장 흔함 (낮은 [E]에서 F_∞,obs < F_∞,true). 또는 hydrogel/confined 기질에서 접근성 차이.
                    - **대응**: (1) **Shared F_∞** (사이드바)로 α = (F_t−F₀)/(F_∞,max−F₀) 사용. (2) **Kinetics는 k_obs vs [E]** 사용 — v₀ vs [E]보다 안정적.
                    - **논문 문구 예**: *"The final fluorescence intensity increased with enzyme concentration, suggesting that complete substrate cleavage was not reached within the experimental time window at lower enzyme levels."*
                    """)
                use_three_param = st.checkbox("Use 3-parameter fit F₀ + (F∞ − F₀)(1 − e^(−k_obs·t))", value=True, key="kobs_three_param")
                try:
                    progress_results, linear_result = fit_progress_curves_pseudo_first_order(
                        df_raw,
                        conc_col=conc_col_raw,
                        time_col='time_min',
                        fluor_col=flu_col,
                        enzyme_mw_kda=float(enzyme_mw),
                        use_three_param=use_three_param,
                    )
                except Exception as e:
                    progress_results = []
                    linear_result = None
                    st.warning(f"Progress curve fitting failed: {e}")
                if len(progress_results) >= 2 and linear_result is not None:
                    concs_display = [r.conc for r in progress_results]
                    k_obs_s = [r.k_obs_per_s for r in progress_results]
                    fig_k = go.Figure()
                    fig_k.add_trace(go.Scatter(x=concs_display, y=k_obs_s, mode='markers', name="k_obs (s⁻¹)", marker=dict(size=12, color='#1f77b4')))
                    E_M = linear_result.E_M
                    k_arr = linear_result.k_obs
                    slope = linear_result.slope_M_inv_s
                    intercept = linear_result.intercept_per_s
                    x_line_ugml = np.linspace(min(concs_display), max(concs_display), 50)
                    E_M_line = (x_line_ugml / enzyme_mw) * 1e-6
                    y_line = slope * E_M_line + intercept
                    fig_k.add_trace(go.Scatter(x=x_line_ugml.tolist(), y=y_line.tolist(), mode='lines', name="Linear fit k_obs vs [E]", line=dict(dash='dash', color='#ff7f0e')))
                    fig_k.update_layout(
                        title="k_obs vs [E] (pseudo-first-order)",
                        xaxis_title="[E] (μg/mL)",
                        yaxis_title="k_obs (s⁻¹)",
                        template="plotly_white",
                        height=400,
                        hovermode="x unified",
                    )
                    st.plotly_chart(fig_k, use_container_width=True)
                    st.markdown(f"**선형성**: R² = {linear_result.r_squared:.4f} — linear에 가까우면 현재 범위에서 **pseudo-first-order kinetics**로 설명 가능.")
                    st.markdown(f"**k_cat/K_M** (slope of k_obs vs [E] in M): **{slope:.3e}** M⁻¹ s⁻¹ (± {linear_result.slope_std:.3e}).")
                    st.markdown("**For reviewers**: *The reaction progress curves were fitted to F(t) = F₀ + (F∞ − F₀)(1 − e^(−k_obs·t)). The observed rate constant k_obs was linear with enzyme concentration (R² = {:.4f}), consistent with pseudo-first-order kinetics and arguing against substrate depletion or diffusion limitation dominating in the time window used.*".format(linear_result.r_squared))
                    table_df = pd.DataFrame({
                        "[E] (μg/mL)": [r.conc for r in progress_results],
                        "k_obs (s⁻¹)": [f"{r.k_obs_per_s:.6f}" for r in progress_results],
                        "k_obs std": [f"{r.k_obs_std:.4e}" for r in progress_results],
                        "F₀": [f"{r.F0:.2f}" for r in progress_results],
                        "F∞": [f"{r.F_inf:.2f}" for r in progress_results],
                        "R²": [f"{r.r_squared:.4f}" for r in progress_results],
                    })
                    st.dataframe(table_df, use_container_width=True, hide_index=True)
                elif len(progress_results) == 1:
                    st.warning("k_obs 피팅에 성공한 농도가 1개뿐이라 k_obs vs [E] 선형 회귀를 할 수 없습니다. 데이터 구간 또는 농도 수를 확인하세요.")
                    table_df = pd.DataFrame({
                        "[E] (μg/mL)": [r.conc for r in progress_results],
                        "k_obs (s⁻¹)": [f"{r.k_obs_per_s:.6f}" for r in progress_results],
                        "F₀": [f"{r.F0:.2f}" for r in progress_results],
                        "F∞": [f"{r.F_inf:.2f}" for r in progress_results],
                        "R²": [f"{r.r_squared:.4f}" for r in progress_results],
                    })
                    st.dataframe(table_df, use_container_width=True, hide_index=True)
                else:
                    st.warning("k_obs 피팅에 성공한 농도가 없습니다. 데이터 구간 또는 농도 수를 확인하세요.")

                st.markdown("---")
                st.subheader("진단 1: Progress curve의 국소 기울기 v(t) = dF/dt")
                st.caption("1) Progress curve smoothing (optional) → 2) Numerical derivative v(t_i) = (F_{i+1} − F_{i-1}) / (t_{i+1} − t_{i-1}) → 3) v(t) vs time. 특히 고농도 [E] 조건을 보세요.")
                use_smooth = st.checkbox("Progress curve smoothing (optional)", value=False, key="v0_window_smooth_curve")
                smooth_window = 5
                if use_smooth:
                    smooth_window = st.slider("Smoothing window (points)", min_value=3, max_value=15, value=5, step=2, key="v0_window_smooth_win")

                colors_vt = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f']
                fig_vt = go.Figure()
                for i, conc in enumerate(concs_sorted):
                    sub = df_raw[df_raw[conc_col_raw] == conc].sort_values('time_min')
                    if len(sub) < 3:
                        continue
                    t = sub['time_min'].values.astype(float)
                    F = sub[flu_col].values.astype(float).astype(np.float64)
                    if use_smooth and len(F) >= smooth_window:
                        win = min(smooth_window if (smooth_window % 2 == 1) else smooth_window - 1, len(F))
                        if win % 2 == 0:
                            win = max(3, win - 1)
                        try:
                            from scipy.signal import savgol_filter
                            F = savgol_filter(F, win, 2)
                        except Exception:
                            k = np.ones(min(smooth_window, len(F))) / min(smooth_window, len(F))
                            F = np.convolve(F, k, mode='same')
                    t_plot, v_t = _central_diff_velocity(t, F)
                    valid = ~np.isnan(v_t)
                    if not np.any(valid):
                        continue
                    fig_vt.add_trace(go.Scatter(
                        x=t_plot[valid],
                        y=v_t[valid],
                        mode='lines+markers',
                        name=f"[E] = {float(conc):.2f} μg/mL",
                        line=dict(width=2, color=colors_vt[i % len(colors_vt)]),
                        marker=dict(size=4, color=colors_vt[i % len(colors_vt)]),
                        legendgroup=str(conc)
                    ))
                fig_vt.update_layout(
                    title="v(t) = dF/dt vs Time (by [E])",
                    xaxis_title="Time (min)",
                    yaxis_title="v(t) (dF/dt, RFU/min)",
                    template="plotly_white",
                    height=450,
                    hovermode="x unified",
                    showlegend=True
                )
                st.plotly_chart(fig_vt, use_container_width=True)
                st.markdown("**해석 포인트 (특히 high [E] 조건)**")
                st.markdown("""
                - **Depletion이면**: v(t) **초반 높음** → v(t) **급격히 감소** (기질 소모로 속도 하락).
                - **Mixing artifact이면**: v(t) **초반 낮음** → v(t) **상승 후 plateau** (혼합 지연 후 안정).
                """)
    
    with tab_kobs:
        st.subheader("📉 k_obs vs [E] (exponential fit → pseudo-first-order)")
        st.caption("기질이 일정하고 pseudo-first-order에 가깝다면, 각 progress curve를 지수 피팅하여 k_obs를 추출한 뒤 **k_obs가 [E]에 비례**하는지 보는 것이 v₀ vs [E]보다 더 자연스럽고 설득력 있습니다.")
        exp_type_kobs = (mm_fit or {}).get('experiment_type', '')
        is_enzyme_var_kobs = exp_type_kobs in ("Enzyme Concentration Variation (Fixed substrate)", "Enzyme Concentration Variation", "Enzyme 농도 변화 (기질 고정)")
        if not is_enzyme_var_kobs:
            st.info("ℹ️ 이 탭은 **효소 농도 변화** 실험(기질 고정, v₀ vs [E])용입니다. 기질 농도 변화(MM) 실험에는 해당하지 않습니다.")
        else:
            r2_full = (mm_fit or {}).get('R_squared')
            conc_vals = np.array((v0_data or {}).get('concentrations', []), dtype=float)
            v0_vals = np.array((v0_data or {}).get('v0_values', []), dtype=float)
            is_supplementary_case = (
                r2_full is not None
                and r2_full < 0.8
                and len(conc_vals) >= 3
                and len(v0_vals) >= 3
                and len(conc_vals) == len(v0_vals)
            )

            if not is_supplementary_case:
                st.info("ℹ️ 이 탭은 v₀ vs [E] 전체 선형 피팅의 R²가 0.8 미만일 때(= supplementary linear fit이 필요한 경우)에만 진행합니다.")
            elif len(conc_vals) < 5:
                st.info("ℹ️ 현재는 low-3 supplementary fit만 가능하고 constant fit(정렬 3–5번째 [E])이 구성되지 않습니다. k_obs vs [E]는 constant-fit [E] 구간이 있을 때만 진행합니다.")
            else:
                sort_idx = np.argsort(conc_vals)
                conc_sorted = conc_vals[sort_idx]
                const_fit_concs = conc_sorted[2:5]
                st.caption(f"Using constant-fit [E] only: {', '.join([f'{c:.4g}' for c in const_fit_concs])} μg/mL")

                conc_col_kobs = 'enzyme_ugml' if 'enzyme_ugml' in df_raw.columns else next((c for c in ['Concentration [ug/mL]', 'Concentration'] if c in df_raw.columns), None)
                if conc_col_kobs is None or 'time_min' not in df_raw.columns:
                    st.warning("Time–fluorescence 곡선이 없습니다. **Time–FLU Interpolated curves**가 포함된 파일을 업로드하거나 Data Load Mode를 실행하세요.")
                else:
                    flu_col_kobs = 'FL_intensity' if 'FL_intensity' in df_raw.columns else 'RFU_Interpolated'
                    if flu_col_kobs not in df_raw.columns and 'RFU_Interpolated' in df_raw.columns:
                        flu_col_kobs = 'RFU_Interpolated'

                    conc_numeric_all = pd.to_numeric(
                        df_raw[conc_col_kobs].astype(str).str.extract(r'([-+]?\d*\.?\d+)')[0],
                        errors='coerce'
                    )
                    tol_abs, tol_rel = 1e-6, 1e-4
                    mask_const = np.zeros(len(df_raw), dtype=bool)
                    conc_array = conc_numeric_all.to_numpy(dtype=float)
                    valid_mask = ~np.isnan(conc_array)
                    for c in const_fit_concs:
                        mask_const |= valid_mask & np.isclose(conc_array, float(c), atol=tol_abs, rtol=tol_rel)
                    df_raw_kobs = df_raw.loc[mask_const].copy()

                    if df_raw_kobs.empty:
                        st.warning("constant-fit [E] 구간(정렬 3–5번째)에 해당하는 time–fluorescence 데이터가 없어 k_obs 분석을 진행할 수 없습니다.")
                    else:
                        use_three_param_kobs = st.checkbox("3-parameter fit: F₀ + (F∞ − F₀)(1 − e^(−k_obs·t))", value=True, key="kobs_tab_three_param")
                        try:
                            progress_results_kobs, linear_result_kobs = fit_progress_curves_pseudo_first_order(
                                df_raw_kobs,
                                conc_col=conc_col_kobs,
                                time_col='time_min',
                                fluor_col=flu_col_kobs,
                                enzyme_mw_kda=float(enzyme_mw),
                                use_three_param=use_three_param_kobs,
                            )
                        except Exception as e:
                            progress_results_kobs = []
                            linear_result_kobs = None
                            st.warning(f"Progress curve 피팅 실패: {e}")

                        if len(progress_results_kobs) >= 2 and linear_result_kobs is not None:
                            # 1) Progress curves with exponential fit overlay (optional expander)
                            with st.expander("📈 각 농도별 progress curve + exponential fit (constant-fit [E] only)", expanded=False):
                                colors_curves = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f']
                                fig_curves = go.Figure()
                                conc_numeric_subset = pd.to_numeric(
                                    df_raw_kobs[conc_col_kobs].astype(str).str.extract(r'([-+]?\d*\.?\d+)')[0],
                                    errors='coerce'
                                ).to_numpy(dtype=float)
                                for i, r in enumerate(progress_results_kobs):
                                    r_conc = float(r.conc)
                                    sub_mask = (~np.isnan(conc_numeric_subset)) & np.isclose(conc_numeric_subset, r_conc, atol=tol_abs, rtol=tol_rel)
                                    sub = df_raw_kobs.loc[sub_mask].sort_values('time_min')
                                    if len(sub) < 2:
                                        continue
                                    t_min = sub['time_min'].values.astype(float)
                                    F = sub[flu_col_kobs].values.astype(float)
                                    c = colors_curves[i % len(colors_curves)]
                                    fig_curves.add_trace(go.Scatter(x=t_min, y=F, mode='markers', name=f"[E] = {r_conc:.2f} μg/mL (data)", marker=dict(size=4, color=c)))
                                    t_s = t_min * 60.0
                                    F_fit = r.F0 + (r.F_inf - r.F0) * (1.0 - np.exp(-r.k_obs_per_s * t_s))
                                    fig_curves.add_trace(go.Scatter(x=t_min, y=F_fit, mode='lines', name=f"[E] = {r_conc:.2f} fit (k_obs={r.k_obs_per_s:.4f} s⁻¹)", line=dict(dash='dash', color=c)))
                                fig_curves.update_layout(
                                    title="Progress curves with exponential fit F(t) = F₀ + (F∞ − F₀)(1 − e^(−k_obs·t))",
                                    xaxis_title="Time (min)", yaxis_title="Fluorescence",
                                    template="plotly_white", height=450, hovermode="x unified",
                                )
                                st.plotly_chart(fig_curves, use_container_width=True)

                            # 2) k_obs vs [E] plot + linear fit
                            concs_display_kobs = [r.conc for r in progress_results_kobs]
                            k_obs_s_kobs = [r.k_obs_per_s for r in progress_results_kobs]
                            fig_kobs = go.Figure()
                            fig_kobs.add_trace(go.Scatter(x=concs_display_kobs, y=k_obs_s_kobs, mode='markers', name="k_obs (s⁻¹)", marker=dict(size=12, color='#1f77b4')))
                            slope_kobs = linear_result_kobs.slope_M_inv_s
                            intercept_kobs = linear_result_kobs.intercept_per_s
                            x_line_ugml = np.linspace(min(concs_display_kobs), max(concs_display_kobs), 50)
                            E_M_line = (x_line_ugml / enzyme_mw) * 1e-6
                            y_line = slope_kobs * E_M_line + intercept_kobs
                            fig_kobs.add_trace(go.Scatter(x=x_line_ugml.tolist(), y=y_line.tolist(), mode='lines', name="Linear fit k_obs vs [E]", line=dict(dash='dash', color='#ff7f0e')))
                            fig_kobs.update_layout(
                                title="k_obs vs [E] (constant-fit [E] subset)",
                                xaxis_title="[E] (μg/mL)",
                                yaxis_title="k_obs (s⁻¹)",
                                template="plotly_white",
                                height=420,
                                hovermode="x unified",
                            )
                            st.plotly_chart(fig_kobs, use_container_width=True)
                            st.markdown(f"**선형성**: R² = {linear_result_kobs.r_squared:.4f} — constant-fit [E] 구간에서의 선형성입니다.")
                            st.markdown(f"**k_cat/K_M** (slope of k_obs vs [E] in M): **{slope_kobs:.3e}** M⁻¹ s⁻¹ (± {linear_result_kobs.slope_std:.3e}).")
                            table_kobs = pd.DataFrame({
                                "[E] (μg/mL)": [r.conc for r in progress_results_kobs],
                                "k_obs (s⁻¹)": [f"{r.k_obs_per_s:.6f}" for r in progress_results_kobs],
                                "k_obs std": [f"{r.k_obs_std:.4e}" for r in progress_results_kobs],
                                "F₀": [f"{r.F0:.2f}" for r in progress_results_kobs],
                                "F∞": [f"{r.F_inf:.2f}" for r in progress_results_kobs],
                                "R²": [f"{r.r_squared:.4f}" for r in progress_results_kobs],
                            })
                            st.dataframe(table_kobs, use_container_width=True, hide_index=True)
                        elif len(progress_results_kobs) == 1:
                            st.warning("k_obs 피팅에 성공한 농도가 1개뿐이라 k_obs vs [E] 선형 회귀를 할 수 없습니다.")
                            st.dataframe(pd.DataFrame({
                                "[E] (μg/mL)": [r.conc for r in progress_results_kobs],
                                "k_obs (s⁻¹)": [f"{r.k_obs_per_s:.6f}" for r in progress_results_kobs],
                                "F₀": [f"{r.F0:.2f}" for r in progress_results_kobs],
                                "F∞": [f"{r.F_inf:.2f}" for r in progress_results_kobs],
                                "R²": [f"{r.r_squared:.4f}" for r in progress_results_kobs],
                            }), use_container_width=True, hide_index=True)
                        else:
                            st.warning("k_obs 피팅에 성공한 농도가 없습니다. 데이터 구간 또는 농도 수를 확인하세요.")

    with tab_alpha:
        st.subheader("📈 Alpha (α) Calculation")

        with st.expander("📖 What is Alpha (α)? — definition and formula", expanded=False):
            st.markdown("""
            **What is Alpha (α)?**
            Normalized cleavage ratio: 0 (no cleavage) ~ 1 (full cleavage).

            **Formula**: α(t) = (F_t − F₀) / (F_∞ − F₀)
            - **F_t (F(t))**: Fluorescence at time t  
              From Data Load: interpolated curve (RFU_Interpolated). Otherwise: raw fluorescence.
            - **F₀**: Initial fluorescence (t = 0)  
              Per concentration. From Data Load (MM Results) if available; otherwise min(F) per concentration.
            - **F_∞**: Complete cleavage fluorescence  
              **Often the same value for all enzyme concentrations** (full cleavage → same signal).  
              Sidebar option **"Use shared F_∞ (same for all [E])"** uses one F_∞ for all conditions (typical in papers).  
              If off, F_∞ per concentration is chosen by:
              - **plateau_avg**: mean fluorescence of the plateau (mean(F_plateau))
              - **exponential_fit**: F∞ from exponential fit (F(t) = F₀ + A·(1 − e^(−k·t)), Fmax = F₀ + A)
              - **fallback_max**: maximum fluorescence (max(F)) when no plateau or insufficient rise
            """)

        # Check if alpha column exists
        if 'alpha' not in df.columns:
            st.error("❌ Alpha was not computed. Data normalization is required.")
            st.info("💡 Data is not normalized. Please check the data load and normalization steps.")
        else:
            # 플롯 타겟 데이터 설명 및 미리보기 (비교용)
            conc_col_alpha = 'enzyme_ugml' if 'enzyme_ugml' in df.columns else (df['conc_col_name'].iloc[0] if 'conc_col_name' in df.columns else None)
            with st.expander("📌 플롯 타겟 데이터 (Plot target data — 비교용)", expanded=False):
                st.markdown("""
                **데이터 소스**  
                플롯과 테이블(Alpha Statistics) 모두 **동일한 DataFrame `df`** 사용.
                - 컬럼: `time_s`(시간), `alpha`(절단비율), 농도 컬럼(예: `enzyme_ugml`)
                - `df` 생성 경로: Fitted Curves(보간 곡선) → 단위 표준화 → 정규화(α 계산)
                """)
                st.markdown("**플롯에 실제로 그려지는 값** (분 단위 점 + 꺾은선 모드)")
                st.markdown("""
                - **x**: 분 단위 그리드 `0, 1, 2, …, max(분)`
                - **y**: 각 그리드 시간에서 `(time_s, alpha)`를 선형 보간한 α
                - 즉 원본 `df`의 모든 `(time_s, alpha)`가 아니라 **정수 분 시점에서만** 보간된 α를 점으로 찍고 직선으로 연결.
                """)
                st.markdown("**테이블(Alpha mean)**")
                st.markdown("""
                - 해당 농도에서 **전체** `(time_s, alpha)` 행에 대한 **평균** = `df[농도].alpha.mean()`
                - 플롯은 “분 단위 31개 점”만 쓰고, 테이블은 “모든 시간점”을 쓰므로, 시간점 개수/분포가 다르면 수치가 달라질 수 있음.
                """)
                if conc_col_alpha:
                    # 분 단위 보간값 테이블 (플롯에 그려지는 (분, α) 일부만 표시)
                    plot_target_rows = []
                    for conc in sorted(df[conc_col_alpha].unique()):
                        subset = df[df[conc_col_alpha] == conc].sort_values('time_s')
                        t = subset['time_s'].values
                        alpha = subset['alpha'].values
                        if len(t) == 0:
                            continue
                        t_max = float(np.nanmax(t))
                        if t_max <= 120:
                            minute_grid = np.arange(0, int(np.ceil(t_max)) + 1, 1, dtype=float)
                        else:
                            n_minutes = int(np.ceil(t_max / 60))
                            minute_grid = np.arange(0, n_minutes + 1, 1) * 60.0
                        alpha_at_minutes = np.interp(minute_grid, t, alpha)
                        for i, (tm, a) in enumerate(zip(minute_grid, alpha_at_minutes)):
                            if i % 5 == 0 or i == len(minute_grid) - 1:  # 5분 간격 + 마지막 점만
                                plot_target_rows.append({
                                    f'Concentration ({conc_unit})': conc,
                                    'Time (plot x)': tm,
                                    'α (plot y)': round(a, 4)
                                })
                    if plot_target_rows:
                        st.caption("플롯에 사용되는 (시간, α) 샘플 (5분 간격 + 마지막)")
                        st.dataframe(pd.DataFrame(plot_target_rows), use_container_width=True, hide_index=True)
            
            # Alpha vs Time Plot
            st.subheader("📊 Normalized Data: α(t) vs Time")
            
            fig_alpha = Visualizer.plot_normalized_data(df, conc_unit, time_label, 
                                                       use_lines=True,
                                                       enzyme_name=enzyme_name,
                                                       substrate_name=substrate_name,
                                                       experiment_type=experiment_type,
                                                       points_every_minute=False)
            # 원본 시간 범위로 xaxis 설정
            original_time_max = st.session_state.get('original_time_max', df['time_s'].max())
            if time_unit == 'min':
                fig_alpha.update_xaxes(range=[0, original_time_max])
            else:
                fig_alpha.update_xaxes(range=[0, original_time_max])
            st.plotly_chart(fig_alpha, use_container_width=True)
            # 색상–농도 매핑 명시 (혼동 방지): Plotly Set1 순서 = 농도 오름차순
            conc_order = sorted(df[conc_col_alpha].unique()) if conc_col_alpha else []
            if len(conc_order) <= 5:
                color_names = ['빨강', '파랑', '초록', '보라', '주황'][:len(conc_order)]
                mapping = ', '.join(f'{color_names[i]}={c}' for i, c in enumerate(conc_order))
                st.caption(f"**색상–농도 매핑** (농도 오름차순): {mapping} {conc_unit}. 범례와 선 색상은 동일 기준입니다.")
            fig_alpha_export = go.Figure(fig_alpha)
            
            # Alpha Statistics
            st.subheader("📋 Alpha Statistics by Concentration")
            
            conc_col = 'enzyme_ugml' if 'enzyme_ugml' in df.columns else df['conc_col_name'].iloc[0] if 'conc_col_name' in df.columns else None
            
            if conc_col:
                alpha_stats = []
                # 실험(raw) 시간점 (분 기준): 0, 1, 5, 10, 15, 20, 25, 30
                exp_times_min = np.array([0, 1, 5, 10, 15, 20, 25, 30], dtype=float)
                # 현재 time 단위에 맞게 변환
                current_time_unit = st.session_state.get('time_unit', 'min')
                if current_time_unit == 's':
                    exp_times_target = exp_times_min * 60.0
                else:
                    exp_times_target = exp_times_min
                
                for conc in sorted(df[conc_col].unique()):
                    subset = df[df[conc_col] == conc].copy()
                    alpha_all = subset['alpha'].dropna().values
                    mean_all = float(alpha_all.mean()) if alpha_all.size > 0 else float('nan')
                    std_all = float(alpha_all.std(ddof=0)) if alpha_all.size > 0 else float('nan')
                    
                    # 8개 실험 시간점에서의 α만 사용해 mean/std 계산
                    times = subset['time_s'].values
                    alpha_exp_points = []
                    for t_target in exp_times_target:
                        if times.size == 0:
                            continue
                        idx = int(np.argmin(np.abs(times - t_target)))
                        alpha_exp_points.append(alpha_all[idx] if idx < alpha_all.size else subset['alpha'].iloc[idx])
                    alpha_exp_points = np.array(alpha_exp_points, dtype=float)
                    if alpha_exp_points.size > 0:
                        mean_exp = float(alpha_exp_points.mean())
                        std_exp = float(alpha_exp_points.std(ddof=0))
                        mean_str = f"{mean_all:.4f} ({mean_exp:.4f})"
                        std_str = f"{std_all:.4f} ({std_exp:.4f})"
                    else:
                        mean_str = f"{mean_all:.4f}"
                        std_str = f"{std_all:.4f}"
                    
                    alpha_stats.append({
                        f'Concentration ({conc_unit})': conc,
                        'Alpha mean': mean_str,
                        'Alpha std': std_str
                    })
                
                df_alpha_stats = pd.DataFrame(alpha_stats)
                st.caption("괄호 안 값은 0, 1, 5, 10, 15, 20, 25, 30분 **8개 실험 시간점만 사용했을 때의 mean/std** 입니다.")
                st.dataframe(df_alpha_stats, use_container_width=True, hide_index=True)
                df_alpha_stats_export = df_alpha_stats.copy()
            
            # F0, Fmax 정보
            st.subheader("🔬 Normalization Parameters (F₀, Fmax)")
            
            # Check if fitted parameters are being used
            fitted_params_used = st.session_state.get('fitted_params', None)
            using_fitted_params = fitted_params_used is not None and len(fitted_params_used) > 0
            
            if use_shared_Finf:
                st.info("📌 **Shared F_∞**: One F_∞ (complete cleavage) for all [E]. α(t) = (F_t − F₀) / (F_∞ − F₀).")
            if using_fitted_params:
                st.success(f"✅ F0, Fmax parameters loaded ({len(fitted_params_used)} concentrations)")
                st.info("💡 F0, Fmax are constants from the normalization exponential in Data Load Mode.")
                st.info("📊 Formula used: F(t) = F₀ + (Fmax - F₀)·[1 - exp(-k_obs·t)]")
            else:
                st.info("ℹ️ Using default normalization (Region-based calculation)")
            
            # F0, Fmax 테이블
            if conc_col and 'F0' in df.columns and 'Fmax' in df.columns:
                st.caption("**Fmax method**: how F_∞ was computed for this concentration — **plateau_avg** = mean of plateau region; **exponential_fit** = F∞ from exponential curve; **fallback_max** = max(F) when no plateau/sufficient rise. **Source**: origin of F₀/Fmax (e.g. Data Load or region-based).")
                f0_fmax_data = []
                _method_display = {
                    'plateau_mean': 'plateau_avg',
                    'exponential_Finf': 'exponential_fit',
                    'fallback_max': 'fallback_max',
                    'fallback_temp': 'fallback_max',
                    'fitted_from_data_load': 'Data Load (exponential)',
                    'shared_Finf': 'shared F_∞',
                    'shared_Finf (from data)': 'shared F_∞ (from data)',
                }
                for conc in sorted(df[conc_col].unique()):
                    subset = df[df[conc_col] == conc]
                    fmax_method_raw = subset['Fmax_method'].iloc[0] if 'Fmax_method' in subset.columns else "N/A"
                    fmax_method = _method_display.get(fmax_method_raw, fmax_method_raw)
                    
                    # F0, Fmax 값의 출처 확인
                    df_F0 = subset['F0'].iloc[0]
                    df_Fmax = subset['Fmax'].iloc[0]
                    
                    # fitted_params에서 값 확인
                    source_info = "Region-based (see Fmax method)"
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
                                if fmax_method_raw == 'fitted_from_data_load':
                                    source_info = "Data Load (normalization exponential)"
                                else:
                                    source_info = "Data Load (normalization exponential)"
                            else:
                                source_info = f"Data Load (F0={fitted_F0:.2f}, Fmax={fitted_Fmax:.2f})"
                    
                    f0_fmax_data.append({
                        f'Concentration ({conc_unit})': conc,
                        'F₀ (initial)': f"{df_F0:.2f}",
                        'Fmax (F_∞)': f"{df_Fmax:.2f}",
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
                       - Fmax: **plateau_avg** if plateau exists; else **exponential_fit** if ≥3 points in rise; else **fallback_max**
                       - α = (F - F0) / (Fmax - F0)
                    
                    **Fmax method (how F_∞ was computed):**
                    - **plateau_avg**: mean fluorescence of the plateau (mean(F_plateau))
                    - **exponential_fit**: F∞ from exponential fit (F(t) = F₀ + A·(1 − e^(−k·t)), Fmax = F₀ + A)
                    - **fallback_max**: maximum fluorescence (max(F)) when no plateau or insufficient exponential rise
                    """)
            
            # Download Alpha data
            st.subheader("💾 Download Alpha Data")
            
            # Alpha 데이터 준비
            alpha_download_df = df[['time_s', conc_col, 'alpha', 'F0', 'Fmax']].copy() if conc_col else df[['time_s', 'alpha', 'F0', 'Fmax']].copy()
            alpha_download_df = alpha_download_df.sort_values(['time_s', conc_col] if conc_col else 'time_s')
            df_alpha_download_export = alpha_download_df.copy()
            
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
                fig_ev_export = go.Figure(fig_ev)
                df_ev_export = agg.copy()
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

    with tab_export:
        st.subheader("📤 Export Tables & Plots")
        st.caption("Save analysis tables as one XLSX workbook (multiple sheets) and plots as PNG (individually or bundled ZIP).")

        table_items = []
        if df_v0_table_export is not None and len(df_v0_table_export) > 0:
            table_items.append(("v0_fit_table", df_v0_table_export))
        if df_alpha_stats_export is not None and len(df_alpha_stats_export) > 0:
            table_items.append(("alpha_stats_by_concentration", df_alpha_stats_export))
        if df_alpha_download_export is not None and len(df_alpha_download_export) > 0:
            table_items.append(("alpha_time_series", df_alpha_download_export))
        if df_ev_export is not None and len(df_ev_export) > 0:
            table_items.append(("E_vs_alpha_summary", df_ev_export))

        plot_items = []
        if fig_v0_export is not None:
            plot_items.append(("v0_fit_plot", fig_v0_export))
        if fig_v0_supp_export is not None:
            plot_items.append(("v0_supplementary_plot", fig_v0_supp_export))
        if fig_alpha_export is not None:
            plot_items.append(("alpha_vs_time_plot", fig_alpha_export))
        if fig_ev_export is not None:
            plot_items.append(("E_vs_alpha_plot", fig_ev_export))

        st.markdown(f"- Tables: **{len(table_items)}**")
        st.markdown(f"- Plots: **{len(plot_items)}**")

        if table_items:
            st.markdown("### Table Export (Single XLSX)")
            xlsx_output = BytesIO()
            with pd.ExcelWriter(xlsx_output, engine="openpyxl") as writer:
                for table_name, table_df in table_items:
                    sheet_map = {
                        "v0_fit_table": "v0 fit table",
                        "alpha_stats_by_concentration": "alpha stats",
                        "alpha_time_series": "alpha time series",
                        "E_vs_alpha_summary": "E vs alpha summary",
                    }
                    sheet_name = sheet_map.get(table_name, table_name)[:31]
                    table_df.to_excel(writer, sheet_name=sheet_name, index=False)
            xlsx_bytes = xlsx_output.getvalue()
            st.download_button(
                label=f"📥 Save All Tables (XLSX, {len(table_items)} sheets)",
                data=xlsx_bytes,
                file_name="model_simulation_tables.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                key="general_export_tables_xlsx",
            )
            for table_name, table_df in table_items:
                with st.expander(f"📋 {table_name}", expanded=False):
                    st.dataframe(table_df, use_container_width=True, hide_index=True, height=260)
        else:
            st.info("No table data available for export.")

        if plot_items:
            st.markdown("### Plot Export (PNG)")
            save_all_png_placeholder = st.empty()
            save_all_png_placeholder.caption("Rendering PNG...")
            png_items = []
            failed_plot_count = 0
            ordered_plot_items = []
            for idx, (plot_name, fig) in enumerate(plot_items, start=1):
                fig_for_export = _dl_apply_display_layout_like_export(go.Figure(fig))
                fig_for_export.update_layout(
                    paper_bgcolor="rgba(0,0,0,0)",
                    plot_bgcolor="rgba(0,0,0,0)",
                )
                ordered_plot_name = f"{idx:02d}_{plot_name}"
                ordered_plot_items.append((ordered_plot_name, fig_for_export))
                with st.expander(f"📈 {plot_name}", expanded=False):
                    st.plotly_chart(fig_for_export, use_container_width=True)
                    png_bytes = _dl_export_fig_to_png_bytes(fig_for_export, plot_name=plot_name)
                    if png_bytes:
                        st.download_button(
                            label=f"💾 Save PNG ({plot_name})",
                            data=png_bytes,
                            file_name=f"{ordered_plot_name}.png",
                            mime="image/png",
                            key=f"general_export_png_{plot_name}",
                        )
                        png_items.append((ordered_plot_name, png_bytes))
                    else:
                        failed_plot_count += 1
                        _dl_render_client_side_png_download(fig_for_export, ordered_plot_name, iframe_height=780)

            if png_items:
                if failed_plot_count == 0:
                    zip_buffer = BytesIO()
                    with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zf:
                        for name, png_data in png_items:
                            zf.writestr(f"{name}.png", png_data)
                    zip_buffer.seek(0)
                    with save_all_png_placeholder.container():
                        st.download_button(
                            label=f"📦 Save All PNG ({len(png_items)})",
                            data=zip_buffer.getvalue(),
                            file_name="model_simulation_plots.zip",
                            mime="application/zip",
                            key="general_export_all_png_zip",
                        )
                else:
                    with save_all_png_placeholder.container():
                        _dl_render_client_side_download_all(
                            ordered_plot_items,
                            iframe_height=74,
                            compact=True,
                            button_label=f"📥 Save All PNG ({len(ordered_plot_items)})",
                        )
            else:
                with save_all_png_placeholder.container():
                    _dl_render_client_side_download_all(
                        ordered_plot_items,
                        iframe_height=74,
                        compact=True,
                        button_label=f"📥 Save All PNG ({len(ordered_plot_items)})",
                    )
        else:
            st.info("No plot data available for export.")


