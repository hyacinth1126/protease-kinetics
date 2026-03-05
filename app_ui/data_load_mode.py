import os
import re
import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
from PIL import Image

try:
    import cv2
    CV2_AVAILABLE = True
except ImportError:
    CV2_AVAILABLE = False

try:
    import pytesseract
    from pytesseract import Output
    TESSERACT_AVAILABLE = True
except ImportError:
    TESSERACT_AVAILABLE = False

try:
    import easyocr
    EASYOCR_AVAILABLE = True
except ImportError:
    EASYOCR_AVAILABLE = False

from mode_prep_raw_data.prep import (
    read_raw_data,
    fit_time_course,
    fit_calibration_curve,
    michaelis_menten_calibration,
    calculate_initial_velocity
)
from data_interpolation_mode.interpolate_prism import (
    exponential_association,
    create_prism_interpolation_range
)
from scipy.optimize import curve_fit


def detect_lines_and_points(image_array):
    """
    이미지에서 선과 점을 감지하는 함수
    """
    if not CV2_AVAILABLE:
        return None, None
    
    try:
        # 그레이스케일 변환
        if len(image_array.shape) == 3:
            gray = cv2.cvtColor(image_array, cv2.COLOR_RGB2GRAY)
        else:
            gray = image_array
        
        # 이진화
        _, binary = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)
        
        # 선 감지 (HoughLinesP)
        lines = cv2.HoughLinesP(binary, 1, np.pi/180, threshold=50, minLineLength=30, maxLineGap=10)
        
        # 점 감지 (contour 기반)
        contours, _ = cv2.findContours(binary, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        points = []
        for contour in contours:
            area = cv2.contourArea(contour)
            if 5 < area < 100:  # 점 크기 범위
                M = cv2.moments(contour)
                if M["m00"] != 0:
                    cx = int(M["m10"] / M["m00"])
                    cy = int(M["m01"] / M["m00"])
                    points.append((cx, cy))
        
        return lines, points
    except Exception as e:
        st.warning(f"자동 감지 오류: {e}")
        return None, None


def extract_line_data_from_image(image_file, lines):
    """
    이미지에서 선 데이터를 추출하고 exponential association 모델로 fitting
    """
    try:
        image = Image.open(image_file)
        img_array = np.array(image)
        
        if lines is None or len(lines) == 0:
            return None
        
        # 선에서 데이터 포인트 추출 (간단한 예시)
        # 실제로는 좌표 변환 및 축 스케일 추출이 필요
        st.info("💡 선 데이터 추출: Exponential Association 모델로 fitting합니다.")
        
        # 여기서는 수동 입력으로 대체
        return None
        
    except Exception as e:
        st.error(f"선 데이터 추출 오류: {e}")
        return None


def extract_legend_text(image_array):
    """
    이미지에서 범례 텍스트를 추출하는 함수 (OCR 사용)
    """
    legend_texts = []
    
    if TESSERACT_AVAILABLE:
        try:
            # 범례 영역은 보통 이미지의 오른쪽 상단 또는 하단에 위치
            # 전체 이미지에서 텍스트 추출 시도
            data = pytesseract.image_to_data(image_array, output_type=Output.DICT, lang='eng')
            
            # 텍스트가 있는 영역 찾기
            n_boxes = len(data['text'])
            for i in range(n_boxes):
                text = data['text'][i].strip()
                conf = int(data['conf'][i])
                if text and conf > 30:  # 신뢰도 30 이상
                    # 농도 관련 텍스트 패턴 찾기 (예: "5 ug/mL", "0.5", "10μM" 등)
                    import re
                    # 숫자와 단위가 포함된 텍스트 찾기
                    if re.search(r'\d+\.?\d*\s*(ug/mL|μg/mL|μM|uM|mg/mL|mM|%)', text, re.IGNORECASE) or \
                       re.search(r'^\d+\.?\d*$', text):
                        legend_texts.append(text)
        except Exception as e:
            st.warning(f"Tesseract OCR 오류: {e}")
    
    if EASYOCR_AVAILABLE and len(legend_texts) == 0:
        try:
            reader = easyocr.Reader(['en'], gpu=False)
            results = reader.readtext(image_array)
            
            for (bbox, text, confidence) in results:
                if confidence > 0.5:  # 신뢰도 0.5 이상
                    import re
                    # 농도 관련 텍스트 패턴 찾기
                    if re.search(r'\d+\.?\d*\s*(ug/mL|μg/mL|μM|uM|mg/mL|mM|%)', text, re.IGNORECASE) or \
                       re.search(r'^\d+\.?\d*$', text):
                        legend_texts.append(text)
        except Exception as e:
            st.warning(f"EasyOCR 오류: {e}")
    
    return legend_texts


def convert_image_coords_to_data(x_img, y_img, img_width, img_height, 
                                  x_min, x_max, y_min, y_max,
                                  plot_x_min, plot_x_max, plot_y_min, plot_y_max):
    """
    이미지 좌표를 실제 데이터 좌표로 변환
    
    Args:
        x_img, y_img: 이미지 상의 픽셀 좌표
        img_width, img_height: 이미지 전체 크기
        x_min, x_max, y_min, y_max: 그래프 축의 실제 데이터 범위
        plot_x_min, plot_x_max, plot_y_min, plot_y_max: 그래프 영역의 픽셀 좌표
    """
    # Y축은 이미지 좌표계에서 위가 0이므로 반전 필요
    y_img_flipped = img_height - y_img
    
    # 그래프 영역 내에서의 상대 위치 계산
    x_relative = (x_img - plot_x_min) / (plot_x_max - plot_x_min)
    y_relative = (y_img_flipped - plot_y_min) / (plot_y_max - plot_y_min)
    
    # 실제 데이터 좌표로 변환
    x_data = x_min + x_relative * (x_max - x_min)
    y_data = y_min + y_relative * (y_max - y_min)
    
    return x_data, y_data


def extract_point_data_from_image(image_file, points):
    """
    이미지에서 점 데이터를 추출
    """
    try:
        image = Image.open(image_file)
        img_array = np.array(image)
        
        if points is None or len(points) == 0:
            return None
        
        # 점에서 데이터 포인트 추출 (간단한 예시)
        # 실제로는 좌표 변환 및 축 스케일 추출이 필요
        st.info("💡 점 데이터 추출: Prism 스타일 interpolation을 수행합니다.")
        
        # 여기서는 수동 입력으로 대체
        return None
        
    except Exception as e:
        st.error(f"점 데이터 추출 오류: {e}")
        return None


def exponential_fit_simple(t, F_max, k_obs):
    """
    Exponential fit 모델: F(t) = F_max(1 - e^(-k_obs*t))
    이미지에서 참고한 모델
    """
    return F_max * (1 - np.exp(-k_obs * t))


def _conc_to_export_label(conc_value, exp_type):
    """Export용 파일명 라벨 (예: 0.3125ugmL, 10uM)."""
    if exp_type == "Substrate Concentration Variation (Standard MM)":
        return re.sub(r'[^\w.]', '', f"{conc_value}uM")
    return re.sub(r'[^\w.]', '', f"{conc_value}ugmL")


def _build_experimental_fig(results):
    """Experimental Results 탭과 동일한 figure 생성 (export용)."""
    exp_type = results.get('mm_fit_results', {}).get('experiment_type', 'Substrate Concentration Variation (Standard MM)')
    fig = go.Figure()
    colors = ['blue', 'red', 'orange', 'green', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']
    conc_col = None
    for col in ['Concentration [μM]', 'Concentration [ug/mL]']:
        if col in results['mm_results_df'].columns:
            conc_col = col
            break
    conc_order = results['mm_results_df'].sort_values(conc_col)['Concentration'].tolist() if conc_col else results['mm_results_df']['Concentration'].tolist()
    for idx, conc_name in enumerate(conc_order):
        color = colors[idx % len(colors)]
        conc_match = re.search(r'(\d+\.?\d*)', conc_name)
        legend_name = f"{float(conc_match.group(1))} μM" if conc_match and "Substrate" in exp_type else (f"{float(conc_match.group(1))} μg/mL" if conc_match else conc_name)
        if 'raw_data' in results and conc_name in results['raw_data']:
            raw_conc_data = results['raw_data'][conc_name]
            fig.add_trace(go.Scatter(
                x=raw_conc_data['time'], y=raw_conc_data['value'],
                mode='markers', name=legend_name,
                marker=dict(size=8, color=color, symbol='circle', line=dict(width=1, color='white')),
                legendgroup=conc_name, showlegend=True
            ))
    fig.update_layout(
        xaxis_title='Time (min)', yaxis_title='RFU', height=600, template='plotly_white',
        plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)', hovermode='x unified',
        legend=dict(orientation="v", yanchor="bottom", y=0.05, xanchor="right", x=0.99,
                   bgcolor="rgba(0,0,0,0)", bordercolor="rgba(0,0,0,0)", borderwidth=0, font=dict(color="black"))
    )
    fig.update_xaxes(range=[results['x_data_min'], results['x_data_max']])
    fig.update_yaxes(rangemode='tozero')
    return fig


def _build_norm_fig_full(norm_data, exp_type):
    """정규화 곡선 figure (전체 시간 범위, export용)."""
    t_min = norm_data['times'].min()
    t_max = norm_data['times'].max()
    v0 = norm_data['k_obs'] * (norm_data['Fmax'] - norm_data['F0']) if norm_data.get('k_obs') else None
    v0_label = f"Initial linear region (v₀={v0:.2f} RFU/min)" if v0 is not None else "Initial linear region"
    fig = go.Figure()
    if norm_data.get('k_obs') and norm_data['k_obs'] > 0:
        t_fit = np.linspace(t_min, t_max, 200)
        fit_curve = exponential_fit_simple(t_fit, 1.0, norm_data['k_obs'])
        fig.add_trace(go.Scatter(x=t_fit, y=fit_curve, mode='lines', name='Exponential increase (Full kinetics)', line=dict(color='orange', width=2.5)))
        linear_curve = norm_data['k_obs'] * t_fit
        fig.add_trace(go.Scatter(x=t_fit, y=linear_curve, mode='lines', name=v0_label, line=dict(color='lightblue', width=2.5, dash='dash')))
        tau = norm_data.get('tau')
        if tau and not np.isinf(tau) and tau > 0:
            t_initial, t_exponential_end, t_plateau = 0.1 * tau, 3.0 * tau, 3.0 * tau
            if t_initial <= t_max:
                fig.add_vline(x=t_initial, line_dash="dash", line_color="orange")
            if t_exponential_end <= t_max:
                fig.add_vline(x=0.1 * tau, line_dash="dash", line_color="purple")
                fig.add_vline(x=t_exponential_end, line_dash="dash", line_color="purple")
            if t_plateau <= t_max:
                fig.add_vline(x=t_plateau, line_dash="dash", line_color="brown")
    fig.add_trace(go.Scatter(x=norm_data['times'], y=norm_data['normalized_values'], mode='markers', name='Normalized data', marker=dict(size=6, color='blue')))
    fig.update_layout(
        xaxis_title='Time (min)', yaxis_title='Fluorescence intensity (a.u.)',
        title='Enzyme-quenched peptide fluorescence kinetics', height=600, template='plotly_white',
        plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)', hovermode='x unified',
        yaxis=dict(range=[0, 1.20]), xaxis=dict(range=[t_min, t_max]),
        legend=dict(orientation="v", yanchor="bottom", y=0.05, xanchor="right", x=0.99, bgcolor="rgba(0,0,0,0)", bordercolor="rgba(0,0,0,0)", borderwidth=0)
    )
    return fig


def _build_norm_fig_uptoy1(norm_data, exp_type):
    """정규화 곡선 figure (y=1 도달 구간만, export용). 조건 미충족 시 None."""
    if not norm_data.get('k_obs') or norm_data['k_obs'] <= 0:
        return None
    tau_scaled = norm_data.get('tau')
    if not tau_scaled or np.isinf(tau_scaled) or tau_scaled <= 0:
        return None
    t_min = norm_data['times'].min()
    t_max = norm_data['times'].max()
    normalized_values = norm_data['normalized_values']
    times_norm = norm_data['times']
    target_value = 0.99
    t_y1 = None
    for i, val in enumerate(normalized_values):
        if val >= target_value:
            t_y1 = times_norm[i]
            break
    if t_y1 is None:
        t_y1 = -np.log(1 - target_value) / norm_data['k_obs']
    t_display_max = t_y1
    v0 = norm_data['k_obs'] * (norm_data['Fmax'] - norm_data['F0'])
    v0_label = f"Initial linear region (v₀={v0:.2f} RFU/min)" if v0 else "Initial linear region"
    fig = go.Figure()
    t_fit_scaled = np.linspace(t_min, t_display_max, 200)
    fit_curve_scaled = exponential_fit_simple(t_fit_scaled, 1.0, norm_data['k_obs'])
    fig.add_trace(go.Scatter(x=t_fit_scaled, y=fit_curve_scaled, mode='lines', name='Exponential increase (Full kinetics)', line=dict(color='orange', width=2.5)))
    linear_curve_scaled = norm_data['k_obs'] * t_fit_scaled
    fig.add_trace(go.Scatter(x=t_fit_scaled, y=linear_curve_scaled, mode='lines', name=v0_label, line=dict(color='lightblue', width=2.5, dash='dash')))
    t_initial_scaled, t_exponential_end_scaled, t_plateau_scaled = 0.1 * tau_scaled, 3.0 * tau_scaled, 3.0 * tau_scaled
    if t_initial_scaled <= t_display_max:
        fig.add_vline(x=t_initial_scaled, line_dash="dash", line_color="orange")
    if t_exponential_end_scaled <= t_display_max:
        fig.add_vline(x=0.1 * tau_scaled, line_dash="dash", line_color="purple")
        fig.add_vline(x=t_exponential_end_scaled, line_dash="dash", line_color="purple")
    if t_plateau_scaled <= t_display_max:
        fig.add_vline(x=t_plateau_scaled, line_dash="dash", line_color="brown")
    fig.add_trace(go.Scatter(x=norm_data['times'], y=norm_data['normalized_values'], mode='markers', name='Normalized data', marker=dict(size=6, color='blue')))
    fig.update_layout(
        xaxis_title='Time (min)', yaxis_title='Fluorescence intensity (a.u.)',
        title='Enzyme-quenched peptide fluorescence kinetics (up to y=1)', height=600, template='plotly_white',
        plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)', hovermode='x unified',
        yaxis=dict(range=[0, 1.20]), xaxis=dict(range=[t_min, t_display_max]),
        legend=dict(orientation="v", yanchor="bottom", y=0.05, xanchor="right", x=0.99, bgcolor="rgba(0,0,0,0)", bordercolor="rgba(0,0,0,0)", borderwidth=0)
    )
    return fig


def _build_v0_fig(results):
    """v₀ vs 농도 / Linear fit figure (export용)."""
    v0_data = results['v0_vs_concentration']
    mm_fit = results['mm_fit_results']
    exp_type = mm_fit.get('experiment_type', 'Substrate Concentration Variation (Standard MM)')
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=v0_data['concentrations'], y=v0_data['v0_values'],
        mode='markers', name='Experimental v₀',
        marker=dict(size=10, color='red', line=dict(width=2, color='black'))
    ))
    if exp_type == "Substrate Concentration Variation (Standard MM)" and mm_fit.get('fit_success') and mm_fit.get('Vmax') is not None and mm_fit.get('Km') is not None:
        conc_min, conc_max = min(v0_data['concentrations']), max(v0_data['concentrations'])
        conc_range = np.linspace(conc_min * 0.5, conc_max * 1.5, 200)
        v0_fitted = michaelis_menten_calibration(conc_range, mm_fit['Vmax'], mm_fit['Km'])
        fig.add_trace(go.Scatter(x=conc_range, y=v0_fitted, mode='lines', name=mm_fit.get('equation', ''), line=dict(width=2.5, color='blue')))
        fig.update_layout(
            title='Initial Velocity (v₀) vs Substrate Concentration [S]',
            xaxis_title='[S] (μM)', yaxis_title='Initial Velocity v₀ (Fluorescence Units / Time)',
            template='plotly_white', height=600, hovermode='x unified',
            plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)'
        )
    else:
        if mm_fit.get('fit_success') and mm_fit.get('slope') is not None:
            conc_min, conc_max = min(v0_data['concentrations']), max(v0_data['concentrations'])
            conc_range = np.linspace(conc_min * 0.5, conc_max * 1.5, 200)
            slope = mm_fit['slope']
            intercept = mm_fit.get('intercept', 0)
            v0_fitted = slope * conc_range + intercept
            fig.add_trace(go.Scatter(x=conc_range, y=v0_fitted, mode='lines', name=f'Linear Fit: {mm_fit.get("equation", "")}', line=dict(width=2.5, color='blue', dash='dash')))
        fig.update_layout(
            title='Initial Velocity (v₀) vs Enzyme Concentration [E] (Constant Substrate)',
            xaxis_title='[E] (μg/mL)', yaxis_title='Initial Velocity v₀ (Fluorescence Units / Time)',
            template='plotly_white', height=600, hovermode='x unified',
            plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)'
        )
    return fig


def _build_all_export_figures(results):
    """분석에서 생성되는 모든 figure를 (파일명_접미사제외, fig) 리스트로 반환."""
    out = []
    exp_type = results.get('mm_fit_results', {}).get('experiment_type', 'Substrate Concentration Variation (Standard MM)')
    out.append(("Experimental_Results", _build_experimental_fig(results)))
    if 'normalization_results' in results and results['normalization_results']:
        norm_results = results['normalization_results']
        conc_order = sorted(norm_results.keys(), key=lambda x: norm_results[x]['concentration'])
        for conc_name in conc_order:
            norm_data = norm_results[conc_name]
            conc_label = _conc_to_export_label(norm_data['concentration'], exp_type)
            out.append((f"Normalization_{conc_label}", _build_norm_fig_full(norm_data, exp_type)))
            fig_scaled = _build_norm_fig_uptoy1(norm_data, exp_type)
            if fig_scaled is not None:
                out.append((f"Normalization_to_plateau_{conc_label}", fig_scaled))
    if 'v0_vs_concentration' in results and results.get('mm_fit_results', {}).get('fit_success'):
        if exp_type == "Substrate Concentration Variation (Standard MM)":
            out.append(("v0_vs_S_Fit", _build_v0_fig(results)))
        else:
            out.append(("Linear_fit", _build_v0_fig(results)))
    return out


def normalize_iterative(times, values, num_iterations=2):
    """
    반복 정규화 수행
    
    1차 정규화: F0 = value at time 0, Fmax = max(F)
    2차 정규화: Exponential fit F(t) = F_max(1 - e^(-k_obs*t))
    
    Parameters:
    - times: 시간 배열
    - values: 형광값 배열
    - num_iterations: 반복 횟수 (최소 2번)
    
    Returns:
    - normalized_times: 정규화된 시간 배열
    - normalized_values: 정규화된 형광값 배열
    - F0: 최종 F0 값
    - Fmax: 최종 Fmax 값
    - k_obs: 최종 k_obs 값
    - tau: 최종 τ = 1/k_obs 값
    - r_squared: 최종 R² 값
    - equation: 방정식 문자열
    """
    times = np.array(times)
    values = np.array(values)
    
    # 정렬 (시간 순서대로)
    sort_idx = np.argsort(times)
    times = times[sort_idx]
    values = values[sort_idx]
    
    # 초기값
    current_values = values.copy()
    F0 = None
    Fmax = None
    k_obs = None
    tau = None
    r_squared = 0
    equation = ""
    
    # 반복 정규화 (최소 2번)
    for iteration in range(max(2, num_iterations)):
        # 1차 정규화: F0 = value at time 0, Fmax = max(F)
        F0 = current_values[0]  # time이 0일 때의 값
        Fmax = np.max(current_values)  # max(F)
        
        # 정규화: (F - F0) / (Fmax - F0)
        if Fmax > F0:
            normalized = (current_values - F0) / (Fmax - F0)
        else:
            normalized = current_values - F0
        
        # 2차 정규화: Exponential fit
        # F(t) = F_max(1 - e^(-k_obs*t))
        # 정규화된 데이터에 대해 피팅
        try:
            # 초기값 추정
            F_max_init = 1.0  # 정규화된 값이므로 1.0
            k_obs_init = 0.1  # 초기 추정값
            
            # Exponential fit
            popt, pcov = curve_fit(
                exponential_fit_simple,
                times,
                normalized,
                p0=[F_max_init, k_obs_init],
                bounds=([0.1, 0.001], [2.0, 10.0]),
                maxfev=5000
            )
            
            F_max_fit, k_obs_fit = popt
            
            # 피팅된 값 계산
            fit_values = exponential_fit_simple(times, F_max_fit, k_obs_fit)
            
            # R² 계산
            ss_res = np.sum((normalized - fit_values) ** 2)
            ss_tot = np.sum((normalized - np.mean(normalized)) ** 2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            
            # τ = 1/k_obs
            tau = 1.0 / k_obs_fit if k_obs_fit > 0 else np.inf
            
            k_obs = k_obs_fit
            
            # 방정식 생성
            equation = f"F(t) = {F_max_fit:.4f}(1 - e^(-{k_obs_fit:.4f}*t)), τ = {tau:.4f}"
            
            # 다음 반복을 위해 피팅된 값을 역정규화하여 사용 (정규화 개선)
            if iteration < max(2, num_iterations) - 1:  # 마지막 반복이 아니면
                # 역정규화: fit_values * (Fmax - F0) + F0
                if Fmax > F0:
                    current_values = fit_values * (Fmax - F0) + F0
                else:
                    current_values = fit_values + F0
            
        except Exception as e:
            # 피팅 실패 시 정규화된 값 유지
            if iteration == 0:
                # 첫 반복에서 실패하면 기본값 사용
                k_obs = 0.1
                tau = 10.0
                equation = f"F(t) = {Fmax:.2f}(1 - e^(-{k_obs:.4f}*t)), τ = {tau:.4f} (피팅 실패)"
    
    return times, normalized, F0, Fmax, k_obs, tau, r_squared, equation


def manual_data_entry(data_type="점"):
    """
    수동으로 데이터 포인트를 입력받는 함수
    data_type: "점" 또는 "선"
    """
    st.subheader(f"📝 수동 데이터 입력 ({data_type} 데이터)")
    
    num_curves = st.number_input("곡선 개수 (농도 조건 수)", min_value=1, max_value=20, value=1)
    
    all_curves_data = {}
    
    for curve_idx in range(num_curves):
        with st.expander(f"곡선 {curve_idx + 1} (농도 조건)", expanded=(curve_idx == 0)):
            conc_name = st.text_input(f"농도 이름 {curve_idx + 1}", value=f"{curve_idx + 1} ug/mL", key=f"conc_{curve_idx}")
            conc_value = st.number_input(f"농도 값 (ug/mL) {curve_idx + 1}", value=float(curve_idx + 1), step=0.1, key=f"conc_val_{curve_idx}")
            
            num_points = st.number_input(f"데이터 포인트 개수 {curve_idx + 1}", min_value=2, max_value=100, value=10, key=f"num_{curve_idx}")
            
            data_points = []
            cols = st.columns(2)
            
            with cols[0]:
                st.write("**시간 (min)**")
            with cols[1]:
                st.write("**RFU 값**")
            
            for i in range(num_points):
                cols = st.columns(2)
                with cols[0]:
                    time_val = st.number_input(f"시간 {i+1}", key=f"time_{curve_idx}_{i}", value=float(i*5), step=0.1)
                with cols[1]:
                    rfu_val = st.number_input(f"RFU {i+1}", key=f"rfu_{curve_idx}_{i}", value=float(100+i*10), step=0.1)
                
                data_points.append({'Time_min': time_val, 'RFU': rfu_val})
            
            all_curves_data[conc_name] = {
                'concentration': conc_value,
                'data': data_points
            }
    
    if st.button("데이터 확인", key="confirm_data"):
        return all_curves_data
    
    return None


def data_load_mode(st):
    """Data Load 모드 - CSV 파일 업로드 또는 이미지에서 데이터 추출"""
    
    # 폴더 구조 생성
    os.makedirs("prep_raw_data_mode", exist_ok=True)
    os.makedirs("prep_raw_data_mode/results", exist_ok=True)
    os.makedirs("data_interpolation_mode/results", exist_ok=True)
    
    st.header("📥 Data Load Mode")
    st.markdown("---")
    
    # Sidebar configuration
    st.sidebar.title("⚙️ Data Load Settings")
    
    # Experiment condition selection (before file upload)
    st.sidebar.subheader("🔬 Experiment Condition")
    experiment_type = st.sidebar.radio(
        "Experiment Type",
        ["Substrate Concentration Variation (Standard MM)", "Enzyme Concentration Variation (Fixed substrate)"],
        help="Substrate Concentration Variation: Standard MM applicable | Enzyme Concentration Variation: Standard MM not applicable, linear relationship"
    )
    
    # Determine sample file path based on experiment type
    if experiment_type == "Substrate Concentration Variation (Standard MM)":
        sample_file_path = "raw/raw_substrate.csv"
        sample_file_name = "raw_substrate_sample.xlsx"
        sample_file_label = "Download Sample raw_substrate.xlsx"
    else:  # Enzyme Concentration Variation (Fixed substrate)
        sample_file_path = "raw/raw_enzyme.csv"
        sample_file_name = "raw_enzyme_sample.xlsx"
        sample_file_label = "Download Sample raw_enzyme.xlsx"
    
    # CSV/XLSX file upload
    st.sidebar.subheader("📁 Data File Upload")
    uploaded_file = st.sidebar.file_uploader(
        "Upload Prep Raw Data File (CSV or XLSX)",
        type=['csv', 'xlsx'],
        help="prep_raw.csv/xlsx format: Time, concentration values, SD, replicates (3 columns each)"
    )
    
    # 샘플 데이터 다운로드 (XLSX 형태로 제공, 실험 타입에 따라 다른 파일)
    try:
        import io
        df_sample = pd.read_csv(sample_file_path, header=None, sep='\t')
        output = io.BytesIO()
        with pd.ExcelWriter(output, engine='openpyxl') as writer:
            df_sample.to_excel(writer, index=False, header=False, sheet_name='Sheet1')
        sample_bytes = output.getvalue()
        st.sidebar.download_button(
            label=sample_file_label,
            data=sample_bytes,
            file_name=sample_file_name,
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )
    except Exception:
        pass
    
    # 데이터 로드
    if uploaded_file is not None:
        # 업로드된 파일을 임시로 저장하고 읽기
        import tempfile
        
        # 파일 확장자 확인
        file_extension = uploaded_file.name.split('.')[-1].lower()
        suffix = f'.{file_extension}'
        
        with tempfile.NamedTemporaryFile(delete=False, suffix=suffix, mode='wb') as tmp_file:
            tmp_file.write(uploaded_file.getbuffer())
            tmp_path = tmp_file.name
        
        try:
            raw_data = read_raw_data(tmp_path)
            os.unlink(tmp_path)
        except Exception as e:
            st.error(f"File reading error: {e}")
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)
            return
    else:
        # 기본 샘플 데이터 사용 (실험 타입에 따라 다른 파일)
        from pathlib import Path
        
        # Determine sample file path based on experiment type
        if experiment_type == "Substrate Concentration Variation (Standard MM)":
            default_sample_paths = [
                'raw/raw_substrate.csv',  # Current working directory
                str(Path(__file__).parent.parent / 'raw' / 'raw_substrate.csv'),  # Script directory
            ]
            default_sample_name = "raw/raw_substrate.csv"
        else:  # Enzyme Concentration Variation (Fixed substrate)
            default_sample_paths = [
                'raw/raw_enzyme.csv',  # Current working directory
                str(Path(__file__).parent.parent / 'raw' / 'raw_enzyme.csv'),  # Script directory
            ]
            default_sample_name = "raw/raw_enzyme.csv"
        
        raw_data = None
        used_path = None
        
        for path in default_sample_paths:
            try:
                if os.path.exists(path):
                    raw_data = read_raw_data(path)
                    used_path = path
                    break
            except Exception:
                continue
        
        if raw_data is None:
            # Last attempt: search in current working directory
            try:
                raw_data = read_raw_data(default_sample_name)
                st.sidebar.info(f"Using {default_sample_name}")
            except Exception as e:
                st.error(f"Data file not found. Please upload a CSV or XLSX file.\nError: {str(e)}")
                st.stop()
        else:
            st.sidebar.info(f"Using {used_path}")
    
    # Data preview
    st.subheader("📋 Data Preview")
    
    # 반응 시간 계산 (최대값)
    all_times = [time_val for data in raw_data.values() for time_val in data['time']] if raw_data else []
    reaction_time = f"{max(all_times):.0f} min" if all_times else "N/A"
    
    # N 값 읽기
    try:
        if uploaded_file is not None:
            uploaded_file.seek(0)
            first_line = uploaded_file.readline().decode('utf-8')
            second_line = uploaded_file.readline().decode('utf-8')
            third_line = uploaded_file.readline().decode('utf-8')
            n_value = int(third_line.split('\t')[3]) if len(third_line.split('\t')) > 3 else 50
            uploaded_file.seek(0)
        else:
            # 실험 타입에 따라 다른 샘플 파일 사용
            if experiment_type == "Substrate Concentration Variation (Standard MM)":
                default_n_file = 'raw/raw_substrate.csv'
            else:  # Enzyme 농도 변화 (Substrate 고정)
                default_n_file = 'raw/raw_enzyme.csv'
            
            with open(default_n_file, 'r', encoding='utf-8') as f:
                f.readline()
                f.readline()
                third_line = f.readline()
                n_value = int(third_line.split('\t')[3]) if len(third_line.split('\t')) > 3 else 50
    except:
        n_value = 50
    
    # Show error message if raw_data is not available
    if not raw_data:
        st.error("Unable to load data. Please upload a CSV or XLSX file.")
        return
    
    # Calculate number of data points per concentration (same for all concentrations)
    sorted_conc = sorted(raw_data.items(), key=lambda x: x[1]['concentration'])
    num_data_points = len(sorted_conc[0][1]['time']) if len(sorted_conc) > 0 else 0
    
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Number of Concentrations", len(raw_data))
    with col2:
        st.metric("Data Points per Concentration", num_data_points)
    with col3:
        st.metric("Reaction Time", reaction_time)
    with col4:
        st.metric("N (Number of Replicates)", n_value)
    
    # Display concentration-specific information
    with st.expander("Concentration Data Information", expanded=False):
        if len(sorted_conc) > 0:
            first_data = sorted_conc[0][1]
            times = first_data['time']
            
            detail_data = {'time_min': times}
            for conc_name, data in sorted_conc:
                conc_label = f"{data['concentration']}"
                detail_data[f'{conc_label}_mean'] = data['value']
                if data.get('SD') is not None:
                    detail_data[f'{conc_label}_SD'] = data['SD']
            
            detail_df = pd.DataFrame(detail_data)
            st.dataframe(detail_df, use_container_width=True, hide_index=True, height=400)
        else:
            st.info("No data available.")
    
    if experiment_type == "Enzyme Concentration Variation (Fixed substrate)":
        st.sidebar.warning("""
        ⚠️ **Constant Substrate + Enzyme Concentration Variation**
        
        - v is **linearly** related to [E]
        - **Km cannot be determined** (substrate concentration gradient required)
        - **Vmax cannot be determined** (standard MM requires fixed [E])
        - Obtainable parameters:
          - **slope = kcat × [S] / (Km + [S])**
          - At very low [S]: **slope ≈ kcat/Km × [S]**
        """)
    
    # Enzyme 농도 입력 (kcat 계산용, Substrate 농도 변화 실험에서만 필요)
    if experiment_type == "Substrate Concentration Variation (Standard MM)":
        st.sidebar.subheader("🧪 Enzyme concentration setting (for kcat calculation)")
        enzyme_conc_input = st.sidebar.number_input(
            "Enzyme concentration [E] (μg/mL)",
            min_value=0.0,
            value=51.43,
            step=0.1,
            help="Required for kcat = Vmax / [E]_T. Enter the enzyme concentration used in the experiment."
        )
    else:
        enzyme_conc_input = None
    
    # Michaelis-Menten 모델 실행 버튼
    if st.button("🚀 Run Michaelis-Menten Model", type="primary"):
            with st.spinner("Running Michaelis-Menten model fitting..."):
                progress_bar = st.progress(0)
                status_text = st.empty()
                
                # 1. 각 농도별 시간 경과 곡선 피팅
                status_text.text("1️⃣ Fitting time-course curves per concentration...")
                progress_bar.progress(0.2)
                
                mm_results = {}
                all_fit_data = []
                
                for conc_name, data in raw_data.items():
                    times = data['time']
                    values = data['value']
                    
                    # 초기 속도 계산 (최적화된 방법: (F∞-F0)의 5-10% 범위에서 R² 최대 구간 선택)
                    params, fit_values, r_sq = fit_time_course(times, values, model='linear', use_optimized=True)
                    
                    # 선형 구간 데이터는 params에서 가져오기 (linear_times, linear_values는 나중에 저장됨)
                    # 최적화된 방법에서는 calculate_initial_velocity_optimized가 이미 호출됨
                    from mode_prep_raw_data.prep import calculate_initial_velocity_optimized
                    v0_calc, F0_calc, r_sq_calc, linear_times, linear_values, conversion_used = calculate_initial_velocity_optimized(times, values)
                    optimal_percent = conversion_used * 100 if conversion_used is not None else None
                    
                    # 초기 속도 파라미터 추출
                    v0 = params['v0']  # 초기 속도
                    F0 = params['F0']  # 초기 형광값
                    Fmax = params['Fmax']  # 최대 형광값
                    
                    mm_results[conc_name] = {
                        'concentration': data['concentration'],
                        'v0': v0,
                        'F0': F0,
                        'Fmax': Fmax,
                        'R_squared': r_sq,
                        'linear_fraction': params['linear_fraction'],
                        'optimal_percent': optimal_percent,  # 최적화된 퍼센트
                        'linear_times': linear_times,  # 초기속도 탭용
                        'linear_values': linear_values,  # 초기속도 탭용
                        'times': times,  # 원본 시간 데이터
                        'values': values  # 원본 형광 데이터
                    }
                    
                    # Fit curve 데이터 저장 (선형 구간만)
                    valid_mask = ~np.isnan(fit_values)
                    # 농도 단위 결정: 실험 타입에 따라
                    # Substrate 농도 변화: uM (몰농도)
                    # Enzyme 농도 변화: ug/mL (질량 농도)
                    if experiment_type == "Substrate Concentration Variation (Standard MM)":
                        conc_unit_col = 'Concentration [μM]'
                    else:  # Enzyme 농도 변화
                        conc_unit_col = 'Concentration [ug/mL]'
                    
                    for t, val, fit_val in zip(times[valid_mask], values[valid_mask], fit_values[valid_mask]):
                        fit_row = {
                            'Concentration': conc_name,
                            'Time_min': t,
                            'Observed_Value': val,
                            'Fit_Value': fit_val,
                            'Residual': val - fit_val
                        }
                        fit_row[conc_unit_col] = data['concentration']
                        all_fit_data.append(fit_row)
                
                progress_bar.progress(0.4)
                
                # 2. Calculate interpolation range
                status_text.text("2️⃣ Calculating interpolation range...")
                
                all_times = [time_val for data in raw_data.values() for time_val in data['time']]
                x_data_min = min(all_times)
                x_data_max = max(all_times)
                # 원본 데이터 범위만 사용 (Prism 확장 범위 사용 안 함)
                x_range_min = x_data_min
                x_range_max = x_data_max
                
                # 보간 포인트 개수 설정 (고정값 사용)
                n_points = 1000  # 기본값으로 고정
                
                # 고밀도 보간 포인트 생성
                x_interp = np.linspace(x_range_min, x_range_max, n_points + 1)
                
                progress_bar.progress(0.6)
                
                # 3. Perform interpolation (using normalization results)
                status_text.text("3️⃣ Generating interpolation curves...")
                
                all_interp_data = []
                # 농도 단위 결정: 실험 타입에 따라
                # Substrate 농도 변화: uM (몰농도)
                # Enzyme 농도 변화: ug/mL (질량 농도)
                if experiment_type == "Substrate Concentration Variation (Standard MM)":
                    conc_unit_col = 'Concentration [μM]'
                else:  # Enzyme 농도 변화
                    conc_unit_col = 'Concentration [ug/mL]'
                
                # 정규화 결과가 있는지 확인 (나중에 생성되므로 임시로 mm_results 사용)
                # 정규화 결과는 나중에 생성되므로, 여기서는 mm_results를 사용하되
                # 정규화 결과가 있으면 그것을 사용하도록 수정 필요
                # 일단 mm_results를 사용하되, 정규화 결과가 생성된 후에는 그것을 사용
                
                for conc_name, params in mm_results.items():
                    v0 = params['v0']
                    F0 = params['F0']
                    Fmax = params['Fmax']
                    
                    # 실험 타입에 따라 다른 보간 방법 사용
                    if experiment_type == "Substrate Concentration Variation (Standard MM)":
                        # 기질 농도 변화: Michaelis-Menten kinetics에 따른 지수 결합 모델 사용
                        # F(t) = F0 + (Fmax - F0) * [1 - exp(-k*t)]
                        # where k = v0 / (Fmax - F0)
                        if Fmax > F0 and v0 > 0:
                            k = v0 / (Fmax - F0)
                            y_interp = F0 + (Fmax - F0) * (1 - np.exp(-k * x_interp))
                        else:
                            # Fallback to linear if invalid parameters
                            y_interp = F0 + v0 * x_interp
                    else:
                        # Enzyme 농도 변화: 선형 관계 사용 (기질 고정)
                        # F(t) = F0 + v0 * t
                        y_interp = F0 + v0 * x_interp
                    
                    # Fmax를 넘지 않도록 제한
                    y_interp = np.clip(y_interp, F0, Fmax)
                    
                    for x, y in zip(x_interp, y_interp):
                        interp_row = {
                            'Concentration': conc_name,
                            'Time_min': x,
                            'RFU_Interpolated': y
                        }
                        interp_row[conc_unit_col] = params['concentration']
                        all_interp_data.append(interp_row)
                
                interp_df = pd.DataFrame(all_interp_data)
                
                progress_bar.progress(0.7)
                
                # 4. Fit v₀ vs concentration (varies by experiment condition)
                if experiment_type == "Substrate Concentration Variation (Standard MM)":
                    status_text.text("4️⃣ Fitting v₀ vs [S] Michaelis-Menten...")
                    
                    # 농도와 초기 속도 데이터 수집
                    concentrations = [params['concentration'] for params in sorted(mm_results.values(), 
                                                                                  key=lambda x: x['concentration'])]
                    v0_values = [params['v0'] for params in sorted(mm_results.values(), 
                                                                  key=lambda x: x['concentration'])]
                    
                    # MM calibration curve 피팅: v₀ = Vmax * [S] / (Km + [S])
                    if len(concentrations) >= 2 and len(v0_values) >= 2:
                        try:
                            cal_params, cal_fit_values, cal_equation = fit_calibration_curve(concentrations, v0_values)
                            Vmax = cal_params['Vmax_cal']
                            Km = cal_params['Km_cal']
                            mm_r_squared = cal_params['R_squared']
                            
                            # kcat 계산 (enzyme 농도 필요)
                            kcat = None
                            mm_fit_success = True
                        except Exception as e:
                            st.warning(f"⚠️ MM fitting failed: {e}")
                            Vmax = None
                            Km = None
                            kcat = None
                            mm_r_squared = 0
                            cal_equation = "Fitting failed"
                            mm_fit_success = False
                    else:
                        Vmax = None
                        Km = None
                        kcat = None
                        mm_r_squared = 0
                        cal_equation = "Insufficient data (at least 2 concentrations required)"
                        mm_fit_success = False
                
                else:  # Enzyme Concentration Variation (Fixed substrate)
                    status_text.text("4️⃣ Fitting v₀ vs [E] linear... (not standard MM)")
                    
                    # 농도와 초기 속도 데이터 수집
                    concentrations = [params['concentration'] for params in sorted(mm_results.values(), 
                                                                                  key=lambda x: x['concentration'])]
                    v0_values = [params['v0'] for params in sorted(mm_results.values(), 
                                                                  key=lambda x: x['concentration'])]
                    
                    # 선형 피팅: v = kcat * [E] * [S] / (Km + [S])
                    # Substrate 고정이므로 slope = kcat * [S] / (Km + [S])
                    if len(concentrations) >= 2 and len(v0_values) >= 2:
                        try:
                            # 선형 회귀
                            coeffs = np.polyfit(concentrations, v0_values, 1)
                            slope = coeffs[0]  # kcat * [S] / (Km + [S])
                            intercept = coeffs[1]
                            
                            # 피팅된 값
                            v0_fitted = np.polyval(coeffs, concentrations)
                            
                            # R² 계산
                            ss_res = np.sum((v0_values - v0_fitted) ** 2)
                            ss_tot = np.sum((v0_values - np.mean(v0_values)) ** 2)
                            mm_r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
                            
                            # Substrate 고정 조건에서는 Km을 구할 수 없음
                            Vmax = None  # Vmax는 [E] 고정 조건에서만 정의됨
                            Km = None  # 구할 수 없음
                            kcat = None  # 단독으로 구할 수 없음 (kcat/Km만 가능)
                            
                            # slope = kcat * [S] / (Km + [S])
                            # Substrate 농도가 알려져 있으면 kcat/Km을 추정할 수 있음 (희석 조건)
                            cal_equation = f"v₀ = {slope:.4f} * [E] + {intercept:.4f} (Linear)"
                            mm_fit_success = True
                        except Exception as e:
                            st.warning(f"⚠️ Linear fitting failed: {e}")
                            Vmax = None
                            Km = None
                            kcat = None
                            mm_r_squared = 0
                            cal_equation = "Fitting failed"
                            mm_fit_success = False
                    else:
                        Vmax = None
                        Km = None
                        kcat = None
                        mm_r_squared = 0
                        cal_equation = "Insufficient data (at least 2 concentrations required)"
                        mm_fit_success = False
                
                progress_bar.progress(0.85)
                
                # 5. 결과 저장
                status_text.text("5️⃣ Saving results...")
                
                # 초기 속도 Results 저장 (MM 파라미터 포함)
                # 농도 단위 결정: 실험 타입에 따라
                # Substrate 농도 변화: uM (몰농도)
                # Enzyme 농도 변화: ug/mL (질량 농도)
                if experiment_type == "Substrate Concentration Variation (Standard MM)":
                    conc_unit_col = 'Concentration [μM]'
                else:  # Enzyme 농도 변화
                    conc_unit_col = 'Concentration [ug/mL]'
                
                results_data = []
                for conc_name, params in sorted(mm_results.items(), key=lambda x: x[1]['concentration']):
                    eq = f"v0 = {params['v0']:.2f} (선형 구간 기울기)"
                    result_row = {
                        'Concentration': conc_name,
                        'v0': params['v0'],
                        'F0': params['F0'],
                        'Fmax': params['Fmax'],
                        'R_squared': params['R_squared'],
                        'linear_fraction': params['linear_fraction'],
                        'Equation': eq
                    }
                    result_row[conc_unit_col] = params['concentration']
                    results_data.append(result_row)
                
                mm_results_df = pd.DataFrame(results_data)
                
                # Enzyme 농도 가져오기 (kcat 계산용)
                # 우선순위: 1) 사용자 입력값, 2) xlsx 파일에서 읽기
                enzyme_conc = None
                
                # 1) 사용자 입력값 확인
                if experiment_type == "Substrate Concentration Variation (Standard MM)" and enzyme_conc_input is not None and enzyme_conc_input > 0:
                    enzyme_conc = enzyme_conc_input
                
                # 2) xlsx 파일에서 읽기 시도 (사용자 입력값이 없을 때만)
                if enzyme_conc is None:
                    try:
                        xlsx_path = 'Michaelis-Menten_calibration_results.xlsx'
                        if os.path.exists(xlsx_path):
                            df_mm_read = pd.read_excel(xlsx_path, sheet_name='MM Results', engine='openpyxl')
                            # enzyme 농도 컬럼 찾기 (다양한 이름 시도)
                            enzyme_conc_col = None
                            for col in ['Enzyme [ug/mL]', 'Enzyme_ug/mL', 'enzyme_ug/mL', '[E] (ug/mL)', 'E_conc', 'Enzyme']:
                                if col in df_mm_read.columns:
                                    enzyme_conc_col = col
                                    break
                            
                            if enzyme_conc_col is not None:
                                # 첫 번째 유효한 enzyme 농도 값 사용
                                enzyme_conc_values = df_mm_read[enzyme_conc_col].dropna()
                                if len(enzyme_conc_values) > 0:
                                    enzyme_conc = float(enzyme_conc_values.iloc[0])
                    except Exception as e:
                        # enzyme 농도 읽기 실패해도 계속 진행
                        pass
                
                # kcat 계산: kcat = Vmax / [E]_T
                if mm_fit_success and Vmax is not None and enzyme_conc is not None and enzyme_conc > 0:
                    kcat = Vmax / enzyme_conc
                else:
                    kcat = None
                    if experiment_type == "Substrate Concentration Variation (Standard MM)" and mm_fit_success and Vmax is not None:
                        st.sidebar.warning("⚠️ Please enter enzyme concentration for kcat calculation.")
                
                # MM 피팅 결과를 별도로 저장
                mm_fit_results = {
                    'Vmax': Vmax,
                    'Km': Km,
                    'kcat': kcat,
                    'enzyme_conc': enzyme_conc,
                    'R_squared': mm_r_squared,
                    'equation': cal_equation,
                    'fit_success': mm_fit_success,
                    'experiment_type': experiment_type,
                    'slope': None  # Enzyme 농도 변화인 경우 slope 저장
                }
                
                # Enzyme 농도 변화인 경우 slope 저장
                if experiment_type == "Enzyme Concentration Variation (Fixed substrate)" and mm_fit_success:
                    concentrations = [params['concentration'] for params in sorted(mm_results.values(), 
                                                                                  key=lambda x: x['concentration'])]
                    v0_values = [params['v0'] for params in sorted(mm_results.values(), 
                                                                  key=lambda x: x['concentration'])]
                    if len(concentrations) >= 2:
                        coeffs = np.polyfit(concentrations, v0_values, 1)
                        mm_fit_results['slope'] = coeffs[0]
                        mm_fit_results['intercept'] = coeffs[1]
                
                try:
                    # Save interpolated curves (CSV)
                    interp_df.to_csv('data_interpolation_mode/results/MM_interpolated_curves.csv', index=False)
                    
                    # Save MM results (CSV)
                    mm_results_df.to_csv('prep_raw_data_mode/results/MM_results_detailed.csv', index=False)
                    
                    st.sidebar.success("✅ Result files saved!")
                except Exception as e:
                    st.sidebar.warning(f"⚠️ Error saving files: {e}")
                
                progress_bar.progress(0.9)
                
                # 6. Perform normalization
                status_text.text("6️⃣ Normalizing...")
                
                normalization_results = {}
                # 정규화 기반 v0 값들을 저장할 딕셔너리
                norm_v0_values = {}
                
                for conc_name, data in raw_data.items():
                    times = data['time']
                    values = data['value']
                    
                    # 반복 정규화 수행 (최소 2번)
                    norm_times, norm_values, F0, Fmax, k_obs, tau, r_sq, equation = normalize_iterative(
                        times, values, num_iterations=2
                    )
                    
                    # 정규화 기반 v0 계산: v0 = k_obs * (Fmax - F0)
                    v0_norm = k_obs * (Fmax - F0) if k_obs is not None and k_obs > 0 else 0
                    norm_v0_values[conc_name] = v0_norm
                    
                    normalization_results[conc_name] = {
                        'concentration': data['concentration'],
                        'times': norm_times,
                        'normalized_values': norm_values,
                        'F0': F0,
                        'Fmax': Fmax,
                        'k_obs': k_obs,
                        'tau': tau,
                        'R_squared': r_sq,
                        'equation': equation,
                        'original_times': times,
                        'original_values': values,
                        'v0': v0_norm  # 정규화 기반 v0 추가
                    }
                
                # 정규화 기반 v0으로 MM fit 다시 수행
                status_text.text("7️⃣ Re-running MM fitting with normalized v₀...")
                
                # 정규화 기반 v0 값들로 농도와 v0 데이터 수집
                norm_concentrations = []
                norm_v0_list = []
                
                for conc_name in sorted(normalization_results.keys(), 
                                       key=lambda x: normalization_results[x]['concentration']):
                    norm_concentrations.append(normalization_results[conc_name]['concentration'])
                    norm_v0_list.append(normalization_results[conc_name]['v0'])
                
                # Enzyme 농도 가져오기 (kcat 계산용)
                # 우선순위: 1) 사용자 입력값, 2) xlsx 파일에서 읽기
                norm_enzyme_conc = None
                
                # 1) 사용자 입력값 확인
                if experiment_type == "Substrate Concentration Variation (Standard MM)" and enzyme_conc_input is not None and enzyme_conc_input > 0:
                    norm_enzyme_conc = enzyme_conc_input
                
                # 2) xlsx 파일에서 읽기 시도 (사용자 입력값이 없을 때만)
                if norm_enzyme_conc is None:
                    try:
                        xlsx_path = 'Michaelis-Menten_calibration_results.xlsx'
                        if os.path.exists(xlsx_path):
                            df_mm_read = pd.read_excel(xlsx_path, sheet_name='MM Results', engine='openpyxl')
                            # enzyme 농도 컬럼 찾기 (다양한 이름 시도)
                            enzyme_conc_col = None
                            for col in ['Enzyme [ug/mL]', 'Enzyme_ug/mL', 'enzyme_ug/mL', '[E] (ug/mL)', 'E_conc', 'Enzyme']:
                                if col in df_mm_read.columns:
                                    enzyme_conc_col = col
                                    break
                            
                            if enzyme_conc_col is not None:
                                # 첫 번째 유효한 enzyme 농도 값 사용
                                enzyme_conc_values = df_mm_read[enzyme_conc_col].dropna()
                                if len(enzyme_conc_values) > 0:
                                    norm_enzyme_conc = float(enzyme_conc_values.iloc[0])
                    except Exception as e:
                        # enzyme 농도 읽기 실패해도 계속 진행
                        pass
                
                # MM fit 재수행 (정규화 기반 v0 사용)
                if experiment_type == "Substrate Concentration Variation (Standard MM)":
                    if len(norm_concentrations) >= 2 and len(norm_v0_list) >= 2:
                        try:
                            cal_params, cal_fit_values, cal_equation = fit_calibration_curve(norm_concentrations, norm_v0_list)
                            Vmax = cal_params['Vmax_cal']
                            Km = cal_params['Km_cal']
                            mm_r_squared = cal_params['R_squared']
                            
                            # kcat 계산: kcat = Vmax / [E]_T
                            if Vmax is not None and norm_enzyme_conc is not None and norm_enzyme_conc > 0:
                                kcat = Vmax / norm_enzyme_conc
                            else:
                                kcat = None
                                if Vmax is not None:
                                    st.sidebar.warning("⚠️ Please enter enzyme concentration for kcat calculation.")
                            
                            mm_fit_success = True
                            
                            # mm_fit_results 업데이트
                            mm_fit_results = {
                                'Vmax': Vmax,
                                'Km': Km,
                                'kcat': kcat,
                                'enzyme_conc': norm_enzyme_conc,
                                'R_squared': mm_r_squared,
                                'equation': cal_equation,
                                'fit_success': mm_fit_success,
                                'experiment_type': experiment_type,
                                'slope': None
                            }
                        except Exception as e:
                            st.warning(f"⚠️ Normalized MM fitting failed: {e}")
                            mm_fit_success = False
                    else:
                        mm_fit_success = False
                else:  # Enzyme 농도 변화
                    if len(norm_concentrations) >= 2 and len(norm_v0_list) >= 2:
                        try:
                            # 선형 회귀
                            coeffs = np.polyfit(norm_concentrations, norm_v0_list, 1)
                            slope = coeffs[0]
                            intercept = coeffs[1]
                            
                            # 피팅된 값
                            v0_fitted = np.polyval(coeffs, norm_concentrations)
                            
                            # R² 계산
                            ss_res = np.sum((norm_v0_list - v0_fitted) ** 2)
                            ss_tot = np.sum((norm_v0_list - np.mean(norm_v0_list)) ** 2)
                            mm_r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
                            
                            cal_equation = f"v₀ = {slope:.4f} * [E] + {intercept:.4f} (Linear)"
                            mm_fit_success = True
                            
                            # mm_fit_results 업데이트
                            mm_fit_results = {
                                'Vmax': None,
                                'Km': None,
                                'kcat': None,
                                'enzyme_conc': None,
                                'R_squared': mm_r_squared,
                                'equation': cal_equation,
                                'fit_success': mm_fit_success,
                                'experiment_type': experiment_type,
                                'slope': slope,
                                'intercept': intercept
                            }
                        except Exception as e:
                            st.warning(f"⚠️ Normalization-based linear fitting failed: {e}")
                            mm_fit_success = False
                    else:
                        mm_fit_success = False
                
                # 정규화 결과를 사용하여 보간 곡선 재생성 (Exponential 식 사용)
                status_text.text("8️⃣ Regenerating interpolation curves from normalization...")
                
                if 'normalization_results' in locals() and normalization_results:
                    all_interp_data_new = []
                    
                    for conc_name in sorted(normalization_results.keys(), 
                                           key=lambda x: normalization_results[x]['concentration']):
                        n_data = normalization_results[conc_name]
                        F0 = n_data['F0']
                        Fmax = n_data['Fmax']
                        k_obs = n_data.get('k_obs', None)
                        
                        if k_obs is not None and k_obs > 0:
                            # Exponential 식 사용: F(t) = F0 + (Fmax - F0) * [1 - exp(-k_obs * t)]
                            # normalize_iterative에서 times는 원본 시간 단위(보통 분)를 사용하므로
                            # k_obs는 분^-1 단위입니다. x_interp도 분 단위이므로 그대로 사용
                            y_interp = F0 + (Fmax - F0) * (1 - np.exp(-k_obs * x_interp))
                            # Fmax를 넘지 않도록 제한
                            y_interp = np.clip(y_interp, F0, Fmax)
                        else:
                            # k_obs가 없으면 선형 보간 사용 (fallback)
                            v0 = n_data.get('v0', 0)
                            y_interp = F0 + v0 * x_interp
                            y_interp = np.clip(y_interp, F0, Fmax)
                        
                        conc_value = n_data['concentration']
                        
                        for x, y in zip(x_interp, y_interp):
                            interp_row = {
                                'Concentration': conc_name,
                                'Time_min': x,
                                'RFU_Interpolated': y
                            }
                            interp_row[conc_unit_col] = conc_value
                            all_interp_data_new.append(interp_row)
                    
                    # 새로운 보간 데이터로 업데이트
                    interp_df = pd.DataFrame(all_interp_data_new)
                
                progress_bar.progress(1.0)
                status_text.text("✅ Michaelis-Menten model fitting and normalization complete!")
                
            # Session state에 저장 (정규화 기반 v0 사용)
            st.session_state['interpolation_results'] = {
                'interp_df': interp_df,
                'mm_results_df': mm_results_df,
                'mm_results': mm_results,  # 초기속도 탭용 (원본 v0 유지)
                'mm_fit_results': mm_fit_results,  # 정규화 기반 MM fit 결과
                'x_range_min': x_range_min,
                'x_range_max': x_range_max,
                'x_data_min': x_data_min,
                'x_data_max': x_data_max,
                'raw_data': raw_data,
                'v0_vs_concentration': {
                    'concentrations': norm_concentrations,  # 정규화 기반 농도
                    'v0_values': norm_v0_list  # 정규화 기반 v0
                },
                'experiment_type': experiment_type,
                'normalization_results': normalization_results,  # 정규화 결과 추가
                'uploaded_filename': uploaded_file.name if uploaded_file is not None else None  # 다운로드 파일명용
            }
            
            # 결과 적용 플래그 설정
            st.session_state['mm_data_ready'] = True
    
    # 결과 표시
    if 'interpolation_results' in st.session_state:
            results = st.session_state['interpolation_results']
            
            st.markdown("---")
            st.subheader("📊 Michaelis-Menten Model Results")
            
            # Display MM fitting results (varies by experiment condition)
            if 'mm_fit_results' in results and results['mm_fit_results']['fit_success']:
                mm_fit = results['mm_fit_results']
                exp_type = mm_fit.get('experiment_type', 'Substrate Concentration Variation (Standard MM)')
                
                if exp_type == "Substrate Concentration Variation (Standard MM)":
                    # 표준 MM 결과 표시 (Substrate는 μM 단위)
                    col1, col2, col3, col4 = st.columns(4)
                    with col1:
                        st.metric("Vmax", f"{mm_fit['Vmax']:.2f}" if mm_fit['Vmax'] is not None else "N/A")
                    with col2:
                        st.metric("Km (μM)", f"{mm_fit['Km']:.4f}" if mm_fit['Km'] is not None else "N/A")
                    with col3:
                        st.metric("kcat", f"{mm_fit['kcat']:.2f}" if mm_fit['kcat'] is not None else "N/A")
                    with col4:
                        st.metric("R²", f"{mm_fit['R_squared']:.4f}")
                    
                    st.info(f"**MM Equation:** {mm_fit['equation']}")
                else:
                    # Display Enzyme concentration variation results
                    st.warning("⚠️ **Fixed substrate + Enzyme Concentration Variation Experiment** (not standard MM)")
                    
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        slope = mm_fit.get('slope', None)
                        st.metric("Slope (v₀ vs [E])", f"{slope:.4f}" if slope is not None else "N/A",
                                 help="slope = kcat × [S] / (Km + [S])")
                    with col2:
                        intercept = mm_fit.get('intercept', None)
                        st.metric("Intercept", f"{intercept:.4f}" if intercept is not None else "N/A")
                    with col3:
                        st.metric("R²", f"{mm_fit['R_squared']:.4f}")
                    
                    st.info(f"**Linear Equation:** {mm_fit['equation']}")
                    st.info("""
                    📌 **Experiment Characteristics:**
                    - v has a **linear** relationship with [E]
                    
                    📌 **Available Parameters:**
                    - **Slope**: kcat × [S] / (Km + [S])
                    - If substrate concentration is very low: slope ≈ kcat/Km × [S]
                    
                    ❌ **Unavailable Parameters:**
                    - **Km**: Substrate concentration gradient required
                    - **Vmax**: Standard MM definition requires [E] fixed
                    - **kcat**: Cannot be determined alone (only kcat/Km possible)
                    """)
            elif 'mm_fit_results' in results:
                st.warning("⚠️ MM fitting failed or insufficient data")
            
            # 탭 구성 (st.tabs 대신 st.radio를 사용하여 상태 제어)
            exp_type = results.get('mm_fit_results', {}).get('experiment_type', 'Substrate Concentration Variation (Standard MM)')
            
            tab_titles = []
            if exp_type == "Substrate Concentration Variation (Standard MM)":
                tab_titles = ["📊 Experimental Results", "🔄 Normalization", "📊 v₀ vs [S] Fit", "📋 Data Table"]
            else:
                tab_titles = ["📊 Experimental Results", "🔄 Normalization", "📊 v₀ vs [E] Linear Fit", "📋 Data Table"]
            
            # 탭 상태 초기화
            if 'current_data_load_tab' not in st.session_state:
                st.session_state['current_data_load_tab'] = tab_titles[0]
            
            # 탭 메뉴 (라디오 버튼으로 구현하여 상태 제어)
            # 위쪽 여백 추가 (margin-top: 24px)
            st.markdown(
                """
                <style>
                div[data-testid="stRadio"] {
                    margin-top: 24px;
                }
                </style>
                """,
                unsafe_allow_html=True
            )
            
            st.markdown("---")

            selected_tab = st.radio(
                "Tabs",
                tab_titles,
                index=tab_titles.index(st.session_state.get('current_data_load_tab', tab_titles[0])) if st.session_state.get('current_data_load_tab') in tab_titles else 0,
                horizontal=True,
                label_visibility="collapsed",
                key="data_load_tab_radio"
            )
            
            # 선택된 탭 상태 업데이트
            st.session_state['current_data_load_tab'] = selected_tab
            
            # Tab 1: 실험결과 (점만 표시)
            if selected_tab == tab_titles[0]:
                # with tab_objects[0]: 대신 직접 렌더링
                st.subheader("Experimental Results")
                
                fig = go.Figure()
                colors = ['blue', 'red', 'orange', 'green', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']
                
                # 농도 순서대로 정렬
                conc_col = None
                for col in ['Concentration [μM]', 'Concentration [ug/mL]']:
                    if col in results['mm_results_df'].columns:
                        conc_col = col
                        break
                
                if conc_col:
                    conc_order = results['mm_results_df'].sort_values(conc_col)['Concentration'].tolist()
                else:
                    conc_order = results['mm_results_df']['Concentration'].tolist()
                
                x_data_min = results['x_data_min']
                x_data_max = results['x_data_max']
                
                for idx, conc_name in enumerate(conc_order):
                    color = colors[idx % len(colors)]
                    
                    # 범례에 표시할 농도 이름 (실험 타입에 따라 단위 변환)
                    # 숫자 추출
                    conc_match = re.search(r'(\d+\.?\d*)', conc_name)
                    if conc_match:
                        conc_value = float(conc_match.group(1))
                        # 실험 타입에 따라 단위 결정
                        if exp_type == "Substrate Concentration Variation (Standard MM)":
                            legend_name = f"{conc_value} μM"
                        else:  # Enzyme 농도 변화
                            legend_name = f"{conc_value} μg/mL"
                    else:
                        legend_name = conc_name
                    
                    # 원본 데이터 포인트만 표시 (점만)
                    if 'raw_data' in results and conc_name in results['raw_data']:
                        raw_conc_data = results['raw_data'][conc_name]
                        times_raw = raw_conc_data['time']
                        values_raw = raw_conc_data['value']
                        
                        fig.add_trace(go.Scatter(
                            x=times_raw,
                            y=values_raw,
                            mode='markers',
                            name=legend_name,
                            marker=dict(size=8, color=color, symbol='circle', line=dict(width=1, color='white')),
                            legendgroup=conc_name,
                            showlegend=True
                        ))
                        
                        # Error bars (SD가 있는 경우)
                        # Substrate 조건: SD 표시 안 함 (모두 0이므로)
                        # Enzyme 조건: SD가 0이 아닐 때만 표시
                        if raw_conc_data.get('SD') is not None:
                            sd_values = raw_conc_data['SD']
                            # 실험 타입에 따라 조건부 표시
                            if exp_type == "Enzyme Concentration Variation (Fixed substrate)":
                                # Enzyme 조건: SD가 0이 아닌 값이 하나라도 있으면 표시
                                if isinstance(sd_values, (list, np.ndarray)):
                                    has_nonzero_sd = np.any(np.array(sd_values) > 0)
                                else:
                                    has_nonzero_sd = sd_values > 0 if sd_values is not None else False
                                
                                if has_nonzero_sd:
                                    fig.add_trace(go.Scatter(
                                        x=times_raw,
                                        y=values_raw,
                                        error_y=dict(type='data', array=sd_values, visible=True),
                                        mode='markers',
                                        marker=dict(size=0, opacity=0),
                                        legendgroup=conc_name,
                                        showlegend=False
                                    ))
                            # Substrate 조건에서는 SD 표시 안 함
                
                fig.update_layout(
                    xaxis_title='Time (min)',
                    yaxis_title='RFU',
                    height=700,
                    template='plotly_white',
                    plot_bgcolor='rgba(0,0,0,0)',
                    paper_bgcolor='rgba(0,0,0,0)',
                    hovermode='x unified',
                    legend=dict(
                        orientation="v",
                        yanchor="bottom",
                        y=0.05,
                        xanchor="right",
                        x=0.99,
                        bgcolor="rgba(0,0,0,0)",
                        bordercolor="rgba(0,0,0,0)",
                        borderwidth=0,
                        font=dict(color="white")
                    )
                )
                
                # 원본 데이터 시간 범위로 제한
                fig.update_xaxes(range=[results['x_data_min'], results['x_data_max']])
                fig.update_yaxes(rangemode='tozero')
                
                st.plotly_chart(fig, use_container_width=True)
            
            # Tab 2: 정규화
            norm_tab_idx = 1
            if selected_tab == tab_titles[norm_tab_idx]:
                st.subheader("🔄 Normalization Results")
                
                if 'normalization_results' in results and results['normalization_results']:
                    norm_results = results['normalization_results']
                    
                    # 농도 순서 정렬
                    conc_col = None
                    for col in ['Concentration [μM]', 'Concentration [ug/mL]']:
                        if col in results['mm_results_df'].columns:
                            conc_col = col
                            break
                    
                    if conc_col:
                        conc_order = sorted(norm_results.keys(), 
                                          key=lambda x: norm_results[x]['concentration'])
                    else:
                        conc_order = list(norm_results.keys())
                    
                    # 농도 선택 (옆으로 넘기기)
                    if len(conc_order) > 0:
                        # session_state를 사용하여 선택된 농도 인덱스 저장
                        if 'normalization_selected_conc_idx' not in st.session_state:
                            st.session_state['normalization_selected_conc_idx'] = 0
                        
                        # 실험 타입에 따라 단위 변환 로직
                        exp_type = results.get('mm_fit_results', {}).get('experiment_type', 'Substrate Concentration Variation (Standard MM)')
                        
                        # 미리 포맷팅된 옵션 리스트 생성
                        formatted_options = []
                        for i in range(len(conc_order)):
                            conc_name = conc_order[i]
                            conc_value = norm_results[conc_name]['concentration']
                            if exp_type == "Substrate Concentration Variation (Standard MM)":
                                formatted_options.append(f"{conc_value} μM")
                            else:
                                formatted_options.append(f"{conc_value} μg/mL")
                        
                        # Selectbox 생성
                        selected_option = st.selectbox(
                            "농도 선택",
                            formatted_options,
                            index=0,
                            key="normalization_conc_select_box"
                        )
                        
                        # 선택된 옵션에 해당하는 농도 이름 찾기
                        conc_idx = formatted_options.index(selected_option)
                        selected_conc = conc_order[conc_idx]
                        norm_data = norm_results[selected_conc]
                        
                        # 정규화된 데이터 플롯
                        fig_norm = go.Figure()
                        
                        # 데이터 범위 계산
                        t_min = norm_data['times'].min()
                        t_max = norm_data['times'].max()
                        
                        # v0 계산 (정규화된 데이터에서 initial linear region의 기울기)
                        # 정규화된 데이터: F_linear(t) = k_obs * t
                        # 원본 데이터로 변환: v0 = k_obs * (Fmax - F0)
                        v0 = norm_data['k_obs'] * (norm_data['Fmax'] - norm_data['F0']) if norm_data['k_obs'] is not None else None
                        
                        # Exponential fit 곡선 (Full kinetics) - 주황색 실선
                        if norm_data['k_obs'] is not None and norm_data['k_obs'] > 0:
                            # X축을 데이터 범위로만 제한 (extrapolation 제거)
                            t_fit = np.linspace(t_min, t_max, 200)
                            F_max = 1.0  # 정규화된 값이므로 최종 F_max = 1.0
                            fit_curve = exponential_fit_simple(t_fit, F_max, norm_data['k_obs'])
                            
                            fig_norm.add_trace(go.Scatter(
                                x=t_fit,
                                y=fit_curve,
                                mode='lines',
                                name='Exponential increase (Full kinetics)',
                                line=dict(color='orange', width=2.5)
                            ))
                            
                            # Initial linear region - 파란색 점선
                            # t=0에서의 접선: F_linear(t) = k_obs * t (정규화된 데이터, F0=0)
                            initial_slope = norm_data['k_obs']  # 정규화된 데이터에서의 기울기
                            linear_curve = initial_slope * t_fit
                            
                            # v0 정보를 범례에 포함
                            v0_label = f"Initial linear region (v₀={v0:.2f} RFU/min)" if v0 is not None else "Initial linear region"
                            
                            fig_norm.add_trace(go.Scatter(
                                x=t_fit,
                                y=linear_curve,
                                mode='lines',
                                name=v0_label,
                                line=dict(color='lightblue', width=2.5, dash='dash')
                            ))
                            
                            # 구간별 세로선 표시
                            tau = norm_data['tau']
                            if tau is not None and not np.isinf(tau) and tau > 0:
                                # 초기 구간: t ≤ 0.1τ
                                t_initial = 0.1 * tau
                                # 지수 구간: 0.1τ ≤ t ≤ 3τ
                                t_exponential_start = 0.1 * tau
                                t_exponential_end = 3.0 * tau
                                # Plateau 구간: t ≥ 3τ
                                t_plateau = 3.0 * tau
                                
                                # 세로선 추가 (데이터 범위 내에 있는 경우만)
                                # 각 구간을 화살표와 텍스트로 표시 (배경 없이 흰색 글씨, y 위치를 다르게 설정)
                                if t_initial <= t_max:
                                    # 초기 구간 세로선
                                    fig_norm.add_vline(
                                        x=t_initial,
                                        line_dash="dash",
                                        line_color="orange"
                                    )
                                    # 화살표와 텍스트 annotation 추가 (y=1.05, 가장 아래)
                                    fig_norm.add_annotation(
                                        x=t_initial,
                                        y=1.05,
                                        text="초기 구간 (t ≤ 0.1τ)",
                                        showarrow=True,
                                        arrowhead=2,
                                        arrowsize=1,
                                        arrowwidth=2,
                                        arrowcolor="orange",
                                        ax=0,
                                        ay=-30,
                                        bgcolor="rgba(0,0,0,0)",
                                        bordercolor="rgba(0,0,0,0)",
                                        borderwidth=0,
                                        font=dict(size=11, color="white")
                                    )
                                
                                # 지수 구간: 0.1τ와 3τ 사이의 중간 지점에 annotation 추가
                                if t_exponential_end <= t_max and t_initial <= t_max:
                                    # 지수 구간 시작과 끝 세로선
                                    fig_norm.add_vline(
                                        x=t_exponential_start,
                                        line_dash="dash",
                                        line_color="purple"
                                    )
                                    fig_norm.add_vline(
                                        x=t_exponential_end,
                                        line_dash="dash",
                                        line_color="purple"
                                    )
                                    # 지수 구간 중간 지점 계산 (0.1τ와 3τ의 중간)
                                    t_exponential_mid = (t_exponential_start + t_exponential_end) / 2.0
                                    # 중간 지점이 데이터 범위 내에 있는 경우에만 annotation 추가
                                    if t_exponential_mid <= t_max:
                                        fig_norm.add_annotation(
                                            x=t_exponential_mid,
                                            y=1.10,
                                            text="지수 구간 (0.1τ < t < 3τ)",
                                            showarrow=True,
                                            arrowhead=2,
                                            arrowsize=1,
                                            arrowwidth=2,
                                            arrowcolor="purple",
                                            ax=0,
                                            ay=-30,
                                            bgcolor="rgba(0,0,0,0)",
                                            bordercolor="rgba(0,0,0,0)",
                                            borderwidth=0,
                                            font=dict(size=11, color="white")
                                        )
                                
                                if t_plateau <= t_max:
                                    # plateau 세로선 추가
                                    fig_norm.add_vline(
                                        x=t_plateau,
                                        line_dash="dash",
                                        line_color="brown"
                                    )
                                    # 화살표와 텍스트 annotation 추가 (y=1.15, 가장 위)
                                    # t_exponential_end와 같은 위치인지 확인하여 x 위치 조정
                                    if abs(t_plateau - t_exponential_end) < 0.001:
                                        # 같은 위치이면 x 방향으로 약간 이동하여 겹침 방지
                                        fig_norm.add_annotation(
                                            x=t_plateau,
                                            y=1.15,
                                            text="Plateau 구간 (t ≥ 3τ)",
                                            showarrow=True,
                                            arrowhead=2,
                                            arrowsize=1,
                                            arrowwidth=2,
                                            arrowcolor="brown",
                                            ax=20,
                                            ay=-30,
                                            bgcolor="rgba(0,0,0,0)",
                                            bordercolor="rgba(0,0,0,0)",
                                            borderwidth=0,
                                            font=dict(size=11, color="white")
                                        )
                                    else:
                                        fig_norm.add_annotation(
                                            x=t_plateau,
                                            y=1.15,
                                            text="Plateau 구간 (t ≥ 3τ)",
                                            showarrow=True,
                                            arrowhead=2,
                                            arrowsize=1,
                                            arrowwidth=2,
                                            arrowcolor="brown",
                                            ax=0,
                                            ay=-30,
                                            bgcolor="rgba(0,0,0,0)",
                                            bordercolor="rgba(0,0,0,0)",
                                            borderwidth=0,
                                            font=dict(size=11, color="white")
                                        )
                        
                        fig_norm.update_layout(
                            xaxis_title='Time (min)',
                            yaxis_title='Fluorescence intensity (a.u.)',
                            title='Enzyme-quenched peptide fluorescence kinetics',
                            height=600,
                            template='plotly_white',
                            plot_bgcolor='rgba(0,0,0,0)',
                            paper_bgcolor='rgba(0,0,0,0)',
                            hovermode='x unified',
                            # Y축 범위를 0-1로 고정 (정규화된 데이터)
                            # annotation이 잘 보이도록 범위 확장
                            yaxis=dict(range=[0, 1.20]),  # annotation 공간 확보
                            # X축 범위를 데이터 범위로 제한
                            xaxis=dict(range=[t_min, t_max]),
                            legend=dict(
                                orientation="v",
                                yanchor="bottom",
                                y=0.05,
                                xanchor="right",
                                x=0.99,
                                bgcolor="rgba(0,0,0,0)",
                                bordercolor="rgba(0,0,0,0)",
                                borderwidth=0
                            )
                        )
                        
                        st.plotly_chart(fig_norm, use_container_width=True)
                        
                        # y=1에 도달하는 시간까지만 표시하는 플롯 추가
                        if norm_data['k_obs'] is not None and norm_data['k_obs'] > 0:
                            tau_scaled = norm_data['tau']
                            if tau_scaled is not None and not np.isinf(tau_scaled) and tau_scaled > 0:
                                # 변수 재정의 (스코프 문제 해결)
                                F_max_scaled = 1.0  # 정규화된 값이므로 최종 F_max = 1.0
                                initial_slope_scaled = norm_data['k_obs']  # 정규화된 데이터에서의 기울기
                                v0_label_scaled = f"Initial linear region (v₀={v0:.2f} RFU/min)" if v0 is not None else "Initial linear region"
                                
                                # y=1에 도달하는 시간 찾기 (exponential 식 사용)
                                # F(t) = 1.0 * (1 - exp(-k_obs * t)) = target_value
                                # target_value = 0.99 (99%에 도달하는 시간)
                                target_value = 0.99
                                # 1 - exp(-k_obs * t) = target_value
                                # exp(-k_obs * t) = 1 - target_value
                                # -k_obs * t = ln(1 - target_value)
                                # t = -ln(1 - target_value) / k_obs
                                
                                normalized_values = norm_data['normalized_values']
                                times_norm = norm_data['times']
                                
                                # 실제 데이터에서 정규화된 값이 target_value 이상이 되는 첫 번째 시간 찾기
                                t_y1 = None
                                for i, val in enumerate(normalized_values):
                                    if val >= target_value:
                                        t_y1 = times_norm[i]
                                        break
                                
                                # 찾지 못하면 exponential 식으로 계산 (실제 데이터 범위를 넘어서도 계산)
                                if t_y1 is None:
                                    # t = -ln(1 - target_value) / k_obs
                                    # target_value = 0.99일 때: t = -ln(0.01) / k_obs ≈ 4.6 / k_obs ≈ 4.6 * tau
                                    t_y1 = -np.log(1 - target_value) / norm_data['k_obs']
                                
                                # t_display_max는 t_y1 사용 (t_max로 제한하지 않음 - 곡선 확장)
                                t_display_max = t_y1
                                
                                # 구간 계산
                                t_initial_scaled = 0.1 * tau_scaled
                                t_exponential_start_scaled = 0.1 * tau_scaled
                                t_exponential_end_scaled = 3.0 * tau_scaled
                                t_plateau_scaled = 3.0 * tau_scaled
                                t_exponential_mid_scaled = (t_exponential_start_scaled + t_exponential_end_scaled) / 2.0
                                
                                # y=1에 도달하는 시간까지만 플롯 생성 (원본 시간 스케일 유지)
                                fig_norm_scaled = go.Figure()
                                
                                # Exponential fit 곡선 (t_display_max까지만)
                                t_fit_scaled = np.linspace(t_min, t_display_max, 200)
                                fit_curve_scaled = exponential_fit_simple(t_fit_scaled, F_max_scaled, norm_data['k_obs'])
                                
                                fig_norm_scaled.add_trace(go.Scatter(
                                    x=t_fit_scaled,
                                    y=fit_curve_scaled,
                                    mode='lines',
                                    name='Exponential increase (Full kinetics)',
                                    line=dict(color='orange', width=2.5)
                                ))
                                
                                # Initial linear region
                                linear_curve_scaled = initial_slope_scaled * t_fit_scaled
                                
                                fig_norm_scaled.add_trace(go.Scatter(
                                    x=t_fit_scaled,
                                    y=linear_curve_scaled,
                                    mode='lines',
                                    name=v0_label_scaled,
                                    line=dict(color='lightblue', width=2.5, dash='dash')
                                ))
                                
                                # 구간별 세로선 표시
                                if t_initial_scaled <= t_display_max:
                                    fig_norm_scaled.add_vline(
                                        x=t_initial_scaled,
                                        line_dash="dash",
                                        line_color="orange"
                                    )
                                    fig_norm_scaled.add_annotation(
                                        x=t_initial_scaled,
                                        y=1.05,
                                        text="초기 구간 (t ≤ 0.1τ)",
                                        showarrow=True,
                                        arrowhead=2,
                                        arrowsize=1,
                                        arrowwidth=2,
                                        arrowcolor="orange",
                                        ax=0,
                                        ay=-30,
                                        bgcolor="rgba(0,0,0,0)",
                                        bordercolor="rgba(0,0,0,0)",
                                        borderwidth=0,
                                        font=dict(size=11, color="white")
                                    )
                                
                                if t_exponential_end_scaled <= t_display_max and t_initial_scaled <= t_display_max:
                                    fig_norm_scaled.add_vline(
                                        x=t_exponential_start_scaled,
                                        line_dash="dash",
                                        line_color="purple"
                                    )
                                    fig_norm_scaled.add_vline(
                                        x=t_exponential_end_scaled,
                                        line_dash="dash",
                                        line_color="purple"
                                    )
                                    if t_exponential_mid_scaled <= t_display_max:
                                        fig_norm_scaled.add_annotation(
                                            x=t_exponential_mid_scaled,
                                            y=1.10,
                                            text="지수 구간 (0.1τ < t < 3τ)",
                                            showarrow=True,
                                            arrowhead=2,
                                            arrowsize=1,
                                            arrowwidth=2,
                                            arrowcolor="purple",
                                            ax=0,
                                            ay=-30,
                                            bgcolor="rgba(0,0,0,0)",
                                            bordercolor="rgba(0,0,0,0)",
                                            borderwidth=0,
                                            font=dict(size=11, color="white")
                                        )
                                
                                if t_plateau_scaled <= t_display_max:
                                    fig_norm_scaled.add_vline(
                                        x=t_plateau_scaled,
                                        line_dash="dash",
                                        line_color="brown"
                                    )
                                    if abs(t_plateau_scaled - t_exponential_end_scaled) < 0.001:
                                        fig_norm_scaled.add_annotation(
                                            x=t_plateau_scaled,
                                            y=1.15,
                                            text="Plateau 구간 (t ≥ 3τ)",
                                            showarrow=True,
                                            arrowhead=2,
                                            arrowsize=1,
                                            arrowwidth=2,
                                            arrowcolor="brown",
                                            ax=20,
                                            ay=-30,
                                            bgcolor="rgba(0,0,0,0)",
                                            bordercolor="rgba(0,0,0,0)",
                                            borderwidth=0,
                                            font=dict(size=11, color="white")
                                        )
                                    else:
                                        fig_norm_scaled.add_annotation(
                                            x=t_plateau_scaled,
                                            y=1.15,
                                            text="Plateau 구간 (t ≥ 3τ)",
                                            showarrow=True,
                                            arrowhead=2,
                                            arrowsize=1,
                                            arrowwidth=2,
                                            arrowcolor="brown",
                                            ax=0,
                                            ay=-30,
                                            bgcolor="rgba(0,0,0,0)",
                                            bordercolor="rgba(0,0,0,0)",
                                            borderwidth=0,
                                            font=dict(size=11, color="white")
                                        )
                                
                                fig_norm_scaled.update_layout(
                                    xaxis_title='Time (min)',
                                    yaxis_title='Fluorescence intensity (a.u.)',
                                    title='Enzyme-quenched peptide fluorescence kinetics (up to y=1)',
                                    height=600,
                                    template='plotly_white',
                                    plot_bgcolor='rgba(0,0,0,0)',
                                    paper_bgcolor='rgba(0,0,0,0)',
                                    hovermode='x unified',
                                    yaxis=dict(range=[0, 1.20]),
                                    xaxis=dict(range=[t_min, t_display_max]),
                                    legend=dict(
                                        orientation="v",
                                        yanchor="bottom",
                                        y=0.05,
                                        xanchor="right",
                                        x=0.99,
                                        bgcolor="rgba(0,0,0,0)",
                                        bordercolor="rgba(0,0,0,0)",
                                        borderwidth=0
                                    )
                                )
                                
                                st.plotly_chart(fig_norm_scaled, use_container_width=True)
                        
                        # 정규화 방법 설명
                        with st.expander("📖 정규화에 사용된 식", expanded=False):
                            st.markdown("""
                            **정규화 과정:**
                            
                            1. **1차 정규화 (임시 정규화)**
                               - F₀ = 시간 0에서의 형광값 (F(t=0))
                               - Fmax = 최대 형광값 (max(F))
                               - 정규화: α_temp = (F - F₀) / (Fmax - F₀)
                            
                            2. **2차 정규화 (Exponential 피팅)**
                               - 정규화된 데이터에 대해 지수 함수 피팅
                               - **정규화된 데이터 식**: F_norm(t) = F_max · [1 - exp(-k_obs · t)]
                                 - F_norm(t): 정규화된 형광값 (0~1 범위)
                                 - F_max: 정규화된 최대값 (보통 1.0)
                                 - k_obs: 관찰된 반응 속도 상수 (분⁻¹)
                                 - t: 시간 (분)
                            
                            3. **원본 데이터로 변환**
                               - **원본 데이터 식**: F(t) = F₀ + (Fmax - F₀) · [1 - exp(-k_obs · t)]
                                 - F(t): 시간 t에서의 형광값
                                 - F₀: 초기 형광값
                                 - Fmax: 최대 형광값
                                 - k_obs: 관찰된 반응 속도 상수 (분⁻¹)
                                 - t: 시간 (분)
                            
                            4. **초기 속도 (v₀) 계산**
                               - v₀ = k_obs · (Fmax - F₀)
                               - 단위: RFU/min (형광 단위/분)
                            
                            5. **특성 시간 (τ)**
                               - τ = 1 / k_obs
                               - 반응이 63.2% 완료되는 시간
                            
                            **반복 정규화:**
                            - 위 과정을 최소 2번 반복하여 정규화를 개선
                            - 각 반복에서 피팅된 값을 역정규화하여 다음 반복에 사용
                            """)
                        
                        # 방정식 및 R² 테이블
                        st.subheader("Normalization Parameters")
                        # 실험 타입에 따라 농도 단위 결정
                        exp_type = results.get('mm_fit_results', {}).get('experiment_type', 'Substrate Concentration Variation (Standard MM)')
                        conc_value = norm_data['concentration']
                        if exp_type == "Substrate Concentration Variation (Standard MM)":
                            conc_display = f"{conc_value} μM"
                        else:  # Enzyme 농도 변화
                            conc_display = f"{conc_value} μg/mL"
                        
                        # 구간 경계 시간 계산
                        tau = norm_data.get('tau', None)
                        if tau is not None and not np.isinf(tau) and tau > 0:
                            t_initial_exponential = 0.1 * tau  # 초기-exponential 경계
                            t_exponential_plateau = 3.0 * tau  # exponential-plateau 경계
                            t_initial_exponential_str = f"{t_initial_exponential:.4f}"
                            t_exponential_plateau_str = f"{t_exponential_plateau:.4f}"
                        else:
                            t_initial_exponential_str = "N/A"
                            t_exponential_plateau_str = "N/A"
                        
                        param_data = {
                            '농도': [conc_display],
                            'F₀': [f"{norm_data['F0']:.4f}"],
                            'F_max': [f"{norm_data['Fmax']:.4f}"],
                            'k_obs': [f"{norm_data['k_obs']:.4f}" if norm_data['k_obs'] is not None else "N/A"],
                            'τ (1/k_obs)': [f"{norm_data['tau']:.4f}" if norm_data['tau'] is not None and not np.isinf(norm_data['tau']) else "N/A"],
                            '초기→지수 (t=0.1τ)': [t_initial_exponential_str],
                            '지수→Plateau (t=3τ)': [t_exponential_plateau_str],
                            'v₀ (RFU/min)': [f"{v0:.2f}" if v0 is not None else "N/A"],
                            'R²': [f"{norm_data['R_squared']:.4f}"],
                            '방정식': [norm_data['equation']]
                        }
                        param_df = pd.DataFrame(param_data)
                        st.dataframe(param_df, use_container_width=True, hide_index=True)
                        
                        # 모든 농도 요약 테이블
                        st.subheader("모든 농도 정규화 요약")
                        summary_data = []
                        for conc_name in conc_order:
                            n_data = norm_results[conc_name]
                            # v0 계산
                            v0_conc = n_data['k_obs'] * (n_data['Fmax'] - n_data['F0']) if n_data['k_obs'] is not None else None
                            # 농도 표시 (실험 타입에 따라 단위 변환)
                            conc_value = n_data['concentration']
                            if exp_type == "Substrate Concentration Variation (Standard MM)":
                                conc_display = f"{conc_value} μM"
                            else:  # Enzyme 농도 변화
                                conc_display = f"{conc_value} μg/mL"
                            
                            # α (절단비율) 계산: normalized_values = (F - F0) / (Fmax - F0)
                            normalized_values = n_data.get('normalized_values', [])
                            if len(normalized_values) > 0:
                                alpha_min = np.min(normalized_values)
                                alpha_max = np.max(normalized_values)
                                alpha_mean = np.mean(normalized_values)
                                alpha_range_str = f"{alpha_min:.3f}-{alpha_max:.3f}"
                                alpha_mean_str = f"{alpha_mean:.3f}"
                            else:
                                alpha_range_str = "N/A"
                                alpha_mean_str = "N/A"
                            
                            # 구간 경계 시간 계산
                            tau_conc = n_data.get('tau', None)
                            if tau_conc is not None and not np.isinf(tau_conc) and tau_conc > 0:
                                t_initial_exponential_conc = 0.1 * tau_conc  # 초기-exponential 경계
                                t_exponential_plateau_conc = 3.0 * tau_conc  # exponential-plateau 경계
                                t_initial_exponential_str_conc = f"{t_initial_exponential_conc:.4f}"
                                t_exponential_plateau_str_conc = f"{t_exponential_plateau_conc:.4f}"
                            else:
                                t_initial_exponential_str_conc = "N/A"
                                t_exponential_plateau_str_conc = "N/A"
                            
                            summary_data.append({
                                '농도': conc_display,
                                'F₀': f"{n_data['F0']:.4f}",
                                'F_max': f"{n_data['Fmax']:.4f}",
                                'k_obs': f"{n_data['k_obs']:.4f}" if n_data['k_obs'] is not None else "N/A",
                                'τ': f"{n_data['tau']:.4f}" if n_data['tau'] is not None and not np.isinf(n_data['tau']) else "N/A",
                                '초기→지수 (t=0.1τ)': t_initial_exponential_str_conc,
                                '지수→Plateau (t=3τ)': t_exponential_plateau_str_conc,
                                'v₀ (RFU/min)': f"{v0_conc:.2f}" if v0_conc is not None else "N/A",
                                'α 범위': alpha_range_str,
                                'α 평균': alpha_mean_str,
                                'R²': f"{n_data['R_squared']:.4f}",
                                '방정식': n_data['equation'][:50] + "..." if len(n_data['equation']) > 50 else n_data['equation']
                            })
                        summary_df = pd.DataFrame(summary_data)
                        st.dataframe(summary_df, use_container_width=True, hide_index=True)
                        
                else:
                    st.info("No normalization results. Please run the 'Run Michaelis-Menten Model' button first.")
            
            # Tab 2: v₀ vs 농도 그래프 (실험 조건에 따라 다름)
            v0_tab_idx = 2 if exp_type == "Substrate Concentration Variation (Standard MM)" else 2
            if selected_tab == tab_titles[v0_tab_idx]:
                if 'v0_vs_concentration' in results and 'mm_fit_results' in results:
                    v0_data = results['v0_vs_concentration']
                    mm_fit = results['mm_fit_results']
                    exp_type = mm_fit.get('experiment_type', 'Substrate Concentration Variation (Standard MM)')
                    
                    fig_v0 = go.Figure()
                    
                    # 실험 데이터 포인트
                    fig_v0.add_trace(go.Scatter(
                        x=v0_data['concentrations'],
                        y=v0_data['v0_values'],
                        mode='markers',
                        name='Experimental v₀',
                        marker=dict(size=10, color='red', line=dict(width=2, color='black'))
                    ))
                    
                    if exp_type == "Substrate Concentration Variation (Standard MM)":
                        st.subheader("v₀ vs [S] Michaelis-Menten Fit")
                        
                        # 실험 데이터 테이블 (expander)
                        with st.expander("📋 Preview Experimental Data", expanded=False):
                            # 농도와 v0 데이터를 테이블로 표시
                            exp_data = {
                                '[S] (μM)': v0_data['concentrations'],
                                'v₀ (RFU/min)': v0_data['v0_values']
                            }
                            exp_df = pd.DataFrame(exp_data)
                            # 농도 순서대로 정렬
                            exp_df = exp_df.sort_values('[S] (μM)')
                            exp_df = exp_df.reset_index(drop=True)
                            st.dataframe(exp_df, use_container_width=True, hide_index=True)
                        
                        # MM 피팅 곡선
                        if mm_fit['fit_success'] and mm_fit['Vmax'] is not None and mm_fit['Km'] is not None:
                            conc_min = min(v0_data['concentrations'])
                            conc_max = max(v0_data['concentrations'])
                            conc_range = np.linspace(conc_min * 0.5, conc_max * 1.5, 200)
                            v0_fitted = michaelis_menten_calibration(conc_range, mm_fit['Vmax'], mm_fit['Km'])
                            
                            fig_v0.add_trace(go.Scatter(
                                x=conc_range,
                                y=v0_fitted,
                                mode='lines',
                                name=f'{mm_fit["equation"]}',
                                line=dict(width=2.5, color='blue')
                            ))
                            
                            # 통계 정보 (Substrate는 μM 단위)
                            stats_text = f"Vmax = {mm_fit['Vmax']:.2f}<br>"
                            stats_text += f"Km = {mm_fit['Km']:.4f} μM<br>"
                            stats_text += f"R² = {mm_fit['R_squared']:.4f}"
                            
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
                            title='Initial Velocity (v₀) vs Substrate Concentration [S]',
                            xaxis_title='[S] (μM)',
                            yaxis_title='Initial Velocity v₀ (Fluorescence Units / Time)',
                            template='plotly_white',
                            height=600,
                            hovermode='x unified'
                        )
                    else:
                        st.subheader("v₀ vs [E] Linear Fit (Constant Substrate)")
                        st.warning("⚠️ This is not the standard Michaelis-Menten model. v is linearly related to [E].")
                        
                        # 실험 데이터 테이블 (expander)
                        with st.expander("📋 Preview Experimental Data", expanded=False):
                            # 농도와 v0 데이터를 테이블로 표시
                            exp_data = {
                                '[E] (μg/mL)': v0_data['concentrations'],
                                'v₀ (RFU/min)': v0_data['v0_values']
                            }
                            exp_df = pd.DataFrame(exp_data)
                            # 농도 순서대로 정렬
                            exp_df = exp_df.sort_values('[E] (μg/mL)')
                            exp_df = exp_df.reset_index(drop=True)
                            st.dataframe(exp_df, use_container_width=True, hide_index=True)
                        
                        # 선형 피팅 곡선
                        if mm_fit['fit_success'] and mm_fit.get('slope') is not None:
                            conc_min = min(v0_data['concentrations'])
                            conc_max = max(v0_data['concentrations'])
                            conc_range = np.linspace(conc_min * 0.5, conc_max * 1.5, 200)
                            slope = mm_fit['slope']
                            intercept = mm_fit.get('intercept', 0)
                            v0_fitted = slope * conc_range + intercept
                            
                            fig_v0.add_trace(go.Scatter(
                                x=conc_range,
                                y=v0_fitted,
                                mode='lines',
                                name=f'Linear Fit: {mm_fit["equation"]}',
                                line=dict(width=2.5, color='blue', dash='dash')
                            ))
                            
                            # 통계 정보
                            stats_text = f"Slope = {slope:.4f}<br>"
                            stats_text += f"Intercept = {intercept:.4f}<br>"
                            stats_text += f"R² = {mm_fit['R_squared']:.4f}<br>"
                            stats_text += "<br><b>⚠️ Cannot calculate Km</b>"
                            
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
                            title='Initial Velocity (v₀) vs Enzyme Concentration [E] (Constant Substrate)',
                            xaxis_title='[E] (μg/mL)',
                            yaxis_title='Initial Velocity v₀ (Fluorescence Units / Time)',
                            template='plotly_white',
                            height=600,
                            hovermode='x unified'
                        )
                    
                    st.plotly_chart(fig_v0, use_container_width=True)
                else:
                    st.warning("No v₀ vs concentration data available.")
            
            # 마지막 탭: 데이터 테이블
            data_tab_idx = 3 if exp_type == "Substrate Concentration Variation (Standard MM)" else 3
            if selected_tab == tab_titles[data_tab_idx]:
                st.subheader("Detailed Parameters")
                
                # 상세 파라미터 테이블 (정규화 기반 v0 사용)
                # 정규화 결과에서 v0 가져오기
                if 'normalization_results' in results and results['normalization_results']:
                    norm_results = results['normalization_results']
                    
                    # 정규화 기반 v0으로 업데이트된 데이터프레임 생성
                    detail_data = []
                    for conc_name in sorted(norm_results.keys(), key=lambda x: norm_results[x]['concentration']):
                        norm_data = norm_results[conc_name]
                        conc_value = norm_data['concentration']
                        
                        # 정규화 기반 v0 계산
                        v0_norm = norm_data.get('v0', 0)
                        if v0_norm == 0 and norm_data.get('k_obs') is not None:
                            v0_norm = norm_data['k_obs'] * (norm_data['Fmax'] - norm_data['F0'])
                        
                        # 실험 타입에 따라 농도 단위 결정
                        if exp_type == "Substrate Concentration Variation (Standard MM)":
                            conc_col_name = 'Concentration [μM]'
                        else:
                            conc_col_name = 'Concentration [ug/mL]'
                        
                        # mm_results에서 해당 농도 찾기
                        mm_data = None
                        for mm_conc_name, mm_params in results.get('mm_results', {}).items():
                            if mm_params.get('concentration') == conc_value:
                                mm_data = mm_params
                                break
                        
                        row_data = {
                            conc_col_name: conc_value,
                            'v0': v0_norm,  # 정규화 기반 v0
                            'F0': norm_data['F0'],
                            'Fmax': norm_data['Fmax'],
                            'R_squared': norm_data['R_squared'],
                            'k_obs': norm_data.get('k_obs', None),
                            'τ': norm_data.get('tau', None),
                            '방정식': norm_data['equation']
                        }
                        
                        detail_data.append(row_data)
                    
                    detail_df = pd.DataFrame(detail_data)
                    st.dataframe(detail_df, use_container_width=True, hide_index=True)
                else:
                    # 정규화 결과가 없으면 기존 방식 사용
                    detail_cols = ['Concentration [μM]', 'Concentration [ug/mL]', 'v0', 'F0', 'Fmax', 'R_squared', 'Equation']
                    available_cols = [col for col in detail_cols if col in results['mm_results_df'].columns]
                    st.dataframe(results['mm_results_df'][available_cols], use_container_width=True, hide_index=True)
                
                # MM Fit 결과 표시
                st.markdown("---")
                st.subheader("MM Fit Results")
                if 'mm_fit_results' in results and results['mm_fit_results'].get('fit_success'):
                    mm_fit = results['mm_fit_results']
                    if exp_type == "Substrate Concentration Variation (Standard MM)":
                        mm_fit_data = {
                            'Parameter': ['Vmax', 'Km (μM)', 'kcat', 'R²'],
                            'Value': [
                                f"{mm_fit['Vmax']:.2f}" if mm_fit['Vmax'] is not None else "N/A",
                                f"{mm_fit['Km']:.4f}" if mm_fit['Km'] is not None else "N/A",
                                f"{mm_fit['kcat']:.2f}" if mm_fit['kcat'] is not None else "N/A",
                                f"{mm_fit['R_squared']:.4f}"
                            ]
                        }
                    else:
                        mm_fit_data = {
                            'Parameter': ['Slope', 'Intercept', 'R²'],
                            'Value': [
                                f"{mm_fit.get('slope', 0):.4f}" if mm_fit.get('slope') is not None else "N/A",
                                f"{mm_fit.get('intercept', 0):.4f}" if mm_fit.get('intercept') is not None else "N/A",
                                f"{mm_fit['R_squared']:.4f}"
                            ]
                        }
                    mm_fit_df = pd.DataFrame(mm_fit_data)
                    st.dataframe(mm_fit_df, use_container_width=True, hide_index=True)
                else:
                    st.warning("No MM Fit results available.")
                
                # 파일 다운로드 버튼
                st.markdown("---")
                col1, col2 = st.columns(2)
                
                # 다운로드용 데이터프레임 결정
                if 'normalization_results' in results and results['normalization_results']:
                    download_df = detail_df  # 정규화 결과 포함
                else:
                    detail_cols = ['Concentration [μM]', 'Concentration [ug/mL]', 'v0', 'F0', 'Fmax', 'R_squared', 'linear_fraction', 'Equation']
                    available_cols = [col for col in detail_cols if col in results['mm_results_df'].columns]
                    download_df = results['mm_results_df'][available_cols]
                
                # XLSX 다운로드 버튼 및 자동 저장
                with col1:
                    try:
                        from io import BytesIO
                        output = BytesIO()
                        with pd.ExcelWriter(output, engine='openpyxl') as writer:
                            download_df.to_excel(writer, sheet_name='MM Results', index=False)
                            
                            # Michaelis-Menten Curves 시트: Concentration 컬럼명 수정
                            interp_df_copy = results['interp_df'].copy()
                            # 'Concentration' 컬럼이 있으면 제거 (conc_unit_col 사용)
                            if 'Concentration' in interp_df_copy.columns:
                                interp_df_copy = interp_df_copy.drop(columns=['Concentration'])
                            
                            interp_df_copy.to_excel(writer, sheet_name='Michaelis-Menten Curves', index=False)
                            
                            # 정규화 결과 시트 추가
                            if 'normalization_results' in results and results['normalization_results']:
                                norm_results = results['normalization_results']
                                norm_summary_data = []
                                for conc_name in sorted(norm_results.keys(), key=lambda x: norm_results[x]['concentration']):
                                    n_data = norm_results[conc_name]
                                    conc_value = n_data['concentration']
                                    v0_conc = n_data.get('v0', 0)
                                    if v0_conc == 0 and n_data.get('k_obs') is not None:
                                        v0_conc = n_data['k_obs'] * (n_data['Fmax'] - n_data['F0'])
                                    
                                    if exp_type == "Substrate Concentration Variation (Standard MM)":
                                        conc_display = f"{conc_value} μM"
                                    else:
                                        conc_display = f"{conc_value} μg/mL"
                                    
                                    norm_summary_data.append({
                                        'Concentration': conc_display,
                                        'F₀': n_data['F0'],
                                        'F_max': n_data['Fmax'],
                                        'k_obs': n_data.get('k_obs', None),
                                        'τ': n_data.get('tau', None),
                                        'v₀ (RFU/min)': v0_conc,
                                        'R²': n_data['R_squared'],
                                        'Equation': n_data['equation']
                                    })
                                
                                if norm_summary_data:
                                    norm_summary_df = pd.DataFrame(norm_summary_data)
                                    norm_summary_df.to_excel(writer, sheet_name='Normalization Results', index=False)
                            
                            # MM Fit 결과 시트 추가
                            if 'mm_fit_results' in results and results['mm_fit_results'].get('fit_success'):
                                mm_fit = results['mm_fit_results']
                                if exp_type == "Substrate Concentration Variation (Standard MM)":
                                    mm_fit_data = {
                                        'Parameter': ['Vmax', 'Km (μM)', 'kcat', 'R²', 'Equation'],
                                        'Value': [
                                            mm_fit['Vmax'] if mm_fit['Vmax'] is not None else "N/A",
                                            mm_fit['Km'] if mm_fit['Km'] is not None else "N/A",
                                            mm_fit['kcat'] if mm_fit['kcat'] is not None else "N/A",
                                            mm_fit['R_squared'],
                                            mm_fit.get('equation', 'N/A')
                                        ]
                                    }
                                else:
                                    mm_fit_data = {
                                        'Parameter': ['Slope', 'Intercept', 'R²', 'Equation'],
                                        'Value': [
                                            mm_fit.get('slope', None),
                                            mm_fit.get('intercept', None),
                                            mm_fit['R_squared'],
                                            mm_fit.get('equation', 'N/A')
                                        ]
                                    }
                                mm_fit_df = pd.DataFrame(mm_fit_data)
                                mm_fit_df.to_excel(writer, sheet_name='MM Fit Results', index=False)
                        
                        output.seek(0)
                        xlsx_data = output.getvalue()
                        
                        # 다운로드 파일명: 업로드 파일이 raw_* 이면 result_*.xlsx 로 저장
                        uploaded_filename = results.get('uploaded_filename') or ''
                        if uploaded_filename.startswith('raw_'):
                            base = 'result_' + uploaded_filename[4:].rsplit('.', 1)[0]
                            xlsx_download_name = base + '.xlsx'
                        else:
                            xlsx_download_name = 'Michaelis-Menten_calibration_results.xlsx'
                        
                        # XLSX 파일 자동 저장 (Analysis 모드에서 자동 로드용)
                        try:
                            with open('Michaelis-Menten_calibration_results.xlsx', 'wb') as f:
                                f.write(xlsx_data)
                        except Exception as save_err:
                            st.sidebar.warning(f"⚠️ XLSX auto-save failed: {save_err}")
                        
                        st.download_button(
                            label="📥 Download All Results (XLSX)",
                            data=xlsx_data,
                            file_name=xlsx_download_name,
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                            use_container_width=True,
                            help="Full Excel file including MM Results, Normalization Results, MM Fit Results, and Michaelis-Menten Curves sheets."
                        )
                    except Exception as e:
                        st.warning(f"Error preparing XLSX download: {e}")
                
                # 모든 분석 이미지 ZIP 다운로드 (투명 PNG)
                with col2:
                    try:
                        import zipfile
                        from io import BytesIO
                        fig_list = _build_all_export_figures(results)
                        zip_buffer = BytesIO()
                        with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zf:
                            for name, fig in fig_list:
                                try:
                                    img_bytes = fig.to_image(format="png", scale=2, engine="kaleido")
                                    zf.writestr(f"{name}.png", img_bytes)
                                except Exception as img_err:
                                    err_msg = str(img_err)
                                    if "kaleido" in err_msg.lower() or "chrome" in err_msg.lower():
                                        err_msg += " → Chrome 설치 불필요. 터미널에서 'pip install -U kaleido' 실행 후 Streamlit 앱을 재시작해주세요. (kaleido 0.2+는 Chrome 없이 동작합니다.)"
                                    st.warning(f"Image generation failed ({name}): {err_msg}")
                        zip_buffer.seek(0)
                        zip_bytes = zip_buffer.getvalue()
                        if zip_bytes:
                            st.download_button(
                                label="📥 모든 이미지 저장 (ZIP)",
                                data=zip_bytes,
                                file_name="analysis_figures.zip",
                                mime="application/zip",
                                use_container_width=True,
                                help="Experimental Results, Normalization, Linear fit 등 분석에서 생성된 모든 그래프를 PNG로 저장한 ZIP 파일입니다."
                            )
                        else:
                            st.caption("이미지 생성에 실패했습니다. kaleido 설치 후 Streamlit 앱을 재시작해주세요. (pip install -U kaleido)")
                    except Exception as zip_err:
                        st.warning(f"Error preparing image ZIP: {zip_err}")

