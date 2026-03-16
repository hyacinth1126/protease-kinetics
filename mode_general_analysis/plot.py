#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Hydrogel FRET Advanced Kinetic Analysis - Visualization Tools
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from typing import List
from .analysis import ModelResults


class Visualizer:
    """Visualization tools for model comparison"""
    
    @staticmethod
    def _detect_conc_col(df: pd.DataFrame) -> str:
        """Detect concentration column name robustly."""
        if 'conc_col_name' in df.columns:
            name = df['conc_col_name'].iloc[0]
            if isinstance(name, str) and name in df.columns:
                return name
        # Fallback candidates
        candidates = [c for c in df.columns
                      if (c.endswith('uM') or c.endswith('nM') or '_uM' in c or '_nM' in c)
                      and 'time' not in c.lower()
                      and 'alpha' not in c.lower()
                      and 'region' not in c.lower()
                      and 'FL_' not in c and 'intensity' not in c.lower()]
        if candidates:
            return candidates[0]
        # Last resort
        return 'enzyme_ugml'
    
    @staticmethod
    def plot_raw_data(df: pd.DataFrame, conc_unit: str = 'μg/mL', time_label: str = '시간 (초)', use_lines: bool = False, 
                     enzyme_name: str = 'enzyme', substrate_name: str = 'substrate'):
        """Plot raw fluorescence data with exponential fits and asymptotes
        
        Args:
            df: DataFrame with time_s, FL_intensity columns
            conc_unit: Concentration unit for labels
            time_label: Time axis label
            use_lines: If True, use lines instead of markers (for fitted curves)
            enzyme_name: Custom name for enzyme (default: 'enzyme')
            substrate_name: Custom name for substrate (default: 'substrate')
        """
        fig = go.Figure()
        
        colors = px.colors.qualitative.Set1
        
        conc_col = Visualizer._detect_conc_col(df)
        base = conc_col.split('_')[0].lower()
        entity_name = substrate_name if base.startswith('pep') or base.startswith('sub') else enzyme_name
        for idx, conc in enumerate(sorted(df[conc_col].unique())):
            subset = df[df[conc_col] == conc]
            color = colors[idx % len(colors)]
            
            # Plot experimental data (markers or lines)
            plot_mode = 'lines' if use_lines else 'markers'
            marker_dict = dict(size=8, color=color) if not use_lines else None
            line_dict = dict(color=color, width=2.5) if use_lines else None
            
            fig.add_trace(go.Scatter(
                x=subset['time_s'],
                y=subset['FL_intensity'],
                mode=plot_mode,
                name=f'{entity_name} {conc} {conc_unit}',
                marker=marker_dict if marker_dict else None,
                line=line_dict if line_dict else None,
                error_y=dict(type='data', array=subset['SD'], visible=True) if ('SD' in subset.columns and not use_lines) else None
            ))
            
            # Plot exponential fit (if available)
            if 't_fit' in subset.columns and 'F_fit' in subset.columns:
                t_fit = subset['t_fit'].iloc[0]
                F_fit = subset['F_fit'].iloc[0]
                
                fig.add_trace(go.Scatter(
                    x=t_fit,
                    y=F_fit,
                    mode='lines',
                    name=f'{conc} {conc_unit} (지수 피팅)',
                    line=dict(color=color, dash='dash', width=2),
                    showlegend=True
                ))
                
                # Add asymptote line (Fmax)
                Fmax = subset['Fmax'].iloc[0]
                fig.add_trace(go.Scatter(
                    x=[t_fit[0], t_fit[-1]],
                    y=[Fmax, Fmax],
                    mode='lines',
                    name=f'{conc} {conc_unit} (점근선 Fmax)',
                    line=dict(color=color, dash='dot', width=1.5),
                    showlegend=True,
                    hovertemplate=f'Fmax = {Fmax:.1f}<extra></extra>'
                ))
        
        fig.update_layout(
            title='Raw data: 시간-형광 그래프',
            xaxis_title=time_label,
            yaxis_title='형광 강도 (RFU)',
            template='plotly_white',
            hovermode='x unified',
            height=600
        )
        
        return fig
    
    @staticmethod
    def plot_normalized_data(df: pd.DataFrame, conc_unit: str = 'μg/mL', time_label: str = '시간 (초)', use_lines: bool = False,
                            enzyme_name: str = 'enzyme', substrate_name: str = 'substrate', experiment_type: str = None,
                            points_every_minute: bool = False):
        """Plot normalized data (fraction cleaved)
        
        Args:
            df: DataFrame with time_s, alpha columns
            conc_unit: Concentration unit for labels
            time_label: Time axis label
            use_lines: If True, use lines instead of markers (for fitted curves)
            enzyme_name: Custom name for enzyme (default: 'enzyme')
            substrate_name: Custom name for substrate (default: 'substrate')
            experiment_type: Experiment type ("Substrate 농도 변화 (표준 MM)" or "Enzyme 농도 변화")
            points_every_minute: If True, resample to integer minutes and plot as points + line segments (꺾은선)
        """
        fig = go.Figure()
        
        colors = px.colors.qualitative.Set1
        conc_col = Visualizer._detect_conc_col(df)
        base = conc_col.split('_')[0].lower()
        
        # 실험 타입에 따라 entity_name 결정
        if experiment_type == "Substrate 농도 변화 (표준 MM)":
            entity_name = substrate_name
        elif experiment_type == "Enzyme 농도 변화" or experiment_type == "Enzyme 농도 변화 (Substrate 고정)":
            entity_name = enzyme_name
        else:
            # 실험 타입이 없으면 컬럼 이름으로 판단
            entity_name = substrate_name if base.startswith('pep') or base.startswith('sub') else enzyme_name
        
        concentrations = sorted(df[conc_col].unique())
        
        for idx, conc in enumerate(concentrations):
            subset = df[df[conc_col] == conc].sort_values('time_s')
            color = colors[idx % len(colors)]
            t = subset['time_s'].values
            alpha = subset['alpha'].values
            
            if points_every_minute and len(t) > 0:
                # 분 단위 점만 사용, 직선으로 연결 (꺾은선, 스무스 곡선 아님)
                t_max = float(np.nanmax(t))
                if t_max <= 120:
                    minute_grid = np.arange(0, int(np.ceil(t_max)) + 1, 1, dtype=float)
                else:
                    n_minutes = int(np.ceil(t_max / 60))
                    minute_grid = np.arange(0, n_minutes + 1, 1) * 60.0
                alpha_at_minutes = np.interp(minute_grid, t, alpha)
                x_plot = minute_grid
                y_plot = alpha_at_minutes
                # 선: 꺾은선만
                fig.add_trace(go.Scatter(
                    x=x_plot,
                    y=y_plot,
                    mode='lines',
                    name=f'{entity_name} {conc} {conc_unit}',
                    line=dict(color=color, width=2, shape='linear'),
                    legendgroup=f'group{idx}'
                ))
                # 점: 분 단위 마커 (나중에 그려서 선 위에 보이게)
                fig.add_trace(go.Scatter(
                    x=x_plot,
                    y=y_plot,
                    mode='markers',
                    name=f'{entity_name} {conc} {conc_unit}',
                    marker=dict(
                        size=10,
                        color=color,
                        line=dict(width=2, color='white'),
                        symbol='circle'
                    ),
                    legendgroup=f'group{idx}',
                    showlegend=False
                ))
            else:
                x_plot = t
                y_plot = alpha
                plot_mode = 'lines' if use_lines else 'markers'
                marker_dict = dict(size=8, color=color) if not use_lines else None
                line_dict = dict(color=color, width=2.5, shape='linear') if use_lines else None
                fig.add_trace(go.Scatter(
                    x=x_plot,
                    y=y_plot,
                    mode=plot_mode,
                    name=f'{entity_name} {conc} {conc_unit}',
                    marker=marker_dict if marker_dict is not None else None,
                    line=line_dict if line_dict is not None else None,
                    legendgroup=f'group{idx}'
                ))
        
        fig.update_layout(
            title='정규화 데이터: 절단 비율 α(t)',
            xaxis_title=time_label,
            yaxis_title='절단 비율 α',
            template='plotly_white',
            hovermode='x unified',
            height=600
        )
        
        return fig
    
    @staticmethod
    def plot_time_concentration_heatmap(df: pd.DataFrame, conc_unit: str = 'μg/mL', time_label: str = '시간 (초)',
                                       enzyme_name: str = 'enzyme', substrate_name: str = 'substrate'):
        """Plot time-concentration heatmap for normalized data
        
        Args:
            df: DataFrame with time_s, alpha columns
            conc_unit: Concentration unit for labels
            time_label: Time axis label
            enzyme_name: Custom name for enzyme (default: 'enzyme')
            substrate_name: Custom name for substrate (default: 'substrate')
        """
        conc_col = Visualizer._detect_conc_col(df)
        base = conc_col.split('_')[0].lower()
        entity_name = substrate_name if base.startswith('pep') or base.startswith('sub') else enzyme_name
        
        # Get unique time points and concentrations
        times = sorted(df['time_s'].unique())
        concentrations = sorted(df[conc_col].unique())
        
        # Create a matrix for heatmap: rows = concentrations, cols = times
        heatmap_data = []
        for conc in concentrations:
            row = []
            subset = df[df[conc_col] == conc]
            for time in times:
                time_subset = subset[subset['time_s'] == time]
                if len(time_subset) > 0:
                    # Use mean alpha if multiple values exist for same time
                    alpha_val = time_subset['alpha'].mean()
                else:
                    # Interpolate if missing
                    alpha_val = np.nan
                row.append(alpha_val)
            heatmap_data.append(row)
        
        # Create heatmap
        fig = go.Figure(data=go.Heatmap(
            z=heatmap_data,
            x=times,
            y=[f'{entity_name} {c} {conc_unit}' for c in concentrations],
            colorscale='Viridis',
            colorbar=dict(title='절단 비율 α'),
            hovertemplate='시간: %{x:.2f}<br>농도: %{y}<br>α: %{z:.3f}<extra></extra>'
        ))
        
        fig.update_layout(
            title='정규화 데이터: 시간-농도 그래프',
            xaxis_title=time_label,
            yaxis_title=f'농도 ({conc_unit})',
            template='plotly_white',
            height=400 + len(concentrations) * 30,  # Dynamic height based on number of concentrations
            yaxis=dict(autorange='reversed')  # Reverse y-axis so highest concentration is at top
        )
        
        return fig
    
    @staticmethod
    def plot_model_fits(df: pd.DataFrame, results: List[ModelResults], 
                       conc_unit: str = 'μg/mL', time_label: str = '시간 (초)',
                       enzyme_name: str = 'enzyme', substrate_name: str = 'substrate'):
        """Plot all model fits together
        
        Args:
            df: DataFrame with time_s, alpha columns
            results: List of ModelResults objects
            conc_unit: Concentration unit for labels
            time_label: Time axis label
            enzyme_name: Custom name for enzyme (default: 'enzyme')
            substrate_name: Custom name for substrate (default: 'substrate')
        """
        fig = make_subplots(
            rows=2, cols=1,
            subplot_titles=('모델 피팅', '잔차'),
            vertical_spacing=0.15,
            row_heights=[0.7, 0.3]
        )
        
        colors = px.colors.qualitative.Set1
        
        # Plot experimental data
        conc_col = Visualizer._detect_conc_col(df)
        base = conc_col.split('_')[0].lower()
        entity_name = substrate_name if base.startswith('pep') or base.startswith('sub') else enzyme_name
        for idx, conc in enumerate(sorted(df[conc_col].unique())):
            subset = df[df[conc_col] == conc]
            color = colors[idx % len(colors)]
            
            fig.add_trace(go.Scatter(
                x=subset['time_s'],
                y=subset['alpha'],
                mode='markers',
                name=f'{entity_name} {conc} {conc_unit}',
                marker=dict(size=6, color=color),
                showlegend=(idx == 0)
            ), row=1, col=1)
        
        # Plot model predictions
        line_styles = ['solid', 'dash', 'dot']
        for model_idx, result in enumerate(results):
            if result is None:
                continue
            
            for idx, conc in enumerate(sorted(df[conc_col].unique())):
                subset = df[df[conc_col] == conc]
                indices = subset.index
                
                fig.add_trace(go.Scatter(
                    x=subset['time_s'],
                    y=result.predictions[indices],
                    mode='lines',
                    name=result.name if idx == 0 else None,
                    line=dict(width=2, dash=line_styles[model_idx % 3]),
                    showlegend=(idx == 0),
                    legendgroup=result.name
                ), row=1, col=1)
        
        # Plot residuals
        for model_idx, result in enumerate(results):
            if result is None:
                continue
            
            fig.add_trace(go.Scatter(
                x=df['time_s'],
                y=result.residuals,
                mode='markers',
                name=f'{result.name} 잔차',
                marker=dict(size=4),
                showlegend=False
            ), row=2, col=1)
        
        fig.add_hline(y=0, line_dash="dash", line_color="gray", row=2, col=1)
        
        fig.update_xaxes(title_text=time_label, row=2, col=1)
        fig.update_yaxes(title_text="α", row=1, col=1)
        fig.update_yaxes(title_text="잔차", row=2, col=1)
        
        fig.update_layout(
            height=800,
            template='plotly_white',
            hovermode='x unified'
        )
        
        return fig
    
    @staticmethod
    def plot_initial_rates(df: pd.DataFrame, conc_unit: str = 'μg/mL', time_unit: str = 's'):
        """Plot initial rates v0 vs [E] to check linearity"""
        # Calculate initial slopes (0-2 time units)
        cutoff_time = 2 if time_unit == 's' else 0.5  # 2 seconds or 0.5 minutes
        initial_data = df[df['time_s'] <= cutoff_time].copy()
        
        rates = []
        concentrations = []
        
        conc_col = Visualizer._detect_conc_col(df)
        base = conc_col.split('_')[0].lower()
        entity_name = 'substrate' if base.startswith('pep') or base.startswith('sub') else 'enzyme'
        for conc in sorted(df[conc_col].unique()):
            subset = initial_data[initial_data[conc_col] == conc]
            if len(subset) >= 2:
                # Linear fit to get slope
                try:
                    # Check for valid data
                    if subset['alpha'].isnull().any() or np.isinf(subset['alpha']).any():
                        continue
                        
                    coeffs = np.polyfit(subset['time_s'], subset['alpha'], 1)
                    v0 = coeffs[0]  # slope = dα/dt
                    rates.append(v0)
                    concentrations.append(conc)
                except (np.linalg.LinAlgError, ValueError, TypeError):
                    # Skip this concentration if fit fails
                    continue
        
        fig = go.Figure()
        
        fig.add_trace(go.Scatter(
            x=concentrations,
            y=rates,
            mode='markers+lines',
            name='초기 속도',
            marker=dict(size=10)
        ))
        
        # Linear fit
        if len(concentrations) >= 2:
            coeffs = np.polyfit(concentrations, rates, 1)
            x_fit = np.linspace(min(concentrations), max(concentrations), 100)
            y_fit = coeffs[0] * x_fit + coeffs[1]
            
            fig.add_trace(go.Scatter(
                x=x_fit,
                y=y_fit,
                mode='lines',
                name=f'선형 피팅 (R² = {np.corrcoef(concentrations, rates)[0,1]**2:.4f})',
                line=dict(dash='dash')
            ))
        
        time_unit_label = "초" if time_unit == 's' else "분"
        time_unit_abbr = "s⁻¹" if time_unit == 's' else "min⁻¹"
        fig.update_layout(
            title=f'초기 속도 분석 (0-{cutoff_time}{time_unit_label})',
            xaxis_title=f'[{"기질" if entity_name=="substrate" else "효소"}] ({conc_unit})',
            yaxis_title=f'v₀ ({time_unit_abbr})',
            template='plotly_white'
        )
        
        return fig
    
    @staticmethod
    def plot_individual_model(df: pd.DataFrame, result: ModelResults, 
                             conc_unit: str = 'μg/mL', time_label: str = '시간 (초)',
                             model_color: str = '#FF6B6B'):
        """Plot a single model fit with experimental data"""
        fig = make_subplots(
            rows=2, cols=1,
            subplot_titles=(f'{result.name} 피팅 결과', '잔차'),
            vertical_spacing=0.15,
            row_heights=[0.7, 0.3]
        )
        
        colors = px.colors.qualitative.Set1
        
        # Plot experimental data
        conc_col = Visualizer._detect_conc_col(df)
        base = conc_col.split('_')[0].lower()
        entity_name = 'substrate' if base.startswith('pep') or base.startswith('sub') else 'enzyme'
        for idx, conc in enumerate(sorted(df[conc_col].unique())):
            subset = df[df[conc_col] == conc]
            color = colors[idx % len(colors)]
            
            fig.add_trace(go.Scatter(
                x=subset['time_s'],
                y=subset['alpha'],
                mode='markers',
                name=f'{entity_name} {conc} {conc_unit}',
                marker=dict(size=8, color=color),
                legendgroup=f'conc_{idx}'
            ), row=1, col=1)
            
            # Plot model fit for this concentration
            indices = subset.index
            fig.add_trace(go.Scatter(
                x=subset['time_s'],
                y=result.predictions[indices],
                mode='lines',
                name=f'{conc} {conc_unit} (피팅)',
                line=dict(width=3, color=color),
                showlegend=False,
                legendgroup=f'conc_{idx}'
            ), row=1, col=1)
            
            # Plot residuals for this concentration
            fig.add_trace(go.Scatter(
                x=subset['time_s'],
                y=result.residuals[indices],
                mode='markers',
                marker=dict(size=5, color=color),
                showlegend=False,
                legendgroup=f'conc_{idx}'
            ), row=2, col=1)
        
        fig.add_hline(y=0, line_dash="dash", line_color="gray", row=2, col=1)
        
        fig.update_xaxes(title_text=time_label, row=2, col=1)
        fig.update_yaxes(title_text="α", row=1, col=1)
        fig.update_yaxes(title_text="잔차", row=2, col=1)
        
        # Add statistics annotation
        stats_text = f"<b>피팅 통계</b><br>R² = {result.r_squared:.4f}<br>RMSE = {result.rmse:.4f}"
        fig.add_annotation(
            xref="paper", yref="paper",
            x=0.02, y=0.98,
            xanchor='left', yanchor='top',
            text=stats_text,
            showarrow=False,
            bgcolor="rgba(255,255,255,0.8)",
            bordercolor=model_color,
            borderwidth=2,
            font=dict(size=11),
            row=1, col=1
        )
        
        fig.update_layout(
            height=800,
            template='plotly_white',
            hovermode='x unified',
            showlegend=True
        )
        
        return fig
    
    @staticmethod
    def create_comparison_table(results: List[ModelResults]) -> pd.DataFrame:
        """Create comparison table for models"""
        data = []
        
        # Define relevant parameters for each model
        model_relevant_params = {
            "Model A: Substrate Depletion": ['kcat_KM'],
            "Model B: Enzyme Deactivation": ['kcat_KM', 'kd'],
            "Model C: Mass-Transfer Limitation": ['kcat_KM', 'km', 'Gamma0'],
            "Model D: Concentration-Dependent Fmax": ['kcat_KM', 'alpha_inf', 'k_access'],
            "Model E: Product Inhibition": ['kcat_KM', 'Ki_eff'],
            "Model F: Enzyme Surface Sequestration": ['kcat_KM', 'k_ads', 'K_ads']
        }
        
        # Master list of all parameters
        all_params = ['kcat_KM', 'kd', 'km', 'Gamma0', 'alpha_inf', 'k_access', 'Ki_eff', 'k_ads', 'K_ads']
        
        for result in results:
            if result is None:
                continue
            
            row = {
                '모델': result.name,
                'R²': f"{result.r_squared:.4f}",
                'RMSE': f"{result.rmse:.4f}",
                'AIC': f"{result.aic:.2f}",
                'BIC': f"{result.bic:.2f}",
                '파라미터 수': len(result.params)
            }
            
            # Get relevant params for this model
            # If model name matches partially (for robustness), find the key
            relevant = []
            for model_name, params in model_relevant_params.items():
                if result.name == model_name:
                    relevant = params
                    break
            
            if not relevant and result.name in model_relevant_params:
                 relevant = model_relevant_params[result.name]

            # Populate columns
            for param_key in all_params:
                if param_key in relevant:
                    # Relevant parameter
                    if param_key in result.params:
                        value = result.params[param_key]
                        # Format value
                        if param_key == 'kcat_KM':
                            row[param_key] = f"{value:.2e} M⁻¹s⁻¹"
                        elif param_key == 'kd':
                            row[param_key] = f"{value:.4f} s⁻¹"
                        elif param_key == 'km':
                            row[param_key] = f"{value:.2e} m/s"
                        elif param_key == 'Gamma0':
                            row[param_key] = f"{value:.2f} pmol/cm²"
                        elif param_key == 'alpha_inf':
                            row[param_key] = f"{value:.4f}"
                        elif param_key == 'k_access':
                            row[param_key] = f"{value:.2e} M⁻¹"
                        elif param_key == 'Ki_eff':
                            row[param_key] = f"{value:.4f}"
                        elif param_key == 'k_ads':
                            row[param_key] = f"{value:.4f} s⁻¹"
                        elif param_key == 'K_ads':
                            row[param_key] = f"{value:.2e} M⁻¹"
                        else:
                            row[param_key] = f"{value:.4e}"
                    else:
                        # Relevant but missing
                        row[param_key] = None
                else:
                    # Irrelevant parameter
                    row[param_key] = '-'
            
            data.append(row)
        
        return pd.DataFrame(data)

