#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# .\.venv\Scripts\python.exe -m streamlit run app.py .\.venv\Scripts\python.exe -m streamlit run app.py
# author: hyacinth1126
"""
Hydrogel FRET Advanced Kinetic Analysis - Streamlit Application
"""

import streamlit as st
from app_ui.footer import render_footer

# 무거운 모듈은 각 모드 선택 시 로드 (Cloud spawn error 방지 + 실패 시 에러 메시지 확인 가능)
def _configure_plotting():
    import matplotlib.pyplot as plt
    import seaborn as sns
    plt.rcParams['font.family'] = 'DejaVu Sans'
    sns.set_style("whitegrid")


def main():
    """Main Streamlit app"""
    st.set_page_config(
        page_title="Hydrogel FRET Advanced Analysis",
        page_icon="🔬",
        layout="wide"
    )
    
    st.title("🔬  FRET Protease Simulation")
    st.markdown("---")
    
    # Mode selection
    analysis_mode = st.sidebar.radio(
        "Select Analysis Mode",
        ["Data Load Mode", "Model Simulation Mode"],
        help="Data Load Mode: Upload CSV file or extract data from image | Model Simulation Mode: Standard FRET analysis"
    )
    # Always render footer at bottom
    render_footer()
    
    # Data Load Mode
    if analysis_mode == "Data Load Mode":
        try:
            _configure_plotting()
            from app_ui.data_load_mode import data_load_mode
            data_load_mode(st)
        except Exception as e:
            st.error("Data Load Mode 로드 중 오류가 발생했습니다.")
            st.code(str(e), language="text")
            st.exception(e)
        return
    
    # Model Simulation Mode
    try:
        _configure_plotting()
        from app_ui.general_analysis_mode import general_analysis_mode
        general_analysis_mode(st)
    except Exception as e:
        st.error("Model Simulation Mode 로드 중 오류가 발생했습니다.")
        st.code(str(e), language="text")
        st.exception(e)


if __name__ == "__main__":
    main()
