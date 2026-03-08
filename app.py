#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# .\.venv\Scripts\python.exe -m streamlit run app.py .\.venv\Scripts\python.exe -m streamlit run app.py
# author: hyacinth1126
"""
Hydrogel FRET Advanced Kinetic Analysis - Streamlit Application
"""
import sys
import traceback
from datetime import datetime

# Cloud 로그용: stderr에 출력 후 flush (Community Cloud 로그에만 보임)
def _cloud_log(msg: str) -> None:
    try:
        ts = datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")
        sys.stderr.write(f"[{ts}] [app] {msg}\n")
        sys.stderr.flush()
    except Exception:
        pass

def _cloud_log_exception() -> None:
    if sys.exc_info()[0] is not None:
        sys.stderr.write(traceback.format_exc())
        sys.stderr.flush()

_cloud_log("app.py: loading (imports starting)")

try:
    import streamlit as st
    _cloud_log("app.py: streamlit imported")
    from app_ui.footer import render_footer
    _cloud_log("app.py: footer imported")
except Exception:
    _cloud_log("app.py: import failed")
    _cloud_log_exception()
    raise

# 무거운 모듈은 각 모드 선택 시 로드 (Cloud spawn error 방지 + 실패 시 에러 메시지 확인 가능)
def _configure_plotting():
    _cloud_log("_configure_plotting(): importing matplotlib, seaborn")
    import matplotlib.pyplot as plt
    import seaborn as sns
    plt.rcParams['font.family'] = 'DejaVu Sans'
    sns.set_style("whitegrid")
    _cloud_log("_configure_plotting(): done")


def main():
    """Main Streamlit app"""
    _cloud_log("main(): entered")
    st.set_page_config(
        page_title="Hydrogel FRET Advanced Analysis",
        page_icon="🔬",
        layout="wide"
    )
    _cloud_log("main(): set_page_config done")
    st.title("🔬  FRET Protease Simulation")
    st.markdown("---")
    
    # Mode selection
    analysis_mode = st.sidebar.radio(
        "Select Analysis Mode",
        ["Data Load Mode", "Model Simulation Mode"],
        help="Data Load Mode: Upload CSV file or extract data from image | Model Simulation Mode: Standard FRET analysis"
    )
    _cloud_log(f"main(): mode selected = {analysis_mode!r}")
    # Always render footer at bottom
    render_footer()
    
    # Data Load Mode
    if analysis_mode == "Data Load Mode":
        try:
            _cloud_log("main(): configuring plotting")
            _configure_plotting()
            _cloud_log("main(): importing data_load_mode")
            from app_ui.data_load_mode import data_load_mode
            _cloud_log("main(): calling data_load_mode(st)")
            data_load_mode(st)
            _cloud_log("main(): data_load_mode(st) returned")
        except Exception as e:
            _cloud_log("main(): Data Load Mode failed")
            _cloud_log_exception()
            st.error("Data Load Mode 로드 중 오류가 발생했습니다.")
            st.code(str(e), language="text")
            st.exception(e)
        return
    
    # Model Simulation Mode
    try:
        _cloud_log("main(): configuring plotting (Model Simulation)")
        _configure_plotting()
        _cloud_log("main(): importing general_analysis_mode")
        from app_ui.general_analysis_mode import general_analysis_mode
        _cloud_log("main(): calling general_analysis_mode(st)")
        general_analysis_mode(st)
        _cloud_log("main(): general_analysis_mode(st) returned")
    except Exception as e:
        _cloud_log("main(): Model Simulation Mode failed")
        _cloud_log_exception()
        st.error("An error occurred while loading Model Simulation Mode.")
        st.code(str(e), language="text")
        st.exception(e)


if __name__ == "__main__":
    _cloud_log("__main__: calling main()")
    main()
    _cloud_log("__main__: main() returned")
