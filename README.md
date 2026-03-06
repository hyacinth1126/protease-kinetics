# Protease Simulator 3

## 빠른 실행 (venv)

이 프로젝트는 가상환경(venv) 사용을 권장합니다. 다른 프로젝트와 패키지 충돌을 피할 수 있습니다.

**1. venv 활성화 (PowerShell)**

```powershell
cd c:\hyacinth1126\6-protease-simulator-3
.\venv\Scripts\Activate.ps1
```

활성화되면 프롬프트 앞에 `(venv)`가 붙습니다.

**2. 앱 실행 (Streamlit)**

```powershell
streamlit run app.py
```

브라우저에서 **http://localhost:8501** 로 접속합니다.

**3. 종료**

터미널에서 `deactivate` 입력 시 venv를 빠져나옵니다.

---

## 다른 PC에서 셋업 (clone / pull 후 venv 만들기)

`venv` 폴더는 `.gitignore`에 있어서 저장소에 포함되지 않습니다. 다른 컴퓨터에서 `git clone` 또는 `git pull` 한 뒤에는 **가상환경을 새로 만들고** 의존성을 설치해야 합니다.

**1. 저장소 받기**

```bash
# 처음 받을 때
git clone <저장소-URL>
cd 6-protease-simulator-3

# 이미 클론한 뒤 최신만 받을 때
git pull
```

**2. 가상환경 생성**

프로젝트 루트에서 실행합니다.

```bash
python -m venv venv
```

**3. 가상환경 활성화**

- **Windows (PowerShell)**  
  ```powershell
  .\venv\Scripts\Activate.ps1
  ```
- **Windows (CMD)**  
  ```cmd
  venv\Scripts\activate.bat
  ```
- **macOS / Linux**  
  ```bash
  source venv/bin/activate
  ```

프롬프트 앞에 `(venv)` 가 붙으면 활성화된 것입니다.

**4. 의존성 설치**

```bash
pip install -r requirements.txt
```

**5. 앱 실행**

```bash
streamlit run app.py
```

이후에는 이 PC에서도 **빠른 실행 (venv)** 절의 1–3단계(활성화 → 실행 → 종료)만 반복하면 됩니다.

---

## Streamlit Community Cloud 배포 시 참고

- **venv는 로컬 전용입니다.** Cloud는 저장소에 올라간 코드만 사용하며, `.gitignore`에 있는 `venv/`·`.venv/`는 배포되지 않습니다. Cloud 측에서 자체 가상환경에 `requirements.txt`로 의존성을 설치합니다.
- **"Oh no. Error running app" / `spawn error`** 가 나고 로그에 에러가 안 보일 때:
  1. **진입 파일을 최소 앱으로 바꿔 보기**  
     Cloud 앱 설정에서 **Main file**을 `app_cloud_debug.py`로 변경한 뒤 재배포하세요.  
     - **"Cloud OK" 화면이 뜨면** → 문제는 `app.py` 또는 의존성 로드 쪽입니다. Main file을 다시 `app.py`로 되돌리고, **Advanced settings**에서 **Python 버전**을 **3.12**로 선택한 뒤 앱을 삭제하고 같은 설정으로 재배포해 보세요.  
     - **최소 앱도 안 뜨면** → Streamlit/Cloud 환경 문제 가능성이 큽니다. Python 3.12로 재배포하거나 [Community Cloud 상태](https://www.streamlitstatus.com/)를 확인하세요.
  2. **서버 설정**  
     `.streamlit/config.toml`에 `server.address = "0.0.0.0"` 이 있어야 health check가 통과합니다. (이미 설정됨)
  3. **의존성**  
     `requirements.txt`에 `protobuf>=3.20,<6` 제한이 있어야 합니다. (이미 추가됨)
  4. **Export plots 탭에서 PNG 다운로드**  
     Kaleido 1.x는 Chrome/Chromium이 필요합니다. Cloud에서는 저장소 루트의 `packages.txt`에 `chromium` 한 줄이 있으면 자동 설치됩니다. (이미 추가됨)
- 로컬에서는 `streamlit run app.py` 만 실행하면 됩니다.

---

단백질 분해 효소(Protease) 반응 시뮬레이션 및 데이터 분석 도구입니다.

---

## 프로젝트 구조 및 역할

| 폴더 / 파일 | 역할 |
|-------------|------|
| **`app.py`** | Streamlit 앱 진입점. 모드 선택 후 해당 UI로 분기합니다. |
| **`app_cloud_debug.py`** | Cloud 디버깅용 최소 앱. spawn error 시 Main file을 이걸로 바꿔 환경 정상 여부를 확인할 수 있습니다. |
| **`packages.txt`** | Cloud 배포 시 apt로 설치할 패키지 목록. `chromium`으로 Export plots PNG 내보내기가 동작합니다. |
| **`app_ui/`** | 웹 UI 코드. Data Load / Model Simulation 모드 화면과 로직을 담당합니다. |
| **`mode_prep_raw_data/`** | 원본(raw) 데이터 읽기, 시간–곡선 피팅, MM/선형 보정, 초기 속도(v₀) 계산 등 **Data Load 모드의 핵심 연산**을 수행합니다. |
| **`data_interpolation_mode/`** | Prism 스타일 보간(Exponential Association 등). 보간 곡선 생성 및 결과 저장을 담당합니다. |
| **`mode_general_analysis/`** | Model Simulation 모드의 **핵심 로직**: 단위 표준화, 정규화, 구간 분할, 6가지 반응 모델 피팅 및 비교 분석입니다. |
| **`flowchart/`** | 전체 워크플로우·플로우차트 문서 및 관련 스크립트입니다. |
| **`convert_to_csv/`** | `raw/` 폴더의 XLSX를 prep_raw 형식 CSV로 일괄 변환하는 스크립트가 들어 있습니다. |
| **`raw/`** | 실험 원본 데이터(XLSX/CSV). 앱의 기본 샘플 및 변환 스크립트의 입력으로 사용됩니다. |
| **`csv/`** | `convert_to_csv` 스크립트가 생성한 CSV 결과가 저장되는 폴더입니다. |
| **`prep_raw_data_mode/results/`** | Data Load 모드에서 생성되는 상세 MM 결과(예: `MM_results_detailed.csv`)가 저장됩니다. |
| **`data_interpolation_mode/results/`** | Data Load 모드에서 생성되는 보간 곡선 CSV(예: `MM_interpolated_curves.csv`)가 저장됩니다. |

---

## 주요 코드 설명

### 1. 앱 실행 및 모드 분기
- **`app.py`**  
  - `python -m streamlit run app.py` 로 실행합니다.  
  - **Data Load Mode**: `app_ui.data_load_mode.data_load_mode(st)` 호출 → 실험 데이터 업로드/로드, 피팅, 보간, MM 보정, 정규화, Excel 저장까지 수행.  
  - **Model Simulation Mode**: `app_ui.general_analysis_mode.general_analysis_mode(st)` 호출 → Data Load 결과 또는 업로드 파일로 모델 시뮬레이션·비교 분석.

### 2. Data Load 모드 (데이터 로드 → MM 보정 → 보간 → 정규화)
- **`app_ui/data_load_mode.py`**  
  - 실험 타입 선택(Substrate/Enzyme 농도 변화), 파일 업로드 또는 `raw/` 기본 샘플 사용.  
  - `mode_prep_raw_data.prep`의 `read_raw_data`, `fit_time_course`, `fit_calibration_curve`, `calculate_initial_velocity` 등을 사용해 v₀ 계산 및 MM/선형 피팅.  
  - `data_interpolation_mode.interpolate_prism`으로 보간 범위·곡선 생성.  
  - 정규화 후 `Michaelis-Menten_calibration_results.xlsx` 생성 및 세션에 결과 저장.

- **`mode_prep_raw_data/prep.py`**  
  - **어떤 것**: raw CSV/XLSX 읽기, 농도별 시간–곡선 피팅, 초기 속도(v₀) 계산, Michaelis–Menten/선형 보정 곡선 피팅.  
  - **어떻게**: `read_raw_data()`로 원본 로드 → `fit_time_course()`로 곡선 피팅 → `calculate_initial_velocity()`로 v₀ → `fit_calibration_curve()`로 Vmax, Km, kcat 등 계산.

- **`data_interpolation_mode/interpolate_prism.py`**  
  - **어떤 것**: Prism 스타일 보간(Exponential Association 등), 보간 범위 계산, 곡선 생성.  
  - **어떻게**: `exponential_association()`, `create_prism_interpolation_range()` 등으로 보간 포인트 생성 후 CSV/세션에 저장.

### 3. Model Simulation 모드 (모델 피팅 및 비교)
- **`app_ui/general_analysis_mode.py`**  
  - Data Load 결과(세션 또는 업로드 CSV/XLSX)를 불러와 단위 표준화·정규화·구간 분할 후, 6가지 모델로 피팅하고 결과를 시각화·다운로드.

- **`mode_general_analysis/analysis.py`**  
  - **어떤 것**: 단위 표준화(UnitStandardizer), 정규화(DataNormalizer), 구간 분할(RegionDivider), 6가지 반응 모델(Substrate Depletion, Enzyme Deactivation, Mass-Transfer 등) 피팅 및 AIC/BIC/R² 비교.  
  - **어떻게**: 입력 데이터에 표준화·정규화 적용 → 구간 나눔 → 각 모델로 글로벌 피팅 → 최적 모델 선택 및 파라미터 출력.

- **`mode_general_analysis/plot.py`**  
  - 피팅 결과 시각화(Visualizer).

### 4. XLSX → CSV 변환 (배치)
- **`convert_to_csv/convert_xlsx_to_csv.py`**  
  - **어떤 것**: `raw/` 폴더 안의 모든 XLSX를 prep_raw 형식의 CSV로 변환.  
  - **어떻게**: 인자 없이 실행하면 `raw/*.xlsx`를 읽어 각각 `csv/<동일파일명>.csv`로 저장. 단일 파일 변환 시에는 `python convert_to_csv/convert_xlsx_to_csv.py raw/파일명.xlsx` 형태로 실행 가능.  
  - **실행 예**:  
    ```bash
    python convert_to_csv/convert_xlsx_to_csv.py
    ```

---

## 주요 기능 요약

*   **Data Load 모드**: 실험 데이터(CSV/XLSX) 로드 및 전처리  
    *   Substrate 농도 변화 / Enzyme 농도 변화 실험 지원  
    *   Michaelis–Menten 피팅 및 초기 속도(v₀) 계산  
    *   보간(Interpolation) 및 정규화(Normalization)  
    *   결과 → `Michaelis-Menten_calibration_results.xlsx` 및 results 폴더에 CSV 저장  
*   **Model Simulation 모드**: 분석된 데이터로 6가지 반응 모델 피팅·비교 및 시각화  
*   **convert_to_csv**: `raw/` 내 XLSX → `csv/` 폴더에 prep_raw 형식 CSV 일괄 생성  

---

## 설치 방법

1.  Python 3.8 이상 설치  
2.  의존성 설치:  
    ```bash
    pip install -r requirements.txt
    ```

---

## 실행 방법

```bash
python -m streamlit run app.py
```

---

## 데이터 형식 (prep_raw.csv / prep_raw.xlsx)

*   **1행**: 농도 값 (각 농도마다 mean, SD, N 3열 반복)  
*   **2행**: 컬럼 헤더 (time_min, mean, SD, N, …)  
*   **3행 이후**: 시간(min), 농도별 mean, SD, N  
*   앱에서 **Download Sample**로 `raw_substrate.xlsx` / `raw_enzyme.xlsx` 형식 샘플을 받을 수 있음  
