# α(t) 플롯과 Alpha Statistics 테이블의 계산 방식 비교

이 문서는 Model Simulation 모드에서 **정규화 데이터: 절단 비율 α(t)** 플롯과 **Alpha Statistics by Concentration** 테이블이 각각 α를 어떻게 계산·집계하는지 비교합니다.

---

## 1. 공통 기반: 절단 비율 α(t)의 정의

두 출력 모두 **같은 기본 식**에서 출발합니다.

$$
\alpha(t) = \frac{F_t - F_0}{F_\infty - F_0}
$$

- **\(F_t\)**: 시간 \(t\)에서의 형광 강도
- **\(F_0\)**: 해당 농도에서의 초기 형광 (t = 0)
- **\(F_\infty\)**: 완전 절단 시 형광 (농도별 하나의 값, 또는 "Use shared F_∞" 시 모든 농도 공통)

계산 후 α는 다음 범위로 제한됩니다.

$$
\alpha(t) \in [0,\ 1]
$$

---

## 2. 플롯: 정규화 데이터 α(t) vs 시간

### 목적

각 **농도**별로 **시간에 따른 α(t)** 를 곡선으로 그려, 반응 진행을 시각화합니다.

### 계산 방식

- **x축**: 시간 \(t\) (초 단위, `time_s`)
- **y축**: 위 식으로 구한 **α(t)** (각 시간점마다 하나의 값)

즉, **한 농도·한 시간**에 대해 α를 한 번씩 계산합니다.

$$
\text{플롯의 한 점:}\quad \bigl(t,\ \alpha(t)\bigr) = \left(t,\ \frac{F_t - F_0}{F_\infty - F_0}\right)
$$

농도 \(c\)에 대해 시간점 \(t_1, t_2, \ldots, t_{N_c}\)마다 α를 계산해 이어서 곡선을 만듭니다. **시간에 대한 일대일 대응**이 핵심입니다.

---

## 3. 테이블: Alpha Statistics by Concentration

### 목적

각 **농도**별로, 그 농도에서 나온 **모든 α(t) 값**을 요약해 **평균(Alpha mean)** 과 **표준편차(Alpha std)** 를 표로 보여줍니다.

### 계산 방식

같은 농도 \(c\)에 속하는 모든 \((t_i, \alpha(t_i))\)를 모아, **통계량**만 계산합니다.

**Alpha mean (농도 \(c\)에 대한 α의 평균):**

$$
\overline{\alpha}_c = \frac{1}{N_c} \sum_{i=1}^{N_c} \alpha(t_i^{(c)})
$$

**Alpha std (농도 \(c\)에 대한 α의 표준편차):**

$$
\sigma_{\alpha,c} = \sqrt{\frac{1}{N_c} \sum_{i=1}^{N_c} \bigl( \alpha(t_i^{(c)}) - \overline{\alpha}_c \bigr)^2}
$$

- \(N_c\): 농도 \(c\)에서의 (시간점) 데이터 개수
- \(t_i^{(c)}\): 농도 \(c\)에서의 \(i\)번째 시간
- \(\alpha(t_i^{(c)})\): 1절의 식으로 구한 값

즉, **플롯에서 쓰는 α(t) 값들을 같은 농도 안에서 모아 평균·표준편차를 낸 것**이 테이블입니다.

### 3.1 평균이 “곡선 높이” 순서와 다르게 나올 수 있는 이유

Alpha mean은 **해당 농도에서의 모든 \((t_i, \alpha(t_i))\)에 대한 단순 산술평균**입니다. 따라서:

- 농도 A가 농도 B보다 **모든 시간에서** α가 더 크면 (곡선이 위에 있으면), 같은 시간 그리드에서는 **반드시** mean(A) > mean(B)입니다.
- 하지만 **농도마다 시간점 개수·범위가 다르면** 순서가 뒤바뀔 수 있습니다.
  - 예: 농도 1.25는 0~30분 전 구간에 점이 많이 있고, 농도 2.5는 5~30분 구간에만 점이 있다면, 1.25의 평균에는 “초반 저α” 점이 많이 포함되어 평균이 낮아질 수 있습니다.
  - 이 경우 **플롯**에서는 “같은 시간에서 1.25가 2.5보다 위”인데, **테이블**에서는 mean(2.5) > mean(1.25)처럼 보일 수 있습니다.

즉, **임의로 넣은 데이터가 아니라**, “전체 시간의 모든 데이터 점”으로 평균을 내는 방식 때문에, **시간점 분포가 농도마다 다르면** 곡선 순서와 평균 순서가 일치하지 않을 수 있습니다. 농도 간 비교가 필요하면 **동일한 시간 그리드**에서만 α를 모아 평균을 내는 방식이 더 적합할 수 있습니다.

---

## 4. 비교 요약

| 구분 | 플롯 (α(t) vs 시간) | 테이블 (Alpha Statistics) |
|------|----------------------|----------------------------|
| **입력** | 시간 \(t\), 형광 \(F_t\), \(F_0\), \(F_\infty\) | 플롯과 동일한 α(t) 값들 (같은 \(df\)) |
| **α 계산** | \(\alpha(t) = \dfrac{F_t - F_0}{F_\infty - F_0}\) (점별) | 동일 식으로 이미 계산된 \(\alpha(t_i)\) 사용 |
| **출력** | 농도별 **시계열** \((t,\, \alpha(t))\) | 농도별 **요약 통계** \(\overline{\alpha}_c\), \(\sigma_{\alpha,c}\) |
| **단위** | 한 점 = 한 시간에서의 α | 한 행 = 한 농도에서 전체 시간에 대한 평균·표준편차 |

**관계**: 테이블의 Alpha mean·std는 **플롯에 그려진 α(t) 곡선(같은 농도) 위의 점들**에 대해 평균과 표준편차를 낸 값입니다. 식은 하나(\(\alpha(t)\) 정의)이고, **플롯은 시계열로**, **테이블은 농도별 집계로** 사용하는 차이입니다.

### 4.1 플롯과 테이블이 사용하는 시간점·해상도

- **플롯**과 **테이블**은 **같은 DataFrame \(df\)** 를 사용합니다. 따라서 사용하는 시간점은 **완전히 동일**하고, resolution 차이는 없습니다. 테이블의 평균·표준편차는 “플롯에 그려진 그 점들”에 대한 것입니다.

### 4.2 농도별 시간점·해상도 (동일 시간 범위에서 resolution 차이)

- **농도 간** 시간점·해상도는 **데이터 소스**(Fitted Curves / 보간 곡선, 즉 \(df\_fitted\))가 어떻게 구성되어 있는지에 따라 결정됩니다.
- 코드에서는 `df_fitted`의 **전체** `Time_min` 값들을 모은 `unique_times`를 기준으로, 각 시간 \(t\)에 대해 “그 시간이 존재하는 (농도, RFU)” 조합만 행으로 넣어 \(df\_raw\) → 정규화 후 \(df\)를 만듭니다.  
  따라서 **소스에 어떤 (시간, 농도) 조합이 있느냐**에 따라 농도별로 갖는 시간점 집합이 달라질 수 있습니다.

| 데이터 소스 형태 | 농도별 시간점·해상도 |
|------------------|----------------------|
| **모든 농도가 같은 시간 그리드** (예: long format 보간 곡선에서 농도마다 동일한 `Time_min` 목록) | **동일**. 모든 농도가 같은 \(N\)개 시간점, 같은 시간 범위, 같은 resolution(예: \(\Delta t = 0.03\) min)을 가짐. |
| **농도마다 다른 시간 그리드** (예: 농도별로 측정 시점이 다르거나, 구간이 다른 경우) | **다를 수 있음**. 농도 \(c\)는 “그 농도에 대해 존재하는 \((t, F)\) 행”만 가지므로, 농도별로 \(N_c\), 시간 범위, 실질적 resolution이 달라질 수 있음. |

- 프로젝트에 포함된 **샘플** `data_interpolation_mode/results/MM_interpolated_curves.csv`에서는 모든 농도가 **동일한 1001개 시간점**(0~30 min, \(\Delta t \approx 0.03\) min)을 사용하므로, 농도 간 시간점·해상도 차이는 없습니다.

---

## 5. Raw data 루트 (둘이 사용하는 데이터 소스)

플롯과 테이블은 **같은 DataFrame \(df\)** 를 사용합니다. 그 \(df\)에 이르기까지의 raw data 경로는 아래와 같습니다.

### 5.1 데이터 소스 우선순위 (Model Simulation 진입 시)

**1단계: Fitted Curves 로드 → `df_fitted`**

| 우선순위 | 소스 | 설명 |
|----------|------|------|
| 0 | **Session (Data Load 결과)** | `st.session_state['interpolation_results']['interp_df']` — Data Load 모드에서 Run 후 메모리에 있는 보간 곡선 |
| 1 | **업로드 파일** | 사이드바에서 업로드한 CSV 또는 XLSX (시트: "Time–FLU Interpolated curves" 등) |
| 2 | **로컬 파일** | `Michaelis-Menten_calibration_results.xlsx`, `data_interpolation_mode/results/MM_interpolated_curves.csv` 등 |
| 3 | **내장 샘플** | 위 파일이 없을 때 사용하는 예시 데이터 |

`df_fitted` 컬럼 예: `Time_min`, 농도 컬럼(예: `Concentration [ug/mL]`), `RFU_Interpolated`(또는 `RFU_Calculated`).

**2단계: 표준 형식으로 변환 → `df_raw`**

- `df_fitted`를 **한 행 = (시간, 농도, 형광)** 형태로 변환.
- 출력 컬럼: `time_min`, `enzyme_ugml`, `FL_intensity`, `SD`(보간 데이터는 0).
- 코드: `app_ui/general_analysis_mode.py` 353–381라인 (루프로 행 생성 후 `pd.DataFrame(df_raw_converted)`).

**3단계: 단위 표준화 → `df_standardized`**

$$
\text{df\_raw}\ \xrightarrow{\ \texttt{UnitStandardizer.standardize}\ }\ \text{df\_standardized}
$$

- 시간 단위 통일, 농도/형광 컬럼 정리. (`app_ui/general_analysis_mode.py` 615–616라인)

**4단계: 정규화 및 α 계산 → `df_current` = \(df\)**

$$
\text{df\_standardized}\ \xrightarrow{\ \texttt{normalize\_temporary}\ }\ \xrightarrow{\ \texttt{divide\_regions}\ }\ \xrightarrow{\ \texttt{normalize\_final}\ }\ \text{df}
$$

- `normalize_final` 단계에서 **α 열 생성** (공식: \(\alpha = (F_t - F_0)/(F_\infty - F_0)\), 0–1 clip).
- 코드: `mode_general_analysis/analysis.py` 560–563라인; 호출은 `app_ui/general_analysis_mode.py` 627–635라인.

### 5.2 한 번에 보는 루트

```
[Raw source]
  → Session interp_df / 업로드 CSV·XLSX / 로컬 xlsx·csv / 샘플
       ↓
  df_fitted  (Time_min, 농도, RFU_Interpolated)
       ↓  변환 (353–381)
  df_raw  (time_min, enzyme_ugml, FL_intensity)
       ↓  UnitStandardizer.standardize (615–616)
  df_standardized
       ↓  normalize_temporary → divide_regions → normalize_final (627–633)
  df  (time_s, enzyme_ugml, alpha, F0, Fmax, ...)
       ↓
  ┌─────────────────────────────────────────────────────────┐
  │  플롯: x = df['time_s'],  y = df['alpha'] (농도별)     │
  │  테이블: 농도별 df['alpha'].mean(), df['alpha'].std()   │
  └─────────────────────────────────────────────────────────┘
```

즉, **플롯**과 **Alpha Statistics 테이블**은 같은 raw 루트에서 나온 **같은 \(df\)** 의 `time_s`·`alpha`·농도 컬럼만 각각 시계열 플롯용·집계용으로 사용합니다.
