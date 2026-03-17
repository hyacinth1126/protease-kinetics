# LOD / LOQ 계산 방식 정리

본 문서는 protease-simulator-3에서 **검출 한계(LOD)** 와 **정량 한계(LOQ)** 를 어떻게 계산하는지 정리한 것입니다.  
**Enzyme 농도 변화 실험**(기질 고정, v₀ vs [E] 선형)에서만 농도 단위 LOD/LOQ가 계산됩니다.

---

## 1. 정의

| 용어 | 의미 |
|------|------|
| **LOD** (Limit of Detection) | blank와 통계적으로 구별 가능한 최소 신호/농도. 일반적으로 blank mean + 3σ 사용. |
| **LOQ** (Limit of Quantification) | 정량적으로 신뢰할 수 있는 최소 신호/농도. 일반적으로 blank mean + 10σ 사용. |

본 코드에서는 **3σ (LOD)** 와 **10σ (LOQ)** 계수를 기본으로 사용합니다.

---

## 2. Blank (PBS) 데이터 입력

LOD/LOQ 계산을 위해 **선택적으로** blank(PBS) 데이터를 raw 데이터에 넣을 수 있습니다.

### 2.1 새 형식 CSV

- **PBS** 또는 **Blank** 컬럼: blank 신호 값(평균 또는 복수 측정값).
- **PBS_SD** 또는 **Blank_SD** 컬럼 (선택): 한 지점에서 준 **표준편차**.
  - 한 시점(예: 0, 1, 5, 10, 15, 20, 25, 30 min 중 하나)에서 **mean**과 **SD**를 함께 주는 경우, 해당 (mean, SD) 쌍을 그대로 사용합니다.
  - SD 컬럼이 없으면: PBS 컬럼의 모든 비-NaN 값으로 **평균**과 **표준편차**(표본 SD, ddof=1)를 계산합니다.

### 2.2 구 형식 (탭 구분)

- 마지막에 **1개 컬럼**이 더 있으면 그 값을 blank로 사용합니다 (여러 행이 있으면 값들의 평균·표준편차 계산).

### 2.3 내부 저장

- Blank가 있으면 `raw_data['_blank']` 에 다음이 저장됩니다.
  - `mean`: blank 평균
  - `sd`: blank 표준편차 (한 지점 mean/SD를 준 경우 해당 SD, 아니면 값들로부터 계산)
  - `n`: 사용된 값 개수(또는 한 지점이면 1)
  - `values`: 원본 값 배열
- `_blank` 는 농도 조건으로 취급되지 않으며, 피팅/보정 곡선 계산 시 제외됩니다.

---

## 3. 신호 공간 LOD/LOQ (Signal)

Blank mean과 SD가 있으면 **신호 공간**에서 LOD·LOQ를 다음 식으로 계산합니다.

$$
\mathrm{LOD}_{\mathrm{signal}} = \mu_{\mathrm{blank}} + k_{\mathrm{LOD}} \cdot \sigma_{\mathrm{blank}}
$$

$$
\mathrm{LOQ}_{\mathrm{signal}} = \mu_{\mathrm{blank}} + k_{\mathrm{LOQ}} \cdot \sigma_{\mathrm{blank}}
$$

- **기본값**: \(k_{\mathrm{LOD}} = 3\), \(k_{\mathrm{LOQ}} = 10\).
- **Blank SD = 0** (한 점만 있거나 복제 1개)인 경우:  
  \(\mathrm{LOD}_{\mathrm{signal}} = \mathrm{LOQ}_{\mathrm{signal}} = \mu_{\mathrm{blank}}\) 로 둡니다.

**코드**: `mode_prep_raw_data.prep.compute_lod_loq_signal(blank_mean, blank_sd, k_lod=3, k_loq=10)`

**단위**: LOD(signal), LOQ(signal)은 **신호와 같은 단위**입니다. 본 앱에서 v₀(초기 속도)가 RFU/min이면 **RFU/min**입니다.

---

## 4. 농도 공간 LOD/LOQ (Concentration)

Enzyme 농도 변화 실험에서는 v₀ vs [E] 가 선형으로 피팅됩니다:

$$
v_0 = \mathrm{slope} \times [E] + \mathrm{intercept}
$$

이 관계를 이용해 **신호 LOD/LOQ**를 **농도 LOD/LOQ**로 변환합니다.

### 4.1 Blank 기반 (선형 보정 곡선 사용)

신호 공간 LOD/LOQ를 위 선형 식의 **역변환**으로 농도로 바꿉니다.

$$
[E]_{\mathrm{LOD}} = \frac{\mathrm{LOD}_{\mathrm{signal}} - \mathrm{intercept}}{\mathrm{slope}}
,\qquad
[E]_{\mathrm{LOQ}} = \frac{\mathrm{LOQ}_{\mathrm{signal}} - \mathrm{intercept}}{\mathrm{slope}}
$$

- **slope ≤ 0** 이면 농도 LOD/LOQ는 계산하지 않습니다.
- 음수 농도는 0으로 클리핑합니다.

**단위**: LOD([E]), LOQ([E])는 **농도 단위**입니다. 효소 농도가 μg/mL이면 **μg/mL**입니다.

**LOD(signal)과 LOD([E]) 숫자가 다른 이유**: 같은 “검출 한계”를 **신호 단위**(예: 2637.70 RFU/min)로 쓰면 LOD(signal), **농도 단위**(예: 0.5843 μg/mL)로 쓰면 LOD([E])입니다. 보정 곡선으로 서로 변환한 값이라 단위가 다르고, 따라서 **숫자도 당연히 다릅니다.**

**코드**: `compute_lod_loq_concentration_from_linear(lod_signal, loq_signal, slope, intercept)`

- Blank(PBS) 컬럼이 있고, 위에서 계산한 \(\mathrm{LOD}_{\mathrm{signal}}\), \(\mathrm{LOQ}_{\mathrm{signal}}\) 이 있을 때만 이 방식으로 **농도 LOD/LOQ**가 채워집니다.

### 4.2 보정 곡선 잔차 기반 (Residual)

Blank의 SD를 쓰지 않고, **v₀ vs [E] 선형 회귀의 잔차 표준편차**로 농도 LOD/LOQ를 추정할 수 있습니다.  
Blank가 한 점만 있거나 SD가 없을 때 유용합니다.

**왜 Blank mean/SD를 쓰지 않나?**  
두 방식은 **정의가 다릅니다.**  
- **Blank 기반 (4.1)**: “blank 신호의 변동(σ_blank)”을 기준으로, blank와 구별 가능한 최소 신호 → 농도로 환산.  
- **잔차 기반 (4.2)**: “보정 곡선 주변의 측정 산포(잔차의 σ)”를 기준으로, **농도 축에서의 불확실성**만으로 LOD/LOQ를 추정.  

즉, 잔차 방식은 “각 농도에서 v₀가 피팅 직선에서 얼마나 흩어져 있는지”만 사용하므로, **Blank mean/SD는 이 공식에 들어가지 않습니다.**  
(Blank 데이터가 있어도 “From blank (PBS)” 결과에만 쓰이고, “From calibration residual” 계산에는 사용하지 않습니다.)

- 잔차 제곱합: \(\mathrm{SS}_{\mathrm{res}} = \sum_i (v_{0,i} - \hat{v}_{0,i})^2\)
- 자유도: \(n - 2\) (선형 회귀)
- 잔차 표준편차: \(s_y = \sqrt{\mathrm{SS}_{\mathrm{res}} / (n-2)}\)

$$
[E]_{\mathrm{LOD}} = k_{\mathrm{LOD}} \cdot \frac{s_y}{\mathrm{slope}}
,\qquad
[E]_{\mathrm{LOQ}} = k_{\mathrm{LOQ}} \cdot \frac{s_y}{\mathrm{slope}}
$$

- **기본값**: \(k_{\mathrm{LOD}} = 3\), \(k_{\mathrm{LOQ}} = 10\).
- slope > 0 이고 \(s_y\) 가 계산 가능할 때만 사용합니다.

**코드**: `compute_lod_loq_from_residual(residual_std, slope, k_lod=3, k_loq=10)`

- **Data Load** 실행 시 v₀ vs [E] 선형 피팅이 성공하면, Blank 유무와 관계없이 **농도 LOD/LOQ (from residual)** 이 계산되면 결과에 포함됩니다.

---

## 5. 적용 조건 요약

| 항목 | 조건 |
|------|------|
| 실험 유형 | **Enzyme 농도 변화** (기질 고정, v₀ vs [E] 선형) |
| 신호 LOD/LOQ | Raw 데이터에 **Blank(PBS)** 컬럼이 있을 때 |
| 농도 LOD/LOQ (blank 기반) | Blank 컬럼 + **v₀ vs [E] 선형 피팅** (slope > 0) |
| 농도 LOD/LOQ (residual 기반) | v₀ vs [E] 선형 피팅 성공 시 (Blank 없어도 계산) |

Substrate 농도 변화(표준 Michaelis–Menten) 실험에서는 현재 **농도 LOD/LOQ는 계산하지 않습니다**.

---

## 6. 결과 표시 및 내보내기

- **Data Load 모드**  
  - Data Preview: Blank 컬럼이 있으면 blank n, mean, SD를 요약 표시.  
  - MM 결과 아래: **LOD/LOQ** 섹션에서  
    - **From blank (PBS)**: Blank mean/SD, LOD(signal), LOQ(signal), LOD([E]), LOQ([E]) (μg/mL)  
    - **From calibration residual**: LOD([E]), LOQ([E]) (μg/mL)
- **엑셀 내보내기**  
  - **Fit results** 시트에 Blank mean, Blank SD, Blank n, LOD(signal), LOQ(signal), LOD([E]), LOQ([E]), residual 기반 LOD/LOQ가 포함됩니다.

---

## 7. 참고 문헌·관례

- IUPAC 등에서 널리 쓰는 **blank + 3σ (LOD), 10σ (LOQ)** 관례를 따릅니다.
- 보정 곡선 잔차를 이용한 LOD/LOQ 추정은 선형 보정에서 흔히 사용되는 방식입니다.
