# 플롯별 사용 식 정리 (Enzyme 농도 변화 기준)

저장되는 각 플롯에서 쓰이는 수식을 정리한 문서입니다.  
**본 문서는 효소 농도 변화(Enzyme Concentration Variation, 기질 고정) 실험만을 기준으로 합니다.** 기질 농도 변화(Substrate 농도 변화, 표준 MM) 실험에 대한 식은 포함하지 않습니다.

---

## 파라미터 정의 (논문용)

본 문서의 식에 등장하는 기호는 아래와 같이 정의하며, Methods 또는 기호 설명란에 그대로 인용할 수 있다.

| 기호 | 의미 (논문용 설명) |
|------|---------------------|
| $F(t)$ | **시간 $t$에서의 형광 강도** (RFU). 반응 진행에 따라 변하는 형광 신호. |
| $F_0$ | **초기 형광** (baseline fluorescence). $t=0$에서의 형광 강도. |
| $F_{\max}$ | **포화 형광** (plateau / maximum fluorescence). 반응이 완료에 가깝게 진행된 후 도달하는 최대 형광 강도. |
| $F_\infty$ | **완전 절단 시 형광** (fluorescence at full cleavage). 농도별 하나의 값이거나, 공통 $F_\infty$를 쓸 수 있다. |
| $F_t$ | **시간 $t$에서의 측정 형광**. $\alpha(t)$ 계산 시 해당 시점의 관측값. |
| $F_{\mathrm{norm}}$ | **정규화 형광**. $(F - F_0)/(F_{\max} - F_0)$로, 0~1 범위의 비율(진행도). |
| $k_{\mathrm{obs}}$ | **관측 속도 상수** (observed rate constant, min$^{-1}$). progress curve 피팅에서 얻은 지수 상수. |
| $\tau$ | **시간 상수**. $\tau = 1/k_{\mathrm{obs}}$. |
| $t$ | **시간** (min 또는 s). 반응 시작($t=0$)으로부터의 경과 시간. |
| $v_0$ | **초기 속도** (initial velocity). 반응 초기 구간에서의 형광 증가율 (RFU/min 등). |
| $[E]$ | **효소 농도** (enzyme concentration, μg/mL 등). |
| $\alpha(t)$ | **절단 비율** (fraction cleaved). $\alpha = (F_t - F_0)/(F_\infty - F_0)$, 0~1. |
| $\Gamma$, $\Gamma_0$, $\Gamma_t$ | **기질(또는 잔여 기질) 양**. Model A 전역 피팅에서 사용; $\alpha = 1 - \Gamma_t/\Gamma_0$. |
| $k_{\mathrm{cat}}/K_M$ | **특이성 상수**. 단순 1차 근사에서 $k_{\mathrm{obs}} \propto (k_{\mathrm{cat}}/K_M)[E]$. |

**Enzyme 농도 변화 실험**  
기질 농도는 고정하고 효소 농도 $[E]$만 변할 때, progress curve는 지수 포화 또는 초기 구간에서 선형 $F(t) = F_0 + v_0 t$로 근사할 수 있으며, $v_0$는 $[E]$에 비례한다.

---

## 기질 농도 고정 시 v₀ vs [E] 선형 관계 유도 (Michaelis–Menten)

효소 농도 변화 실험에서는 **기질 농도 [S]를 고정**하고 효소 농도 [E]만 바꾼다. 이때 초기 속도 $v_0$가 [E]에 비례하는 것은 Michaelis–Menten 식에서 다음과 같이 유도된다.

**1. Michaelis–Menten 식 (기질 농도에 따른 초기 속도)**

$$
v_0 = \frac{V_{\max} \cdot [S]}{K_M + [S]}
$$

- $[S]$: 기질 농도  
- $V_{\max}$: 최대 속도 (기질이 포화될 때의 속도)  
- $K_M$: Michaelis 상수  

**2. $V_{\max}$와 효소 농도 [E]의 관계**

총 효소 농도를 $[E]_T$라 하면, $V_{\max} = k_{\mathrm{cat}} \cdot [E]_T$ 이다. (단위에 맞게 [E]_T를 M 등으로 쓰면 $k_{\mathrm{cat}}$은 turnover number.)

**3. [S]를 상수로 두면**

[S]를 고정한 실험에서는 $[S]$가 변하지 않으므로

$$
\frac{V_{\max} \cdot [S]}{K_M + [S]} = \frac{k_{\mathrm{cat}} \cdot [S]}{K_M + [S]} \cdot [E]_T
$$

에서 $\dfrac{k_{\mathrm{cat}} \cdot [S]}{K_M + [S]}$ 는 [E]에 무관한 **상수**이다. 이를 기울기로 두면

$$
v_0 = \underbrace{\frac{k_{\mathrm{cat}} \cdot [S]}{K_M + [S]}}_{\text{slope}} \cdot [E]_T
$$

**4. 정리**

기질 농도 [S]가 고정일 때 **$v_0$는 [E]에 비례**한다. 따라서 v₀ vs [E] 플롯은 **직선**이 되며, 기울기 = $k_{\mathrm{cat}}[S]/(K_M + [S])$ 이다.  
(실제 피팅에서는 절편을 넣은 $v_0 = \mathrm{slope} \cdot [E] + \mathrm{intercept}$ 로 쓸 수 있다.)

---

## 1. Data Load 모드 — Export 플롯 (Enzyme 농도 변화)

### 1.1 Experimental_Results

- **설명**: 원시 데이터만 표시 (Time vs RFU). 곡선 식 없음.
- **축**: x = Time (min), y = RFU

---

### 1.2 Time_Fluorescence_Interpolated_Curves

- **설명**: 농도별(효소 농도 [E]별) Time–RFU 보간 곡선.
- **곡선 식 (Enzyme 농도 변화)**

  Data Load 모드에서 곡선을 만드는 방식이 두 가지 있다.

  **1) 선형 구간 사용 시**  
  반응 **초기**에 형광이 거의 직선으로 올라가는 구간만 쓰는 경우. 이 구간의 기울기를 초기 속도 $v_0$로 두고, 그 직선을 곡선 전체의 대표로 쓸 때 다음 식을 쓴다.
$$
F(t) = F_0 + v_0 \, t
$$

  **2) 지수 포화 사용 시**  
  시간에 따른 형광 증가를 **지수 포화 곡선**으로 피팅하는 경우. 전체 시간 구간에 대해 다음 식을 쓴다.
$$
F(t) = F_0 + (F_{\max} - F_0)(1 - e^{-k t})
$$

- **기호**: $F_0$ = t=0 형광, $F_{\max}$ = 포화 형광, $k$ = 속도 상수, $v_0$ = 초기 속도, $t$ = 시간 (min).

---

### 1.3 Normalized_Time_Fluorescence_exponential_curves_3tau_max_range

- **설명**: 정규화된 형광 곡선 (0~1), 지수 구간 강조, 최대 3τ 범위.
- **정규화**

$$
F_{\mathrm{norm}} = \frac{F - F_0}{F_{\max} - F_0}
$$

- **피팅 곡선 (정규화 공간, $F_{\max}=1$)**
$$
F_{\mathrm{norm}}(t) = 1 - e^{-k_{\mathrm{obs}} t}
$$

- **초기 선형 구간 (대시선)**
$$
y = k_{\mathrm{obs}} \, t
$$

- **시간 상수**: $\tau = 1 / k_{\mathrm{obs}}$.

---

### 1.4 Linear_fit (v₀ vs [E])

- **설명**: 기질 고정 시 v₀ vs [E] 선형 피팅.
- **곡선 식**

$$
v_0 = \mathrm{slope} \cdot [E] + \mathrm{intercept}
$$

- **기호**: $[E]$ = 효소 농도 (μg/mL 등), slope·intercept = 피팅 파라미터.

---

### 1.5 Supplementary_low3_const_fit

- **설명**: 낮은 [E] 3점만 사용한 보조 선형 피팅 (상수항 포함).
- **곡선 식**: $v_0 = a \cdot [E] + b$.

---

### 1.6 Normalization_{농도} (Enzyme-quenched peptide fluorescence kinetics)

- **설명**: [E]별 정규화 데이터 + 지수 피팅 + 초기 선형 구간.
- **정규화 데이터**: $F_{\mathrm{norm}} = (F - F_0)/(F_{\max} - F_0)$.
- **지수 곡선 (주황)**

$$
F_{\mathrm{norm}}(t) = F_{\max}(1 - e^{-k_{\mathrm{obs}} t}), \quad F_{\max}=1
$$

- **초기 선형 (점선)**
$$
y = k_{\mathrm{obs}} \, t
$$

- **v₀ (원시 단위)**: $v_0 = k_{\mathrm{obs}}(F_{\max} - F_0)$.

---

### 1.7 Normalization_to_plateau_{농도}

- **설명**: 위와 동일 식, y≈1(plateau) 도달 구간까지만 표시.

---

## 2. Data Load 모드 — 정규화 피팅 (내부, Enzyme 농도 변화)

- **정규화**: $F_{\mathrm{norm}} = (F - F_0)/(F_{\max} - F_0)$, $F_0$ = 첫 시간점 값, $F_{\max}$ = 최대값.
- **정규화 공간 지수 피팅**

$$
F_{\mathrm{norm}}(t) = F_{\max}(1 - e^{-k_{\mathrm{obs}} t})
$$

($F_{\max}$ 피팅, 보통 ≈1.)
- **역정규화**: $F(t) = F_{\mathrm{norm}}(t)(F_{\max} - F_0) + F_0$.

---

## 3. Model Simulation (General Analysis) — α 및 모델 (Enzyme 농도 변화)

### 3.1 절단 비율 α(t)

$$
\alpha(t) = \frac{F_t - F_0}{F_\infty - F_0}
$$

- $F_t$: 시간 $t$에서 형광, $F_0$: 초기 형광, $F_\infty$: 완전 절단 시 형광. $\alpha \in [0, 1]$로 클리핑.
- **Shared F_∞** 옵션: 모든 [E]에 하나의 $F_\infty$ 사용 (논문에서 흔한 방식).

---

### 3.2 Progress curve 피팅 (형광 F vs t)

- **3-파라미터**

$$
F(t) = F_0 + (F_\infty - F_0)(1 - e^{-k_{\mathrm{obs}} t})
$$

- **2-파라미터** ($F_0=0$ 가정)

$$
F(t) = F_\infty(1 - e^{-k_{\mathrm{obs}} t})
$$

---

### 3.3 Model A: Substrate Depletion (단순 1차, [E] 변수)

$$
\alpha(t) = 1 - e^{-k_{\mathrm{obs}} t}, \qquad k_{\mathrm{obs}} = (k_{\mathrm{cat}}/K_M) \cdot [E]
$$

- $[E]$: 효소 농도, $k_{\mathrm{cat}}/K_M$ 피팅.

---

### 3.4 Model A: Full (Michaelis–Menten, 수치 적분)

- $d\Gamma/dt = -v$, $v = k_{\mathrm{cat}} [E] \Gamma / (K_M + \Gamma)$.
- $\alpha = 1 - \Gamma_t / \Gamma_0$, $\alpha \in [0, 1]$ 클리핑.

---

### 3.5 [E] vs α 플롯

- **x축**: $[E]$ (효소 농도, μg/mL 등)
- **y축**: $\alpha$ (절단 비율, 평균 또는 시점값)
- **피팅**: Exponential·Hyperbolic(MM형) 등으로 R²·AIC 비교. (Enzyme 농도 변화 실험에서 α mean vs [E] 곡선.)

---

## 요약 표 (Enzyme 농도 변화만)

| 플롯/출력 | 사용 식 |
|-----------|---------|
| Time-Fluorescence Interpolated | $F(t) = F_0 + v_0 t$ 또는 $F_0 + (F_{\max}-F_0)(1 - e^{-kt})$ |
| Normalized exponential curves | $F_{\mathrm{norm}} = 1 - e^{-k_{\mathrm{obs}} t}$, $y = k_{\mathrm{obs}} t$ |
| v₀ vs [E] | $v_0 = \mathrm{slope}\cdot [E] + \mathrm{intercept}$ |
| Normalization 탭 | $F_{\mathrm{norm}} = (F-F_0)/(F_{\max}-F_0)$, $1 - e^{-k_{\mathrm{obs}} t}$ |
| α(t) (Model Simulation) | $\alpha = (F_t - F_0)/(F_\infty - F_0)$ |
| Model A 단순 | $\alpha = 1 - e^{-k_{\mathrm{obs}} t}$, $k_{\mathrm{obs}} \propto [E]$ |
