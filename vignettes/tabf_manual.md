# `tabf` 패키지 사용자 가이드북 (AMH 프레임워크 전용)

`tabf`는 AMH(AI-driven epidemiological Analysis and Modeling) 프레임워크에 최적화된 **역학 연구용 기초 통계 및 표 생성 R 패키지**입니다. 복잡한 코딩 없이 논문에 들어갈 `Table 1` (기초통계표)과 다변량 로지스틱 회귀분석 결과표(OR, 95% CI)를 즉시 추출하고, Forest Plot 시각화까지 도와줍니다.

---

## 1. 패키지 설치 및 로드표준 세팅

`tabf` 환경을 구성하려면 아래의 코드를 사용하세요. 이미 설치되어 있으면 로드만 수행하고, 설치되어 있지 않으면 GitHub에서 자동으로 다운로드 받습니다.

```r
# 필수 패키지 설치 및 로드
if(!require("devtools")) install.packages("devtools")
library(devtools)

# tabf (v2) 설치 및 로드
if(!require("tabf")) install_github("jinha-dspub/tabf")
library(tabf)

library(dplyr)
library(ggplot2)
```

---

## 2. 데이터 탐색 및 전처리 (Phase 2)

### 2.1 결측치 파악 (`missingf`)
분석 전 데이터의 변수별 결측치 갯수와 비율을 쉽게 파악합니다.

```r
# 사용법
missing_summary <- missingf(my_data)
print(missing_summary)

# 결과 예시
# Variables   Total N   Missing N  Missing %
# wtime_con1     5000          50       1.0%
# age            5000          25       0.5%
```

### 2.2 집단 구분 없는 전체 특성 요약 (Univariate Table 1)
종속 변수(`stratas`)를 넣지 않고, **연구 대상자 전체의 기초 특성(단일 집단)**을 한눈에 파악할 때 유용합니다. (Total N, %, 평균±SD 제공 / p-value 제외)

```r
tab_univ <- tabf(
  dat1 = my_data, 
  stratas = NULL,                     # 단일 집단 모드
  catVars = c("gender", "education"), # 범주형 변수
  conVars = c("age", "bmi")           # 연속형 변수
)

# HTML 표로 출력 시
library(htmlTable)
tab_univ %>% htmlTable()
```

---

## 3. 기초 통계표 생성 (Table 1) (Phase 3)

논문 맨 첫 페이지에 들어가는 **특성별 집단 비교 표(Table 1)**를 생성합니다.

### 3.1 집단 간 특성 비교 (`tabf`)
가장 널리 쓰이는 형태이며, **열(Column) 기준의 백분율(%)**을 계산합니다.
- `stratas` 변수에 들어가는 집단의 갯수를 자동으로 인식하여 연속형 변수에 대해 **2개 그룹이면 T-test, 3개 이상이면 ANOVA**를 자체적으로 판단하여 p-value를 도출합니다.

```r
tab1 <- tabf(
  dat1 = my_data,
  stratas = "disease_status",         # 그룹을 나눌 기준 (ex: 0=정상, 1=질병)
  catVars = c("gender", "smoking"),   # 카이제곱 검정 수행
  conVars = c("age", "bmi")           # T-test 또는 ANOVA 수행
)

tab1 %>% htmlTable(caption = "Table 1. Baseline Characteristics")
```

### 3.2 행(Row) 기준 백분율 계산 (`tabf2`)
`tabf()`와 인자는 동일하지만, 백분율을 층화 변수별이 아닌 **변수 내부 수준별(Row)로 계산**합니다. 유병률/발생률 표기 시 유용합니다.

```r
tab2_row <- tabf2(dat1 = my_data, stratas = "disease_status", catVars = c("gender"))
```

---

## 4. 로지스틱 회귀분석 요약표 (Phase 4)

각종 변수를 보정한(Adjusted) 여러 모델을 묶어서 바로 논문에 복붙 가능한 오즈비(OR) 테이블로 만들어줍니다.

### 4.1 통합 OR 테이블 구성 (`oddsTabf`)
Model 1(Crude), Model 2(부분 보정), Model 3(완전 보정)의 회귀분석 결과를 묶어서 보여줍니다.

```r
# 모델 생성
mod1 <- glm(disease_status ~ exposure, data = my_data, family = binomial())
mod2 <- glm(disease_status ~ exposure + age + gender, data = my_data, family = binomial())
mod3 <- glm(disease_status ~ exposure + age + gender + bmi, data = my_data, family = binomial())

# 통합 테이블 추출
or_table <- oddsTabf(mod1, mod2, mod3)
or_table %>% htmlTable(caption = "Table 2. Odds Ratios for Disease Status")

# 출력 예시
# Variables  Values  Crude            Age&Gender      Fully Adjusted
# exposure        1  1.00 (ref)       1.00 (ref)      1.00 (ref)
# exposure        2  1.52 (1.1-2.0)   1.48 (1.0-1.9)  1.32 (0.9-1.8)
```

---

## 5. 출판용 시각화 (Phase 5)

### 5.1 Forest Plot 생성 (`forestf`)
로지스틱 회귀 모델을 통째로 넣어 시각적 인사이트(숲 그림)를 바로 뽑아냅니다. p-value < 0.05에 따라 색상이 자동으로 강조됩니다.

```r
# 모델을 집어넣어 바로 플롯 저장
plt <- forestf(mod3, title = "Adjusted Odds Ratios for Disease Status")

# 이미지 파일로 저장 (ggplot2 객체)
ggsave("forest_plot.png", plot = plt, width = 7, height = 5)
```

---

## 💡 AI 챗봇(Agent)을 위한 프롬프트 가이드
분석을 수행할 때 AI에게 다음과 같이 지시하면 최고의 결과를 얻을 수 있습니다:
> *"Phase 2: `tabf` 패키지의 `missingf()`를 사용해 결측치 테이블을 보여줘. 그리고 `stratas=NULL` 조건으로 `tabf()`를 돌려서 전체 대상자의 특성표를 만들어봐."*
> *"Phase 3: `tabf()` 함수를 써서 질병 유무(disease_status)에 따른 Table 1을 HTML 테이블로 출력해줘. 연속형은 age, 범주형은 gender야."*
> *"Phase 4 및 5: `oddsTabf(mod1, mod2, mod3)`을 실행해서 OR 요약표를 만들어줘. 그리고 완전 보정 모델인 `mod3`를 `forestf()`에 넣어서 Forest Plot으로 시각화해줘."*
