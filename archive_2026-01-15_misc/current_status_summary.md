# 현재 진행상황 및 문제점 정리

## 작업 목표
신규 분리 미생물의 GSM 모델에서 FBA가 작동하지 않는 문제 해결
- 레퍼런스 모델 (Stenotrophomonas 폴더)에서는 FBA가 성공
- 신규 모델 (Stenotrophomonas-causal AI 폴더)에서는 FBA가 실패

---

## 지금까지 진행된 작업

### 1. 누락된 반응 분석 (초기 작업)
- 레퍼런스 모델과 신규 모델 비교
- 레퍼런스 모델에 있지만 신규 모델에 없는 반응: **114개** 발견
- 이전에 114개 반응을 추가했을 때: **성장률 1.410437로 성공**
- 최소 필요 반응: **99개** (9개 active + 90개 medium priority)

### 2. 핵심 반응 추가 (최근 작업)
- **4개 핵심 반응 추가**:
  1. **ACS_ADP**: Acetyl-CoA synthetase (ADP-forming)
     - 반응식: `ac_c + atp_c + coa_c <=> accoa_c + adp_c + pi_c`
     - Acetate를 Acetyl-CoA로 전환 (ATP → ADP, 더 효율적)
  2. **SUCDi**: Succinate dehydrogenase (irreversible)
     - 반응식: `q8_c + succ_c --> fum_c + q8h2_c`
     - TCA cycle에서 Succinate를 Fumarate로 전환 (비가역)
  3. **PEPCK_ATP**: Phosphoenolpyruvate carboxykinase (ATP)
     - 반응식: `atp_c + oaa_c --> adp_c + co2_c + pep_c`
     - Gluconeogenesis 경로 (OAA → PEP)
  4. **ACtexi**: Acetate transport (exchange)
     - 반응식: `ac_e <=> ac_c`
     - Acetate를 세포 외부에서 내부로 운반

### 3. Exchange bounds 수정
- 레퍼런스 모델과 동일하게 설정
- **14개 Exchange bounds 수정**:
  - EX_ac_e: [-19, -19] → [-1000, 1000]
  - EX_o2_e: [-100, 1000] → [-1000, 1000]
  - EX_hco3_e: 추가됨 [-1000, 1000]
  - 기타 Exchange들도 레퍼런스와 동일하게 수정

### 4. ATPM bounds 수정
- 레퍼런스 모델: [0.0, 1000.0]
- 신규 모델 (수정 전): [10.0, 1000.0]
- **수정 후**: [0.0, 1000.0] (레퍼런스와 동일)

### 5. 모델 파일 저장
- **BaseModel_with_ACtexi.xml**: 모든 수정사항 반영
  - 4개 반응 추가 완료
  - Exchange bounds 수정 완료
  - ATPM bounds 수정 완료

---

## 현재 상황

### FBA 실행 결과 (BaseModel_with_ACtexi.xml)

#### ATPM=0일 때
- **성장률**: 0.000000
- **상태**: optimal
- **주요 반응 플럭스**: 모두 0
  - CS: 0.000000
  - ACS_ADP: 0.000000
  - SUCDi: 0.000000
  - ICL: 0.000000
  - MALS: 0.000000

#### ATPM>=5일 때
- **성장률**: 0.000000 (여전히 성장하지 않음)
- **상태**: optimal
- **경로는 작동함**:
  - ATPM=5: CS=0.909091, SUCDi=0.909091
  - ATPM=10: CS=1.818182, SUCDi=1.818182
  - ATPM=20: CS=3.636364, SUCDi=3.636364
- **ACS_ADP**: 여전히 0 (작동하지 않음)

---

## 발견된 문제점

### 1. ATPM=0일 때 경로가 작동하지 않음
**원인 분석**:
- ATPM 반응식: `atp_c + h2o_c --> adp_c + h_c + pi_c`
- ATPM은 **ADP를 생성**함
- ATPM=0이면 ADP가 생성되지 않음
- ATP 생성 경로 (ATPS4rpp, PYK 등)는 **ADP를 필요로 함**
- **결과**: ADP가 없어서 ATP 생성 경로가 작동하지 않음 → 전체 경로 막힘

**참고사항**:
- 이전에 성장했을 때 (성장률 1.410437)도 ATPM 설정 확인 필요
- 레퍼런스 FBA 결과 파일 (`fba_flux_gradient_acid.csv`)에서는 ATPM=0일 때도 경로가 작동함 (ACS_ADP=0.943448)

### 2. ATPM>=5일 때 경로는 작동하지만 성장하지 않음
**현상**:
- CS, SUCDi 등 일부 경로는 작동함
- 하지만 **ACS_ADP는 여전히 작동하지 않음**
- 성장률은 0

**가능한 원인**:
- ACS_ADP 반응이 추가되었지만 플럭스가 0
- 다른 제약 조건이나 경로 문제
- CoA 생성 문제 (이전에 확인했던 문제)

### 3. 이전 성공 사례와의 차이
**이전에 성장했을 때**:
- 114개 반응 추가 시: 성장률 1.410437
- 미디어 설정: setup_acetate_medium_growing (EX_ac_e=-10 등)
- ATPM 설정: 기본값 (확인 필요)

**현재 상황**:
- 4개 핵심 반응만 추가
- Exchange bounds는 레퍼런스와 동일하게 설정
- ATPM=0에서 작동하지 않음

**차이점**:
- 이전에는 **114개 반응**을 모두 추가했음
- 현재는 **4개 핵심 반응**만 추가
- 나머지 110개 반응이 필요할 수 있음

---

## 핵심 문제 요약

### 문제 1: ATPM=0에서 경로 작동 불가
- **원인**: ATPM이 ADP를 생성하는 반응이므로, ATPM=0이면 ADP 부족
- **현상**: ATP 생성 경로가 작동하지 않아 전체 경로 막힘
- **참고**: 레퍼런스 FBA 결과에서는 ATPM=0일 때도 작동함 (다른 설정일 수 있음)

### 문제 2: ATPM>=5일 때 경로는 작동하지만 ACS_ADP가 작동하지 않음
- **현상**: CS, SUCDi는 작동하지만 ACS_ADP는 0
- **가능한 원인**: 
  - CoA 생성 문제 (이전에 확인)
  - 다른 경로 문제
  - 추가 반응 필요

### 문제 3: 성장률이 0
- **현상**: 경로가 작동해도 성장률이 0
- **가능한 원인**:
  - 미디어 설정 문제
  - 추가 반응 필요 (나머지 110개 반응)
  - Biomass 반응 문제

---

## 다음 단계 제안

### 옵션 1: 이전 성공 설정 재현
- 114개 반응을 모두 추가
- 이전에 성공했을 때의 미디어 설정 사용
- ATPM 설정 확인

### 옵션 2: ACS_ADP 작동 문제 해결
- ACS_ADP가 작동하지 않는 원인 확인
- CoA 생성 경로 확인
- 필요한 전구체 확인

### 옵션 3: 레퍼런스 FBA 결과 재현
- 레퍼런스 FBA 결과 파일의 정확한 설정 확인
- 그 설정을 재현

---

## 현재 모델 파일
- **BaseModel_with_ACtexi.xml**: 
  - 4개 반응 추가 완료 (ACS_ADP, SUCDi, PEPCK_ATP, ACtexi)
  - Exchange bounds 수정 완료
  - ATPM bounds 수정 완료 ([0.0, 1000.0])

---

## FBA 관점에서의 핵심 해석 (최신 분석)

### 현재 관찰된 패턴의 의미

**핵심 해석**:
- "현재 제약조건에서 biomass 반응은 애초에 flux가 날 수 없다(= 최대 성장률이 0)"
- ATPM=0에서 all-zero flux는 "ADP가 없어서"라기보다, **현재 조건에서 biomass가 infeasible(최대 성장률=0)**이기 때문에 FBA가 **trivial 해(0 flux)**를 반환하는 전형적인 상황
- ATPM을 걸면 유지비를 충족해야 하니 catabolism/호흡은 도는데, biomass 경로가 막혀 성장은 여전히 0

**정상적인 ATPM 거동**:
- ATPM=0에서 성장률이 가장 높고
- ATPM을 올릴수록 성장률이 단조 감소
- 어떤 ATPM 이상에서 성장 0(maintenance-only)로 전환
- 그 구간에서 TCA/호흡 flux와 CO₂가 증가하는 패턴

### 원인 (가장 흔한 4가지)

1. **Biomass 전구체/보조인자 생산 경로가 빠져 있음** (가장 흔함)
   - 114개 누락 반응 발견
   - 과거에 114개 반응 추가 시 성장률 1.41로 성공
   - BCAA, Proline, NAD 계열 경로가 1순위 후보

2. **인실리코 배지가 실험 배지와 다름** (yeast extract/비타민)
   - 실험 배지에 yeast extract가 있다면 외부 공급 필요

3. **Adenylate/인산 연결(ATP–ADP–AMP)이 끊김**
   - ADK1 같은 연결 반응 누락 가능성

4. **수송/compartment/방향성 문제**
   - 필수 무기물 수송 제약

### 해결을 위한 진단 순서

**Step 4-1**: 배지 조건을 강제로 고정해서 성장 가능성 판정 ✅ (완료)
- 결과: 성장률 0 (biomass가 막혀 있음)

**Step 4-2**: Biomass가 요구하는 성분 중 뭐가 안 만들어지는지 찾기 (다음 단계)

**Step 4-3**: 막힌 전구체 주변을 우선 복구

**Step 4-4**: ATPM-ADP 원인 확인 ✅ (완료)
- 결과: "ADP Deadlock" 가설은 성립하지 않음
- ATPM=0일 때도 성장 가능 (성장률: 3.29)
- ADP는 다양한 경로에서 충분히 생성됨 (265개 이상의 반응)
- 실제 원인: Biosynthetic Gaps (BCAA, NAD/NADP, 이온 수송 경로 누락)

**Step 4-5**: 조건부 gap-filling 수행

### 해결 루트

**루트 A**: 빠르게 성장 복구(실무형)
- 과거 114개 반응 세트 재적용 → 성장 가능한 baseline 복원 → 반응 삭제로 최소화

**루트 B**: 원인까지 깔끔하게 규명(연구형)
- Biomass 성분 demand 테스트 → 막힌 전구체 pinpoint → targeted gapfill

상세 내용은 `fba_interpretation_and_solution.md` 파일 참조
