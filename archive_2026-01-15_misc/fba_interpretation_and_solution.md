# FBA 관점에서의 ATPM=0 문제 해석 및 해결 방안

## 핵심 해석: ATPM=0에서 성장이 0인 현상

### 현재 관찰된 패턴

**ATPM=0일 때**:
- 성장률: 0
- CS, SUCDi 등 주요 반응 flux: 모두 0
- Solver status: optimal (해는 존재함)

**ATPM≥5일 때**:
- 성장률: 여전히 0
- CS, SUCDi 같은 TCA/호흡 관련 flux: ATPM에 비례하여 증가
  - ATPM=5: CS=0.909091, SUCDi=0.909091
  - ATPM=10: CS=1.818182, SUCDi=1.818182

### 핵심 해석

**"현재 제약조건에서 biomass 반응은 애초에 flux가 날 수 없다(= 최대 성장률이 0)"**

1. **ATPM=0일 때**:
   - Biomass=0을 만족하는 해가 무수히 많음
   - Solver가 가장 "쉬운" 해(대개 all-zero flux)를 반환
   - → "전부 0"이 나옴
   - **이건 ADP가 '없어서'라기보다, 성장(또는 ATP 요구)이 전혀 없으니 굳이 돌릴 이유가 없어서 0이 되는 것**

2. **ATPM을 양(+)으로 걸 때**:
   - "ATPM flux를 만족시키려면 ATP를 만들어서 계속 태워야 함"이라는 추가 요구조건이 생김
   - Biomass는 여전히 0이지만, ATP를 만들기 위해 acetate 산화 + 호흡(TCA/ETC) 쪽 flux가 강제로 돌아감
   - → CS, SUCDi 등이 ATPM에 비례해서 커지는 현상

**✅ 결론**: 
- "ATPM이 있어야 돌아간다"는 건 보통
  - (A) 성장(바이오매스) 자체가 막혀서 최적 성장=0인 상태에서
  - (B) ATPM이 '유지비'를 강제하니 catabolism만 돌아가는 상태로 이해하는 게 가장 자연스럽습니다.

---

## 왜 biomass가 막혔나? (가장 흔한 원인 4가지)

### 원인 1: Biomass 전구체/보조인자 생산 경로가 빠져 있음 (가장 흔함)

**현황**:
- Reference 대비 누락 반응: 114개 발견
- 과거에 114개 반응을 추가했을 때: 성장률 1.41로 "성공"
- 첨부 초안에서도 acetate minimal growth를 위해 특정 경로(BCAA, proline, NAD/니코티네이트, 수송반응 등)를 gap-filling로 추가했다고 언급

**결론**: "새 분리주 초안 모델이 최소배지에서 성장 못 하는 건 흔한 일"입니다.

### 원인 2: 인실리코 배지가 실험 배지와 다름 (특히 yeast extract/비타민)

**현황**:
- 실험 배지에 yeast extract가 들어갔다면, 실제로는 외부에서 비타민/아미노산/조효소 전구체를 얻고 있을 수 있음
- 모델은 그것을 공급하지 않으면 auxotrophy처럼 성장 불가

**결론**: Reference strain 모델은 내부 합성 경로가 갖춰져서 성장하지만, 새 분리주 모델은 그 경로를 빠뜨려 성장 못 하는 경우가 잦습니다.

### 원인 3: Adenylate/인산 연결(ATP–ADP–AMP)이 끊겨 있음

**가능성**:
- ADP를 생성하는 반응(ATP→ADP)이 사실상 없거나
- Acetate 활성화가 AMP-forming인데 adenylate kinase(ADK1: AMP+ATP ↔ 2ADP) 같은 연결이 누락
- → ATP 생성(ATPS, PYK 등)이 막히며 성장도 막힐 수 있음

**결론**: ATPM을 억지로 걸면(그리고 호흡으로 ATP를 만들 수 있으면) ATPM ↔ ATPS 형태의 유지비 루프가 돌아가면서 "돌아가는 것처럼" 보일 수 있습니다.

### 원인 4: 수송/compartment/방향성 문제

**가능성**:
- Acetate는 들어오는데 NH4, Pi, SO4, 금속이온, H2O 등 필수 무기물이 닫혀 있으면 biomass는 못 만듦
- Periplasm/cytosol 구획 모델에서 proton/oxygen 관련 수송 제약이 꼬이면, respiration은 되는데 특정 전구체 합성이 막히는 식의 이상이 생김

---

## 해결을 위한 "진단 순서"

### Step 4-1: 배지 조건을 "강제로 고정"해서 성장 가능성부터 판정

**작업**:
1. Acetate uptake를 고정
   - 예: `EX_ac_e.lower_bound = EX_ac_e.upper_bound = -19` (또는 실험값)
2. O₂ uptake도 제한/고정
   - 예: `EX_o2_e.lower_bound = -100` (초안과 유사한 설정)
3. 나머지 필수 무기물(NH4, Pi, SO4, Mg, K, Na 등)은 충분히 열어두기

**판정 기준**:
- 이렇게 "탄소/전자수용체가 확실히 들어오게" 만든 상태에서
- `maximize biomass` 했는데도 μ=0이면, **진짜로 biomass가 막혀 있는 것**
- ("optional uptake라서 안 먹은" 문제가 아닙니다)

### Step 4-2: "Biomass가 요구하는 성분 중 뭐가 안 만들어지는지" 찾기

**가장 빠른 방법**:
- Biomass 반응이 소비(–)하는 각 metabolite m_i에 대해 DM_mi (demand)를 만들고
- `max DM_mi`가 0이면 그 전구체/보조인자 생산 경로가 끊긴 것

**결과**:
- 이 작업을 하면 "114개 중 뭘 넣어야 하나?"가 바로 핵심 후보 몇 개로 줄어듭니다.

### Step 4-3: 막힌 전구체가 나오면, 그 주변을 우선 복구

**보통 acetate minimal growth에서 자주 막히는 모듈**:
1. Glyoxylate shunt (ICL, MALS)
2. Gluconeogenesis/PPP (PEP↔OAA, FBPase 등)
3. BCAA (Val/Leu/Ile), Proline
4. NAD/니코티네이트, CoA(판토텐산), Quinone(유비퀴논/메나퀴논)
5. 지질/세포벽 전구체 (ACP, UDP-당 등)

**참고**: 첨부 초안에서도 실제로 acetate growth를 위해 수송 + BCAA/Proline + NAD 계열 경로를 gap-filling로 보완했다고 서술되어 있어, 이 모델에서도 이쪽이 1순위 후보입니다.

### Step 4-4: "ATPM 때문에 ADP가 생긴다"가 진짜 원인인지 확인하는 체크

**1분 안에 확인 가능**:
1. 모델에서 ADK1(adenylate kinase), NDPK, PPase 같은 연결 반응이 있는지 확인
2. Acetate 활성화가 AMP-forming(Acs) 위주면 ADK1이 없을 때 adenylate pool이 끊기는 경우가 많음
3. Biomass 반응 자체가 ATP를 소비하면서 ADP를 생성하지 않는 비정상 stoichiometry인지 확인

### Step 4-5: 마지막으로 gap-filling을 "조건부(acetate minimal + biomass)"로 수행

**전략**:
- 목표: acetate minimal 조건에서 biomass > 0
- 최소 추가 반응 수: 가능한 작게
- 조건부 gapfill을 돌려 "정말 필요한 반응만" 뽑기

**팁**:
- 한 번에 114개를 다 넣기보다
- Gapfill이 제안한 최소 집합 → 성장 확인
- 그 집합에서 reaction essentiality(단일 반응 knockout)로 "진짜 핵심"만 남기기
- → 모델 해석력이 확 좋아집니다.

---

## 지금 현상을 어떻게 바로 해결할 수 있나? (권장 루트 2개)

### 루트 A: 빠르게 성장 복구(실무형)

**전략**:
1. 과거에 성장시켰던 99개/114개 반응 세트를 다시 적용해서
2. "성장 가능한 baseline"을 먼저 복원
3. 그 상태에서 **반응 삭제(리덕션)**를 하며 "성장에 꼭 필요한 반응"만 남기기

**장점**: 
- 이미 "9개 active + 90개 medium"까지 좁힌 히스토리가 있으니, 이게 가장 빠름

### 루트 B: 원인까지 깔끔하게 규명(연구형)

**전략**:
1. Step 4-2의 "biomass 성분 demand 테스트"로 막힌 전구체를 pinpoint
2. 그 전구체 주변 경로만 targeted gapfill/수기 큐레이션
3. 실험 배지에 yeast extract/비타민이 있으면 모델에도 해당 uptake를 열어
4. "실험 조건 그대로" 먼저 맞춘 뒤, 이후 완전 최소배지로 단계적으로 줄이기

---

## 최종 검증: 정상화되면 ATPM 스윕이 이렇게 나와야 합니다

**모델이 제대로 연결되면**:
- ATPM=0에서 성장률이 가장 높고
- ATPM을 올릴수록 성장률이 단조 감소
- 어떤 ATPM 이상에서 성장 0(maintenance-only)로 전환
- 그 구간에서 TCA/호흡 flux와 CO₂가 증가하는 패턴

**이게 첨부 초안에서 기술한 "energetic stress에 따른 decoupling" 그림과 일치하는 정상적인 거동입니다.**

---

## (짧게) 지금 로그에 기반한 "가장 가능성 높은 결론"

1. **ATPM=0에서 all-zero flux는 "ADP가 없어서"라기보다**
   - **현재 조건에서 biomass가 infeasible(최대 성장률=0)**이기 때문에 FBA가 **trivial 해(0 flux)**를 반환하는 전형적인 상황일 가능성이 큼

2. **ATPM을 걸면 유지비를 충족해야 하니 catabolism/호흡은 도는데, biomass 경로가 막혀 성장은 여전히 0**

3. **해결책**:
   - Biomass 전구체 생산 경로 복구 (114개 반응 중 필요한 것만 추가)
   - 특히 BCAA, Proline, NAD 계열 경로가 1순위 후보
