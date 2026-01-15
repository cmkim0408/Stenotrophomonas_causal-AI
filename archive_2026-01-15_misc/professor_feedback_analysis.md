# 교수님 피드백 분석 및 해결 방안

## 문제 원인 분석

### 1. ADP Deadlock (ADP 고갈 문제)

**원인**:
- ATPM 반응식: `ATP + H₂O → ADP + Pi + H⁺`
- ATPM은 단순히 에너지를 소모하는 것이 아니라, **ADP를 재생산(Recycling)**하는 핵심 역할
- ATPM=0으로 강제하면, 시스템 내에서 ADP가 재생산되지 않음
- 산화적 인산화(ATP Synthase)가 작동하려면 기질로 ADP가 필요한데, ADP 공급원이 차단되어 전체 에너지 대사가 멈춤

**Reference 모델이 ATPM=0에서도 작동한 이유**:
- Adenylate Kinase (ADK1) (`ATP + AMP ↔ 2ADP`) 같은 우회 경로가 살아있었을 가능성
- 현재 신규 모델에서는 이 경로가 막혀있거나(AMP 생성 경로 부재), 초기 조건에서 ADP가 부족한 상태

---

### 2. Biosynthetic Gaps (바이오매스 합성 경로의 결실)

**현재 상황**:
- 114개 누락 반응 중 '핵심 4개(수송, TCA 등)'만 추가됨
- "연료(Acetate)를 태우는 엔진"은 복구했지만, **"세포를 만드는 조립 라인"은 여전히 끊겨 있음**

**문제점**:
- 바이오매스 식(Biomass Objective Function)은 모든 아미노산을 필요로 함
- 에너지가 아무리 많아도, 단 하나의 필수 아미노산(예: Valine)이라도 합성되지 않으면 생장률은 수학적으로 0

**논문의 언급**:
> "초기 시뮬레이션 결과, 발린(Valine), 류신(Leucine), 프롤린(Proline) 및 필수 조효소(NAD/NADP)의 생합성 경로가 불완전하여 바이오매스 합성을 막고 있음이 확인되었다."

---

## 해결 방안

논문의 [Materials and Methods 2.3] 섹션에 언급된 **19개 필수 반응**을 모두 추가해야 합니다.

### 1. 분지쇄 아미노산(BCAA) 합성 (Valine/Leucine 등)

#### ALS (Acetolactate synthase)
- 반응식: `2 pyruvate → 2-acetolactate + CO₂`
- 유전자: ilvB, ilvI

#### IPMS (Isopropylmalate synthase)
- 반응식: `2-oxoisovalerate + acetyl-CoA → 3-isopropylmalate + CoA`
- 유전자: leuA

#### IPMD (3-Isopropylmalate dehydrogenase)
- 반응식: `3-isopropylmalate + NAD⁺ → 2-oxoisocaproate + CO₂ + NADH`
- 유전자: leuB
- 참고: comprehensive_missing_reactions.csv에는 IPMDH로 표기됨 (`3ippm_c + nad_c --> 4mop_c + co2_c + nadh_c`)

#### KARI (Ketol-acid reductoisomerase) / 2OXOVISO (2-Oxoisovalerate synthase)
- 반응식: `2-acetolactate + NADPH + H⁺ → 2,3-dihydroxyisovalerate + NADP⁺`
- 유전자: ilvC
- 참고: comprehensive_missing_reactions.csv에는 KARI로 표기됨 (`alac__S_c + h_c + nadph_c --> dhiv_c + nadp_c`)

#### DHAD (Dihydroxyacid dehydratase)
- 반응식: `2,3-dihydroxyisovalerate → 2-oxoisovalerate + H₂O`
- 유전자: ilvD
- 반응식: `dhiv_c --> 2kiv_c + h2o_c` (comprehensive_missing_reactions.csv)

#### IPMI (Isopropylmalate isomerase)
- 반응식: `2-isopropylmalate ⇄ 3-isopropylmalate`
- 유전자: leuC, leuD
- 반응식: `2ippm_c <=> 3ippm_c` (comprehensive_missing_reactions.csv)

#### BCAT_VAL, BCAT_LEU (Branched-chain amino acid transaminase)
- BCAT_VAL: `2-oxoisovalerate + glutamate ⇄ valine + 2-oxoglutarate`
- BCAT_LEU: `2-oxoisocaproate + glutamate ⇄ leucine + 2-oxoglutarate`
- 유전자: ilvE, bcaT
- 반응식: 
  - BCAT_VAL: `2kiv_c + glu__L_c <=> akg_c + val__L_c`
  - BCAT_LEU: `4mop_c + glu__L_c <=> akg_c + leu__L_c` (comprehensive_missing_reactions.csv)

---

### 2. 프롤린(Proline) 합성

#### P5CS (Pyrroline-5-carboxylate Synthase)
- 반응식: `glutamate + ATP + NADPH → pyrroline-5-carboxylate + ADP + Pi + NADP⁺`
- 유전자: proB, proA

#### P5CD / P5CR (Pyrroline-5-carboxylate Reductase)
- 반응식: `pyrroline-5-carboxylate + NADPH → proline + NADP⁺`
- 유전자: proC
- 참고: P5CD와 P5CR은 같은 반응을 지칭함 (gap_filling_paper_method.py에 P5CD로 표기됨)

---

### 3. 조효소 밸런싱 (NAD/NADP)

#### THD2pp (NADP transhydrogenase)
- 반응식: `NADH + NADP⁺ + H⁺_periplasm → NAD⁺ + NADPH + H⁺_cytosol`
- Acetate 대사로 넘쳐나는 NADH를 생장(Anabolism)에 필요한 NADPH로 변환해주는 필수 효소
- 참고: comprehensive_missing_reactions.csv에는 직접적으로 없지만, NADTRHD나 NADH16_MQ 같은 transhydrogenase 관련 반응이 있을 수 있음
- 레퍼런스 모델에서 확인 필요

---

## 현재 누락된 반응 상태

### comprehensive_missing_reactions.csv에 포함된 반응

#### BCAA 관련 (LOW 우선순위)
- **BCAT_LEU**: `4mop_c + glu__L_c <=> akg_c + leu__L_c` (LOW)
- **BCAT_VAL**: `2kiv_c + glu__L_c <=> akg_c + val__L_c` (LOW)
- **IPMDH**: `3ippm_c + nad_c --> 4mop_c + co2_c + nadh_c` (MEDIUM, TCA Cycle 경로로 분류)
- **DHAD**: `dhiv_c --> 2kiv_c + h2o_c` (LOW)
- **IPMI**: `2ippm_c <=> 3ippm_c` (LOW)
- **KARI**: `alac__S_c + h_c + nadph_c --> dhiv_c + nadp_c` (LOW)

#### Proline 관련
- comprehensive_missing_reactions.csv에 **직접적으로 없음**
- gap_filling_paper_method.py에 P5CS, P5CD 추가 코드가 있음

#### Transhydrogenase 관련
- comprehensive_missing_reactions.csv에 THD2pp가 직접적으로 없음
- NADH16_MQ, MKRED, MQN8RD 등은 있지만 이들은 다른 기능

---

## 추가 확인 사항

### 1. 레퍼런스 모델 확인
- 레퍼런스 모델에서 ALS, IPMS, IPMD, P5CD/P5CR, THD2pp 등의 반응이 실제로 있는지 확인 필요

### 2. 반응 ID 확인
- IPMD vs IPMDH: 레퍼런스 모델에서 어떤 ID를 사용하는지 확인
- P5CD vs P5CR: 레퍼런스 모델에서 어떤 ID를 사용하는지 확인
- 2OXOVISO vs KARI: 실제 경로에 필요한 반응 확인

### 3. 논문의 19개 필수 반응
- 논문의 Materials and Methods 2.3 섹션의 정확한 19개 반응 목록 확인 필요
- 현재 파악된 반응들만으로는 19개를 채우지 못함

---

## 다음 단계

1. **레퍼런스 모델 확인**: 레퍼런스 모델에서 생합성 경로 관련 반응 확인
2. **누락 반응 추가**: 확인된 반응들을 모델에 추가
3. **FBA 테스트**: 추가 후 FBA 실행하여 성장률 확인
4. **ADP Deadlock 해결**: ADK1 경로 확인 및 필요시 추가
