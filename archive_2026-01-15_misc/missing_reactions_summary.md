# 레퍼런스 모델과 신규 모델 비교 분석 결과

## 요약

- **레퍼런스 모델 반응 수**: 2,198개
- **신규 모델 반응 수**: 2,092개  
- **레퍼런스에만 있는 반응**: 114개
- **신규 모델에만 있는 반응**: 8개
- **레퍼런스에만 있는 유전자**: 0개 (공통 유전자: 931개)

## 누락된 반응 우선순위 분류

### HIGH 우선순위 (FBA 성공에 핵심) - 5개

1. **PC** (Pyruvate carboxylase)
   - 반응식: ATP + HCO3 + Pyruvate → ADP + H + OAA + Pi
   - 설명: Pyruvate에서 OAA 생산 (Gluconeogenesis 핵심)
   - 이유: PEP 생산에 필요, OAA 생산 경로

2. **PPDK** (Pyruvate phosphate dikinase)
   - 반응식: ATP + Pi + Pyruvate → AMP + PEP + PPI
   - 설명: Pyruvate에서 PEP 생산 (대체 경로)
   - 이유: PEP 생산 대체 경로

3. **PEPCK** (PEP carboxykinase, GTP)
   - 반응식: GTP + OAA → CO2 + GDP + PEP
   - 설명: OAA에서 PEP 생산 (GTP 사용)
   - 이유: PEP 생산 핵심 반응

4. **PEPCK_ATP** (PEP carboxykinase, ATP)
   - 반응식: ATP + OAA → ADP + CO2 + PEP
   - 설명: OAA에서 PEP 생산 (ATP 사용)
   - 이유: PEP 생산 대체 반응

5. **ACS_ADP** (Acetate-CoA ligase, ADP-forming)
   - 반응식: Acetate + ATP + CoA ⇄ AcCoA + ADP + Pi
   - 설명: Acetate에서 AcCoA 생산
   - 이유: Acetate에서 AcCoA 생산 핵심 반응

### MEDIUM 우선순위 (FBA 성공에 중요) - 11개

1. **SUCDi** (Succinate dehydrogenase)
   - 반응식: Succinate + Q8 → Fumarate + Q8H2
   - 설명: TCA 사이클, 전자 전달 사슬

2. **FRD7** (Fumarate reductase)
   - 반응식: Fumarate + Q8H2 → Succinate + Q8
   - 설명: 무산소 조건에서 사용 가능

3. **MKRED** (Menaquinone-8 reduction by NADH)
   - 반응식: H + MQN8 + NADH → MQL8 + NAD
   - 설명: 전자 전달 사슬 관련

4. **MQN8RD** (Menaquinone-8 reduction by NADH, alternative)
   - 반응식: H + MQN8 + NADH → MQL8 + NAD
   - 설명: 전자 전달 사슬 관련

5. **MQN8r_NADH** (Menaquinone-8 reduction by NADH)
   - 반응식: 2 H + MQN8 + NADH → MQL8 + NAD
   - 설명: 전자 전달 사슬 관련

6. **MQN8red** (Menaquinone-8 NADH dehydrogenase, lumped)
   - 반응식: H + MQN8 + NADH → MQL8 + NAD
   - 설명: 전자 전달 사슬 관련

7. **NADH16_MQ** (NDH-1 to MQN8)
   - 반응식: 4 H + MQN8 + NADH → 3 H_p + MQN8H2 + NAD
   - 설명: Complex I 관련 전자 전달

8. **MKOX** (Menaquinol-8 oxidation)
   - 반응식: MQL8 + NAD → H + MQN8 + NADH
   - 설명: 전자 전달 사슬 관련

9. **FADS** (FAD synthetase)
   - 반응식: ATP + FMN → FAD + PPI
   - 설명: FAD 보조인자 합성

10. **RIBFLVKin** (Riboflavin kinase)
    - 반응식: ATP + Riboflavin → ADP + FMN + H
    - 설명: FMN 생합성 경로

11. **CA** (Carbonic anhydrase)
    - 반응식: CO2 + H2O ⇄ H + HCO3
    - 설명: HCO3 생산/소비 균형

### LOW 우선순위 - 4개

1. **BCAT_LEU** (Branched-chain amino acid transaminase, Leu)
2. **BCAT_VAL** (Branched-chain amino acid transaminase, Val)
3. **SUPPLY_accoa_c** (AcCoA supply reaction, bootstrap)
4. **SUPPLY_oaa_c** (OAA supply reaction, bootstrap)

### 기타 (Exchange/Transport) - 94개

- Exchange 반응: 12개 (주로 보조인자 교환)
- Transport 반응: 75개 (아미노산, 이온, 보조인자 운송)
- 기타 반응: 7개

## 핵심 분석 결과

### 1. Gluconeogenesis 경로 문제

레퍼런스 모델에는 PEP 생산을 위한 여러 경로가 있습니다:
- **PC**: Pyruvate → OAA
- **PEPCK/PEPCK_ATP**: OAA → PEP
- **PPDK**: Pyruvate → PEP (대체 경로)

신규 모델에는 이러한 반응들이 누락되어 있어 PEP 생산이 불가능합니다.

### 2. Acetyl-CoA 생산 경로

- **ACS_ADP**: Acetate + ATP + CoA → AcCoA (누락)
- 레퍼런스 모델에는 Acetate에서 AcCoA 생산을 위한 반응이 있으나 신규 모델에는 없습니다.

### 3. 전자 전달 사슬 (ETC)

레퍼런스 모델에는 Menaquinone 관련 반응들이 다수 존재:
- Menaquinone reduction (MKRED, MQN8RD, MQN8r_NADH, MQN8red, NADH16_MQ)
- Menaquinol oxidation (MKOX)
- Succinate dehydrogenase (SUCDi)
- Fumarate reductase (FRD7)

신규 모델에는 이러한 반응들이 누락되어 전자 전달 사슬이 불완전합니다.

### 4. 보조인자 합성

- **FADS**: FMN → FAD
- **RIBFLVKin**: Riboflavin → FMN

이러한 보조인자 생합성 경로가 누락되어 있습니다.

## 권장 사항

### 우선 추가해야 할 반응 (HIGH 우선순위)

1. **PEPCK** 또는 **PEPCK_ATP**: OAA → PEP 생산
2. **PC**: Pyruvate → OAA 생산 (PEPCK과 함께 사용)
3. **ACS_ADP**: Acetate → AcCoA 생산
4. **PPDK**: Pyruvate → PEP 생산 (대체 경로, 선택적)

### 다음 단계로 추가할 반응 (MEDIUM 우선순위)

1. **SUCDi**: Succinate dehydrogenase (TCA 사이클)
2. **MKRED** 또는 **NADH16_MQ**: Menaquinone 관련 전자 전달
3. **FADS**: FAD 합성
4. **CA**: Carbonic anhydrase (HCO3 균형)

### 유전자 확인 작업

다음 단계로 신규 분리 미생물의 genome에서 다음 효소들의 유전자를 확인해야 합니다:

1. **PEPCK** (PEP carboxykinase) - pckA/pck
2. **PC** (Pyruvate carboxylase) - pyc/pycA
3. **ACS** (Acetate-CoA ligase) - acsA
4. **PPDK** (Pyruvate phosphate dikinase) - ppdK
5. **SUCDi** (Succinate dehydrogenase) - sdhABCD
6. **Menaquinone 관련 유전자** - menA, menB, menC, menD, menE 등
7. **FADS** (FAD synthetase) - ribF
8. **RIBFLVKin** (Riboflavin kinase) - ribF
