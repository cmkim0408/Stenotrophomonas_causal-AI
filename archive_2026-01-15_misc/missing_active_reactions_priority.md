# 레퍼런스 모델 실제 FBA 결과 기반 누락 반응 분석

## 분석 개요

레퍼런스 모델(`model_YE0p5.xml`)의 실제 FBA 플럭스 결과(`fba_flux_gradient_acid.csv`)를 기반으로 분석한 결과입니다.

- **레퍼런스 모델에서 실제 사용된 반응**: 77개 (|flux| > 1e-6)
- **그 중 신규 모델에 없는 반응**: 9개

## 누락된 실제 사용 반응 (9개)

### HIGH 우선순위 (핵심 대사 경로)

#### 1. **ACS_ADP** (Acetate-CoA ligase, ADP-forming)
- **반응식**: `ac_c + atp_c + coa_c <=> accoa_c + adp_c + pi_c`
- **경로**: Acetate Metabolism
- **중요성**: Acetate에서 Acetyl-CoA 생산 핵심 반응
- **유전자**: 확인 필요 (acsA 관련)
- **비고**: 레퍼런스 모델에서 실제로 사용된 Acetate 동화 경로

#### 2. **PEPCK_ATP** (PEP carboxykinase, ATP)
- **반응식**: `atp_c + oaa_c --> adp_c + co2_c + pep_c`
- **경로**: Gluconeogenesis
- **중요성**: OAA에서 PEP 생산 (Gluconeogenesis 핵심)
- **유전자**: pckA (PEP carboxykinase)
- **비고**: PEPCK (GTP 버전)는 레퍼런스 모델에도 있지만 실제 FBA에서는 사용되지 않음 (플럭스 0)

#### 3. **SUCDi** (Succinate dehydrogenase)
- **반응식**: `q8_c + succ_c --> fum_c + q8h2_c`
- **경로**: TCA Cycle / ETC
- **중요성**: TCA 사이클, 전자 전달 사슬 연결
- **유전자**: sdhABCD (succinate dehydrogenase)
- **비고**: TCA 사이클에서 fumarate 생산 및 전자 전달

### MEDIUM 우선순위 (Transport 및 Exchange)

#### 4. **ACtexi** (Acetate transport)
- **반응식**: `ac_e <=> ac_c`
- **경로**: Transport
- **중요성**: Acetate 세포 내 유입
- **비고**: Acetate 운반체

#### 5. **EX_hco3_e** (Bicarbonate exchange)
- **반응식**: `hco3_e <=>`
- **경로**: Exchange
- **중요성**: HCO3 교환 (CO2 고정 관련)
- **비고**: Bicarbonate 공급

#### 6. **T_hco3_e_to_c** (Bicarbonate transport)
- **반응식**: `hco3_e <=> hco3_c`
- **경로**: Transport
- **중요성**: HCO3 세포 내 유입
- **비고**: EX_hco3_e와 함께 작동

#### 7. **T_o2_e_to_o2_c** (Oxygen transport)
- **반응식**: `o2_e <=> o2_c`
- **경로**: Transport
- **중요성**: O2 세포 내 유입
- **비고**: 호기성 대사에 필요

#### 8. **T_nh4_e_to_nh4_c** (Ammonium transport)
- **반응식**: `nh4_e <=> nh4_c`
- **경로**: Transport
- **중요성**: NH4 세포 내 유입
- **비고**: 질소 공급

#### 9. **T_fe3_e_to_fe3_c** (Fe3+ transport)
- **반응식**: `fe3_e <=> fe3_c`
- **경로**: Transport
- **중요성**: Fe3+ 세포 내 유입
- **비고**: 미량원소 공급

## 중요 발견 사항

1. **PEPCK vs PEPCK_ATP**
   - 레퍼런스 모델에는 PEPCK (GTP 버전)와 PEPCK_ATP (ATP 버전) 둘 다 존재
   - 실제 FBA 결과: PEPCK (GTP) 플럭스 = 0, PEPCK_ATP 플럭스 > 0
   - **결론**: PEPCK_ATP만 실제로 사용됨, PEPCK (GTP)는 필요 없음

2. **ACS_ADP**
   - 레퍼런스 모델에서 Acetate 동화에 실제로 사용된 반응
   - 신규 모델에 없으면 Acetate 대사 불가

3. **SUCDi**
   - TCA 사이클과 전자 전달 사슬 연결에 필수
   - 신규 모델에 없으면 TCA 사이클이 불완전

4. **Transport 반응**
   - 기본적인 영양소 운반체들이 누락
   - 특히 Acetate, HCO3, O2, NH4, Fe3+ 운반체

## 우선순위별 권장 사항

### 최우선 추가 반응 (HIGH)

1. **ACS_ADP** - Acetate 대사 핵심
2. **PEPCK_ATP** - Gluconeogenesis 핵심
3. **SUCDi** - TCA 사이클/ETC 연결

### 다음 단계 추가 반응 (MEDIUM)

4. **ACtexi** - Acetate 운반
5. **EX_hco3_e** + **T_hco3_e_to_c** - HCO3 공급/운반 (PEPCK_ATP와 함께 필요)
6. **T_o2_e_to_o2_c** - O2 운반 (호기성 대사)
7. **T_nh4_e_to_nh4_c** - NH4 운반 (질소 공급)
8. **T_fe3_e_to_fe3_c** - Fe3+ 운반 (미량원소)

## 유전자 확인 작업

다음 효소들의 유전자를 신규 분리 미생물 genome에서 확인:

1. **ACS_ADP** → acsA (Acetate-CoA ligase)
2. **PEPCK_ATP** → pckA (PEP carboxykinase)
3. **SUCDi** → sdhABCD (Succinate dehydrogenase complex)
4. **ACtexi** → acetate transport 유전자
5. **T_hco3_e_to_c** → bicarbonate transport 유전자
6. **T_o2_e_to_o2_c** → O2 transport 관련
7. **T_nh4_e_to_nh4_c** → ammonium transport 유전자
8. **T_fe3_e_to_fe3_c** → Fe3+ transport 유전자
