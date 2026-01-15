# 레퍼런스 모델 vs 신규 모델 전체 비교 분석 결과

## 분석 개요

레퍼런스 모델(`model_YE0p5.xml`)과 신규 모델(`BaseModel.xml`)의 전체 비교 분석 결과입니다.

- **레퍼런스 모델 반응 수**: 2,198개
- **신규 모델 반응 수**: 2,092개
- **레퍼런스에만 있는 반응**: 114개

## 실제 FBA 결과 분석

레퍼런스 모델의 실제 FBA 플럭스 결과(`fba_flux_gradient_acid.csv`)를 기반으로 분석:

- **레퍼런스 모델에서 실제 사용된 반응**: 77개 (|flux| > 1e-6)
  - **신규 모델에 있음**: 68개
  - **신규 모델에 없음**: 9개 ⭐

## 누락된 반응 분류

### 실제 사용된 반응 중 누락 (9개) - HIGH 우선순위

이 9개 반응은 레퍼런스 모델의 FBA에서 실제로 사용되었습니다.

#### 1. **SUCDi** (Succinate dehydrogenase)
- **반응식**: `q8_c + succ_c --> fum_c + q8h2_c`
- **경로**: TCA Cycle
- **최대 플럭스**: 1.0
- **유전자**: sdhABCD (succinate dehydrogenase)

#### 2. **ACS_ADP** (Acetate-CoA ligase, ADP-forming)
- **반응식**: `ac_c + atp_c + coa_c <=> accoa_c + adp_c + pi_c`
- **경로**: Acetate Metabolism
- **최대 플럭스**: 0.94
- **유전자**: acsA

#### 3. **PEPCK_ATP** (PEP carboxykinase, ATP)
- **반응식**: `atp_c + oaa_c --> adp_c + co2_c + pep_c`
- **경로**: Gluconeogenesis / Glycolysis
- **최대 플럭스**: 0.079
- **유전자**: pckA
- **비고**: PEPCK (GTP 버전)는 레퍼런스 모델에도 있지만 실제 FBA에서 사용되지 않음 (플럭스 0)

#### 4-9. **Transport 반응들**
- **T_o2_e_to_o2_c**: O2 transport (최대 플럭스: 2.0)
- **EX_hco3_e**: Bicarbonate exchange (최대 플럭스: 1.0)
- **T_hco3_e_to_c**: Bicarbonate transport (최대 플럭스: 1.0)
- **ACtexi**: Acetate transport (최대 플럭스: 1.0)
- **T_nh4_e_to_nh4_c**: NH4 transport (최대 플럭스: 0.10)
- **T_fe3_e_to_fe3_c**: Fe3+ transport (최대 플럭스: 0.018)

### 실제 사용되지 않은 반응 중 누락 (105개)

레퍼런스 모델에는 있지만 실제 FBA에서 사용되지 않은 반응들입니다.

#### 우선순위 분류

- **HIGH**: 0개
- **MEDIUM**: 82개
  - Transport: 69개
  - Exchange: 11개
  - Other: 2개
- **LOW**: 15개

#### 주요 카테고리

1. **Transport 반응 (69개)**
   - 대부분 아미노산, 이온, 보조인자 운송 반응
   - 선택적으로 필요할 수 있음

2. **Exchange 반응 (11개)**
   - 보조인자 교환 (biotin, CoA, folate, lipoate, menaquinone, nicotinamide, pantothenate, ubiquinone, riboflavin, thiamine 등)

3. **대사 경로 반응 (약 25개)**
   - PEPCK (GTP 버전) - 사용 안 됨
   - PC (Pyruvate carboxylase) - 사용 안 됨
   - PPDK (Pyruvate phosphate dikinase) - 사용 안 됨
   - ACKr (Acetate kinase reverse) - 사용 안 됨
   - Menaquinone 관련 반응들 (MKRED, MQN8RD, MQN8red 등) - 사용 안 됨
   - FRD7 (Fumarate reductase) - 사용 안 됨
   - 기타 반응들

## 결론 및 권장사항

### 최우선 추가 반응 (실제 사용된 9개)

이 9개 반응은 레퍼런스 모델의 FBA에서 실제로 사용되었으므로 **최우선으로 추가**해야 합니다:

1. **SUCDi** - TCA 사이클/ETC 연결
2. **ACS_ADP** - Acetate 대사
3. **PEPCK_ATP** - Gluconeogenesis
4. **ACtexi** - Acetate transport
5. **EX_hco3_e** + **T_hco3_e_to_c** - Bicarbonate 공급/운반
6. **T_o2_e_to_o2_c** - O2 운반
7. **T_nh4_e_to_nh4_c** - NH4 운반
8. **T_fe3_e_to_fe3_c** - Fe3+ 운반

### 추가 확인 필요 (실제 사용 안 된 반응)

실제 FBA에서 사용되지 않았지만, 신규 모델의 구조나 조건에 따라 필요할 수 있는 반응들:

- Transport 반응들 (선택적)
- Exchange 반응들 (보조인자 공급용)
- 대체 경로 반응들 (PC, PPDK, PEPCK GTP 등)

### 테스트 결과

9개 반응을 모두 추가했을 때 FBA 테스트 결과:
- **FBA 상태**: FAIL (성장률 = 0.0)
- **원인**: 추가 반응이 더 필요하거나, 대사물질 연결 문제

### 다음 단계

1. **9개 반응 추가 후 테스트**
   - 신규 모델에 9개 반응 추가
   - FBA 실행 및 결과 확인

2. **추가 진단 필요 시**
   - Blocked reactions 분석
   - Biomass 구성 요소 생산 경로 확인
   - 대사물질 연결 확인

3. **유전자 확인**
   - 9개 반응의 유전자 확인
   - 신규 분리 미생물 genome에서 검색

## 생성된 파일

1. **`comprehensive_missing_reactions.csv`** - 레퍼런스에만 있는 전체 반응 114개
2. **`missing_active_reactions_detailed.csv`** - 실제 사용된 반응 중 누락 9개 (상세)
3. **`all_missing_active_reactions.csv`** - 실제 사용된 반응 중 누락 9개 (간단)
4. **`missing_active_reactions.csv`** - 실제 사용된 반응 중 누락 9개 (경로 분류)
