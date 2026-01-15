# 최종 모델 요약

## 모델 파일
**BaseModel_final_cleaned.xml**

## 최종 성장률
- **ATPM=0**: 0.396812
- **상태**: optimal

## 제거된 반응 (아티팩트)

### 1. Acetate 관련
- **ACt2rpp**: Acetate 재순환 방지

### 2. Acetyl-CoA 관련
- **ACS_ADP, SUCOAACTr**: 비현실적 경로

### 3. PEPCK 관련
- **PEPCK_ATP**: PEPCK 제거 (maeB로 대체)

### 4. 에너지 회수 루프
- **ACCOAL, APAT_1, PACPT_1**: 에너지 회수 루프

### 5. PEP 루프
- **GALpts, GALt2, A6PAG**: 비현실적 PEP 순환

### 6. 레독스 셔틀
- **ACOAD2, ACOAD2f, SUCD**: FADH2 → NADH 변환 셔틀
- **SUCDi 사용** (Q8 기반, 정상)

### 7. 기타
- **DM_coa_c**: CoA 생합성 경로 활성화
- **sink_4hba_c**: 비현실적 sink

## 추가된 반응

### 1. BCAA 합성 (6개)
- KARI, DHAD, IPMI, IPMDH, BCAT_VAL, BCAT_LEU

### 2. 이온 수송 (3개)
- T_cl_e_to_cl_c, T_cu2_e_to_cu2_c, T_cobalt2_e_to_cobalt2_c

### 3. NAD/NADP (3개)
- T_nac_e_to_nac_c, EX_nac_e, EX_ncam_e

### 4. 기타
- **maeB**: NADP-dependent malic enzyme
- **EX_pnto__R_e, T_pnto__R_e_to_c**: 판토텐산 (yeast extract 가정)
- **EX_4hba_e, T_4hba_e_to_c**: 4hba 배출

## 수정된 반응
- **PPA_1pp**: 비활성화 (H+-translocating PPase)

## 유지된 반응 (작은 플럭스)
- **MQN8t**: 메나퀴논 생합성 경로 불완전, 플럭스 매우 작음 (0.000040)

## 배지 가정
- **EX_pnto__R_e**: 판토텐산 exchange (yeast extract 가정)
- **EX_nac_e, EX_ncam_e**: NAD 전구체 exchange (yeast extract 가정)
- 실제 배지 조건에 맞게 조정 가능

## 주요 Flux (ATPM=0)

### 탄소 전환
- **EX_ac_e**: -19.0
- **ACS**: 19.06
- **CS**: 13.47 (TCA cycle)
- **ICL**: 5.16 (Glyoxylate shunt)
- **MALS**: 5.16 (Glyoxylate shunt)
- **SUCDi**: 13.02 (Q8 기반)

### 에너지
- **ATPS4rpp**: 64.98 (ATP 생성)
- **Growth**: 21.48 (ATP 소비)
- **ACS**: 19.06 (ATP 소비)

## Export된 파일
- **flux_analysis_final_complete.csv**: 주요 반응 25개
