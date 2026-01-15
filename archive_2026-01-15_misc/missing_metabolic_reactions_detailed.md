# 레퍼런스 모델에 있지만 신규 모델에 없는 대사 경로 반응 리스트
## (Transport 및 Exchange 제외)

총 **27개** 반응

---

## 실제 대사 반응: 25개

### Pseudo 반응 (제외): 2개

- **SUPPLY_accoa_c**:  --> accoa_c (Pseudo reaction)
- **SUPPLY_oaa_c**:  --> oaa_c (Pseudo reaction)

---

## HIGH 우선순위 반응 (실제 사용됨)

### 1. ACS_ADP
- **반응식**: `ac_c + atp_c + coa_c <=> accoa_c + adp_c + pi_c`
- **경로**: Acetate Metabolism
- **실제 사용됨**: Yes (플럭스: 0.943448)

### 2. PEPCK_ATP
- **이름**: PEP carboxykinase (ATP)
- **반응식**: `atp_c + oaa_c --> adp_c + co2_c + pep_c`
- **경로**: Gluconeogenesis / Glycolysis
- **실제 사용됨**: Yes (플럭스: 0.078621)

### 3. SUCDi
- **이름**: Succinate dehydrogenase
- **반응식**: `q8_c + succ_c --> fum_c + q8h2_c`
- **경로**: TCA Cycle
- **실제 사용됨**: Yes (플럭스: 1.000000)

## MEDIUM 우선순위 반응

### 1. ACKr
- **반응식**: `ac_c + atp_c <=> actp_c + adp_c`
- **경로**: Acetate Metabolism

### 2. MKRED
- **이름**: menaquinone-8 reduction by NADH
- **반응식**: `h_c + mqn8_c + nadh_c --> mql8_c + nad_c`
- **경로**: ETC / Electron Transport

### 3. MQN8RD
- **이름**: menaquinone-8 reduction by NADH
- **반응식**: `h_c + mqn8_c + nadh_c --> mql8_c + nad_c`
- **경로**: ETC / Electron Transport

### 4. MQN8red
- **이름**: Menaquinone-8 NADH dehydrogenase (lumped)
- **반응식**: `h_c + mqn8_c + nadh_c --> mql8_c + nad_c`
- **경로**: ETC / Electron Transport

### 5. PC
- **이름**: pyruvate carboxylase
- **반응식**: `atp_c + hco3_c + pyr_c --> adp_c + h_c + oaa_c + pi_c`
- **경로**: Gluconeogenesis / Glycolysis

### 6. PEPCK
- **이름**: phosphoenolpyruvate carboxykinase (GTP)
- **반응식**: `gtp_c + oaa_c --> co2_c + gdp_c + pep_c`
- **경로**: Gluconeogenesis / Glycolysis

### 7. PPDK
- **이름**: Pyruvate phosphate dikinase
- **반응식**: `atp_c + pi_c + pyr_c --> amp_c + pep_c + ppi_c`
- **경로**: Gluconeogenesis / Glycolysis

### 8. FRD7
- **이름**: Fumarate reductase
- **반응식**: `fum_c + q8h2_c --> q8_c + succ_c`
- **경로**: TCA Cycle

### 9. IPMDH
- **반응식**: `3ippm_c + nad_c --> 4mop_c + co2_c + nadh_c`
- **경로**: TCA Cycle

## LOW 우선순위 반응

### 1. FADS
- **반응식**: `atp_c + fmn_c --> fad_c + ppi_c`
- **경로**: Cofactor Metabolism

### 2. MQN8r_NADH
- **반응식**: `2.0 h_c + mqn8_c + nadh_c --> mql8_c + nad_c`
- **경로**: Cofactor Metabolism

### 3. NADH16_MQ
- **이름**: NDH-1 to MQN8
- **반응식**: `4.0 h_c + mqn8_c + nadh_c --> 3.0 h_p + mqn8h2_c + nad_c`
- **경로**: Cofactor Metabolism

### 4. BCAT_LEU
- **반응식**: `4mop_c + glu__L_c <=> akg_c + leu__L_c`
- **경로**: Other

### 5. BCAT_VAL
- **반응식**: `2kiv_c + glu__L_c <=> akg_c + val__L_c`
- **경로**: Other

### 6. CA
- **이름**: carbonic anhydrase
- **반응식**: `co2_c + h2o_c <=> h_c + hco3_c`
- **경로**: Other

### 7. DHAD
- **반응식**: `dhiv_c --> 2kiv_c + h2o_c`
- **경로**: Other

### 8. IPMI
- **반응식**: `2ippm_c <=> 3ippm_c`
- **경로**: Other

### 9. KARI
- **반응식**: `alac__S_c + h_c + nadph_c --> dhiv_c + nadp_c`
- **경로**: Other

### 10. MKOX
- **이름**: menaquinol-8 oxidation (pseudo)
- **반응식**: `mql8_c + nad_c --> h_c + mqn8_c + nadh_c`
- **경로**: Other

### 11. PDXK
- **반응식**: `atp_c + pydxn_c --> adp_c + h_c + pydx5p_c`
- **경로**: Other

### 12. RIBFLVKin
- **반응식**: `atp_c + ribflv_c --> adp_c + fmn_c + h_c`
- **경로**: Other

### 13. THIK
- **반응식**: `atp_c + thm_c --> adp_c + h_c + thmpp_c`
- **경로**: Other

