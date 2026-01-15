# 레퍼런스 모델에 있지만 신규 모델에 없는 대사 경로 반응 리스트
## (Transport, Exchange 및 Pseudo 반응 제외)

총 **25개** 실제 대사 반응
- Pseudo 반응 2개 제외 (SUPPLY_accoa_c, SUPPLY_oaa_c)

---

## 요약

- **총 실제 대사 반응**: 25개
- **Pseudo 반응 (제외)**: 2개
  - SUPPLY_accoa_c: ` --> accoa_c`
  - SUPPLY_oaa_c: ` --> oaa_c`

### 경로별 분류
- **Other**: 11개
- **Cofactor Metabolism**: 4개
- **Gluconeogenesis / Glycolysis**: 4개
- **ETC / Electron Transport**: 3개
- **TCA Cycle**: 3개
- **Acetate Metabolism**: 2개

### 우선순위별 분류
- **HIGH**: 3개 (실제 사용됨)
- **MEDIUM**: 12개
- **LOW**: 10개

---

## HIGH 우선순위 반응 (3개) - 실제 사용됨

### 1. ACS_ADP (Acetate-CoA ligase, ADP-forming)
- **반응 ID**: ACS_ADP
- **반응식**: `ac_c + atp_c + coa_c <=> accoa_c + adp_c + pi_c`
- **경로**: Acetate Metabolism
- **실제 사용됨**: Yes (플럭스: 0.943448)
- **가능한 유전자**: acsA (Acetate-CoA synthetase)

### 2. PEPCK_ATP (PEP carboxykinase, ATP)
- **반응 ID**: PEPCK_ATP
- **반응식**: `atp_c + oaa_c --> adp_c + co2_c + pep_c`
- **경로**: Gluconeogenesis / Glycolysis
- **실제 사용됨**: Yes (플럭스: 0.078621)
- **가능한 유전자**: pckA (Phosphoenolpyruvate carboxykinase)

### 3. SUCDi (Succinate dehydrogenase)
- **반응 ID**: SUCDi
- **반응식**: `q8_c + succ_c --> fum_c + q8h2_c`
- **경로**: TCA Cycle
- **실제 사용됨**: Yes (플럭스: 1.000000)
- **가능한 유전자**: sdhABCD (Succinate dehydrogenase complex)

---

## MEDIUM 우선순위 반응 (12개)

### Acetate Metabolism (1개)

#### 4. ACKr (Acetate kinase, reverse)
- **반응식**: `ac_c + atp_c <=> actp_c + adp_c`
- **가능한 유전자**: ackA (Acetate kinase)

### ETC / Electron Transport (3개)

#### 5. MKRED (Menaquinone-8 reduction by NADH)
- **반응식**: `h_c + mqn8_c + nadh_c --> mql8_c + nad_c`
- **가능한 유전자**: menA, menB, menC, menD, menE (Menaquinone biosynthesis)

#### 6. MQN8RD (Menaquinone-8 reduction by NADH, alternative)
- **반응식**: `h_c + mqn8_c + nadh_c --> mql8_c + nad_c`
- **가능한 유전자**: menaquinone 관련

#### 7. MQN8red (Menaquinone-8 NADH dehydrogenase, lumped)
- **반응식**: `h_c + mqn8_c + nadh_c --> mql8_c + nad_c`
- **가능한 유전자**: menaquinone 관련

### Gluconeogenesis / Glycolysis (3개)

#### 8. PC (Pyruvate carboxylase)
- **반응식**: `atp_c + hco3_c + pyr_c --> adp_c + h_c + oaa_c + pi_c`
- **가능한 유전자**: pyc/pycA (Pyruvate carboxylase)

#### 9. PEPCK (Phosphoenolpyruvate carboxykinase, GTP)
- **반응식**: `gtp_c + oaa_c --> co2_c + gdp_c + pep_c`
- **가능한 유전자**: pckA (Phosphoenolpyruvate carboxykinase)
- **비고**: PEPCK_ATP이 실제 사용됨

#### 10. PPDK (Pyruvate phosphate dikinase)
- **반응식**: `atp_c + pi_c + pyr_c --> amp_c + pep_c + ppi_c`
- **가능한 유전자**: ppdK (Pyruvate phosphate dikinase)

### TCA Cycle (2개)

#### 11. FRD7 (Fumarate reductase)
- **반응식**: `fum_c + q8h2_c --> q8_c + succ_c`
- **가능한 유전자**: frdABCD (Fumarate reductase complex)

#### 12. IPMDH (Isopropylmalate dehydrogenase)
- **반응식**: `3ippm_c + nad_c --> 4mop_c + co2_c + nadh_c`
- **가능한 유전자**: leuB (Isopropylmalate dehydrogenase)

---

## LOW 우선순위 반응 (10개)

### Cofactor Metabolism (3개)

#### 13. FADS (FAD synthetase)
- **반응식**: `atp_c + fmn_c --> fad_c + ppi_c`
- **가능한 유전자**: FAD synthase 관련

#### 14. MQN8r_NADH (Menaquinone-8 reduction by NADH, alternative)
- **반응식**: `2.0 h_c + mqn8_c + nadh_c --> mql8_c + nad_c`
- **가능한 유전자**: menaquinone 관련

#### 15. NADH16_MQ (NDH-1 to MQN8)
- **반응식**: `4.0 h_c + mqn8_c + nadh_c --> 3.0 h_p + mqn8h2_c + nad_c`
- **가능한 유전자**: NADH dehydrogenase 관련

### Other (7개)

#### 16. BCAT_LEU (Branched-chain amino acid aminotransferase, Leucine)
- **반응식**: `4mop_c + glu__L_c <=> akg_c + leu__L_c`
- **가능한 유전자**: BCAT (branched-chain aminotransferase)

#### 17. BCAT_VAL (Branched-chain amino acid aminotransferase, Valine)
- **반응식**: `2kiv_c + glu__L_c <=> akg_c + val__L_c`
- **가능한 유전자**: BCAT (branched-chain aminotransferase)

#### 18. CA (Carbonic anhydrase)
- **반응식**: `co2_c + h2o_c <=> h_c + hco3_c`
- **가능한 유전자**: can (carbonic anhydrase)

#### 19. DHAD (Dihydroxyacid dehydratase)
- **반응식**: `dhiv_c --> 2kiv_c + h2o_c`
- **가능한 유전자**: ilvD (dihydroxyacid dehydratase)

#### 20. IPMI (Isopropylmalate isomerase)
- **반응식**: `2ippm_c <=> 3ippm_c`
- **가능한 유전자**: leuC, leuD (isopropylmalate isomerase)

#### 21. KARI (Ketol-acid reductoisomerase)
- **반응식**: `alac__S_c + h_c + nadph_c --> dhiv_c + nadp_c`
- **가능한 유전자**: ilvC (ketol-acid reductoisomerase)

#### 22. MKOX (Menaquinol-8 oxidation, pseudo)
- **반응식**: `mql8_c + nad_c --> h_c + mqn8_c + nadh_c`
- **비고**: Pseudo reaction (유전자 없음)

#### 23. PDXK (Pyridoxal kinase)
- **반응식**: `atp_c + pydxn_c --> adp_c + h_c + pydx5p_c`
- **가능한 유전자**: pdxK (pyridoxal kinase)

#### 24. RIBFLVKin (Riboflavin kinase)
- **반응식**: `atp_c + ribflv_c --> adp_c + fmn_c + h_c`
- **가능한 유전자**: riboflavin kinase 관련

#### 25. THIK (Thiamine kinase)
- **반응식**: `atp_c + thm_c --> adp_c + h_c + thmpp_c`
- **가능한 유전자**: thiK (thiamine kinase)

---

## 제외된 Pseudo 반응 (유전자 없음)

### 1. SUPPLY_accoa_c
- **반응식**: ` --> accoa_c`
- **비고**: 공급 반응 (pseudo reaction) - 유전자 없음, 추가 불필요

### 2. SUPPLY_oaa_c
- **반응식**: ` --> oaa_c`
- **비고**: 공급 반응 (pseudo reaction) - 유전자 없음, 추가 불필요

---

## 주요 확인 대상 유전자 (우선순위별)

### HIGH 우선순위 (필수 확인) - 3개

1. **acsA** - ACS_ADP (Acetate-CoA synthetase)
2. **pckA** - PEPCK_ATP (Phosphoenolpyruvate carboxykinase)
3. **sdhABCD** - SUCDi (Succinate dehydrogenase complex)

### MEDIUM 우선순위 (추가 확인) - 12개

4. **ackA** - ACKr (Acetate kinase)
5. **menA, menB, menC, menD, menE** - Menaquinone 관련 (MKRED, MQN8RD, MQN8red)
6. **pyc/pycA** - PC (Pyruvate carboxylase)
7. **ppdK** - PPDK (Pyruvate phosphate dikinase)
8. **frdABCD** - FRD7 (Fumarate reductase complex)
9. **leuB** - IPMDH (Isopropylmalate dehydrogenase)

### LOW 우선순위 (선택 확인) - 10개

10. **ilvD** - DHAD (Dihydroxyacid dehydratase)
11. **leuC, leuD** - IPMI (Isopropylmalate isomerase)
12. **ilvC** - KARI (Ketol-acid reductoisomerase)
13. **pdxK** - PDXK (Pyridoxal kinase)
14. **thiK** - THIK (Thiamine kinase)
15. **can** - CA (Carbonic anhydrase)
16. 기타 cofactor 관련 유전자

---

## 다음 단계

이 **25개 실제 대사 반응**에 해당하는 유전자를 신규 분리 미생물의 whole genome DNA에서 찾아야 합니다.

특히 **HIGH 우선순위 3개 반응** (ACS_ADP, PEPCK_ATP, SUCDi)은 레퍼런스 모델의 FBA에서 실제로 사용되었으므로, 신규 분리 미생물에서도 이 반응들이 작동하려면 해당 유전자들이 필수입니다.

**참고**: SUPPLY_accoa_c와 SUPPLY_oaa_c는 pseudo reaction이므로 유전자를 찾을 필요가 없습니다. 이들은 모델링을 위한 공급 반응일 뿐 실제 대사 반응이 아닙니다.
