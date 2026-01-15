# 레퍼런스 모델에 있지만 신규 모델에 없는 대사 경로 반응 리스트
## (Transport 및 Exchange 제외)

총 **27개** 반응

---

## HIGH 우선순위 (3개)

### 1. SUCDi (Succinate dehydrogenase)
- **반응식**: `q8_c + succ_c --> fum_c + q8h2_c`
- **경로**: TCA Cycle
- **실제 사용됨**: Yes (플럭스: 1.0)
- **유전자**: sdhABCD (succinate dehydrogenase)

### 2. ACS_ADP (Acetate-CoA ligase, ADP-forming)
- **반응식**: `ac_c + atp_c + coa_c <=> accoa_c + adp_c + pi_c`
- **경로**: Acetate Metabolism
- **실제 사용됨**: Yes (플럭스: 0.94)
- **유전자**: acsA (Acetate-CoA ligase)

### 3. PEPCK_ATP (PEP carboxykinase, ATP)
- **반응식**: `atp_c + oaa_c --> adp_c + co2_c + pep_c`
- **경로**: Gluconeogenesis / Glycolysis
- **실제 사용됨**: Yes (플럭스: 0.079)
- **유전자**: pckA (PEP carboxykinase)
- **비고**: PEPCK (GTP 버전)는 레퍼런스 모델에도 있지만 실제 FBA에서 사용되지 않음

---

## MEDIUM 우선순위 (24개)

### Gluconeogenesis / Glycolysis (3개)

4. **PC** (Pyruvate carboxylase)
   - 반응식: `atp_c + hco3_c + pyr_c --> adp_c + h_c + oaa_c + pi_c`
   - 실제 사용됨: No
   - 유전자: pyc/pycA (Pyruvate carboxylase)

5. **PEPCK** (PEP carboxykinase, GTP)
   - 반응식: `gtp_c + oaa_c --> co2_c + gdp_c + pep_c`
   - 실제 사용됨: No (PEPCK_ATP 사용됨)
   - 유전자: pckA (PEP carboxykinase)

6. **PPDK** (Pyruvate phosphate dikinase)
   - 반응식: `atp_c + pi_c + pyr_c --> amp_c + pep_c + ppi_c`
   - 실제 사용됨: No
   - 유전자: ppdK (Pyruvate phosphate dikinase)

### TCA Cycle (2개)

7. **FRD7** (Fumarate reductase)
   - 반응식: `fum_c + q8h2_c --> q8_c + succ_c`
   - 실제 사용됨: No
   - 유전자: frdABCD (Fumarate reductase)

8. **IPMDH** (Isopropylmalate dehydrogenase)
   - 반응식: `3ippm_c + nad_c --> 4mop_c + co2_c + nadh_c`
   - 실제 사용됨: No
   - 유전자: leuB (Isopropylmalate dehydrogenase)

### Acetate Metabolism (1개)

9. **ACKr** (Acetate kinase, reverse)
   - 반응식: `ac_c + atp_c <=> actp_c + adp_c`
   - 실제 사용됨: No
   - 유전자: ackA (Acetate kinase)

### ETC / Electron Transport (3개)

10. **MKRED** (Menaquinone-8 reduction by NADH)
    - 반응식: `h_c + mqn8_c + nadh_c --> mql8_c + nad_c`
    - 실제 사용됨: No
    - 유전자: menaquinone 관련

11. **MQN8RD** (Menaquinone-8 reduction by NADH, alternative)
    - 반응식: `h_c + mqn8_c + nadh_c --> mql8_c + nad_c`
    - 실제 사용됨: No
    - 유전자: menaquinone 관련

12. **MQN8red** (Menaquinone-8 NADH dehydrogenase, lumped)
    - 반응식: `h_c + mqn8_c + nadh_c --> mql8_c + nad_c`
    - 실제 사용됨: No
    - 유전자: menaquinone 관련

### Cofactor Metabolism (추가 확인 필요)
- 해당 반응 없음

### Other (15개)

13. **DHAD** (Dihydroxyacid dehydratase)
    - 반응식: `dhiv_c --> 2kiv_c + h2o_c`
    - 실제 사용됨: No

14. **IPMI** (Isopropylmalate isomerase)
    - 반응식: `2ippm_c <=> 3ippm_c`
    - 실제 사용됨: No

15. **KARI** (Ketol-acid reductoisomerase)
    - 반응식: `alac__S_c + h_c + nadph_c --> dhiv_c + nadp_c`
    - 실제 사용됨: No

16. **PDXK** (Pyridoxal kinase)
    - 반응식: `atp_c + pydxn_c --> adp_c + h_c + pydx5p_c`
    - 실제 사용됨: No

17. **THIK** (Thiamine kinase)
    - 반응식: `atp_c + thm_c --> adp_c + h_c + thmpp_c`
    - 실제 사용됨: No

18-27. 기타 반응들 (CSV 파일 참조)

---

## 요약

### 경로별 분류 (Transport/Exchange 제외)

- **Gluconeogenesis / Glycolysis**: 4개 (PC, PEPCK, PPDK, PEPCK_ATP)
- **TCA Cycle**: 3개 (SUCDi, FRD7, IPMDH)
- **Acetate Metabolism**: 2개 (ACS_ADP, ACKr)
- **ETC / Electron Transport**: 3개 (MKRED, MQN8RD, MQN8red)
- **Other**: 15개

### 우선순위별 분류

- **HIGH**: 3개 (모두 실제 사용됨)
- **MEDIUM**: 24개 (실제 사용 안 됨)
- **LOW**: 0개

---

## 유전자 확인 작업

신규 분리 미생물 genome에서 확인할 주요 유전자:

### HIGH 우선순위 (필수)

1. **acsA** - ACS_ADP
2. **pckA** - PEPCK_ATP
3. **sdhABCD** - SUCDi

### MEDIUM 우선순위 (추가 확인)

4. **pyc/pycA** - PC
5. **ppdK** - PPDK
6. **frdABCD** - FRD7
7. **leuB** - IPMDH
8. **ackA** - ACKr
9. **menA, menB, menC, menD, menE** - Menaquinone 관련
10. **기타** - DHAD, IPMI, KARI, PDXK, THIK 등
