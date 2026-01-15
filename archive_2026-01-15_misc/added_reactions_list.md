# 추가된 반응 목록

**원본 모델**: BaseModel.xml (2092개 반응)
**최종 모델**: BaseModel_with_BCAA_cofactors_ions_nad_transport.xml (2109개 반응)
**추가된 반응**: 17개

---

## Acetate Metabolism (1개)

### ACS_ADP

**반응식**: `ac_c + atp_c + coa_c <=> accoa_c + adp_c + pi_c`

**bounds**: [-1000.0, 1000.0]

**가역성**: True

---

## Acetate Transport (1개)

### ACtexi

**이름**: Acetate transport (e<->c)

**반응식**: `ac_e <=> ac_c`

**bounds**: [-1000.0, 1000.0]

**가역성**: True

---

## BCAA Synthesis (6개)

### BCAT_LEU

**이름**: BCAT_LEU

**반응식**: `4mop_c + glu__L_c <=> akg_c + leu__L_c`

**bounds**: [-1000.0, 1000.0]

**가역성**: True

---

### BCAT_VAL

**이름**: BCAT_VAL

**반응식**: `2kiv_c + glu__L_c <=> akg_c + val__L_c`

**bounds**: [-1000.0, 1000.0]

**가역성**: True

---

### DHAD

**이름**: DHAD

**반응식**: `dhiv_c --> 2kiv_c + h2o_c`

**bounds**: [0.0, 1000.0]

**가역성**: False

---

### IPMDH

**이름**: IPMDH

**반응식**: `3ippm_c + nad_c --> 4mop_c + co2_c + nadh_c`

**bounds**: [0.0, 1000.0]

**가역성**: False

---

### IPMI

**이름**: IPMI

**반응식**: `2ippm_c <=> 3ippm_c`

**bounds**: [-1000.0, 1000.0]

**가역성**: True

---

### KARI

**이름**: KARI

**반응식**: `alac__S_c + h_c + nadph_c --> dhiv_c + nadp_c`

**bounds**: [0.0, 1000.0]

**가역성**: False

---

## Exchange (1개)

### EX_hco3_e

**이름**: Bicarbonate exchange

**반응식**: `hco3_e <=> `

**bounds**: [-1000.0, 1000.0]

**가역성**: True

---

## Gluconeogenesis (1개)

### PEPCK_ATP

**이름**: PEP carboxykinase (ATP)

**반응식**: `atp_c + oaa_c --> adp_c + co2_c + pep_c`

**bounds**: [0.0, 1000.0]

**가역성**: False

---

## Ion Transport (1개)

### T_nac_e_to_nac_c

**이름**: Nicotinate transport (e<->c)

**반응식**: `nac_e <=> nac_c`

**bounds**: [-1000.0, 1000.0]

**가역성**: True

---

## Ion Transport (Chloride) (1개)

### T_cl_e_to_cl_c

**이름**: T_cl_e_to_cl_c

**반응식**: `cl_e <=> cl_c`

**bounds**: [-1000.0, 1000.0]

**가역성**: True

---

## Ion Transport (Cobalt) (1개)

### T_cobalt2_e_to_cobalt2_c

**이름**: T_cobalt2_e_to_cobalt2_c

**반응식**: `cobalt2_e <=> cobalt2_c`

**bounds**: [-1000.0, 1000.0]

**가역성**: True

---

## Ion Transport (Copper) (1개)

### T_cu2_e_to_cu2_c

**이름**: T_cu2_e_to_cu2_c

**반응식**: `cu2_e <=> cu2_c`

**bounds**: [-1000.0, 1000.0]

**가역성**: True

---

## NAD/NADP Exchange (2개)

### EX_nac_e

**이름**: nac_e exchange

**반응식**: `nac_e <=> `

**bounds**: [-1000.0, 1000.0]

**가역성**: True

---

### EX_ncam_e

**이름**: ncam_e exchange

**반응식**: `ncam_e <=> `

**bounds**: [-1000.0, 1000.0]

**가역성**: True

---

## TCA Cycle (1개)

### SUCDi

**이름**: Succinate dehydrogenase

**반응식**: `q8_c + succ_c --> fum_c + q8h2_c`

**bounds**: [0.0, 1000.0]

**가역성**: False

---

