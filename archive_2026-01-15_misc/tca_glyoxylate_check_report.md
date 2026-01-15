# TCA Cycle 및 Glyoxylate Shunt 점검 보고서

## 📊 점검 결과 요약

### ✅ 전체 상태: 정상

- **TCA Cycle**: 9/9 반응 존재
- **Glyoxylate Shunt**: 2/2 반응 존재  
- **경로 연결성**: OK
- **블록된 반응**: 0개

---

## 🔄 TCA Cycle 반응 상세

### 1. CS (Citrate Synthase) - ✅ 존재
- **반응식**: `accoa_c + h2o_c + oaa_c --> cit_c + coa_c + h_c`
- **가역성**: No
- **유전자**: FAGFNPBA_03340 (1개)
- **역할**: TCA cycle 시작점, Acetyl-CoA와 Oxaloacetate 결합

### 2. ACONT (Aconitase) - ✅ 존재
- **반응식**: `cit_c <=> icit_c`
- **가역성**: Yes
- **유전자**: FAGFNPBA_01989, FAGFNPBA_03113, FAGFNPBA_01993 (3개)
- **역할**: Citrate를 Isocitrate로 이성질화

### 3. ICDHx (Isocitrate Dehydrogenase, NAD+) - ✅ 존재
- **반응식**: `icit_c + nad_c <=> akg_c + co2_c + nadh_c`
- **가역성**: Yes
- **유전자**: FAGFNPBA_00882 (1개)
- **역할**: Isocitrate 탈수소화, 첫 번째 NADH 생성

### 4. ICDHyr (Isocitrate Dehydrogenase, NADP+) - ✅ 존재
- **반응식**: `icit_c + nadp_c <=> akg_c + co2_c + nadph_c`
- **가역성**: Yes
- **유전자**: FAGFNPBA_00882, FAGFNPBA_03758 (2개)
- **역할**: NADP+ 의존성 대체 경로

### 5. AKGDH (α-Ketoglutarate Dehydrogenase) - ✅ 존재
- **반응식**: `akg_c + coa_c + nad_c --> co2_c + nadh_c + succoa_c`
- **가역성**: No
- **유전자**: FAGFNPBA_02721, FAGFNPBA_02719, FAGFNPBA_02720, FAGFNPBA_03607 (4개)
- **역할**: 복합체 반응, 두 번째 NADH 생성

### 6. SUCOAS (Succinyl-CoA Synthetase) - ✅ 존재
- **반응식**: `atp_c + coa_c + succ_c <=> adp_c + pi_c + succoa_c`
- **가역성**: Yes
- **유전자**: FAGFNPBA_03254, FAGFNPBA_03255, FAGFNPBA_00400 (3개)
- **역할**: Substrate-level phosphorylation, ATP 생성

### 7. SUCD (Succinate Dehydrogenase) - ✅ 추가 완료
- **반응식**: `succ_c + fad_c --> fum_c + fadh2_c`
- **가역성**: No
- **유전자**: Smlt1796 (1개) - succinate dehydrogenase cytochrome b-556 subunit
- **역할**: Succinate를 Fumarate로 변환, FADH2 생성, 전자전달계 복합체 II
- **상태**: ✅ 유전자 기반으로 반응 추가 완료

### 8. FUM (Fumarase) - ✅ 존재
- **반응식**: `fum_c + h2o_c <=> mal__L_c`
- **가역성**: Yes
- **유전자**: FAGFNPBA_02715, FAGFNPBA_02744 (2개)
- **역할**: Fumarate를 Malate로 가수화

### 9. MDH (Malate Dehydrogenase) - ✅ 존재
- **반응식**: `mal__L_c + nad_c <=> h_c + nadh_c + oaa_c`
- **가역성**: Yes
- **유전자**: FAGFNPBA_00845 (1개)
- **역할**: 마지막 탈수소화, OAA 재생성, 세 번째 NADH 생성

---

## 🔁 Glyoxylate Shunt 반응 상세

### 1. ICL (Isocitrate Lyase) - ✅ 존재
- **반응식**: `icit_c --> glx_c + succ_c`
- **가역성**: No
- **유전자**: FAGFNPBA_00226 (1개)
- **역할**: Isocitrate를 분해하여 Glyoxylate와 Succinate 생성

### 2. MALS (Malate Synthase) - ✅ 존재
- **반응식**: `accoa_c + glx_c + h2o_c --> coa_c + h_c + mal__L_c`
- **가역성**: No
- **유전자**: FAGFNPBA_00225 (1개)
- **역할**: Glyoxylate와 Acetyl-CoA로부터 Malate 합성

---

## 🔗 경로 연결성 확인

### TCA Cycle 경로 연결성
모든 단계가 정상적으로 연결되어 있습니다:

1. ✅ `accoa_c` → `cit_c` (1개 반응)
2. ✅ `cit_c` → `icit_c` (1개 반응)
3. ✅ `icit_c` → `akg_c` (2개 반응: ICDHx, ICDHyr)
4. ✅ `akg_c` → `succoa_c` (3개 반응)
5. ✅ `succoa_c` → `succ_c` (5개 반응)
6. ✅ `succ_c` → `fum_c` (2개 반응, SUCFUMt 포함)
7. ✅ `fum_c` → `mal__L_c` (1개 반응)
8. ✅ `mal__L_c` → `oaa_c` (4개 반응)

### Glyoxylate Shunt 경로 연결성
모든 단계가 정상적으로 연결되어 있습니다:

1. ✅ `icit_c` → `glx_c` (1개 반응: ICL)
2. ✅ `glx_c` → `mal__L_c` (1개 반응: MALS)

---

## 📦 주요 대사물질 확인

모든 필수 대사물질이 정상적으로 존재합니다:

- ✅ Acetyl-CoA (accoa_c): 63개 반응 참여
- ✅ Citrate (cit_c): 10개 반응 참여
- ✅ Isocitrate (icit_c): 6개 반응 참여
- ✅ α-Ketoglutarate (akg_c): 40개 반응 참여
- ✅ Succinyl-CoA (succoa_c): 11개 반응 참여
- ✅ Succinate (succ_c): 16개 반응 참여
- ✅ Fumarate (fum_c): 10개 반응 참여
- ✅ Malate (mal__L_c): 10개 반응 참여
- ✅ Oxaloacetate (oaa_c): 12개 반응 참여
- ✅ Glyoxylate (glx_c): 4개 반응 참여

---

## ✅ 최근 업데이트

### SUCD 반응 추가 완료
- **문제**: 초기 모델에 SUCD 반응이 없었습니다.
- **해결**: Smlt1796 유전자(sdhC, succinate dehydrogenase cytochrome b-556 subunit)를 확인하고 SUCD 반응을 추가했습니다.
- **반응식**: `succ_c + fad_c --> fum_c + fadh2_c`
- **유전자 연결**: Smlt1796
- **결과**: TCA cycle이 완전히 구성되었습니다 (9/9 반응).

---

## ✅ 최종 결론

### TCA Cycle
- **상태**: ✅ 정상 (SUCD 반응 추가 완료)
- **반응 존재율**: 9/9 (100%)
- **경로 연결성**: ✅ 완전
- **블록된 반응**: 0개
- **최근 업데이트**: SUCD 반응 추가 (Smlt1796 유전자 연결)

### Glyoxylate Shunt
- **상태**: ✅ 정상
- **반응 존재율**: 2/2 (100%)
- **경로 연결성**: ✅ 완전
- **블록된 반응**: 0개

### 전체 평가
**TCA Cycle과 Glyoxylate Shunt가 완전히 구성되어 있으며, 모든 경로가 정상적으로 연결되어 있습니다.**

Acetate 기반 성장에 필요한 모든 대사 경로가 준비되어 있습니다:
- ✅ Acetate → Acetyl-CoA
- ✅ Acetyl-CoA → TCA cycle
- ✅ Acetyl-CoA → Glyoxylate shunt → Malate
- ✅ Malate → Oxaloacetate 재생성

---

**점검일**: 2024년
**모델**: BaseModel.xml
