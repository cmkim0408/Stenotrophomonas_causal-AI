# NAD/NADP 생산 경로 분석

## 현재 상황

### 막힌 Metabolite (Step 4-2 재실행 결과)
- nad_c (NAD): 계수 0.001831
- nadp_c (NADP): 계수 0.000447

### 확인된 반응들

**이미 모델에 있는 NAD/NADP 관련 반응**:
- ✅ NADK: `atp_c + nad_c --> adp_c + h_c + nadp_c` (NAD → NADP)
- ✅ NADS1: `atp_c + dnad_c + nh4_c --> amp_c + h_c + nad_c + ppi_c` (dnad_c → nad_c)
- ✅ NADS2: `atp_c + dnad_c + gln__L_c + h2o_c --> amp_c + glu__L_c + h_c + nad_c + ppi_c` (dnad_c → nad_c)
- ✅ THD2pp: `2.0 h_p + nadh_c + nadp_c --> 2.0 h_c + nad_c + nadph_c` (Transhydrogenase)
- ✅ NADTRHD: `nad_c + nadph_c --> nadh_c + nadp_c` (Transhydrogenase)
- ✅ NAPRT: `h_c + nac_c + prpp_c <=> nicrnt_c + ppi_c` (Nicotinate → Nicotinate ribonucleotide)
- ✅ NMNAT: `atp_c + h_c + nmn_c --> nad_c + ppi_c` (NMN → NAD)
- ✅ NNATr: `atp_c + h_c + nicrnt_c <=> dnad_c + ppi_c` (nicrnt_c → dnad_c)

**추가된 Exchange 반응**:
- ✅ EX_nac_e: Nicotinate exchange
- ✅ EX_ncam_e: Nicotinamide exchange

---

## NAD 생산 경로

### 경로 1: Nicotinate 경로 (De novo)
1. **EX_nac_e** → nac_c (Nicotinate) ✅
2. **NAPRT**: nac_c + prpp_c → nicrnt_c + ppi_c ✅
3. **NNATr**: nicrnt_c + atp_c + h_c → dnad_c + ppi_c ✅
4. **NADS1/NADS2**: dnad_c → nad_c ✅
5. **NADK**: nad_c → nadp_c ✅

### 경로 2: Nicotinamide 경로 (Salvage)
1. **EX_ncam_e** → ncam_c (Nicotinamide) ✅
2. **NAMPT**: ncam_c + prpp_c → nmn_c + ppi_c (❌ 반응 없음)
3. **NMNAT**: nmn_c + atp_c + h_c → nad_c + ppi_c ✅

---

## 문제점 분석

### 경로 1 문제
- 모든 반응이 있는데도 생성되지 않음
- 가능한 원인:
  1. **PRPP 생성 문제**: PRPP가 생성되지 않으면 NAPRT가 작동하지 않음
  2. **nicrnt_c 생성 문제**: NAPRT가 작동하지 않으면 nicrnt_c가 생성되지 않음
  3. **경로 연결 문제**: 반응들이 서로 연결되지 않음

### 경로 2 문제
- **NAMPT 반응이 없음**: ncam_c → nmn_c 경로가 막혀 있음
- 레퍼런스 모델에도 NAMPT가 없는 것으로 확인됨

---

## 다음 단계

1. **PRPP 생산 확인**: PRPP가 생성되는지 확인
2. **nicrnt_c 생산 확인**: NAPRT가 작동하는지 확인
3. **전체 경로 연결 확인**: 경로가 완전히 연결되어 있는지 확인
4. **NAMPT 반응 추가 검토**: 레퍼런스 모델에 있는지 확인
