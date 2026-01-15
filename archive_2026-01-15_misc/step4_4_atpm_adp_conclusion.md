# Step 4-4: ATPM-ADP 원인 확인 결과

## 1. ATPM 반응식 확인

**ATPM 반응식:**
```
atp_c + h2o_c --> adp_c + h_c + pi_c
```

**핵심 발견:**
- ATPM은 ATP를 소비하고 **ADP를 생성**합니다 (ADP 계수: +1.0)
- ATPM이 ADP 재생산 경로 중 하나입니다

## 2. ATPM=0일 때 ADP/ATP 균형

**FBA 결과:**
- 상태: optimal
- **성장률: 3.293964** ✅ (성장 가능!)

**ADP 생성 경로 (ATPM=0일 때):**
1. ACS_ADP: 플럭스 110.51, ADP 생성 +1.0
2. PEPCK_ATP: 플럭스 19.26, ADP 생성 +1.0
3. PGK: 플럭스 11.55, ADP 생성 +1.0
4. ADK1: 플럭스 -9.47, ADP 생성 +2.0 (AMP + ATP → 2ADP)
5. GLNS, ASPK, CDPMEK, CYTK1, NDPK3 등: 여러 반응에서 ADP 생성
6. Growth: 플럭스 3.29, ADP 생성 +53.95

**ADP 소비 경로:**
1. ATPS4rpp: 플럭스 345.75, ADP 소비 -1.0 (ATP 생성에 사용)
2. RNDR1b: 플럭스 0.09, ADP 소비 -1.0

**ATP 생성:**
- ATPS4rpp: 플럭스 345.75, ATP 생성 +1.0

## 3. ATP 생성 경로에 ADP 필요 여부

**ATP 생성 반응 중 ADP를 소비하는 반응:**
1. **ATPS4rpp**: `adp_c + 4.0 h_p + pi_c <=> atp_c + h2o_c + 3.0 h_c`
   - ATP 생성 계수: +1.00
   - ADP 소비 계수: -1.00
2. PPAKr: `adp_c + ppap_c <=> atp_c + ppa_c`
3. PYK: `adp_c + h_c + pep_c --> atp_c + pyr_c`

**결론:**
- 주요 ATP 생성 경로(ATPS4rpp)는 ADP를 필요로 합니다
- 하지만 ADP는 ATPM 외에도 **다양한 경로에서 충분히 생성**됩니다

## 4. ATPM=0 vs ATPM=5 비교

### ATPM=0
- 성장률: **3.293964**
- ATPM 플럭스: 0.0
- ATPS4rpp 플럭스: 345.75

### ATPM=5
- 성장률: **3.247537** (약간 낮음)
- ATPM 플럭스: 5.0
- ATPS4rpp 플럭스: 345.88

**관찰:**
- ATPM=0일 때도 성장이 가능합니다
- ATPM이 높을수록 ATP 소비가 증가하므로 성장률이 약간 낮아집니다
- ATPM은 필수 경로가 아닙니다 (다른 경로에서 ADP 충분히 생성)

## 5. ADP 생성 경로 (총 265개)

주요 ADP 생성 경로:
- **ATPM**: ATP → ADP + Pi
- **ACS_ADP**: Acetate + ATP → Acetyl-CoA + ADP + Pi
- **PEPCK_ATP**: OAA + ATP → PEP + CO2 + ADP
- **PGK**: 1,3-BPG + ADP → 3PG + ATP (역방향)
- **ADK1**: AMP + ATP → 2ADP (Adenylate Kinase)
- **Growth**: Biomass 합성 과정에서 ADP 생성
- 기타 260개 이상의 반응에서 ADP 생성

## 결론

### "ADP Deadlock" 가설 검증 결과:

**원래 가설:**
- ATPM=0일 때 ADP가 재생산되지 않아 전체 에너지 대사가 멈춤
- ATP 생성 경로(ATPS4rpp)가 ADP를 필요로 하므로, ADP 부족 시 ATP 생성 불가

**실제 관찰:**
1. ✅ **ATPM=0일 때도 성장 가능** (성장률: 3.29)
2. ✅ **ADP는 다양한 경로에서 생성됨**:
   - ACS_ADP, PEPCK_ATP, PGK, ADK1, Growth 등
   - 총 265개 이상의 반응에서 ADP 생성
3. ✅ **ATP 생성 경로(ATPS4rpp)는 ADP를 필요로 하나**, 충분히 공급됨
4. ✅ **ATPM은 ADP 재생산 경로 중 하나일 뿐**, 필수 경로가 아님

### 최종 결론:

**"ADP Deadlock" 가설은 현재 모델에서는 성립하지 않습니다.**

ATPM=0일 때 성장이 불가능했던 원인은:
- ❌ ADP 부족 (X)
- ✅ **Biosynthetic Gaps** (O)
  - BCAA 합성 경로 누락
  - NAD/NADP 합성 경로 누락 (transport 반응 포함)
  - 이온 수송 경로 누락

**현재 상태:**
- 모든 biosynthetic gaps가 해결됨
- ATPM=0일 때도 성장 가능 (성장률: 3.29)
- ATPM은 선택적 경로 (있어도 되고 없어도 됨)
