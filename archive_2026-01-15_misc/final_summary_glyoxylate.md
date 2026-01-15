# Glyoxylate Shunt vs TCA Cycle 최종 정리

## 레퍼런스 모델 결과 (Stenotrophomonas/fba_flux_gradient_acid.csv)

| ATPM | ICL | MALS | ICDHx | 경로 |
|------|-----|------|-------|------|
| 0 | 0.367 | 0.367 | 0.0 | **Glyoxylate shunt만** 사용 |
| 5 | 0.074 | 0.074 | 0.793 | 혼합 (Glyoxylate + TCA) |
| 10 | 0.0 | 0.0 | 1.0 | **TCA cycle만** 사용 |

**핵심**: ATPM이 낮을수록 → Glyoxylate shunt 활성 (탄소 보존)  
**핵심**: ATPM이 높을수록 → TCA cycle 활성 (에너지 생성)

---

## 신규 모델 결과

| ATPM | Growth | ICL | MALS | ICDHx | CS | 경로 |
|------|--------|-----|------|-------|----|------|
| 0 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | **성장 불가** |
| 5 | ~0 | 0.062 | 0.062 | 0.743 | 0.805 | **혼합** (Glyoxylate + TCA) |
| 10 | 0.000 | 0.124 | 0.124 | 1.486 | 1.610 | **혼합** (Glyoxylate + TCA) |

---

## 핵심 차이점

### 1. 레퍼런스 모델
- **ATPM=0**: Glyoxylate shunt만 사용 (ICL=0.367, ICDHx=0.0)
- **ATPM=10**: TCA cycle만 사용 (ICL=0.0, ICDHx=1.0)
- **명확한 전환**: ATPM 값에 따라 한 경로만 선택

### 2. 신규 모델
- **ATPM=0**: 성장 불가 (모든 플럭스 0) ❌
- **ATPM=5, 10**: **Glyoxylate shunt + TCA cycle 동시 사용**
  - ICL/MALS도 활성 (0.062~0.124)
  - ICDHx도 활성 (0.743~1.486)
- **혼합 경로**: 두 경로를 동시에 사용

---

## 사용자 지적사항

> "ATPM에서 ATP 소모가 크지 않다면 glyoxylate shunt가 메인이 되어야 할거 같은데?"

**맞습니다!**

- 레퍼런스 모델: ATPM=0일 때 Glyoxylate shunt가 메인 경로 (ICL=0.367)
- 신규 모델: ATPM=5, 10일 때 Glyoxylate shunt도 활성화되어 있음 (ICL=0.062~0.124)
  - 하지만 TCA cycle도 동시에 사용 (ICDHx=0.743~1.486)
  - 따라서 Glyoxylate shunt가 "메인"이라기보다는 "혼합 사용"

---

## 문제점

### 1. ATPM=0에서 성장 불가 ❌
- 레퍼런스 모델은 ATPM=0에서도 ICL/MALS가 활성화되어 탄소 경로가 작동
- 신규 모델은 ATPM=0에서 모든 플럭스가 0 (성장 불가)
- **원인 불명**: 모델 제약 조건이나 누락된 반응 가능성

### 2. 레퍼런스와 다른 패턴
- 레퍼런스: ATPM 값에 따라 **한 경로만** 선택 (명확한 전환)
- 신규: ATPM 5, 10에서 **두 경로 동시 사용** (혼합 경로)

---

## 결론

1. **신규 모델에서 Glyoxylate shunt는 활성화되어 있습니다**
   - ATPM=5, 10에서 ICL=0.062~0.124, MALS=0.062~0.124

2. **하지만 레퍼런스 모델과는 다른 패턴입니다**
   - 레퍼런스: ATPM 높음 → Glyoxylate shunt 비활성, TCA만
   - 신규: ATPM 5, 10 → Glyoxylate shunt + TCA 혼합

3. **문제점**
   - ATPM=0에서 성장 불가 (레퍼런스 모델과 다름)
   - 추가 분석 필요
