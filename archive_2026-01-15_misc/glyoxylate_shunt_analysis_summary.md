# Glyoxylate Shunt 분석 요약

## 레퍼런스 모델 결과

레퍼런스 모델 (`Stenotrophomonas/fba_flux_gradient_acid.csv`)에서:

| ATPM | ICL | MALS | ICDHx | 경로 |
|------|-----|------|-------|------|
| 0 | 0.367 | 0.367 | 0.0 | **Glyoxylate shunt만** |
| 5 | 0.074 | 0.074 | 0.793 | 혼합 |
| 10 | 0.0 | 0.0 | 1.0 | **TCA cycle만** |

**결론**: 
- ATPM이 낮을수록 → Glyoxylate shunt 활성 (탄소 보존)
- ATPM이 높을수록 → TCA cycle 활성 (에너지 생성)

---

## 신규 모델 결과

신규 모델에서:

| ATPM | ICL | MALS | ICDHx | CS | 성장률 | 경로 |
|------|-----|------|-------|----|--------|------|
| 0 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 성장 불가 |
| 5 | 0.062 | 0.062 | 0.743 | 0.805 | ~0 | **혼합** (Glyoxylate + TCA) |
| 10 | 0.124 | 0.124 | 1.486 | 1.610 | 0.000 | **혼합** (Glyoxylate + TCA) |

**결론**:
- ATPM=0: 성장 불가 (모든 플럭스 0)
- ATPM=5, 10: **Glyoxylate shunt와 TCA cycle을 동시에 사용** (혼합 경로)
- 레퍼런스 모델과 달리 ATPM이 높아도 Glyoxylate shunt가 비활성화되지 않음

---

## 차이점 분석

### 레퍼런스 모델
- **ATPM 낮음** → Glyoxylate shunt만 사용 (탄소 효율 최대화)
- **ATPM 높음** → TCA cycle만 사용 (에너지 효율 최대화)
- **명확한 전환**: ATPM 값에 따라 한 경로만 선택

### 신규 모델
- **ATPM=0** → 성장 불가
- **ATPM=5, 10** → Glyoxylate shunt + TCA cycle **동시 사용**
- **혼합 경로**: 두 경로를 동시에 사용하여 탄소와 에너지 균형 유지

---

## 의미

1. **신규 모델이 더 유연함**: 
   - ATPM 값에 관계없이 Glyoxylate shunt와 TCA cycle을 동시에 사용 가능
   - 레퍼런스 모델처럼 한 경로만 선택하는 것이 아님

2. **문제점**:
   - ATPM=0일 때 성장 불가 (레퍼런스는 성장 가능)
   - 이는 모델에 다른 제약이나 누락된 반응이 있을 수 있음을 시사

3. **긍정적 측면**:
   - ATPM=5, 10에서 Glyoxylate shunt가 활성화되어 탄소 보존 경로 사용
   - TCA cycle과 함께 사용하여 에너지 생성도 유지

---

## 사용자 지적사항

> "ATPM에서 ATP 소모가 크지 않다면 glyoxylate shunt가 메인이 되어야 할거 같은데?"

**맞습니다!** 

- 레퍼런스 모델에서는 ATPM이 낮을수록 Glyoxylate shunt가 메인 경로입니다.
- 신규 모델에서도 ATPM=5, 10에서 Glyoxylate shunt가 활성화되어 있습니다 (플럭스: 0.062~0.124).
- 하지만 TCA cycle도 동시에 활성화되어 있어서 (ICDHx 플럭스: 0.743~1.486), Glyoxylate shunt가 "메인"이라고 하기보다는 "혼합 경로"로 보입니다.

---

## 결론

1. **신규 모델에서 Glyoxylate shunt는 활성화되어 있습니다.**
   - ATPM=5, 10에서 ICL과 MALS 플럭스가 0이 아님
   
2. **하지만 레퍼런스 모델과는 다른 패턴입니다.**
   - 레퍼런스: ATPM 높음 → TCA만, ATPM 낮음 → Glyoxylate만
   - 신규: ATPM 5, 10 → Glyoxylate + TCA 동시 사용

3. **추가 확인 필요**:
   - ATPM=0일 때 성장 불가 원인
   - 왜 레퍼런스 모델처럼 한 경로만 선택하지 않는지
   - 모델 제약 조건이나 반응 누락 여부
