# 신규 모델에서 아세트산(acetate) 전환 경로 분석

## 요약

현재 신규 모델에서는 아세트산이 다음과 같은 경로로 전환됩니다:

### 1. 직접 전환 경로 (2개)

#### ACS (Acetyl-CoA synthetase)
- **반응식**: `ac_c + atp_c + coa_c --> accoa_c + amp_c + ppi_c`
- **경로**: Acetate Metabolism
- **비고**: ATP를 AMP + PPi로 분해 (ADP 형성 아님)
- **레퍼런스 모델과의 차이**: 레퍼런스 모델에는 ACS_ADP 반응이 있음 (ADP 형성)

#### SUCOAACTr (Succinyl-CoA:acetate CoA transferase)
- **반응식**: `ac_c + succoa_c <=> accoa_c + succ_c`
- **경로**: Acetate Metabolism
- **비고**: Succinyl-CoA와 CoA를 교환하여 Acetyl-CoA 형성
- **특징**: ATP 소모 없음, Succinyl-CoA가 필요

---

## 레퍼런스 모델과의 차이점

### 레퍼런스 모델에 있지만 신규 모델에 없는 반응

#### ACS_ADP (Acetate-CoA ligase, ADP-forming)
- **반응식**: `ac_c + atp_c + coa_c <=> accoa_c + adp_c + pi_c`
- **차이점**: 
  - ATP → ADP + Pi (에너지 효율이 더 좋음)
  - Reversible (양방향 반응)
- **레퍼런스 모델 FBA에서 실제 사용됨**: Yes (플럭스: 0.943448)

---

## 현재 상황에서의 Acetate 전환 경로

### 경로 1: ACS (주 경로)
```
ac_c + atp_c + coa_c --> accoa_c + amp_c + ppi_c
```
- **에너지 소모**: ATP 1개 → AMP + PPi
- **에너지 비용**: 높음 (AMP는 ATP로 재생성하는데 2 ATP 필요)

### 경로 2: SUCOAACTr (대체 경로)
```
ac_c + succoa_c <=> accoa_c + succ_c
```
- **에너지 소모**: 없음 (CoA 전이 반응)
- **필수 조건**: Succinyl-CoA 필요
- **특징**: TCA cycle이 작동해야 가능

---

## 문제점 및 의의

1. **ACS vs ACS_ADP 차이**
   - 신규 모델: ACS (ATP → AMP + PPi)
   - 레퍼런스 모델: ACS_ADP (ATP → ADP + Pi)
   - ACS가 에너지 효율이 낮음 (AMP → ATP 재생성에 추가 에너지 필요)

2. **레퍼런스 모델에서 ACS_ADP가 실제 사용됨**
   - 레퍼런스 모델의 FBA 결과에서 ACS_ADP가 실제로 사용됨 (플럭스 0.943448)
   - 신규 모델에는 이 반응이 없어서 효율성이 낮을 수 있음

3. **대체 경로 (SUCOAACTr)의 제약**
   - Succinyl-CoA가 필요하므로 TCA cycle이 선행되어야 함
   - 초기 단계에서는 활용하기 어려울 수 있음

---

## 결론

현재 신규 모델에서는:
- ✅ ACS 반응을 통해 Acetate → Acetyl-CoA 전환이 가능
- ✅ SUCOAACTr 반응을 통한 대체 경로도 존재
- ⚠️ 하지만 레퍼런스 모델의 ACS_ADP 반응이 없어서 에너지 효율이 낮을 수 있음
- 🔍 ACS_ADP 반응 추가 시 에너지 효율성 향상 가능
