# 신규 모델 탄소 플럭스 경로 분석 요약

## FBA 결과

**성장률**: 모델이 성장 가능 (FBA 최적화 성공)

---

## 주요 발견사항

### 1. TCA Cycle - ✅ 활성화됨

TCA cycle이 완전히 작동하고 있으며, 모든 핵심 반응이 플럭스 **1.600000**으로 활성화되어 있습니다:

1. **CS** (Citrate Synthase): 1.600000
   - `accoa_c + oaa_c --> cit_c`
   - Acetyl-CoA와 Oxaloacetate 결합

2. **ACONT** (Aconitase): 1.600000
   - `cit_c <=> icit_c`
   - Citrate를 Isocitrate로 이성질화

3. **ICDHx** (Isocitrate Dehydrogenase, NAD+): 1.600000
   - `icit_c + nad_c --> akg_c + co2_c + nadh_c`
   - 첫 번째 NADH 생성, CO2 방출

4. **AKGDH** (α-Ketoglutarate Dehydrogenase): 1.600000
   - `akg_c --> succoa_c + co2_c + nadh_c`
   - 두 번째 NADH 생성, CO2 방출

5. **SUCD** (Succinate Dehydrogenase): 1.600000
   - `succ_c --> fum_c`
   - FADH2 생성 (전자 전달)

6. **FUM** (Fumarase): 1.600000
   - `fum_c --> mal__L_c`
   - Fumarate를 Malate로 전환

7. **MDH** (Malate Dehydrogenase): 1.600000
   - `mal__L_c + nad_c --> oaa_c + nadh_c`
   - 세 번째 NADH 생성, Oxaloacetate 재생성

**결론**: TCA cycle이 정상적으로 작동하여 에너지(NADH, FADH2)를 생성하고 있습니다.

---

### 2. Glyoxylate Shunt - ❌ 비활성화됨

Glyoxylate shunt의 두 반응 모두 플럭스가 **0.000000**으로 비활성화되어 있습니다:

1. **ICL** (Isocitrate Lyase): 0.000000
   - `icit_c --> glx_c + succ_c`
   - 비활성

2. **MALS** (Malate Synthase): 0.000000
   - `accoa_c + glx_c --> mal__L_c`
   - 비활성

**결론**: 
- Glyoxylate shunt가 사용되지 않고 있습니다.
- Acetate를 탄소원으로 사용할 때는 일반적으로 Glyoxylate shunt가 활성화되어 탄소를 보존합니다.
- 현재 모델에서는 TCA cycle만 사용하고 있어 탄소가 CO2로 손실됩니다.

---

### 3. Acetate 전환 경로

#### 주 경로: SUCOAACTr (활성화됨)
- **반응**: `ac_c + succoa_c <=> accoa_c + succ_c`
- **플럭스**: 1.600000
- **특징**: 
  - ATP 소모 없음
  - Succinyl-CoA와 CoA 전이 반응
  - TCA cycle의 Succinyl-CoA를 활용

#### 대체 경로: ACS (비활성)
- **반응**: `ac_c + atp_c + coa_c --> accoa_c + amp_c + ppi_c`
- **플럭스**: 0.000000
- **비고**: 사용되지 않음 (에너지 비용이 높음)

**결론**: SUCOAACTr를 통해 Acetate가 Acetyl-CoA로 전환되고 있습니다.

---

## 탄소 플럭스 흐름 경로

### 실제 사용되는 경로:

```
Acetate (ac_c)
  ↓ [SUCOAACTr: 1.600000]
Acetyl-CoA (accoa_c)
  ↓ [CS: 1.600000]
Citrate (cit_c)
  ↓ [ACONT: 1.600000]
Isocitrate (icit_c)
  ↓ [ICDHx: 1.600000]  ← CO2 방출, NADH 생성
α-Ketoglutarate (akg_c)
  ↓ [AKGDH: 1.600000]  ← CO2 방출, NADH 생성
Succinyl-CoA (succoa_c)
  ↓ [SUCOAS]
Succinate (succ_c)
  ↓ [SUCD: 1.600000]  ← FADH2 생성
Fumarate (fum_c)
  ↓ [FUM: 1.600000]
Malate (mal__L_c)
  ↓ [MDH: 1.600000]  ← NADH 생성
Oxaloacetate (oaa_c)
  ↓ [CS: 1.600000]
Citrate (cit_c)  ← TCA cycle 완료
```

### 사용되지 않는 경로:

```
Isocitrate (icit_c)
  ↓ [ICL: 0.000000]  ← 비활성
Glyoxylate (glx_c) + Succinate (succ_c)
  ↓ [MALS: 0.000000]  ← 비활성
Malate (mal__L_c)  ← 탄소 보존 경로 비활성
```

---

## 핵심 요약

1. **TCA Cycle**: ✅ 완전히 활성화되어 정상 작동
2. **Glyoxylate Shunt**: ❌ 비활성화 (탄소 보존 경로 미사용)
3. **Acetate 전환**: SUCOAACTr를 통한 CoA 전이 반응 사용
4. **탄소 손실**: CO2로 인한 탄소 손실 발생 (ICDHx, AKGDH에서)
5. **에너지 생성**: NADH 3개, FADH2 1개 생성 (TCA cycle 1회전당)

---

## 레퍼런스 모델과의 차이

레퍼런스 모델의 경우:
- ACS_ADP 반응 사용 (ATP → ADP + Pi, 더 효율적)
- Glyoxylate shunt 활성화 가능성 (추가 분석 필요)

현재 신규 모델:
- SUCOAACTr 사용 (ATP 소모 없음, Succinyl-CoA 필요)
- Glyoxylate shunt 비활성

---

## 결론 및 권장사항

1. **현재 상태**: 모델이 TCA cycle을 통해 성장 가능
2. **문제점**: 
   - Glyoxylate shunt가 비활성화되어 탄소 효율이 낮을 수 있음
   - CO2로 인한 탄소 손실 발생
3. **개선 방안**:
   - Glyoxylate shunt 활성화 조건 확인
   - 레퍼런스 모델과의 비교 분석
   - ACS_ADP 반응 추가 고려 (에너지 효율 향상)
