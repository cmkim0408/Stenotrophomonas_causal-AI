# 포도당 생장 문제 해결 방안 가이드

## 문제 요약

**질문**: 포도당을 완전 산화시켜서 TCA 반응까지 유도한 뒤, 이걸로 NADH를 만들어 전자전달계로 ATP를 만들 수 있나요?

**답변**: **이론적으로는 가능하지만, 실제로는 초기 부트스트랩 문제로 작동하지 않습니다.**

## 근본 원인

### 1. 순환 의존성 문제
- 포도당 transport → ATP 필요 (GLCabc) 또는 PEP 필요 (GLCpts)
- ATP 생산 → 포도당 필요
- PEP 생산 → ATP 필요
- **결과**: 초기 순환 의존성 발생

### 2. HEX1 단계 문제
- Hexokinase (HEX1): `ATP + Glucose → G6P + ADP`
- 포도당이 세포 내에 있어도 ATP 없이는 G6P 생성 불가능

### 3. 다른 생합성 경로 불완전
- NAD+ 생산 경로 불완전
- CoA 생산 경로 문제 가능성
- 뉴클레오티드 생합성 경로 불완전
- 아미노산 생합성 경로 불완전

## 해결 방안

### 방안 1: 실험적 부트스트랩 (권장)

**원리**: 초기 소량의 ATP/NAD+/CoA를 외부에서 공급하여 순환 의존성 해결

**방법**:
```python
# FBA 실행 시 부트스트랩 추가
demand_atp = cobra.Reaction('DM_atp_c_bootstrap')
demand_atp.lower_bound = -0.01  # 매우 소량만
demand_atp.add_metabolites({atp_c: -1})

demand_nad = cobra.Reaction('DM_nad_c_bootstrap')
demand_nad.lower_bound = -0.01
demand_nad.add_metabolites({nad_c: -1})

demand_coa = cobra.Reaction('DM_coa_c_bootstrap')
demand_coa.lower_bound = -0.001  # 더욱 소량
demand_coa.add_metabolites({coa_c: -1})

model.add_reactions([demand_atp, demand_nad, demand_coa])
```

**실험적 검증**:
- 부유한 배지에서 사전 배양 후 포도당 최소 배지로 전환
- 또는 초기 소량의 ATP/NAD+/CoA를 배지에 추가

### 방안 2: Transport 경로 수정

**방법 A**: ATP 필요 없는 transport 경로 추가
```python
# GLCt: glc__D_e <=> glc__D_c (단순 확산)
glc_transport = cobra.Reaction('GLCt')
glc_transport.lower_bound = -1000  # 가역적
glc_transport.upper_bound = 1000
glc_transport.add_metabolites({
    model.metabolites.get_by_id('glc__D_e'): -1,
    model.metabolites.get_by_id('glc__D_c'): 1
})
model.add_reactions([glc_transport])
```

**방법 B**: PTS 경로 활성화 (PEP 필요하지만 ATP보다 덜 제한적)
- GLCpts: `glc__D_e + PEP → G6P + Pyruvate`
- PEP 부트스트랩으로 시작

### 방안 3: HEX1 가역성 확보

HEX1을 가역적으로 만들어 ATP 없이도 작동 가능하게:
```python
hex1 = model.reactions.get_by_id('HEX1')
hex1.lower_bound = -1000  # 역방향 허용
```

**주의**: 열역학적으로는 역방향이 불리하지만, 모델에서는 작동 가능

### 방안 4: 누락된 반응 추가 (Gap-filling)

**필요한 반응**:
1. **NAD+ 생산 경로**:
   - Nicotinate 생산 경로
   - Salvage pathway

2. **CoA 생산 경로**:
   - Pantothenate 생산
   - CoA 생합성 경로

3. **뉴클레오티드 생산 경로**:
   - Purine 생합성
   - Pyrimidine 생합성
   - NDPK 경로

4. **아미노산 생산 경로**:
   - Serine, Glycine 생합성
   - Aromatic 아미노산 생합성
   - 기타 필수 아미노산

**Gap-filling 도구 사용**:
- Meneco
- GapFind/GapFill
- ModelSEED

## 실용적 해결 방법 (우선순위)

### 즉시 시도 가능한 방법

1. **부트스트랩 추가** (가장 빠름)
   - 소량의 ATP/NAD+/CoA demand reaction 추가
   - 모델 수정 최소화

2. **Transport 경로 추가**
   - GLCt 반응 추가 (ATP 필요 없음)
   - HEX1 가역성 확보

### 장기적 해결 방법

3. **Gap-filling 수행**
   - 누락된 반응 식별
   - 반응 추가 및 검증

4. **실험적 검증**
   - 부트스트랩이 실제로 필요한지 확인
   - 어떤 구성 요소가 실제로 문제인지 확인

## 적용 스크립트

`fix_glucose_growth.py`를 실행하여 다음을 수행:
1. GLCt transport 경로 추가
2. HEX1 가역성 확보
3. 최소 부트스트랩 반응 추가

```bash
python fix_glucose_growth.py
```

## 결론

**포도당 완전 산화 경로는 이론적으로 작동 가능합니다:**
- 포도당 → Glycolysis → Pyruvate
- Pyruvate → PDH → Acetyl-CoA
- Acetyl-CoA → TCA Cycle → NADH
- NADH → ETC → h_p → ATP

**하지만 초기 부트스트랩 문제로 작동하지 않으므로:**
1. 소량의 ATP/NAD+/CoA 부트스트랩 추가
2. 또는 Transport 경로 수정
3. 또는 누락된 반응 추가 (gap-filling)

**가장 실용적인 방법**: 부트스트랩 추가 + Transport 경로 수정
