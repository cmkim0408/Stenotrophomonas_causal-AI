# 포도당 생장 문제 해결 방안 요약

## 핵심 문제

포도당 완전 산화 경로(포도당 → TCA → NADH → ETC → ATP)는 **이론적으로 존재하고 작동 가능**하지만, **초기 부트스트랩 문제**로 작동하지 않습니다.

## 즉시 적용 가능한 해결 방법

### 방법 1: 부트스트랩 추가 (가장 빠름)

```python
import cobra

model = cobra.io.read_sbml_model("BaseModel.xml")

# 소량의 ATP/NAD+/CoA 부트스트랩 추가
atp_c = model.metabolites.get_by_id('atp_c')
dm_atp = cobra.Reaction('DM_atp_c_bootstrap')
dm_atp.lower_bound = -0.01
dm_atp.add_metabolites({atp_c: -1})
model.add_reactions([dm_atp])

nad_c = model.metabolites.get_by_id('nad_c')
dm_nad = cobra.Reaction('DM_nad_c_bootstrap')
dm_nad.lower_bound = -0.01
dm_nad.add_metabolites({nad_c: -1})
model.add_reactions([dm_nad])

coa_c = model.metabolites.get_by_id('coa_c')
dm_coa = cobra.Reaction('DM_coa_c_bootstrap')
dm_coa.lower_bound = -0.001
dm_coa.add_metabolites({coa_c: -1})
model.add_reactions([dm_coa])
```

### 방법 2: Transport 경로 수정

```python
# ATP 필요 없는 transport 추가
glc__D_e = model.metabolites.get_by_id('glc__D_e')
glc__D_c = model.metabolites.get_by_id('glc__D_c')
glc_transport = cobra.Reaction('GLCt')
glc_transport.lower_bound = -1000
glc_transport.upper_bound = 1000
glc_transport.add_metabolites({glc__D_e: -1, glc__D_c: 1})
model.add_reactions([glc_transport])

# HEX1 가역성 확보
hex1 = model.reactions.get_by_id('HEX1')
hex1.lower_bound = -1000
```

## 전체 해결 스크립트

`fix_glucose_growth.py` 실행:
- Transport 경로 수정
- HEX1 가역성 확보  
- 최소 부트스트랩 추가

**주의**: 현재 모델에서는 부트스트랩만으로는 부족하고, 추가 gap-filling이 필요할 수 있습니다.

## 다음 단계

1. **실험적 검증**: 부트스트랩이 실제로 필요한지 확인
2. **Gap-filling**: 누락된 반응 추가
3. **모델 검증**: 포도당으로 실제 생장 가능한지 확인
