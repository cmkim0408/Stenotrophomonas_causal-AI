# 성장 실패 원인 분석 리포트

## 확인된 사실
- ✓ Exchange reaction (EX_ac_e) 추가 완료
- ✓ Transport reaction (ACt) 추가 완료  
- ✓ Acetate 경로 (ACS, TCA, Glyoxylate) 반응 존재
- ✓ Biomass reaction (Growth) 존재
- ✓ 에너지 루프 없음 (정상)

## 발견된 문제들

### 1. ATP 생산 불가
- **원인**: Acetate만으로는 ATP를 생성할 수 없음
- **영향**: 모든 에너지 의존 반응 (ACS 포함) 차단
- **증상**: ATP 생산 테스트 실패 (objective = 0)

### 2. CoA 생산 불가  
- **원인**: Acetyl-CoA가 없어서 CoA를 재생할 수 없음
- **영향**: ACS 반응 작동 불가
- **증상**: CoA 생산 테스트 실패 (objective = 0)

### 3. OAA Bootstrap 문제
- **원인**: 초기 OAA 생성 경로 없음
- **순환 의존성**:
  - CS는 OAA 필요
  - OAA는 MDH (Malate → OAA)로 생성 가능
  - Malate는 MALS (Glyoxylate shunt)로 생성
  - Glyoxylate는 ICL (Isocitrate → Glyoxylate)로 생성
  - Isocitrate는 ACONT (Citrate → Isocitrate)로 생성
  - Citrate는 CS (Acetyl-CoA + OAA → Citrate)로 생성
  - → **순환! 초기 OAA 없이는 시작 불가**

### 4. PEP 생산 불가
- **원인**: ATP가 없어서 PEP 생산 경로 작동 불가
- **영향**: PPC (PEP → OAA) 경로 사용 불가

## 근본 원인: Bootstrap Problem

모델이 '닭과 달걀' 문제에 직면:
- ATP가 필요하지만 ATP 생산에 Acetyl-CoA 필요
- Acetyl-CoA 생산에 ATP 필요 (ACS 반응)
- 순환 의존성으로 인해 경로 시작 불가

## 해결 방안

### 방안 1: 초기 Bootstrap 대사물질 제공
- 소량의 ATP 제공 (예: 0.1 mmol/gDW/h)
- 또는 소량의 CoA 제공
- 또는 소량의 OAA 제공
- → 이렇게 하면 경로가 시작될 수 있음

### 방안 2: 대체 경로 추가
- ATP 생산을 위한 다른 경로
- OAA 초기화 경로 (예: Pyruvate carboxylase 활성화)
- CoA 초기화 경로

### 방안 3: 반응 방향성 조정
- 일부 반응의 가역성 확인
- 경계 조건 조정

## 권장 조치사항
1. Bootstrap 테스트 수행 - 소량의 ATP/CoA/OAA 제공
2. Gap analysis 상세 수행 - 누락된 반응 식별
3. 모델 수정 - 필요한 반응 추가, 반응 방향성 수정

