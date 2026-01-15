#!/usr/bin/env python
"""
종합 진단 보고서 생성
1, 2, 3단계 모든 결과를 요약
"""

import cobra
import pandas as pd

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    return model

def find_biomass_reaction(model):
    for name in ['Growth', 'BIOMASS']:
        try:
            return model.reactions.get_by_id(name)
        except KeyError:
            continue
    return None

def create_diagnosis_report():
    """종합 진단 보고서 생성"""
    
    print("="*70)
    print("모델 종합 진단 보고서")
    print("="*70)
    
    print("\n[진단 완료된 항목]")
    print("\n1. Biomass 반응 검증")
    print("  - Biomass 반응: Growth (존재 확인)")
    print("  - 총 구성 요소: 57개")
    print("  - 주요 구성 요소:")
    print("    * ATP: -54.12 (대량 필요, 생산 필요)")
    print("    * 아미노산: 20종 (생산 필요)")
    print("    * 기타 뉴클레오티드: GTP, UTP, CTP, dNTP (생산 필요)")
    print("    * 보조인자: CoA, NAD+, NADP+, FAD 등")
    
    print("\n2. 기본 대사 경로 연결성")
    print("  - Glycolysis: 10/10 반응 존재 [OK]")
    print("  - TCA Cycle: 8/8 반응 존재 [OK]")
    print("  - Glyoxylate Shunt: 2/2 반응 존재 [OK]")
    print("  - 주요 대사물질: 모두 존재 및 연결됨 [OK]")
    
    print("\n3. 모델 구조 검토")
    print("  - 무제한 영양소: 성장 가능 (Biomass flux: 63.37 1/h) [OK]")
    print("  - 포도당만 허용: infeasible [FAIL]")
    print("  - Acetate만 허용: optimal이지만 Biomass flux = 0 [FAIL]")
    print("  - Acetate + 부트스트랩: 여전히 Biomass flux = 0 [FAIL]")
    
    print("\n[핵심 문제 발견]")
    print("\n문제 1: 뉴클레오티드 생산 경로")
    print("  - ATP (계수: -54.12) 생산 경로 작동 안함")
    print("  - GTP, UTP, CTP, dNTP 생산 경로 문제")
    print("  - 포도당만으로는 뉴클레오티드 생산 불가")
    print("  - 부트스트랩 추가해도 해결 안됨")
    
    print("\n문제 2: 부트스트랩 문제")
    print("  - CoA 생산 경로 작동 안함")
    print("  - Acetate -> Acetyl-CoA 경로가 CoA 필요")
    print("  - CoA 생산 경로가 Acetyl-CoA 필요 (순환 의존성)")
    
    print("\n문제 3: 뉴클레오티드 생합성 경로 불완전 가능성")
    print("  - 포도당으로도 성장하지 못함")
    print("  - 무제한 영양소로는 성장 가능")
    print("  - -> 탄소원 처리 경로는 정상")
    print("  - -> 뉴클레오티드/아미노산 생합성 경로 문제 가능성")
    
    print("\n[추가 조사 필요 항목]")
    print("\n1. 뉴클레오티드 생합성 경로 확인")
    print("  - De novo purine synthesis 경로")
    print("  - De novo pyrimidine synthesis 경로")
    print("  - Nucleotide salvage pathway")
    print("  - dNTP 합성 경로")
    
    print("\n2. 아미노산 생합성 경로 확인")
    print("  - 20종 아미노산 생합성 경로")
    print("  - 특히 필수 아미노산 경로")
    
    print("\n3. 보조인자 생산 경로 확인")
    print("  - CoA 생산 경로 (Pantothenate -> CoA)")
    print("  - NAD+/NADP+ 생산 경로")
    print("  - FAD 생산 경로")
    print("  - Folate 관련 경로")
    
    print("\n[현재 상태 요약]")
    print("  - 모델 구조: 정상 (무제한 영양소로 성장 가능)")
    print("  - TCA/Glyoxylate 경로: 완전")
    print("  - 문제 영역: 뉴클레오티드 및 아미노산 생합성 경로")
    print("  - 부트스트랩: 일부 해결되지만 여전히 성장 불가")
    
    print("\n" + "="*70)
    print("진단 완료")
    print("="*70)

def main():
    print("="*70)
    print("모델 종합 진단 보고서")
    print("="*70)
    
    create_diagnosis_report()
    
    # 파일로 저장할 수 있도록 요약
    summary = {
        '항목': [
            '모델 구조',
            'Glycolysis',
            'TCA Cycle',
            'Glyoxylate Shunt',
            '무제한 영양소 성장',
            '포도당 기반 성장',
            'Acetate 기반 성장',
            '뉴클레오티드 생산',
            'CoA 생산'
        ],
        '상태': [
            '정상',
            '10/10 반응 존재',
            '8/8 반응 존재',
            '2/2 반응 존재',
            '가능 (63.37 1/h)',
            '불가 (infeasible)',
            '불가 (flux=0)',
            '경로 문제',
            '부트스트랩 필요'
        ]
    }
    
    df_summary = pd.DataFrame(summary)
    df_summary.to_csv('model_diagnosis_summary.csv', index=False)
    print(f"\n[OK] 요약 저장: model_diagnosis_summary.csv")

if __name__ == "__main__":
    main()
