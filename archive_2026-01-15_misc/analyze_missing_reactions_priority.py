#!/usr/bin/env python
"""
누락된 반응 우선순위 분석 및 정리
레퍼런스 모델에 있지만 신규 모델에 없는 반응 중
FBA 성공에 중요한 반응들을 우선순위로 정리합니다.
"""

import pandas as pd
from pathlib import Path

def analyze_missing_reactions():
    """누락된 반응 우선순위 분석"""
    
    # CSV 파일 읽기
    script_dir = Path(__file__).parent
    csv_file = script_dir / "missing_reactions_vs_reference.csv"
    df = pd.read_csv(csv_file)
    
    # 우선순위 정의
    high_priority_reactions = {
        # Central Carbon Metabolism - PEP 생산
        'PC': {
            'priority': 'HIGH',
            'category': 'Central Carbon Metabolism',
            'description': 'Pyruvate carboxylase - Pyruvate + HCO3 → OAA (Gluconeogenesis 핵심)',
            'reason': 'PEP 생산에 필요, OAA 생산 경로'
        },
        'PPDK': {
            'priority': 'HIGH',
            'category': 'Central Carbon Metabolism',
            'description': 'Pyruvate phosphate dikinase - Pyruvate → PEP (대체 경로)',
            'reason': 'PEP 생산 대체 경로'
        },
        'PEPCK': {
            'priority': 'HIGH',
            'category': 'Gluconeogenesis',
            'description': 'PEP carboxykinase (GTP) - OAA + GTP → PEP + CO2 + GDP',
            'reason': 'PEP 생산 핵심 반응'
        },
        'PEPCK_ATP': {
            'priority': 'HIGH',
            'category': 'Gluconeogenesis',
            'description': 'PEP carboxykinase (ATP) - OAA + ATP → PEP + CO2 + ADP',
            'reason': 'PEP 생산 대체 반응 (ATP 사용)'
        },
        
        # Acetyl-CoA 생산
        'ACS_ADP': {
            'priority': 'HIGH',
            'category': 'Acetyl-CoA Synthesis',
            'description': 'Acetate-CoA ligase (ADP-forming) - Acetate + ATP + CoA ⇄ AcCoA + ADP + Pi',
            'reason': 'Acetate에서 AcCoA 생산 핵심 반응'
        },
        
        # Cofactor - Menaquinone 관련 (ETC)
        'MKRED': {
            'priority': 'MEDIUM',
            'category': 'ETC / Electron Transport',
            'description': 'Menaquinone-8 reduction by NADH',
            'reason': '전자 전달 사슬 관련'
        },
        'MQN8RD': {
            'priority': 'MEDIUM',
            'category': 'ETC / Electron Transport',
            'description': 'Menaquinone-8 reduction by NADH (alternative)',
            'reason': '전자 전달 사슬 관련'
        },
        'MQN8r_NADH': {
            'priority': 'MEDIUM',
            'category': 'ETC / Electron Transport',
            'description': 'Menaquinone-8 reduction by NADH',
            'reason': '전자 전달 사슬 관련'
        },
        'MQN8red': {
            'priority': 'MEDIUM',
            'category': 'ETC / Electron Transport',
            'description': 'Menaquinone-8 NADH dehydrogenase (lumped)',
            'reason': '전자 전달 사슬 관련'
        },
        'NADH16_MQ': {
            'priority': 'MEDIUM',
            'category': 'ETC / Electron Transport',
            'description': 'NDH-1 to MQN8 - Complex I to menaquinone',
            'reason': 'Complex I 관련 전자 전달'
        },
        'MKOX': {
            'priority': 'MEDIUM',
            'category': 'ETC / Electron Transport',
            'description': 'Menaquinol-8 oxidation',
            'reason': '전자 전달 사슬 관련'
        },
        
        # Succinate 관련
        'SUCDi': {
            'priority': 'MEDIUM',
            'category': 'TCA Cycle',
            'description': 'Succinate dehydrogenase - Succinate + Q8 → Fumarate + Q8H2',
            'reason': 'TCA 사이클, 전자 전달 사슬'
        },
        'FRD7': {
            'priority': 'MEDIUM',
            'category': 'TCA Cycle',
            'description': 'Fumarate reductase - Fumarate + Q8H2 → Succinate + Q8',
            'reason': '무산소 조건에서 사용 가능'
        },
        
        # Cofactor synthesis
        'FADS': {
            'priority': 'MEDIUM',
            'category': 'Cofactor Synthesis',
            'description': 'FAD synthetase - ATP + FMN → FAD + PPI',
            'reason': 'FAD 보조인자 합성'
        },
        'RIBFLVKin': {
            'priority': 'MEDIUM',
            'category': 'Cofactor Synthesis',
            'description': 'Riboflavin kinase - ATP + Riboflavin → ADP + FMN + H',
            'reason': 'FMN 생합성 경로'
        },
        
        # Carbon anhydrase
        'CA': {
            'priority': 'MEDIUM',
            'category': 'Carbon Metabolism',
            'description': 'Carbonic anhydrase - CO2 + H2O ⇄ H + HCO3',
            'reason': 'HCO3 생산/소비 균형'
        },
        
        # Amino acid metabolism
        'BCAT_LEU': {
            'priority': 'LOW',
            'category': 'Amino Acid Metabolism',
            'description': 'Branched-chain amino acid transaminase (Leu)',
            'reason': 'Leucine 대사'
        },
        'BCAT_VAL': {
            'priority': 'LOW',
            'category': 'Amino Acid Metabolism',
            'description': 'Branched-chain amino acid transaminase (Val)',
            'reason': 'Valine 대사'
        },
        
        # Supply reactions (bootstrap)
        'SUPPLY_accoa_c': {
            'priority': 'LOW',
            'category': 'Bootstrap',
            'description': 'AcCoA supply reaction (bootstrap)',
            'reason': '부트스트랩 반응, 실제 경로가 필요'
        },
        'SUPPLY_oaa_c': {
            'priority': 'LOW',
            'category': 'Bootstrap',
            'description': 'OAA supply reaction (bootstrap)',
            'reason': '부트스트랩 반응, 실제 경로가 필요'
        },
    }
    
    # 결과 데이터프레임 생성
    results = []
    
    for _, row in df.iterrows():
        rxn_id = row['reaction_id']
        priority_info = high_priority_reactions.get(rxn_id, {})
        
        result = {
            'reaction_id': rxn_id,
            'name': row['name'],
            'equation': row['equation'],
            'genes': row['genes'],
            'pathway': row['pathway'],
            'priority': priority_info.get('priority', 'UNKNOWN'),
            'category': priority_info.get('category', row['pathway']),
            'description': priority_info.get('description', ''),
            'reason': priority_info.get('reason', ''),
            'lower_bound': row['lower_bound'],
            'upper_bound': row['upper_bound'],
            'reversible': row['reversible']
        }
        results.append(result)
    
    # 우선순위별 정렬
    priority_order = {'HIGH': 1, 'MEDIUM': 2, 'LOW': 3, 'UNKNOWN': 4}
    df_results = pd.DataFrame(results)
    df_results['priority_order'] = df_results['priority'].map(priority_order)
    df_results = df_results.sort_values(['priority_order', 'category', 'reaction_id'])
    df_results = df_results.drop('priority_order', axis=1)
    
    # CSV 저장
    script_dir = Path(__file__).parent
    output_file = script_dir / "missing_reactions_prioritized.csv"
    df_results.to_csv(output_file, index=False, encoding='utf-8-sig')
    print(f"[OK] 우선순위별 반응 리스트 저장: {output_file}")
    
    # 요약 출력
    print("\n" + "="*70)
    print("우선순위별 누락된 반응 요약")
    print("="*70)
    
    priority_summary = df_results.groupby('priority').size()
    print("\n우선순위별 개수:")
    for priority, count in priority_summary.items():
        print(f"  {priority}: {count}개")
    
    print("\n" + "="*70)
    print("HIGH 우선순위 반응 (FBA 성공에 핵심)")
    print("="*70)
    
    high_priority = df_results[df_results['priority'] == 'HIGH']
    print(f"\nHIGH 우선순위 반응 수: {len(high_priority)}개")
    for i, (_, row) in enumerate(high_priority.iterrows(), 1):
        print(f"\n{i}. {row['reaction_id']}")
        if pd.notna(row['name']):
            print(f"   Name: {row['name']}")
        print(f"   Equation: {row['equation']}")
        if row['description']:
            print(f"   Description: {row['description']}")
        if row['reason']:
            print(f"   Reason: {row['reason']}")
        if pd.notna(row['genes']) and row['genes']:
            print(f"   Genes: {row['genes']}")
    
    print("\n" + "="*70)
    print("MEDIUM 우선순위 반응 (FBA 성공에 중요)")
    print("="*70)
    
    medium_priority = df_results[df_results['priority'] == 'MEDIUM']
    print(f"\nMEDIUM 우선순위 반응 수: {len(medium_priority)}개")
    for i, (_, row) in enumerate(medium_priority.iterrows(), 1):
        print(f"\n{i}. {row['reaction_id']}")
        if pd.notna(row['name']):
            print(f"   Name: {row['name']}")
        print(f"   Equation: {row['equation']}")
        if row['description']:
            print(f"   Description: {row['description']}")
        if pd.notna(row['genes']) and row['genes']:
            print(f"   Genes: {row['genes']}")
    
    # Exchange 및 Transport 반응 제외한 핵심 반응만 필터링
    core_reactions = df_results[
        ~df_results['pathway'].isin(['Exchange', 'Transport']) & 
        (df_results['priority'].isin(['HIGH', 'MEDIUM']))
    ]
    
    core_output_file = script_dir / "missing_core_reactions_priority.csv"
    core_reactions.to_csv(core_output_file, index=False, encoding='utf-8-sig')
    print(f"\n[OK] 핵심 반응만 필터링한 리스트 저장: {core_output_file}")
    print(f"  (Exchange/Transport 제외, HIGH/MEDIUM 우선순위만)")
    print(f"  총 {len(core_reactions)}개 반응")
    
    return df_results, core_reactions

if __name__ == "__main__":
    df_all, df_core = analyze_missing_reactions()
    
    print("\n" + "="*70)
    print("분석 완료")
    print("="*70)
