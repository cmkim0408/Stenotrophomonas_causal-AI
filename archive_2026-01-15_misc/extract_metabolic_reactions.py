#!/usr/bin/env python
"""
Transport와 Exchange 제외한 대사 경로 반응만 추출
"""

import pandas as pd
from pathlib import Path

def main():
    base_path = Path(__file__).parent.parent
    input_file = base_path / "Stenotrophomonas-causal AI" / "comprehensive_missing_reactions.csv"
    output_file = base_path / "Stenotrophomonas-causal AI" / "missing_metabolic_reactions_only.csv"
    
    # CSV 읽기
    df = pd.read_csv(input_file)
    
    # Transport와 Exchange 제외
    df_metabolic = df[~df['pathway'].isin(['Transport', 'Exchange'])].copy()
    
    # 정렬 (우선순위, 경로, 반응 ID)
    priority_order = {'HIGH': 1, 'MEDIUM': 2, 'LOW': 3}
    df_metabolic['priority_order'] = df_metabolic['priority'].map(priority_order)
    df_metabolic = df_metabolic.sort_values(['priority_order', 'pathway', 'reaction_id'])
    df_metabolic = df_metabolic.drop('priority_order', axis=1)
    
    # CSV 저장
    df_metabolic.to_csv(output_file, index=False, encoding='utf-8-sig')
    
    print("="*70)
    print("Transport/Exchange 제외한 대사 경로 반응 리스트")
    print("="*70)
    print(f"\n총 {len(df_metabolic)}개 반응\n")
    
    # 경로별 분류
    print("경로별 분류:")
    pathway_summary = df_metabolic.groupby('pathway').size().sort_values(ascending=False)
    for pathway, count in pathway_summary.items():
        print(f"  {pathway}: {count}개")
    
    # 우선순위별 분류
    print("\n우선순위별 분류:")
    priority_summary = df_metabolic.groupby('priority').size()
    for priority, count in priority_summary.items():
        print(f"  {priority}: {count}개")
    
    # 실제 사용 여부
    print("\n실제 사용 여부:")
    active_summary = df_metabolic.groupby('is_active').size()
    print(f"  실제 사용됨: {active_summary.get(True, 0)}개")
    print(f"  실제 사용 안 됨: {active_summary.get(False, 0)}개")
    
    # 전체 반응 리스트
    print("\n" + "="*70)
    print("전체 반응 리스트")
    print("="*70)
    
    for i, (_, row) in enumerate(df_metabolic.iterrows(), 1):
        print(f"\n{i}. {row['reaction_id']} ({row['pathway']}, {row['priority']})")
        print(f"   반응식: {row['equation']}")
        if pd.notna(row['name']) and row['name']:
            print(f"   이름: {row['name']}")
        if row['is_active']:
            print(f"   실제 사용됨: Yes (플럭스: {row['max_flux']:.6f})")
        else:
            print(f"   실제 사용됨: No")
        if pd.notna(row['genes']) and row['genes']:
            print(f"   유전자: {row['genes']}")
    
    print(f"\n[OK] 결과 저장: {output_file}")
    print(f"총 {len(df_metabolic)}개 대사 경로 반응")

if __name__ == "__main__":
    main()
