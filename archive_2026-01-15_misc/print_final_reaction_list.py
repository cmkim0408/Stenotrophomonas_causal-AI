#!/usr/bin/env python
"""
최종 반응 리스트 출력 (이름 포함, Pseudo 반응 제외)
"""

import pandas as pd
from pathlib import Path

def main():
    base_path = Path(__file__).parent.parent
    input_csv = base_path / "Stenotrophomonas-causal AI" / "missing_metabolic_reactions_detailed.csv"
    
    # CSV 읽기
    df = pd.read_csv(input_csv)
    
    # Pseudo 반응 제외
    df_real = df[~df['is_pseudo']].copy()
    
    # 정렬
    priority_order = {'HIGH': 1, 'MEDIUM': 2, 'LOW': 3}
    df_real['priority_num'] = df_real['priority'].map(priority_order)
    df_real = df_real.sort_values(['priority_num', 'pathway', 'reaction_id'])
    df_real = df_real.drop('priority_num', axis=1)
    
    print("="*80)
    print("레퍼런스 모델에 있지만 신규 모델에 없는 대사 경로 반응 리스트")
    print("(Transport, Exchange 및 Pseudo 반응 제외)")
    print("="*80)
    print(f"\n총 {len(df_real)}개 실제 대사 반응\n")
    
    # Pseudo 반응 리스트
    df_pseudo = df[df['is_pseudo']]
    if len(df_pseudo) > 0:
        print(f"제외된 Pseudo 반응 ({len(df_pseudo)}개):")
        for _, row in df_pseudo.iterrows():
            print(f"  - {row['reaction_id']}: {row['equation']}")
        print()
    
    # HIGH 우선순위
    df_high = df_real[df_real['priority'] == 'HIGH']
    if len(df_high) > 0:
        print("="*80)
        print("HIGH 우선순위 반응 (실제 사용됨) - 3개")
        print("="*80)
        for i, (_, row) in enumerate(df_high.iterrows(), 1):
            print(f"\n{i}. {row['reaction_id']}")
            if pd.notna(row['name']) and row['name']:
                print(f"   이름: {row['name']}")
            print(f"   반응식: {row['equation']}")
            print(f"   경로: {row['pathway']}")
            if row['is_active']:
                print(f"   실제 사용됨: Yes (플럭스: {row['max_flux']:.6f})")
            print()
    
    # MEDIUM 우선순위
    df_medium = df_real[df_real['priority'] == 'MEDIUM']
    if len(df_medium) > 0:
        print("\n" + "="*80)
        print("MEDIUM 우선순위 반응 - 12개")
        print("="*80)
        for i, (_, row) in enumerate(df_medium.iterrows(), 1):
            print(f"\n{i}. {row['reaction_id']}")
            if pd.notna(row['name']) and row['name']:
                print(f"   이름: {row['name']}")
            print(f"   반응식: {row['equation']}")
            print(f"   경로: {row['pathway']}")
            print()
    
    # LOW 우선순위
    df_low = df_real[df_real['priority'] == 'LOW']
    if len(df_low) > 0:
        print("\n" + "="*80)
        print("LOW 우선순위 반응 - 10개")
        print("="*80)
        for i, (_, row) in enumerate(df_low.iterrows(), 1):
            print(f"\n{i}. {row['reaction_id']}")
            if pd.notna(row['name']) and row['name']:
                print(f"   이름: {row['name']}")
            print(f"   반응식: {row['equation']}")
            print(f"   경로: {row['pathway']}")
            print()
    
    print("\n" + "="*80)
    print("요약")
    print("="*80)
    print(f"\n총 실제 대사 반응: {len(df_real)}개")
    print(f"  HIGH: {len(df_high)}개")
    print(f"  MEDIUM: {len(df_medium)}개")
    print(f"  LOW: {len(df_low)}개")
    print(f"\n제외된 Pseudo 반응: {len(df_pseudo)}개")

if __name__ == "__main__":
    main()
