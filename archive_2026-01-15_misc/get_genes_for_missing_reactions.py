#!/usr/bin/env python
"""
레퍼런스 모델에서 누락된 반응들의 유전자 정보 추출
"""

import cobra
import pandas as pd
from pathlib import Path

def load_model(model_path):
    print(f"모델 로드 중: {model_path}")
    model = cobra.io.read_sbml_model(model_path)
    print(f"[OK] 모델 로드 완료: {model.id}")
    return model

def get_reaction_genes(model, rxn_id):
    """반응의 유전자 정보 가져오기"""
    try:
        rxn = model.reactions.get_by_id(rxn_id)
        genes = [g.id for g in rxn.genes]
        gpr = str(rxn.gene_reaction_rule) if hasattr(rxn, 'gene_reaction_rule') else ''
        return genes, gpr
    except KeyError:
        return [], ''

def main():
    base_path = Path(__file__).parent.parent
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    input_csv = base_path / "Stenotrophomonas-causal AI" / "missing_metabolic_reactions_only.csv"
    output_csv = base_path / "Stenotrophomonas-causal AI" / "missing_metabolic_reactions_with_genes.csv"
    
    # 모델 로드
    ref_model = load_model(str(ref_model_path))
    
    # CSV 읽기
    df = pd.read_csv(input_csv)
    
    # 유전자 정보 추가
    genes_list = []
    gpr_list = []
    
    for _, row in df.iterrows():
        rxn_id = row['reaction_id']
        genes, gpr = get_reaction_genes(ref_model, rxn_id)
        genes_list.append(', '.join(genes) if genes else '')
        gpr_list.append(gpr)
    
    df['genes'] = genes_list
    df['gene_reaction_rule'] = gpr_list
    
    # CSV 저장
    df.to_csv(output_csv, index=False, encoding='utf-8-sig')
    print(f"\n[OK] 결과 저장: {output_csv}")
    
    # 요약 출력
    print("\n" + "="*70)
    print("유전자 정보 요약")
    print("="*70)
    
    total = len(df)
    with_genes = len(df[df['genes'] != ''])
    without_genes = total - with_genes
    
    print(f"\n총 {total}개 반응")
    print(f"  유전자 정보 있음: {with_genes}개")
    print(f"  유전자 정보 없음: {without_genes}개")
    
    print("\n유전자 정보가 있는 반응:")
    for _, row in df[df['genes'] != ''].iterrows():
        print(f"\n  {row['reaction_id']} ({row['pathway']}, {row['priority']})")
        print(f"    반응식: {row['equation']}")
        print(f"    유전자: {row['genes']}")
        if row['gene_reaction_rule']:
            print(f"    GPR: {row['gene_reaction_rule']}")

if __name__ == "__main__":
    main()
