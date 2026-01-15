#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
원래 annotation된 모델에서 추가된 모든 반응 리스트업

비교:
- 원본: BaseModel.xml
- 최종: BaseModel_with_BCAA_cofactors_ions_nad_transport.xml
"""

import cobra
from pathlib import Path
import pandas as pd
import sys

def load_model(model_path):
    try:
        model = cobra.io.read_sbml_model(str(model_path))
        print(f"[OK] 모델 로드 완료: {model_path.name}")
        print(f"    반응 수: {len(model.reactions)}")
        return model
    except Exception as e:
        print(f"[ERROR] 모델 로드 실패: {e}")
        sys.exit(1)

def get_reaction_info(rxn):
    """반응 정보 추출"""
    return {
        'reaction_id': rxn.id,
        'reaction_name': rxn.name if rxn.name else '',
        'reaction_equation': rxn.reaction,
        'lower_bound': rxn.lower_bound,
        'upper_bound': rxn.upper_bound,
        'reversible': rxn.reversibility,
    }

def categorize_reaction(rxn_id, rxn_name, rxn_equation):
    """반응 카테고리 분류"""
    category = "Other"
    
    # BCAA 관련
    if any(x in rxn_id.upper() for x in ['KARI', 'DHAD', 'IPMI', 'IPMDH', 'BCAT_VAL', 'BCAT_LEU']):
        category = "BCAA Synthesis"
    elif 'BCAA' in rxn_name.upper() or 'VAL' in rxn_id.upper() or 'LEU' in rxn_id.upper():
        category = "BCAA Synthesis"
    
    # 이온 수송
    elif rxn_id.startswith('T_') and ('_e_to_' in rxn_id or '_c_to_' in rxn_id):
        if 'cl' in rxn_id.lower():
            category = "Ion Transport (Chloride)"
        elif 'cu' in rxn_id.lower():
            category = "Ion Transport (Copper)"
        elif 'cobalt' in rxn_id.lower():
            category = "Ion Transport (Cobalt)"
        else:
            category = "Ion Transport"
    
    # NAD/NADP 관련
    elif rxn_id.startswith('EX_nac') or rxn_id.startswith('EX_ncam'):
        category = "NAD/NADP Exchange"
    elif 'T_nac' in rxn_id.lower() or 'nac_e_to_nac_c' in rxn_id.lower():
        category = "NAD/NADP Transport"
    elif 'NAD' in rxn_id or 'NADP' in rxn_id:
        category = "NAD/NADP Synthesis"
    
    # Acetate 관련
    elif rxn_id == 'ACS_ADP':
        category = "Acetate Metabolism"
    elif rxn_id == 'ACtexi':
        category = "Acetate Transport"
    
    # TCA cycle
    elif rxn_id == 'SUCDi':
        category = "TCA Cycle"
    
    # Gluconeogenesis
    elif rxn_id == 'PEPCK_ATP':
        category = "Gluconeogenesis"
    
    # Exchange 반응
    elif rxn_id.startswith('EX_'):
        category = "Exchange"
    
    # Transport 반응
    elif rxn_id.startswith('T_'):
        category = "Transport"
    
    return category

def main():
    base_path = Path(__file__).parent.parent
    original_model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    final_model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_with_BCAA_cofactors_ions_nad_transport.xml"
    output_dir = base_path / "Stenotrophomonas-causal AI"
    
    print("="*70)
    print("추가된 반응 리스트업")
    print("="*70)
    
    # 모델 로드
    print("\n[모델 로드]")
    original_model = load_model(original_model_path)
    final_model = load_model(final_model_path)
    
    # 반응 ID 비교
    original_rxn_ids = set(r.id for r in original_model.reactions)
    final_rxn_ids = set(r.id for r in final_model.reactions)
    
    added_rxn_ids = final_rxn_ids - original_rxn_ids
    
    print(f"\n[반응 비교]")
    print(f"  원본 모델 반응 수: {len(original_rxn_ids)}")
    print(f"  최종 모델 반응 수: {len(final_rxn_ids)}")
    print(f"  추가된 반응 수: {len(added_rxn_ids)}")
    
    # 추가된 반응 정보 수집
    added_reactions = []
    for rxn_id in sorted(added_rxn_ids):
        rxn = final_model.reactions.get_by_id(rxn_id)
        info = get_reaction_info(rxn)
        info['category'] = categorize_reaction(rxn_id, rxn.name if rxn.name else '', rxn.reaction)
        added_reactions.append(info)
    
    # 카테고리별로 정리
    categories = {}
    for rxn_info in added_reactions:
        cat = rxn_info['category']
        if cat not in categories:
            categories[cat] = []
        categories[cat].append(rxn_info)
    
    # 출력
    print("\n" + "="*70)
    print("추가된 반응 목록 (카테고리별)")
    print("="*70)
    
    for category in sorted(categories.keys()):
        rxns = categories[category]
        print(f"\n[{category}] ({len(rxns)}개)")
        print("-" * 70)
        for rxn_info in rxns:
            print(f"\n  ID: {rxn_info['reaction_id']}")
            if rxn_info['reaction_name']:
                print(f"  이름: {rxn_info['reaction_name']}")
            print(f"  반응식: {rxn_info['reaction_equation']}")
            print(f"  bounds: [{rxn_info['lower_bound']}, {rxn_info['upper_bound']}]")
            print(f"  가역성: {rxn_info['reversible']}")
    
    # DataFrame 생성 및 CSV 저장
    df = pd.DataFrame(added_reactions)
    df = df.sort_values(['category', 'reaction_id'])
    
    output_file = output_dir / "added_reactions_list.csv"
    df.to_csv(output_file, index=False, encoding='utf-8-sig')
    print(f"\n\n[CSV 저장] {output_file}")
    print(f"  총 {len(df)}개 반응")
    
    # 카테고리별 요약
    print("\n" + "="*70)
    print("카테고리별 요약")
    print("="*70)
    category_summary = df.groupby('category').size().sort_values(ascending=False)
    for category, count in category_summary.items():
        print(f"  {category}: {count}개")
    
    # 상세 표시용 마크다운 파일도 생성
    md_file = output_dir / "added_reactions_list.md"
    with open(md_file, 'w', encoding='utf-8') as f:
        f.write("# 추가된 반응 목록\n\n")
        f.write(f"**원본 모델**: BaseModel.xml ({len(original_rxn_ids)}개 반응)\n")
        f.write(f"**최종 모델**: BaseModel_with_BCAA_cofactors_ions_nad_transport.xml ({len(final_rxn_ids)}개 반응)\n")
        f.write(f"**추가된 반응**: {len(added_rxn_ids)}개\n\n")
        f.write("---\n\n")
        
        for category in sorted(categories.keys()):
            rxns = categories[category]
            f.write(f"## {category} ({len(rxns)}개)\n\n")
            for rxn_info in rxns:
                f.write(f"### {rxn_info['reaction_id']}\n\n")
                if rxn_info['reaction_name']:
                    f.write(f"**이름**: {rxn_info['reaction_name']}\n\n")
                f.write(f"**반응식**: `{rxn_info['reaction_equation']}`\n\n")
                f.write(f"**bounds**: [{rxn_info['lower_bound']}, {rxn_info['upper_bound']}]\n\n")
                f.write(f"**가역성**: {rxn_info['reversible']}\n\n")
                f.write("---\n\n")
    
    print(f"[Markdown 저장] {md_file}")
    
    return df, categories

if __name__ == "__main__":
    df, categories = main()
