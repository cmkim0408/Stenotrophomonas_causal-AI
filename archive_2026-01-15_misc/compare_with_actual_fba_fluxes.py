#!/usr/bin/env python
"""
레퍼런스 모델의 실제 FBA 플럭스를 기반으로 비교
레퍼런스 모델의 FBA 결과에서 실제로 사용된 반응들을 확인하고,
신규 모델과 비교하여 누락된 반응을 재분석합니다.
"""

import cobra
from pathlib import Path
import pandas as pd
import numpy as np

def load_model(model_path):
    """모델 로드"""
    print(f"모델 로드 중: {model_path}")
    model = cobra.io.read_sbml_model(str(model_path))
    print(f"[OK] 모델 로드 완료: {model.id}")
    return model

def get_active_reactions_from_fba(flux_file, threshold=1e-6):
    """FBA 플럭스 파일에서 실제로 사용된 반응 추출"""
    print(f"\nFBA 플럭스 파일 분석: {flux_file}")
    
    df = pd.read_csv(flux_file, index_col=0)
    
    # 모든 컬럼에 대해 절댓값이 threshold보다 큰 반응 찾기
    active_mask = (df.abs() > threshold).any(axis=1)
    active_reactions = df[active_mask].index.tolist()
    
    print(f"  총 반응 수: {len(df)}")
    print(f"  실제 사용된 반응 수 (|flux| > {threshold}): {len(active_reactions)}")
    
    return set(active_reactions)

def compare_with_active_reactions(ref_model, new_model, active_rxn_ids):
    """실제 사용된 반응만 비교"""
    print("\n" + "="*70)
    print("실제 FBA에서 사용된 반응 비교")
    print("="*70)
    
    # 레퍼런스 모델에 있는 반응 중 실제 사용된 것만
    ref_active_in_model = active_rxn_ids & set(ref_model.reactions.list_attr('id'))
    new_rxns = set(new_model.reactions.list_attr('id'))
    
    # 실제 사용되었지만 신규 모델에 없는 반응
    missing_active = ref_active_in_model - new_rxns
    
    print(f"\n레퍼런스 모델에서 실제 사용된 반응: {len(ref_active_in_model)}개")
    print(f"그 중 신규 모델에 없는 반응: {len(missing_active)}개")
    
    return missing_active

def analyze_missing_active_reactions(ref_model, missing_rxn_ids):
    """누락된 실제 사용 반응 상세 분석"""
    print("\n" + "="*70)
    print("누락된 실제 사용 반응 상세 분석")
    print("="*70)
    
    missing_reactions = []
    
    for rxn_id in sorted(missing_rxn_ids):
        if rxn_id in ref_model.reactions:
            rxn = ref_model.reactions.get_by_id(rxn_id)
            
            # 경로 분류
            pathway = classify_pathway(rxn_id, rxn)
            
            missing_reactions.append({
                'reaction_id': rxn_id,
                'name': rxn.name,
                'equation': rxn.reaction,
                'genes': ', '.join([g.id for g in rxn.genes]),
                'pathway': pathway,
                'lower_bound': rxn.lower_bound,
                'upper_bound': rxn.upper_bound,
                'reversible': rxn.lower_bound < 0
            })
    
    return missing_reactions

def classify_pathway(rxn_id, reaction):
    """반응의 경로 분류"""
    rxn_id_lower = rxn_id.lower()
    name_lower = reaction.name.lower() if reaction.name else ""
    
    # Exchange reactions
    if rxn_id.startswith('EX_'):
        return 'Exchange'
    
    # Transport reactions
    if 'transport' in name_lower or '_t' in rxn_id_lower or '_pp' in rxn_id_lower or '_p' in rxn_id_lower:
        return 'Transport'
    
    # Acetate pathway
    if any(keyword in rxn_id_lower or keyword in name_lower for keyword in ['acs', 'ack', 'pta', 'aceti']):
        return 'Acetate Metabolism'
    
    # TCA cycle
    if any(keyword in rxn_id_lower or keyword in name_lower for keyword in ['citrate', 'aconit', 'icdh', 'akgdh', 'succ', 'fumar', 'malate', 'mdh']):
        return 'TCA Cycle'
    
    # Glyoxylate shunt
    if any(keyword in rxn_id_lower or keyword in name_lower for keyword in ['icl', 'mals', 'glyoxylate']):
        return 'Glyoxylate Shunt'
    
    # Gluconeogenesis
    if any(keyword in rxn_id_lower or keyword in name_lower for keyword in ['pepck', 'ppdk', 'pc', 'fructose', 'glucose', 'gluconeogenesis']):
        return 'Gluconeogenesis'
    
    # Glycolysis
    if any(keyword in rxn_id_lower or keyword in name_lower for keyword in ['glycolysis', 'hexokinase', 'phosphofructo', 'aldolase', 'glyceraldehyde']):
        return 'Glycolysis'
    
    # PPP
    if any(keyword in rxn_id_lower or keyword in name_lower for keyword in ['pentose', 'transketolase', 'transaldolase', 'ribose', 'xylulose']):
        return 'PPP'
    
    # ETC / Electron transport
    if any(keyword in rxn_id_lower or keyword in name_lower for keyword in ['nadl', 'cytbo', 'atps', 'succinate dehydrogenase', 'fumarate reductase']):
        return 'ETC / Electron Transport'
    
    # Amino acid
    aa_keywords = ['serine', 'glycine', 'tyrosine', 'phenylalanine', 'tryptophan', 
                   'leucine', 'isoleucine', 'valine', 'methionine', 'cysteine',
                   'aspartate', 'asparagine', 'glutamate', 'glutamine', 'lysine',
                   'arginine', 'histidine', 'proline', 'threonine', 'alanine']
    if any(keyword in name_lower for keyword in aa_keywords):
        return 'Amino Acid Synthesis'
    
    # Nucleotide
    if any(keyword in rxn_id_lower or keyword in name_lower for keyword in ['nucleotide', 'atp', 'gtp', 'utp', 'ctp', 'adp', 'gdp', 'udp', 'cdp']):
        return 'Nucleotide Synthesis'
    
    # Cofactor
    if any(keyword in rxn_id_lower or keyword in name_lower for keyword in ['coa', 'nad', 'nadp', 'fad', 'fmn', 'thf', 'folate', 'biotin', 'pantothenate', 'menaquinone', 'ubiquinone']):
        return 'Cofactor Synthesis'
    
    return 'Other'

def main():
    # 경로 설정
    base_path = Path(__file__).parent.parent
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    new_model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    flux_file = base_path / "Stenotrophomonas" / "fba_flux_gradient_acid.csv"
    
    # 모델 로드
    print("="*70)
    print("레퍼런스 모델 FBA 결과 기반 비교")
    print("="*70)
    
    ref_model = load_model(ref_model_path)
    new_model = load_model(new_model_path)
    
    # 실제 사용된 반응 추출
    active_rxn_ids = get_active_reactions_from_fba(flux_file, threshold=1e-6)
    
    # 비교
    missing_active = compare_with_active_reactions(ref_model, new_model, active_rxn_ids)
    
    # 상세 분석
    missing_reactions = analyze_missing_active_reactions(ref_model, missing_active)
    
    # 결과를 DataFrame으로 변환
    df_missing = pd.DataFrame(missing_reactions)
    
    # CSV 파일로 저장
    script_dir = Path(__file__).parent
    output_file = script_dir / "missing_active_reactions.csv"
    df_missing.to_csv(output_file, index=False, encoding='utf-8-sig')
    print(f"\n[OK] 누락된 실제 사용 반응 리스트 저장: {output_file}")
    
    # 경로별 요약
    print("\n" + "="*70)
    print("경로별 누락된 실제 사용 반응 요약")
    print("="*70)
    
    if len(df_missing) > 0:
        pathway_summary = df_missing.groupby('pathway').size().sort_values(ascending=False)
        for pathway, count in pathway_summary.items():
            print(f"  {pathway}: {count}개")
        
        # 상위 30개 반응 출력
        print("\n" + "="*70)
        print("누락된 실제 사용 반응 상위 30개")
        print("="*70)
        
        for i, (_, row) in enumerate(df_missing.head(30).iterrows(), 1):
            print(f"\n{i}. {row['reaction_id']} ({row['pathway']})")
            if pd.notna(row['name']):
                print(f"   이름: {row['name']}")
            print(f"   반응식: {row['equation']}")
            if pd.notna(row['genes']) and row['genes']:
                print(f"   유전자: {row['genes']}")
    
    # 전체 요약
    print("\n" + "="*70)
    print("비교 요약")
    print("="*70)
    print(f"\n레퍼런스 모델에서 실제 사용된 반응: {len(active_rxn_ids)}개")
    print(f"그 중 신규 모델에 없는 반응: {len(missing_active)}개")
    
    print("\n[OK] 비교 완료!")

if __name__ == "__main__":
    main()
