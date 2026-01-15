#!/usr/bin/env python
"""
레퍼런스 모델의 실제 FBA 플럭스에서 사용된 모든 반응 확인
신규 모델과 비교하여 누락된 모든 반응 찾기
"""

import cobra
from pathlib import Path
import pandas as pd
import numpy as np
from collections import defaultdict

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
    
    # 각 반응의 최대 플럭스 값도 저장
    max_fluxes = {}
    for rxn_id in active_reactions:
        max_flux = df.loc[rxn_id].abs().max()
        max_fluxes[rxn_id] = max_flux
    
    print(f"  총 반응 수: {len(df)}")
    print(f"  실제 사용된 반응 수 (|flux| > {threshold}): {len(active_reactions)}")
    
    return set(active_reactions), max_fluxes

def classify_pathway(rxn_id, reaction, max_flux=None):
    """반응의 경로 분류"""
    rxn_id_lower = rxn_id.lower()
    name_lower = reaction.name.lower() if reaction.name else ""
    
    # Exchange reactions
    if rxn_id.startswith('EX_'):
        return 'Exchange'
    
    # Transport reactions
    if rxn_id.startswith('T_') or 'transport' in name_lower or '_pp' in rxn_id_lower or '_p' in rxn_id_lower:
        return 'Transport'
    
    # Acetate pathway
    if any(keyword in rxn_id_lower or keyword in name_lower for keyword in ['acs', 'ack', 'pta', 'aceti']):
        return 'Acetate Metabolism'
    
    # TCA cycle
    tca_keywords = ['citrate', 'aconit', 'icdh', 'akgdh', 'succ', 'fumar', 'malate', 'mdh']
    if any(keyword in rxn_id_lower or keyword in name_lower for keyword in tca_keywords):
        return 'TCA Cycle'
    
    # Glyoxylate shunt
    if any(keyword in rxn_id_lower or keyword in name_lower for keyword in ['icl', 'mals', 'glyoxylate']):
        return 'Glyoxylate Shunt'
    
    # Gluconeogenesis
    if any(keyword in rxn_id_lower or keyword in name_lower for keyword in ['pepck', 'ppdk', 'pc', 'fructose', 'glucose', 'gluconeogenesis', 'fba', 'fbp', 'pgi', 'tpi']):
        return 'Gluconeogenesis / Glycolysis'
    
    # PPP
    if any(keyword in rxn_id_lower or keyword in name_lower for keyword in ['pentose', 'transketolase', 'transaldolase', 'ribose', 'xylulose', 'g6pdh', 'gnd', 'rpi', 'rpe']):
        return 'PPP'
    
    # ETC / Electron transport
    etc_keywords = ['nadl', 'nad16', 'cytbo', 'atps', 'succinate dehydrogenase', 'fumarate reductase', 'menaquinone', 'ubiquinone', 'q8']
    if any(keyword in rxn_id_lower or keyword in name_lower for keyword in etc_keywords):
        return 'ETC / Electron Transport'
    
    # Amino acid
    aa_keywords = ['serine', 'glycine', 'tyrosine', 'phenylalanine', 'tryptophan', 
                   'leucine', 'isoleucine', 'valine', 'methionine', 'cysteine',
                   'aspartate', 'asparagine', 'glutamate', 'glutamine', 'lysine',
                   'arginine', 'histidine', 'proline', 'threonine', 'alanine']
    if any(keyword in name_lower for keyword in aa_keywords):
        return 'Amino Acid Metabolism'
    
    # Nucleotide
    if any(keyword in rxn_id_lower or keyword in name_lower for keyword in ['nucleotide', 'atp', 'gtp', 'utp', 'ctp', 'adp', 'gdp', 'udp', 'cdp', 'dntp', 'dctp', 'datp', 'dgtp', 'dttp']):
        return 'Nucleotide Metabolism'
    
    # Cofactor
    cofactor_keywords = ['coa', 'nad', 'nadp', 'fad', 'fmn', 'thf', 'folate', 'biotin', 'pantothenate', 'riboflavin', 'thiamine', 'pyridoxine']
    if any(keyword in rxn_id_lower or keyword in name_lower for keyword in cofactor_keywords):
        return 'Cofactor Metabolism'
    
    # Lipid
    if any(keyword in rxn_id_lower or keyword in name_lower for keyword in ['fatty', 'lipid', 'acyl', 'glycerol']):
        return 'Lipid Metabolism'
    
    # Biomass
    if 'biomass' in rxn_id_lower or 'growth' in rxn_id_lower:
        return 'Biomass'
    
    return 'Other'

def compare_active_reactions(ref_model, new_model, active_rxn_ids, max_fluxes):
    """실제 사용된 반응 비교 및 상세 분석"""
    print("\n" + "="*70)
    print("실제 사용된 반응 비교 및 분석")
    print("="*70)
    
    ref_rxns = set(ref_model.reactions.list_attr('id'))
    new_rxns = set(new_model.reactions.list_attr('id'))
    
    # 레퍼런스 모델에 있고 실제 사용된 반응
    ref_active_in_model = active_rxn_ids & ref_rxns
    
    # 그 중 신규 모델에 없는 반응
    missing_active = ref_active_in_model - new_rxns
    
    print(f"\n레퍼런스 모델에서 실제 사용된 반응: {len(ref_active_in_model)}개")
    print(f"그 중 신규 모델에 없는 반응: {len(missing_active)}개")
    
    # 누락된 반응 상세 분석
    missing_reactions = []
    
    for rxn_id in sorted(missing_active):
        if rxn_id in ref_model.reactions:
            rxn = ref_model.reactions.get_by_id(rxn_id)
            max_flux = max_fluxes.get(rxn_id, 0)
            pathway = classify_pathway(rxn_id, rxn, max_flux)
            
            missing_reactions.append({
                'reaction_id': rxn_id,
                'name': rxn.name if rxn.name else '',
                'equation': rxn.reaction,
                'genes': ', '.join([g.id for g in rxn.genes]),
                'pathway': pathway,
                'max_flux': max_flux,
                'lower_bound': rxn.lower_bound,
                'upper_bound': rxn.upper_bound,
                'reversible': rxn.lower_bound < 0
            })
    
    return missing_reactions

def prioritize_reactions(missing_reactions):
    """반응 우선순위 정리"""
    # 플럭스 크기와 경로 기반 우선순위
    high_priority_pathways = [
        'Acetate Metabolism',
        'TCA Cycle',
        'Glyoxylate Shunt',
        'Gluconeogenesis / Glycolysis',
        'ETC / Electron Transport'
    ]
    
    medium_priority_pathways = [
        'PPP',
        'Amino Acid Metabolism',
        'Nucleotide Metabolism',
        'Cofactor Metabolism'
    ]
    
    for rxn in missing_reactions:
        pathway = rxn['pathway']
        max_flux = rxn['max_flux']
        
        if pathway in high_priority_pathways:
            if max_flux > 0.01:
                rxn['priority'] = 'HIGH'
            else:
                rxn['priority'] = 'MEDIUM'
        elif pathway in medium_priority_pathways:
            rxn['priority'] = 'MEDIUM'
        elif pathway in ['Exchange', 'Transport']:
            rxn['priority'] = 'MEDIUM'
        else:
            rxn['priority'] = 'LOW'
    
    return missing_reactions

def main():
    # 경로 설정
    base_path = Path(__file__).parent.parent
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    new_model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    flux_file = base_path / "Stenotrophomonas" / "fba_flux_gradient_acid.csv"
    
    print("="*70)
    print("레퍼런스 모델 FBA 결과 기반 전체 누락 반응 분석")
    print("="*70)
    
    # 모델 로드
    ref_model = load_model(ref_model_path)
    new_model = load_model(new_model_path)
    
    # 실제 사용된 반응 추출
    active_rxn_ids, max_fluxes = get_active_reactions_from_fba(flux_file, threshold=1e-6)
    
    # 비교 및 분석
    missing_reactions = compare_active_reactions(ref_model, new_model, active_rxn_ids, max_fluxes)
    
    # 우선순위 정리
    missing_reactions = prioritize_reactions(missing_reactions)
    
    # DataFrame으로 변환
    df_missing = pd.DataFrame(missing_reactions)
    
    # 우선순위별 정렬
    priority_order = {'HIGH': 1, 'MEDIUM': 2, 'LOW': 3}
    df_missing['priority_order'] = df_missing['priority'].map(priority_order)
    df_missing = df_missing.sort_values(['priority_order', 'max_flux', 'pathway', 'reaction_id'], ascending=[True, False, True, True])
    df_missing = df_missing.drop('priority_order', axis=1)
    
    # CSV 파일로 저장
    script_dir = Path(__file__).parent
    output_file = script_dir / "all_missing_active_reactions.csv"
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
        
        # 우선순위별 요약
        print("\n" + "="*70)
        print("우선순위별 요약")
        print("="*70)
        priority_summary = df_missing.groupby('priority').size()
        for priority, count in priority_summary.items():
            print(f"  {priority}: {count}개")
        
        # HIGH 우선순위 반응 출력
        print("\n" + "="*70)
        print("HIGH 우선순위 누락 반응")
        print("="*70)
        
        high_priority = df_missing[df_missing['priority'] == 'HIGH']
        for i, (_, row) in enumerate(high_priority.iterrows(), 1):
            print(f"\n{i}. {row['reaction_id']} ({row['pathway']})")
            if pd.notna(row['name']) and row['name']:
                print(f"   이름: {row['name']}")
            print(f"   반응식: {row['equation']}")
            print(f"   최대 플럭스: {row['max_flux']:.6f}")
            if pd.notna(row['genes']) and row['genes']:
                print(f"   유전자: {row['genes']}")
        
        # MEDIUM 우선순위 일부 출력
        print("\n" + "="*70)
        print("MEDIUM 우선순위 누락 반응 (상위 20개)")
        print("="*70)
        
        medium_priority = df_missing[df_missing['priority'] == 'MEDIUM'].head(20)
        for i, (_, row) in enumerate(medium_priority.iterrows(), 1):
            print(f"\n{i}. {row['reaction_id']} ({row['pathway']})")
            if pd.notna(row['name']) and row['name']:
                print(f"   이름: {row['name']}")
            print(f"   반응식: {row['equation']}")
            print(f"   최대 플럭스: {row['max_flux']:.6f}")
    
    # 전체 요약
    print("\n" + "="*70)
    print("전체 요약")
    print("="*70)
    print(f"레퍼런스 모델에서 실제 사용된 반응: {len(active_rxn_ids)}개")
    print(f"그 중 신규 모델에 없는 반응: {len(missing_reactions)}개")
    
    print("\n[OK] 분석 완료!")

if __name__ == "__main__":
    main()
