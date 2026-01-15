#!/usr/bin/env python
"""
레퍼런스 모델과 신규 모델의 전체 비교
레퍼런스 모델의 실제 FBA 결과를 기반으로 누락된 반응을 찾되,
FBA 성공에 필요한 모든 반응을 포함하여 분석
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
    df = pd.read_csv(flux_file, index_col=0)
    active_mask = (df.abs() > threshold).any(axis=1)
    active_reactions = df[active_mask].index.tolist()
    max_fluxes = {}
    for rxn_id in active_reactions:
        max_flux = df.loc[rxn_id].abs().max()
        max_fluxes[rxn_id] = max_flux
    return set(active_reactions), max_fluxes

def compare_all_reactions(ref_model, new_model, active_rxn_ids, max_fluxes):
    """전체 반응 비교 및 분석"""
    print("\n" + "="*70)
    print("전체 반응 비교")
    print("="*70)
    
    ref_rxns = set(ref_model.reactions.list_attr('id'))
    new_rxns = set(new_model.reactions.list_attr('id'))
    
    # 레퍼런스에만 있는 반응
    ref_only = ref_rxns - new_rxns
    
    # 실제 사용된 반응 중 누락된 것
    active_missing = (active_rxn_ids & ref_rxns) - new_rxns
    
    print(f"\n레퍼런스 모델 반응 수: {len(ref_rxns)}")
    print(f"신규 모델 반응 수: {len(new_rxns)}")
    print(f"레퍼런스에만 있는 반응: {len(ref_only)}개")
    print(f"실제 사용된 반응 중 누락: {len(active_missing)}개")
    
    # 실제 사용된 반응들 중 신규 모델에 있는 것
    active_present = (active_rxn_ids & ref_rxns) & new_rxns
    print(f"실제 사용된 반응 중 신규 모델에 있음: {len(active_present)}개")
    
    return ref_only, active_missing, active_present

def analyze_reactions_in_detail(ref_model, missing_rxn_ids, active_rxn_ids, max_fluxes):
    """누락된 반응 상세 분석"""
    missing_reactions = []
    
    for rxn_id in sorted(missing_rxn_ids):
        if rxn_id in ref_model.reactions:
            rxn = ref_model.reactions.get_by_id(rxn_id)
            is_active = rxn_id in active_rxn_ids
            max_flux = max_fluxes.get(rxn_id, 0) if is_active else 0
            
            # 경로 분류
            pathway = classify_pathway(rxn_id, rxn)
            
            missing_reactions.append({
                'reaction_id': rxn_id,
                'name': rxn.name if rxn.name else '',
                'equation': rxn.reaction,
                'genes': ', '.join([g.id for g in rxn.genes]),
                'pathway': pathway,
                'is_active': is_active,
                'max_flux': max_flux,
                'lower_bound': rxn.lower_bound,
                'upper_bound': rxn.upper_bound,
                'reversible': rxn.lower_bound < 0
            })
    
    return missing_reactions

def classify_pathway(rxn_id, reaction):
    """반응의 경로 분류"""
    rxn_id_lower = rxn_id.lower()
    name_lower = reaction.name.lower() if reaction.name else ""
    
    if rxn_id.startswith('EX_'):
        return 'Exchange'
    if rxn_id.startswith('T_') or 'transport' in name_lower or '_pp' in rxn_id_lower:
        return 'Transport'
    if any(k in rxn_id_lower or k in name_lower for k in ['acs', 'ack', 'pta']):
        return 'Acetate Metabolism'
    if any(k in rxn_id_lower or k in name_lower for k in ['citrate', 'aconit', 'icdh', 'akgdh', 'succ', 'fumar', 'malate', 'mdh']):
        return 'TCA Cycle'
    if any(k in rxn_id_lower or k in name_lower for k in ['icl', 'mals', 'glyoxylate']):
        return 'Glyoxylate Shunt'
    if any(k in rxn_id_lower or k in name_lower for k in ['pepck', 'ppdk', 'pc', 'fructose', 'glucose', 'fba', 'fbp', 'pgi', 'tpi']):
        return 'Gluconeogenesis / Glycolysis'
    if any(k in rxn_id_lower or k in name_lower for k in ['pentose', 'transketolase', 'transaldolase', 'ribose', 'g6pdh', 'gnd']):
        return 'PPP'
    if any(k in rxn_id_lower or k in name_lower for k in ['nadl', 'cytbo', 'atps', 'menaquinone', 'ubiquinone', 'q8']):
        return 'ETC / Electron Transport'
    if any(k in name_lower for k in ['serine', 'glycine', 'tyrosine', 'phenylalanine', 'tryptophan', 'leucine', 'isoleucine', 'valine', 'methionine', 'cysteine', 'aspartate', 'asparagine', 'glutamate', 'glutamine', 'lysine', 'arginine', 'histidine', 'proline', 'threonine', 'alanine']):
        return 'Amino Acid Metabolism'
    if any(k in rxn_id_lower or k in name_lower for k in ['nucleotide', 'atp', 'gtp', 'utp', 'ctp']):
        return 'Nucleotide Metabolism'
    if any(k in rxn_id_lower or k in name_lower for k in ['coa', 'nad', 'nadp', 'fad', 'fmn', 'thf', 'folate', 'riboflavin', 'thiamine']):
        return 'Cofactor Metabolism'
    return 'Other'

def prioritize_missing_reactions(missing_reactions):
    """누락된 반응 우선순위 정리"""
    high_pathways = ['Acetate Metabolism', 'TCA Cycle', 'Glyoxylate Shunt', 
                     'Gluconeogenesis / Glycolysis', 'ETC / Electron Transport']
    
    for rxn in missing_reactions:
        if rxn['is_active']:
            if rxn['pathway'] in high_pathways and rxn['max_flux'] > 0.01:
                rxn['priority'] = 'HIGH'
            elif rxn['max_flux'] > 0.1:
                rxn['priority'] = 'HIGH'
            else:
                rxn['priority'] = 'MEDIUM'
        else:
            if rxn['pathway'] in high_pathways:
                rxn['priority'] = 'MEDIUM'
            elif rxn['pathway'] in ['Exchange', 'Transport']:
                rxn['priority'] = 'MEDIUM'
            else:
                rxn['priority'] = 'LOW'
    
    return missing_reactions

def main():
    base_path = Path(__file__).parent.parent
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    new_model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    flux_file = base_path / "Stenotrophomonas" / "fba_flux_gradient_acid.csv"
    
    print("="*70)
    print("레퍼런스 모델 vs 신규 모델 전체 비교 분석")
    print("="*70)
    
    # 모델 로드
    ref_model = load_model(ref_model_path)
    new_model = load_model(new_model_path)
    
    # 실제 사용된 반응 추출
    active_rxn_ids, max_fluxes = get_active_reactions_from_fba(flux_file, threshold=1e-6)
    
    # 전체 비교
    ref_only, active_missing, active_present = compare_all_reactions(
        ref_model, new_model, active_rxn_ids, max_fluxes)
    
    # 누락된 반응 상세 분석
    missing_reactions = analyze_reactions_in_detail(
        ref_model, ref_only, active_rxn_ids, max_fluxes)
    
    # 우선순위 정리
    missing_reactions = prioritize_missing_reactions(missing_reactions)
    
    # DataFrame 생성
    df_missing = pd.DataFrame(missing_reactions)
    
    # 우선순위별 정렬
    priority_order = {'HIGH': 1, 'MEDIUM': 2, 'LOW': 3}
    df_missing['priority_order'] = df_missing['priority'].map(priority_order)
    df_missing = df_missing.sort_values(
        ['is_active', 'priority_order', 'max_flux', 'pathway', 'reaction_id'],
        ascending=[False, True, False, True, True])
    df_missing = df_missing.drop('priority_order', axis=1)
    
    # CSV 저장
    script_dir = Path(__file__).parent
    output_file = script_dir / "comprehensive_missing_reactions.csv"
    df_missing.to_csv(output_file, index=False, encoding='utf-8-sig')
    print(f"\n[OK] 전체 누락 반응 리스트 저장: {output_file}")
    
    # 실제 사용된 반응 중 누락된 것만 별도 저장
    df_active_missing = df_missing[df_missing['is_active'] == True]
    active_output_file = script_dir / "missing_active_reactions_detailed.csv"
    df_active_missing.to_csv(active_output_file, index=False, encoding='utf-8-sig')
    print(f"[OK] 실제 사용된 반응 중 누락: {active_output_file} ({len(df_active_missing)}개)")
    
    # 요약 출력
    print("\n" + "="*70)
    print("요약")
    print("="*70)
    
    print(f"\n레퍼런스 모델에서 실제 사용된 반응: {len(active_rxn_ids)}개")
    print(f"  - 신규 모델에 있음: {len(active_present)}개")
    print(f"  - 신규 모델에 없음: {len(active_missing)}개")
    
    print(f"\n레퍼런스에만 있는 전체 반응: {len(ref_only)}개")
    print(f"  - 실제 사용됨: {len(active_missing)}개")
    print(f"  - 실제 사용 안 됨: {len(ref_only) - len(active_missing)}개")
    
    print("\n우선순위별 분류:")
    priority_summary = df_missing.groupby('priority').size()
    for priority, count in priority_summary.items():
        print(f"  {priority}: {count}개")
    
    print("\n경로별 분류 (상위 10개):")
    pathway_summary = df_missing.groupby('pathway').size().sort_values(ascending=False).head(10)
    for pathway, count in pathway_summary.items():
        print(f"  {pathway}: {count}개")
    
    print("\n[OK] 분석 완료!")

if __name__ == "__main__":
    main()
