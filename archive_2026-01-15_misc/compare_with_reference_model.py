#!/usr/bin/env python
"""
레퍼런스 모델과 신규 모델 비교 스크립트
레퍼런스 모델(Stenotrophomonas/scenarios/YE0p5_clean/model_YE0p5.xml)과 
신규 모델(Stenotrophomonas-causal AI/BaseModel.xml)을 비교하여
차이점을 분석하고 누락된 반응을 리스트업합니다.
"""

import cobra
from pathlib import Path
import pandas as pd
from collections import defaultdict

def load_model(model_path):
    """모델 로드"""
    print(f"모델 로드 중: {model_path}")
    model = cobra.io.read_sbml_model(str(model_path))
    print(f"[OK] 모델 로드 완료: {model.id}")
    print(f"  - 반응 수: {len(model.reactions)}")
    print(f"  - 대사물질 수: {len(model.metabolites)}")
    print(f"  - 유전자 수: {len(model.genes)}")
    return model

def get_reaction_info(reaction):
    """반응 정보 추출"""
    info = {
        'id': reaction.id,
        'name': reaction.name,
        'equation': reaction.reaction,
        'genes': [g.id for g in reaction.genes],
        'lower_bound': reaction.lower_bound,
        'upper_bound': reaction.upper_bound,
        'reversible': reaction.lower_bound < 0
    }
    return info

def compare_reactions(ref_model, new_model):
    """반응 비교"""
    print("\n" + "="*70)
    print("반응 비교")
    print("="*70)
    
    ref_rxns = set(ref_model.reactions.list_attr('id'))
    new_rxns = set(new_model.reactions.list_attr('id'))
    
    common_rxns = ref_rxns & new_rxns
    ref_only_rxns = ref_rxns - new_rxns
    new_only_rxns = new_rxns - ref_rxns
    
    print(f"\n레퍼런스 모델 반응 수: {len(ref_rxns)}")
    print(f"신규 모델 반응 수: {len(new_rxns)}")
    print(f"공통 반응 수: {len(common_rxns)}")
    print(f"레퍼런스에만 있는 반응 수: {len(ref_only_rxns)}")
    print(f"신규 모델에만 있는 반응 수: {len(new_only_rxns)}")
    
    return ref_only_rxns, new_only_rxns, common_rxns

def analyze_missing_reactions(ref_model, new_model, missing_rxn_ids):
    """누락된 반응 상세 분석"""
    print("\n" + "="*70)
    print("레퍼런스에만 있는 반응 상세 분석")
    print("="*70)
    
    missing_reactions = []
    
    for rxn_id in sorted(missing_rxn_ids):
        if rxn_id in ref_model.reactions:
            rxn = ref_model.reactions.get_by_id(rxn_id)
            info = get_reaction_info(rxn)
            
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
    equation_lower = reaction.reaction.lower()
    
    # Exchange reactions
    if rxn_id.startswith('EX_'):
        return 'Exchange'
    
    # Transport reactions
    if 'transport' in name_lower or '_t' in rxn_id_lower or '_pp' in rxn_id_lower or '_p' in rxn_id_lower:
        return 'Transport'
    
    # Amino acid synthesis
    aa_keywords = ['serine', 'glycine', 'tyrosine', 'phenylalanine', 'tryptophan', 
                   'leucine', 'isoleucine', 'valine', 'methionine', 'cysteine',
                   'aspartate', 'asparagine', 'glutamate', 'glutamine', 'lysine',
                   'arginine', 'histidine', 'proline', 'threonine', 'alanine']
    if any(keyword in name_lower or keyword in equation_lower for keyword in aa_keywords):
        return 'Amino Acid Synthesis'
    
    # Nucleotide synthesis
    if any(keyword in rxn_id_lower or keyword in name_lower for keyword in ['nucleotide', 'atp', 'gtp', 'utp', 'ctp', 'adp', 'gdp', 'udp', 'cdp']):
        return 'Nucleotide Synthesis'
    
    # Coenzyme/cofactor synthesis
    if any(keyword in rxn_id_lower or keyword in name_lower for keyword in ['coa', 'nad', 'nadp', 'fad', 'fmn', 'thf', 'folate', 'biotin', 'pantothenate']):
        return 'Cofactor Synthesis'
    
    # Central carbon metabolism
    if any(keyword in rxn_id_lower or keyword in name_lower for keyword in ['glycolysis', 'tca', 'citrate', 'pyruvate', 'glucose', 'pentose', 'ppp']):
        return 'Central Carbon Metabolism'
    
    # Lipid metabolism
    if any(keyword in rxn_id_lower or keyword in name_lower for keyword in ['fatty', 'lipid', 'acyl', 'glycerol']):
        return 'Lipid Metabolism'
    
    # Biomass
    if 'biomass' in rxn_id_lower or 'growth' in rxn_id_lower:
        return 'Biomass'
    
    return 'Other'

def compare_genes(ref_model, new_model):
    """유전자 비교"""
    print("\n" + "="*70)
    print("유전자 비교")
    print("="*70)
    
    ref_genes = set(ref_model.genes.list_attr('id'))
    new_genes = set(new_model.genes.list_attr('id'))
    
    common_genes = ref_genes & new_genes
    ref_only_genes = ref_genes - new_genes
    new_only_genes = new_genes - ref_genes
    
    print(f"\n레퍼런스 모델 유전자 수: {len(ref_genes)}")
    print(f"신규 모델 유전자 수: {len(new_genes)}")
    print(f"공통 유전자 수: {len(common_genes)}")
    print(f"레퍼런스에만 있는 유전자 수: {len(ref_only_genes)}")
    print(f"신규 모델에만 있는 유전자 수: {len(new_only_genes)}")
    
    return ref_only_genes

def get_reactions_by_genes(ref_model, missing_genes):
    """누락된 유전자가 관련된 반응 찾기"""
    gene_to_reactions = defaultdict(list)
    
    for gene_id in missing_genes:
        if gene_id in ref_model.genes:
            gene = ref_model.genes.get_by_id(gene_id)
            for reaction in gene.reactions:
                gene_to_reactions[gene_id].append(reaction.id)
    
    return gene_to_reactions

def prioritize_missing_reactions(missing_reactions):
    """누락된 반응 우선순위 정리"""
    print("\n" + "="*70)
    print("누락된 반응 우선순위 분석")
    print("="*70)
    
    # 경로별 그룹화
    pathway_groups = defaultdict(list)
    for rxn in missing_reactions:
        pathway_groups[rxn['pathway']].append(rxn)
    
    # 우선순위 정의
    priority_order = [
        'Exchange',
        'Transport',
        'Central Carbon Metabolism',
        'Amino Acid Synthesis',
        'Nucleotide Synthesis',
        'Cofactor Synthesis',
        'Lipid Metabolism',
        'Biomass',
        'Other'
    ]
    
    prioritized = []
    for pathway in priority_order:
        if pathway in pathway_groups:
            prioritized.extend(pathway_groups[pathway])
    
    return prioritized

def main():
    # 경로 설정
    base_path = Path(__file__).parent.parent
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    new_model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    
    # 모델 로드
    print("="*70)
    print("레퍼런스 모델과 신규 모델 비교")
    print("="*70)
    
    ref_model = load_model(ref_model_path)
    print()
    new_model = load_model(new_model_path)
    
    # 반응 비교
    ref_only_rxns, new_only_rxns, common_rxns = compare_reactions(ref_model, new_model)
    
    # 누락된 반응 분석
    missing_reactions = analyze_missing_reactions(ref_model, new_model, ref_only_rxns)
    
    # 유전자 비교
    ref_only_genes = compare_genes(ref_model, new_model)
    
    # 누락된 반응 우선순위 정리
    prioritized_reactions = prioritize_missing_reactions(missing_reactions)
    
    # 결과를 DataFrame으로 변환
    df_missing = pd.DataFrame(prioritized_reactions)
    
    # CSV 파일로 저장
    output_file = Path(__file__).parent / "missing_reactions_vs_reference.csv"
    df_missing.to_csv(output_file, index=False, encoding='utf-8-sig')
    print(f"\n[OK] 누락된 반응 리스트 저장: {output_file}")
    
    # 경로별 요약
    print("\n" + "="*70)
    print("경로별 누락된 반응 요약")
    print("="*70)
    
    pathway_summary = df_missing.groupby('pathway').size().sort_values(ascending=False)
    for pathway, count in pathway_summary.items():
        print(f"  {pathway}: {count}개")
    
    # 상위 20개 반응 출력
    print("\n" + "="*70)
    print("우선순위 상위 20개 누락된 반응")
    print("="*70)
    
    for i, rxn in enumerate(prioritized_reactions[:20], 1):
        print(f"\n{i}. {rxn['reaction_id']} ({rxn['pathway']})")
        print(f"   이름: {rxn['name']}")
        print(f"   반응식: {rxn['equation']}")
        if rxn['genes']:
            print(f"   유전자: {rxn['genes']}")
    
    # 전체 요약
    print("\n" + "="*70)
    print("비교 요약")
    print("="*70)
    print(f"\n레퍼런스 모델에만 있는 반응: {len(ref_only_rxns)}개")
    print(f"신규 모델에만 있는 반응: {len(new_only_rxns)}개")
    print(f"레퍼런스 모델에만 있는 유전자: {len(ref_only_genes)}개")
    
    print("\n[OK] 비교 완료!")

if __name__ == "__main__":
    main()
