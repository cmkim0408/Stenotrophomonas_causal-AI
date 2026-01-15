#!/usr/bin/env python
"""
누락된 반응 리스트업 및 검색
진단 결과에서 발견한 누락된 반응들을 체계적으로 정리하고 모델에서 검색
"""

import cobra
import pandas as pd

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    print(f"[OK] 모델 로드: {model.id}")
    return model

def search_reaction_in_model(model, rxn_id, search_patterns=None):
    """모델에서 반응 검색 (정확한 ID, 유사한 ID, 이름으로 검색)"""
    results = {
        'exact_id': None,
        'similar_ids': [],
        'similar_names': [],
        'exists': False
    }
    
    # 1. 정확한 ID로 검색
    try:
        rxn = model.reactions.get_by_id(rxn_id)
        results['exact_id'] = rxn.id
        results['exists'] = True
        return results
    except KeyError:
        pass
    
    # 2. 유사한 ID 패턴 검색
    if search_patterns:
        for pattern in search_patterns:
            for rxn in model.reactions:
                if pattern.lower() in rxn.id.lower():
                    if rxn.id not in results['similar_ids']:
                        results['similar_ids'].append(rxn.id)
    
    # 3. 이름으로 검색
    rxn_id_base = rxn_id.replace('_', '').replace('-', '').lower()
    for rxn in model.reactions:
        rxn_name_lower = rxn.name.lower().replace('_', '').replace('-', '') if rxn.name else ''
        if rxn_id_base in rxn_name_lower or (search_patterns and any(p.lower() in rxn_name_lower for p in search_patterns)):
            if rxn.id not in results['similar_names']:
                results['similar_names'].append(rxn.id)
    
    return results

def list_missing_coa_pathway(model):
    """CoA 생산 경로 누락 반응 리스트"""
    print("\n" + "="*70)
    print("[1] CoA 생산 경로 누락 반응")
    print("="*70)
    
    coa_reactions = [
        {
            'id': 'PNTO',
            'name': 'Pantothenate kinase',
            'equation': 'atp_c + pnto_c --> adp_c + 4ppan_c',
            'pathway': 'Pantothenate -> 4\'-Phosphopantothenate',
            'search_patterns': ['PNTO', 'pantothenate', 'kinase', 'pnto']
        },
        {
            'id': 'PPCS',
            'name': '4-phosphopantothenoylcysteine synthetase',
            'equation': 'atp_c + 4ppan_c + cys__L_c --> adp_c + pi_c + pppcs_c',
            'pathway': '4\'-Phosphopantothenate -> 4\'-Phosphopantothenoylcysteine',
            'search_patterns': ['PPCS', 'phosphopantothenoyl', 'cysteine', 'synthetase']
        },
        {
            'id': 'PPAT',
            'name': 'Phosphopantetheine adenylyltransferase',
            'equation': 'atp_c + pppan_c --> ppi_c + dcoa_c',
            'pathway': '4\'-Phosphopantetheine -> Dephospho-CoA',
            'search_patterns': ['PPAT', 'phosphopantetheine', 'adenylyltransferase', 'dephospho', 'coa']
        }
    ]
    
    missing_list = []
    
    for rxn_info in coa_reactions:
        print(f"\n반응 ID: {rxn_info['id']}")
        print(f"  이름: {rxn_info['name']}")
        print(f"  경로: {rxn_info['pathway']}")
        print(f"  반응식: {rxn_info['equation']}")
        
        # 모델에서 검색
        results = search_reaction_in_model(model, rxn_info['id'], rxn_info['search_patterns'])
        
        if results['exists']:
            print(f"  상태: [EXISTS] {results['exact_id']}")
            rxn = model.reactions.get_by_id(results['exact_id'])
            print(f"    실제 반응식: {rxn.reaction}")
        else:
            print(f"  상태: [MISSING] 모델에 없음")
            if results['similar_ids']:
                print(f"    유사한 ID: {', '.join(results['similar_ids'][:5])}")
            if results['similar_names']:
                print(f"    유사한 이름: {', '.join(results['similar_names'][:5])}")
            
            missing_list.append({
                'Category': 'CoA Pathway',
                'Reaction_ID': rxn_info['id'],
                'Name': rxn_info['name'],
                'Equation': rxn_info['equation'],
                'Pathway': rxn_info['pathway'],
                'Status': 'MISSING',
                'Similar_IDs': ', '.join(results['similar_ids'][:3]),
                'Similar_Names': ', '.join(results['similar_names'][:3])
            })
    
    return missing_list

def list_missing_amino_acid_pathway(model):
    """아미노산 생합성 경로 누락 반응 리스트"""
    print("\n" + "="*70)
    print("[2] 아미노산 생합성 경로 누락 반응")
    print("="*70)
    
    aa_reactions = [
        # Serine pathway
        {
            'id': 'SER',
            'name': 'Serine synthase / Phosphoserine phosphatase',
            'equation': 'pser__L_c + h2o_c --> ser__L_c + pi_c',
            'pathway': 'Serine pathway (Phosphoserine -> Serine)',
            'search_patterns': ['SER', 'serine', 'phosphatase', 'phosphoserine']
        },
        {
            'id': 'SERD',
            'name': 'Serine dehydratase',
            'equation': 'ser__L_c --> nh3_c + pyr_c + h2o_c',
            'pathway': 'Serine pathway (Serine -> Pyruvate)',
            'search_patterns': ['SERD', 'serine', 'dehydratase', 'dehydrase']
        },
        # Glycine pathway
        {
            'id': 'GCY',
            'name': 'Glycine synthase / Glycine cleavage system',
            'equation': 'co2_c + nh4_c + 5fthf_c + nad_c --> gly_c + nadh_c + mlthf_c',
            'pathway': 'Glycine pathway (CO2 + NH4 -> Glycine)',
            'search_patterns': ['GCY', 'glycine', 'synthase', 'cleavage', 'gly']
        },
        {
            'id': 'SHMT',
            'name': 'Serine hydroxymethyltransferase',
            'equation': 'ser__L_c + thf_c --> gly_c + mlthf_c + h2o_c',
            'pathway': 'Glycine pathway (Serine -> Glycine)',
            'search_patterns': ['SHMT', 'serine', 'hydroxymethyltransferase', 'glycine']
        },
        # Aromatic AA
        {
            'id': 'DAHPS',
            'name': '3-deoxy-7-phosphoheptulonate synthase',
            'equation': 'e4p_c + pep_c --> dah7p_c + pi_c + h2o_c',
            'pathway': 'Aromatic AA (Erythrose-4-P + PEP -> DAHP)',
            'search_patterns': ['DAHPS', 'deoxy', 'phosphoheptulonate', 'dahp', 'aro']
        },
        {
            'id': 'TRPS',
            'name': 'Tryptophan synthase',
            'equation': 'indole_c + ser__L_c --> trp__L_c + h2o_c',
            'pathway': 'Aromatic AA (Indole + Serine -> Tryptophan)',
            'search_patterns': ['TRPS', 'tryptophan', 'synthase', 'indole']
        },
        {
            'id': 'TYRS',
            'name': 'Tyrosine synthase',
            'equation': 'prephenate_c --> tyr__L_c + co2_c + h2o_c',
            'pathway': 'Aromatic AA (Prephenate -> Tyrosine)',
            'search_patterns': ['TYRS', 'tyrosine', 'synthase', 'prephenate']
        },
        {
            'id': 'PHES',
            'name': 'Phenylalanine synthase',
            'equation': 'prephenate_c --> phe__L_c + co2_c + h2o_c',
            'pathway': 'Aromatic AA (Prephenate -> Phenylalanine)',
            'search_patterns': ['PHES', 'phenylalanine', 'synthase', 'prephenate']
        },
        # Aspartate family
        {
            'id': 'HOM',
            'name': 'Homoserine dehydrogenase',
            'equation': 'aspartate_semialdehyde_c + nadph_c + h_c --> hom__L_c + nadp_c + h2o_c',
            'pathway': 'Aspartate family (Aspartate semialdehyde -> Homoserine)',
            'search_patterns': ['HOM', 'homoserine', 'dehydrogenase', 'hom']
        },
        {
            'id': 'LYSS',
            'name': 'Lysine synthase / Diaminopimelate pathway',
            'equation': 'meso_dap_c --> lys__L_c + succ_c',
            'pathway': 'Aspartate family (meso-Diaminopimelate -> Lysine)',
            'search_patterns': ['LYSS', 'lysine', 'synthase', 'diaminopimelate', 'dap']
        }
    ]
    
    missing_list = []
    
    for rxn_info in aa_reactions:
        print(f"\n반응 ID: {rxn_info['id']}")
        print(f"  이름: {rxn_info['name']}")
        print(f"  경로: {rxn_info['pathway']}")
        print(f"  반응식: {rxn_info['equation']}")
        
        # 모델에서 검색
        results = search_reaction_in_model(model, rxn_info['id'], rxn_info['search_patterns'])
        
        if results['exists']:
            print(f"  상태: [EXISTS] {results['exact_id']}")
            rxn = model.reactions.get_by_id(results['exact_id'])
            print(f"    실제 반응식: {rxn.reaction}")
        else:
            print(f"  상태: [MISSING] 모델에 없음")
            if results['similar_ids']:
                print(f"    유사한 ID: {', '.join(results['similar_ids'][:5])}")
            if results['similar_names']:
                print(f"    유사한 이름: {', '.join(results['similar_names'][:5])}")
            
            missing_list.append({
                'Category': 'Amino Acid Pathway',
                'Reaction_ID': rxn_info['id'],
                'Name': rxn_info['name'],
                'Equation': rxn_info['equation'],
                'Pathway': rxn_info['pathway'],
                'Status': 'MISSING',
                'Similar_IDs': ', '.join(results['similar_ids'][:3]),
                'Similar_Names': ', '.join(results['similar_names'][:3])
            })
    
    return missing_list

def list_missing_nucleotide_pathway(model):
    """뉴클레오티드 생합성 경로 누락 반응 리스트"""
    print("\n" + "="*70)
    print("[3] 뉴클레오티드 생합성 경로 누락 반응")
    print("="*70)
    
    # ATP 생산 경로 문제 분석
    print("\nATP 생산 경로 분석:")
    try:
        atp_c = model.metabolites.get_by_id('atp_c')
        atp_producing = [r for r in atp_c.reactions if atp_c in r.products]
        
        print(f"  ATP 생산 반응: {len(atp_producing)}개")
        print("\n  주요 ATP 생산 반응:")
        for rxn in atp_producing[:5]:
            print(f"    - {rxn.id}: {rxn.name}")
            print(f"      {rxn.reaction}")
        
        if len(atp_producing) < 10:
            print("\n  [WARNING] ATP 생산 반응이 적음")
            print("    - NDPK (Nucleoside diphosphate kinase) 경로 확인 필요")
            print("    - Substrate level phosphorylation 경로 확인 필요")
            print("    - Oxidative phosphorylation 경로 확인 필요")
    except KeyError:
        print("  [ERROR] ATP metabolite 없음")
    
    # NDPK 관련 반응 확인
    print("\nNucleoside Diphosphate Kinase (NDPK) 경로:")
    ndpk_patterns = ['NDPK', 'NDK', 'nucleoside', 'diphosphate', 'kinase']
    
    ndpk_found = []
    for pattern in ndpk_patterns:
        for rxn in model.reactions:
            if pattern.lower() in rxn.id.lower() or (rxn.name and pattern.lower() in rxn.name.lower()):
                if 'diphosphate' in rxn.name.lower() and 'kinase' in rxn.name.lower():
                    if rxn not in ndpk_found:
                        ndpk_found.append(rxn)
    
    if ndpk_found:
        print(f"  발견된 NDPK 관련 반응: {len(ndpk_found)}개")
        for rxn in ndpk_found[:5]:
            print(f"    - {rxn.id}: {rxn.name}")
            print(f"      {rxn.reaction}")
    else:
        print("  [WARNING] NDPK 관련 반응 없음")
        print("    - NDPK1 (ATP:GDP) 필요할 수 있음")
        print("    - NDPK2 (ATP:UDP) 필요할 수 있음")
        print("    - NDPK3 (ATP:CDP) 필요할 수 있음")
    
    # PRPP 생산 경로 확인
    print("\nPRPP 생산 경로:")
    try:
        prpp_c = model.metabolites.get_by_id('prpp_c')
        prpp_producing = [r for r in prpp_c.reactions if prpp_c in r.products]
        
        if prpp_producing:
            print(f"  PRPP 생산 반응: {len(prpp_producing)}개")
            for rxn in prpp_producing[:3]:
                print(f"    - {rxn.id}: {rxn.name}")
                print(f"      {rxn.reaction}")
        else:
            print("  [WARNING] PRPP 생산 반응 없음")
            print("    - PRPPS (Phosphoribosyl pyrophosphate synthetase) 필요")
    except KeyError:
        print("  [WARNING] PRPP metabolite 없음")
        print("    - PRPPS 반응 확인 필요")
    
    missing_list = []
    
    # 주요 누락 가능 반응들
    nucleotide_reactions = [
        {
            'id': 'PRPPS',
            'name': 'Phosphoribosyl pyrophosphate synthetase',
            'equation': 'atp_c + r5p_c --> prpp_c + amp_c',
            'pathway': 'Purine/Pyrimidine synthesis (R5P -> PRPP)',
            'search_patterns': ['PRPPS', 'phosphoribosyl', 'pyrophosphate', 'synthetase', 'prpp']
        }
    ]
    
    for rxn_info in nucleotide_reactions:
        print(f"\n반응 ID: {rxn_info['id']}")
        print(f"  이름: {rxn_info['name']}")
        print(f"  경로: {rxn_info['pathway']}")
        print(f"  반응식: {rxn_info['equation']}")
        
        results = search_reaction_in_model(model, rxn_info['id'], rxn_info['search_patterns'])
        
        if results['exists']:
            print(f"  상태: [EXISTS] {results['exact_id']}")
        else:
            print(f"  상태: [MISSING] 모델에 없음")
            if results['similar_ids'] or results['similar_names']:
                print(f"    유사한 반응: {', '.join(results['similar_ids'][:3] + results['similar_names'][:3])}")
            
            missing_list.append({
                'Category': 'Nucleotide Pathway',
                'Reaction_ID': rxn_info['id'],
                'Name': rxn_info['name'],
                'Equation': rxn_info['equation'],
                'Pathway': rxn_info['pathway'],
                'Status': 'MISSING',
                'Similar_IDs': ', '.join(results['similar_ids'][:3]),
                'Similar_Names': ', '.join(results['similar_names'][:3])
            })
    
    return missing_list

def check_existing_similar_reactions(model, reaction_list):
    """유사한 반응이 실제로 존재하는지 상세 확인"""
    print("\n" + "="*70)
    print("[4] 유사한 반응 상세 확인")
    print("="*70)
    
    for rxn_info in reaction_list:
        if not rxn_info.get('Similar_IDs') and not rxn_info.get('Similar_Names'):
            continue
        
        print(f"\n반응 ID: {rxn_info['Reaction_ID']}")
        print(f"  찾은 유사 반응:")
        
        # Similar_IDs 확인
        if rxn_info.get('Similar_IDs'):
            for sim_id in rxn_info['Similar_IDs'].split(', '):
                if sim_id:
                    try:
                        sim_rxn = model.reactions.get_by_id(sim_id.strip())
                        print(f"    - {sim_id}: {sim_rxn.name}")
                        print(f"      {sim_rxn.reaction}")
                    except KeyError:
                        pass
        
        # Similar_Names 확인
        if rxn_info.get('Similar_Names'):
            for sim_id in rxn_info['Similar_Names'].split(', '):
                if sim_id:
                    try:
                        sim_rxn = model.reactions.get_by_id(sim_id.strip())
                        print(f"    - {sim_id}: {sim_rxn.name}")
                        print(f"      {sim_rxn.reaction}")
                    except KeyError:
                        pass

def main():
    print("="*70)
    print("누락된 반응 리스트업 및 검색")
    print("="*70)
    
    # 모델 로드
    model = load_model("BaseModel.xml")
    
    print(f"\n총 반응 수: {len(model.reactions)}")
    print(f"총 대사물질 수: {len(model.metabolites)}")
    
    # 1. CoA 경로 누락 반응
    coa_missing = list_missing_coa_pathway(model)
    
    # 2. 아미노산 경로 누락 반응
    aa_missing = list_missing_amino_acid_pathway(model)
    
    # 3. 뉴클레오티드 경로 누락 반응
    nucleotide_missing = list_missing_nucleotide_pathway(model)
    
    # 모든 누락 반응 통합
    all_missing = coa_missing + aa_missing + nucleotide_missing
    
    # 유사한 반응 상세 확인
    check_existing_similar_reactions(model, all_missing)
    
    # 결과 저장
    if all_missing:
        df_missing = pd.DataFrame(all_missing)
        df_missing.to_csv('missing_reactions_list.csv', index=False)
        print(f"\n[OK] 누락 반응 리스트 저장: missing_reactions_list.csv")
        
        # 카테고리별 요약
        print("\n" + "="*70)
        print("누락 반응 요약")
        print("="*70)
        
        by_category = df_missing.groupby('Category').size()
        print("\n카테고리별 누락 반응 수:")
        for category, count in by_category.items():
            print(f"  {category}: {count}개")
        
        print(f"\n총 누락 반응 수: {len(all_missing)}개")
        
        print("\n상세 리스트:")
        for _, row in df_missing.iterrows():
            print(f"\n  [{row['Category']}] {row['Reaction_ID']}: {row['Name']}")
            print(f"    경로: {row['Pathway']}")
            print(f"    반응식: {row['Equation']}")
    
    print("\n" + "="*70)
    print("리스트업 완료")
    print("="*70)

if __name__ == "__main__":
    main()
