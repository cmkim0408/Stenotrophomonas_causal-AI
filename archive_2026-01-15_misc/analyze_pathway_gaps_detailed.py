#!/usr/bin/env python
"""
경로 Gap 상세 분석
차단된 대사물질의 생산 경로를 추적하여 누락된 반응 식별
"""

import cobra
import pandas as pd
from collections import deque

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    return model

def find_biomass_reaction(model):
    for name in ['Growth', 'BIOMASS']:
        try:
            return model.reactions.get_by_id(name)
        except KeyError:
            continue
    return None

def setup_test_environment(model, substrate_type='glucose'):
    """테스트 환경 설정"""
    # 모든 exchange 초기화
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    # 기질 설정
    if substrate_type == 'glucose':
        try:
            ex_sub = model.reactions.get_by_id('EX_glc__D_e')
            ex_sub.lower_bound = -100
            ex_sub.upper_bound = 1000
        except KeyError:
            pass
    elif substrate_type == 'acetate':
        try:
            ex_sub = model.reactions.get_by_id('EX_ac_e')
            ex_sub.lower_bound = -100
            ex_sub.upper_bound = 1000
        except KeyError:
            pass
    
    # 필수 영양소
    essentials = ['EX_nh4_e', 'EX_h2o_e', 'EX_h_e', 'EX_pi_e', 'EX_so4_e',
                  'EX_k_e', 'EX_na1_e', 'EX_mg2_e', 'EX_ca2_e', 'EX_fe2_e',
                  'EX_mn2_e', 'EX_zn2_e', 'EX_co2_e', 'EX_o2_e']
    
    for ex_id in essentials:
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            if ex_id in ['EX_co2_e', 'EX_o2_e']:
                ex_rxn.lower_bound = -1000
                ex_rxn.upper_bound = 1000
            else:
                ex_rxn.lower_bound = -1000
        except KeyError:
            pass
    
    # Transport 수정
    try:
        if substrate_type == 'glucose':
            try:
                model.reactions.get_by_id('GLCt')
            except KeyError:
                glc__D_e = model.metabolites.get_by_id('glc__D_e')
                glc__D_c = model.metabolites.get_by_id('glc__D_c')
                glc_transport = cobra.Reaction('GLCt')
                glc_transport.lower_bound = -1000
                glc_transport.upper_bound = 1000
                glc_transport.add_metabolites({glc__D_e: -1, glc__D_c: 1})
                model.add_reactions([glc_transport])
        
        hex1 = model.reactions.get_by_id('HEX1')
        hex1.lower_bound = -1000
    except:
        pass

def test_metabolite_production(model, met_id):
    """대사물질 생산 가능 여부 테스트"""
    try:
        met = model.metabolites.get_by_id(met_id)
        test_rxn = cobra.Reaction(f'TEST_{met_id}')
        test_rxn.add_metabolites({met: 1})
        test_rxn.lower_bound = 0
        test_rxn.upper_bound = 1000
        model.add_reactions([test_rxn])
        model.objective = test_rxn.id
        
        solution = model.optimize()
        model.remove_reactions([test_rxn])
        
        return solution.status == 'optimal' and solution.objective_value > 1e-6
    except:
        return False

def trace_production_pathway(model, target_met_id, depth=3, visited=None):
    """대사물질 생산 경로 추적"""
    if visited is None:
        visited = set()
    
    if target_met_id in visited or depth <= 0:
        return []
    
    visited.add(target_met_id)
    
    try:
        met = model.metabolites.get_by_id(target_met_id)
    except KeyError:
        return []
    
    pathways = []
    
    # 이 대사물질을 생성하는 반응 찾기
    for rxn in met.reactions:
        if met in rxn.products:
            # 반응물 확인
            reactants = [m for m in rxn.reactants if m != met]
            
            pathway = {
                'Reaction_ID': rxn.id,
                'Reaction_Name': rxn.name,
                'Equation': rxn.reaction,
                'Reactants': [m.id for m in reactants],
                'GPR': getattr(rxn, 'gene_reaction_rule', 'N/A'),
                'Depth': 3 - depth
            }
            
            # 반응물 중 생산 불가능한 것 확인
            blocked_reactants = []
            for reactant in reactants:
                if not test_metabolite_production(model, reactant.id):
                    blocked_reactants.append(reactant.id)
                    # 재귀적으로 추적
                    sub_pathways = trace_production_pathway(model, reactant.id, depth-1, visited.copy())
                    pathway['Sub_pathways'] = sub_pathways
            
            pathway['Blocked_Reactants'] = blocked_reactants
            
            if blocked_reactants:
                pathways.append(pathway)
    
    return pathways

def analyze_specific_pathways(model):
    """주요 경로 상세 분석"""
    print("="*70)
    print("주요 경로 상세 분석")
    print("="*70)
    
    setup_test_environment(model, 'glucose')
    
    key_pathways = {
        'atp_c': {
            'name': 'ATP',
            'priority_pathways': [
                'ATPS4rpp',  # ATP synthase (ETC)
                'PYK',       # Pyruvate kinase
                'PGK',       # Phosphoglycerate kinase
            ]
        },
        'nad_c': {
            'name': 'NAD+',
            'priority_pathways': [
                'NAD synthetase',
                'NAD salvage',
            ]
        },
        'coa_c': {
            'name': 'CoA',
            'priority_pathways': [
                'PNTK',      # Pantothenate kinase
                'PPCDC',     # Phosphopantothenoylcysteine decarboxylase
                'DPCOAK',    # Dephospho-CoA kinase
            ]
        },
        'pep_c': {
            'name': 'PEP',
            'priority_pathways': [
                'PEPCK',     # PEP carboxykinase
                'PPS',       # PEP synthase
            ]
        },
    }
    
    pathway_analysis = []
    
    for met_id, pathway_info in key_pathways.items():
        print(f"\n[{pathway_info['name']} ({met_id}) 생산 경로 분석]")
        print("-" * 70)
        
        # 생산 가능 여부 확인
        can_produce = test_metabolite_production(model, met_id)
        
        if not can_produce:
            print(f"  [BLOCKED] {pathway_info['name']} 생산 불가능")
            
            # 우선순위 경로 확인
            print("\n  우선순위 경로 확인:")
            for pathway_name in pathway_info['priority_pathways']:
                try:
                    rxn = model.reactions.get_by_id(pathway_name)
                    print(f"    [OK] {pathway_name}: {rxn.name}")
                    print(f"      반응식: {rxn.reaction}")
                    if hasattr(rxn, 'gene_reaction_rule'):
                        print(f"      GPR: {rxn.gene_reaction_rule}")
                    
                    # 반응물 생산 가능 여부 확인
                    reactants_blocked = []
                    for reactant in rxn.reactants:
                        reactant_id = reactant.id
                        if not test_metabolite_production(model, reactant_id):
                            reactants_blocked.append(reactant_id)
                    
                    if reactants_blocked:
                        print(f"      [GAP] 차단된 반응물: {', '.join(reactants_blocked)}")
                        pathway_analysis.append({
                            'Target_Metabolite': pathway_info['name'],
                            'Target_ID': met_id,
                            'Reaction_ID': pathway_name,
                            'Reaction_Name': rxn.name,
                            'Equation': rxn.reaction,
                            'GPR': getattr(rxn, 'gene_reaction_rule', 'N/A'),
                            'Blocked_Reactants': ', '.join(reactants_blocked),
                            'Pathway_Type': 'Priority'
                        })
                    else:
                        print(f"      [OK] 모든 반응물 생산 가능")
                except KeyError:
                    print(f"    [MISSING] {pathway_name}: 반응 없음")
                    pathway_analysis.append({
                        'Target_Metabolite': pathway_info['name'],
                        'Target_ID': met_id,
                        'Reaction_ID': pathway_name,
                        'Reaction_Name': 'NOT FOUND',
                        'Equation': 'N/A',
                        'GPR': 'N/A',
                        'Blocked_Reactants': 'Reaction missing',
                        'Pathway_Type': 'Missing'
                    })
            
            # 경로 추적
            print("\n  생산 경로 추적:")
            pathways = trace_production_pathway(model, met_id, depth=3)
            if pathways:
                for pathway in pathways[:3]:  # 처음 3개만
                    print(f"    [{pathway['Reaction_ID']}]: {pathway['Reaction_Name']}")
                    if pathway['Blocked_Reactants']:
                        print(f"      차단된 반응물: {', '.join(pathway['Blocked_Reactants'])}")
        else:
            print(f"  [OK] {pathway_info['name']} 생산 가능")
    
    return pathway_analysis

def find_ec_numbers_for_reactions(reaction_ids):
    """반응 ID에 대한 EC 번호 찾기 (KEGG 정보 기반)"""
    print("\n" + "="*70)
    print("EC 번호 및 효소 정보 찾기")
    print("="*70)
    
    # 일반적인 반응-EC 매핑 (KEGG 기반)
    reaction_ec_mapping = {
        'ATPS4rpp': ['EC:7.1.2.2'],  # ATP synthase
        'PYK': ['EC:2.7.1.40'],      # Pyruvate kinase
        'PGK': ['EC:2.7.2.3'],       # Phosphoglycerate kinase
        'PEPCK': ['EC:4.1.1.32', 'EC:4.1.1.49'],  # PEP carboxykinase
        'PPS': ['EC:2.7.9.2'],       # PEP synthase
        'PNTK': ['EC:2.7.1.33'],     # Pantothenate kinase
        'PPCDC': ['EC:4.1.1.36'],    # Phosphopantothenoylcysteine decarboxylase
        'DPCOAK': ['EC:2.7.1.24'],   # Dephospho-CoA kinase
    }
    
    enzyme_info = []
    
    for rxn_id in reaction_ids:
        print(f"\n[{rxn_id}]")
        print("-" * 70)
        
        if rxn_id in reaction_ec_mapping:
            ec_numbers = reaction_ec_mapping[rxn_id]
            print(f"  EC 번호: {', '.join(ec_numbers)}")
            
            for ec in ec_numbers:
                enzyme_info.append({
                    'Reaction_ID': rxn_id,
                    'EC_Number': ec,
                    'Source': 'Manual mapping'
                })
        else:
            print(f"  [WARNING] EC 번호 매핑 없음")
            print(f"  → KEGG 또는 UniProt에서 수동으로 확인 필요")
            enzyme_info.append({
                'Reaction_ID': rxn_id,
                'EC_Number': 'Unknown',
                'Source': 'Not found'
            })
    
    return enzyme_info

def main():
    print("="*70)
    print("경로 Gap 상세 분석 및 효소 찾기")
    print("="*70)
    
    model = load_model("BaseModel.xml")
    
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass: {biomass_rxn.id}")
    
    # 주요 경로 분석
    pathway_analysis = analyze_specific_pathways(model)
    
    if pathway_analysis:
        # 결과 저장
        df_pathways = pd.DataFrame(pathway_analysis)
        df_pathways.to_csv('pathway_gaps_detailed.csv', index=False, encoding='utf-8-sig')
        print(f"\n[OK] 경로 Gap 상세 분석 저장: pathway_gaps_detailed.csv")
        
        # 누락된 반응 ID 추출
        missing_reaction_ids = [p['Reaction_ID'] for p in pathway_analysis if p['Pathway_Type'] == 'Missing']
        gap_reaction_ids = [p['Reaction_ID'] for p in pathway_analysis if p['Blocked_Reactants'] and p['Pathway_Type'] == 'Priority']
        
        all_reaction_ids = list(set(missing_reaction_ids + gap_reaction_ids))
        
        if all_reaction_ids:
            # EC 번호 찾기
            enzyme_info = find_ec_numbers_for_reactions(all_reaction_ids)
            
            if enzyme_info:
                df_enzymes = pd.DataFrame(enzyme_info)
                df_enzymes.to_csv('reaction_ec_mapping.csv', index=False, encoding='utf-8-sig')
                print(f"[OK] 반응-EC 매핑 저장: reaction_ec_mapping.csv")
    
    print("\n" + "="*70)
    print("다음 단계")
    print("="*70)
    print("1. pathway_gaps_detailed.csv 확인")
    print("2. reaction_ec_mapping.csv 확인")
    print("3. 누락된 반응의 EC 번호를 기반으로 KEGG/UniProt에서 효소 찾기")
    print("4. Stenotrophomonas maltophilia K279a 게놈에서 해당 유전자 확인")
    print("="*70)

if __name__ == "__main__":
    main()
