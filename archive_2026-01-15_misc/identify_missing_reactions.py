#!/usr/bin/env python
"""
누락된 반응 식별
포도당/아세트산 생장 실패 원인을 정확히 파악하여 필요한 반응 식별
"""

import cobra
import pandas as pd

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

def analyze_blocked_reactions_with_substrate(model, biomass_rxn, substrate_type='glucose'):
    """특정 기질에서 차단된 반응 분석"""
    print("="*70)
    print(f"{substrate_type} 기반 차단 반응 분석")
    print("="*70)
    
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
            print(f"[ERROR] EX_glc__D_e 없음")
            return None
    elif substrate_type == 'acetate':
        try:
            ex_sub = model.reactions.get_by_id('EX_ac_e')
            ex_sub.lower_bound = -100
            ex_sub.upper_bound = 1000
        except KeyError:
            print(f"[ERROR] EX_ac_e 없음")
            return None
    
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
    
    # Transport 수정 (이미 적용되었다고 가정)
    try:
        if substrate_type == 'glucose':
            try:
                glc_transport = model.reactions.get_by_id('GLCt')
            except KeyError:
                glc__D_e = model.metabolites.get_by_id('glc__D_e')
                glc__D_c = model.metabolites.get_by_id('glc__D_c')
                glc_transport = cobra.Reaction('GLCt')
                glc_transport.lower_bound = -1000
                glc_transport.upper_bound = 1000
                glc_transport.add_metabolites({glc__D_e: -1, glc__D_c: 1})
                model.add_reactions([glc_transport])
    except:
        pass
    
    try:
        hex1 = model.reactions.get_by_id('HEX1')
        hex1.lower_bound = -1000
    except:
        pass
    
    # Biomass 최적화 시도
    model.objective = biomass_rxn.id
    solution = model.optimize()
    
    if solution.status == 'infeasible':
        print(f"\n[상태] {substrate_type}에서 infeasible")
        print("\n차단된 대사물질 식별 중...")
        
        # 주요 대사물질 생산 가능 여부 테스트
        key_metabolites = {
            'atp_c': 'ATP',
            'nad_c': 'NAD+',
            'nadh_c': 'NADH',
            'coa_c': 'CoA',
            'pep_c': 'PEP',
            'g6p_c': 'Glucose-6-phosphate',
            'pyr_c': 'Pyruvate',
            'accoa_c': 'Acetyl-CoA',
            'h_p': 'Proton motive force',
            'imp_c': 'IMP (Purine)',
            'ump_c': 'UMP (Pyrimidine)',
            'ala__L_c': 'Alanine',
            'gly_c': 'Glycine',
            'ser__L_c': 'Serine',
        }
        
        blocked_metabolites = []
        
        for met_id, met_name in key_metabolites.items():
            try:
                met = model.metabolites.get_by_id(met_id)
                
                test_rxn = cobra.Reaction(f'TEST_{met_id}')
                test_rxn.add_metabolites({met: 1})
                test_rxn.lower_bound = 0
                test_rxn.upper_bound = 1000
                model.add_reactions([test_rxn])
                model.objective = test_rxn.id
                
                test_solution = model.optimize()
                model.remove_reactions([test_rxn])
                
                if test_solution.status != 'optimal' or test_solution.objective_value < 1e-6:
                    blocked_metabolites.append({
                        'Metabolite_ID': met_id,
                        'Name': met_name,
                        'Status': test_solution.status,
                        'Max_Flux': test_solution.objective_value if test_solution.status == 'optimal' else 0
                    })
                    print(f"  [BLOCKED] {met_name} ({met_id}): {test_solution.status}")
                else:
                    print(f"  [OK] {met_name} ({met_id}): 생산 가능")
            except KeyError:
                print(f"  [MISSING] {met_name} ({met_id}): metabolite 없음")
        
        return blocked_metabolites
    
    elif solution.status == 'optimal':
        biomass_flux = solution.fluxes.get(biomass_rxn.id, 0)
        if biomass_flux > 1e-6:
            print(f"\n[OK] {substrate_type}에서 생장 가능!")
            return []
        else:
            print(f"\n[FAIL] {substrate_type}에서 flux = 0")
            return []
    
    return None

def identify_missing_reactions_for_metabolite(model, met_id):
    """특정 대사물질 생산을 위해 필요한 반응 식별"""
    print(f"\n[{met_id} 생산을 위한 필요한 반응 식별]")
    print("-" * 70)
    
    try:
        met = model.metabolites.get_by_id(met_id)
        
        # 이 대사물질을 생성하는 반응 찾기
        producing_reactions = []
        
        for rxn in met.reactions:
            if met in rxn.products:
                # 반응이 활성화될 수 있는지 확인
                # (간단한 방법: 반응의 반응물이 생산 가능한지)
                
                reactants = [m for m in rxn.reactants if m != met]
                
                producing_reactions.append({
                    'Reaction_ID': rxn.id,
                    'Reaction_Name': rxn.name,
                    'Equation': rxn.reaction,
                    'Reactants': [m.id for m in reactants],
                    'GPR': rxn.gene_reaction_rule if hasattr(rxn, 'gene_reaction_rule') else 'N/A'
                })
        
        if producing_reactions:
            print(f"  발견된 생산 반응: {len(producing_reactions)}개")
            for rxn_info in producing_reactions[:5]:  # 처음 5개만
                print(f"\n    [{rxn_info['Reaction_ID']}]: {rxn_info['Reaction_Name']}")
                print(f"      반응식: {rxn_info['Equation']}")
                print(f"      반응물: {', '.join(rxn_info['Reactants'])}")
                if rxn_info['GPR'] != 'N/A':
                    print(f"      GPR: {rxn_info['GPR']}")
        else:
            print(f"  [WARNING] {met_id}을 생산하는 반응 없음")
        
        return producing_reactions
    
    except KeyError:
        print(f"  [ERROR] {met_id} metabolite 없음")
        return []

def find_pathway_gaps(model, blocked_metabolites):
    """차단된 대사물질들의 생산 경로 gap 찾기"""
    print("\n" + "="*70)
    print("경로 Gap 분석")
    print("="*70)
    
    missing_reactions_summary = []
    
    for met_info in blocked_metabolites:
        met_id = met_info['Metabolite_ID']
        met_name = met_info['Name']
        
        print(f"\n[{met_name} ({met_id})]")
        print("-" * 70)
        
        # 생산 반응 찾기
        producing_rxns = identify_missing_reactions_for_metabolite(model, met_id)
        
        # 반응물 생산 가능 여부 확인
        for rxn_info in producing_rxns:
            reactants_blocked = []
            
            for reactant_id in rxn_info['Reactants']:
                try:
                    reactant = model.metabolites.get_by_id(reactant_id)
                    
                    # 반응물 생산 테스트
                    test_rxn = cobra.Reaction(f'TEST_{reactant_id}')
                    test_rxn.add_metabolites({reactant: 1})
                    test_rxn.lower_bound = 0
                    test_rxn.upper_bound = 1000
                    model.add_reactions([test_rxn])
                    model.objective = test_rxn.id
                    
                    test_solution = model.optimize()
                    model.remove_reactions([test_rxn])
                    
                    if test_solution.status != 'optimal' or test_solution.objective_value < 1e-6:
                        reactants_blocked.append(reactant_id)
                except:
                    pass
            
            if reactants_blocked:
                missing_reactions_summary.append({
                    'Target_Metabolite': met_name,
                    'Target_ID': met_id,
                    'Reaction_ID': rxn_info['Reaction_ID'],
                    'Reaction_Name': rxn_info['Reaction_Name'],
                    'Equation': rxn_info['Equation'],
                    'Blocked_Reactants': ', '.join(reactants_blocked),
                    'GPR': rxn_info['GPR']
                })
                
                print(f"\n    [GAP 발견] {rxn_info['Reaction_ID']}")
                print(f"      차단된 반응물: {', '.join(reactants_blocked)}")
    
    return missing_reactions_summary

def main():
    print("="*70)
    print("누락된 반응 식별")
    print("="*70)
    
    model = load_model("BaseModel.xml")
    
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass: {biomass_rxn.id}")
    
    # 포도당에서 차단된 대사물질 식별
    print("\n" + "="*70)
    glucose_blocked = analyze_blocked_reactions_with_substrate(model, biomass_rxn, 'glucose')
    
    # 아세트산에서 차단된 대사물질 식별
    print("\n" + "="*70)
    acetate_blocked = analyze_blocked_reactions_with_substrate(model, biomass_rxn, 'acetate')
    
    # 공통 차단 대사물질 식별
    if glucose_blocked and acetate_blocked:
        glucose_ids = {m['Metabolite_ID'] for m in glucose_blocked}
        acetate_ids = {m['Metabolite_ID'] for m in acetate_blocked}
        common_blocked = glucose_ids & acetate_ids
        
        print("\n" + "="*70)
        print("공통 차단 대사물질")
        print("="*70)
        for met_id in common_blocked:
            glc_info = next(m for m in glucose_blocked if m['Metabolite_ID'] == met_id)
            print(f"  - {glc_info['Name']} ({met_id})")
        
        # 공통 차단 대사물질에 대한 gap 분석
        common_blocked_list = [m for m in glucose_blocked if m['Metabolite_ID'] in common_blocked]
        missing_reactions = find_pathway_gaps(model, common_blocked_list)
        
        # 결과 저장
        if missing_reactions:
            df_missing = pd.DataFrame(missing_reactions)
            df_missing.to_csv('missing_reactions_identified.csv', index=False, encoding='utf-8-sig')
            print(f"\n[OK] 누락된 반응 목록 저장: missing_reactions_identified.csv")
    
    # 개별 결과 저장
    if glucose_blocked:
        df_glucose = pd.DataFrame(glucose_blocked)
        df_glucose.to_csv('glucose_blocked_metabolites.csv', index=False, encoding='utf-8-sig')
        print(f"[OK] 포도당 차단 대사물질 저장: glucose_blocked_metabolites.csv")
    
    if acetate_blocked:
        df_acetate = pd.DataFrame(acetate_blocked)
        df_acetate.to_csv('acetate_blocked_metabolites.csv', index=False, encoding='utf-8-sig')
        print(f"[OK] 아세트산 차단 대사물질 저장: acetate_blocked_metabolites.csv")
    
    print("\n" + "="*70)
    print("다음 단계: missing_reactions_identified.csv를 확인하여")
    print("각 반응을 수행하는 효소를 찾아보세요.")
    print("="*70)

if __name__ == "__main__":
    main()
