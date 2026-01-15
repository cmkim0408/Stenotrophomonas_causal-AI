#!/usr/bin/env python
"""
9개 반응 추가 후 FBA 실패 원인 진단
"""

import cobra
from pathlib import Path
from cobra.flux_analysis import find_essential_reactions, find_blocked_reactions

def load_model(model_path):
    """모델 로드"""
    model = cobra.io.read_sbml_model(str(model_path))
    return model

def add_reaction_to_model(model, ref_model, rxn_id):
    """레퍼런스 모델에서 반응을 신규 모델에 추가"""
    if rxn_id in model.reactions:
        return False
    
    if rxn_id not in ref_model.reactions:
        return False
    
    ref_rxn = ref_model.reactions.get_by_id(rxn_id)
    
    # 대사물질 먼저 확인/추가
    metabolites_dict = {}
    for met in ref_rxn.metabolites:
        if met.id not in model.metabolites:
            new_met = cobra.Metabolite(met.id, formula=met.formula, name=met.name, 
                                      compartment=met.compartment, charge=met.charge)
            model.add_metabolites([new_met])
        metabolites_dict[model.metabolites.get_by_id(met.id)] = ref_rxn.metabolites[met]
    
    # 반응 생성 및 추가
    new_rxn = cobra.Reaction(rxn_id)
    new_rxn.name = ref_rxn.name
    new_rxn.lower_bound = ref_rxn.lower_bound
    new_rxn.upper_bound = ref_rxn.upper_bound
    new_rxn.add_metabolites(metabolites_dict)
    
    model.add_reactions([new_rxn])
    
    # 유전자 연결
    for gene in ref_rxn.genes:
        if gene.id in model.genes:
            model.genes.get_by_id(gene.id).reactions.add(new_rxn)
    
    return True

def setup_acetate_medium(model):
    """Acetate 기반 미디어 설정"""
    exchanges = {
        'EX_ac_e': (-10, 1000),
        'EX_o2_e': (-1000, 1000),
        'EX_h2o_e': (-1000, 1000),
        'EX_h_e': (-1000, 1000),
        'EX_nh4_e': (-10, 1000),
        'EX_pi_e': (-10, 1000),
        'EX_so4_e': (-10, 1000),
        'EX_k_e': (-1000, 1000),
        'EX_na1_e': (-1000, 1000),
        'EX_mg2_e': (-1000, 1000),
        'EX_ca2_e': (-1000, 1000),
        'EX_fe2_e': (-10, 1000),
        'EX_fe3_e': (-10, 1000),
        'EX_hco3_e': (-1000, 1000),
        'EX_co2_e': (-1000, 1000),
    }
    
    for ex_id, bounds in exchanges.items():
        if ex_id in model.reactions:
            model.reactions.get_by_id(ex_id).bounds = bounds
    
    return model

def diagnose_growth_failure(model):
    """성장 실패 원인 진단"""
    print("\n" + "="*70)
    print("성장 실패 원인 진단")
    print("="*70)
    
    # Biomass 반응 찾기
    biomass_rxns = [r for r in model.reactions if 'growth' in r.id.lower() or 'biomass' in r.id.lower()]
    if not biomass_rxns:
        print("[ERROR] Biomass 반응을 찾을 수 없습니다")
        return
    
    biomass_rxn = biomass_rxns[0]
    model.objective = biomass_rxn.id
    
    # 1. Biomass 구성 요소 확인
    print("\n1. Biomass 반응 구성 요소 확인:")
    biomass_components = []
    for met, coeff in biomass_rxn.metabolites.items():
        biomass_components.append((met.id, coeff))
        print(f"  {met.id}: {coeff:.4f}")
    
    # 2. 각 구성 요소 생산 가능 여부 확인
    print("\n2. Biomass 구성 요소 생산 가능 여부:")
    from cobra.flux_analysis import flux_variability_analysis
    
    for met_id, coeff in biomass_components[:20]:  # 처음 20개만
        if met_id in model.metabolites:
            met = model.metabolites.get_by_id(met_id)
            # 해당 대사물질 생산 가능 여부 확인
            try:
                # 생산 반응 찾기
                producers = [r for r in met.reactions if met in r.reactions and r.reactions[met] > 0]
                if producers:
                    # 생산 가능성 테스트
                    test_rxn = cobra.Reaction(f"TEST_{met_id}")
                    test_rxn.add_metabolites({met: 1})
                    test_rxn.lower_bound = 0
                    test_rxn.upper_bound = 1000
                    model.add_reactions([test_rxn])
                    model.objective = test_rxn.id
                    test_sol = model.optimize()
                    model.remove_reactions([test_rxn])
                    model.objective = biomass_rxn.id
                    
                    if test_sol.status == 'optimal' and test_sol.objective_value > 1e-6:
                        print(f"  {met_id}: 생산 가능")
                    else:
                        print(f"  {met_id}: 생산 불가능 (blocked)")
                else:
                    print(f"  {met_id}: 생산 반응 없음")
            except:
                print(f"  {met_id}: 확인 실패")
    
    # 3. Blocked reactions 확인
    print("\n3. Blocked reactions 확인:")
    try:
        blocked = find_blocked_reactions(model)
        print(f"  Blocked reactions 수: {len(blocked)}")
        if len(blocked) < 50:
            print("  일부 blocked reactions:")
            for rxn in list(blocked)[:20]:
                print(f"    {rxn.id}")
    except Exception as e:
        print(f"  Blocked reactions 분석 실패: {e}")
    
    # 4. 핵심 경로 확인
    print("\n4. 핵심 경로 반응 확인:")
    key_reactions = ['ACS_ADP', 'PEPCK_ATP', 'SUCDi', 'CS', 'MDH', 'ICL', 'MALS', 'ACONT']
    for rxn_id in key_reactions:
        if rxn_id in model.reactions:
            rxn = model.reactions.get_by_id(rxn_id)
            print(f"  {rxn_id}: 존재 (bounds: [{rxn.lower_bound}, {rxn.upper_bound}])")
        else:
            print(f"  {rxn_id}: 누락")

def main():
    base_path = Path(__file__).parent.parent
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    new_model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    
    print("="*70)
    print("FBA 실패 원인 진단")
    print("="*70)
    
    # 모델 로드
    ref_model = load_model(ref_model_path)
    new_model = load_model(new_model_path)
    
    # 9개 반응 추가
    missing_reactions = [
        'ACS_ADP', 'PEPCK_ATP', 'SUCDi', 'ACtexi', 'EX_hco3_e',
        'T_hco3_e_to_c', 'T_o2_e_to_o2_c', 'T_nh4_e_to_nh4_c', 'T_fe3_e_to_fe3_c'
    ]
    
    print("\n반응 추가 중...")
    for rxn_id in missing_reactions:
        add_reaction_to_model(new_model, ref_model, rxn_id)
    
    # 미디어 설정
    new_model = setup_acetate_medium(new_model)
    
    # 진단
    diagnose_growth_failure(new_model)
    
    print("\n" + "="*70)
    print("결론")
    print("="*70)
    print("9개 반응만으로는 부족합니다.")
    print("추가로 필요한 것들:")
    print("  - 더 많은 대사 경로 반응")
    print("  - 보조인자 생합성 경로")
    print("  - 아미노산/뉴클레오티드 생합성 경로")
    print("  - 대사물질 연결 확인")

if __name__ == "__main__":
    main()
