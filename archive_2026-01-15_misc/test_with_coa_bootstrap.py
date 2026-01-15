#!/usr/bin/env python
"""
CoA를 초기 조건으로 제공하여 아세트산으로부터 ATP 생산 테스트
실제 생물학적으로는 세포에 초기 CoA가 존재함
"""

import cobra
from cobra import Reaction

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    return model

def setup_acetate_medium(model, provide_coa=False, coa_amount=0.01):
    """Acetate medium 설정 (선택적으로 CoA 제공)"""
    for rxn in model.exchanges:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    
    essentials = {
        'EX_ac_e': (-1000, 1000),
        'EX_nh4_e': (-1000, 0),
        'EX_h2o_e': (-1000, 0),
        'EX_h_e': (-1000, 0),
        'EX_pi_e': (-1000, 0),
        'EX_so4_e': (-1000, 0),
        'EX_k_e': (-1000, 0),
        'EX_mg2_e': (-1000, 0),
        'EX_fe2_e': (-1000, 0),
        'EX_mn2_e': (-1000, 0),
        'EX_zn2_e': (-1000, 0),
        'EX_co2_e': (-1000, 1000),
        'EX_o2_e': (-1000, 1000)
    }
    
    for ex_id, bounds in essentials.items():
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            ex_rxn.lower_bound = bounds[0]
            ex_rxn.upper_bound = bounds[1]
        except KeyError:
            pass
    
    # CoA 제공 (선택적)
    if provide_coa:
        coa_c = model.metabolites.get_by_id('coa_c')
        ex_coa = Reaction('EX_coa_boot')
        ex_coa.lower_bound = -coa_amount
        ex_coa.upper_bound = 0
        ex_coa.add_metabolites({coa_c: -1})
        model.add_reactions([ex_coa])
        return model, ex_coa
    else:
        return model, None

def test_atp_production(model, provide_coa=False):
    """ATP 생산 테스트"""
    print("="*70)
    if provide_coa:
        print("ATP 생산 테스트 (CoA 초기 제공)")
    else:
        print("ATP 생산 테스트 (CoA 없음)")
    print("="*70)
    
    model, ex_coa = setup_acetate_medium(model, provide_coa=provide_coa)
    
    # ATP 생산을 objective로
    atp_c = model.metabolites.get_by_id('atp_c')
    test_atp = Reaction('TEST_atp')
    test_atp.lower_bound = 0
    test_atp.upper_bound = 1000
    test_atp.add_metabolites({atp_c: -1})
    model.add_reactions([test_atp])
    model.objective = 'TEST_atp'
    
    solution = model.optimize()
    
    print(f"\n결과:")
    print(f"  상태: {solution.status}")
    print(f"  ATP 생산: {solution.objective_value:.6f}")
    
    if solution.objective_value > 1e-6:
        print(f"  [SUCCESS] ATP 생산 가능!")
        print(f"\n  주요 경로 플럭스:")
        key_rxns = ['EX_ac_e', 'ACt', 'ACS', 'CS', 'ICDHx', 'AKGDH', 'SUCOAS', 
                   'SUCD', 'FUM', 'MDH', 'ATPS4rpp']
        for rxn_id in key_rxns:
            flux = solution.fluxes.get(rxn_id, 0)
            if abs(flux) > 1e-6:
                print(f"    {rxn_id}: {flux:.4f}")
        
        if ex_coa:
            print(f"    EX_coa_boot: {solution.fluxes.get('EX_coa_boot', 0):.4f}")
    else:
        print(f"  [FAIL] ATP 생산 불가")
    
    model.remove_reactions([test_atp])
    if ex_coa:
        model.remove_reactions([ex_coa])
    
    return solution

def test_growth(model, provide_coa=False):
    """성장 테스트"""
    print("\n" + "="*70)
    if provide_coa:
        print("성장 테스트 (CoA 초기 제공)")
    else:
        print("성장 테스트 (CoA 없음)")
    print("="*70)
    
    model, ex_coa = setup_acetate_medium(model, provide_coa=provide_coa)
    model.objective = 'Growth'
    
    solution = model.optimize()
    
    print(f"\n결과:")
    print(f"  상태: {solution.status}")
    print(f"  Biomass: {solution.objective_value:.6f}")
    
    if solution.objective_value > 1e-6:
        print(f"  [SUCCESS] 성장 가능!")
        print(f"\n  주요 플럭스:")
        for rxn_id in ['EX_ac_e', 'ACt', 'ACS', 'CS', 'ICL', 'MALS', 'ME1', 'ME2', 'Growth']:
            flux = solution.fluxes.get(rxn_id, 0)
            if abs(flux) > 1e-6:
                print(f"    {rxn_id}: {flux:.4f}")
        if ex_coa:
            print(f"    EX_coa_boot: {solution.fluxes.get('EX_coa_boot', 0):.4f}")
    else:
        print(f"  [FAIL] 성장 불가")
    
    if ex_coa:
        model.remove_reactions([ex_coa])
    
    return solution

def main():
    model = load_model("BaseModel.xml")
    
    print("\n" + "="*70)
    print("아세트산만으로 ATP 생산 가능 여부 테스트")
    print("="*70)
    print("\n[중요] 실제 생물학적으로는 세포에 초기 CoA가 존재합니다.")
    print("모델에서도 CoA를 초기 조건으로 제공하면 아세트산으로부터")
    print("ATP 생산 및 성장이 가능한지 테스트합니다.\n")
    
    # Test 1: CoA 없이
    test_atp_production(model, provide_coa=False)
    
    # Test 2: CoA 제공
    test_atp_production(model, provide_coa=True, coa_amount=0.01)
    
    # Test 3: 성장 테스트 (CoA 제공)
    test_growth(model, provide_coa=True, coa_amount=0.01)
    
    print("\n" + "="*70)
    print("결론:")
    print("="*70)
    print("아세트산만으로 ATP를 만들 수 있는지?")
    print("→ CoA가 초기에 존재한다면 가능해야 합니다.")
    print("→ ACS 반응: ac_c + ATP + CoA → Acetyl-CoA")
    print("→ TCA cycle → NADH/FADH2 → 전자전달계 → ATP")
    print("→ 순환이 시작되면 지속적으로 ATP 생성 가능")
    print("="*70)

if __name__ == "__main__":
    main()

