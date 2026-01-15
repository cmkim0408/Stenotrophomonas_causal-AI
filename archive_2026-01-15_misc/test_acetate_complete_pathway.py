#!/usr/bin/env python
"""
아세트산만으로 ATP 생산 - 완전한 경로 분석
모든 필요한 초기 조건을 하나씩 추가하며 테스트
"""

import cobra
from cobra import Reaction

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    return model

def setup_medium(model):
    """기본 medium 설정"""
    for rxn in model.exchanges:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    
    for ex_id, bounds in [
        ('EX_ac_e', (-1000, 1000)),
        ('EX_nh4_e', (-1000, 0)),
        ('EX_h2o_e', (-1000, 0)),
        ('EX_h_e', (-1000, 0)),
        ('EX_pi_e', (-1000, 0)),
        ('EX_so4_e', (-1000, 0)),
        ('EX_k_e', (-1000, 0)),
        ('EX_mg2_e', (-1000, 0)),
        ('EX_fe2_e', (-1000, 0)),
        ('EX_o2_e', (-1000, 1000))
    ]:
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            ex_rxn.lower_bound = bounds[0]
            ex_rxn.upper_bound = bounds[1]
        except KeyError:
            pass
    
    return model

def test_acs_requirements(model):
    """ACS 작동에 필요한 것들 확인"""
    print("="*70)
    print("ACS 작동 요구사항 분석")
    print("="*70)
    
    model = setup_medium(model)
    
    try:
        acs = model.reactions.get_by_id('ACS')
        print(f"ACS 반응: {acs.reaction}")
        print(f"\n필요한 대사물질:")
        
        required = {}
        for met, coeff in acs.metabolites.items():
            if coeff < 0:
                required[met.id] = -coeff
                print(f"  {met.id}: {-coeff}")
        
        # 각각이 제공 가능한지 테스트
        print(f"\n각 대사물질 제공 가능 여부:")
        
        # Test 1: ac_c
        ac_c = model.metabolites.get_by_id('ac_c')
        test_ac = Reaction('TEST_ac')
        test_ac.lower_bound = 0
        test_ac.upper_bound = 1000
        test_ac.add_metabolites({ac_c: -1})
        model.add_reactions([test_ac])
        model.objective = 'TEST_ac'
        sol_ac = model.optimize()
        print(f"  ac_c: {sol_ac.objective_value:.6f} {'[OK]' if sol_ac.objective_value > 1e-6 else '[FAIL]'}")
        model.remove_reactions([test_ac])
        
        # Test 2: atp_c
        atp_c = model.metabolites.get_by_id('atp_c')
        test_atp = Reaction('TEST_atp')
        test_atp.lower_bound = 0
        test_atp.upper_bound = 1000
        test_atp.add_metabolites({atp_c: -1})
        model.add_reactions([test_atp])
        model.objective = 'TEST_atp'
        sol_atp = model.optimize()
        print(f"  atp_c: {sol_atp.objective_value:.6f} {'[OK]' if sol_atp.objective_value > 1e-6 else '[FAIL]'}")
        model.remove_reactions([test_atp])
        
        # Test 3: coa_c
        coa_c = model.metabolites.get_by_id('coa_c')
        test_coa = Reaction('TEST_coa')
        test_coa.lower_bound = 0
        test_coa.upper_bound = 1000
        test_coa.add_metabolites({coa_c: -1})
        model.add_reactions([test_coa])
        model.objective = 'TEST_coa'
        sol_coa = model.optimize()
        print(f"  coa_c: {sol_coa.objective_value:.6f} {'[OK]' if sol_coa.objective_value > 1e-6 else '[FAIL]'}")
        model.remove_reactions([test_coa])
        
    except KeyError as e:
        print(f"[ERROR] {e}")

def test_tca_cycle_start(model):
    """TCA cycle 시작을 위한 요구사항"""
    print("\n" + "="*70)
    print("TCA Cycle 시작 요구사항 분석")
    print("="*70)
    
    model = setup_medium(model)
    
    try:
        cs = model.reactions.get_by_id('CS')
        print(f"CS 반응 (Citrate Synthase): {cs.reaction}")
        print(f"\n필요한 대사물질:")
        
        for met, coeff in cs.metabolites.items():
            if coeff < 0:
                print(f"  {met.id}: {-coeff}")
                
                # OAA 확인
                if met.id == 'oaa_c':
                    print(f"    [중요] OAA는 TCA cycle의 시작점!")
                    
                    # OAA 생산 경로 확인
                    oaa_c = met
                    oaa_producing = [r for r in oaa_c.reactions 
                                   if oaa_c in r.products or (r.reversibility and oaa_c in r.reactants)]
                    print(f"    OAA 생산 반응: {len(oaa_producing)}개")
                    print(f"    주요 반응:")
                    for rxn in oaa_producing[:5]:
                        print(f"      {rxn.id}: {rxn.reaction}")
                    
                    # OAA 생산 테스트
                    test_oaa = Reaction('TEST_oaa')
                    test_oaa.lower_bound = 0
                    test_oaa.upper_bound = 1000
                    test_oaa.add_metabolites({oaa_c: -1})
                    model.add_reactions([test_oaa])
                    model.objective = 'TEST_oaa'
                    sol_oaa = model.optimize()
                    print(f"    OAA 생산 가능 여부: {sol_oaa.objective_value:.6f} {'[OK]' if sol_oaa.objective_value > 1e-6 else '[FAIL]'}")
                    
                    if sol_oaa.objective_value < 1e-6:
                        print(f"    [ROOT CAUSE] OAA를 생성할 수 없습니다!")
                        print(f"      → TCA cycle이 시작되지 않음")
                        print(f"      → 아세트산만으로는 OAA를 만들 수 없을 수 있음")
                    
                    model.remove_reactions([test_oaa])
        
    except KeyError as e:
        print(f"[ERROR] {e}")

def test_with_all_bootstrap(model):
    """모든 필요한 초기 조건 제공"""
    print("\n" + "="*70)
    print("완전한 Bootstrap 테스트 (CoA + OAA 제공)")
    print("="*70)
    
    model = setup_medium(model)
    
    # Bootstrap 대사물질 제공
    coa_c = model.metabolites.get_by_id('coa_c')
    oaa_c = model.metabolites.get_by_id('oaa_c')
    
    ex_coa = Reaction('EX_coa')
    ex_coa.lower_bound = -0.01
    ex_coa.add_metabolites({coa_c: -1})
    
    ex_oaa = Reaction('EX_oaa')
    ex_oaa.lower_bound = -0.01
    ex_oaa.add_metabolites({oaa_c: -1})
    
    model.add_reactions([ex_coa, ex_oaa])
    
    # ATP 생산 테스트
    atp_c = model.metabolites.get_by_id('atp_c')
    test_atp = Reaction('TEST_atp')
    test_atp.lower_bound = 0
    test_atp.upper_bound = 1000
    test_atp.add_metabolites({atp_c: -1})
    model.add_reactions([test_atp])
    model.objective = 'TEST_atp'
    
    solution = model.optimize()
    
    print(f"CoA + OAA 제공 후 ATP 생산:")
    print(f"  상태: {solution.status}")
    print(f"  ATP 생산: {solution.objective_value:.6f}")
    
    if solution.objective_value > 0.01:
        print(f"  [SUCCESS] Bootstrap으로 ATP 생산 가능!")
        print(f"    → 순환이 시작되면 지속 가능")
        
        print(f"\n  주요 플럭스:")
        for rxn_id in ['EX_ac_e', 'ACt', 'ACS', 'CS', 'ICDHx', 'AKGDH', 'SUCOAS', 'MDH', 'ATPS4rpp']:
            flux = solution.fluxes.get(rxn_id, 0)
            if abs(flux) > 1e-6:
                print(f"    {rxn_id}: {flux:.4f}")
    else:
        print(f"  [FAIL] Bootstrap으로도 ATP 생산 불가")
    
    # Biomass 테스트
    model.objective = 'Growth'
    sol_growth = model.optimize()
    print(f"\n  Biomass: {sol_growth.objective_value:.6f}")
    
    model.remove_reactions([ex_coa, ex_oaa, test_atp])

def main():
    model = load_model("BaseModel.xml")
    
    test_acs_requirements(model)
    test_tca_cycle_start(model)
    test_with_all_bootstrap(model)
    
    print("\n" + "="*70)
    print("결론")
    print("="*70)
    print("아세트산만으로 ATP를 만들 수 있는가?")
    print("→ 이론적으로는 가능하지만, 모델에서는:")
    print("  1. CoA 초기화 필요")
    print("  2. OAA 초기화 필요 (또는 OAA 생산 경로 필요)")
    print("→ 실제 생물학적으로는 세포에 초기 대사물질이 존재")
    print("→ 모델에서는 'bootstrap' 대사물질을 제공해야 할 수 있음")
    print("="*70)

if __name__ == "__main__":
    main()

