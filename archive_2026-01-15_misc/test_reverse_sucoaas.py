#!/usr/bin/env python
"""
SUCOAS 역방향 테스트: Succinyl-CoA → ATP + CoA
이것이 CoA 초기화 경로가 될 수 있음
"""

import cobra
from cobra import Reaction

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    return model

def setup_medium(model):
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

def test_sucoaas_reverse(model):
    """SUCOAS 역방향 테스트"""
    print("="*70)
    print("SUCOAS 역방향 테스트 (Succinyl-CoA → ATP + CoA)")
    print("="*70)
    
    model = setup_medium(model)
    
    try:
        sucoaas = model.reactions.get_by_id('SUCOAS')
        print(f"SUCOAS 반응: {sucoaas.reaction}")
        print(f"  가역성: {sucoaas.reversibility}")
        print(f"  정방향: atp_c + coa_c + succ_c → adp_c + pi_c + succoa_c")
        print(f"  역방향: adp_c + pi_c + succoa_c → atp_c + coa_c + succ_c")
        
        # Succinyl-CoA 제공
        succoa_c = model.metabolites.get_by_id('succoa_c')
        succoa_boot = Reaction('EX_succoa')
        succoa_boot.lower_bound = -0.1
        succoa_boot.upper_bound = 0
        succoa_boot.add_metabolites({succoa_c: -1})
        model.add_reactions([succoa_boot])
        
        # ATP + CoA 생산 테스트
        atp_c = model.metabolites.get_by_id('atp_c')
        coa_c = model.metabolites.get_by_id('coa_c')
        
        test_atp_coa = Reaction('TEST_atp_coa')
        test_atp_coa.lower_bound = 0
        test_atp_coa.upper_bound = 1000
        test_atp_coa.add_metabolites({atp_c: -1, coa_c: -1})
        model.add_reactions([test_atp_coa])
        model.objective = 'TEST_atp_coa'
        
        solution = model.optimize()
        
        print(f"\nSUCOAS 역방향 테스트 결과:")
        print(f"  상태: {solution.status}")
        print(f"  ATP + CoA 생산: {solution.objective_value:.6f}")
        print(f"  SUCOAS 플럭스: {solution.fluxes.get('SUCOAS', 0):.6f}")
        
        if solution.objective_value > 1e-6:
            print(f"  [SUCCESS] SUCOAS 역방향으로 ATP + CoA 생산 가능!")
            print(f"    → 이것이 CoA 초기화 경로가 될 수 있음!")
        else:
            print(f"  [FAIL] SUCOAS 역방향으로도 생산 불가")
        
        model.remove_reactions([succoa_boot, test_atp_coa])
        
    except KeyError as e:
        print(f"[ERROR] {e}")

def test_succinyl_coa_production(model):
    """Succinyl-CoA 생산 경로 확인"""
    print("\n" + "="*70)
    print("Succinyl-CoA 생산 경로 확인")
    print("="*70)
    
    model = setup_medium(model)
    
    try:
        succoa_c = model.metabolites.get_by_id('succoa_c')
        
        # Succinyl-CoA 생산 반응
        producing = [r for r in succoa_c.reactions 
                    if succoa_c in r.products or (r.reversibility and succoa_c in r.reactants)]
        
        print(f"Succinyl-CoA 생산 반응: {len(producing)}개\n")
        
        # 각 반응이 작동 가능한지 확인
        for rxn in producing[:10]:
            print(f"{rxn.id}: {rxn.reaction}")
            print(f"  가역성: {rxn.reversibility}")
            
            # 필요한 대사물질 확인
            required = {m.id: -coeff for m, coeff in rxn.metabolites.items() if coeff < 0}
            if required:
                print(f"  필요한 대사물질: {list(required.keys())}")
        
        # Succinyl-CoA 생산 테스트
        test_succoa = Reaction('TEST_succoa')
        test_succoa.lower_bound = 0
        test_succoa.upper_bound = 1000
        test_succoa.add_metabolites({succoa_c: -1})
        model.add_reactions([test_succoa])
        model.objective = 'TEST_succoa'
        
        solution = model.optimize()
        print(f"\nSuccinyl-CoA 생산 테스트:")
        print(f"  상태: {solution.status}")
        print(f"  Objective: {solution.objective_value:.6f}")
        
        if solution.objective_value > 1e-6:
            print(f"  [OK] Succinyl-CoA 생산 가능!")
            # 어떤 경로로 생산되는지
            for rxn in producing:
                flux = solution.fluxes.get(rxn.id, 0)
                if abs(flux) > 1e-6:
                    print(f"    {rxn.id}: {flux:.4f}")
        else:
            print(f"  [FAIL] Succinyl-CoA 생산 불가")
        
        model.remove_reactions([test_succoa])
        
    except KeyError:
        print("[ERROR] succoa_c 없음")

def test_complete_cycle(model):
    """완전한 순환 테스트"""
    print("\n" + "="*70)
    print("완전한 순환 테스트")
    print("="*70)
    
    model = setup_medium(model)
    
    # 매우 작은 Succinyl-CoA 제공 (bootstrap)
    try:
        succoa_c = model.metabolites.get_by_id('succoa_c')
        succoa_boot = Reaction('EX_succoa_boot')
        succoa_boot.lower_bound = -0.01
        succoa_boot.upper_bound = 0
        succoa_boot.add_metabolites({succoa_c: -1})
        model.add_reactions([succoa_boot])
        
        # Biomass objective
        model.objective = 'Growth'
        
        solution = model.optimize()
        
        print(f"Succinyl-CoA 0.01 제공 후 성장 테스트:")
        print(f"  상태: {solution.status}")
        print(f"  Biomass: {solution.objective_value:.6f}")
        
        if solution.objective_value > 1e-6:
            print(f"  [SUCCESS] 성장 가능!")
            print(f"\n  주요 플럭스:")
            for rxn_id in ['EX_ac_e', 'ACt', 'ACS', 'SUCOAS', 'SUCOAACTr', 'CS', 'Growth']:
                flux = solution.fluxes.get(rxn_id, 0)
                if abs(flux) > 1e-6:
                    print(f"    {rxn_id}: {flux:.4f}")
        else:
            print(f"  [FAIL] 성장 불가")
        
        model.remove_reactions([succoa_boot])
        
    except Exception as e:
        print(f"[ERROR] {e}")

def main():
    model = load_model("BaseModel.xml")
    
    test_sucoaas_reverse(model)
    test_succinyl_coa_production(model)
    test_complete_cycle(model)
    
    print("\n" + "="*70)
    print("SUCOAS 역방향 경로 분석 완료")
    print("="*70)

if __name__ == "__main__":
    main()

