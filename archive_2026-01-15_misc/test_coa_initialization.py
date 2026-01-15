#!/usr/bin/env python
"""
CoA 초기화 경로 테스트
Acetyl-CoA 없이 CoA를 생성할 수 있는지 확인
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

def test_coa_production_without_accoa(model):
    """Acetyl-CoA 없이 CoA 생산 가능 여부"""
    print("="*70)
    print("CoA 생산 테스트 (Acetyl-CoA 없이)")
    print("="*70)
    
    model = setup_medium(model)
    
    # Acetyl-CoA 관련 반응 차단 (CoA 생산 경로 테스트용)
    try:
        accoa_c = model.metabolites.get_by_id('accoa_c')
        # Acetyl-CoA를 사용하는 반응들 확인
        accoa_rxns = list(accoa_c.reactions)
        print(f"Acetyl-CoA 관련 반응: {len(accoa_rxns)}개")
        
        # 일단 차단하지 말고, CoA 생산 테스트
        coa_c = model.metabolites.get_by_id('coa_c')
        
        # CoA 생산 테스트
        test_coa = Reaction('TEST_coa')
        test_coa.lower_bound = 0
        test_coa.upper_bound = 1000
        test_coa.add_metabolites({coa_c: -1})
        model.add_reactions([test_coa])
        model.objective = 'TEST_coa'
        
        solution = model.optimize()
        
        print(f"\nCoA 생산 테스트:")
        print(f"  상태: {solution.status}")
        print(f"  Objective: {solution.objective_value:.6f}")
        
        if solution.status == 'optimal' and solution.objective_value > 1e-6:
            print(f"  [OK] CoA 생산 가능!")
            
            # 어떤 경로로 생산되는지 확인
            print(f"\n  활성화된 CoA 생산 반응:")
            coa_producing = [r for r in coa_c.reactions 
                           if coa_c in r.products or (r.reversibility and coa_c in r.reactants)]
            for rxn in coa_producing:
                flux = solution.fluxes.get(rxn.id, 0)
                if abs(flux) > 1e-6:
                    print(f"    {rxn.id}: {flux:.4f} | {rxn.reaction}")
        else:
            print(f"  [FAIL] CoA 생산 불가")
            print(f"  [원인] Acetyl-CoA 없이는 CoA를 생성할 수 없음")
        
        model.remove_reactions([test_coa])
        
    except KeyError as e:
        print(f"[ERROR] {e}")

def test_complete_pathway_with_minimal_atp(model):
    """최소 ATP로 완전한 경로 테스트"""
    print("\n" + "="*70)
    print("완전한 경로 테스트 (최소 ATP 제공)")
    print("="*70)
    
    model = setup_medium(model)
    
    # 매우 작은 ATP 제공
    atp_c = model.metabolites.get_by_id('atp_c')
    atp_boot = Reaction('EX_atp_min')
    atp_boot.lower_bound = -0.01  # 매우 작은 양
    atp_boot.upper_bound = 0
    atp_boot.add_metabolites({atp_c: -1})
    model.add_reactions([atp_boot])
    
    # Biomass objective
    model.objective = 'Growth'
    
    solution = model.optimize()
    
    print(f"ATP 0.01 제공 후 성장 테스트:")
    print(f"  상태: {solution.status}")
    print(f"  Biomass: {solution.objective_value:.6f}")
    
    if solution.objective_value > 1e-6:
        print(f"  [SUCCESS] 성장 가능!")
        print(f"\n  주요 플럭스:")
        for rxn_id in ['EX_ac_e', 'ACt', 'ACS', 'CS', 'ICL', 'MALS', 'MDH', 'Growth']:
            flux = solution.fluxes.get(rxn_id, 0)
            if abs(flux) > 1e-6:
                print(f"    {rxn_id}: {flux:.4f}")
    else:
        print(f"  [FAIL] 성장 불가")
        print(f"\n  원인 확인:")
        print(f"    EX_ac_e: {solution.fluxes.get('EX_ac_e', 0):.6f}")
        print(f"    ACt: {solution.fluxes.get('ACt', 0):.6f}")
        print(f"    ACS: {solution.fluxes.get('ACS', 0):.6f}")
        
        # CoA 상태 확인
        try:
            coa_c = model.metabolites.get_by_id('coa_c')
            coa_producing = [r for r in coa_c.reactions 
                           if coa_c in r.products or (r.reversibility and coa_c in r.reactants)]
            active_coa = False
            for rxn in coa_producing:
                flux = solution.fluxes.get(rxn.id, 0)
                if abs(flux) > 1e-6:
                    active_coa = True
                    print(f"    CoA 생산 반응 활성: {rxn.id} ({flux:.4f})")
                    break
            if not active_coa:
                print(f"    [문제] CoA 생산 반응이 활성화되지 않음")
        except:
            pass
    
    model.remove_reactions([atp_boot])

def check_succoa_coa_transferase(model):
    """SUCOAACTr 반응 확인 (CoA 초기화 가능성)"""
    print("\n" + "="*70)
    print("SUCOAACTr 반응 확인 (CoA 초기화 가능성)")
    print("="*70)
    
    try:
        sucoaacr = model.reactions.get_by_id('SUCOAACTr')
        print(f"[OK] SUCOAACTr 발견:")
        print(f"  반응식: {sucoaacr.reaction}")
        print(f"  가역성: {sucoaacr.reversibility}")
        print(f"  하한: {sucoaacr.lower_bound}, 상한: {sucoaacr.upper_bound}")
        
        # 이 반응은 ac_c + succoa_c <=> accoa_c + succ_c
        # 역방향으로: accoa_c + succ_c <=> ac_c + succoa_c
        # 즉, Succinyl-CoA가 있으면 Acetyl-CoA를 만들 수 있음!
        
        print(f"\n  [분석]")
        print(f"    정방향: ac_c + succoa_c → accoa_c + succ_c")
        print(f"    역방향: accoa_c + succ_c → ac_c + succoa_c")
        print(f"    → Succinyl-CoA가 있으면 CoA transfer 가능")
        
        # Succinyl-CoA 생산 경로 확인
        try:
            succoa_c = model.metabolites.get_by_id('succoa_c')
            succoa_producing = [r for r in succoa_c.reactions 
                              if succoa_c in r.products or (r.reversibility and succoa_c in r.reactants)]
            print(f"\n  Succinyl-CoA 생산 반응: {len(succoa_producing)}개")
            for rxn in succoa_producing[:5]:
                print(f"    {rxn.id}: {rxn.reaction}")
        except KeyError:
            pass
            
    except KeyError:
        print("[NOT FOUND] SUCOAACTr 반응 없음")

def main():
    model = load_model("BaseModel.xml")
    
    test_coa_production_without_accoa(model)
    check_succoa_coa_transferase(model)
    test_complete_pathway_with_minimal_atp(model)
    
    print("\n" + "="*70)
    print("CoA 초기화 경로 분석 완료")
    print("="*70)

if __name__ == "__main__":
    main()

