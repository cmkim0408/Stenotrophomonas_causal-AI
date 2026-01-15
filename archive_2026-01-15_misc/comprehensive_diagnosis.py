#!/usr/bin/env python
"""
종합 진단: ATP, CoA, PEP 생산 경로 확인
"""

import cobra

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    return model

def setup_medium(model):
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
    
    return model

def check_atp_production(model):
    """ATP 생산 경로 확인"""
    print("="*70)
    print("ATP 생산 경로 분석")
    print("="*70)
    
    model = setup_medium(model)
    
    try:
        atp_c = model.metabolites.get_by_id('atp_c')
        
        # ATP 생산 반응
        producing = [r for r in atp_c.reactions if atp_c in r.products or (r.reversibility and atp_c in r.reactants)]
        print(f"ATP 생산 반응: {len(producing)}개\n")
        print("주요 ATP 생산 반응 (상위 10개):")
        for rxn in producing[:10]:
            print(f"  {rxn.id}: {rxn.reaction}")
            print(f"    가역성: {rxn.reversibility}")
        
        # ATP 생산 테스트
        from cobra import Reaction
        test_atp = Reaction('TEST_atp')
        test_atp.lower_bound = 0
        test_atp.upper_bound = 1000
        test_atp.add_metabolites({atp_c: -1})
        model.add_reactions([test_atp])
        model.objective = 'TEST_atp'
        
        solution = model.optimize()
        print(f"\nATP 생산 테스트:")
        print(f"  상태: {solution.status}")
        print(f"  Objective: {solution.objective_value:.6f}")
        
        if solution.status == 'optimal' and solution.objective_value > 1e-6:
            print(f"  [OK] ATP 생산 가능!")
            # 주요 플럭스 확인
            for rxn_id in ['EX_ac_e', 'ACt', 'ACS']:
                flux = solution.fluxes.get(rxn_id, 0)
                if abs(flux) > 1e-6:
                    print(f"    {rxn_id}: {flux:.4f}")
        else:
            print(f"  [FAIL] ATP 생산 불가")
            print(f"  [ROOT CAUSE] ATP가 생성되지 않아 모든 에너지 의존 반응이 작동하지 않습니다!")
        
        model.remove_reactions([test_atp])
        
    except KeyError:
        print("[ERROR] atp_c 없음")

def check_coa_production(model):
    """CoA 생산 경로 확인"""
    print("\n" + "="*70)
    print("CoA 생산 경로 분석")
    print("="*70)
    
    model = setup_medium(model)
    
    try:
        coa_c = model.metabolites.get_by_id('coa_c')
        
        # CoA 생산 반응
        producing = [r for r in coa_c.reactions if coa_c in r.products or (r.reversibility and coa_c in r.reactants)]
        print(f"CoA 생산 반응: {len(producing)}개\n")
        print("주요 CoA 생산 반응 (상위 10개):")
        for rxn in producing[:10]:
            print(f"  {rxn.id}: {rxn.reaction}")
        
        # CoA 생산 테스트
        from cobra import Reaction
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
        else:
            print(f"  [FAIL] CoA 생산 불가")
            print(f"  [ROOT CAUSE] CoA가 생성되지 않아 ACS가 작동하지 않습니다!")
        
        model.remove_reactions([test_coa])
        
    except KeyError:
        print("[ERROR] coa_c 없음")

def check_pep_production(model):
    """PEP 생산 경로 확인"""
    print("\n" + "="*70)
    print("PEP 생산 경로 분석 (OAA bootstrap용)")
    print("="*70)
    
    model = setup_medium(model)
    
    try:
        pep_c = model.metabolites.get_by_id('pep_c')
        
        # PEP 생산 반응
        producing = [r for r in pep_c.reactions if pep_c in r.products or (r.reversibility and pep_c in r.reactants)]
        print(f"PEP 생산 반응: {len(producing)}개\n")
        print("주요 PEP 생산 반응 (상위 10개):")
        for rxn in producing[:10]:
            print(f"  {rxn.id}: {rxn.reaction}")
            print(f"    가역성: {rxn.reversibility}")
        
        # PEP 생산 테스트
        from cobra import Reaction
        test_pep = Reaction('TEST_pep')
        test_pep.lower_bound = 0
        test_pep.upper_bound = 1000
        test_pep.add_metabolites({pep_c: -1})
        model.add_reactions([test_pep])
        model.objective = 'TEST_pep'
        
        solution = model.optimize()
        print(f"\nPEP 생산 테스트:")
        print(f"  상태: {solution.status}")
        print(f"  Objective: {solution.objective_value:.6f}")
        
        if solution.status == 'optimal' and solution.objective_value > 1e-6:
            print(f"  [OK] PEP 생산 가능!")
            # PPC를 통한 OAA 생산 테스트
            try:
                ppc = model.reactions.get_by_id('PPC')
                print(f"\n  PPC 확인: {ppc.reaction}")
                
                # OAA 생산 테스트 (PPC 사용)
                oaa_c = model.metabolites.get_by_id('oaa_c')
                test_oaa = Reaction('TEST_oaa_ppc')
                test_oaa.lower_bound = 0
                test_oaa.upper_bound = 1000
                test_oaa.add_metabolites({oaa_c: -1})
                model.add_reactions([test_oaa])
                model.objective = 'TEST_oaa_ppc'
                
                solution2 = model.optimize()
                if solution2.status == 'optimal' and solution2.objective_value > 1e-6:
                    print(f"  [OK] PEP → OAA (via PPC) 가능: {solution2.objective_value:.6f}")
                else:
                    print(f"  [FAIL] PEP → OAA 불가")
                
                model.remove_reactions([test_oaa])
                
            except KeyError:
                print(f"  [ERROR] PPC 반응 없음")
        else:
            print(f"  [FAIL] PEP 생산 불가")
        
        model.remove_reactions([test_pep])
        
    except KeyError:
        print("[ERROR] pep_c 없음")

def summary_and_recommendations(model):
    """요약 및 권장사항"""
    print("\n" + "="*70)
    print("진단 요약 및 권장사항")
    print("="*70)
    
    print("\n[발견된 문제들]")
    print("  1. ATP 생산 불가 → 모든 에너지 의존 반응 차단")
    print("  2. CoA 생산 불가 → ACS 및 CoA 의존 반응 차단")
    print("  3. OAA bootstrap 문제 → TCA/Glyoxylate cycle 시작 불가")
    
    print("\n[가능한 해결책]")
    print("  1. ATP 생산 경로 확인 및 수정 필요")
    print("     - 아마도 Acetyl-CoA가 없어서 에너지 생성 불가")
    print("     - 순환 의존성 문제")
    print("  2. 초기 bootstrap 대사물질 제공")
    print("     - 소량의 ATP, CoA, 또는 OAA를 제공")
    print("  3. 또는 다른 초기화 경로 추가")

def main():
    model = load_model("BaseModel.xml")
    
    check_atp_production(model)
    check_coa_production(model)
    check_pep_production(model)
    summary_and_recommendations(model)
    
    print("\n" + "="*70)
    print("종합 진단 완료")
    print("="*70)

if __name__ == "__main__":
    main()

