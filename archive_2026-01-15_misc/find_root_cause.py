#!/usr/bin/env python
"""
근본 원인 찾기: 왜 ACS와 OAA 생산이 안 되는가?
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

def check_acs_requirements(model):
    """ACS 반응에 필요한 대사물질 확인"""
    print("="*70)
    print("ACS 반응 요구사항 분석")
    print("="*70)
    
    model = setup_medium(model)
    
    try:
        acs = model.reactions.get_by_id('ACS')
        print(f"ACS 반응: {acs.reaction}\n")
        
        # 필요한 대사물질들
        required = {m.id: -coeff for m, coeff in acs.metabolites.items() if coeff < 0}
        print("필요한 대사물질:")
        for met_id, req_coeff in required.items():
            print(f"  {met_id}: {req_coeff}")
            
            # 각 대사물질 생산 가능 여부 확인
            try:
                met = model.metabolites.get_by_id(met_id)
                from cobra import Reaction
                
                test_rxn = Reaction(f'TEST_{met_id}')
                test_rxn.lower_bound = 0
                test_rxn.upper_bound = 1000
                test_rxn.add_metabolites({met: -1})
                model.add_reactions([test_rxn])
                model.objective = test_rxn.id
                
                solution = model.optimize()
                if solution.status == 'optimal' and solution.objective_value > 1e-6:
                    print(f"    [OK] {met_id} 생산 가능: {solution.objective_value:.6f}")
                else:
                    print(f"    [FAIL] {met_id} 생산 불가")
                    # 어떤 반응이 필요한지 확인
                    producing = [r for r in met.reactions if met in r.products or (r.reversibility and met in r.reactants)]
                    print(f"      생산 반응 수: {len(producing)}")
                    for r in producing[:3]:
                        print(f"        {r.id}: {r.reaction}")
                
                model.remove_reactions([test_rxn])
                
            except KeyError:
                print(f"    [NOT FOUND] {met_id}")
        
    except KeyError:
        print("[ERROR] ACS 반응 없음")

def check_oaa_cycle_problem(model):
    """OAA 순환 문제 확인"""
    print("\n" + "="*70)
    print("OAA 순환 문제 분석")
    print("="*70)
    
    model = setup_medium(model)
    
    try:
        oaa_c = model.metabolites.get_by_id('oaa_c')
        print(f"OAA 생산 반응:")
        producing = [r for r in oaa_c.reactions if oaa_c in r.products or (r.reversibility and oaa_c in r.reactants)]
        
        for rxn in producing:
            print(f"  {rxn.id}: {rxn.reaction}")
            print(f"    가역성: {rxn.reversibility}")
        
        # MDH로 OAA 생산 테스트 (Malate에서)
        print(f"\nMDH를 통한 OAA 생산 테스트:")
        try:
            mdh = model.reactions.get_by_id('MDH')
            print(f"  MDH: {mdh.reaction}")
            print(f"  가역성: {mdh.reversibility}")
            
            # Malate가 생성 가능한지 먼저 확인
            mal_c = model.metabolites.get_by_id('mal__L_c')
            from cobra import Reaction
            
            test_mal = Reaction('TEST_mal')
            test_mal.lower_bound = 0
            test_mal.upper_bound = 1000
            test_mal.add_metabolites({mal_c: -1})
            model.add_reactions([test_mal])
            model.objective = 'TEST_mal'
            
            solution = model.optimize()
            if solution.status == 'optimal' and solution.objective_value > 1e-6:
                print(f"  [OK] Malate 생산 가능: {solution.objective_value:.6f}")
            else:
                print(f"  [FAIL] Malate 생산 불가 - 이것이 OAA 생산을 막고 있음!")
            
            model.remove_reactions([test_mal])
            
        except KeyError:
            print("  [ERROR] MDH 반응 없음")
            
    except KeyError:
        print("[ERROR] oaa_c 없음")

def check_bootstrap_problem(model):
    """Bootstrap 문제: 초기 대사물질이 없어서 경로가 시작되지 않는 문제"""
    print("\n" + "="*70)
    print("Bootstrap 문제 확인")
    print("="*70)
    
    model = setup_medium(model)
    
    # 핵심 순환 확인
    print("핵심 순환 확인:")
    print("  1. CS는 OAA가 필요함")
    print("  2. OAA는 MDH (Malate → OAA)로 생성 가능")
    print("  3. Malate는 MALS (Glyoxylate shunt)로 생성 가능")
    print("  4. Glyoxylate는 ICL (Isocitrate → Glyoxylate)로 생성")
    print("  5. Isocitrate는 ACONT (Citrate → Isocitrate)로 생성")
    print("  6. Citrate는 CS (Acetyl-CoA + OAA → Citrate)로 생성")
    print("\n  → 순환! 초기 OAA가 없으면 시작 불가")
    
    print("\n해결책 확인:")
    print("  - OAA를 생산하는 다른 경로가 있는지 확인")
    
    try:
        oaa_c = model.metabolites.get_by_id('oaa_c')
        # OAA 생산 반응 중 MDH 외 다른 것 확인
        producing = [r for r in oaa_c.reactions if oaa_c in r.products]
        
        print(f"\n  OAA를 직접 생산하는 반응 (비순환):")
        non_mdh = [r for r in producing if 'MDH' not in r.id and r.id != 'MDH']
        if non_mdh:
            for rxn in non_mdh:
                print(f"    {rxn.id}: {rxn.reaction}")
                # 이 반응들이 활성화 가능한지 확인
                # 필요한 대사물질이 생성 가능한지 확인
        else:
            print(f"    없음 - 이것이 문제!")
            print(f"\n  [ROOT CAUSE] OAA를 초기화할 수 있는 경로가 없습니다!")
            print(f"    → Pyruvate carboxylase나 다른 OAA 생산 경로가 필요할 수 있습니다")
            
    except KeyError:
        print("[ERROR] oaa_c 없음")

def check_pyruvate_carboxylase(model):
    """Pyruvate carboxylase (Pyruvate → OAA) 확인"""
    print("\n" + "="*70)
    print("Pyruvate Carboxylase 확인")
    print("="*70)
    
    # Pyruvate → OAA 반응 찾기
    try:
        pyr_c = model.metabolites.get_by_id('pyr_c')
        oaa_c = model.metabolites.get_by_id('oaa_c')
        
        # Pyruvate와 OAA를 모두 포함하는 반응
        pyr_oaa_rxns = [r for r in model.reactions 
                       if pyr_c in r.metabolites and oaa_c in r.metabolites]
        
        if pyr_oaa_rxns:
            print(f"[OK] Pyruvate → OAA 반응 발견: {len(pyr_oaa_rxns)}개")
            for rxn in pyr_oaa_rxns:
                print(f"  {rxn.id}: {rxn.reaction}")
                if pyr_c in rxn.reactants and oaa_c in rxn.products:
                    print(f"    [OK] Pyruvate → OAA 방향 가능")
        else:
            print(f"[FAIL] Pyruvate → OAA 반응 없음")
            print(f"  [ROOT CAUSE] Pyruvate carboxylase가 필요합니다!")
            
    except KeyError:
        print("[ERROR] pyr_c 또는 oaa_c 없음")

def main():
    model = load_model("BaseModel.xml")
    
    check_acs_requirements(model)
    check_oaa_cycle_problem(model)
    check_bootstrap_problem(model)
    check_pyruvate_carboxylase(model)
    
    print("\n" + "="*70)
    print("근본 원인 분석 완료")
    print("="*70)
    print("\n[결론]")
    print("  성장 실패의 주요 원인:")
    print("    1. OAA bootstrap 문제 - 초기 OAA 생성 경로 없음")
    print("    2. Pyruvate carboxylase (Pyruvate → OAA) 반응이 없을 가능성")
    print("    3. 또는 다른 OAA 초기화 경로 필요")

if __name__ == "__main__":
    main()

