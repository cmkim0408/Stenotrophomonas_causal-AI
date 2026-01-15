#!/usr/bin/env python
"""
Acetate uptake가 작동하지 않는 정확한 원인 찾기
"""

import cobra
from cobra import Reaction

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    return model

def test_exchange_directly(model):
    """Exchange 반응 직접 테스트"""
    print("="*70)
    print("Exchange 반응 직접 테스트")
    print("="*70)
    
    # 모든 exchange 차단
    for rxn in model.exchanges:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    
    # EX_ac_e만 활성화
    try:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.lower_bound = -10
        ex_ac.upper_bound = 10
        
        print(f"[OK] EX_ac_e 설정:")
        print(f"  반응식: {ex_ac.reaction}")
        print(f"  하한: {ex_ac.lower_bound}, 상한: {ex_ac.upper_bound}")
        
        # ac_e 확인
        ac_e = model.metabolites.get_by_id('ac_e')
        print(f"\n  ac_e metabolite:")
        print(f"    ID: {ac_e.id}")
        print(f"    구획: {ac_e.compartment}")
        print(f"    연결된 반응: {len(ac_e.reactions)}개")
        for rxn in ac_e.reactions:
            coeff = rxn.metabolites[ac_e]
            print(f"      {rxn.id}: 계수={coeff}, {rxn.reaction}")
        
        # ac_e 생산 테스트
        test_ac_e = Reaction('TEST_ac_e')
        test_ac_e.lower_bound = 0
        test_ac_e.upper_bound = 1000
        test_ac_e.add_metabolites({ac_e: -1})
        model.add_reactions([test_ac_e])
        model.objective = 'TEST_ac_e'
        
        solution = model.optimize()
        print(f"\n  ac_e 생산 테스트:")
        print(f"    상태: {solution.status}")
        print(f"    Objective: {solution.objective_value:.6f}")
        print(f"    EX_ac_e 플럭스: {solution.fluxes.get('EX_ac_e', 0):.6f}")
        
        if solution.objective_value > 1e-6:
            print(f"    [OK] ac_e 생산 가능!")
        else:
            print(f"    [FAIL] ac_e 생산 불가")
            print(f"    [원인] EX_ac_e가 작동하지 않음")
        
        model.remove_reactions([test_ac_e])
        
    except KeyError as e:
        print(f"[ERROR] {e}")

def test_transport_directly(model):
    """Transport 반응 직접 테스트"""
    print("\n" + "="*70)
    print("Transport 반응 직접 테스트")
    print("="*70)
    
    # Medium 설정
    for rxn in model.exchanges:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    
    # EX_ac_e 활성화
    ex_ac = model.reactions.get_by_id('EX_ac_e')
    ex_ac.lower_bound = -10
    
    # ACt 활성화
    ac_transport = model.reactions.get_by_id('ACt')
    ac_transport.lower_bound = -1000
    ac_transport.upper_bound = 1000
    
    # ac_c 생산 테스트
    ac_c = model.metabolites.get_by_id('ac_c')
    test_ac_c = Reaction('TEST_ac_c')
    test_ac_c.lower_bound = 0
    test_ac_c.upper_bound = 1000
    test_ac_c.add_metabolites({ac_c: -1})
    model.add_reactions([test_ac_c])
    model.objective = 'TEST_ac_c'
    
    solution = model.optimize()
    
    print(f"ac_c 생산 테스트 (EX_ac_e + ACt 활성화):")
    print(f"  상태: {solution.status}")
    print(f"  Objective: {solution.objective_value:.6f}")
    print(f"  EX_ac_e 플럭스: {solution.fluxes.get('EX_ac_e', 0):.6f}")
    print(f"  ACt 플럭스: {solution.fluxes.get('ACt', 0):.6f}")
    
    if solution.objective_value > 1e-6:
        print(f"  [OK] ac_c 생산 가능!")
    else:
        print(f"  [FAIL] ac_c 생산 불가")
        print(f"  [원인] Transport 경로가 작동하지 않음")
    
    model.remove_reactions([test_ac_c])

def check_compartment_issue(model):
    """구획 문제 확인"""
    print("\n" + "="*70)
    print("구획 일관성 확인")
    print("="*70)
    
    try:
        ac_e = model.metabolites.get_by_id('ac_e')
        ac_c = model.metabolites.get_by_id('ac_c')
        
        print(f"ac_e:")
        print(f"  ID: {ac_e.id}")
        print(f"  구획: {ac_e.compartment}")
        print(f"  반응 수: {len(ac_e.reactions)}")
        
        print(f"\nac_c:")
        print(f"  ID: {ac_c.id}")
        print(f"  구획: {ac_c.compartment}")
        print(f"  반응 수: {len(ac_c.reactions)}")
        
        # ACt 확인
        ac_transport = model.reactions.get_by_id('ACt')
        print(f"\nACt transport:")
        print(f"  반응식: {ac_transport.reaction}")
        print(f"  대사물질:")
        for met, coeff in ac_transport.metabolites.items():
            print(f"    {met.id} ({met.compartment}): {coeff}")
        
        # 구획이 일치하는지 확인
        if ac_e.compartment == 'C_e' and ac_c.compartment == 'C_c':
            print(f"\n  [OK] 구획이 올바름 (ac_e: C_e, ac_c: C_c)")
        else:
            print(f"\n  [WARNING] 구획이 예상과 다름")
        
    except Exception as e:
        print(f"[ERROR] {e}")

def test_simple_loop(model):
    """간단한 순환 테스트: ac_e → ac_c → (소비)"""
    print("\n" + "="*70)
    print("간단한 순환 테스트")
    print("="*70)
    
    # Medium 설정
    for rxn in model.exchanges:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    
    ex_ac = model.reactions.get_by_id('EX_ac_e')
    ex_ac.lower_bound = -10
    
    # ac_c를 소비하는 반응 찾기
    ac_c = model.metabolites.get_by_id('ac_c')
    ac_consuming = [r for r in ac_c.reactions if ac_c in r.reactants]
    
    print(f"ac_c를 소비하는 반응: {len(ac_consuming)}개")
    print(f"주요 소비 반응:")
    for rxn in ac_consuming[:10]:
        print(f"  {rxn.id}: {rxn.reaction}")
    
    # ACS로 테스트
    try:
        acs = model.reactions.get_by_id('ACS')
        print(f"\nACS 반응으로 테스트:")
        print(f"  {acs.reaction}")
        
        # ACS에 필요한 것들 확인
        print(f"  필요한 대사물질:")
        for met, coeff in acs.metabolites.items():
            if coeff < 0:
                print(f"    {met.id}: {-coeff}")
                # 이 대사물질이 생성 가능한지
                if met.id == 'ac_c':
                    print(f"      → ac_c는 ACt를 통해 생성 가능해야 함")
                elif met.id == 'atp_c':
                    print(f"      → atp_c는 초기에 제공 필요")
                elif met.id == 'coa_c':
                    print(f"      → coa_c는 초기에 제공 필요 또는 다른 경로")
        
    except KeyError:
        print("[ERROR] ACS 없음")

def main():
    model = load_model("BaseModel.xml")
    
    test_exchange_directly(model)
    test_transport_directly(model)
    check_compartment_issue(model)
    test_simple_loop(model)
    
    print("\n" + "="*70)
    print("Acetate uptake 디버깅 완료")
    print("="*70)

if __name__ == "__main__":
    main()

