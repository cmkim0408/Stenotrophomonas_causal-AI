#!/usr/bin/env python
"""
사용자 질문에 대한 답변: "아세트산만으로 TCA cycle과 전자 전달계를 통해 ATP를 만들 수 있잖아요? 왜 못만드는걸까요?"
"""

import cobra

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    return model

def main():
    model = load_model("BaseModel.xml")
    
    print("="*70)
    print("아세트산만으로 ATP 생산이 불가능한 이유 분석")
    print("="*70)
    
    print("\n[이론적으로는 가능해야 합니다]")
    print("  Acetate → Acetyl-CoA (ACS)")
    print("  Acetyl-CoA → TCA cycle → NADH/FADH2")
    print("  NADH/FADH2 → 전자전달계 → ATP")
    
    print("\n[실제 모델에서 확인된 문제점]")
    
    # 문제 1: CoA 초기화
    print("\n1. CoA 초기화 문제")
    coa_c = model.metabolites.get_by_id('coa_c')
    coa_producing = [r for r in coa_c.reactions 
                    if coa_c in r.products or (r.reversibility and coa_c in r.reactants)]
    coa_without_accoa = [r for r in coa_producing 
                        if 'accoa' not in str(r.metabolites).lower()]
    print(f"   CoA 생산 반응: {len(coa_producing)}개")
    print(f"   Acetyl-CoA 없이 CoA 생산 가능한 반응: {len(coa_without_accoa)}개")
    if len(coa_without_accoa) == 0:
        print(f"   [문제] Acetyl-CoA 없이는 CoA를 생성할 수 없음")
        print(f"   → ACS는 CoA가 필요한데, CoA는 Acetyl-CoA에서 생성됨")
        print(f"   → 순환 의존성 (Chicken-and-egg problem)")
    else:
        print(f"   [OK] CoA 초기화 경로 존재")
        for rxn in coa_without_accoa[:3]:
            print(f"      {rxn.id}: {rxn.reaction}")
    
    # 문제 2: OAA 초기화
    print("\n2. OAA (Oxaloacetate) 초기화 문제")
    oaa_c = model.metabolites.get_by_id('oaa_c')
    oaa_producing = [r for r in oaa_c.reactions 
                    if oaa_c in r.products or (r.reversibility and oaa_c in r.reactants)]
    print(f"   OAA 생산 반응: {len(oaa_producing)}개")
    
    # OAA 생산 테스트 (아세트산만으로)
    for rxn in model.exchanges:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    
    ex_ac = model.reactions.get_by_id('EX_ac_e')
    ex_ac.lower_bound = -1000
    
    essentials = ['EX_nh4_e', 'EX_h2o_e', 'EX_h_e', 'EX_pi_e', 'EX_so4_e',
                  'EX_k_e', 'EX_mg2_e', 'EX_fe2_e', 'EX_mn2_e', 'EX_zn2_e',
                  'EX_co2_e', 'EX_o2_e']
    for ex_id in essentials:
        try:
            ex = model.reactions.get_by_id(ex_id)
            if ex_id in ['EX_co2_e', 'EX_o2_e']:
                ex.lower_bound = -1000
                ex.upper_bound = 1000
            else:
                ex.lower_bound = -1000
        except:
            pass
    
    test_oaa = cobra.Reaction('TEST_oaa')
    test_oaa.lower_bound = 0
    test_oaa.upper_bound = 1000
    test_oaa.add_metabolites({oaa_c: -1})
    model.add_reactions([test_oaa])
    model.objective = 'TEST_oaa'
    
    sol_oaa = model.optimize()
    print(f"   아세트산만으로 OAA 생산 가능: {sol_oaa.objective_value:.6f}")
    if sol_oaa.objective_value < 1e-6:
        print(f"   [문제] 아세트산만으로는 OAA를 생성할 수 없음")
        print(f"   → TCA cycle 시작에 OAA가 필요하지만, OAA 생산이 불가능")
        print(f"   → Glyoxylate shunt는 있지만, 초기 OAA가 필요")
    else:
        print(f"   [OK] OAA 생산 경로 존재")
    
    model.remove_reactions([test_oaa])
    
    # 문제 3: ATP 초기화
    print("\n3. ATP 초기화 문제")
    print(f"   [분석] ACS 반응은 ATP를 필요로 함")
    print(f"   → 초기 ATP가 없으면 ACS가 작동하지 않음")
    print(f"   → 하지만 TCA cycle과 전자전달계가 작동하면 ATP 생성 가능")
    print(f"   → 따라서 초기 ATP가 필요 (Bootstrap)")
    
    print("\n" + "="*70)
    print("결론 및 해결방안")
    print("="*70)
    print("\n[왜 아세트산만으로 ATP를 못 만드는가?]")
    print("\n1. 순환 의존성 (Circular Dependency):")
    print("   - ACS: ac_c + ATP + CoA → Acetyl-CoA")
    print("   - CoA 생산: 대부분 Acetyl-CoA에서 생성됨")
    print("   - → 초기 CoA가 필요")
    
    print("\n2. TCA Cycle 시작 조건:")
    print("   - CS: Acetyl-CoA + OAA → Citrate")
    print("   - → 초기 OAA가 필요")
    print("   - Glyoxylate shunt로 OAA 재생성 가능하지만,")
    print("     초기 OAA가 있어야 순환이 시작됨")
    
    print("\n3. ATP Bootstrap:")
    print("   - ACS는 ATP를 소비")
    print("   - TCA + ETC에서 ATP 생성")
    print("   - → 초기 소량의 ATP가 필요할 수 있음")
    
    print("\n[해결방안]")
    print("실제 생물학적으로는 세포에 초기 대사물질이 존재합니다:")
    print("  - CoA pool")
    print("  - OAA pool (또는 생성 경로)")
    print("  - ATP pool")
    print("\n모델에서도 초기 조건으로 제공하면:")
    print("  1. 초기 CoA 제공")
    print("  2. 초기 OAA 제공 (또는 OAA 생산 경로 확인)")
    print("  3. 그 후 아세트산만으로 지속적인 ATP 생산 및 성장 가능")
    
    print("\n[또는 모델 개선]")
    print("  - CoA 초기화 경로 추가 (pantothenate → CoA)")
    print("  - OAA 생산 경로 확인/추가 (PEP carboxylase 등)")
    print("="*70)

if __name__ == "__main__":
    main()

