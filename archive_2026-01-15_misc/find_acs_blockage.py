#!/usr/bin/env python
"""
ACS가 작동하지 않는 정확한 원인 찾기
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

def test_acs_with_bootstrap(model):
    """ACS 작동 테스트: 필요한 것들을 하나씩 제공"""
    print("="*70)
    print("ACS 작동 조건 테스트")
    print("="*70)
    
    model = setup_medium(model)
    
    # Test 1: ATP만 제공
    print("\n[Test 1] ATP만 제공 (0.1)")
    atp_c = model.metabolites.get_by_id('atp_c')
    atp_boot = Reaction('EX_atp')
    atp_boot.lower_bound = -0.1
    atp_boot.add_metabolites({atp_c: -1})
    model.add_reactions([atp_boot])
    
    # ACS 활성화 테스트
    acs = model.reactions.get_by_id('ACS')
    test_accoa = Reaction('TEST_accoa')
    test_accoa.add_metabolites({model.metabolites.get_by_id('accoa_c'): -1})
    test_accoa.lower_bound = 0
    test_accoa.upper_bound = 1000
    model.add_reactions([test_accoa])
    model.objective = 'TEST_accoa'
    
    sol1 = model.optimize()
    print(f"  상태: {sol1.status}, Acetyl-CoA: {sol1.objective_value:.6f}")
    print(f"  ACS 플럭스: {sol1.fluxes.get('ACS', 0):.6f}")
    print(f"  ac_c 플럭스: {sol1.fluxes.get('ACt', 0):.6f}")
    print(f"  atp_c 플럭스: {sol1.fluxes.get('EX_atp', 0):.6f}")
    print(f"  coa_c 관련 플럭스:")
    coa_rxns = [r for r in model.reactions if 'coa_c' in str(r.metabolites)]
    for rxn in coa_rxns[:5]:
        flux = sol1.fluxes.get(rxn.id, 0)
        if abs(flux) > 1e-6:
            print(f"    {rxn.id}: {flux:.4f}")
    
    model.remove_reactions([atp_boot, test_accoa])
    
    # Test 2: ATP + CoA 제공
    print("\n[Test 2] ATP + CoA 제공")
    atp_boot = Reaction('EX_atp')
    atp_boot.lower_bound = -0.1
    atp_boot.add_metabolites({atp_c: -1})
    
    coa_c = model.metabolites.get_by_id('coa_c')
    coa_boot = Reaction('EX_coa')
    coa_boot.lower_bound = -0.1
    coa_boot.add_metabolites({coa_c: -1})
    
    model.add_reactions([atp_boot, coa_boot])
    
    test_accoa = Reaction('TEST_accoa')
    test_accoa.add_metabolites({model.metabolites.get_by_id('accoa_c'): -1})
    test_accoa.lower_bound = 0
    test_accoa.upper_bound = 1000
    model.add_reactions([test_accoa])
    model.objective = 'TEST_accoa'
    
    sol2 = model.optimize()
    print(f"  상태: {sol2.status}, Acetyl-CoA: {sol2.objective_value:.6f}")
    print(f"  ACS 플럭스: {sol2.fluxes.get('ACS', 0):.6f}")
    
    if sol2.objective_value > 0.1:
        print(f"  [SUCCESS] ATP + CoA 제공으로 ACS 작동 가능!")
    
    model.remove_reactions([atp_boot, coa_boot, test_accoa])
    
    # Test 3: ATP + CoA + OAA 제공 (완전한 bootstrap)
    print("\n[Test 3] ATP + CoA + OAA 제공 (완전 bootstrap)")
    atp_boot = Reaction('EX_atp')
    atp_boot.lower_bound = -0.1
    atp_boot.add_metabolites({atp_c: -1})
    
    coa_boot = Reaction('EX_coa')
    coa_boot.lower_bound = -0.1
    coa_boot.add_metabolites({coa_c: -1})
    
    oaa_c = model.metabolites.get_by_id('oaa_c')
    oaa_boot = Reaction('EX_oaa')
    oaa_boot.lower_bound = -0.1
    oaa_boot.add_metabolites({oaa_c: -1})
    
    model.add_reactions([atp_boot, coa_boot, oaa_boot])
    
    model.objective = 'Growth'
    
    sol3 = model.optimize()
    print(f"  상태: {sol3.status}, Biomass: {sol3.objective_value:.6f}")
    
    if sol3.objective_value > 1e-6:
        print(f"  [SUCCESS] 완전 bootstrap으로 성장 가능!")
        print(f"  주요 플럭스:")
        for rxn_id in ['EX_ac_e', 'ACt', 'ACS', 'CS', 'ICL', 'MALS', 'Growth']:
            flux = sol3.fluxes.get(rxn_id, 0)
            if abs(flux) > 1e-6:
                print(f"    {rxn_id}: {flux:.4f}")
    else:
        print(f"  [FAIL] 완전 bootstrap으로도 성장 불가")
    
    model.remove_reactions([atp_boot, coa_boot, oaa_boot])

def check_coa_requirement(model):
    """CoA 요구사항 확인"""
    print("\n" + "="*70)
    print("CoA 생산 경로 상세 분석")
    print("="*70)
    
    model = setup_medium(model)
    
    try:
        coa_c = model.metabolites.get_by_id('coa_c')
        
        # CoA를 생성하는 반응 중 Acetyl-CoA를 사용하는 것들
        accoa_producing_coa = []
        for rxn in coa_c.reactions:
            if 'accoa' in str(rxn.metabolites).lower():
                accoa_producing_coa.append(rxn)
        
        print(f"Acetyl-CoA를 사용하여 CoA를 생성하는 반응: {len(accoa_producing_coa)}개")
        for rxn in accoa_producing_coa[:10]:
            print(f"  {rxn.id}: {rxn.reaction}")
        
        # CoA를 직접 생성하는 반응 (Acetyl-CoA 없이)
        direct_coa = [r for r in coa_c.reactions 
                     if 'accoa' not in str(r.metabolites).lower() 
                     and coa_c in r.products]
        
        print(f"\nAcetyl-CoA 없이 CoA를 생성하는 반응: {len(direct_coa)}개")
        for rxn in direct_coa[:10]:
            print(f"  {rxn.id}: {rxn.reaction}")
        
        if len(direct_coa) == 0:
            print(f"  [ROOT CAUSE] CoA를 초기화할 수 있는 경로가 없습니다!")
            print(f"    → CoA는 Acetyl-CoA에서만 생성되지만,")
            print(f"    → Acetyl-CoA 생산에 CoA가 필요 (순환 의존성)")
        
    except KeyError:
        print("[ERROR] coa_c 없음")

def main():
    model = load_model("BaseModel.xml")
    
    test_acs_with_bootstrap(model)
    check_coa_requirement(model)
    
    print("\n" + "="*70)
    print("ACS 차단 원인 분석 완료")
    print("="*70)

if __name__ == "__main__":
    main()

