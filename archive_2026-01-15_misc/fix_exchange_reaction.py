#!/usr/bin/env python
"""
Exchange reaction 문제 수정
ac_e가 boundary metabolite인지 확인 및 수정
"""

import cobra

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    return model

def check_boundary_metabolites(model):
    """Boundary metabolites 확인"""
    print("="*70)
    print("Boundary Metabolites 확인")
    print("="*70)
    
    # Exchange reaction과 연결된 metabolites 확인
    boundary_mets = set()
    for rxn in model.exchanges:
        for met in rxn.metabolites:
            boundary_mets.add(met.id)
    
    print(f"Exchange reaction과 연결된 대사물질: {len(boundary_mets)}개")
    
    # ac_e 확인
    try:
        ac_e = model.metabolites.get_by_id('ac_e')
        if ac_e.id in boundary_mets:
            print(f"[OK] ac_e는 boundary metabolite입니다")
        else:
            print(f"[WARNING] ac_e가 boundary metabolite 목록에 없습니다")
    except KeyError:
        print(f"[ERROR] ac_e metabolite가 없습니다")
    
    # 다른 예시 확인 (glucose)
    try:
        glc_e = model.metabolites.get_by_id('glc__D_e')
        if glc_e.id in boundary_mets:
            print(f"[OK] glc__D_e는 boundary metabolite입니다 (참고)")
    except KeyError:
        pass

def verify_exchange_reaction(model):
    """EX_ac_e 반응 검증"""
    print("\n" + "="*70)
    print("EX_ac_e 반응 검증")
    print("="*70)
    
    try:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        print(f"[OK] EX_ac_e 존재")
        print(f"  반응식: {ex_ac.reaction}")
        print(f"  대사물질:")
        for met, coeff in ex_ac.metabolites.items():
            print(f"    {met.id} ({met.compartment}): {coeff}")
        
        # Exchange reaction은 보통 하나의 metabolite만 있어야 함
        if len(ex_ac.metabolites) == 1:
            met = list(ex_ac.metabolites.keys())[0]
            coeff = ex_ac.metabolites[met]
            print(f"\n  [OK] Exchange 반응 형식 정상 (1개 대사물질)")
            print(f"    대사물질: {met.id}, 계수: {coeff}")
            
            # 계수가 -1이어야 uptake 가능
            if coeff == -1:
                print(f"    [OK] 계수가 -1입니다 (uptake 가능)")
            else:
                print(f"    [WARNING] 계수가 {coeff}입니다. -1이어야 합니다!")
        else:
            print(f"\n  [WARNING] Exchange 반응에 대사물질이 {len(ex_ac.metabolites)}개 있습니다")
            
    except KeyError:
        print("[ERROR] EX_ac_e 반응이 없습니다!")

def test_simple_uptake(model):
    """단순 uptake 테스트"""
    print("\n" + "="*70)
    print("단순 Uptake 테스트")
    print("="*70)
    
    # Medium 설정
    for rxn in model.exchanges:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    
    # EX_ac_e만 활성화
    try:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.lower_bound = -10  # Uptake 허용
        ex_ac.upper_bound = 10
        
        # ACt 활성화
        ac_transport = model.reactions.get_by_id('ACt')
        ac_transport.lower_bound = -1000
        ac_transport.upper_bound = 1000
        
        # 필수 영양소도 활성화
        for ex_id in ['EX_h2o_e', 'EX_h_e', 'EX_pi_e', 'EX_o2_e']:
            try:
                ex_rxn = model.reactions.get_by_id(ex_id)
                if ex_id == 'EX_o2_e':
                    ex_rxn.lower_bound = -1000
                    ex_rxn.upper_bound = 1000
                else:
                    ex_rxn.lower_bound = -1000
            except KeyError:
                pass
        
        # Objective를 ac_c 생산으로 변경 (간단한 테스트)
        print("테스트: ac_c 생산 가능 여부 확인")
        
        # 임시 반응 생성: ac_c -> (테스트용)
        try:
            test_rxn = model.reactions.get_by_id('TEST_ac_c')
        except KeyError:
            from cobra import Reaction
            ac_c = model.metabolites.get_by_id('ac_c')
            test_rxn = Reaction('TEST_ac_c')
            test_rxn.name = 'Test ac_c production'
            test_rxn.lower_bound = 0
            test_rxn.upper_bound = 1000
            test_rxn.add_metabolites({ac_c: -1})
            model.add_reactions([test_rxn])
            print(f"[OK] 테스트 반응 생성: {test_rxn.reaction}")
        
        model.objective = 'TEST_ac_c'
        
        solution = model.optimize()
        
        print(f"\n결과:")
        print(f"  상태: {solution.status}")
        print(f"  Objective (ac_c 생산): {solution.objective_value:.6f}")
        
        if solution.status == 'optimal':
            ex_ac_flux = solution.fluxes.get('EX_ac_e', 0)
            ac_transport_flux = solution.fluxes.get('ACt', 0)
            test_flux = solution.fluxes.get('TEST_ac_c', 0)
            
            print(f"  EX_ac_e 플럭스: {ex_ac_flux:.6f}")
            print(f"  ACt 플럭스: {ac_transport_flux:.6f}")
            print(f"  TEST_ac_c 플럭스: {test_flux:.6f}")
            
            if abs(ex_ac_flux) > 1e-6:
                print(f"  [OK] EX_ac_e 활성화 가능!")
            else:
                print(f"  [FAIL] EX_ac_e 활성화 불가")
            
            if abs(test_flux) > 1e-6:
                print(f"  [OK] ac_c 생산 가능!")
            else:
                print(f"  [FAIL] ac_c 생산 불가")
        
        # 테스트 반응 제거
        model.remove_reactions([test_rxn])
        
    except Exception as e:
        print(f"[ERROR] 테스트 중 오류: {e}")
        import traceback
        traceback.print_exc()

def check_compartment_consistency(model):
    """구획 일관성 확인"""
    print("\n" + "="*70)
    print("구획 일관성 확인")
    print("="*70)
    
    try:
        ac_e = model.metabolites.get_by_id('ac_e')
        ac_c = model.metabolites.get_by_id('ac_c')
        
        print(f"ac_e 구획: {ac_e.compartment}")
        print(f"ac_c 구획: {ac_c.compartment}")
        
        # EX_ac_e가 ac_e와 연결되어 있는지
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        if ac_e in ex_ac.metabolites:
            print(f"[OK] EX_ac_e가 ac_e와 연결됨")
        else:
            print(f"[ERROR] EX_ac_e가 ac_e와 연결되지 않음!")
        
        # ACt가 ac_e와 ac_c를 모두 연결하는지
        ac_transport = model.reactions.get_by_id('ACt')
        if ac_e in ac_transport.metabolites and ac_c in ac_transport.metabolites:
            print(f"[OK] ACt가 ac_e와 ac_c를 연결함")
        else:
            print(f"[ERROR] ACt가 ac_e와 ac_c를 연결하지 않음!")
            
    except Exception as e:
        print(f"[ERROR] 확인 중 오류: {e}")

def main():
    model = load_model("BaseModel.xml")
    
    check_boundary_metabolites(model)
    verify_exchange_reaction(model)
    check_compartment_consistency(model)
    test_simple_uptake(model)
    
    print("\n" + "="*70)
    print("검증 완료")
    print("="*70)

if __name__ == "__main__":
    main()

