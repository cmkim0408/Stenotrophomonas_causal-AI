#!/usr/bin/env python
"""
성장 실패 원인 심층 분석
EX_ac_e와 ACt가 blocked인 이유 확인
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
    
    model.objective = 'Growth'
    return model

def check_metabolite_connectivity(model):
    """ac_e metabolite의 연결성 확인"""
    print("="*70)
    print("ac_e Metabolite 연결성 확인")
    print("="*70)
    
    try:
        ac_e = model.metabolites.get_by_id('ac_e')
        print(f"[OK] ac_e metabolite 존재")
        print(f"  ID: {ac_e.id}")
        print(f"  이름: {ac_e.name}")
        print(f"  구획: {ac_e.compartment}")
        print(f"  연결된 반응 수: {len(ac_e.reactions)}")
        
        print(f"\n  연결된 반응:")
        for rxn in ac_e.reactions:
            coeff = rxn.metabolites[ac_e]
            direction = "소비" if coeff < 0 else "생성"
            print(f"    {rxn.id}: {rxn.reaction} (계수: {coeff}, {direction})")
        
        # ac_e를 소비하는 반응이 있는지
        consuming = [rxn for rxn in ac_e.reactions if rxn.metabolites[ac_e] < 0]
        producing = [rxn for rxn in ac_e.reactions if rxn.metabolites[ac_e] > 0]
        
        print(f"\n  ac_e를 소비하는 반응: {len(consuming)}개")
        print(f"  ac_e를 생성하는 반응: {len(producing)}개")
        
        if len(consuming) == 0:
            print(f"  [WARNING] ac_e를 소비하는 반응이 없습니다!")
        
    except KeyError:
        print("[ERROR] ac_e metabolite가 없습니다!")

def check_transport_reaction(model):
    """ACt transport reaction 상세 확인"""
    print("\n" + "="*70)
    print("ACt Transport Reaction 상세 분석")
    print("="*70)
    
    try:
        ac_transport = model.reactions.get_by_id('ACt')
        print(f"[OK] ACt 반응 존재")
        print(f"  반응식: {ac_transport.reaction}")
        print(f"  하한: {ac_transport.lower_bound}")
        print(f"  상한: {ac_transport.upper_bound}")
        print(f"  가역성: {ac_transport.reversibility}")
        
        # 대사물질 확인
        print(f"\n  반응에 참여하는 대사물질:")
        for met, coeff in ac_transport.metabolites.items():
            print(f"    {met.id} ({met.compartment}): {coeff}")
        
        # ac_e와 ac_c가 모두 있는지 확인
        ac_e = model.metabolites.get_by_id('ac_e')
        ac_c = model.metabolites.get_by_id('ac_c')
        
        has_ac_e = ac_e in ac_transport.metabolites
        has_ac_c = ac_c in ac_transport.metabolites
        
        print(f"\n  ac_e 포함: {has_ac_e}")
        print(f"  ac_c 포함: {has_ac_c}")
        
        if not has_ac_e or not has_ac_c:
            print(f"  [ERROR] Transport 반응에 ac_e 또는 ac_c가 없습니다!")
        
    except KeyError:
        print("[ERROR] ACt transport 반응이 없습니다!")

def test_individual_reactions(model):
    """개별 반응 테스트"""
    print("\n" + "="*70)
    print("개별 반응 활성화 테스트")
    print("="*70)
    
    # 각 반응을 개별적으로 활성화해보기
    test_reactions = ['EX_ac_e', 'ACt', 'ACS']
    
    for rxn_id in test_reactions:
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            print(f"\n{rxn_id} 테스트:")
            
            # 원래 경계 저장
            orig_lb = rxn.lower_bound
            orig_ub = rxn.upper_bound
            
            # 강제로 활성화 시도
            if 'EX' in rxn_id:
                rxn.lower_bound = -10  # Uptake 허용
                rxn.upper_bound = 10
            else:
                rxn.lower_bound = 0
                rxn.upper_bound = 10
            
            # FBA 실행
            solution = model.optimize()
            
            if solution.status == 'optimal':
                flux = solution.fluxes.get(rxn_id, 0)
                print(f"  최적화 상태: {solution.status}")
                print(f"  {rxn_id} 플럭스: {flux:.6f}")
                print(f"  Objective: {solution.objective_value:.6f}")
                
                if abs(flux) > 1e-6:
                    print(f"  [OK] 활성화 가능")
                else:
                    print(f"  [FAIL] 활성화 불가 (플럭스 0)")
            else:
                print(f"  [FAIL] 최적화 실패: {solution.status}")
            
            # 원래 경계 복원
            rxn.lower_bound = orig_lb
            rxn.upper_bound = orig_ub
            
        except KeyError:
            print(f"\n{rxn_id}: [NOT FOUND]")
        except Exception as e:
            print(f"\n{rxn_id}: [ERROR] {e}")

def check_mass_balance(model):
    """질량 균형 확인"""
    print("\n" + "="*70)
    print("질량 균형 확인 (ac_e)")
    print("="*70)
    
    try:
        ac_e = model.metabolites.get_by_id('ac_e')
        
        # ac_e의 생산/소비 균형
        production = 0
        consumption = 0
        
        for rxn in ac_e.reactions:
            coeff = rxn.metabolites[ac_e]
            if coeff > 0:
                production += abs(coeff)
                print(f"  생산: {rxn.id} (계수: {coeff})")
            elif coeff < 0:
                consumption += abs(coeff)
                print(f"  소비: {rxn.id} (계수: {coeff})")
        
        print(f"\n  총 생산 계수: {production}")
        print(f"  총 소비 계수: {consumption}")
        
        if production == 0 and consumption > 0:
            print(f"  [WARNING] ac_e를 생산하는 경로가 없습니다!")
        if consumption == 0 and production > 0:
            print(f"  [WARNING] ac_e를 소비하는 경로가 없습니다!")
            
    except KeyError:
        print("[ERROR] ac_e metabolite가 없습니다!")

def find_blocking_constraints(model):
    """막히는 제약 조건 찾기"""
    print("\n" + "="*70)
    print("Blocking 제약 조건 찾기")
    print("="*70)
    
    # EX_ac_e를 강제로 활성화
    try:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.lower_bound = -10
        
        # ACt도 활성화
        ac_transport = model.reactions.get_by_id('ACt')
        ac_transport.lower_bound = -10
        ac_transport.upper_bound = 10
        
        # ACS도 활성화
        acs = model.reactions.get_by_id('ACS')
        acs.upper_bound = 10
        
        # FBA 실행
        solution = model.optimize()
        
        print(f"강제 활성화 후 FBA 결과:")
        print(f"  상태: {solution.status}")
        print(f"  Objective: {solution.objective_value:.6f}")
        
        if solution.status == 'optimal':
            print(f"\n  주요 반응 플럭스:")
            key_rxns = ['EX_ac_e', 'ACt', 'ACS', 'CS', 'ICL', 'MALS', 'Growth']
            for rxn_id in key_rxns:
                try:
                    flux = solution.fluxes.get(rxn_id, 0)
                    if abs(flux) > 1e-8:
                        print(f"    {rxn_id}: {flux:.6f}")
                except:
                    pass
        else:
            print(f"  [ERROR] 최적화 실패 - 제약 조건 충돌 가능")
            
    except Exception as e:
        print(f"[ERROR] 테스트 중 오류: {e}")
        import traceback
        traceback.print_exc()

def main():
    model = load_model("BaseModel.xml")
    model = setup_medium(model)
    
    check_metabolite_connectivity(model)
    check_transport_reaction(model)
    test_individual_reactions(model)
    check_mass_balance(model)
    find_blocking_constraints(model)
    
    print("\n" + "="*70)
    print("심층 진단 완료")
    print("="*70)

if __name__ == "__main__":
    main()

