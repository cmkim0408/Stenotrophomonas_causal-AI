#!/usr/bin/env python
"""
최종 아세트산 성장 테스트 - 정확한 반응 ID와 bounds 확인
"""

import cobra

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    return model

def main():
    model = load_model("BaseModel.xml")
    
    print("="*70)
    print("모델에서 Acetate 관련 반응 확인")
    print("="*70)
    
    # 모든 반응에서 acetate 찾기
    acetate_rxns = [r for r in model.reactions if 'ac' in r.id.lower() or 'acetate' in r.name.lower()]
    print(f"\nAcetate 관련 반응: {len(acetate_rxns)}개")
    for rxn in acetate_rxns[:10]:
        print(f"  {rxn.id}: {rxn.name}")
        print(f"    반응식: {rxn.reaction}")
        print(f"    하한: {rxn.lower_bound}, 상한: {rxn.upper_bound}")
        print(f"    Exchange? {rxn in model.exchanges}")
    
    # EX_ac_e 정확히 찾기
    print("\n" + "="*70)
    print("Exchange 반응 확인")
    print("="*70)
    
    ex_rxns = [r for r in model.exchanges if 'ac' in r.id.lower()]
    print(f"Acetate exchange 반응: {len(ex_rxns)}개")
    for rxn in ex_rxns:
        print(f"  ID: {rxn.id}")
        print(f"  이름: {rxn.name}")
        print(f"  반응식: {rxn.reaction}")
        print(f"  하한: {rxn.lower_bound}, 상한: {rxn.upper_bound}")
    
    # Medium 설정
    print("\n" + "="*70)
    print("Acetate Medium 설정 및 성장 테스트")
    print("="*70)
    
    # 모든 exchange 차단
    for rxn in model.exchanges:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    
    # Acetate exchange 찾아서 활성화
    ex_ac = None
    for rxn in model.exchanges:
        if 'ac' in rxn.id.lower() and ('ex' in rxn.id.lower() or rxn.id.startswith('EX_')):
            ex_ac = rxn
            break
    
    if ex_ac:
        ex_ac.lower_bound = -1000
        ex_ac.upper_bound = 1000
        print(f"[OK] {ex_ac.id} 활성화")
    else:
        print("[ERROR] Acetate exchange 반응을 찾을 수 없음")
        return
    
    # 필수 영양소 활성화
    essential = ['EX_nh4_e', 'EX_h2o_e', 'EX_h_e', 'EX_pi_e', 'EX_so4_e', 
                'EX_k_e', 'EX_mg2_e', 'EX_fe2_e', 'EX_mn2_e', 'EX_zn2_e',
                'EX_co2_e', 'EX_o2_e']
    
    for ex_id in essential:
        try:
            ex = model.reactions.get_by_id(ex_id)
            if ex_id in ['EX_co2_e', 'EX_o2_e']:
                ex.lower_bound = -1000
                ex.upper_bound = 1000
            else:
                ex.lower_bound = -1000
            print(f"[OK] {ex_id} 활성화")
        except KeyError:
            print(f"[WARN] {ex_id} 없음")
    
    # Biomass objective
    try:
        model.objective = 'Growth'
        print(f"\n[OK] Objective: Growth")
    except:
        print("[ERROR] Growth 반응 없음")
        return
    
    # FBA
    print("\nFBA 실행...")
    solution = model.optimize()
    
    print(f"\n결과:")
    print(f"  상태: {solution.status}")
    print(f"  Biomass: {solution.objective_value:.6f}")
    
    if solution.objective_value > 1e-6:
        print(f"  [SUCCESS] 성장 가능!")
    else:
        print(f"  [FAIL] 성장 불가")
        print(f"\n  주요 반응 플럭스:")
        for rxn_id in ['EX_ac_e', 'R_EX_ac_e', 'ACt', 'R_ACt', 'ACS', 'R_ACS', 'CS', 'Growth']:
            try:
                flux = solution.fluxes.get(rxn_id, 0)
                if abs(flux) > 1e-8:
                    print(f"    {rxn_id}: {flux:.6f}")
            except:
                pass
        
        # 모든 활성 반응 확인
        print(f"\n  활성 exchange 반응:")
        for rxn in model.exchanges:
            flux = solution.fluxes.get(rxn.id, 0)
            if abs(flux) > 1e-6:
                print(f"    {rxn.id}: {flux:.4f}")

if __name__ == "__main__":
    main()

