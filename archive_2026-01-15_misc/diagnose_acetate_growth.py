#!/usr/bin/env python
"""
Acetate 성장 실패 원인 진단
"""

import cobra

def load_model(model_path="BaseModel.xml"):
    """모델 로드"""
    model = cobra.io.read_sbml_model(model_path)
    print(f"[OK] 모델 로드 완료: {model.id}\n")
    return model

def diagnose_acetate_pathway(model):
    """Acetate 경로 진단"""
    print("="*70)
    print("Acetate 경로 진단")
    print("="*70)
    
    # Medium 설정
    for rxn in model.exchanges:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    
    # Acetate 및 필수 영양소만 허용
    essentials = ['EX_ac_e', 'EX_nh4_e', 'EX_h2o_e', 'EX_h_e', 'EX_pi_e', 
                  'EX_so4_e', 'EX_k_e', 'EX_mg2_e', 'EX_fe2_e', 'EX_mn2_e', 
                  'EX_zn2_e', 'EX_co2_e', 'EX_o2_e']
    
    for ex_id in essentials:
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            if ex_id in ['EX_co2_e', 'EX_o2_e']:
                ex_rxn.lower_bound = -1000
                ex_rxn.upper_bound = 1000
            else:
                ex_rxn.lower_bound = -1000
        except KeyError:
            pass
    
    # Biomass objective
    model.objective = 'Growth'
    
    # FBA 실행
    print("\nFBA 최적화 중...")
    solution = model.optimize()
    
    print(f"\n최적화 상태: {solution.status}")
    print(f"Objective value: {solution.objective_value:.6f}")
    
    # 주요 반응 플럭스 확인
    print("\n주요 반응 플럭스:")
    key_reactions = ['EX_ac_e', 'ACt', 'ACS', 'CS', 'ICL', 'MALS', 'ME1', 'ME2', 'Growth']
    
    for rxn_id in key_reactions:
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            flux = solution.fluxes.get(rxn_id, 0)
            if abs(flux) > 1e-8:
                print(f"  {rxn_id}: {flux:.6f} | {rxn.reaction}")
            else:
                print(f"  {rxn_id}: {flux:.6f} (inactive)")
        except KeyError:
            print(f"  {rxn_id}: [NOT FOUND]")
    
    # Exchange 플럭스 확인
    print("\nExchange 플럭스:")
    for rxn in model.exchanges:
        flux = solution.fluxes.get(rxn.id, 0)
        if abs(flux) > 1e-6:
            print(f"  {rxn.id}: {flux:.4f}")
    
    return solution

def main():
    model = load_model("BaseModel.xml")
    solution = diagnose_acetate_pathway(model)
    
    if solution.objective_value > 0:
        print("\n[SUCCESS] Acetate에서 성장 가능!")
    else:
        print("\n[FAIL] Acetate에서 성장 불가")
        print("\n다음 단계:")
        print("  1. Gap analysis 수행 (필요한 반응 확인)")
        print("  2. Biomass reaction 구성 요소 확인")
        print("  3. 필수 대사 경로 완성도 확인")

if __name__ == "__main__":
    main()


