#!/usr/bin/env python
"""
BaseModel_with_ACtexi.xml로 FBA 실행
"""

import cobra
from pathlib import Path

def load_model(model_path):
    model = cobra.io.read_sbml_model(model_path)
    return model

def setup_acetate_medium_growing(model):
    """성장했을 때의 미디어 설정"""
    exchanges = {
        'EX_ac_e': (-10, 1000),
        'EX_o2_e': (-1000, 1000),
        'EX_h2o_e': (-1000, 1000),
        'EX_h_e': (-1000, 1000),
        'EX_nh4_e': (-10, 1000),
        'EX_pi_e': (-10, 1000),
        'EX_so4_e': (-10, 1000),
        'EX_k_e': (-1000, 1000),
        'EX_na1_e': (-1000, 1000),
        'EX_mg2_e': (-1000, 1000),
        'EX_ca2_e': (-1000, 1000),
        'EX_fe2_e': (-10, 1000),
        'EX_fe3_e': (-10, 1000),
        'EX_hco3_e': (-1000, 1000),
        'EX_co2_e': (-1000, 1000),
    }
    
    for ex_id, bounds in exchanges.items():
        if ex_id in model.reactions:
            model.reactions.get_by_id(ex_id).bounds = bounds
    
    return model

def run_fba(model):
    """FBA 실행"""
    print("="*80)
    print("FBA 실행")
    print("="*80)
    
    # ATPM 설정 확인
    atpm_rxn = model.reactions.get_by_id('ATPM')
    print(f"\n[ATPM 설정]")
    print(f"  bounds: [{atpm_rxn.lower_bound}, {atpm_rxn.upper_bound}]")
    
    # ATPM=0으로 설정
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    model.objective = 'Growth'
    solution = model.optimize()
    
    print(f"\n[FBA 결과]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    
    # 주요 반응 플럭스
    key_reactions = [
        'ACS_ADP', 'ACS', 'SUCDi', 'PEPCK_ATP', 'ACtexi',
        'CS', 'ICL', 'MALS', 
        'ATPS4rpp', 'NADH16pp', 'CYTBO3_4pp',
        'ADK1', 'EX_ac_e', 'EX_o2_e'
    ]
    
    print(f"\n[주요 반응 플럭스]")
    active_fluxes = {}
    for rxn_id in key_reactions:
        if rxn_id in model.reactions:
            flux = solution.fluxes.get(rxn_id, 0.0)
            if abs(flux) > 1e-6:
                print(f"  {rxn_id:20s}: {flux:>12.6f}")
                active_fluxes[rxn_id] = flux
    
    if len(active_fluxes) == 0:
        print(f"  [문제] 모든 주요 반응 플럭스가 0")
    else:
        print(f"\n[작동하는 반응] {len(active_fluxes)}개")
    
    return solution, active_fluxes

def test_multiple_atpm_values(model):
    """여러 ATPM 값에서 테스트"""
    print("\n" + "="*80)
    print("여러 ATPM 값에서 테스트")
    print("="*80)
    
    atpm_values = [0, 5, 10, 20]
    
    print(f"\n{'ATPM':<10} {'성장률':<15} {'상태':<15} {'CS':<12} {'ACS_ADP':<12} {'SUCDi':<12}")
    print("-" * 75)
    
    for atpm_val in atpm_values:
        model_test = model.copy()
        atpm_rxn = model_test.reactions.get_by_id('ATPM')
        atpm_rxn.lower_bound = atpm_val
        atpm_rxn.upper_bound = 1000
        
        model_test.objective = 'Growth'
        solution = model_test.optimize()
        
        growth = solution.objective_value if solution.status == 'optimal' else 0.0
        cs_flux = solution.fluxes.get('CS', 0.0) if solution.status == 'optimal' else 0.0
        acs_adp_flux = solution.fluxes.get('ACS_ADP', 0.0) if solution.status == 'optimal' else 0.0
        sucdi_flux = solution.fluxes.get('SUCDi', 0.0) if solution.status == 'optimal' else 0.0
        
        print(f"{atpm_val:<10.1f} {growth:<15.6f} {solution.status:<15} {cs_flux:<12.6f} {acs_adp_flux:<12.6f} {sucdi_flux:<12.6f}")

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_with_ACtexi.xml"
    
    print("="*80)
    print("BaseModel_with_ACtexi.xml로 FBA 실행")
    print("="*80)
    
    # 모델 로드
    model = load_model(str(model_path))
    
    print(f"\n[모델 정보]")
    print(f"  반응 수: {len(model.reactions)}")
    print(f"  메타볼라이트 수: {len(model.metabolites)}")
    
    # 추가된 반응 확인
    key_reactions = ['ACS_ADP', 'SUCDi', 'PEPCK_ATP', 'ACtexi']
    print(f"\n[추가된 반응 확인]")
    for rxn_id in key_reactions:
        if rxn_id in model.reactions:
            print(f"  {rxn_id}: 존재")
        else:
            print(f"  {rxn_id}: 없음")
    
    # 미디어 설정
    model = setup_acetate_medium_growing(model)
    
    # FBA 실행
    solution, active_fluxes = run_fba(model)
    
    # 여러 ATPM 값에서 테스트
    test_multiple_atpm_values(model)
    
    print("\n" + "="*80)
    print("최종 결과")
    print("="*80)
    
    if solution.status == 'optimal' and solution.objective_value > 1e-6:
        print(f"\n[성공] FBA 성공! 성장률: {solution.objective_value:.6f}")
    elif len(active_fluxes) > 0:
        print(f"\n[부분 성공] 경로는 작동하지만 성장률이 0")
        print(f"  -> {len(active_fluxes)}개 반응이 작동함")
    else:
        print(f"\n[실패] 경로가 작동하지 않음")
        print(f"  -> 추가 조사 필요")

if __name__ == "__main__":
    main()
