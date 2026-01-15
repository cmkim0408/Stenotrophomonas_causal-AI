#!/usr/bin/env python
"""
9개 누락된 반응을 신규 모델에 추가하고 FBA 테스트
실제 사용된 9개 반응을 추가하여 FBA 성공 여부 확인
"""

import cobra
from pathlib import Path
import pandas as pd
from cobra.flux_analysis import find_blocked_reactions

def load_model(model_path):
    """모델 로드"""
    print(f"모델 로드 중: {model_path}")
    model = cobra.io.read_sbml_model(str(model_path))
    print(f"[OK] 모델 로드 완료: {model.id}")
    print(f"  반응 수: {len(model.reactions)}, 대사물질 수: {len(model.metabolites)}")
    return model

def add_reaction_to_model(model, ref_model, rxn_id):
    """레퍼런스 모델에서 반응을 신규 모델에 추가"""
    if rxn_id in model.reactions:
        print(f"  [SKIP] {rxn_id}: 이미 모델에 존재")
        return False
    
    if rxn_id not in ref_model.reactions:
        print(f"  [ERROR] {rxn_id}: 레퍼런스 모델에 없음")
        return False
    
    ref_rxn = ref_model.reactions.get_by_id(rxn_id)
    
    # 대사물질 먼저 확인/추가
    metabolites_dict = {}
    for met in ref_rxn.metabolites:
        if met.id not in model.metabolites:
            new_met = cobra.Metabolite(met.id, formula=met.formula, name=met.name, 
                                      compartment=met.compartment, charge=met.charge)
            model.add_metabolites([new_met])
        metabolites_dict[model.metabolites.get_by_id(met.id)] = ref_rxn.metabolites[met]
    
    # 반응 생성 및 추가
    new_rxn = cobra.Reaction(rxn_id)
    new_rxn.name = ref_rxn.name
    new_rxn.lower_bound = ref_rxn.lower_bound
    new_rxn.upper_bound = ref_rxn.upper_bound
    new_rxn.add_metabolites(metabolites_dict)
    
    model.add_reactions([new_rxn])
    
    # 유전자 연결
    for gene in ref_rxn.genes:
        if gene.id in model.genes:
            model.genes.get_by_id(gene.id).reactions.add(new_rxn)
    
    print(f"  [ADDED] {rxn_id}: {ref_rxn.reaction}")
    return True

def setup_acetate_medium(model):
    """Acetate 기반 미디어 설정 (레퍼런스와 유사하게)"""
    print("\n미디어 설정 중...")
    
    exchanges = {
        'EX_ac_e': (-10, 1000),      # Acetate uptake
        'EX_o2_e': (-1000, 1000),    # Oxygen
        'EX_h2o_e': (-1000, 1000),   # Water
        'EX_h_e': (-1000, 1000),     # H+
        'EX_nh4_e': (-10, 1000),     # Ammonium
        'EX_pi_e': (-10, 1000),      # Phosphate
        'EX_so4_e': (-10, 1000),     # Sulfate
        'EX_k_e': (-1000, 1000),     # Potassium
        'EX_na1_e': (-1000, 1000),   # Sodium
        'EX_mg2_e': (-1000, 1000),   # Magnesium
        'EX_ca2_e': (-1000, 1000),   # Calcium
        'EX_fe2_e': (-10, 1000),     # Fe2+
        'EX_fe3_e': (-10, 1000),     # Fe3+
        'EX_hco3_e': (-1000, 1000),  # Bicarbonate
        'EX_co2_e': (-1000, 1000),   # CO2
    }
    
    applied = 0
    for ex_id, bounds in exchanges.items():
        if ex_id in model.reactions:
            model.reactions.get_by_id(ex_id).bounds = bounds
            applied += 1
        else:
            print(f"  [WARN] {ex_id}: 모델에 없음")
    
    print(f"  {applied}개 exchange 반응 설정 완료")
    return model

def test_fba(model):
    """FBA 테스트 및 진단"""
    print("\n" + "="*70)
    print("FBA 테스트")
    print("="*70)
    
    # Biomass 반응 찾기
    biomass_rxns = [r for r in model.reactions if 'growth' in r.id.lower() or 'biomass' in r.id.lower()]
    if not biomass_rxns:
        print("[ERROR] Biomass 반응을 찾을 수 없습니다")
        return None
    
    biomass_rxn = biomass_rxns[0]
    model.objective = biomass_rxn.id
    print(f"목적함수: {biomass_rxn.id}")
    
    # FBA 실행
    try:
        solution = model.optimize()
        print(f"\nFBA 상태: {solution.status}")
        print(f"목적함수 값 (성장률): {solution.objective_value}")
        
        if solution.status == 'optimal' and solution.objective_value > 1e-6:
            print("\n" + "="*70)
            print("[SUCCESS] FBA 성공! 모델이 성장할 수 있습니다.")
            print("="*70)
            
            # 추가된 9개 반응의 플럭스 확인
            print("\n추가된 9개 반응의 플럭스:")
            key_reactions = ['ACS_ADP', 'PEPCK_ATP', 'SUCDi', 'ACtexi', 
                            'EX_hco3_e', 'T_hco3_e_to_c', 'T_o2_e_to_o2_c', 
                            'T_nh4_e_to_nh4_c', 'T_fe3_e_to_fe3_c']
            for rxn_id in key_reactions:
                if rxn_id in model.reactions:
                    flux = solution.fluxes.get(rxn_id, 0)
                    print(f"  {rxn_id:20s}: {flux:>12.6f}")
                else:
                    print(f"  {rxn_id:20s}: [누락]")
            
            # Exchange 플럭스 확인
            print("\n주요 Exchange 플럭스:")
            ex_reactions = ['EX_ac_e', 'EX_o2_e', 'EX_hco3_e', 'EX_nh4_e', 'EX_co2_e']
            for ex_id in ex_reactions:
                if ex_id in model.reactions:
                    flux = solution.fluxes.get(ex_id, 0)
                    if abs(flux) > 1e-6:
                        print(f"  {ex_id:20s}: {flux:>12.6f}")
            
            return solution
        else:
            print("\n" + "="*70)
            print("[FAIL] FBA 실패 또는 성장 불가능")
            print("="*70)
            
            # 진단 정보
            if solution.status == 'optimal' and solution.objective_value <= 1e-6:
                print("\n진단: 최적해는 있지만 성장률이 0입니다.")
                
                # Blocked reactions 확인
                try:
                    blocked = find_blocked_reactions(model)
                    print(f"  Blocked reactions 수: {len(blocked)}")
                except:
                    pass
                
                # 추가된 반응들의 플럭스 확인
                print("\n추가된 9개 반응의 플럭스 (0일 수 있음):")
                key_reactions = ['ACS_ADP', 'PEPCK_ATP', 'SUCDi', 'ACtexi', 
                                'EX_hco3_e', 'T_hco3_e_to_c', 'T_o2_e_to_o2_c', 
                                'T_nh4_e_to_nh4_c', 'T_fe3_e_to_fe3_c']
                for rxn_id in key_reactions:
                    if rxn_id in model.reactions:
                        flux = solution.fluxes.get(rxn_id, 0)
                        print(f"  {rxn_id:20s}: {flux:>12.6f}")
            
            return solution
            
    except Exception as e:
        print(f"\n[ERROR] FBA 실행 중 오류: {e}")
        import traceback
        traceback.print_exc()
        return None

def main():
    # 경로 설정
    base_path = Path(__file__).parent.parent
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    new_model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    
    print("="*70)
    print("9개 누락된 반응 추가 및 FBA 테스트")
    print("="*70)
    
    # 모델 로드
    ref_model = load_model(ref_model_path)
    new_model = load_model(new_model_path)
    
    # 누락된 9개 반응 리스트 (실제 사용된 반응)
    missing_reactions = [
        'ACS_ADP',
        'PEPCK_ATP',
        'SUCDi',
        'ACtexi',
        'EX_hco3_e',
        'T_hco3_e_to_c',
        'T_o2_e_to_o2_c',
        'T_nh4_e_to_nh4_c',
        'T_fe3_e_to_fe3_c'
    ]
    
    print("\n" + "="*70)
    print("누락된 9개 반응 추가")
    print("="*70)
    
    added_count = 0
    for rxn_id in missing_reactions:
        if add_reaction_to_model(new_model, ref_model, rxn_id):
            added_count += 1
    
    print(f"\n총 {added_count}/{len(missing_reactions)}개 반응 추가 완료")
    
    # 미디어 설정
    new_model = setup_acetate_medium(new_model)
    
    # FBA 테스트
    solution = test_fba(new_model)
    
    # 최종 결과 요약
    print("\n" + "="*70)
    print("최종 결과 요약")
    print("="*70)
    print(f"추가된 반응 수: {added_count}/{len(missing_reactions)}")
    
    if solution and solution.status == 'optimal' and solution.objective_value > 1e-6:
        print(f"FBA 상태: SUCCESS")
        print(f"성장률: {solution.objective_value:.6f} 1/h")
        print("\n결론: 9개 반응 추가로 FBA가 성공했습니다!")
        print("-> 이 9개 반응의 유전자를 신규 분리 미생물에서 확인하세요.")
    else:
        print("FBA 상태: FAIL")
        if solution:
            print(f"성장률: {solution.objective_value:.6f} 1/h")
        else:
            print("성장률: N/A")
        print("\n결론: 9개 반응만으로는 부족합니다.")
        print("-> 추가 반응이 필요할 수 있습니다.")
        print("-> 대사물질 연결 문제일 수도 있습니다.")
    
    print("\n" + "="*70)

if __name__ == "__main__":
    main()
