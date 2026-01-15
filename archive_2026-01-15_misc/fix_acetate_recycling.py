#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Acetate 재순환 문제 해결

문제:
- ACt2rpp: h_p <=> ac_c + h_c (Acetate를 비현실적으로 재생성)
- AADb/AADa 경로가 비정상적으로 높은 플럭스
- 실제 Acetate uptake(19)와 Acetyl-CoA 생산이 맞지 않음

해결 방법:
1. ACt2rpp 반응 제거 또는 제한
2. ACS 반응 사용 확인
3. Acetate balance 확인
"""

import cobra
from pathlib import Path
import sys

def load_model(model_path):
    try:
        model = cobra.io.read_sbml_model(str(model_path))
        return model
    except Exception as e:
        print(f"[ERROR] 모델 로드 실패: {e}")
        sys.exit(1)

def setup_media_forced(model):
    """배지 조건을 강제로 고정"""
    if 'EX_ac_e' in model.reactions:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.lower_bound = -19.0
        ex_ac.upper_bound = -19.0
    
    if 'EX_o2_e' in model.reactions:
        ex_o2 = model.reactions.get_by_id('EX_o2_e')
        ex_o2.lower_bound = -100.0
        ex_o2.upper_bound = 1000.0
    
    essential_exchanges = {
        'EX_nh4_e': (-1000.0, 1000.0),
        'EX_pi_e': (-1000.0, 1000.0),
        'EX_so4_e': (-1000.0, 1000.0),
        'EX_mg2_e': (-1000.0, 1000.0),
        'EX_k_e': (-1000.0, 1000.0),
        'EX_na1_e': (-1000.0, 1000.0),
        'EX_fe2_e': (-1000.0, 1000.0),
        'EX_fe3_e': (-1000.0, 1000.0),
        'EX_h2o_e': (-1000.0, 1000.0),
        'EX_h_e': (-1000.0, 1000.0),
        'EX_co2_e': (-1000.0, 1000.0),
        'EX_hco3_e': (-1000.0, 1000.0),
        'EX_nac_e': (-1000.0, 1000.0),
        'EX_ncam_e': (-1000.0, 1000.0),
    }
    
    for ex_id, (lb, ub) in essential_exchanges.items():
        if ex_id in model.reactions:
            ex_rxn = model.reactions.get_by_id(ex_id)
            ex_rxn.lower_bound = lb
            ex_rxn.upper_bound = ub
    
    return model

def check_act2rpp_reaction(model, solution):
    """ACt2rpp 반응 확인"""
    print("\n" + "="*70)
    print("ACt2rpp 반응 확인")
    print("="*70)
    
    if 'ACt2rpp' not in model.reactions:
        print("[INFO] ACt2rpp 반응이 모델에 없습니다.")
        return None
    
    act2rpp = model.reactions.get_by_id('ACt2rpp')
    act2rpp_flux = solution.fluxes.get('ACt2rpp', 0.0)
    
    print(f"\n[ACt2rpp 반응]")
    print(f"  반응식: {act2rpp.reaction}")
    print(f"  플럭스: {act2rpp_flux:.6f}")
    print(f"  bounds: [{act2rpp.lower_bound}, {act2rpp.upper_bound}]")
    print(f"  설명: h_p <=> ac_c + h_c (Acetate를 재생성)")
    
    if act2rpp_flux > 1e-6:
        print(f"\n[문제] ACt2rpp가 Acetate를 {act2rpp_flux:.6f}만큼 생성하고 있습니다!")
        print(f"       이것이 Acetate 재순환의 원인입니다.")
    
    return act2rpp

def test_remove_act2rpp(model):
    """ACt2rpp 제거 테스트"""
    print("\n" + "="*70)
    print("ACt2rpp 제거 테스트")
    print("="*70)
    
    # ATPM=0 설정
    if 'ATPM' in model.reactions:
        atpm = model.reactions.get_by_id('ATPM')
        atpm.lower_bound = 0.0
        atpm.upper_bound = 0.0
    
    # 제거 전
    print("\n[제거 전]")
    solution_before = model.optimize()
    print(f"  상태: {solution_before.status}")
    print(f"  성장률: {solution_before.objective_value:.6f}")
    
    if 'ACt2rpp' in model.reactions:
        act2rpp_flux_before = solution_before.fluxes.get('ACt2rpp', 0.0)
        print(f"  ACt2rpp 플럭스: {act2rpp_flux_before:.6f}")
    
    if 'AADb' in model.reactions:
        aadb_flux_before = solution_before.fluxes.get('AADb', 0.0)
        print(f"  AADb 플럭스: {aadb_flux_before:.6f}")
    
    if 'ACS' in model.reactions:
        acs_flux_before = solution_before.fluxes.get('ACS', 0.0)
        print(f"  ACS 플럭스: {acs_flux_before:.6f}")
    
    # ACt2rpp 제거
    if 'ACt2rpp' in model.reactions:
        act2rpp = model.reactions.get_by_id('ACt2rpp')
        print(f"\n[ACt2rpp 제거]")
        print(f"  반응식: {act2rpp.reaction}")
        model.remove_reactions(['ACt2rpp'])
        print(f"  [OK] ACt2rpp 제거 완료")
    else:
        print(f"\n[INFO] ACt2rpp가 이미 모델에 없습니다.")
    
    # 제거 후
    print("\n[제거 후]")
    solution_after = model.optimize()
    print(f"  상태: {solution_after.status}")
    print(f"  성장률: {solution_after.objective_value:.6f}")
    
    if 'AADb' in model.reactions:
        aadb_flux_after = solution_after.fluxes.get('AADb', 0.0)
        print(f"  AADb 플럭스: {aadb_flux_after:.6f}")
    
    if 'ACS' in model.reactions:
        acs_flux_after = solution_after.fluxes.get('ACS', 0.0)
        print(f"  ACS 플럭스: {acs_flux_after:.6f}")
    
    # Acetate balance 확인
    if 'EX_ac_e' in model.reactions:
        ex_ac_flux = solution_after.fluxes.get('EX_ac_e', 0.0)
        print(f"  EX_ac_e 플럭스: {ex_ac_flux:.6f}")
    
    # Acetate → Acetyl-CoA 전환 확인
    print("\n[Acetate → Acetyl-CoA 전환 확인]")
    ac_to_accoa_rxns = ['ACS', 'AADb']
    total_ac_consumed = 0.0
    
    for rxn_id in ac_to_accoa_rxns:
        if rxn_id in model.reactions:
            rxn = model.reactions.get_by_id(rxn_id)
            flux = solution_after.fluxes.get(rxn_id, 0.0)
            if abs(flux) > 1e-6:
                if 'ac_c' in [m.id for m in rxn.metabolites]:
                    ac_coeff = rxn.metabolites.get(model.metabolites.get_by_id('ac_c'), 0)
                    if ac_coeff < 0:  # 소비
                        ac_consumed = abs(flux * ac_coeff)
                        total_ac_consumed += ac_consumed
                        print(f"  {rxn_id:10s}: 플럭스 {flux:10.6f}, Acetate 소비 {ac_consumed:.6f}")
    
    if 'EX_ac_e' in model.reactions:
        ex_ac_flux = solution_after.fluxes.get('EX_ac_e', 0.0)
        acetate_uptake = abs(ex_ac_flux)
        print(f"\n  Acetate uptake: {acetate_uptake:.6f}")
        print(f"  Acetate 소비 (Acetyl-CoA 생성을 위해): {total_ac_consumed:.6f}")
        
        if abs(acetate_uptake - total_ac_consumed) < 1e-3:
            print(f"  [OK] Acetate balance가 맞습니다!")
        else:
            print(f"  [주의] Acetate balance 차이: {abs(acetate_uptake - total_ac_consumed):.6f}")
    
    return solution_before, solution_after

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_without_ACS_ADP_SUCOAACTr.xml"
    output_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_fixed_acetate_recycling.xml"
    
    print("="*70)
    print("Acetate 재순환 문제 해결")
    print("="*70)
    
    # 모델 로드
    model = load_model(model_path)
    model = setup_media_forced(model)
    
    # ACt2rpp 확인
    solution = model.optimize()
    act2rpp = check_act2rpp_reaction(model, solution)
    
    # ACt2rpp 제거 테스트
    sol_before, sol_after = test_remove_act2rpp(model)
    
    # 결과 요약
    print("\n" + "="*70)
    print("결과 요약")
    print("="*70)
    
    growth_before = sol_before.objective_value
    growth_after = sol_after.objective_value
    
    print(f"\n[성장률 비교]")
    print(f"  ACt2rpp 있을 때: {growth_before:.6f}")
    print(f"  ACt2rpp 없을 때: {growth_after:.6f}")
    print(f"  차이: {growth_after - growth_before:.6f}")
    
    if growth_after > 1e-6:
        print(f"\n[결론] ACt2rpp 제거 후에도 성장 가능 (성장률: {growth_after:.6f})")
        
        # Acetate balance 확인
        if 'EX_ac_e' in model.reactions:
            ex_ac_flux = sol_after.fluxes.get('EX_ac_e', 0.0)
            acetate_uptake = abs(ex_ac_flux)
            
            # Acetate 소비량 계산
            total_ac_consumed = 0.0
            for rxn in model.reactions:
                if 'ac_c' in [m.id for m in rxn.metabolites]:
                    ac_coeff = rxn.metabolites.get(model.metabolites.get_by_id('ac_c'), 0)
                    if ac_coeff < 0:  # 소비
                        flux = sol_after.fluxes.get(rxn.id, 0.0)
                        if abs(flux) > 1e-6:
                            total_ac_consumed += abs(flux * ac_coeff)
            
            print(f"\n[Acetate Balance]")
            print(f"  Uptake: {acetate_uptake:.6f}")
            print(f"  소비량: {total_ac_consumed:.6f}")
            
            if abs(acetate_uptake - total_ac_consumed) < 1e-3:
                print(f"  [OK] Acetate balance가 맞습니다!")
            else:
                print(f"  [주의] Acetate balance 차이: {abs(acetate_uptake - total_ac_consumed):.6f}")
        
        # 모델 저장
        cobra.io.write_sbml_model(model, str(output_path))
        print(f"\n[모델 저장] {output_path}")
    else:
        print(f"\n[결론] ACt2rpp 제거 후 성장 불가 (성장률: {growth_after:.6f})")
        print(f"  -> ACt2rpp가 필수일 수 있습니다.")
    
    return model, sol_before, sol_after

if __name__ == "__main__":
    model, sol_before, sol_after = main()
