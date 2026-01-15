#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 4-1: 배지 조건을 "강제로 고정"해서 성장 가능성부터 판정

작업:
1. Acetate uptake를 고정 (EX_ac_e = -19)
2. O₂ uptake도 제한/고정 (EX_o2_e.lower_bound = -100)
3. 나머지 필수 무기물(NH4, Pi, SO4, Mg, K, Na 등)은 충분히 열어두기
4. 이 상태에서 maximize biomass 했는데도 μ=0이면, 진짜로 biomass가 막혀 있는 것
"""

import cobra
from pathlib import Path
import sys

def setup_media_forced(model):
    """배지 조건을 강제로 고정"""
    print("\n" + "="*70)
    print("Step 4-1: 배지 조건 강제 고정")
    print("="*70)
    
    # 1. Acetate uptake 고정
    if 'EX_ac_e' in model.reactions:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.lower_bound = -19.0
        ex_ac.upper_bound = -19.0
        print(f"[OK] EX_ac_e 고정: [{ex_ac.lower_bound}, {ex_ac.upper_bound}]")
    else:
        print("[WARNING] EX_ac_e 반응이 없습니다.")
    
    # 2. O₂ uptake 제한/고정
    if 'EX_o2_e' in model.reactions:
        ex_o2 = model.reactions.get_by_id('EX_o2_e')
        ex_o2.lower_bound = -100.0
        ex_o2.upper_bound = 1000.0
        print(f"[OK] EX_o2_e 제한: [{ex_o2.lower_bound}, {ex_o2.upper_bound}]")
    else:
        print("[WARNING] EX_o2_e 반응이 없습니다.")
    
    # 3. 필수 무기물 열어두기
    essential_exchanges = {
        'EX_nh4_e': (-1000.0, 1000.0),  # NH4
        'EX_pi_e': (-1000.0, 1000.0),   # Pi
        'EX_so4_e': (-1000.0, 1000.0),  # SO4
        'EX_mg2_e': (-1000.0, 1000.0),  # Mg
        'EX_k_e': (-1000.0, 1000.0),    # K
        'EX_na1_e': (-1000.0, 1000.0),  # Na
        'EX_fe2_e': (-1000.0, 1000.0),  # Fe2
        'EX_fe3_e': (-1000.0, 1000.0),  # Fe3
        'EX_h2o_e': (-1000.0, 1000.0),  # H2O
        'EX_h_e': (-1000.0, 1000.0),    # H+
        'EX_co2_e': (-1000.0, 1000.0),  # CO2
        'EX_hco3_e': (-1000.0, 1000.0), # HCO3
    }
    
    opened_count = 0
    for ex_id, (lb, ub) in essential_exchanges.items():
        if ex_id in model.reactions:
            ex_rxn = model.reactions.get_by_id(ex_id)
            ex_rxn.lower_bound = lb
            ex_rxn.upper_bound = ub
            opened_count += 1
        else:
            print(f"[WARNING] {ex_id} 반응이 없습니다.")
    
    print(f"[OK] 필수 무기물 {opened_count}개 열어둠")
    
    return model

def find_biomass_reaction(model):
    """Biomass 반응 찾기"""
    biomass_keywords = ['biomass', 'growth', 'BIOMASS', 'Growth']
    for rxn in model.reactions:
        if any(keyword in rxn.id for keyword in biomass_keywords):
            return rxn
    return None

def test_biomass_growth(model):
    """Biomass 성장 가능성 테스트"""
    print("\n" + "="*70)
    print("Biomass 성장 가능성 테스트")
    print("="*70)
    
    # Biomass 반응 찾기
    biomass_rxn = find_biomass_reaction(model)
    if biomass_rxn is None:
        print("[ERROR] Biomass 반응을 찾을 수 없습니다.")
        return None
    
    print(f"\n[Biomass 반응] {biomass_rxn.id}")
    print(f"  반응식: {biomass_rxn.reaction}")
    
    # Biomass를 objective로 설정
    model.objective = biomass_rxn.id
    model.objective_direction = 'max'
    
    # FBA 실행
    solution = model.optimize()
    
    print(f"\n[FBA 결과]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    
    if solution.status == 'optimal':
        if solution.objective_value > 1e-6:
            print(f"\n[결론] [OK] 성장 가능함 (성장률: {solution.objective_value:.6f})")
            return True
        else:
            print(f"\n[결론] [NO] 성장 불가 (성장률: 0)")
            print("  -> 진짜로 biomass가 막혀 있는 것")
            return False
    else:
        print(f"\n[결론] [ERROR] FBA 최적화 실패 (상태: {solution.status})")
        return False

def check_key_fluxes(model, solution):
    """주요 경로 플럭스 확인"""
    print("\n[주요 경로 플럭스]")
    key_reactions = ['ACS_ADP', 'CS', 'SUCDi', 'ICL', 'MALS', 'PEPCK_ATP']
    
    for rxn_id in key_reactions:
        if rxn_id in model.reactions:
            flux = solution.fluxes.get(rxn_id, 0.0)
            if abs(flux) > 1e-6:
                print(f"  {rxn_id:12s}: {flux:10.6f}")
            else:
                print(f"  {rxn_id:12s}: {flux:10.6f} (0)")

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_with_BCAA_cofactors_ions_nad_transport.xml"
    
    print("="*70)
    print("Step 4-1: 배지 조건 강제 고정 및 성장 가능성 판정")
    print("="*70)
    print(f"\n모델: {model_path}")
    
    # 모델 로드
    if not model_path.exists():
        print(f"[ERROR] 모델 파일이 없습니다: {model_path}")
        return
    
    model = cobra.io.read_sbml_model(str(model_path))
    print(f"[OK] 모델 로드 완료 (반응 수: {len(model.reactions)})")
    
    # 배지 조건 강제 고정
    model = setup_media_forced(model)
    
    # Biomass 성장 가능성 테스트
    can_grow = test_biomass_growth(model)
    
    if can_grow is False:
        solution = model.optimize()
        check_key_fluxes(model, solution)
        print("\n" + "="*70)
        print("다음 단계: Step 4-2 - Biomass가 요구하는 성분 중 뭐가 안 만들어지는지 찾기")
        print("="*70)
    
    return model, can_grow

if __name__ == "__main__":
    model, can_grow = main()
