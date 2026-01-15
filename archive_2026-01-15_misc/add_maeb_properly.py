#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
maeB 반응 명시적으로 추가

maeB: NADP-dependent malic enzyme
반응식: mal__L_c + nadp_c --> co2_c + nadph_c + pyr_c
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

def add_maeb_reaction(model):
    """maeB 반응 추가"""
    print("\n" + "="*70)
    print("maeB 반응 추가")
    print("="*70)
    
    # maeB 반응 ID
    maeB_id = 'MAEB'
    
    # 이미 있는지 확인
    if maeB_id in model.reactions:
        rxn = model.reactions.get_by_id(maeB_id)
        print(f"[INFO] {maeB_id} 반응이 이미 모델에 있습니다.")
        print(f"  반응식: {rxn.reaction}")
        return True
    
    # 필요한 metabolite 확인
    required_mets = {
        'mal__L_c': 'L-Malate',
        'nadp_c': 'NADP+',
        'co2_c': 'CO2',
        'nadph_c': 'NADPH',
        'pyr_c': 'Pyruvate',
    }
    
    missing_mets = []
    for met_id, met_name in required_mets.items():
        if met_id not in model.metabolites:
            missing_mets.append(f"{met_id} ({met_name})")
    
    if missing_mets:
        print(f"[ERROR] 필요한 metabolite가 모델에 없습니다:")
        for met in missing_mets:
            print(f"  {met}")
        return False
    
    # maeB 반응 생성
    maeB_rxn = cobra.Reaction(maeB_id)
    maeB_rxn.name = "NADP-dependent malic enzyme (maeB)"
    
    # 반응식: mal__L_c + nadp_c --> co2_c + nadph_c + pyr_c
    maeB_rxn.add_metabolites({
        model.metabolites.get_by_id('mal__L_c'): -1,
        model.metabolites.get_by_id('nadp_c'): -1,
        model.metabolites.get_by_id('co2_c'): 1,
        model.metabolites.get_by_id('nadph_c'): 1,
        model.metabolites.get_by_id('pyr_c'): 1,
    })
    
    maeB_rxn.lower_bound = -1000.0  # 가역 반응
    maeB_rxn.upper_bound = 1000.0
    
    model.add_reactions([maeB_rxn])
    
    print(f"[OK] {maeB_id} 반응 추가 완료")
    print(f"  반응식: {maeB_rxn.reaction}")
    print(f"  bounds: [{maeB_rxn.lower_bound}, {maeB_rxn.upper_bound}]")
    print(f"  설명: malate + NADP+ → pyruvate + CO2 + NADPH")
    print(f"  EC: 1.1.1.40")
    print(f"  Gene: maeB (Smlt3940)")
    
    return True

def test_model(model):
    """모델 테스트"""
    print("\n" + "="*70)
    print("모델 테스트")
    print("="*70)
    
    # 배지 조건 설정
    if 'EX_ac_e' in model.reactions:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.lower_bound = -19.0
        ex_ac.upper_bound = -19.0
    
    if 'EX_o2_e' in model.reactions:
        ex_o2 = model.reactions.get_by_id('EX_o2_e')
        ex_o2.lower_bound = -100.0
        ex_o2.upper_bound = 1000.0
    
    # ATPM=0 설정
    if 'ATPM' in model.reactions:
        atpm = model.reactions.get_by_id('ATPM')
        atpm.lower_bound = 0.0
        atpm.upper_bound = 0.0
    
    solution = model.optimize()
    
    print(f"\n[FBA 결과]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    
    # maeB 플럭스 확인
    if 'MAEB' in model.reactions:
        maeb_flux = solution.fluxes.get('MAEB', 0.0)
        print(f"  MAEB 플럭스: {maeb_flux:.6f}")
        if abs(maeb_flux) > 1e-6:
            print(f"  [OK] MAEB가 사용되고 있습니다.")
        else:
            print(f"  [INFO] MAEB 플럭스가 0입니다 (사용되지 않음)")
    
    # PEPCK 확인
    pepck_rxns = ['PEPCK_ATP', 'PCK', 'PYCK']
    pepck_found = []
    for rxn_id in pepck_rxns:
        if rxn_id in model.reactions:
            pepck_found.append(rxn_id)
    
    if pepck_found:
        print(f"\n[주의] PEPCK 반응이 여전히 모델에 있습니다: {pepck_found}")
    else:
        print(f"\n[OK] PEPCK 반응이 제거되었습니다.")
    
    # 주요 경로 플럭스
    print(f"\n[주요 경로 플럭스]")
    key_rxns = ['EX_ac_e', 'ACS', 'CS', 'ICL', 'MALS', 'MDH', 'MAEB', 'Growth']
    for rxn_id in key_rxns:
        if rxn_id in model.reactions:
            flux = solution.fluxes.get(rxn_id, 0.0)
            print(f"  {rxn_id:15s}: {flux:10.6f}")
    
    return solution

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_with_maeB_no_PEPCK.xml"
    output_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_with_maeB_no_PEPCK_final.xml"
    
    print("="*70)
    print("maeB 반응 명시적으로 추가")
    print("="*70)
    
    # 모델 로드
    model = load_model(model_path)
    
    # maeB 추가
    success = add_maeb_reaction(model)
    
    if not success:
        print("[ERROR] maeB 반응 추가 실패")
        return
    
    # 모델 테스트
    solution = test_model(model)
    
    # 모델 저장
    if solution.status == 'optimal':
        cobra.io.write_sbml_model(model, str(output_path))
        print(f"\n[모델 저장] {output_path}")
    else:
        print(f"\n[주의] FBA 최적화 실패. 저장하지 않습니다.")
    
    return model, solution

if __name__ == "__main__":
    model, solution = main()
