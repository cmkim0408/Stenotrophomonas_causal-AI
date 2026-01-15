#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
maeB (NADP-dependent malic enzyme) 확인 및 추가

1. 모델에 maeB 관련 반응이 있는지 확인
2. PEPCK 제거
3. maeB 반응 추가 (없는 경우)
"""

import cobra
from pathlib import Path
import sys

def load_model(model_path):
    try:
        model = cobra.io.read_sbml_model(str(model_path))
        print(f"[OK] 모델 로드 완료 (반응 수: {len(model.reactions)})")
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

def check_maeb_reactions(model):
    """maeB 관련 반응 확인"""
    print("\n" + "="*70)
    print("maeB (NADP-dependent malic enzyme) 반응 확인")
    print("="*70)
    
    # maeB 관련 반응 ID 패턴
    maeB_patterns = ['MAE', 'MAEB', 'ME', 'MALME', 'maeB']
    
    found_rxns = []
    
    for rxn in model.reactions:
        rxn_id_upper = rxn.id.upper()
        rxn_name_upper = (rxn.name if rxn.name else '').upper()
        
        # 반응 ID나 이름에 maeB 관련 키워드가 있는지 확인
        if any(pattern in rxn_id_upper or pattern in rxn_name_upper for pattern in maeB_patterns):
            found_rxns.append(rxn)
        
        # 반응식에서 malate와 pyruvate, NADP가 모두 있는지 확인
        mets = [m.id for m in rxn.metabolites]
        if 'mal__L_c' in mets and 'pyr_c' in mets and 'nadp_c' in mets:
            # malate를 소비하고 pyruvate를 생성하는지 확인
            mal_coeff = rxn.metabolites.get(model.metabolites.get_by_id('mal__L_c'), 0)
            pyr_coeff = rxn.metabolites.get(model.metabolites.get_by_id('pyr_c'), 0)
            nadp_coeff = rxn.metabolites.get(model.metabolites.get_by_id('nadp_c'), 0)
            
            if mal_coeff < 0 and pyr_coeff > 0 and nadp_coeff < 0:  # malate 소비, pyruvate 생성, NADP 소비
                found_rxns.append(rxn)
    
    if found_rxns:
        print(f"\n[찾은 maeB 관련 반응] (총 {len(found_rxns)}개)")
        for rxn in found_rxns:
            print(f"\n  {rxn.id}")
            if rxn.name:
                print(f"    이름: {rxn.name}")
            print(f"    반응식: {rxn.reaction}")
            print(f"    bounds: [{rxn.lower_bound}, {rxn.upper_bound}]")
        return found_rxns[0] if found_rxns else None
    else:
        print(f"\n[결과] maeB 관련 반응이 모델에 없습니다.")
        return None

def check_pepck_reactions(model):
    """PEPCK 관련 반응 확인"""
    print("\n" + "="*70)
    print("PEPCK 관련 반응 확인")
    print("="*70)
    
    pepck_patterns = ['PEPCK', 'PCK', 'PYCK']
    
    found_rxns = []
    
    for rxn in model.reactions:
        rxn_id_upper = rxn.id.upper()
        if any(pattern in rxn_id_upper for pattern in pepck_patterns):
            found_rxns.append(rxn)
    
    if found_rxns:
        print(f"\n[찾은 PEPCK 관련 반응] (총 {len(found_rxns)}개)")
        for rxn in found_rxns:
            print(f"\n  {rxn.id}")
            if rxn.name:
                print(f"    이름: {rxn.name}")
            print(f"    반응식: {rxn.reaction}")
            print(f"    bounds: [{rxn.lower_bound}, {rxn.upper_bound}]")
        return found_rxns
    else:
        print(f"\n[결과] PEPCK 관련 반응이 모델에 없습니다.")
        return []

def add_maeb_reaction(model):
    """maeB 반응 추가"""
    print("\n" + "="*70)
    print("maeB 반응 추가")
    print("="*70)
    
    # maeB 반응식: malate + NADP+ → pyruvate + CO2 + NADPH
    # 또는: mal__L_c + nadp_c --> co2_c + nadph_c + pyr_c
    
    # 필요한 metabolite 확인
    required_mets = ['mal__L_c', 'nadp_c', 'co2_c', 'nadph_c', 'pyr_c']
    
    for met_id in required_mets:
        if met_id not in model.metabolites:
            print(f"[ERROR] {met_id}가 모델에 없습니다.")
            return False
    
    # maeB 반응 ID
    maeB_id = 'MAEB'
    
    if maeB_id in model.reactions:
        print(f"[INFO] {maeB_id} 반응이 이미 모델에 있습니다.")
        return True
    
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
    
    return True

def remove_pepck_reactions(model, pepck_rxns):
    """PEPCK 반응 제거"""
    print("\n" + "="*70)
    print("PEPCK 반응 제거")
    print("="*70)
    
    if not pepck_rxns:
        print("[INFO] 제거할 PEPCK 반응이 없습니다.")
        return
    
    rxn_ids = [rxn.id for rxn in pepck_rxns]
    print(f"\n[제거할 반응]")
    for rxn_id in rxn_ids:
        print(f"  {rxn_id}")
    
    model.remove_reactions(rxn_ids)
    print(f"\n[OK] {len(rxn_ids)}개 PEPCK 반응 제거 완료")

def test_model(model):
    """모델 테스트"""
    print("\n" + "="*70)
    print("모델 테스트")
    print("="*70)
    
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
    
    # PEPCK 플럭스 확인 (제거되었는지)
    pepck_rxns = ['PEPCK_ATP', 'PCK', 'PYCK']
    for rxn_id in pepck_rxns:
        if rxn_id in model.reactions:
            print(f"  [주의] {rxn_id}가 여전히 모델에 있습니다.")
    
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
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_fixed_acetate_recycling.xml"
    output_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_with_maeB_no_PEPCK.xml"
    
    print("="*70)
    print("maeB 확인 및 추가, PEPCK 제거")
    print("="*70)
    
    # 모델 로드
    model = load_model(model_path)
    model = setup_media_forced(model)
    
    # maeB 확인
    maeb_rxn = check_maeb_reactions(model)
    
    # PEPCK 확인
    pepck_rxns = check_pepck_reactions(model)
    
    # PEPCK 제거
    remove_pepck_reactions(model, pepck_rxns)
    
    # maeB 추가 (없는 경우)
    if maeb_rxn is None:
        success = add_maeb_reaction(model)
        if not success:
            print("[ERROR] maeB 반응 추가 실패")
            return
    else:
        print(f"\n[INFO] maeB 반응이 이미 모델에 있습니다: {maeb_rxn.id}")
    
    # 모델 테스트
    solution = test_model(model)
    
    # 모델 저장
    if solution.status == 'optimal' and solution.objective_value > 1e-6:
        cobra.io.write_sbml_model(model, str(output_path))
        print(f"\n[모델 저장] {output_path}")
    else:
        print(f"\n[주의] 모델이 성장하지 않습니다. 저장하지 않습니다.")
    
    return model, solution

if __name__ == "__main__":
    model, solution = main()
