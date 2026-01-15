#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
레독스 셔틀 및 PPase 문제 수정

문제점:
1. SUCD + ACOAD2/ACOAD2f: FAD/FADH2를 자유 대사체로 사용하는 레독스 셔틀
2. PPA_1pp: H⁺-translocating PPase가 에너지 회수 아티팩트 생성 가능
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

def check_redox_shuttle(model, solution):
    """레독스 셔틀 확인"""
    print("\n" + "="*70)
    print("레독스 셔틀 확인 (SUCD + ACOAD2/ACOAD2f)")
    print("="*70)
    
    # SUCD 확인
    if 'SUCD' in model.reactions:
        sucd = model.reactions.get_by_id('SUCD')
        sucd_flux = solution.fluxes.get('SUCD', 0.0)
        print(f"\n[SUCD]")
        print(f"  반응식: {sucd.reaction}")
        print(f"  플럭스: {sucd_flux:.6f}")
        print(f"  bounds: [{sucd.lower_bound}, {sucd.upper_bound}]")
        print(f"  [문제] free FAD/FADH2를 사용 (Q8 기반이 아님)")
    else:
        print(f"\n[INFO] SUCD가 모델에 없습니다.")
    
    # ACOAD2, ACOAD2f 확인
    for rxn_id in ['ACOAD2', 'ACOAD2f']:
        if rxn_id in model.reactions:
            rxn = model.reactions.get_by_id(rxn_id)
            flux = solution.fluxes.get(rxn_id, 0.0)
            print(f"\n[{rxn_id}]")
            print(f"  반응식: {rxn.reaction}")
            print(f"  플럭스: {flux:.6f}")
            print(f"  bounds: [{rxn.lower_bound}, {rxn.upper_bound}]")
    
    # 루프 확인
    if all(rxn_id in model.reactions for rxn_id in ['SUCD', 'ACOAD2', 'ACOAD2f']):
        sucd_flux = solution.fluxes.get('SUCD', 0.0)
        acoad2_flux = solution.fluxes.get('ACOAD2', 0.0)
        acoad2f_flux = solution.fluxes.get('ACOAD2f', 0.0)
        
        print(f"\n[레독스 셔틀 확인]")
        print(f"  SUCD: {sucd_flux:.6f}")
        print(f"  ACOAD2: {acoad2_flux:.6f}")
        print(f"  ACOAD2f: {acoad2f_flux:.6f}")
        
        if abs(sucd_flux) > 1e-6 and abs(acoad2_flux) > 1e-6 and abs(acoad2f_flux) > 1e-6:
            if abs(abs(sucd_flux) - abs(acoad2_flux)) < 1e-3 and abs(abs(sucd_flux) - abs(acoad2f_flux)) < 1e-3:
                print(f"  [문제] 세 반응이 거의 같은 크기로 맞물려 레독스 셔틀을 형성합니다!")
                print(f"         FADH2 → NADH 변환으로 ATP 수율이 비현실적으로 좋아질 수 있습니다.")

def check_ppase(model, solution):
    """PPase 확인"""
    print("\n" + "="*70)
    print("PPase 확인 (PPA_1pp)")
    print("="*70)
    
    if 'PPA_1pp' in model.reactions:
        ppa = model.reactions.get_by_id('PPA_1pp')
        ppa_flux = solution.fluxes.get('PPA_1pp', 0.0)
        print(f"\n[PPA_1pp]")
        print(f"  반응식: {ppa.reaction}")
        print(f"  플럭스: {ppa_flux:.6f}")
        print(f"  bounds: [{ppa.lower_bound}, {ppa.upper_bound}]")
        print(f"  [주의] H+-translocating PPase (pmf 생성)")
        print(f"         PPi → pmf → ATP로 에너지 회수 아티팩트 가능성")
    else:
        print(f"\n[INFO] PPA_1pp가 모델에 없습니다.")
    
    # 일반 cytosolic PPase 확인
    cytosolic_ppase = ['PPA', 'PPA2', 'PPA3']
    for rxn_id in cytosolic_ppase:
        if rxn_id in model.reactions:
            rxn = model.reactions.get_by_id(rxn_id)
            flux = solution.fluxes.get(rxn_id, 0.0)
            if abs(flux) > 1e-6:
                print(f"\n[{rxn_id}] (cytosolic PPase)")
                print(f"  반응식: {rxn.reaction}")
                print(f"  플럭스: {flux:.6f}")

def check_sucdi_reaction(model):
    """SUCDi 반응 확인 (Q8 기반)"""
    print("\n" + "="*70)
    print("SUCDi 반응 확인 (Q8 기반)")
    print("="*70)
    
    if 'SUCDi' in model.reactions:
        sucdi = model.reactions.get_by_id('SUCDi')
        print(f"\n[SUCDi]")
        print(f"  반응식: {sucdi.reaction}")
        print(f"  bounds: [{sucdi.lower_bound}, {sucdi.upper_bound}]")
        print(f"  [OK] Q8 기반 SUCDi 반응이 있습니다.")
        return True
    else:
        print(f"\n[INFO] SUCDi가 모델에 없습니다.")
        return False

def fix_redox_shuttle(model):
    """레독스 셔틀 수정"""
    print("\n" + "="*70)
    print("레독스 셔틀 수정")
    print("="*70)
    
    reactions_to_remove = []
    
    # ACOAD2, ACOAD2f 제거
    for rxn_id in ['ACOAD2', 'ACOAD2f']:
        if rxn_id in model.reactions:
            rxn = model.reactions.get_by_id(rxn_id)
            print(f"\n[제거] {rxn_id}")
            print(f"  반응식: {rxn.reaction}")
            print(f"  이유: 레독스 셔틀 형성 (FADH2 → NADH 변환)")
            reactions_to_remove.append(rxn_id)
    
    # SUCD 확인 및 제거 (SUCDi가 있는 경우)
    if 'SUCDi' in model.reactions:
        if 'SUCD' in model.reactions:
            sucd = model.reactions.get_by_id('SUCD')
            print(f"\n[제거] SUCD")
            print(f"  반응식: {sucd.reaction}")
            print(f"  이유: free FAD/FADH2 사용 (SUCDi가 Q8 기반으로 있음)")
            reactions_to_remove.append('SUCD')
    else:
        print(f"\n[주의] SUCDi가 없어서 SUCD를 유지합니다.")
        print(f"       SUCDi를 추가하는 것을 권장합니다.")
    
    if reactions_to_remove:
        model.remove_reactions(reactions_to_remove)
        print(f"\n[OK] {len(reactions_to_remove)}개 반응 제거 완료")
    else:
        print(f"\n[INFO] 제거할 반응이 없습니다.")
    
    return reactions_to_remove

def fix_ppase(model):
    """PPase 수정"""
    print("\n" + "="*70)
    print("PPase 수정")
    print("="*70)
    
    # PPA_1pp를 0으로 막기 (또는 제거)
    if 'PPA_1pp' in model.reactions:
        ppa = model.reactions.get_by_id('PPA_1pp')
        print(f"\n[PPA_1pp 수정]")
        print(f"  현재 반응식: {ppa.reaction}")
        print(f"  현재 bounds: [{ppa.lower_bound}, {ppa.upper_bound}]")
        print(f"  이유: H+-translocating PPase가 에너지 회수 아티팩트 생성 가능")
        
        # 옵션 1: 제거
        # model.remove_reactions(['PPA_1pp'])
        # print(f"  [제거] PPA_1pp 제거 완료")
        
        # 옵션 2: 0으로 막기 (비가역만 유지)
        ppa.lower_bound = 0.0
        ppa.upper_bound = 0.0
        print(f"  [차단] PPA_1pp bounds를 [0.0, 0.0]으로 설정 (비활성화)")
        print(f"         일반 cytosolic PPase가 있으면 그것이 사용됩니다.")
        
        return True
    else:
        print(f"\n[INFO] PPA_1pp가 모델에 없습니다.")
        return False

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
    
    # 주요 경로 플럭스
    print(f"\n[주요 경로 플럭스]")
    key_rxns = ['EX_ac_e', 'ACS', 'CS', 'ICL', 'MALS', 'SUCDi', 'SUCD', 'MDH', 'Growth', 'ATPS4rpp']
    for rxn_id in key_rxns:
        if rxn_id in model.reactions:
            flux = solution.fluxes.get(rxn_id, 0.0)
            if abs(flux) > 1e-6:
                print(f"  {rxn_id:15s}: {flux:10.6f}")
    
    # 제거된 반응 확인
    removed_rxns = ['ACOAD2', 'ACOAD2f', 'SUCD']
    still_present = []
    for rxn_id in removed_rxns:
        if rxn_id in model.reactions:
            still_present.append(rxn_id)
    
    if still_present:
        print(f"\n[주의] 다음 반응이 여전히 모델에 있습니다: {still_present}")
    else:
        print(f"\n[OK] 레독스 셔틀 반응이 제거되었습니다.")
    
    # PPA_1pp 확인
    if 'PPA_1pp' in model.reactions:
        ppa_flux = solution.fluxes.get('PPA_1pp', 0.0)
        if abs(ppa_flux) < 1e-6:
            print(f"\n[OK] PPA_1pp가 비활성화되었습니다 (플럭스: {ppa_flux:.6f})")
        else:
            print(f"\n[주의] PPA_1pp가 여전히 작동 중입니다 (플럭스: {ppa_flux:.6f})")
    
    return solution

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_no_artificial_loops.xml"
    output_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_no_redox_shuttle.xml"
    
    print("="*70)
    print("레독스 셔틀 및 PPase 문제 수정")
    print("="*70)
    
    # 모델 로드
    model = load_model(model_path)
    model = setup_media_forced(model)
    
    # ATPM=0 설정
    if 'ATPM' in model.reactions:
        atpm = model.reactions.get_by_id('ATPM')
        atpm.lower_bound = 0.0
        atpm.upper_bound = 0.0
    
    # 제거 전 FBA
    solution_before = model.optimize()
    print(f"\n[제거 전 FBA]")
    print(f"  상태: {solution_before.status}")
    print(f"  성장률: {solution_before.objective_value:.6f}")
    
    # 문제 확인
    check_redox_shuttle(model, solution_before)
    check_ppase(model, solution_before)
    check_sucdi_reaction(model)
    
    # 수정
    removed_redox = fix_redox_shuttle(model)
    fixed_ppase = fix_ppase(model)
    
    # 제거 후 FBA
    solution_after = test_model(model)
    
    # 결과 요약
    print("\n" + "="*70)
    print("결과 요약")
    print("="*70)
    
    growth_before = solution_before.objective_value
    growth_after = solution_after.objective_value
    
    print(f"\n[성장률 비교]")
    print(f"  수정 전: {growth_before:.6f}")
    print(f"  수정 후: {growth_after:.6f}")
    print(f"  차이: {growth_after - growth_before:.6f}")
    
    print(f"\n[제거/수정된 반응]")
    for rxn_id in removed_redox:
        print(f"  {rxn_id} (제거)")
    if fixed_ppase:
        print(f"  PPA_1pp (비활성화)")
    
    if growth_after > 1e-6:
        print(f"\n[결론] 레독스 셔틀 및 PPase 수정 후에도 성장 가능 (성장률: {growth_after:.6f})")
        
        # 모델 저장
        cobra.io.write_sbml_model(model, str(output_path))
        print(f"\n[모델 저장] {output_path}")
    else:
        print(f"\n[주의] 수정 후 성장 불가 (성장률: {growth_after:.6f})")
    
    return model, solution_before, solution_after

if __name__ == "__main__":
    model, sol_before, sol_after = main()
