#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
인공 루프 제거

문제점:
1. ACCOAL: 이름은 "Acetate-CoA ligase"인데 반응식은 ppa(propionate) 사용
2. ACCOAL + APAT_1 + PACPT_1: 에너지 회수 루프 (수학적 아티팩트)
3. GALpts 관련 PEP 루프: 비현실적인 순환
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

def check_problematic_reactions(model, solution):
    """문제가 되는 반응 확인"""
    print("\n" + "="*70)
    print("문제가 되는 반응 확인")
    print("="*70)
    
    # 1. ACCOAL 확인
    print("\n[1. ACCOAL 반응 확인]")
    if 'ACCOAL' in model.reactions:
        accoal = model.reactions.get_by_id('ACCOAL')
        accoal_flux = solution.fluxes.get('ACCOAL', 0.0)
        print(f"  ID: {accoal.id}")
        print(f"  이름: {accoal.name if accoal.name else ''}")
        print(f"  반응식: {accoal.reaction}")
        print(f"  플럭스: {accoal_flux:.6f}")
        print(f"  [문제] 이름은 'Acetate-CoA ligase'인데 반응식은 ppa(propionate) 사용")
    else:
        print(f"  [INFO] ACCOAL이 모델에 없습니다.")
    
    # 2. APAT_1, PACPT_1 확인
    print("\n[2. APAT_1, PACPT_1 확인]")
    problematic_rxns = ['APAT_1', 'PACPT_1']
    for rxn_id in problematic_rxns:
        if rxn_id in model.reactions:
            rxn = model.reactions.get_by_id(rxn_id)
            flux = solution.fluxes.get(rxn_id, 0.0)
            print(f"  {rxn_id:10s}: 플럭스 {flux:10.6f}")
            print(f"             반응식: {rxn.reaction}")
            print(f"             bounds: [{rxn.lower_bound}, {rxn.upper_bound}]")
        else:
            print(f"  {rxn_id:10s}: 모델에 없음")
    
    # 3. GALpts 관련 확인
    print("\n[3. GALpts 관련 반응 확인]")
    gal_rxns = ['GALpts', 'GALt2', 'A6PAG']
    for rxn_id in gal_rxns:
        if rxn_id in model.reactions:
            rxn = model.reactions.get_by_id(rxn_id)
            flux = solution.fluxes.get(rxn_id, 0.0)
            if abs(flux) > 1e-6:
                print(f"  {rxn_id:10s}: 플럭스 {flux:10.6f}")
                print(f"             반응식: {rxn.reaction}")
                print(f"             bounds: [{rxn.lower_bound}, {rxn.upper_bound}]")
        else:
            print(f"  {rxn_id:10s}: 모델에 없음")
    
    # 4. 루프 확인
    print("\n[4. 에너지 회수 루프 확인]")
    if all(rxn_id in model.reactions for rxn_id in ['ACCOAL', 'APAT_1', 'PACPT_1']):
        accoal_flux = solution.fluxes.get('ACCOAL', 0.0)
        apat_flux = solution.fluxes.get('APAT_1', 0.0)
        pacpt_flux = solution.fluxes.get('PACPT_1', 0.0)
        
        print(f"  ACCOAL: {accoal_flux:.6f}")
        print(f"  APAT_1: {apat_flux:.6f}")
        print(f"  PACPT_1: {pacpt_flux:.6f}")
        
        if abs(accoal_flux) > 1e-6 and abs(apat_flux) > 1e-6 and abs(pacpt_flux) > 1e-6:
            if abs(abs(accoal_flux) - abs(apat_flux)) < 1e-3 and abs(abs(accoal_flux) - abs(pacpt_flux)) < 1e-3:
                print(f"  [문제] 세 반응이 거의 같은 크기로 맞물려 에너지 회수 루프를 형성합니다!")

def remove_problematic_reactions(model):
    """문제가 되는 반응 제거"""
    print("\n" + "="*70)
    print("문제가 되는 반응 제거")
    print("="*70)
    
    reactions_to_remove = []
    
    # 1. ACCOAL 제거 (이름과 반응식 불일치)
    if 'ACCOAL' in model.reactions:
        accoal = model.reactions.get_by_id('ACCOAL')
        print(f"\n[제거] ACCOAL")
        print(f"  이름: {accoal.name if accoal.name else ''}")
        print(f"  반응식: {accoal.reaction}")
        print(f"  이유: 이름은 'Acetate-CoA ligase'인데 반응식은 ppa(propionate) 사용")
        reactions_to_remove.append('ACCOAL')
    
    # 2. APAT_1, PACPT_1 제거 (에너지 회수 루프)
    for rxn_id in ['APAT_1', 'PACPT_1']:
        if rxn_id in model.reactions:
            rxn = model.reactions.get_by_id(rxn_id)
            print(f"\n[제거] {rxn_id}")
            print(f"  반응식: {rxn.reaction}")
            print(f"  이유: 에너지 회수 루프 형성 (수학적 아티팩트)")
            reactions_to_remove.append(rxn_id)
    
    # 3. GALpts 관련 루프 제거
    gal_rxns_to_remove = []
    if 'GALpts' in model.reactions:
        galpts = model.reactions.get_by_id('GALpts')
        print(f"\n[제거] GALpts")
        print(f"  반응식: {galpts.reaction}")
        print(f"  이유: PEP 루프 형성 (비현실적 순환)")
        gal_rxns_to_remove.append('GALpts')
    
    if 'GALt2' in model.reactions:
        galt2 = model.reactions.get_by_id('GALt2')
        print(f"\n[제거] GALt2")
        print(f"  반응식: {galt2.reaction}")
        print(f"  이유: PEP 루프 형성")
        gal_rxns_to_remove.append('GALt2')
    
    if 'A6PAG' in model.reactions:
        a6pag = model.reactions.get_by_id('A6PAG')
        print(f"\n[제거] A6PAG")
        print(f"  반응식: {a6pag.reaction}")
        print(f"  이유: PEP 루프 형성")
        gal_rxns_to_remove.append('A6PAG')
    
    reactions_to_remove.extend(gal_rxns_to_remove)
    
    # 반응 제거
    if reactions_to_remove:
        print(f"\n[제거할 반응 목록]")
        for rxn_id in reactions_to_remove:
            print(f"  {rxn_id}")
        
        model.remove_reactions(reactions_to_remove)
        print(f"\n[OK] {len(reactions_to_remove)}개 반응 제거 완료")
    else:
        print(f"\n[INFO] 제거할 반응이 없습니다.")
    
    return reactions_to_remove

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
    key_rxns = ['EX_ac_e', 'ACS', 'CS', 'ICL', 'MALS', 'MDH', 'MAEB', 'Growth']
    for rxn_id in key_rxns:
        if rxn_id in model.reactions:
            flux = solution.fluxes.get(rxn_id, 0.0)
            print(f"  {rxn_id:15s}: {flux:10.6f}")
    
    # 제거된 반응이 여전히 있는지 확인
    removed_rxns = ['ACCOAL', 'APAT_1', 'PACPT_1', 'GALpts', 'GALt2', 'A6PAG']
    still_present = []
    for rxn_id in removed_rxns:
        if rxn_id in model.reactions:
            still_present.append(rxn_id)
    
    if still_present:
        print(f"\n[주의] 다음 반응이 여전히 모델에 있습니다: {still_present}")
    else:
        print(f"\n[OK] 모든 문제 반응이 제거되었습니다.")
    
    return solution

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_with_maeB_no_PEPCK_final.xml"
    output_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_no_artificial_loops.xml"
    
    print("="*70)
    print("인공 루프 제거")
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
    
    # 문제 반응 확인
    check_problematic_reactions(model, solution_before)
    
    # 문제 반응 제거
    removed_rxns = remove_problematic_reactions(model)
    
    # 제거 후 FBA
    solution_after = test_model(model)
    
    # 결과 요약
    print("\n" + "="*70)
    print("결과 요약")
    print("="*70)
    
    growth_before = solution_before.objective_value
    growth_after = solution_after.objective_value
    
    print(f"\n[성장률 비교]")
    print(f"  제거 전: {growth_before:.6f}")
    print(f"  제거 후: {growth_after:.6f}")
    print(f"  차이: {growth_after - growth_before:.6f}")
    
    print(f"\n[제거된 반응]")
    for rxn_id in removed_rxns:
        print(f"  {rxn_id}")
    
    if growth_after > 1e-6:
        print(f"\n[결론] 인공 루프 제거 후에도 성장 가능 (성장률: {growth_after:.6f})")
        
        # 모델 저장
        cobra.io.write_sbml_model(model, str(output_path))
        print(f"\n[모델 저장] {output_path}")
    else:
        print(f"\n[주의] 인공 루프 제거 후 성장 불가 (성장률: {growth_after:.6f})")
    
    return model, solution_before, solution_after

if __name__ == "__main__":
    model, sol_before, sol_after = main()
