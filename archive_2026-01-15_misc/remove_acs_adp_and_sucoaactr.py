#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ACS_ADP와 SUCOAACTr 제거 테스트

두 반응을 모두 제거하고 FBA를 실행하여:
1. 성장률 확인
2. Acetate → Acetyl-CoA 전환 경로 확인
3. ACS 반응 사용 여부 확인
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

def test_removal(model, reactions_to_remove):
    """반응 제거 및 테스트"""
    print("\n" + "="*70)
    print("반응 제거 테스트")
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
    
    for rxn_id in reactions_to_remove:
        if rxn_id in model.reactions:
            flux = solution_before.fluxes.get(rxn_id, 0.0)
            print(f"  {rxn_id}: 플럭스 {flux:.6f}")
    
    # 반응 제거
    print(f"\n[제거할 반응]")
    reactions_to_remove_actual = []
    for rxn_id in reactions_to_remove:
        if rxn_id in model.reactions:
            rxn = model.reactions.get_by_id(rxn_id)
            print(f"  {rxn_id}: {rxn.reaction}")
            reactions_to_remove_actual.append(rxn_id)
        else:
            print(f"  {rxn_id}: 모델에 없음")
    
    if reactions_to_remove_actual:
        model.remove_reactions(reactions_to_remove_actual)
        print(f"\n[OK] {len(reactions_to_remove_actual)}개 반응 제거 완료")
    
    # 제거 후
    print("\n[제거 후]")
    solution_after = model.optimize()
    print(f"  상태: {solution_after.status}")
    print(f"  성장률: {solution_after.objective_value:.6f}")
    
    # Acetate → Acetyl-CoA 전환 경로 확인
    print("\n[Acetate → Acetyl-CoA 전환 경로 확인]")
    acetate_to_accoa_rxns = {
        'ACS': 'Acetyl-CoA synthetase (AMP-forming)',
        'ACACT1r': 'Acetyl-CoA C-acetyltransferase',
        'PTAr': 'Phosphotransacetylase',
        'ACKr': 'Acetate kinase',
    }
    
    has_pathway = False
    for rxn_id, desc in acetate_to_accoa_rxns.items():
        if rxn_id in model.reactions:
            flux = solution_after.fluxes.get(rxn_id, 0.0)
            if abs(flux) > 1e-6:
                rxn = model.reactions.get_by_id(rxn_id)
                print(f"  [사용 중] {rxn_id:15s}: 플럭스 {flux:10.6f}")
                print(f"             반응식: {rxn.reaction}")
                has_pathway = True
            else:
                print(f"  [미사용] {rxn_id:15s}: 플럭스 {flux:10.6f}")
    
    # Acetyl-CoA 생산 확인
    print("\n[Acetyl-CoA 생산 경로 확인]")
    accoa_producing = []
    for rxn in model.reactions:
        if 'accoa_c' in [m.id for m in rxn.metabolites]:
            accoa_coeff = rxn.metabolites.get(model.metabolites.get_by_id('accoa_c'), 0)
            if accoa_coeff > 0:  # 생산
                flux = solution_after.fluxes.get(rxn.id, 0.0)
                if abs(flux) > 1e-6:
                    accoa_producing.append((rxn.id, rxn.reaction, flux, accoa_coeff))
    
    print(f"  (플럭스 > 1e-6, 총 {len(accoa_producing)}개)")
    for rxn_id, rxn_eq, flux, coeff in accoa_producing[:5]:
        accoa_gen = flux * coeff
        print(f"  {rxn_id:20s}: 플럭스 {flux:10.6f}, Acetyl-CoA 생산 {accoa_gen:.6f}")
        print(f"                     반응식: {rxn_eq}")
    
    # 주요 경로 플럭스
    print("\n[주요 경로 플럭스]")
    key_rxns = ['EX_ac_e', 'CS', 'ICL', 'MALS', 'SUCDi', 'Growth']
    for rxn_id in key_rxns:
        if rxn_id in model.reactions:
            flux = solution_after.fluxes.get(rxn_id, 0.0)
            print(f"  {rxn_id:15s}: {flux:10.6f}")
    
    return solution_before, solution_after, has_pathway

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_with_BCAA_cofactors_ions_nad_transport.xml"
    output_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_without_ACS_ADP_SUCOAACTr.xml"
    
    print("="*70)
    print("ACS_ADP와 SUCOAACTr 제거 테스트")
    print("="*70)
    
    # 모델 로드
    model = load_model(model_path)
    model = setup_media_forced(model)
    
    # 제거할 반응
    reactions_to_remove = ['ACS_ADP', 'SUCOAACTr']
    
    # 제거 테스트
    solution_before, solution_after, has_pathway = test_removal(model, reactions_to_remove)
    
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
    
    if growth_after > 1e-6:
        print(f"\n[결론] 성장 가능 (성장률: {growth_after:.6f})")
        if has_pathway:
            print(f"  -> Acetate → Acetyl-CoA 전환 경로 있음 (ACS 등)")
        else:
            print(f"  -> Acetate → Acetyl-CoA 직접 전환 경로 없음")
            print(f"  -> 다른 경로로 Acetyl-CoA 생성 중")
        
        # 모델 저장
        cobra.io.write_sbml_model(model, str(output_path))
        print(f"\n[모델 저장] {output_path}")
    else:
        print(f"\n[결론] 성장 불가 (성장률: {growth_after:.6f})")
        print(f"  -> ACS_ADP 또는 SUCOAACTr 중 하나는 필수일 수 있음")
    
    return model, solution_before, solution_after

if __name__ == "__main__":
    model, sol_before, sol_after = main()
