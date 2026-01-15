#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ACS_ADP 제거 후 Acetate 경로 분석

ACS_ADP가 없어도 성장률이 동일한 이유 확인
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

def find_acetate_to_accoa_pathways(model):
    """Acetate → Acetyl-CoA 전환 가능한 모든 경로 찾기"""
    print("\n" + "="*70)
    print("Acetate → Acetyl-CoA 전환 경로 찾기")
    print("="*70)
    
    ac_to_accoa_rxns = []
    
    for rxn in model.reactions:
        mets = [m.id for m in rxn.metabolites]
        if 'ac_c' in mets and 'accoa_c' in mets:
            ac_coeff = rxn.metabolites.get(model.metabolites.get_by_id('ac_c'), 0)
            accoa_coeff = rxn.metabolites.get(model.metabolites.get_by_id('accoa_c'), 0)
            
            # ac_c를 소비하고 accoa_c를 생성하는 반응
            if ac_coeff < 0 and accoa_coeff > 0:
                ac_to_accoa_rxns.append((rxn.id, rxn.reaction, ac_coeff, accoa_coeff))
    
    print(f"\n[찾은 반응] (총 {len(ac_to_accoa_rxns)}개)")
    for rxn_id, rxn_eq, ac_coeff, accoa_coeff in ac_to_accoa_rxns:
        print(f"\n  {rxn_id}")
        print(f"    반응식: {rxn_eq}")
        print(f"    ac_c 계수: {ac_coeff:+.2f}")
        print(f"    accoa_c 계수: {accoa_coeff:+.2f}")
    
    return ac_to_accoa_rxns

def check_accoa_production(model, solution):
    """Acetyl-CoA 생산 가능 경로 확인"""
    print("\n" + "="*70)
    print("Acetyl-CoA 생산 경로 확인")
    print("="*70)
    
    # Acetyl-CoA를 생산하는 모든 반응 찾기
    accoa_producing = []
    
    for rxn in model.reactions:
        if 'accoa_c' in [m.id for m in rxn.metabolites]:
            accoa_coeff = rxn.metabolites.get(model.metabolites.get_by_id('accoa_c'), 0)
            if accoa_coeff > 0:  # 생산
                flux = solution.fluxes.get(rxn.id, 0.0)
                if abs(flux) > 1e-6:
                    accoa_producing.append((rxn.id, rxn.reaction, flux, accoa_coeff))
    
    print(f"\n[Acetyl-CoA 생산 반응] (플럭스 > 1e-6, 총 {len(accoa_producing)}개)")
    for rxn_id, rxn_eq, flux, coeff in accoa_producing[:10]:
        accoa_gen = flux * coeff
        print(f"\n  {rxn_id}")
        print(f"    반응식: {rxn_eq}")
        print(f"    플럭스: {flux:.6f}")
        print(f"    accoa_c 계수: {coeff:+.2f}")
        print(f"    Acetyl-CoA 생산량: {accoa_gen:.6f}")
    
    return accoa_producing

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_with_BCAA_cofactors_ions_nad_transport.xml"
    
    print("="*70)
    print("ACS_ADP 제거 후 Acetate 경로 분석")
    print("="*70)
    
    # 모델 로드
    model = load_model(model_path)
    model = setup_media_forced(model)
    
    # ATPM=0 설정
    if 'ATPM' in model.reactions:
        atpm = model.reactions.get_by_id('ATPM')
        atpm.lower_bound = 0.0
        atpm.upper_bound = 0.0
    
    # ACS_ADP 제거 전
    print("\n[ACS_ADP 제거 전]")
    solution_with = model.optimize()
    print(f"  성장률: {solution_with.objective_value:.6f}")
    
    # Acetate → Acetyl-CoA 경로 찾기
    ac_to_accoa_rxns = find_acetate_to_accoa_pathways(model)
    
    # ACS_ADP 제거
    if 'ACS_ADP' in model.reactions:
        model.remove_reactions(['ACS_ADP'])
        print("\n[ACS_ADP 제거 완료]")
    
    # ACS_ADP 제거 후
    print("\n[ACS_ADP 제거 후]")
    solution_without = model.optimize()
    print(f"  성장률: {solution_without.objective_value:.6f}")
    
    # Acetyl-CoA 생산 경로 확인
    accoa_producing = check_accoa_production(model, solution_without)
    
    # 주요 반응 플럭스 비교
    print("\n" + "="*70)
    print("주요 반응 플럭스 비교")
    print("="*70)
    
    key_rxns = ['EX_ac_e', 'CS', 'ICL', 'MALS', 'Growth']
    print(f"\n{'반응':<15} {'ACS_ADP 있을 때':<20} {'ACS_ADP 없을 때':<20} {'차이':<20}")
    print("-" * 75)
    
    for rxn_id in key_rxns:
        if rxn_id in model.reactions:
            flux_with = solution_with.fluxes.get(rxn_id, 0.0) if rxn_id in solution_with.fluxes else 0.0
            flux_without = solution_without.fluxes.get(rxn_id, 0.0)
            diff = flux_without - flux_with
            print(f"{rxn_id:<15} {flux_with:<20.6f} {flux_without:<20.6f} {diff:<20.6f}")
    
    return model, solution_with, solution_without

if __name__ == "__main__":
    model, sol_with, sol_without = main()
