#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Acetate 사용 경로 확인

Acetate가 19만큼 들어갔는데 Acetyl-CoA로 직접 전환되지 않는 이유 확인
- Acetate가 어디로 가는지 확인
- Acetate를 사용하는 모든 반응 확인
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

def analyze_acetate_usage(model, solution):
    """Acetate 사용 경로 분석"""
    print("\n" + "="*70)
    print("Acetate 사용 경로 분석")
    print("="*70)
    
    ac_c = model.metabolites.get_by_id('ac_c')
    
    # Acetate를 소비하는 반응
    consuming_rxns = []
    # Acetate를 생산하는 반응
    producing_rxns = []
    
    for rxn in model.reactions:
        if ac_c in rxn.metabolites:
            coeff = rxn.metabolites[ac_c]
            flux = solution.fluxes.get(rxn.id, 0.0)
            
            if abs(flux) > 1e-6:
                if coeff < 0:  # 소비
                    ac_consumed = abs(flux * coeff)
                    consuming_rxns.append((rxn.id, rxn.reaction, flux, coeff, ac_consumed))
                elif coeff > 0:  # 생산
                    ac_produced = flux * coeff
                    producing_rxns.append((rxn.id, rxn.reaction, flux, coeff, ac_produced))
    
    # 소비 정렬
    consuming_rxns.sort(key=lambda x: x[4], reverse=True)
    # 생산 정렬
    producing_rxns.sort(key=lambda x: x[4], reverse=True)
    
    print(f"\n[Acetate 소비 반응] (플럭스 * 계수 기준)")
    total_consumption = 0.0
    for rxn_id, rxn_eq, flux, coeff, consumption in consuming_rxns:
        print(f"  {rxn_id:20s}: 플럭스 {flux:10.6f}, 계수 {coeff:+.2f}, 소비량 {consumption:10.6f}")
        print(f"                     반응식: {rxn_eq}")
        total_consumption += consumption
    print(f"\n  총 소비량: {total_consumption:.6f}")
    
    print(f"\n[Acetate 생산 반응] (플럭스 * 계수 기준)")
    total_production = 0.0
    for rxn_id, rxn_eq, flux, coeff, production in producing_rxns:
        print(f"  {rxn_id:20s}: 플럭스 {flux:10.6f}, 계수 {coeff:+.2f}, 생산량 {production:10.6f}")
        print(f"                     반응식: {rxn_eq}")
        total_production += production
    print(f"\n  총 생산량: {total_production:.6f}")
    
    # Acetate balance
    if 'EX_ac_e' in model.reactions:
        ex_ac_flux = solution.fluxes.get('EX_ac_e', 0.0)
        acetate_in = abs(ex_ac_flux)  # 음수이므로 절댓값
        net_consumption = total_consumption - total_production
        
        print(f"\n[Acetate 균형]")
        print(f"  Uptake (EX_ac_e): {acetate_in:.6f}")
        print(f"  총 소비량: {total_consumption:.6f}")
        print(f"  총 생산량: {total_production:.6f}")
        print(f"  순 소비량: {net_consumption:.6f}")
        print(f"  균형: {acetate_in - net_consumption:.6f} (0에 가까워야 함)")
    
    return consuming_rxns, producing_rxns, total_consumption, total_production

def check_acetate_transport(model, solution):
    """Acetate 수송 반응 확인"""
    print("\n" + "="*70)
    print("Acetate 수송 반응 확인")
    print("="*70)
    
    transport_rxns = []
    
    for rxn in model.reactions:
        if 'ac_e' in [m.id for m in rxn.metabolites] and 'ac_c' in [m.id for m in rxn.metabolites]:
            flux = solution.fluxes.get(rxn.id, 0.0)
            if abs(flux) > 1e-6:
                ac_e_coeff = rxn.metabolites.get(model.metabolites.get_by_id('ac_e'), 0)
                ac_c_coeff = rxn.metabolites.get(model.metabolites.get_by_id('ac_c'), 0)
                transport_rxns.append((rxn.id, rxn.reaction, flux, ac_e_coeff, ac_c_coeff))
    
    print(f"\n[Acetate 수송 반응]")
    for rxn_id, rxn_eq, flux, ac_e_coeff, ac_c_coeff in transport_rxns:
        print(f"  {rxn_id}")
        print(f"    반응식: {rxn_eq}")
        print(f"    플럭스: {flux:.6f}")
        print(f"    ac_e 계수: {ac_e_coeff:+.2f}")
        print(f"    ac_c 계수: {ac_c_coeff:+.2f}")
    
    return transport_rxns

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_without_ACS_ADP_SUCOAACTr.xml"
    
    print("="*70)
    print("Acetate 사용 경로 확인")
    print("="*70)
    
    # 모델 로드
    model = load_model(model_path)
    model = setup_media_forced(model)
    
    # ATPM=0 설정
    if 'ATPM' in model.reactions:
        atpm = model.reactions.get_by_id('ATPM')
        atpm.lower_bound = 0.0
        atpm.upper_bound = 0.0
    
    # FBA 실행
    solution = model.optimize()
    print(f"\n[FBA 결과]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    
    # Acetate 수송 확인
    transport_rxns = check_acetate_transport(model, solution)
    
    # Acetate 사용 경로 분석
    consuming, producing, total_cons, total_prod = analyze_acetate_usage(model, solution)
    
    # 종합 분석
    print("\n" + "="*70)
    print("종합 분석")
    print("="*70)
    
    if 'EX_ac_e' in model.reactions:
        ex_ac_flux = solution.fluxes.get('EX_ac_e', 0.0)
        acetate_in = abs(ex_ac_flux)
        net_cons = total_cons - total_prod
        
        print(f"\n[Acetate 흐름]")
        print(f"  Uptake: {acetate_in:.6f}")
        print(f"  순 소비량: {net_cons:.6f}")
        
        if abs(acetate_in - net_cons) < 1e-3:
            print(f"\n[OK] Acetate 균형이 맞습니다.")
        else:
            print(f"\n[주의] Acetate 균형이 맞지 않습니다!")
            print(f"       차이: {abs(acetate_in - net_cons):.6f}")
        
        # Acetate → Acetyl-CoA 직접 전환 확인
        ac_to_accoa_consumption = 0.0
        for rxn_id, rxn_eq, flux, coeff, consumption in consuming:
            if 'accoa_c' in [m.id for m in model.reactions.get_by_id(rxn_id).metabolites]:
                ac_to_accoa_consumption += consumption
        
        print(f"\n[Acetate → Acetyl-CoA 직접 전환]")
        print(f"  Acetate 소비량 (Acetyl-CoA 생성을 위해): {ac_to_accoa_consumption:.6f}")
        if ac_to_accoa_consumption < 1e-6:
            print(f"  [문제] Acetate가 Acetyl-CoA로 직접 전환되지 않고 있습니다!")
        else:
            print(f"  [OK] Acetate가 Acetyl-CoA로 전환되고 있습니다.")
    
    return model, solution

if __name__ == "__main__":
    model, solution = main()
