#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Acetyl-CoA 균형 분석

Acetate uptake가 -19.0인데 ACS flux가 123.81로 높은 이유 분석
- Acetyl-CoA의 생성/소비 균형 확인
- 실제 Acetate → Acetyl-CoA 전환량 확인
- Acetyl-CoA 재순환 여부 확인
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

def analyze_acetylcoa_balance(model, solution):
    """Acetyl-CoA 균형 분석"""
    print("\n" + "="*70)
    print("Acetyl-CoA 균형 분석")
    print("="*70)
    
    accoa_c = model.metabolites.get_by_id('accoa_c')
    
    # Acetyl-CoA를 생산하는 반응
    producing_rxns = []
    # Acetyl-CoA를 소비하는 반응
    consuming_rxns = []
    
    for rxn in model.reactions:
        if accoa_c in rxn.metabolites:
            coeff = rxn.metabolites[accoa_c]
            flux = solution.fluxes.get(rxn.id, 0.0)
            
            if abs(flux) > 1e-6:
                if coeff > 0:  # 생산
                    producing_rxns.append((rxn.id, rxn.reaction, flux, coeff, flux * coeff))
                elif coeff < 0:  # 소비
                    consuming_rxns.append((rxn.id, rxn.reaction, flux, coeff, abs(flux * coeff)))
    
    # 생산 정렬 (생산량 기준)
    producing_rxns.sort(key=lambda x: x[4], reverse=True)
    # 소비 정렬 (소비량 기준)
    consuming_rxns.sort(key=lambda x: x[4], reverse=True)
    
    print(f"\n[Acetyl-CoA 생산 반응] (플럭스 * 계수 기준, 상위 10개)")
    total_production = 0.0
    for rxn_id, rxn_eq, flux, coeff, production in producing_rxns[:10]:
        print(f"  {rxn_id:20s}: 플럭스 {flux:10.6f}, 계수 {coeff:+.2f}, 생산량 {production:10.6f}")
        print(f"                     반응식: {rxn_eq}")
        total_production += production
    print(f"\n  총 생산량: {total_production:.6f}")
    
    print(f"\n[Acetyl-CoA 소비 반응] (플럭스 * 계수 기준, 상위 10개)")
    total_consumption = 0.0
    for rxn_id, rxn_eq, flux, coeff, consumption in consuming_rxns[:10]:
        print(f"  {rxn_id:20s}: 플럭스 {flux:10.6f}, 계수 {coeff:+.2f}, 소비량 {consumption:10.6f}")
        print(f"                     반응식: {rxn_eq}")
        total_consumption += consumption
    print(f"\n  총 소비량: {total_consumption:.6f}")
    
    print(f"\n[균형]")
    print(f"  총 생산량: {total_production:.6f}")
    print(f"  총 소비량: {total_consumption:.6f}")
    print(f"  균형: {total_production - total_consumption:.6f} (0에 가까워야 함)")
    
    return producing_rxns, consuming_rxns, total_production, total_consumption

def analyze_acetate_to_accoa(model, solution):
    """Acetate → Acetyl-CoA 직접 전환 분석"""
    print("\n" + "="*70)
    print("Acetate → Acetyl-CoA 직접 전환 분석")
    print("="*70)
    
    # Acetate uptake
    if 'EX_ac_e' in model.reactions:
        ex_ac_flux = solution.fluxes.get('EX_ac_e', 0.0)
        print(f"\n[Acetate Uptake]")
        print(f"  EX_ac_e: {ex_ac_flux:.6f}")
    
    # Acetate → Acetyl-CoA 전환 반응
    acetate_to_accoa = []
    
    for rxn in model.reactions:
        if 'ac_c' in [m.id for m in rxn.metabolites] and 'accoa_c' in [m.id for m in rxn.metabolites]:
            ac_coeff = rxn.metabolites.get(model.metabolites.get_by_id('ac_c'), 0)
            accoa_coeff = rxn.metabolites.get(model.metabolites.get_by_id('accoa_c'), 0)
            
            # ac_c를 소비하고 accoa_c를 생성하는 반응
            if ac_coeff < 0 and accoa_coeff > 0:
                flux = solution.fluxes.get(rxn.id, 0.0)
                if abs(flux) > 1e-6:
                    ac_consumed = abs(flux * ac_coeff)
                    accoa_produced = flux * accoa_coeff
                    acetate_to_accoa.append((rxn.id, rxn.reaction, flux, ac_consumed, accoa_produced))
    
    print(f"\n[Acetate → Acetyl-CoA 전환 반응]")
    total_ac_consumed = 0.0
    total_accoa_produced = 0.0
    
    for rxn_id, rxn_eq, flux, ac_consumed, accoa_produced in acetate_to_accoa:
        print(f"\n  {rxn_id}")
        print(f"    반응식: {rxn_eq}")
        print(f"    플럭스: {flux:.6f}")
        print(f"    Acetate 소비: {ac_consumed:.6f}")
        print(f"    Acetyl-CoA 생성: {accoa_produced:.6f}")
        total_ac_consumed += ac_consumed
        total_accoa_produced += accoa_produced
    
    print(f"\n  [총합]")
    print(f"    총 Acetate 소비: {total_ac_consumed:.6f}")
    print(f"    총 Acetyl-CoA 생성: {total_accoa_produced:.6f}")
    print(f"    Acetate uptake: {abs(ex_ac_flux):.6f}")
    
    if abs(total_ac_consumed - abs(ex_ac_flux)) > 1e-3:
        print(f"\n  [주의] Acetate 소비량 ({total_ac_consumed:.6f})과 uptake ({abs(ex_ac_flux):.6f})가 다릅니다!")
    
    return acetate_to_accoa, total_ac_consumed, total_accoa_produced

def check_acetylcoa_recycling(model, solution):
    """Acetyl-CoA 재순환 확인"""
    print("\n" + "="*70)
    print("Acetyl-CoA 재순환 확인")
    print("="*70)
    
    # Acetyl-CoA를 생성하는 반응 중에서 Acetate를 사용하지 않는 반응
    recycling_rxns = []
    
    for rxn in model.reactions:
        if 'accoa_c' in [m.id for m in rxn.metabolites]:
            accoa_coeff = rxn.metabolites.get(model.metabolites.get_by_id('accoa_c'), 0)
            if accoa_coeff > 0:  # 생산
                # Acetate를 사용하지 않는 경우
                if 'ac_c' not in [m.id for m in rxn.metabolites]:
                    flux = solution.fluxes.get(rxn.id, 0.0)
                    if abs(flux) > 1e-6:
                        accoa_produced = flux * accoa_coeff
                        recycling_rxns.append((rxn.id, rxn.reaction, flux, accoa_produced))
    
    recycling_rxns.sort(key=lambda x: x[3], reverse=True)
    
    print(f"\n[Acetyl-CoA 재순환 반응] (Acetate를 사용하지 않고 Acetyl-CoA 생성)")
    print(f"  (상위 10개)")
    total_recycling = 0.0
    
    for rxn_id, rxn_eq, flux, accoa_produced in recycling_rxns[:10]:
        print(f"  {rxn_id:20s}: 플럭스 {flux:10.6f}, Acetyl-CoA 생산 {accoa_produced:10.6f}")
        print(f"                     반응식: {rxn_eq}")
        total_recycling += accoa_produced
    
    print(f"\n  총 재순환 생산량: {total_recycling:.6f}")
    
    return recycling_rxns, total_recycling

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_without_ACS_ADP_SUCOAACTr.xml"
    
    print("="*70)
    print("Acetyl-CoA 균형 분석")
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
    
    # Acetyl-CoA 균형 분석
    producing, consuming, total_prod, total_cons = analyze_acetylcoa_balance(model, solution)
    
    # Acetate → Acetyl-CoA 직접 전환 분석
    acetate_to_accoa, ac_consumed, accoa_from_ac = analyze_acetate_to_accoa(model, solution)
    
    # Acetyl-CoA 재순환 확인
    recycling_rxns, total_recycling = check_acetylcoa_recycling(model, solution)
    
    # 종합 분석
    print("\n" + "="*70)
    print("종합 분석")
    print("="*70)
    
    print(f"\n[탄소 흐름]")
    if 'EX_ac_e' in model.reactions:
        ex_ac_flux = solution.fluxes.get('EX_ac_e', 0.0)
        print(f"  Acetate uptake: {abs(ex_ac_flux):.6f}")
    print(f"  Acetate → Acetyl-CoA 직접 전환: {accoa_from_ac:.6f}")
    print(f"  Acetyl-CoA 재순환 생산: {total_recycling:.6f}")
    print(f"  총 Acetyl-CoA 생산: {total_prod:.6f}")
    
    if abs(accoa_from_ac - abs(ex_ac_flux)) < 1e-3:
        print(f"\n[결론] Acetate → Acetyl-CoA 직접 전환량이 uptake와 일치합니다.")
    else:
        print(f"\n[주의] Acetate → Acetyl-CoA 직접 전환량 ({accoa_from_ac:.6f})이")
        print(f"       uptake ({abs(ex_ac_flux):.6f})와 다릅니다.")
        print(f"       차이: {abs(accoa_from_ac - abs(ex_ac_flux)):.6f}")
    
    if total_recycling > 1e-3:
        print(f"\n[발견] Acetyl-CoA가 재순환되고 있습니다!")
        print(f"       재순환 생산량: {total_recycling:.6f}")
        print(f"       이로 인해 ACS 플럭스가 Acetate uptake보다 높을 수 있습니다.")
    
    return model, solution

if __name__ == "__main__":
    model, solution = main()
