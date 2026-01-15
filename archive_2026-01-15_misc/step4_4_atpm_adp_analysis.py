#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 4-4: ATPM-ADP 원인 확인

가설: ATPM=0일 때 ADP가 재생산되지 않아서 전체 에너지 대사가 멈춤

확인 사항:
1. ATPM 반응식 확인 (ADP 생성 여부)
2. ATPM=0일 때 ADP 생성량 확인
3. ATP 생성 경로에 ADP 필요 여부 확인
4. ATPM과 ADP/ATP 사이클의 관계 확인
5. ATPM=0 vs ATPM>0일 때 ADP/ATP 비율 확인
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

def check_atpm_reaction(model):
    """ATPM 반응식 확인"""
    print("\n" + "="*70)
    print("1. ATPM 반응식 확인")
    print("="*70)
    
    if 'ATPM' not in model.reactions:
        print("[ERROR] ATPM 반응이 모델에 없습니다.")
        return None
    
    atpm = model.reactions.get_by_id('ATPM')
    print(f"\n[ATPM 반응]")
    print(f"  ID: {atpm.id}")
    print(f"  이름: {atpm.name}")
    print(f"  반응식: {atpm.reaction}")
    print(f"  bounds: [{atpm.lower_bound}, {atpm.upper_bound}]")
    
    # ADP/ATP 포함 여부 확인
    met_ids = [m.id for m in atpm.metabolites]
    print(f"\n  관련 metabolite:")
    for met_id in met_ids:
        met = model.metabolites.get_by_id(met_id)
        coeff = atpm.metabolites[met]
        print(f"    {met_id}: 계수 {coeff:+.6f} ({met.name})")
    
    # ADP 생성 여부 확인
    has_adp = 'adp_c' in met_ids or 'ADP' in str(atpm.reaction)
    has_atp = 'atp_c' in met_ids or 'ATP' in str(atpm.reaction)
    
    print(f"\n  [ADP 포함 여부] {has_adp}")
    print(f"  [ATP 포함 여부] {has_atp}")
    
    if has_adp:
        adp_coeff = atpm.metabolites.get(model.metabolites.get_by_id('adp_c'), 0)
        print(f"  [ADP 계수] {adp_coeff:+.6f} (양수면 생성, 음수면 소비)")
    
    return atpm

def check_adp_atp_balance_at_atpm0(model):
    """ATPM=0일 때 ADP/ATP 균형 확인"""
    print("\n" + "="*70)
    print("2. ATPM=0일 때 ADP/ATP 균형 확인")
    print("="*70)
    
    # ATPM=0 설정
    atpm = model.reactions.get_by_id('ATPM')
    original_lb = atpm.lower_bound
    original_ub = atpm.upper_bound
    
    atpm.lower_bound = 0.0
    atpm.upper_bound = 0.0
    print(f"\n[ATPM 설정] bounds: [{atpm.lower_bound}, {atpm.upper_bound}]")
    
    # FBA 실행
    solution = model.optimize()
    print(f"\n[FBA 결과]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    
    # ADP/ATP 관련 플럭스 확인
    atp_producing_rxns = []
    atp_consuming_rxns = []
    adp_producing_rxns = []
    adp_consuming_rxns = []
    
    for rxn in model.reactions:
        if 'atp_c' in [m.id for m in rxn.metabolites]:
            flux = solution.fluxes.get(rxn.id, 0.0)
            if abs(flux) > 1e-6:
                atp_coeff = rxn.metabolites.get(model.metabolites.get_by_id('atp_c'), 0)
                if atp_coeff > 0:  # ATP 생성
                    atp_producing_rxns.append((rxn.id, flux, atp_coeff))
                elif atp_coeff < 0:  # ATP 소비
                    atp_consuming_rxns.append((rxn.id, flux, atp_coeff))
        
        if 'adp_c' in [m.id for m in rxn.metabolites]:
            flux = solution.fluxes.get(rxn.id, 0.0)
            if abs(flux) > 1e-6:
                adp_coeff = rxn.metabolites.get(model.metabolites.get_by_id('adp_c'), 0)
                if adp_coeff > 0:  # ADP 생성
                    adp_producing_rxns.append((rxn.id, flux, adp_coeff))
                elif adp_coeff < 0:  # ADP 소비
                    adp_consuming_rxns.append((rxn.id, flux, adp_coeff))
    
    print(f"\n[ATP 생성 반응] (플럭스 > 1e-6)")
    if atp_producing_rxns:
        for rxn_id, flux, coeff in sorted(atp_producing_rxns, key=lambda x: abs(x[1]), reverse=True)[:10]:
            print(f"  {rxn_id:20s}: 플럭스 {flux:10.6f}, ATP 계수 {coeff:+.2f}")
    else:
        print(f"  없음")
    
    print(f"\n[ATP 소비 반응] (플럭스 > 1e-6)")
    if atp_consuming_rxns:
        for rxn_id, flux, coeff in sorted(atp_consuming_rxns, key=lambda x: abs(x[1]), reverse=True)[:10]:
            print(f"  {rxn_id:20s}: 플럭스 {flux:10.6f}, ATP 계수 {coeff:+.2f}")
    else:
        print(f"  없음")
    
    print(f"\n[ADP 생성 반응] (플럭스 > 1e-6)")
    if adp_producing_rxns:
        for rxn_id, flux, coeff in sorted(adp_producing_rxns, key=lambda x: abs(x[1]), reverse=True)[:10]:
            print(f"  {rxn_id:20s}: 플럭스 {flux:10.6f}, ADP 계수 {coeff:+.2f}")
    else:
        print(f"  없음")
    
    print(f"\n[ADP 소비 반응] (플럭스 > 1e-6)")
    if adp_consuming_rxns:
        for rxn_id, flux, coeff in sorted(adp_consuming_rxns, key=lambda x: abs(x[1]), reverse=True)[:10]:
            print(f"  {rxn_id:20s}: 플럭스 {flux:10.6f}, ADP 계수 {coeff:+.2f}")
    else:
        print(f"  없음")
    
    # 원래 bounds 복원
    atpm.lower_bound = original_lb
    atpm.upper_bound = original_ub
    
    return {
        'atp_producing': atp_producing_rxns,
        'atp_consuming': atp_consuming_rxns,
        'adp_producing': adp_producing_rxns,
        'adp_consuming': adp_consuming_rxns,
    }

def check_atp_synthase_requirements(model):
    """ATP 생성 경로 (ATP Synthase 등)에 ADP 필요 여부 확인"""
    print("\n" + "="*70)
    print("3. ATP 생성 경로에 ADP 필요 여부 확인")
    print("="*70)
    
    # ATP 생성 주요 반응들
    atp_synthesis_keywords = ['ATPS', 'ATPase', 'ATP synthase', 'ATPase']
    atp_synthesis_rxns = []
    
    for rxn in model.reactions:
        # ATP를 생성하는 반응 찾기
        if 'atp_c' in [m.id for m in rxn.metabolites]:
            atp_coeff = rxn.metabolites.get(model.metabolites.get_by_id('atp_c'), 0)
            if atp_coeff > 0:  # ATP 생성
                adp_coeff = rxn.metabolites.get(model.metabolites.get_by_id('adp_c'), 0)
                if adp_coeff < 0:  # ADP 소비 (ADP가 필요)
                    atp_synthesis_rxns.append((rxn.id, rxn.reaction, atp_coeff, adp_coeff))
    
    print(f"\n[ATP 생성 반응 중 ADP를 소비하는 반응]")
    print(f"  (총 {len(atp_synthesis_rxns)}개)")
    
    # 상위 10개만 출력
    for rxn_id, rxn_str, atp_coeff, adp_coeff in atp_synthesis_rxns[:10]:
        print(f"\n  {rxn_id}")
        print(f"    반응식: {rxn_str}")
        print(f"    ATP 생성 계수: {atp_coeff:+.2f}")
        print(f"    ADP 소비 계수: {adp_coeff:+.2f}")
    
    return atp_synthesis_rxns

def compare_atpm0_vs_atpm5(model):
    """ATPM=0 vs ATPM=5일 때 ADP/ATP 플럭스 비교"""
    print("\n" + "="*70)
    print("4. ATPM=0 vs ATPM=5일 때 ADP/ATP 플럭스 비교")
    print("="*70)
    
    atpm = model.reactions.get_by_id('ATPM')
    original_lb = atpm.lower_bound
    original_ub = atpm.upper_bound
    
    results = {}
    
    for atpm_val in [0, 5]:
        atpm.lower_bound = atpm_val
        atpm.upper_bound = 1000.0
        
        solution = model.optimize()
        
        # ATPM 플럭스
        atpm_flux = solution.fluxes.get('ATPM', 0.0)
        
        # ATP/ADP 관련 주요 반응 플럭스
        key_rxns = ['ATPM', 'ATPS4rpp', 'ATPS', 'Growth']
        key_fluxes = {}
        for rxn_id in key_rxns:
            if rxn_id in model.reactions:
                key_fluxes[rxn_id] = solution.fluxes.get(rxn_id, 0.0)
        
        results[atpm_val] = {
            'status': solution.status,
            'growth': solution.objective_value,
            'atpm_flux': atpm_flux,
            'key_fluxes': key_fluxes,
        }
        
        print(f"\n[ATPM = {atpm_val}]")
        print(f"  상태: {solution.status}")
        print(f"  성장률: {solution.objective_value:.6f}")
        print(f"  ATPM 플럭스: {atpm_flux:.6f}")
        print(f"  주요 반응 플럭스:")
        for rxn_id, flux in key_fluxes.items():
            print(f"    {rxn_id:15s}: {flux:10.6f}")
    
    # 원래 bounds 복원
    atpm.lower_bound = original_lb
    atpm.upper_bound = original_ub
    
    return results

def check_adp_generation_pathways(model):
    """ADP 생성 경로 확인"""
    print("\n" + "="*70)
    print("5. ADP 생성 경로 확인")
    print("="*70)
    
    # ADP를 생성하는 주요 반응들 찾기
    adp_producing_rxns = []
    
    for rxn in model.reactions:
        if 'adp_c' in [m.id for m in rxn.metabolites]:
            adp_coeff = rxn.metabolites.get(model.metabolites.get_by_id('adp_c'), 0)
            if adp_coeff > 0:  # ADP 생성
                adp_producing_rxns.append((rxn.id, rxn.reaction, adp_coeff))
    
    print(f"\n[ADP 생성 반응] (총 {len(adp_producing_rxns)}개)")
    print(f"  상위 20개만 표시:")
    
    for rxn_id, rxn_str, adp_coeff in adp_producing_rxns[:20]:
        # ATPM 반응 강조
        if rxn_id == 'ATPM':
            print(f"\n  *** {rxn_id} (ATPM) ***")
        else:
            print(f"\n  {rxn_id}")
        print(f"    반응식: {rxn_str}")
        print(f"    ADP 생성 계수: {adp_coeff:+.2f}")
    
    return adp_producing_rxns

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_with_BCAA_cofactors_ions_nad_transport.xml"
    
    print("="*70)
    print("Step 4-4: ATPM-ADP 원인 확인")
    print("="*70)
    
    # 모델 로드
    model = load_model(model_path)
    
    # 배지 조건 설정 (고정)
    if 'EX_ac_e' in model.reactions:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.lower_bound = -19.0
        ex_ac.upper_bound = -19.0
    
    if 'EX_o2_e' in model.reactions:
        ex_o2 = model.reactions.get_by_id('EX_o2_e')
        ex_o2.lower_bound = -100.0
        ex_o2.upper_bound = 1000.0
    
    # 1. ATPM 반응식 확인
    atpm = check_atpm_reaction(model)
    
    # 2. ATPM=0일 때 ADP/ATP 균형 확인
    balance_results = check_adp_atp_balance_at_atpm0(model)
    
    # 3. ATP 생성 경로에 ADP 필요 여부 확인
    atp_synth_rxns = check_atp_synthase_requirements(model)
    
    # 4. ATPM=0 vs ATPM=5 비교
    comparison_results = compare_atpm0_vs_atpm5(model)
    
    # 5. ADP 생성 경로 확인
    adp_producing_rxns = check_adp_generation_pathways(model)
    
    print("\n" + "="*70)
    print("결론")
    print("="*70)
    
    print("\n[핵심 발견]")
    print("1. ATPM 반응식 확인 완료")
    print("2. ATPM=0일 때 ADP/ATP 균형 분석 완료")
    print("3. ATP 생성 경로에 ADP 필요 여부 확인 완료")
    print("4. ATPM=0 vs ATPM=5 비교 완료")
    print("5. ADP 생성 경로 확인 완료")
    
    return model, {
        'atpm': atpm,
        'balance': balance_results,
        'atp_synth': atp_synth_rxns,
        'comparison': comparison_results,
        'adp_producing': adp_producing_rxns,
    }

if __name__ == "__main__":
    model, results = main()
