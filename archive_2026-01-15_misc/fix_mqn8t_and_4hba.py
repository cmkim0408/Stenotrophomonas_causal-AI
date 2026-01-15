#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MQN8t 및 4hba 문제 해결

문제:
1. MQN8t 차단 시 성장 불가 - 메나퀴논 생합성 경로 확인 필요
2. sink_4hba_c 제거 시 성장 불가 - 4hba 배출 경로 추가 필요
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

def check_menaquinone_synthesis(model, solution):
    """메나퀴논 생합성 경로 확인"""
    print("\n" + "="*70)
    print("메나퀴논 생합성 경로 확인")
    print("="*70)
    
    # mqn8_c 또는 mql8_c 관련 반응 찾기
    mqn_rxns = []
    for rxn in model.reactions:
        mets = [m.id for m in rxn.metabolites]
        if 'mqn8_c' in mets or 'mql8_c' in mets:
            mqn_rxns.append(rxn)
    
    print(f"\n[메나퀴논 관련 반응] (총 {len(mqn_rxns)}개)")
    for rxn in mqn_rxns[:10]:
        flux = solution.fluxes.get(rxn.id, 0.0)
        if abs(flux) > 1e-6:
            print(f"  {rxn.id:20s}: 플럭스 {flux:10.6f}")
            print(f"                     반응식: {rxn.reaction}")
    
    # mqn8_c 생산 반응 찾기
    mqn8_producing = []
    if 'mqn8_c' in model.metabolites:
        for rxn in model.reactions:
            if 'mqn8_c' in [m.id for m in rxn.metabolites]:
                mqn8_coeff = rxn.metabolites.get(model.metabolites.get_by_id('mqn8_c'), 0)
                if mqn8_coeff > 0:  # mqn8_c 생성
                    flux = solution.fluxes.get(rxn.id, 0.0)
                    if abs(flux) > 1e-6:
                        mqn8_producing.append((rxn.id, rxn.reaction, flux, mqn8_coeff))
    
    if mqn8_producing:
        print(f"\n[mqn8_c 생산 반응] (플럭스 > 1e-6, 총 {len(mqn8_producing)}개)")
        for rxn_id, rxn_eq, flux, coeff in mqn8_producing[:5]:
            mqn8_prod = flux * coeff
            print(f"  {rxn_id:20s}: 플럭스 {flux:10.6f}, mqn8_c 생산 {mqn8_prod:.6f}")
    else:
        print(f"\n[문제] mqn8_c 생산 반응이 작동하지 않습니다!")
    
    return mqn8_producing

def check_biomass_mql8_requirement(model):
    """Biomass의 mql8 요구 확인"""
    print("\n" + "="*70)
    print("Biomass의 mql8 요구 확인")
    print("="*70)
    
    if 'Growth' in model.reactions:
        growth = model.reactions.get_by_id('Growth')
        if 'mql8_c' in [m.id for m in growth.metabolites]:
            mql8_coeff = growth.metabolites.get(model.metabolites.get_by_id('mql8_c'), 0)
            print(f"\n[Growth 반응]")
            print(f"  mql8_c 계수: {mql8_coeff:.6f}")
            print(f"  [주의] Biomass가 mql8_c를 요구합니다.")
            print(f"         mqn8_c → mql8_c 전환 경로가 필요할 수 있습니다.")
            return mql8_coeff
        else:
            print(f"\n[INFO] Growth 반응에 mql8_c가 없습니다.")
            return 0.0
    return 0.0

def check_4hba_exchange_or_pathway(model):
    """4hba 배출 경로 확인"""
    print("\n" + "="*70)
    print("4hba 배출 경로 확인")
    print("="*70)
    
    if '4hba_c' not in model.metabolites:
        print("[INFO] 4hba_c가 모델에 없습니다.")
        return []
    
    # 4hba_c를 소비하는 반응 찾기 (sink 제외)
    hba_consuming = []
    for rxn in model.reactions:
        if '4hba_c' in [m.id for m in rxn.metabolites]:
            if 'sink' not in rxn.id.lower():
                hba_coeff = rxn.metabolites.get(model.metabolites.get_by_id('4hba_c'), 0)
                if hba_coeff < 0:  # 4hba_c 소비
                    hba_consuming.append((rxn.id, rxn.reaction, hba_coeff))
    
    if hba_consuming:
        print(f"\n[4hba_c 소비 반응] (sink 제외, 총 {len(hba_consuming)}개)")
        for rxn_id, rxn_eq, coeff in hba_consuming:
            print(f"  {rxn_id:20s}: 계수 {coeff:+.2f}")
            print(f"                     반응식: {rxn_eq}")
    else:
        print(f"\n[문제] 4hba_c를 소비하는 반응이 없습니다 (sink 제외)")
        print(f"       -> 4hba exchange 또는 다른 배출 경로 필요")
    
    # 4hba exchange 확인
    if 'EX_4hba_e' in model.reactions:
        ex_4hba = model.reactions.get_by_id('EX_4hba_e')
        print(f"\n[4hba Exchange]")
        print(f"  EX_4hba_e: bounds [{ex_4hba.lower_bound}, {ex_4hba.upper_bound}]")
    else:
        print(f"\n[INFO] EX_4hba_e가 모델에 없습니다.")
    
    return hba_consuming

def add_4hba_exchange(model):
    """4hba exchange 추가"""
    print("\n" + "="*70)
    print("4hba Exchange 추가")
    print("="*70)
    
    if '4hba_c' not in model.metabolites:
        print("[ERROR] 4hba_c가 모델에 없습니다.")
        return False
    
    # 4hba_e metabolite 확인/추가
    if '4hba_e' not in model.metabolites:
        hba_e = cobra.Metabolite('4hba_e', name='4-Hydroxybenzoate (extracellular)', compartment='e')
        model.add_metabolites([hba_e])
        print(f"[추가] 4hba_e metabolite 추가")
    
    # EX_4hba_e 추가
    if 'EX_4hba_e' not in model.reactions:
        ex_4hba = cobra.Reaction('EX_4hba_e')
        ex_4hba.name = "4-Hydroxybenzoate exchange"
        ex_4hba.add_metabolites({model.metabolites.get_by_id('4hba_e'): -1})
        ex_4hba.lower_bound = 0.0  # 배출만 허용
        ex_4hba.upper_bound = 1000.0
        model.add_reactions([ex_4hba])
        print(f"[추가] EX_4hba_e 반응 추가")
        print(f"  반응식: {ex_4hba.reaction}")
        print(f"  bounds: [{ex_4hba.lower_bound}, {ex_4hba.upper_bound}] (배출만 허용)")
    
    # Transport 반응 확인/추가
    transport_id = 'T_4hba_e_to_c'
    if transport_id not in model.reactions:
        if '4hba_e' in model.metabolites and '4hba_c' in model.metabolites:
            transport = cobra.Reaction(transport_id)
            transport.name = "4-Hydroxybenzoate transport (e<->c)"
            transport.add_metabolites({
                model.metabolites.get_by_id('4hba_e'): -1,
                model.metabolites.get_by_id('4hba_c'): 1,
            })
            transport.lower_bound = -1000.0
            transport.upper_bound = 1000.0
            model.add_reactions([transport])
            print(f"[추가] {transport_id} 반응 추가")
            print(f"  반응식: {transport.reaction}")
    
    return True

def test_model_without_sink(model):
    """sink 없이 모델 테스트"""
    print("\n" + "="*70)
    print("sink_4hba_c 없이 모델 테스트")
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
    
    # 4hba 관련 플럭스 확인
    if 'EX_4hba_e' in model.reactions:
        ex_4hba_flux = solution.fluxes.get('EX_4hba_e', 0.0)
        print(f"  EX_4hba_e 플럭스: {ex_4hba_flux:.6f}")
    
    if 'sink_4hba_c' in model.reactions:
        sink_flux = solution.fluxes.get('sink_4hba_c', 0.0)
        print(f"  sink_4hba_c 플럭스: {sink_flux:.6f} (여전히 존재)")
    
    return solution

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_coa_synthesis_fixed.xml"
    output_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_final_cleaned.xml"
    
    print("="*70)
    print("MQN8t 및 4hba 문제 해결")
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
    print(f"\n[수정 전 FBA]")
    print(f"  상태: {solution_before.status}")
    print(f"  성장률: {solution_before.objective_value:.6f}")
    
    # 1. 메나퀴논 생합성 경로 확인
    mqn8_producing = check_menaquinone_synthesis(model, solution_before)
    mql8_coeff = check_biomass_mql8_requirement(model)
    
    # 2. 4hba 배출 경로 확인 및 추가
    hba_consuming = check_4hba_exchange_or_pathway(model)
    if not hba_consuming:
        add_4hba_exchange(model)
    
    # 3. sink_4hba_c 제거
    if 'sink_4hba_c' in model.reactions:
        model.remove_reactions(['sink_4hba_c'])
        print(f"\n[제거] sink_4hba_c 제거 완료")
    
    # 4. 수정 후 테스트
    solution_after = test_model_without_sink(model)
    
    # 5. MQN8t 차단 테스트
    if 'MQN8t' in model.reactions:
        mqn8t = model.reactions.get_by_id('MQN8t')
        original_lb = mqn8t.lower_bound
        original_ub = mqn8t.upper_bound
        
        mqn8t.lower_bound = 0.0
        mqn8t.upper_bound = 0.0
        
        solution_mqn_blocked = model.optimize()
        
        print(f"\n[MQN8t 차단 후]")
        print(f"  상태: {solution_mqn_blocked.status}")
        print(f"  성장률: {solution_mqn_blocked.objective_value:.6f}")
        
        can_grow_without_mqn8t = solution_mqn_blocked.objective_value > 1e-6
        
        if can_grow_without_mqn8t:
            print(f"  [OK] MQN8t 없이도 성장 가능")
        else:
            print(f"  [문제] MQN8t 없으면 성장 불가")
            # 원래 bounds 복원
            mqn8t.lower_bound = original_lb
            mqn8t.upper_bound = original_ub
            print(f"  [복원] MQN8t bounds 복원")
    else:
        can_grow_without_mqn8t = True
    
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
    
    print(f"\n[MQN8t]")
    if can_grow_without_mqn8t:
        print(f"  [OK] MQN8t 없이도 성장 가능")
        if 'MQN8t' in model.reactions:
            mqn8t = model.reactions.get_by_id('MQN8t')
            mqn8t.lower_bound = 0.0
            mqn8t.upper_bound = 0.0
            print(f"  [차단] MQN8t bounds를 [0.0, 0.0]으로 설정")
    else:
        print(f"  [주의] MQN8t 없으면 성장 불가 - 메나퀴논 생합성 경로 확인 필요")
    
    print(f"\n[4hba]")
    print(f"  [OK] sink_4hba_c 제거 및 4hba exchange 추가")
    
    if growth_after > 1e-6:
        print(f"\n[결론] 수정 후에도 성장 가능 (성장률: {growth_after:.6f})")
        
        # 모델 저장
        cobra.io.write_sbml_model(model, str(output_path))
        print(f"\n[모델 저장] {output_path}")
    else:
        print(f"\n[주의] 수정 후 성장 불가 (성장률: {growth_after:.6f})")
    
    return model, solution_before, solution_after

if __name__ == "__main__":
    model, sol_before, sol_after = main()
