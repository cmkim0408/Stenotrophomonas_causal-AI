#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
최종 모델 정리

1. MQN8t를 0으로 막고 성장하는지 확인
2. sink_4hba_c 제거 또는 4hba 배출/활용 경로 추가
3. EX_pnto__R_e, EX_nac_e가 실제 배지 가정과 맞는지 확인
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

def test_mqn8t_blocked(model):
    """MQN8t를 0으로 막고 성장하는지 확인"""
    print("\n" + "="*70)
    print("MQN8t 차단 테스트")
    print("="*70)
    
    # ATPM=0 설정
    if 'ATPM' in model.reactions:
        atpm = model.reactions.get_by_id('ATPM')
        atpm.lower_bound = 0.0
        atpm.upper_bound = 0.0
    
    # 차단 전
    print("\n[MQN8t 차단 전]")
    solution_before = model.optimize()
    print(f"  상태: {solution_before.status}")
    print(f"  성장률: {solution_before.objective_value:.6f}")
    
    if 'MQN8t' in model.reactions:
        mqn8t_flux_before = solution_before.fluxes.get('MQN8t', 0.0)
        print(f"  MQN8t 플럭스: {mqn8t_flux_before:.6f}")
    
    # MQN8t 차단
    if 'MQN8t' in model.reactions:
        mqn8t = model.reactions.get_by_id('MQN8t')
        original_lb = mqn8t.lower_bound
        original_ub = mqn8t.upper_bound
        
        mqn8t.lower_bound = 0.0
        mqn8t.upper_bound = 0.0
        print(f"\n[MQN8t 차단] bounds: [{mqn8t.lower_bound}, {mqn8t.upper_bound}]")
        
        # 차단 후
        print("\n[MQN8t 차단 후]")
        solution_after = model.optimize()
        print(f"  상태: {solution_after.status}")
        print(f"  성장률: {solution_after.objective_value:.6f}")
        
        # 메나퀴논 생합성 경로 확인
        print(f"\n[메나퀴논 생합성 경로 확인]")
        mqn_synthesis_rxns = []
        for rxn in model.reactions:
            if 'mqn8_c' in [m.id for m in rxn.metabolites]:
                mqn_coeff = rxn.metabolites.get(model.metabolites.get_by_id('mqn8_c'), 0)
                if mqn_coeff > 0:  # mqn8_c 생성
                    flux = solution_after.fluxes.get(rxn.id, 0.0)
                    if abs(flux) > 1e-6:
                        mqn_synthesis_rxns.append((rxn.id, rxn.reaction, flux, mqn_coeff))
        
        if mqn_synthesis_rxns:
            print(f"  [찾은 반응] (플럭스 > 1e-6, 총 {len(mqn_synthesis_rxns)}개)")
            for rxn_id, rxn_eq, flux, coeff in mqn_synthesis_rxns[:5]:
                mqn_prod = flux * coeff
                print(f"    {rxn_id:20s}: 플럭스 {flux:10.6f}, mqn8_c 생산 {mqn_prod:.6f}")
        else:
            print(f"  [문제] 메나퀴논 생합성 경로가 작동하지 않습니다!")
        
        # 원래 bounds 복원
        mqn8t.lower_bound = original_lb
        mqn8t.upper_bound = original_ub
        
        return solution_before, solution_after, solution_after.objective_value > 1e-6
    else:
        print(f"\n[INFO] MQN8t가 모델에 없습니다.")
        return solution_before, solution_before, True

def check_sink_4hba(model, solution):
    """sink_4hba_c 확인"""
    print("\n" + "="*70)
    print("sink_4hba_c 확인")
    print("="*70)
    
    sink_id = 'sink_4hba_c'
    if sink_id not in model.reactions:
        # 다른 패턴으로 찾기
        for rxn in model.reactions:
            if 'sink' in rxn.id.lower() and '4hba' in rxn.id.lower():
                sink_id = rxn.id
                break
    
    if sink_id in model.reactions:
        sink = model.reactions.get_by_id(sink_id)
        sink_flux = solution.fluxes.get(sink_id, 0.0)
        print(f"\n[{sink_id}]")
        print(f"  반응식: {sink.reaction}")
        print(f"  플럭스: {sink_flux:.6f}")
        print(f"  bounds: [{sink.lower_bound}, {sink.upper_bound}]")
        print(f"  [문제] sink 반응은 대사물질을 무에서 생성/소멸시키는 비현실적 반응")
        return sink_id, sink
    else:
        print(f"\n[INFO] sink_4hba_c가 모델에 없습니다.")
        return None, None

def check_4hba_pathways(model):
    """4hba 배출/활용 경로 확인"""
    print("\n" + "="*70)
    print("4hba 배출/활용 경로 확인")
    print("="*70)
    
    if '4hba_c' not in model.metabolites:
        print("[INFO] 4hba_c가 모델에 없습니다.")
        return []
    
    # 4hba_c를 사용하는 반응 찾기
    hba_rxns = []
    for rxn in model.reactions:
        if '4hba_c' in [m.id for m in rxn.metabolites]:
            hba_coeff = rxn.metabolites.get(model.metabolites.get_by_id('4hba_c'), 0)
            hba_rxns.append((rxn.id, rxn.reaction, hba_coeff))
    
    if hba_rxns:
        print(f"\n[4hba_c 관련 반응] (총 {len(hba_rxns)}개)")
        for rxn_id, rxn_eq, coeff in hba_rxns[:10]:
            print(f"  {rxn_id:20s}: 계수 {coeff:+.2f}")
            print(f"                     반응식: {rxn_eq}")
    else:
        print(f"\n[INFO] 4hba_c 관련 반응을 찾을 수 없습니다.")
    
    return hba_rxns

def check_media_assumptions(model):
    """배지 가정 확인"""
    print("\n" + "="*70)
    print("배지 가정 확인")
    print("="*70)
    
    # 비타민/보조인자 exchange 확인
    vitamin_exchanges = {
        'EX_pnto__R_e': 'Pantothenate (B5)',
        'EX_nac_e': 'Nicotinate (B3)',
        'EX_ncam_e': 'Nicotinamide (B3)',
        'EX_thm_e': 'Thiamine (B1)',
        'EX_ribflv_e': 'Riboflavin (B2)',
        'EX_btn_e': 'Biotin (B7)',
        'EX_fol_e': 'Folate (B9)',
    }
    
    print(f"\n[비타민/보조인자 Exchange 확인]")
    for ex_id, desc in vitamin_exchanges.items():
        if ex_id in model.reactions:
            ex_rxn = model.reactions.get_by_id(ex_id)
            print(f"  {ex_id:20s}: bounds [{ex_rxn.lower_bound:8.1f}, {ex_rxn.upper_bound:8.1f}] ({desc})")
        else:
            print(f"  {ex_id:20s}: 모델에 없음 ({desc})")
    
    print(f"\n[배지 가정 해석]")
    print(f"  - Yeast extract 가정: 비타민 exchange가 열려있어야 함")
    print(f"  - Defined minimal medium: 비타민 생합성 경로가 완전해야 함")
    print(f"  - 현재 상태를 확인하고 적절히 조정 필요")

def remove_sink_4hba(model):
    """sink_4hba_c 제거"""
    print("\n" + "="*70)
    print("sink_4hba_c 제거")
    print("="*70)
    
    sink_id = 'sink_4hba_c'
    if sink_id not in model.reactions:
        # 다른 패턴으로 찾기
        for rxn in model.reactions:
            if 'sink' in rxn.id.lower() and '4hba' in rxn.id.lower():
                sink_id = rxn.id
                break
    
    if sink_id in model.reactions:
        sink = model.reactions.get_by_id(sink_id)
        print(f"\n[제거] {sink_id}")
        print(f"  반응식: {sink.reaction}")
        print(f"  이유: sink 반응은 비현실적 (무에서 생성/소멸)")
        
        model.remove_reactions([sink_id])
        print(f"[OK] {sink_id} 제거 완료")
        return True
    else:
        print(f"\n[INFO] sink_4hba_c가 모델에 없습니다.")
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
    key_rxns = ['EX_ac_e', 'ACS', 'CS', 'ICL', 'MALS', 'SUCDi', 'MDH', 'Growth']
    for rxn_id in key_rxns:
        if rxn_id in model.reactions:
            flux = solution.fluxes.get(rxn_id, 0.0)
            if abs(flux) > 1e-6:
                print(f"  {rxn_id:15s}: {flux:10.6f}")
    
    return solution

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_coa_synthesis_fixed.xml"
    output_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_final_cleaned.xml"
    
    print("="*70)
    print("최종 모델 정리")
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
    
    # 1. MQN8t 차단 테스트
    sol_mqn_before, sol_mqn_after, can_grow_without_mqn8t = test_mqn8t_blocked(model)
    
    # 2. sink_4hba_c 확인 및 제거
    sink_id, sink_rxn = check_sink_4hba(model, solution_before)
    hba_rxns = check_4hba_pathways(model)
    sink_removed = remove_sink_4hba(model)
    
    # 3. 배지 가정 확인
    check_media_assumptions(model)
    
    # 수정 후 FBA
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
    
    print(f"\n[MQN8t 차단 테스트]")
    if can_grow_without_mqn8t:
        print(f"  [OK] MQN8t 없이도 성장 가능 (성장률: {sol_mqn_after.objective_value:.6f})")
        print(f"       -> 메나퀴논 생합성 경로가 작동 중")
    else:
        print(f"  [문제] MQN8t 없으면 성장 불가")
        print(f"         -> (a) 퀴논 생합성 경로 추가 필요")
        print(f"         -> (b) biomass의 mql8 요구 재검토 필요")
    
    print(f"\n[sink_4hba_c]")
    if sink_removed:
        print(f"  [OK] sink_4hba_c 제거 완료")
    else:
        print(f"  [INFO] sink_4hba_c가 없거나 제거되지 않음")
    
    print(f"\n[배지 가정]")
    print(f"  EX_pnto__R_e: 판토텐산 exchange (yeast extract 가정)")
    print(f"  EX_nac_e, EX_ncam_e: NAD 전구체 exchange (yeast extract 가정)")
    print(f"  -> 실제 배지 조건에 맞게 조정 필요")
    
    if growth_after > 1e-6:
        print(f"\n[결론] 최종 정리 후에도 성장 가능 (성장률: {growth_after:.6f})")
        
        # MQN8t 차단 (성장 가능한 경우)
        if can_grow_without_mqn8t and 'MQN8t' in model.reactions:
            mqn8t = model.reactions.get_by_id('MQN8t')
            mqn8t.lower_bound = 0.0
            mqn8t.upper_bound = 0.0
            print(f"\n[MQN8t 차단] bounds를 [0.0, 0.0]으로 설정")
        
        # 모델 저장
        cobra.io.write_sbml_model(model, str(output_path))
        print(f"\n[모델 저장] {output_path}")
    else:
        print(f"\n[주의] 수정 후 성장 불가 (성장률: {growth_after:.6f})")
    
    return model, solution_before, solution_after

if __name__ == "__main__":
    model, sol_before, sol_after = main()
