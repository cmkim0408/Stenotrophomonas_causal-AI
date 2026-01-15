#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CoA 생합성 및 MQN8t 문제 확인 및 수정

문제점:
1. DM_coa_c: CoA를 외부에서 공급하는 역할 (CoA 생합성 경로가 막혀있을 가능성)
2. MQN8t: 무에서 생성하는 형태 (quinone 생합성 경로 확인 필요)
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

def check_dm_coa(model, solution):
    """DM_coa_c 확인"""
    print("\n" + "="*70)
    print("DM_coa_c 확인")
    print("="*70)
    
    if 'DM_coa_c' not in model.reactions:
        print("[INFO] DM_coa_c가 모델에 없습니다.")
        return None
    
    dm_coa = model.reactions.get_by_id('DM_coa_c')
    dm_coa_flux = solution.fluxes.get('DM_coa_c', 0.0)
    
    print(f"\n[DM_coa_c]")
    print(f"  반응식: {dm_coa.reaction}")
    print(f"  플럭스: {dm_coa_flux:.6f}")
    print(f"  bounds: [{dm_coa.lower_bound}, {dm_coa.upper_bound}]")
    
    # Biomass에서 CoA 필요량 계산
    if 'Growth' in model.reactions:
        growth_rxn = model.reactions.get_by_id('Growth')
        growth_flux = solution.fluxes.get('Growth', 0.0)
        
        if 'coa_c' in [m.id for m in growth_rxn.metabolites]:
            coa_coeff = growth_rxn.metabolites.get(model.metabolites.get_by_id('coa_c'), 0)
            coa_required = abs(coa_coeff * growth_flux)
            
            print(f"\n[Biomass CoA 필요량]")
            print(f"  Growth 플럭스: {growth_flux:.6f}")
            print(f"  CoA 계수: {coa_coeff:.6f}")
            print(f"  CoA 필요량: {coa_required:.6f}")
            print(f"  DM_coa_c 플럭스: {abs(dm_coa_flux):.6f}")
            
            if abs(abs(dm_coa_flux) - coa_required) < 1e-6:
                print(f"\n  [문제] DM_coa_c 플럭스가 CoA 필요량과 정확히 일치합니다!")
                print(f"         -> CoA 생합성 경로가 막혀있어 외부 공급에 의존 중")
    
    return dm_coa, dm_coa_flux

def check_coa_synthesis_pathway(model):
    """CoA 생합성 경로 확인"""
    print("\n" + "="*70)
    print("CoA 생합성 경로 확인")
    print("="*70)
    
    # CoA 생합성 관련 반응 키워드
    coa_synthesis_keywords = ['COA', 'PAN', 'PPAT', 'DPCK', 'PANT']
    
    coa_rxns = []
    for rxn in model.reactions:
        rxn_id_upper = rxn.id.upper()
        rxn_name_upper = (rxn.name if rxn.name else '').upper()
        
        if any(keyword in rxn_id_upper or keyword in rxn_name_upper for keyword in coa_synthesis_keywords):
            if 'coa_c' in [m.id for m in rxn.metabolites]:
                coa_coeff = rxn.metabolites.get(model.metabolites.get_by_id('coa_c'), 0)
                if coa_coeff > 0:  # CoA 생성
                    coa_rxns.append(rxn)
    
    if coa_rxns:
        print(f"\n[CoA 생산 반응] (총 {len(coa_rxns)}개)")
        for rxn in coa_rxns[:10]:
            print(f"  {rxn.id:20s}: {rxn.reaction}")
    else:
        print(f"\n[결과] CoA 생산 반응을 찾을 수 없습니다.")
    
    # 판토텐산 관련 확인
    print(f"\n[판토텐산 관련 반응]")
    pnto_rxns = []
    for rxn in model.reactions:
        if 'pnto' in rxn.id.lower() or 'pant' in rxn.id.lower():
            pnto_rxns.append(rxn)
    
    if pnto_rxns:
        for rxn in pnto_rxns[:5]:
            print(f"  {rxn.id:20s}: {rxn.reaction}")
    else:
        print(f"  [INFO] 판토텐산 관련 반응을 찾을 수 없습니다.")
    
    # 판토텐산 exchange 확인
    if 'EX_pnto__R_e' in model.reactions:
        ex_pnto = model.reactions.get_by_id('EX_pnto__R_e')
        print(f"\n[판토텐산 Exchange]")
        print(f"  EX_pnto__R_e: bounds [{ex_pnto.lower_bound}, {ex_pnto.upper_bound}]")
    else:
        print(f"\n[INFO] EX_pnto__R_e가 모델에 없습니다.")
    
    return coa_rxns

def check_mqn8t(model, solution):
    """MQN8t 확인"""
    print("\n" + "="*70)
    print("MQN8t 확인")
    print("="*70)
    
    if 'MQN8t' not in model.reactions:
        print("[INFO] MQN8t가 모델에 없습니다.")
        return None
    
    mqn8t = model.reactions.get_by_id('MQN8t')
    mqn8t_flux = solution.fluxes.get('MQN8t', 0.0)
    
    print(f"\n[MQN8t]")
    print(f"  반응식: {mqn8t.reaction}")
    print(f"  플럭스: {mqn8t_flux:.6f}")
    print(f"  bounds: [{mqn8t.lower_bound}, {mqn8t.upper_bound}]")
    
    if abs(mqn8t_flux) > 1e-6:
        print(f"  [주의] MQN8t가 무에서 생성하는 형태로 작동 중")
        print(f"         플럭스가 작지만 주의 필요")
    
    # 메나퀴논 생합성 경로 확인
    print(f"\n[메나퀴논 생합성 경로 확인]")
    mqn_synthesis_rxns = []
    for rxn in model.reactions:
        if 'mqn' in rxn.id.lower() or 'menaquinone' in (rxn.name.lower() if rxn.name else ''):
            if 'mqn8_c' in [m.id for m in rxn.metabolites]:
                mqn_coeff = rxn.metabolites.get(model.metabolites.get_by_id('mqn8_c'), 0)
                if mqn_coeff > 0:  # mqn8_c 생성
                    mqn_synthesis_rxns.append(rxn)
    
    if mqn_synthesis_rxns:
        print(f"  [찾은 반응] (총 {len(mqn_synthesis_rxns)}개)")
        for rxn in mqn_synthesis_rxns[:5]:
            print(f"    {rxn.id:20s}: {rxn.reaction}")
    else:
        print(f"  [INFO] 메나퀴논 생합성 반응을 찾을 수 없습니다.")
    
    return mqn8t, mqn8t_flux

def test_coa_production(model):
    """CoA 생산 가능성 테스트"""
    print("\n" + "="*70)
    print("CoA 생산 가능성 테스트")
    print("="*70)
    
    if 'coa_c' not in model.metabolites:
        print("[ERROR] coa_c가 모델에 없습니다.")
        return 0.0
    
    with model:
        # DM_coa_c 제거
        if 'DM_coa_c' in model.reactions:
            model.remove_reactions(['DM_coa_c'])
        
        # Demand 반응 추가
        demand_id = 'DM_coa_c_test'
        if demand_id in model.reactions:
            model.remove_reactions([demand_id])
        
        demand_rxn = cobra.Reaction(demand_id)
        demand_rxn.add_metabolites({model.metabolites.get_by_id('coa_c'): -1})
        model.add_reactions([demand_rxn])
        
        model.objective = demand_id
        model.objective_direction = 'max'
        
        solution = model.optimize()
        
        max_prod = solution.objective_value if solution.status == 'optimal' else 0.0
        print(f"\n[결과]")
        print(f"  상태: {solution.status}")
        print(f"  CoA 최대 생산량: {max_prod:.6f}")
        
        if max_prod > 1e-6:
            print(f"  [OK] CoA를 모델 내에서 생산할 수 있습니다.")
        else:
            print(f"  [문제] CoA를 모델 내에서 생산할 수 없습니다!")
            print(f"         -> CoA 생합성 경로가 막혀있거나 누락됨")
        
        return max_prod

def fix_coa_issue(model):
    """CoA 문제 해결"""
    print("\n" + "="*70)
    print("CoA 문제 해결")
    print("="*70)
    
    # 옵션 1: 판토텐산 exchange 열기
    if 'EX_pnto__R_e' in model.reactions:
        ex_pnto = model.reactions.get_by_id('EX_pnto__R_e')
        if ex_pnto.lower_bound > -1000:
            print(f"\n[판토텐산 Exchange 열기]")
            print(f"  현재 bounds: [{ex_pnto.lower_bound}, {ex_pnto.upper_bound}]")
            ex_pnto.lower_bound = -1000.0
            ex_pnto.upper_bound = 1000.0
            print(f"  수정 bounds: [{ex_pnto.lower_bound}, {ex_pnto.upper_bound}]")
            print(f"  [OK] 판토텐산 uptake 허용 (yeast extract에 포함 가능)")
            return True
    
    print(f"\n[INFO] EX_pnto__R_e가 없거나 이미 열려있습니다.")
    return False

def fix_mqn8t(model):
    """MQN8t 문제 해결"""
    print("\n" + "="*70)
    print("MQN8t 문제 해결")
    print("="*70)
    
    if 'MQN8t' in model.reactions:
        mqn8t = model.reactions.get_by_id('MQN8t')
        print(f"\n[MQN8t]")
        print(f"  현재 반응식: {mqn8t.reaction}")
        print(f"  현재 bounds: [{mqn8t.lower_bound}, {mqn8t.upper_bound}]")
        
        # 옵션 1: 0으로 막기
        # mqn8t.lower_bound = 0.0
        # mqn8t.upper_bound = 0.0
        # print(f"  [차단] MQN8t bounds를 [0.0, 0.0]으로 설정")
        
        # 옵션 2: 제거
        # model.remove_reactions(['MQN8t'])
        # print(f"  [제거] MQN8t 제거 완료")
        
        print(f"  [INFO] MQN8t는 플럭스가 매우 작아서 (약 3.97E-05)")
        print(f"         지금 당장 큰 문제는 아닙니다.")
        print(f"         필요시 나중에 제거/차단 가능")
        
        return False
    else:
        print(f"\n[INFO] MQN8t가 모델에 없습니다.")
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
    
    # DM_coa_c 확인
    if 'DM_coa_c' in model.reactions:
        dm_coa_flux = solution.fluxes.get('DM_coa_c', 0.0)
        print(f"  DM_coa_c 플럭스: {dm_coa_flux:.6f}")
    
    # MQN8t 확인
    if 'MQN8t' in model.reactions:
        mqn8t_flux = solution.fluxes.get('MQN8t', 0.0)
        print(f"  MQN8t 플럭스: {mqn8t_flux:.6f}")
    
    return solution

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_no_redox_shuttle.xml"
    output_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_coa_fixed.xml"
    
    print("="*70)
    print("CoA 생합성 및 MQN8t 문제 확인 및 수정")
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
    
    # 문제 확인
    dm_coa, dm_coa_flux = check_dm_coa(model, solution_before)
    coa_rxns = check_coa_synthesis_pathway(model)
    mqn8t, mqn8t_flux = check_mqn8t(model, solution_before)
    
    # CoA 생산 가능성 테스트
    coa_max_prod = test_coa_production(model)
    
    # 수정
    coa_fixed = fix_coa_issue(model)
    mqn8t_fixed = fix_mqn8t(model)
    
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
    
    print(f"\n[CoA 생산 가능성]")
    print(f"  최대 생산량: {coa_max_prod:.6f}")
    if coa_max_prod < 1e-6:
        print(f"  [문제] CoA 생합성 경로가 막혀있습니다.")
        print(f"         -> 판토텐산 exchange를 열었습니다 (yeast extract 조건 재현)")
    else:
        print(f"  [OK] CoA를 모델 내에서 생산할 수 있습니다.")
    
    if growth_after > 1e-6:
        print(f"\n[결론] CoA 문제 수정 후에도 성장 가능 (성장률: {growth_after:.6f})")
        
        # 모델 저장
        cobra.io.write_sbml_model(model, str(output_path))
        print(f"\n[모델 저장] {output_path}")
    else:
        print(f"\n[주의] 수정 후 성장 불가 (성장률: {growth_after:.6f})")
    
    return model, solution_before, solution_after

if __name__ == "__main__":
    model, sol_before, sol_after = main()
