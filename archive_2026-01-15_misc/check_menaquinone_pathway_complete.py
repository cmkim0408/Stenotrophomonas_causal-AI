#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
메나퀴논 생합성 경로 완전성 확인

MQN8t 없이도 성장하려면 메나퀴논 생합성 경로가 완전해야 함
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

def test_mqn8_production_without_mqn8t(model):
    """MQN8t 없이 mqn8_c 생산 가능성 테스트"""
    print("\n" + "="*70)
    print("MQN8t 없이 mqn8_c 생산 가능성 테스트")
    print("="*70)
    
    if 'mqn8_c' not in model.metabolites:
        print("[ERROR] mqn8_c가 모델에 없습니다.")
        return 0.0
    
    # MQN8t 차단
    if 'MQN8t' in model.reactions:
        mqn8t = model.reactions.get_by_id('MQN8t')
        mqn8t.lower_bound = 0.0
        mqn8t.upper_bound = 0.0
    
    with model:
        # Demand 반응 추가
        demand_id = 'DM_mqn8_c'
        if demand_id in model.reactions:
            model.remove_reactions([demand_id])
        
        demand_rxn = cobra.Reaction(demand_id)
        demand_rxn.add_metabolites({model.metabolites.get_by_id('mqn8_c'): -1})
        model.add_reactions([demand_rxn])
        
        model.objective = demand_id
        model.objective_direction = 'max'
        
        solution = model.optimize()
        
        max_prod = solution.objective_value if solution.status == 'optimal' else 0.0
        print(f"\n[결과]")
        print(f"  상태: {solution.status}")
        print(f"  mqn8_c 최대 생산량: {max_prod:.6f}")
        
        if max_prod > 1e-6:
            print(f"  [OK] MQN8t 없이도 mqn8_c를 생산할 수 있습니다.")
        else:
            print(f"  [문제] MQN8t 없이는 mqn8_c를 생산할 수 없습니다!")
            print(f"         -> 메나퀴논 생합성 경로가 불완전함")
        
        return max_prod

def find_menaquinone_synthesis_pathway(model):
    """메나퀴논 생합성 경로 찾기"""
    print("\n" + "="*70)
    print("메나퀴논 생합성 경로 찾기")
    print("="*70)
    
    # 일반적인 메나퀴논 생합성 경로
    # chorismate → isochorismate → 2-succinylbenzoate → ... → menaquinone
    
    pathway_keywords = ['CHOR', 'ICHS', 'MEN', 'DMK', 'UPP', 'OCTAPRENYL']
    
    pathway_rxns = []
    for rxn in model.reactions:
        rxn_id_upper = rxn.id.upper()
        rxn_name_upper = (rxn.name if rxn.name else '').upper()
        
        if any(keyword in rxn_id_upper or keyword in rxn_name_upper for keyword in pathway_keywords):
            if 'mqn' in rxn_id_upper.lower() or 'menaquinone' in rxn_name_upper.lower():
                pathway_rxns.append(rxn)
    
    if pathway_rxns:
        print(f"\n[찾은 메나퀴논 생합성 관련 반응] (총 {len(pathway_rxns)}개)")
        for rxn in pathway_rxns[:10]:
            print(f"  {rxn.id:20s}: {rxn.reaction}")
    else:
        print(f"\n[INFO] 메나퀴논 생합성 관련 반응을 찾을 수 없습니다.")
    
    return pathway_rxns

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_final_cleaned.xml"
    
    print("="*70)
    print("메나퀴논 생합성 경로 완전성 확인")
    print("="*70)
    
    # 모델 로드
    model = load_model(model_path)
    model = setup_media_forced(model)
    
    # ATPM=0 설정
    if 'ATPM' in model.reactions:
        atpm = model.reactions.get_by_id('ATPM')
        atpm.lower_bound = 0.0
        atpm.upper_bound = 0.0
    
    # mqn8_c 생산 가능성 테스트
    max_prod = test_mqn8_production_without_mqn8t(model)
    
    # 메나퀴논 생합성 경로 찾기
    pathway_rxns = find_menaquinone_synthesis_pathway(model)
    
    # 결론
    print("\n" + "="*70)
    print("결론")
    print("="*70)
    
    if max_prod < 1e-6:
        print(f"\n[문제] 메나퀴논 생합성 경로가 불완전합니다.")
        print(f"       MQN8t 없이는 mqn8_c를 생산할 수 없습니다.")
        print(f"\n[권장 사항]")
        print(f"  1. MQN8t를 유지 (플럭스가 매우 작아서 큰 문제는 아님)")
        print(f"  2. 또는 메나퀴논 생합성 경로를 완성")
        print(f"  3. 또는 biomass의 mql8 요구를 재검토")
    else:
        print(f"\n[OK] 메나퀴논 생합성 경로가 작동합니다.")
        print(f"     MQN8t를 차단해도 됩니다.")
    
    return model, max_prod

if __name__ == "__main__":
    model, max_prod = main()
