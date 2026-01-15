#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NAD/NADP 생산 가능성 최종 확인

nac transport 반응 추가 후 NAD/NADP가 생성되는지 확인
"""

import cobra
from pathlib import Path
import sys

def test_metabolite_production(model, metabolite_id):
    """특정 대사물질의 최대 생산량 테스트"""
    if metabolite_id not in model.metabolites:
        return None, "metabolite_not_in_model"
    
    with model:
        demand_id = f'DM_{metabolite_id}'
        if demand_id in model.reactions:
            model.remove_reactions([demand_id])
        
        demand_rxn = cobra.Reaction(demand_id)
        demand_rxn.add_metabolites({model.metabolites.get_by_id(metabolite_id): -1})
        model.add_reactions([demand_rxn])
        
        model.objective = demand_id
        model.objective_direction = 'max'
        
        solution = model.optimize()
        
        if solution.status == 'optimal':
            return solution.objective_value, "optimal"
        else:
            return 0.0, solution.status

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

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_with_BCAA_cofactors_ions_nad_transport.xml"
    
    print("="*70)
    print("NAD/NADP 생산 가능성 최종 확인")
    print("="*70)
    
    # 모델 로드
    print(f"\n[모델 로드]")
    model = cobra.io.read_sbml_model(str(model_path))
    print(f"[OK] 모델 로드 완료")
    
    # 배지 조건 설정
    model = setup_media_forced(model)
    
    # NAD 경로 단계별 확인
    metabolites_to_test = [
        ('nac_c', 'Nicotinate (cytosol)'),
        ('nicrnt_c', 'Nicotinate ribonucleotide'),
        ('dnad_c', 'Deamino-NAD'),
        ('nad_c', 'NAD'),
        ('nadp_c', 'NADP'),
    ]
    
    print("\n" + "="*70)
    print("NAD 생산 경로 단계별 생산 가능성 확인")
    print("="*70)
    print(f"\n{'단계':<20} {'대사물질':<15} {'생산 가능성':<20} {'max_production':<20}")
    print("-" * 75)
    
    results = {}
    for met_id, met_name in metabolites_to_test:
        max_prod, status = test_metabolite_production(model, met_id)
        if max_prod is not None:
            status_str = "[OK]" if max_prod > 1e-6 else "[NO]"
            results[met_id] = (max_prod, status_str)
            print(f"{met_name:<20} {met_id:<15} {status_str:<20} {max_prod:.6f}")
        else:
            results[met_id] = (0.0, "[ERROR]")
            print(f"{met_name:<20} {met_id:<15} {'[ERROR]':<20} {status}")
    
    print("\n" + "="*70)
    print("결과 요약")
    print("="*70)
    
    nad_ok = results.get('nad_c', (0.0, "[NO]"))[1] == "[OK]"
    nadp_ok = results.get('nadp_c', (0.0, "[NO]"))[1] == "[OK]"
    
    if nad_ok and nadp_ok:
        print("\n[SUCCESS] ✅ NAD와 NADP 모두 생성 가능합니다!")
    elif nad_ok:
        print("\n[PARTIAL] ⚠️  NAD는 생성 가능하지만 NADP는 생성 불가합니다.")
    else:
        print("\n[FAILED] ❌ NAD/NADP 생산 경로에 여전히 문제가 있습니다.")
    
    return model, results

if __name__ == "__main__":
    model, results = main()
