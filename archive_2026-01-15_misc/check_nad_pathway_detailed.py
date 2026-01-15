#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NAD/NADP 생산 경로 상세 확인

경로:
1. nac_c + prpp_c → nicrnt_c (NAPRT)
2. nicrnt_c + atp_c + h_c → dnad_c (NNATr)
3. dnad_c → nad_c (NADS1/NADS2)
4. nad_c → nadp_c (NADK)

확인 사항:
- PRPP 생성 가능성
- nicrnt_c 생성 가능성
- dnad_c 생성 가능성
- 경로 연결 확인
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

def check_reaction_flux(model, rxn_id):
    """반응이 플럭스를 가질 수 있는지 확인"""
    if rxn_id not in model.reactions:
        return None, "reaction_not_in_model"
    
    rxn = model.reactions.get_by_id(rxn_id)
    
    # Objective를 이 반응으로 설정
    model.objective = rxn_id
    model.objective_direction = 'max'
    
    solution = model.optimize()
    
    if solution.status == 'optimal':
        max_flux = solution.objective_value
        # 역방향도 확인
        model.objective_direction = 'min'
        solution_min = model.optimize()
        min_flux = solution_min.objective_value if solution_min.status == 'optimal' else 0.0
        return max_flux, min_flux, "optimal"
    else:
        return 0.0, 0.0, solution.status

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

def check_nad_pathway_step_by_step(model):
    """NAD 생산 경로를 단계별로 확인"""
    print("\n" + "="*70)
    print("NAD 생산 경로 단계별 확인")
    print("="*70)
    
    # 경로 단계
    steps = [
        ('PRPP', 'prpp_c', 'PRPP (Phosphoribosyl pyrophosphate)'),
        ('nac_c', 'nac_c', 'Nicotinate (cytosol)'),
        ('NAPRT', 'nicrnt_c', 'Nicotinate ribonucleotide'),
        ('NNATr', 'dnad_c', 'Deamino-NAD'),
        ('NADS1/NADS2', 'nad_c', 'NAD'),
        ('NADK', 'nadp_c', 'NADP'),
    ]
    
    print(f"\n{'단계':<15} {'대사물질':<15} {'이름':<40} {'생산 가능성':<20}")
    print("-" * 90)
    
    for step_name, met_id, met_name in steps:
        if step_name.startswith('NAPRT') or step_name.startswith('NNATr') or step_name.startswith('NADS') or step_name.startswith('NADK'):
            # 반응 ID
            if step_name == 'NAPRT':
                rxn_id = 'NAPRT'
            elif step_name == 'NNATr':
                rxn_id = 'NNATr'
            elif step_name == 'NADS1/NADS2':
                rxn_id = 'NADS1'  # NADS1 확인
            elif step_name == 'NADK':
                rxn_id = 'NADK'
            else:
                rxn_id = None
            
            if rxn_id and rxn_id in model.reactions:
                rxn = model.reactions.get_by_id(rxn_id)
                print(f"{step_name:<15} {met_id:<15} {met_name:<40}")
                print(f"             반응식: {rxn.reaction}")
                print(f"             bounds: [{rxn.lower_bound}, {rxn.upper_bound}]")
                
                # 반응 플럭스 확인
                max_flux, min_flux, status = check_reaction_flux(model, rxn_id)
                if max_flux is not None:
                    if abs(max_flux) > 1e-6 or abs(min_flux) > 1e-6:
                        print(f"             플럭스: [{min_flux:.6f}, {max_flux:.6f}] [OK]")
                    else:
                        print(f"             플럭스: [{min_flux:.6f}, {max_flux:.6f}] [NO]")
            else:
                print(f"{step_name:<15} {met_id:<15} {met_name:<40} [반응 없음]")
        
        # 대사물질 생산 가능성 확인
        if met_id in model.metabolites:
            max_prod, status = test_metabolite_production(model, met_id)
            if max_prod is not None:
                status_str = "[OK]" if max_prod > 1e-6 else "[NO]"
                print(f"             생산 가능성: {status_str} max_production = {max_prod:.6f}")
            else:
                print(f"             생산 가능성: [ERROR] {status}")
        else:
            print(f"             생산 가능성: [NOT IN MODEL]")
        
        print()

def check_prpp_production(model):
    """PRPP 생산 가능성 확인"""
    print("\n" + "="*70)
    print("PRPP 생산 가능성 확인")
    print("="*70)
    
    if 'prpp_c' not in model.metabolites:
        print("[ERROR] prpp_c가 모델에 없습니다.")
        return
    
    # PRPP 생성 반응 확인
    prpp_reactions = ['PRPPS', 'PRPPS1', 'PRPPS2']
    
    print("\n[PRPP 생성 반응]")
    for rxn_id in prpp_reactions:
        if rxn_id in model.reactions:
            rxn = model.reactions.get_by_id(rxn_id)
            print(f"  [OK] {rxn_id}: {rxn.reaction}")
            print(f"       bounds: [{rxn.lower_bound}, {rxn.upper_bound}]")
        else:
            print(f"  [MISSING] {rxn_id}")
    
    # PRPP 생산 테스트
    max_prod, status = test_metabolite_production(model, 'prpp_c')
    if max_prod is not None:
        status_str = "[OK]" if max_prod > 1e-6 else "[NO]"
        print(f"\n[PRPP 생산 가능성]")
        print(f"  {status_str} max_production = {max_prod:.6f}")
    else:
        print(f"\n[PRPP 생산 가능성]")
        print(f"  [ERROR] {status}")

def check_nicrnt_production(model):
    """nicrnt_c 생산 가능성 확인"""
    print("\n" + "="*70)
    print("nicrnt_c 생산 가능성 확인")
    print("="*70)
    
    if 'nicrnt_c' not in model.metabolites:
        print("[ERROR] nicrnt_c가 모델에 없습니다.")
        return
    
    # NAPRT 반응 확인
    if 'NAPRT' in model.reactions:
        rxn = model.reactions.get_by_id('NAPRT')
        print(f"\n[NAPRT 반응]")
        print(f"  반응식: {rxn.reaction}")
        print(f"  bounds: [{rxn.lower_bound}, {rxn.upper_bound}]")
        
        # 필요한 전구체 확인
        print(f"\n  필요한 전구체:")
        for met, coeff in rxn.metabolites.items():
            if coeff < 0:  # 소비되는 metabolite
                print(f"    {met.id}: 계수 {coeff}")
                max_prod, status = test_metabolite_production(model, met.id)
                if max_prod is not None:
                    status_str = "[OK]" if max_prod > 1e-6 else "[NO]"
                    print(f"      생산 가능성: {status_str} max_production = {max_prod:.6f}")
    
    # nicrnt_c 생산 테스트
    max_prod, status = test_metabolite_production(model, 'nicrnt_c')
    if max_prod is not None:
        status_str = "[OK]" if max_prod > 1e-6 else "[NO]"
        print(f"\n[nicrnt_c 생산 가능성]")
        print(f"  {status_str} max_production = {max_prod:.6f}")
    else:
        print(f"\n[nicrnt_c 생산 가능성]")
        print(f"  [ERROR] {status}")

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_with_BCAA_cofactors_ions_nad.xml"
    
    print("="*70)
    print("NAD/NADP 생산 경로 상세 확인")
    print("="*70)
    
    # 모델 로드
    print(f"\n[모델 로드]")
    model = cobra.io.read_sbml_model(str(model_path))
    print(f"[OK] 모델 로드 완료")
    
    # 배지 조건 설정
    model = setup_media_forced(model)
    
    # PRPP 생산 확인
    check_prpp_production(model)
    
    # nicrnt_c 생산 확인
    check_nicrnt_production(model)
    
    # NAD 생산 경로 단계별 확인
    check_nad_pathway_step_by_step(model)
    
    return model

if __name__ == "__main__":
    model = main()
