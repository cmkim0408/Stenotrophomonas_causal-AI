#!/usr/bin/env python
"""
경로 Gap 찾기
각 경로에서 어떤 부분이 끊어졌는지 상세 확인
"""

import cobra
import pandas as pd
from cobra.flux_analysis import find_blocked_reactions

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    return model

def find_gap_in_pep_pathway(model):
    """PEP 생산 경로 Gap 찾기"""
    print("="*70)
    print("PEP 생산 경로 Gap 분석")
    print("="*70)
    
    # Acetate 설정
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    try:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.lower_bound = -100
        ex_ac.upper_bound = 1000
    except KeyError:
        pass
    
    # 필수 영양소
    essentials = ['EX_nh4_e', 'EX_h2o_e', 'EX_h_e', 'EX_pi_e', 'EX_so4_e',
                  'EX_k_e', 'EX_na1_e', 'EX_mg2_e', 'EX_ca2_e', 'EX_fe2_e',
                  'EX_mn2_e', 'EX_zn2_e', 'EX_co2_e', 'EX_o2_e']
    
    for ex_id in essentials:
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            if ex_id in ['EX_co2_e', 'EX_o2_e']:
                ex_rxn.lower_bound = -1000
                ex_rxn.upper_bound = 1000
            else:
                ex_rxn.lower_bound = -1000
        except KeyError:
            pass
    
    # Acetate -> Acetyl-CoA 경로 확인
    print("\n[1] Acetate -> Acetyl-CoA 경로:")
    print("-" * 70)
    
    acetate_pathway = ['EX_ac_e', 'ACt', 'ACt2rpp', 'SUCOAACTr', 'ACS']
    
    for rxn_id in acetate_pathway:
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            print(f"  [OK] {rxn_id}: {rxn.name}")
            print(f"    {rxn.reaction}")
        except KeyError:
            print(f"  [MISSING] {rxn_id}")
    
    # Acetyl-CoA -> Pyruvate 경로 확인
    print("\n[2] Acetyl-CoA -> Pyruvate 경로:")
    print("-" * 70)
    
    # Pyruvate 생산 테스트
    try:
        pyr_c = model.metabolites.get_by_id('pyr_c')
        
        dm_pyr = cobra.Reaction('DM_pyr_c')
        dm_pyr.name = 'Pyruvate demand'
        dm_pyr.lower_bound = 0
        dm_pyr.upper_bound = 1000
        dm_pyr.add_metabolites({pyr_c: -1})
        model.add_reactions([dm_pyr])
        
        model.objective = dm_pyr.id
        solution = model.optimize()
        
        if solution.status == 'optimal' and solution.objective_value > 1e-6:
            print(f"  [OK] Pyruvate 생산 가능: {solution.objective_value:.6f}")
        else:
            print(f"  [FAIL] Pyruvate 생산 불가능 ({solution.status})")
        
        model.remove_reactions([dm_pyr])
        
    except KeyError:
        print("  [ERROR] pyr_c metabolite 없음")
    
    # Pyruvate -> PEP 경로 확인
    print("\n[3] Pyruvate -> PEP 경로:")
    print("-" * 70)
    
    pep_pathway = {
        'PPS': 'Phosphoenolpyruvate synthase (ATP 사용)',
        'PEPCK': 'Phosphoenolpyruvate carboxykinase',
        'PC': 'Pyruvate carboxylase (PEPCK과 함께 사용)'
    }
    
    for rxn_id, description in pep_pathway.items():
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            print(f"  [OK] {rxn_id}: {rxn.name}")
            print(f"    {rxn.reaction}")
            
            # 이 반응의 대사물질 확인
            reactants = [str(m) for m, c in rxn.metabolites.items() if c < 0]
            products = [str(m) for m, c in rxn.metabolites.items() if c > 0]
            
            print(f"    반응물: {reactants}")
            print(f"    생성물: {products}")
            
        except KeyError:
            print(f"  [MISSING] {rxn_id}: {description}")
    
    # PEP 생산 차단 원인 찾기
    print("\n[4] PEP 생산 차단 원인 분석:")
    print("-" * 70)
    
    # PPS가 존재하는데 왜 작동 안하는지 확인
    try:
        pps = model.reactions.get_by_id('PPS')
        print(f"  PPS 반응식: {pps.reaction}")
        print(f"  PPS 경계 조건: LB={pps.lower_bound}, UB={pps.upper_bound}")
        
        # PPS 반응물 확인
        pps_reactants = {m: abs(c) for m, c in pps.metabolites.items() if c < 0}
        print(f"  PPS 필요 반응물: {pps_reactants}")
        
        # 각 반응물 생산 가능 여부 확인
        for met, coeff in pps_reactants.items():
            met_id = met.id
            try:
                # 생산 가능 여부 테스트
                test_rxn = cobra.Reaction(f'TEST_{met_id}')
                test_rxn.add_metabolites({met: 1})
                test_rxn.lower_bound = 0
                test_rxn.upper_bound = 1000
                
                model.add_reactions([test_rxn])
                model.objective = test_rxn.id
                
                sol = model.optimize()
                model.remove_reactions([test_rxn])
                
                can_produce = sol.status == 'optimal' and sol.objective_value > 1e-6
                
                if can_produce:
                    print(f"    [OK] {met_id} 생산 가능")
                else:
                    print(f"    [FAIL] {met_id} 생산 불가능 -> PPS 차단!")
                    
            except Exception as e:
                print(f"    [ERROR] {met_id} 테스트 실패: {e}")
                
    except KeyError:
        print("  PPS 반응 없음")

def find_gap_in_g6p_pathway(model):
    """G6P 생산 경로 Gap 찾기"""
    print("\n" + "="*70)
    print("G6P 생산 경로 Gap 분석")
    print("="*70)
    
    # PEP -> G6P 경로 확인
    gluconeogenesis_steps = [
        ('PEP', 'pep_c'),
        ('2PG', '2pg_c'),
        ('3PG', '3pg_c'),
        ('G3P', 'g3p_c'),
        ('FBP', 'fdp_c'),
        ('F6P', 'f6p_c'),
        ('G6P', 'g6p_c')
    ]
    
    print("\nGluconeogenesis 단계별 대사물질 생산 가능 여부:")
    print("-" * 70)
    
    # 모든 exchange 초기화
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    try:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.lower_bound = -100
        ex_ac.upper_bound = 1000
    except KeyError:
        pass
    
    # 필수 영양소
    essentials = ['EX_nh4_e', 'EX_h2o_e', 'EX_h_e', 'EX_pi_e', 'EX_so4_e',
                  'EX_k_e', 'EX_na1_e', 'EX_mg2_e', 'EX_ca2_e', 'EX_fe2_e',
                  'EX_mn2_e', 'EX_zn2_e', 'EX_co2_e', 'EX_o2_e']
    
    for ex_id in essentials:
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            if ex_id in ['EX_co2_e', 'EX_o2_e']:
                ex_rxn.lower_bound = -1000
                ex_rxn.upper_bound = 1000
            else:
                ex_rxn.lower_bound = -1000
        except KeyError:
            pass
    
    blocked_metabolites = []
    
    for step_name, met_id in gluconeogenesis_steps:
        try:
            met = model.metabolites.get_by_id(met_id)
            
            # 생산 가능 여부 테스트
            test_rxn = cobra.Reaction(f'TEST_{met_id}')
            test_rxn.add_metabolites({met: 1})
            test_rxn.lower_bound = 0
            test_rxn.upper_bound = 1000
            
            model.add_reactions([test_rxn])
            model.objective = test_rxn.id
            
            sol = model.optimize()
            model.remove_reactions([test_rxn])
            
            can_produce = sol.status == 'optimal' and sol.objective_value > 1e-6
            
            if can_produce:
                print(f"  [OK] {step_name} ({met_id}): 생산 가능")
            else:
                print(f"  [FAIL] {step_name} ({met_id}): 생산 불가능 -> Gap!")
                blocked_metabolites.append((step_name, met_id))
                
                # 이전 단계 확인
                if gluconeogenesis_steps.index((step_name, met_id)) > 0:
                    prev_step, prev_met_id = gluconeogenesis_steps[gluconeogenesis_steps.index((step_name, met_id)) - 1]
                    print(f"    이전 단계 {prev_step} ({prev_met_id}) 확인 필요")
        
        except KeyError:
            print(f"  [ERROR] {step_name} ({met_id}): metabolite 없음")
            blocked_metabolites.append((step_name, met_id))
    
    return blocked_metabolites

def find_gap_in_e4p_pathway(model):
    """E4P 생산 경로 Gap 찾기"""
    print("\n" + "="*70)
    print("E4P 생산 경로 Gap 분석 (PPP)")
    print("="*70)
    
    # PPP 경로 단계
    ppp_steps = [
        ('G6P', 'g6p_c'),
        ('6PGL', '6pgl_c'),
        ('6PGC', '6pgc_c'),
        ('Ru5P', 'ru5p__D_c'),
        ('R5P', 'r5p_c'),
        ('Xu5P', 'xu5p__D_c'),
        ('S7P', 's7p_c'),
        ('E4P', 'e4p_c')
    ]
    
    print("\nPPP 단계별 대사물질 생산 가능 여부:")
    print("-" * 70)
    
    # 모든 exchange 초기화
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    try:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.lower_bound = -100
        ex_ac.upper_bound = 1000
    except KeyError:
        pass
    
    # 필수 영양소
    essentials = ['EX_nh4_e', 'EX_h2o_e', 'EX_h_e', 'EX_pi_e', 'EX_so4_e',
                  'EX_k_e', 'EX_na1_e', 'EX_mg2_e', 'EX_ca2_e', 'EX_fe2_e',
                  'EX_mn2_e', 'EX_zn2_e', 'EX_co2_e', 'EX_o2_e']
    
    for ex_id in essentials:
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            if ex_id in ['EX_co2_e', 'EX_o2_e']:
                ex_rxn.lower_bound = -1000
                ex_rxn.upper_bound = 1000
            else:
                ex_rxn.lower_bound = -1000
        except KeyError:
            pass
    
    blocked_metabolites = []
    
    for step_name, met_id in ppp_steps:
        try:
            met = model.metabolites.get_by_id(met_id)
            
            # 생산 가능 여부 테스트
            test_rxn = cobra.Reaction(f'TEST_{met_id}')
            test_rxn.add_metabolites({met: 1})
            test_rxn.lower_bound = 0
            test_rxn.upper_bound = 1000
            
            model.add_reactions([test_rxn])
            model.objective = test_rxn.id
            
            sol = model.optimize()
            model.remove_reactions([test_rxn])
            
            can_produce = sol.status == 'optimal' and sol.objective_value > 1e-6
            
            if can_produce:
                print(f"  [OK] {step_name} ({met_id}): 생산 가능")
            else:
                print(f"  [FAIL] {step_name} ({met_id}): 생산 불가능 -> Gap!")
                blocked_metabolites.append((step_name, met_id))
        
        except KeyError:
            print(f"  [ERROR] {step_name} ({met_id}): metabolite 없음")
            blocked_metabolites.append((step_name, met_id))
    
    return blocked_metabolites

def main():
    print("="*70)
    print("경로 Gap 찾기")
    print("="*70)
    
    model = load_model("BaseModel.xml")
    
    # 1. PEP 경로 Gap
    find_gap_in_pep_pathway(model)
    
    # 2. G6P 경로 Gap
    blocked_g6p = find_gap_in_g6p_pathway(model)
    
    # 3. E4P 경로 Gap
    blocked_e4p = find_gap_in_e4p_pathway(model)
    
    # 요약
    print("\n" + "="*70)
    print("Gap 요약")
    print("="*70)
    
    if blocked_g6p:
        print("\nG6P 경로 차단 지점:")
        for step, met_id in blocked_g6p:
            print(f"  - {step} ({met_id})")
    
    if blocked_e4p:
        print("\nE4P 경로 차단 지점:")
        for step, met_id in blocked_e4p:
            print(f"  - {step} ({met_id})")
    
    print("\n" + "="*70)

if __name__ == "__main__":
    main()
