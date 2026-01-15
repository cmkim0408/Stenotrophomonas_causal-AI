#!/usr/bin/env python
"""
포도당 직접 공급 시 infeasible 원인 진단
포도당이 세포 내에 있어도 생장이 안 되는 이유 찾기
"""

import cobra
import pandas as pd

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    return model

def find_biomass_reaction(model):
    for name in ['Growth', 'BIOMASS']:
        try:
            return model.reactions.get_by_id(name)
        except KeyError:
            continue
    return None

def diagnose_infeasible_with_glucose_direct(model, biomass_rxn):
    """포도당 직접 공급 시 infeasible 원인 진단"""
    print("="*70)
    print("포도당 직접 공급 시 Infeasible 원인 진단")
    print("="*70)
    
    # 모든 exchange 초기화
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
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
    
    # 포도당 직접 공급
    try:
        glc__D_c = model.metabolites.get_by_id('glc__D_c')
        dm_glc = cobra.Reaction('DM_glc__D_c_direct')
        dm_glc.lower_bound = -100
        dm_glc.upper_bound = 1000
        dm_glc.add_metabolites({glc__D_c: -1})
        model.add_reactions([dm_glc])
    except KeyError:
        pass
    
    # 1. 포도당으로 주요 구성 요소 생산 가능 여부 확인
    print("\n[1] 포도당으로 주요 구성 요소 생산 가능 여부:")
    print("-" * 70)
    
    key_components = {
        'atp_c': 'ATP',
        'nad_c': 'NAD+',
        'nadh_c': 'NADH',
        'coa_c': 'CoA',
        'pep_c': 'PEP',
        'g6p_c': 'Glucose-6-phosphate',
        'pyr_c': 'Pyruvate',
        'accoa_c': 'Acetyl-CoA',
        'h_p': 'Proton motive force'
    }
    
    component_status = []
    
    for met_id, met_name in key_components.items():
        try:
            met = model.metabolites.get_by_id(met_id)
            
            test_rxn = cobra.Reaction(f'TEST_{met_id}')
            test_rxn.add_metabolites({met: 1})
            test_rxn.lower_bound = 0
            test_rxn.upper_bound = 1000
            model.add_reactions([test_rxn])
            model.objective = test_rxn.id
            
            solution = model.optimize()
            model.remove_reactions([test_rxn])
            
            can_produce = solution.status == 'optimal' and solution.objective_value > 1e-6
            
            component_status.append({
                'Component': met_name,
                'Metabolite_ID': met_id,
                'Can_Produce': can_produce,
                'Status': solution.status,
                'Max_Flux': solution.objective_value if solution.status == 'optimal' else 0
            })
            
            status_icon = "[OK]" if can_produce else "[FAIL]"
            print(f"  {status_icon} {met_name} ({met_id})")
            
            if not can_produce:
                print(f"      상태: {solution.status}")
                if solution.status == 'optimal':
                    print(f"      최대 플럭스: {solution.objective_value:.6f}")
            
        except KeyError:
            component_status.append({
                'Component': met_name,
                'Metabolite_ID': met_id,
                'Can_Produce': False,
                'Status': 'Metabolite missing',
                'Max_Flux': 0
            })
            print(f"  [MISSING] {met_name} ({met_id})")
    
    # 2. 순차적 부트스트랩 테스트
    print("\n[2] 순차적 부트스트랩 테스트:")
    print("-" * 70)
    
    bootstrap_strategies = [
        {'name': 'ATP만', 'components': {'atp_c': -0.1}},
        {'name': 'NAD+만', 'components': {'nad_c': -0.1}},
        {'name': 'CoA만', 'components': {'coa_c': -0.01}},
        {'name': 'ATP + NAD+', 'components': {'atp_c': -0.1, 'nad_c': -0.1}},
        {'name': 'ATP + CoA', 'components': {'atp_c': -0.1, 'coa_c': -0.01}},
        {'name': 'NAD+ + CoA', 'components': {'nad_c': -0.1, 'coa_c': -0.01}},
        {'name': 'ATP + NAD+ + CoA', 'components': {'atp_c': -0.1, 'nad_c': -0.1, 'coa_c': -0.01}},
    ]
    
    bootstrap_results = []
    
    for strategy in bootstrap_strategies:
        # 부트스트랩 추가
        demand_rxns = []
        for met_id, supply_rate in strategy['components'].items():
            try:
                met = model.metabolites.get_by_id(met_id)
                dm_rxn = cobra.Reaction(f'DM_{met_id}_bs')
                dm_rxn.lower_bound = supply_rate
                dm_rxn.upper_bound = 1000
                dm_rxn.add_metabolites({met: -1})
                model.add_reactions([dm_rxn])
                demand_rxns.append(dm_rxn.id)
            except KeyError:
                pass
        
        # Biomass 최적화
        model.objective = biomass_rxn.id
        solution = model.optimize()
        
        result = {
            'Strategy': strategy['name'],
            'Components': str(strategy['components']),
            'Status': solution.status,
            'Biomass_Flux': solution.objective_value if solution.status == 'optimal' else 0,
            'Can_Grow': solution.status == 'optimal' and solution.objective_value > 1e-6
        }
        
        bootstrap_results.append(result)
        
        status_icon = "[OK]" if result['Can_Grow'] else "[FAIL]"
        print(f"  {status_icon} {strategy['name']}: {result['Status']}")
        
        if result['Can_Grow']:
            print(f"      Biomass flux: {result['Biomass_Flux']:.6f} 1/h")
        
        # 부트스트랩 제거
        model.remove_reactions(demand_rxns)
    
    # 3. 최소 부트스트랩 찾기
    print("\n[3] 최소 부트스트랩 찾기:")
    print("-" * 70)
    
    working_strategies = [r for r in bootstrap_results if r['Can_Grow']]
    
    if working_strategies:
        print(f"  [OK] 작동하는 부트스트랩 전략 {len(working_strategies)}개 발견:")
        for r in working_strategies:
            print(f"    - {r['Strategy']}: Biomass = {r['Biomass_Flux']:.6f} 1/h")
            print(f"      구성 요소: {r['Components']}")
    else:
        print("  [FAIL] 모든 부트스트랩 전략 실패")
        print("  → 더 많은 구성 요소 부트스트랩 필요 또는 다른 문제 존재")
    
    # 결과 저장
    if component_status:
        df_components = pd.DataFrame(component_status)
        df_components.to_csv('glucose_direct_component_status.csv', index=False)
        print(f"\n[OK] 구성 요소 생산 상태 저장: glucose_direct_component_status.csv")
    
    if bootstrap_results:
        df_bootstrap = pd.DataFrame(bootstrap_results)
        df_bootstrap.to_csv('glucose_direct_bootstrap_results.csv', index=False)
        print(f"[OK] 부트스트랩 결과 저장: glucose_direct_bootstrap_results.csv")
    
    return component_status, bootstrap_results

def main():
    print("="*70)
    print("포도당 직접 공급 시 Infeasible 원인 진단")
    print("="*70)
    
    model = load_model("BaseModel.xml")
    
    # Biomass 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass: {biomass_rxn.id}")
    
    # 진단 수행
    component_status, bootstrap_results = diagnose_infeasible_with_glucose_direct(model, biomass_rxn)
    
    # 최종 요약
    print("\n" + "="*70)
    print("최종 요약")
    print("="*70)
    
    if component_status:
        can_produce_count = sum(1 for s in component_status if s['Can_Produce'])
        total_count = len(component_status)
        
        print(f"\n포도당으로 구성 요소 생산 가능 여부:")
        print(f"  생산 가능: {can_produce_count}/{total_count}")
        
        if can_produce_count < total_count:
            print("\n생산 불가능한 구성 요소:")
            for s in component_status:
                if not s['Can_Produce']:
                    print(f"  - {s['Component']} ({s['Metabolite_ID']}): {s['Status']}")
    
    if bootstrap_results:
        working = [r for r in bootstrap_results if r['Can_Grow']]
        
        if working:
            print(f"\n부트스트랩으로 생장 가능한 전략:")
            for r in working:
                print(f"  [OK] {r['Strategy']}: {r['Components']}")
        else:
            print(f"\n부트스트랩으로도 생장 불가능")
            print("  → 포도당이 세포 내에 있어도 다른 근본적인 문제 존재")
    
    print("\n" + "="*70)

if __name__ == "__main__":
    main()
