#!/usr/bin/env python
"""
종합 부트스트랩 분석
무제한 영양소 상태와 비교하여 부족한 요소 확인
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

def compare_unlimited_vs_acetate(model, biomass_rxn):
    """무제한 영양소 vs Acetate 상태 비교"""
    print("="*70)
    print("무제한 영양소 vs Acetate 상태 비교")
    print("="*70)
    
    # 1. 무제한 영양소 상태
    print("\n[1] 무제한 영양소 상태")
    print("-" * 70)
    
    for rxn in model.exchanges:
        rxn.lower_bound = -1000
        rxn.upper_bound = 1000
    
    model.objective = biomass_rxn.id
    solution_unlimited = model.optimize()
    
    if solution_unlimited.status == 'optimal':
        biomass_unlimited = solution_unlimited.objective_value
        print(f"  [SUCCESS] Biomass flux: {biomass_unlimited:.6f} 1/h")
        
        # 활성 exchange 확인
        active_exchanges = []
        for rxn in model.exchanges:
            flux = solution_unlimited.fluxes.get(rxn.id, 0)
            if abs(flux) > 1:
                active_exchanges.append((rxn.id, flux))
        
        print(f"\n  활성 Exchange (절대값 > 1):")
        for rxn_id, flux in sorted(active_exchanges, key=lambda x: abs(x[1]), reverse=True)[:20]:
            direction = "Uptake" if flux < 0 else "Secretion"
            print(f"    {rxn_id}: {flux:.6f} ({direction})")
    
    # 2. Acetate + ATP 부트스트랩 상태
    print("\n[2] Acetate + ATP 부트스트랩 상태")
    print("-" * 70)
    
    # 모든 exchange 초기화
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    # Acetate 설정
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
    
    # ATP 부트스트랩
    try:
        atp_c = model.metabolites.get_by_id('atp_c')
        dm_atp = cobra.Reaction('DM_atp_c')
        dm_atp.name = 'ATP bootstrap'
        dm_atp.lower_bound = -0.1
        dm_atp.upper_bound = 1000
        dm_atp.add_metabolites({atp_c: -1})
        model.add_reactions([dm_atp])
    except KeyError:
        pass
    
    model.objective = biomass_rxn.id
    solution_acetate = model.optimize()
    
    if solution_acetate.status == 'optimal':
        biomass_acetate = solution_acetate.objective_value
        print(f"  Biomass flux: {biomass_acetate:.6f} 1/h")
        
        if biomass_acetate > 1e-6:
            print(f"  [SUCCESS] 성장 가능!")
        else:
            print(f"  [FAIL] 성장 불가능")
            
            # 차단 원인 분석: Biomass 구성 요소 생산 불가능 확인
            print("\n  [분석] Biomass 구성 요소 생산 불가능 확인:")
            
            # 주요 Biomass 구성 요소 확인
            biomass_components = list(biomass_rxn.metabolites.keys())
            sorted_components = sorted(biomass_components, 
                                     key=lambda m: abs(biomass_rxn.metabolites[m]), 
                                     reverse=True)
            
            # 상위 30개 구성 요소 확인
            missing_components = []
            for met in sorted_components[:30]:
                met_coeff = abs(biomass_rxn.metabolites[met])
                if met_coeff < 0.001:  # 너무 작은 계수는 스킵
                    continue
                
                met_id = met.id
                
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
                
                if not can_produce and met_coeff > 0.01:  # 중요한 구성 요소만
                    missing_components.append({
                        'Metabolite_ID': met_id,
                        'Coefficient': met_coeff,
                        'Status': sol.status
                    })
                    print(f"    [FAIL] {met_id} (계수: {met_coeff:.6f}): 생산 불가능 ({sol.status})")
            
            if missing_components:
                print(f"\n  총 {len(missing_components)}개 중요한 구성 요소가 생산 불가능")
    
    # 부트스트랩 제거
    try:
        model.remove_reactions(['DM_atp_c'])
    except:
        pass
    
    return solution_unlimited, solution_acetate

def test_comprehensive_bootstrap(model, biomass_rxn):
    """종합 부트스트랩 테스트 (여러 조합)"""
    print("\n" + "="*70)
    print("종합 부트스트랩 테스트")
    print("="*70)
    
    # 모든 exchange 초기화
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    # Acetate 설정
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
    
    # 다양한 부트스트랩 조합
    bootstrap_combinations = [
        {'atp_c': -0.1},
        {'atp_c': -0.1, 'coa_c': -0.01},
        {'atp_c': -0.1, 'coa_c': -0.01, 'nad_c': -0.01},
        {'atp_c': -0.1, 'coa_c': -0.01, 'nad_c': -0.01, 'gtp_c': -0.01},
        {'atp_c': -0.1, 'coa_c': -0.01, 'nad_c': -0.01, 'gtp_c': -0.01, 'utp_c': -0.01, 'ctp_c': -0.01},
        {'atp_c': -0.5, 'coa_c': -0.01, 'nad_c': -0.1},  # 더 많은 ATP/NAD
    ]
    
    results = []
    
    for i, bootstrap in enumerate(bootstrap_combinations, 1):
        print(f"\n[조합 {i}] {bootstrap}")
        
        # 부트스트랩 demand 추가
        demand_rxns = []
        for met_id, supply_rate in bootstrap.items():
            try:
                met = model.metabolites.get_by_id(met_id)
                dm_rxn = cobra.Reaction(f'DM_{met_id}_bs{i}')
                dm_rxn.name = f'{met_id} bootstrap'
                dm_rxn.lower_bound = supply_rate
                dm_rxn.upper_bound = 1000
                dm_rxn.add_metabolites({met: -1})
                model.add_reactions([dm_rxn])
                demand_rxns.append(dm_rxn.id)
            except KeyError:
                pass
        
        # Biomass 생산 테스트
        model.objective = biomass_rxn.id
        solution = model.optimize()
        
        result = {
            'Combination': str(bootstrap),
            'Status': solution.status,
            'Biomass_Flux': solution.objective_value if solution.status == 'optimal' else 0,
            'Can_Grow': solution.status == 'optimal' and solution.objective_value > 1e-6
        }
        
        if result['Can_Grow']:
            print(f"  [SUCCESS] Biomass flux: {result['Biomass_Flux']:.6f} 1/h")
        else:
            print(f"  [FAIL] Biomass flux: {result['Biomass_Flux']:.6f} 1/h ({solution.status})")
        
        results.append(result)
        
        # 부트스트랩 제거
        model.remove_reactions(demand_rxns)
    
    return results

def main():
    print("="*70)
    print("종합 부트스트랩 분석")
    print("="*70)
    
    model = load_model("BaseModel.xml")
    
    # Biomass 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass: {biomass_rxn.id}")
    
    # 1. 무제한 vs Acetate 비교
    solution_unlimited, solution_acetate = compare_unlimited_vs_acetate(model, biomass_rxn)
    
    # 2. 종합 부트스트랩 테스트
    bootstrap_results = test_comprehensive_bootstrap(model, biomass_rxn)
    
    # 결과 저장
    if bootstrap_results:
        df_results = pd.DataFrame(bootstrap_results)
        df_results.to_csv('comprehensive_bootstrap_results.csv', index=False)
        print(f"\n[OK] 종합 부트스트랩 결과 저장: comprehensive_bootstrap_results.csv")
    
    # 최종 요약
    print("\n" + "="*70)
    print("최종 요약")
    print("="*70)
    
    working_combinations = [r for r in bootstrap_results if r['Can_Grow']]
    
    if working_combinations:
        print(f"\n[SUCCESS] {len(working_combinations)}개 부트스트랩 조합이 작동:")
        for r in working_combinations:
            print(f"  [OK] {r['Combination']}: Biomass = {r['Biomass_Flux']:.6f} 1/h")
    else:
        print("\n[FAIL] 모든 부트스트랩 조합이 실패")
        print("\n문제 분석:")
        print("  1. ATP 부트스트랩만으로는 Biomass 생산 불가능")
        print("  2. 추가 영양소 또는 경로가 필요할 수 있음")
        print("  3. Biomass 구성 요소 생산 경로 불완전 가능성")
    
    print("\n" + "="*70)

if __name__ == "__main__":
    main()
