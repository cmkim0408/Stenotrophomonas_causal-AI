#!/usr/bin/env python
"""
포도당 기반 생장 테스트
Glucose만으로 생장 가능 여부 확인
"""

import cobra
import pandas as pd

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    print(f"[OK] 모델 로드: {model.id}")
    return model

def find_biomass_reaction(model):
    for name in ['Growth', 'BIOMASS']:
        try:
            return model.reactions.get_by_id(name)
        except KeyError:
            continue
    return None

def test_glucose_growth_detailed(model, biomass_rxn):
    """포도당 기반 생장 상세 테스트"""
    print("="*70)
    print("포도당 기반 생장 테스트")
    print("="*70)
    
    # 모든 exchange 초기화
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    # 포도당 설정
    glucose_exchanges = ['EX_glc__D_e', 'EX_glc_e', 'EX_glucose_e']
    glucose_set = False
    
    for ex_id in glucose_exchanges:
        try:
            ex_glc = model.reactions.get_by_id(ex_id)
            ex_glc.lower_bound = -100
            ex_glc.upper_bound = 1000
            print(f"[OK] 포도당 설정: {ex_id} (-100 ~ 1000 mmol/gDCW/h)")
            glucose_set = True
            break
        except KeyError:
            continue
    
    if not glucose_set:
        # 패턴으로 찾기
        for rxn in model.exchanges:
            if 'glc' in rxn.id.lower() or 'glucose' in rxn.id.lower():
                rxn.lower_bound = -100
                rxn.upper_bound = 1000
                print(f"[OK] 포도당 설정: {rxn.id} (-100 ~ 1000 mmol/gDCW/h)")
                glucose_set = True
                break
    
    if not glucose_set:
        print("[ERROR] 포도당 exchange 반응을 찾을 수 없습니다!")
        return None, None
    
    # 필수 영양소 설정
    essentials = ['EX_nh4_e', 'EX_h2o_e', 'EX_h_e', 'EX_pi_e', 'EX_so4_e',
                  'EX_k_e', 'EX_na1_e', 'EX_mg2_e', 'EX_ca2_e', 'EX_fe2_e',
                  'EX_mn2_e', 'EX_zn2_e', 'EX_co2_e', 'EX_o2_e']
    
    print("\n필수 영양소 설정:")
    for ex_id in essentials:
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            if ex_id in ['EX_co2_e', 'EX_o2_e']:
                ex_rxn.lower_bound = -1000
                ex_rxn.upper_bound = 1000
            else:
                ex_rxn.lower_bound = -1000
            print(f"  [OK] {ex_id}")
        except KeyError:
            print(f"  [SKIP] {ex_id} 없음")
    
    # Biomass 최적화
    model.objective = biomass_rxn.id
    
    print("\n최적화 수행 중...")
    solution = model.optimize()
    
    print(f"\n최적화 상태: {solution.status}")
    
    if solution.status == 'optimal':
        biomass_flux = solution.fluxes.get(biomass_rxn.id, 0)
        
        if biomass_flux > 1e-6:
            print(f"[SUCCESS] 포도당으로 생장 가능!")
            print(f"  Biomass flux: {biomass_flux:.6f} 1/h")
            
            # Exchange 플럭스
            print("\nExchange 플럭스 (절대값 > 0.001):")
            exchange_fluxes = []
            for rxn in model.exchanges:
                flux = solution.fluxes.get(rxn.id, 0)
                if abs(flux) > 0.001:
                    exchange_fluxes.append({
                        'Exchange_ID': rxn.id,
                        'Flux': flux,
                        'Direction': 'Uptake' if flux < 0 else 'Secretion'
                    })
                    direction = "Uptake" if flux < 0 else "Secretion"
                    print(f"  {rxn.id}: {flux:.6f} ({direction})")
            
            # Yield 계산
            glucose_uptake = abs(solution.fluxes.get(glucose_exchanges[0] if glucose_set else 'EX_glc__D_e', 0))
            if glucose_uptake > 0:
                yield_biomass = biomass_flux / glucose_uptake
                print(f"\n성장 속도: {biomass_flux:.6f} 1/h")
                print(f"Glucose 섭취: {glucose_uptake:.6f} mmol/gDCW/h")
                print(f"Biomass yield: {yield_biomass:.6f} gDW/mmol glucose")
            
            # 주요 경로 플럭스 확인
            print("\n주요 경로 플럭스:")
            key_reactions = {
                'Glycolysis': ['HEX1', 'PGI', 'PFK', 'FBA', 'GAPD', 'PGK', 'PYK'],
                'TCA': ['PDH', 'CS', 'ACONT', 'ICDHx', 'AKGDH', 'SUCOAS', 'SUCD', 'FUM', 'MDH'],
                'ETC': ['NADH16', 'SUCD', 'QCR', 'CYO3', 'ATPS4rpp'],
                'Biomass': ['Growth']
            }
            
            for category, rxn_ids in key_reactions.items():
                print(f"\n  [{category}]:")
                for rxn_id in rxn_ids:
                    try:
                        flux = solution.fluxes.get(rxn_id, 0)
                        if abs(flux) > 1e-8:
                            print(f"    {rxn_id}: {flux:.6f}")
                    except KeyError:
                        pass
            
            return solution, True
        
        else:
            print(f"[FAIL] 포도당으로 생장 불가능")
            print(f"  Biomass flux: {biomass_flux:.6f} 1/h")
            
            # 왜 생장하지 못하는지 분석
            print("\n차단 원인 분석:")
            
            # ATP 생산 확인
            try:
                atp_c = model.metabolites.get_by_id('atp_c')
                
                # ATP 생산 테스트
                dm_atp = cobra.Reaction('DM_atp_c_test')
                dm_atp.add_metabolites({atp_c: -1})
                dm_atp.lower_bound = 0
                dm_atp.upper_bound = 1000
                model.add_reactions([dm_atp])
                
                model.objective = dm_atp.id
                atp_sol = model.optimize()
                model.remove_reactions([dm_atp])
                
                if atp_sol.status == 'optimal' and atp_sol.objective_value > 1e-6:
                    print("  [OK] ATP 생산 가능")
                else:
                    print("  [FAIL] ATP 생산 불가능 -> 생장 차단")
            except:
                pass
            
            return solution, False
    
    elif solution.status == 'infeasible':
        print(f"[FAIL] 최적화 실패: infeasible")
        print("\n차단 원인 분석:")
        
        # Infeeasible 원인 찾기
        print("  - 모델이 주어진 제약 조건 하에서 해를 찾을 수 없음")
        print("  - 필수 대사물질 생산 불가능 가능성")
        
        # 필수 영양소를 더 추가해보기
        print("\n추가 영양소 테스트:")
        
        # 소량의 뉴클레오티드 부트스트랩
        bootstrap_mets = {
            'atp_c': -0.1,
            'gtp_c': -0.01,
            'utp_c': -0.01,
            'ctp_c': -0.01
        }
        
        for met_id, supply_rate in bootstrap_mets.items():
            try:
                met = model.metabolites.get_by_id(met_id)
                dm_rxn = cobra.Reaction(f'DM_{met_id}_test')
                dm_rxn.add_metabolites({met: -1})
                dm_rxn.lower_bound = supply_rate
                dm_rxn.upper_bound = 1000
                model.add_reactions([dm_rxn])
            except KeyError:
                pass
        
        model.objective = biomass_rxn.id
        solution_bootstrap = model.optimize()
        
        if solution_bootstrap.status == 'optimal':
            biomass_bootstrap = solution_bootstrap.objective_value
            if biomass_bootstrap > 1e-6:
                print(f"  [SUCCESS] 부트스트랩 추가 후 생장 가능: {biomass_bootstrap:.6f} 1/h")
            else:
                print(f"  [FAIL] 부트스트랩 추가 후에도 생장 불가능")
        else:
            print(f"  [FAIL] 부트스트랩 추가 후에도 infeasible")
        
        # 부트스트랩 제거
        for met_id in bootstrap_mets.keys():
            try:
                model.remove_reactions([f'DM_{met_id}_test'])
            except:
                pass
        
        return solution, False
    
    else:
        print(f"[ERROR] 최적화 실패: {solution.status}")
        return solution, False

def compare_with_acetate(model, biomass_rxn):
    """Acetate와 비교"""
    print("\n" + "="*70)
    print("포도당 vs Acetate 비교")
    print("="*70)
    
    results = {
        'Substrate': [],
        'Status': [],
        'Biomass_Flux': [],
        'Can_Grow': []
    }
    
    substrates = [
        {
            'name': 'Glucose',
            'exchange_id': 'EX_glc__D_e',
            'uptake': -100
        },
        {
            'name': 'Acetate',
            'exchange_id': 'EX_ac_e',
            'uptake': -100
        }
    ]
    
    for substrate in substrates:
        # 모든 exchange 초기화
        for rxn in model.exchanges:
            rxn.upper_bound = 0
            rxn.lower_bound = 0
        
        # 해당 탄소원 설정
        try:
            ex_rxn = model.reactions.get_by_id(substrate['exchange_id'])
            ex_rxn.lower_bound = substrate['uptake']
            ex_rxn.upper_bound = 1000
        except KeyError:
            results['Substrate'].append(substrate['name'])
            results['Status'].append('Exchange not found')
            results['Biomass_Flux'].append(0)
            results['Can_Grow'].append(False)
            continue
        
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
        
        # 최적화
        model.objective = biomass_rxn.id
        solution = model.optimize()
        
        results['Substrate'].append(substrate['name'])
        results['Status'].append(solution.status)
        
        if solution.status == 'optimal':
            biomass_flux = solution.fluxes.get(biomass_rxn.id, 0)
            results['Biomass_Flux'].append(biomass_flux)
            results['Can_Grow'].append(biomass_flux > 1e-6)
            
            print(f"\n[{substrate['name']}]:")
            print(f"  상태: {solution.status}")
            print(f"  Biomass flux: {biomass_flux:.6f} 1/h")
            if biomass_flux > 1e-6:
                print(f"  [SUCCESS] 생장 가능")
            else:
                print(f"  [FAIL] 생장 불가능")
        else:
            results['Biomass_Flux'].append(0)
            results['Can_Grow'].append(False)
            print(f"\n[{substrate['name']}]:")
            print(f"  상태: {solution.status}")
            print(f"  [FAIL] 최적화 실패")
    
    return results

def main():
    print("="*70)
    print("포도당 기반 생장 가능 여부 확인")
    print("="*70)
    
    # 모델 로드
    model = load_model("BaseModel.xml")
    
    # Biomass 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass: {biomass_rxn.id}")
    
    # 포도당 생장 테스트
    solution, can_grow = test_glucose_growth_detailed(model, biomass_rxn)
    
    # Acetate와 비교
    comparison_results = compare_with_acetate(model, biomass_rxn)
    
    # 결과 저장
    if comparison_results:
        df_comparison = pd.DataFrame(comparison_results)
        df_comparison.to_csv('glucose_vs_acetate_growth.csv', index=False)
        print(f"\n[OK] 비교 결과 저장: glucose_vs_acetate_growth.csv")
    
    # 최종 결론
    print("\n" + "="*70)
    print("최종 결론")
    print("="*70)
    
    if can_grow:
        print("\n[SUCCESS] 포도당으로 생장 가능합니다!")
        print("  → 모델은 포도당을 탄소원으로 사용 가능")
        print("  → 문제는 Acetate 경로에 특정적으로 존재")
    else:
        print("\n[FAIL] 포도당으로도 생장 불가능합니다")
        print("  → 모델 구조 자체에 문제가 있을 수 있습니다")
        print("  → 또는 부트스트랩 문제가 포도당에서도 존재")
        
        if comparison_results:
            glucose_can_grow = comparison_results['Can_Grow'][0] if comparison_results['Can_Grow'] else False
            acetate_can_grow = comparison_results['Can_Grow'][1] if len(comparison_results['Can_Grow']) > 1 else False
            
            print("\n비교 결과:")
            print(f"  포도당: {'생장 가능' if glucose_can_grow else '생장 불가능'}")
            print(f"  Acetate: {'생장 가능' if acetate_can_grow else '생장 불가능'}")
    
    print("\n" + "="*70)

if __name__ == "__main__":
    main()
