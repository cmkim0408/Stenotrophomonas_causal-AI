#!/usr/bin/env python
"""
Proton Motive Force (PMF) 생성 경로 확인
h_p 생성 경로 및 전자 전달 사슬 확인
"""

import cobra
import pandas as pd

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    return model

def check_proton_motive_force(model):
    """Proton Motive Force 생성 경로 확인"""
    print("="*70)
    print("1. Proton Motive Force (PMF) 생성 경로 확인")
    print("="*70)
    
    # h_p metabolite 확인
    try:
        h_p = model.metabolites.get_by_id('h_p')
        print(f"\n[OK] h_p metabolite 존재")
        print(f"  총 반응 수: {len(h_p.reactions)}개")
        
        # h_p 생성 반응 확인
        h_p_producing = [rxn for rxn in h_p.reactions if h_p in rxn.products]
        h_p_consuming = [rxn for rxn in h_p.reactions if h_p in rxn.reactants]
        
        print(f"  h_p 생성 반응: {len(h_p_producing)}개")
        print(f"  h_p 소비 반응: {len(h_p_consuming)}개")
        
        print("\n주요 h_p 생성 반응:")
        for rxn in h_p_producing[:10]:
            print(f"  - {rxn.id}: {rxn.name}")
            print(f"    {rxn.reaction}")
        
    except KeyError:
        print("\n[ERROR] h_p metabolite 없음")
        return None
    
    # 전자 전달 사슬 (ETC) 반응 확인
    print("\n" + "="*70)
    print("2. 전자 전달 사슬 (ETC) 반응 확인")
    print("="*70)
    
    etc_reactions = {
        'NADH16': 'NADH dehydrogenase (Complex I)',
        'SUCD': 'Succinate dehydrogenase (Complex II)',
        'CYO3': 'Cytochrome c oxidase (Complex IV)',
        'ATPS4rpp': 'ATP synthase (Complex V)',
        'QCR': 'Cytochrome bc1 complex (Complex III)'
    }
    
    print("\nETC Complex 반응:")
    etc_status = []
    
    for rxn_id, rxn_name in etc_reactions.items():
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            etc_status.append({
                'Complex': rxn_name,
                'Reaction_ID': rxn_id,
                'Exists': True,
                'Equation': rxn.reaction,
                'GPR': rxn.gene_reaction_rule
            })
            print(f"  [OK] {rxn_id}: {rxn_name}")
            print(f"    {rxn.reaction}")
            
            # h_p 관련 확인
            if 'h_p' in rxn.reaction:
                h_p_coeff = rxn.metabolites.get(h_p, 0)
                if h_p_coeff > 0:
                    print(f"    [OK] h_p 생성: {h_p_coeff}")
                elif h_p_coeff < 0:
                    print(f"    [OK] h_p 소비: {abs(h_p_coeff)}")
        except KeyError:
            etc_status.append({
                'Complex': rxn_name,
                'Reaction_ID': rxn_id,
                'Exists': False,
                'Equation': '',
                'GPR': ''
            })
            print(f"  [MISSING] {rxn_id}: {rxn_name}")
    
    return etc_status, h_p

def test_etc_pathway(model, h_p):
    """ETC 경로 작동 테스트"""
    print("\n" + "="*70)
    print("3. ETC 경로 작동 테스트")
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
    
    # h_p 생성 테스트
    print("\n[테스트 1] h_p 생성 가능 여부 (Acetate 기반)")
    print("-" * 70)
    
    try:
        # h_p sink 추가
        try:
            dm_h_p = model.reactions.get_by_id('DM_h_p')
            dm_h_p.lower_bound = 0
            dm_h_p.upper_bound = 1000
        except KeyError:
            dm_h_p = cobra.Reaction('DM_h_p')
            dm_h_p.name = 'h_p demand'
            dm_h_p.lower_bound = 0
            dm_h_p.upper_bound = 1000
            dm_h_p.add_metabolites({h_p: -1})
            model.add_reactions([dm_h_p])
        
        model.objective = dm_h_p.id
        solution = model.optimize()
        
        if solution.status == 'optimal':
            h_p_flux = solution.objective_value
            
            if h_p_flux > 1e-6:
                print(f"  [SUCCESS] h_p 생산 가능: {h_p_flux:.6f} mmol/gDCW/h")
                
                # ETC 경로 플럭스 확인
                etc_rxns = ['NADH16', 'SUCD', 'QCR', 'CYO3', 'ATPS4rpp']
                print("\n  ETC 경로 플럭스:")
                for rxn_id in etc_rxns:
                    try:
                        flux = solution.fluxes.get(rxn_id, 0)
                        if abs(flux) > 1e-8:
                            print(f"    {rxn_id}: {flux:.6f}")
                    except KeyError:
                        pass
                
                # NADH/NADH2 경로 확인
                try:
                    nadh_c = model.metabolites.get_by_id('nadh_c')
                    nadh_flux = solution.fluxes.get('NADH16', 0)
                    if abs(nadh_flux) > 1e-8:
                        print(f"\n  NADH16 플럭스: {nadh_flux:.6f}")
                except:
                    pass
                
                model.remove_reactions([dm_h_p])
                return True, solution
            else:
                print(f"  [FAIL] h_p 생산 불가능 (flux: {h_p_flux:.6f})")
                model.remove_reactions([dm_h_p])
                return False, solution
        else:
            print(f"  [FAIL] 최적화 실패: {solution.status}")
            try:
                model.remove_reactions([dm_h_p])
            except:
                pass
            return False, solution
    
    except Exception as e:
        print(f"  [ERROR] 테스트 실패: {e}")
        return False, None

def test_atp_via_etc(model):
    """ETC를 통한 ATP 생산 테스트"""
    print("\n" + "="*70)
    print("4. ETC를 통한 ATP 생산 테스트")
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
    
    # NADH 생산 확인
    print("\n[테스트 1] NADH 생산 가능 여부")
    print("-" * 70)
    
    try:
        nadh_c = model.metabolites.get_by_id('nadh_c')
        
        dm_nadh = cobra.Reaction('DM_nadh_c')
        dm_nadh.name = 'NADH demand'
        dm_nadh.lower_bound = 0
        dm_nadh.upper_bound = 1000
        dm_nadh.add_metabolites({nadh_c: -1})
        model.add_reactions([dm_nadh])
        
        model.objective = dm_nadh.id
        solution = model.optimize()
        
        if solution.status == 'optimal' and solution.objective_value > 1e-6:
            print(f"  [SUCCESS] NADH 생산 가능: {solution.objective_value:.6f}")
            
            # NADH 생산 경로 확인
            nadh_producing = ['ICDHx', 'AKGDH', 'MDH', 'GAPD']
            print("\n  NADH 생산 경로:")
            for rxn_id in nadh_producing:
                try:
                    flux = solution.fluxes.get(rxn_id, 0)
                    if abs(flux) > 1e-8:
                        print(f"    {rxn_id}: {flux:.6f}")
                except KeyError:
                    pass
        else:
            print(f"  [FAIL] NADH 생산 불가능 ({solution.status})")
        
        model.remove_reactions([dm_nadh])
        
    except Exception as e:
        print(f"  [ERROR] NADH 테스트 실패: {e}")
    
    # ATP 생산 (ETC 경로 통한)
    print("\n[테스트 2] ETC 경로를 통한 ATP 생산")
    print("-" * 70)
    
    try:
        atp_c = model.metabolites.get_by_id('atp_c')
        
        dm_atp = cobra.Reaction('DM_atp_c')
        dm_atp.name = 'ATP demand'
        dm_atp.lower_bound = 0
        dm_atp.upper_bound = 1000
        dm_atp.add_metabolites({atp_c: -1})
        model.add_reactions([dm_atp])
        
        model.objective = dm_atp.id
        solution = model.optimize()
        
        if solution.status == 'optimal':
            atp_flux = solution.objective_value
            
            if atp_flux > 1e-6:
                print(f"  [SUCCESS] ATP 생산 가능: {atp_flux:.6f} mmol/gDCW/h")
                
                # ATP 생산 경로 확인
                print("\n  ATP 생산 경로 플럭스:")
                
                # ETC 경로
                etc_rxns = ['NADH16', 'SUCD', 'QCR', 'CYO3', 'ATPS4rpp']
                etc_active = False
                for rxn_id in etc_rxns:
                    try:
                        flux = solution.fluxes.get(rxn_id, 0)
                        if abs(flux) > 1e-8:
                            print(f"    {rxn_id}: {flux:.6f}")
                            etc_active = True
                    except KeyError:
                        pass
                
                if etc_active:
                    print("    [OK] ETC 경로 활성화")
                else:
                    print("    [WARNING] ETC 경로 비활성")
                
                # 다른 ATP 생산 경로
                atp_producing_rxns = ['PYK', 'PPAKr']
                for rxn_id in atp_producing_rxns:
                    try:
                        flux = solution.fluxes.get(rxn_id, 0)
                        if abs(flux) > 1e-8:
                            print(f"    {rxn_id}: {flux:.6f}")
                    except KeyError:
                        pass
                
                model.remove_reactions([dm_atp])
                return True, solution
            else:
                print(f"  [FAIL] ATP 생산 불가능 (flux: {atp_flux:.6f})")
                
                # 왜 안되는지 분석
                print("\n  차단 원인 분석:")
                
                # NADH 생산 가능 여부
                try:
                    nadh_test = cobra.Reaction('TEST_nadh')
                    nadh_test.add_metabolites({nadh_c: -1})
                    nadh_test.lower_bound = 0
                    nadh_test.upper_bound = 1000
                    model.add_reactions([nadh_test])
                    model.objective = nadh_test.id
                    nadh_sol = model.optimize()
                    model.remove_reactions([nadh_test])
                    
                    if nadh_sol.status == 'optimal' and nadh_sol.objective_value > 1e-6:
                        print("    [OK] NADH 생산 가능")
                    else:
                        print("    [FAIL] NADH 생산 불가능 -> ETC 차단")
                except:
                    pass
                
                # h_p 생산 가능 여부
                try:
                    h_p = model.metabolites.get_by_id('h_p')
                    h_p_test = cobra.Reaction('TEST_h_p')
                    h_p_test.add_metabolites({h_p: -1})
                    h_p_test.lower_bound = 0
                    h_p_test.upper_bound = 1000
                    model.add_reactions([h_p_test])
                    model.objective = h_p_test.id
                    h_p_sol = model.optimize()
                    model.remove_reactions([h_p_test])
                    
                    if h_p_sol.status == 'optimal' and h_p_sol.objective_value > 1e-6:
                        print("    [OK] h_p 생산 가능")
                    else:
                        print("    [FAIL] h_p 생산 불가능 -> ATP synthase 차단")
                except:
                    pass
                
                model.remove_reactions([dm_atp])
                return False, solution
        else:
            print(f"  [FAIL] 최적화 실패: {solution.status}")
            try:
                model.remove_reactions([dm_atp])
            except:
                pass
            return False, solution
    
    except Exception as e:
        print(f"  [ERROR] ATP 테스트 실패: {e}")
        return False, None

def analyze_bootstrap_strategy(model):
    """부트스트랩 전략 분석"""
    print("\n" + "="*70)
    print("5. 부트스트랩 전략 분석")
    print("="*70)
    
    # 필요한 부트스트랩 물질들
    bootstrap_candidates = {
        'ATP': {
            'metabolite': 'atp_c',
            'reason': 'PPS, 기타 ATP 필요 반응',
            'suggested_supply': -0.1
        },
        'PEP': {
            'metabolite': 'pep_c',
            'reason': 'PYK (ATP 생산), Gluconeogenesis',
            'suggested_supply': -0.1
        },
        'CoA': {
            'metabolite': 'coa_c',
            'reason': 'Acetate -> Acetyl-CoA',
            'suggested_supply': -0.01
        },
        'NAD+': {
            'metabolite': 'nad_c',
            'reason': 'NADH 생산, ETC 경로',
            'suggested_supply': -0.01
        }
    }
    
    print("\n부트스트랩 후보물질:")
    print("-" * 70)
    
    for name, info in bootstrap_candidates.items():
        print(f"\n{name}:")
        print(f"  대사물질: {info['metabolite']}")
        print(f"  이유: {info['reason']}")
        print(f"  제안 공급량: {info['suggested_supply']} mmol/gDCW/h")
    
    # 최소 부트스트랩 조합 테스트
    print("\n[테스트] 최소 부트스트랩 조합 테스트")
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
    
    # 부트스트랩 조합 테스트
    bootstrap_combinations = [
        {'atp_c': -0.1},
        {'pep_c': -0.1},
        {'atp_c': -0.1, 'pep_c': -0.1},
        {'atp_c': -0.05, 'coa_c': -0.01},
        {'atp_c': -0.05, 'coa_c': -0.01, 'nad_c': -0.01}
    ]
    
    print("\n부트스트랩 조합별 ATP 생산 테스트:")
    results = []
    
    for i, bootstrap in enumerate(bootstrap_combinations, 1):
        print(f"\n[조합 {i}] {bootstrap}")
        
        # 부트스트랩 demand 추가
        demand_rxns = []
        for met_id, supply_rate in bootstrap.items():
            try:
                met = model.metabolites.get_by_id(met_id)
                dm_rxn = cobra.Reaction(f'DM_{met_id}_bootstrap')
                dm_rxn.name = f'{met_id} bootstrap'
                dm_rxn.lower_bound = supply_rate
                dm_rxn.upper_bound = 1000
                dm_rxn.add_metabolites({met: -1})
                model.add_reactions([dm_rxn])
                demand_rxns.append(dm_rxn.id)
            except KeyError:
                print(f"  [SKIP] {met_id} metabolite 없음")
        
        # ATP 생산 테스트
        try:
            atp_c = model.metabolites.get_by_id('atp_c')
            dm_atp = cobra.Reaction('DM_atp_c')
            dm_atp.name = 'ATP demand'
            dm_atp.lower_bound = 0
            dm_atp.upper_bound = 1000
            dm_atp.add_metabolites({atp_c: -1})
            model.add_reactions([dm_atp])
            
            model.objective = dm_atp.id
            solution = model.optimize()
            
            if solution.status == 'optimal':
                atp_flux = solution.objective_value
                can_produce = atp_flux > 1e-6
                
                results.append({
                    'Combination': str(bootstrap),
                    'Status': solution.status,
                    'ATP_Flux': atp_flux,
                    'Can_Produce': can_produce
                })
                
                if can_produce:
                    print(f"  [SUCCESS] ATP 생산 가능: {atp_flux:.6f}")
                else:
                    print(f"  [FAIL] ATP 생산 불가능: {atp_flux:.6f}")
            else:
                results.append({
                    'Combination': str(bootstrap),
                    'Status': solution.status,
                    'ATP_Flux': 0,
                    'Can_Produce': False
                })
                print(f"  [FAIL] 최적화 실패: {solution.status}")
            
            model.remove_reactions([dm_atp.id])
        
        except Exception as e:
            print(f"  [ERROR] 테스트 실패: {e}")
            results.append({
                'Combination': str(bootstrap),
                'Status': 'error',
                'ATP_Flux': 0,
                'Can_Produce': False
            })
        
        # 부트스트랩 제거
        model.remove_reactions(demand_rxns)
    
    # 결과 저장
    if results:
        df_results = pd.DataFrame(results)
        df_results.to_csv('bootstrap_strategy_results.csv', index=False)
        print(f"\n[OK] 부트스트랩 전략 결과 저장: bootstrap_strategy_results.csv")
    
    return results

def main():
    print("="*70)
    print("Proton Motive Force 및 부트스트랩 전략 분석")
    print("="*70)
    
    model = load_model("BaseModel.xml")
    
    # 1. PMF 경로 확인
    etc_status, h_p = check_proton_motive_force(model)
    
    # 2. h_p 생성 테스트
    if h_p:
        can_produce_h_p, h_p_solution = test_etc_pathway(model, h_p)
    
    # 3. ETC를 통한 ATP 생산 테스트
    can_produce_atp_etc, atp_solution = test_atp_via_etc(model)
    
    # 4. 부트스트랩 전략 분석
    bootstrap_results = analyze_bootstrap_strategy(model)
    
    # 최종 요약
    print("\n" + "="*70)
    print("최종 요약")
    print("="*70)
    
    if etc_status:
        existing_complexes = sum(1 for s in etc_status if s['Exists'])
        print(f"\nETC Complex: {existing_complexes}/{len(etc_status)}개 존재")
    
    if h_p:
        print(f"\nPMF (h_p): {'생산 가능' if can_produce_h_p else '생산 불가능'}")
    
    print(f"\nATP (ETC 경로): {'생산 가능' if can_produce_atp_etc else '생산 불가능'}")
    
    if bootstrap_results:
        working_combinations = [r for r in bootstrap_results if r['Can_Produce']]
        if working_combinations:
            print(f"\n부트스트랩: {len(working_combinations)}개 조합이 작동")
            for r in working_combinations:
                print(f"  [OK] {r['Combination']}: ATP flux = {r['ATP_Flux']:.6f}")
        else:
            print("\n부트스트랩: 작동하는 조합 없음")
    
    print("\n" + "="*70)

if __name__ == "__main__":
    main()
