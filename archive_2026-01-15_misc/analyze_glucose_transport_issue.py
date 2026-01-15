#!/usr/bin/env python
"""
포도당 Transport 경로 문제 분석
왜 포도당이 세포 내로 들어가지 못하는지 확인
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

def analyze_glucose_transport(model):
    """포도당 Transport 경로 상세 분석"""
    print("="*70)
    print("포도당 Transport 경로 문제 분석")
    print("="*70)
    
    # 포도당 transport 반응들 확인
    glucose_transports = {
        'GLCabc': 'D-glucose transport via ABC system',
        'GLCabcpp': 'D-glucose transport via ABC system (periplasm)',
        'GLCpts': 'D-glucose transport via PEP:Pyr PTS',
        'GLCt2rpp': 'D-glucose reversible transport via proton symport (periplasm)',
        'GLCtex': 'Glucose transport via diffusion (extracellular to periplasm)'
    }
    
    print("\n포도당 Transport 반응 상세 분석:")
    print("-" * 70)
    
    transport_status = []
    
    for rxn_id, rxn_name in glucose_transports.items():
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            
            print(f"\n[{rxn_id}]: {rxn_name}")
            print(f"  반응식: {rxn.reaction}")
            
            # 반응물 확인
            reactants = {m: abs(c) for m, c in rxn.metabolites.items() if c < 0}
            products = {m: abs(c) for m, c in rxn.metabolites.items() if c > 0}
            
            print(f"  반응물: {[str(m) for m in reactants.keys()]}")
            print(f"  생성물: {[str(m) for m in products.keys()]}")
            
            # ATP 필요 여부 확인
            needs_atp = any('atp' in str(m).lower() for m in reactants.keys())
            needs_pep = any('pep' in str(m).lower() for m in reactants.keys())
            
            if needs_atp:
                print(f"  [WARNING] ATP 필요 -> 순환 의존성 가능")
            if needs_pep:
                print(f"  [NOTE] PEP 필요 (PTS 경로)")
            
            # 경계 조건 확인
            print(f"  경계 조건: LB={rxn.lower_bound}, UB={rxn.upper_bound}")
            print(f"  가역성: {rxn.reversibility}")
            
            transport_status.append({
                'Transport_ID': rxn_id,
                'Name': rxn_name,
                'Equation': rxn.reaction,
                'Needs_ATP': needs_atp,
                'Needs_PEP': needs_pep,
                'Lower_Bound': rxn.lower_bound,
                'Upper_Bound': rxn.upper_bound,
                'Reversible': rxn.reversibility
            })
            
        except KeyError:
            print(f"\n[MISSING] {rxn_id}: {rxn_name}")
            transport_status.append({
                'Transport_ID': rxn_id,
                'Name': rxn_name,
                'Equation': '',
                'Needs_ATP': False,
                'Needs_PEP': False,
                'Lower_Bound': 0,
                'Upper_Bound': 0,
                'Reversible': False
            })
    
    # PTS 경로 테스트 (ATP 필요 없음)
    print("\n[테스트 1] PTS 경로 테스트 (ATP 필요 없음):")
    print("-" * 70)
    
    # 모든 exchange 초기화
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    # 포도당 설정
    try:
        ex_glc = model.reactions.get_by_id('EX_glc__D_e')
        ex_glc.lower_bound = -100
        ex_glc.upper_bound = 1000
    except KeyError:
        pass
    
    # 최소 영양소만 (ATP, PEP 부트스트랩 없이)
    minimal = ['EX_nh4_e', 'EX_h2o_e', 'EX_pi_e', 'EX_o2_e', 'EX_co2_e']
    
    for ex_id in minimal:
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            if ex_id in ['EX_co2_e', 'EX_o2_e']:
                ex_rxn.lower_bound = -1000
                ex_rxn.upper_bound = 1000
            else:
                ex_rxn.lower_bound = -1000
        except KeyError:
            pass
    
    # G6P 생산 테스트 (PTS 경로 통한)
    try:
        g6p_c = model.metabolites.get_by_id('g6p_c')
        
        dm_g6p = cobra.Reaction('DM_g6p_c')
        dm_g6p.add_metabolites({g6p_c: -1})
        dm_g6p.lower_bound = 0
        dm_g6p.upper_bound = 1000
        model.add_reactions([dm_g6p])
        
        model.objective = dm_g6p.id
        solution = model.optimize()
        
        model.remove_reactions([dm_g6p])
        
        if solution.status == 'optimal' and solution.objective_value > 1e-6:
            print(f"  [SUCCESS] G6P 생산 가능 (PTS 경로): {solution.objective_value:.6f}")
            
            # PTS 경로 플럭스 확인
            try:
                glcpts_flux = solution.fluxes.get('GLCpts', 0)
                print(f"    GLCpts 플럭스: {glcpts_flux:.6f}")
            except:
                pass
        else:
            print(f"  [FAIL] G6P 생산 불가능 ({solution.status})")
            
            # PEP 생산 확인
            try:
                pep_c = model.metabolites.get_by_id('pep_c')
                dm_pep = cobra.Reaction('DM_pep_c')
                dm_pep.add_metabolites({pep_c: -1})
                dm_pep.lower_bound = 0
                dm_pep.upper_bound = 1000
                model.add_reactions([dm_pep])
                
                model.objective = dm_pep.id
                pep_sol = model.optimize()
                model.remove_reactions([dm_pep])
                
                if pep_sol.status == 'optimal' and pep_sol.objective_value > 1e-6:
                    print(f"    [OK] PEP 생산 가능: {pep_sol.objective_value:.6f}")
                else:
                    print(f"    [FAIL] PEP 생산 불가능 -> PTS 경로 차단")
            except:
                pass
    
    except KeyError:
        print("  [ERROR] g6p_c metabolite 없음")
    
    # ABC 경로 테스트 (ATP 필요)
    print("\n[테스트 2] ABC 경로 테스트 (ATP 필요):")
    print("-" * 70)
    
    # ATP 부트스트랩 추가
    try:
        atp_c = model.metabolites.get_by_id('atp_c')
        dm_atp = cobra.Reaction('DM_atp_c_test')
        dm_atp.add_metabolites({atp_c: -1})
        dm_atp.lower_bound = -0.1
        dm_atp.upper_bound = 1000
        model.add_reactions([dm_atp])
        print("  [OK] ATP 부트스트랩 추가: -0.1 mmol/gDCW/h")
    except KeyError:
        pass
    
    # G6P 생산 테스트
    try:
        g6p_c = model.metabolites.get_by_id('g6p_c')
        
        dm_g6p = cobra.Reaction('DM_g6p_c_abc')
        dm_g6p.add_metabolites({g6p_c: -1})
        dm_g6p.lower_bound = 0
        dm_g6p.upper_bound = 1000
        model.add_reactions([dm_g6p])
        
        model.objective = dm_g6p.id
        solution = model.optimize()
        
        model.remove_reactions([dm_g6p])
        
        if solution.status == 'optimal' and solution.objective_value > 1e-6:
            print(f"  [SUCCESS] G6P 생산 가능 (ABC 경로): {solution.objective_value:.6f}")
            
            # ABC 경로 플럭스 확인
            abc_rxns = ['GLCabc', 'GLCabcpp', 'GLCtex']
            for rxn_id in abc_rxns:
                try:
                    flux = solution.fluxes.get(rxn_id, 0)
                    if abs(flux) > 1e-8:
                        print(f"    {rxn_id}: {flux:.6f}")
                except:
                    pass
        else:
            print(f"  [FAIL] G6P 생산 불가능 ({solution.status})")
    
    except KeyError:
        pass
    
    # ATP 부트스트랩 제거
    try:
        model.remove_reactions(['DM_atp_c_test'])
    except:
        pass
    
    # 결과 저장
    if transport_status:
        df_transport = pd.DataFrame(transport_status)
        df_transport.to_csv('glucose_transport_analysis.csv', index=False)
        print(f"\n[OK] 포도당 Transport 분석 저장: glucose_transport_analysis.csv")
    
    return transport_status

def test_minimal_bootstrap_for_glucose(model, biomass_rxn):
    """포도당 생장을 위한 최소 부트스트랩 테스트"""
    print("\n" + "="*70)
    print("포도당 생장을 위한 최소 부트스트랩 테스트")
    print("="*70)
    
    # 모든 exchange 초기화
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    # 포도당 설정
    try:
        ex_glc = model.reactions.get_by_id('EX_glc__D_e')
        ex_glc.lower_bound = -100
        ex_glc.upper_bound = 1000
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
        {'name': 'ATP만', 'components': {'atp_c': -0.1}},
        {'name': 'PEP만', 'components': {'pep_c': -0.1}},
        {'name': 'ATP + PEP', 'components': {'atp_c': -0.1, 'pep_c': -0.1}},
        {'name': 'ATP + CoA', 'components': {'atp_c': -0.1, 'coa_c': -0.01}},
        {'name': 'ATP + NAD+', 'components': {'atp_c': -0.1, 'nad_c': -0.01}},
        {'name': 'ATP + CoA + NAD+', 'components': {'atp_c': -0.1, 'coa_c': -0.01, 'nad_c': -0.01}},
    ]
    
    results = []
    
    for strategy in bootstrap_combinations:
        print(f"\n[{strategy['name']}]: {strategy['components']}")
        
        # 부트스트랩 추가
        demand_rxns = []
        for met_id, supply_rate in strategy['components'].items():
            try:
                met = model.metabolites.get_by_id(met_id)
                dm_rxn = cobra.Reaction(f'DM_{met_id}_glc')
                dm_rxn.name = f'{met_id} bootstrap'
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
        
        if result['Can_Grow']:
            print(f"  [SUCCESS] Biomass flux: {result['Biomass_Flux']:.6f} 1/h")
        else:
            print(f"  [FAIL] Biomass flux: {result['Biomass_Flux']:.6f} 1/h ({solution.status})")
        
        results.append(result)
        
        # 부트스트랩 제거
        model.remove_reactions(demand_rxns)
    
    # 결과 저장
    if results:
        df_results = pd.DataFrame(results)
        df_results.to_csv('glucose_bootstrap_results.csv', index=False)
        print(f"\n[OK] 포도당 부트스트랩 결과 저장: glucose_bootstrap_results.csv")
    
    return results

def main():
    print("="*70)
    print("포도당 Transport 및 생장 문제 분석")
    print("="*70)
    
    model = load_model("BaseModel.xml")
    
    # Biomass 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass: {biomass_rxn.id}")
    
    # 1. 포도당 Transport 분석
    transport_status = analyze_glucose_transport(model)
    
    # 2. 포도당 생장을 위한 최소 부트스트랩 테스트
    bootstrap_results = test_minimal_bootstrap_for_glucose(model, biomass_rxn)
    
    # 최종 요약
    print("\n" + "="*70)
    print("최종 요약")
    print("="*70)
    
    print("\n포도당 Transport:")
    abc_needs_atp = sum(1 for t in transport_status if t.get('Needs_ATP', False))
    pts_needs_pep = sum(1 for t in transport_status if t.get('Needs_PEP', False))
    
    print(f"  - ABC 경로 (ATP 필요): {abc_needs_atp}개")
    print(f"  - PTS 경로 (PEP 필요): {pts_needs_pep}개")
    
    if bootstrap_results:
        working_strategies = [r for r in bootstrap_results if r['Can_Grow']]
        
        if working_strategies:
            print(f"\n부트스트랩 성공:")
            for r in working_strategies:
                print(f"  [OK] {r['Strategy']}: Biomass = {r['Biomass_Flux']:.6f} 1/h")
        else:
            print(f"\n부트스트랩 실패:")
            print(f"  모든 부트스트랩 조합이 실패")
            print(f"  → 포도당 생장은 부트스트랩만으로는 해결되지 않음")
    
    print("\n" + "="*70)

if __name__ == "__main__":
    main()
