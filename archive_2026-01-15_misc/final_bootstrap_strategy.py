#!/usr/bin/env python
"""
최종 부트스트랩 전략 수립 및 검증
NADH 생산 경로 문제 확인 포함
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

def analyze_nadh_production_issue(model):
    """NADH 생산 경로 문제 분석"""
    print("="*70)
    print("NADH 생산 경로 문제 분석")
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
    
    print("\n[분석 1] NADH 생산 반응 확인:")
    print("-" * 70)
    
    try:
        nadh_c = model.metabolites.get_by_id('nadh_c')
        nadh_producing_rxns = [rxn for rxn in nadh_c.reactions if nadh_c in rxn.products]
        
        print(f"NADH 생산 반응: {len(nadh_producing_rxns)}개")
        
        # 주요 NADH 생산 반응
        key_nadh_producers = ['ICDHx', 'AKGDH', 'MDH', 'GAPD']
        
        for rxn_id in key_nadh_producers:
            try:
                rxn = model.reactions.get_by_id(rxn_id)
                print(f"\n  {rxn_id}: {rxn.name}")
                print(f"    {rxn.reaction}")
                
                # 반응물 확인
                reactants = {m: abs(c) for m, c in rxn.metabolites.items() if c < 0}
                print(f"    필요 반응물: {[str(m) for m in reactants.keys()]}")
                
                # 각 반응물 생산 가능 여부 간단 확인
                for met in reactants.keys():
                    met_id = met.id
                    # 간단한 패턴 확인 (실제 생산 테스트는 복잡하므로 생략)
                    if met_id in ['nad_c', 'nadp_c']:
                        print(f"      {met_id}: [확인 필요]")
                
            except KeyError:
                print(f"  [MISSING] {rxn_id}")
        
    except KeyError:
        print("[ERROR] nadh_c metabolite 없음")
    
    # NAD+ 생산 가능 여부 확인
    print("\n[분석 2] NAD+ 생산 가능 여부:")
    print("-" * 70)
    
    try:
        nad_c = model.metabolites.get_by_id('nad_c')
        
        dm_nad = cobra.Reaction('DM_nad_c')
        dm_nad.name = 'NAD+ demand'
        dm_nad.lower_bound = 0
        dm_nad.upper_bound = 1000
        dm_nad.add_metabolites({nad_c: -1})
        model.add_reactions([dm_nad])
        
        model.objective = dm_nad.id
        solution = model.optimize()
        
        if solution.status == 'optimal' and solution.objective_value > 1e-6:
            print(f"  [OK] NAD+ 생산 가능: {solution.objective_value:.6f}")
        else:
            print(f"  [FAIL] NAD+ 생산 불가능 ({solution.status})")
            print("    → 이것이 NADH 생산 차단의 원인일 수 있음")
        
        model.remove_reactions([dm_nad])
        
    except KeyError:
        print("  [ERROR] nad_c metabolite 없음")

def test_final_bootstrap_strategy(model, biomass_rxn):
    """최종 부트스트랩 전략 테스트"""
    print("\n" + "="*70)
    print("최종 부트스트랩 전략 테스트")
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
    
    # 최적 부트스트랩 조합 (ATP + CoA)
    print("\n[전략 1] ATP + CoA 부트스트랩 테스트")
    print("-" * 70)
    
    bootstrap_strategies = [
        {
            'name': 'ATP만',
            'components': {'atp_c': -0.1},
            'reason': '가장 간단한 전략'
        },
        {
            'name': 'ATP + CoA',
            'components': {'atp_c': -0.05, 'coa_c': -0.01},
            'reason': 'ATP 생산 + Acetate 활성화'
        },
        {
            'name': 'ATP + CoA + NAD+',
            'components': {'atp_c': -0.05, 'coa_c': -0.01, 'nad_c': -0.01},
            'reason': 'ATP 생산 + Acetate 활성화 + NADH 생산'
        },
        {
            'name': 'ATP + PEP',
            'components': {'atp_c': -0.05, 'pep_c': -0.05},
            'reason': 'ATP 생산 + Gluconeogenesis 시작'
        }
    ]
    
    results = []
    
    for strategy in bootstrap_strategies:
        print(f"\n[{strategy['name']}]")
        print(f"  이유: {strategy['reason']}")
        print(f"  구성: {strategy['components']}")
        
        # 부트스트랩 demand 추가
        demand_rxns = []
        for met_id, supply_rate in strategy['components'].items():
            try:
                met = model.metabolites.get_by_id(met_id)
                dm_rxn = cobra.Reaction(f'DM_{met_id}_bs')
                dm_rxn.name = f'{met_id} bootstrap'
                dm_rxn.lower_bound = supply_rate
                dm_rxn.upper_bound = 1000
                dm_rxn.add_metabolites({met: -1})
                model.add_reactions([dm_rxn])
                demand_rxns.append(dm_rxn.id)
            except KeyError:
                print(f"  [SKIP] {met_id} metabolite 없음")
        
        # Biomass 생산 테스트
        model.objective = biomass_rxn.id
        solution = model.optimize()
        
        if solution.status == 'optimal':
            biomass_flux = solution.objective_value
            
            result = {
                'Strategy': strategy['name'],
                'Components': str(strategy['components']),
                'Status': solution.status,
                'Biomass_Flux': biomass_flux,
                'Can_Grow': biomass_flux > 1e-6
            }
            
            if biomass_flux > 1e-6:
                print(f"  [SUCCESS] 성장 가능: {biomass_flux:.6f} 1/h")
                
                # 주요 경로 플럭스 확인
                key_rxns = ['EX_ac_e', 'ACt', 'SUCOAACTr', 'ACS', 'CS', 'ICDHx', 
                           'AKGDH', 'SUCD', 'ATPS4rpp', 'PPS', 'PYK']
                
                print("\n  주요 경로 플럭스:")
                active_pathways = []
                for rxn_id in key_rxns:
                    try:
                        flux = solution.fluxes.get(rxn_id, 0)
                        if abs(flux) > 1e-8:
                            print(f"    {rxn_id}: {flux:.6f}")
                            active_pathways.append(rxn_id)
                    except KeyError:
                        pass
                
                result['Active_Pathways'] = len(active_pathways)
                
                # Exchange 플럭스
                acetate_uptake = abs(solution.fluxes.get('EX_ac_e', 0))
                if acetate_uptake > 0:
                    yield_biomass = biomass_flux / acetate_uptake
                    print(f"\n  Acetate 섭취: {acetate_uptake:.6f} mmol/gDCW/h")
                    print(f"  Biomass yield: {yield_biomass:.6f} gDW/mmol acetate")
                    result['Acetate_Uptake'] = acetate_uptake
                    result['Biomass_Yield'] = yield_biomass
            else:
                print(f"  [FAIL] 성장 불가능 (Biomass flux: {biomass_flux:.6f})")
                result['Active_Pathways'] = 0
                result['Acetate_Uptake'] = 0
                result['Biomass_Yield'] = 0
            
            results.append(result)
        else:
            print(f"  [FAIL] 최적화 실패: {solution.status}")
            results.append({
                'Strategy': strategy['name'],
                'Components': str(strategy['components']),
                'Status': solution.status,
                'Biomass_Flux': 0,
                'Can_Grow': False,
                'Active_Pathways': 0
            })
        
        # 부트스트랩 제거
        model.remove_reactions(demand_rxns)
    
    return results

def create_bootstrap_reaction_addition_script(model, biomass_rxn):
    """부트스트랩 반응 추가 스크립트 생성"""
    print("\n" + "="*70)
    print("부트스트랩 반응 추가 스크립트 생성")
    print("="*70)
    
    # 최적 전략: ATP + CoA
    optimal_strategy = {
        'atp_c': -0.05,
        'coa_c': -0.01
    }
    
    script_content = f"""#!/usr/bin/env python
\"\"\"
부트스트랩 반응 추가
Acetate 기반 성장을 위한 최소 부트스트랩
\"\"\"

import cobra

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    return model

def add_bootstrap_reactions(model):
    \"\"\"부트스트랩 반응 추가\"\"\"
    print("부트스트랩 반응 추가 중...")
    
    bootstrap_components = {optimal_strategy}
    
    added_reactions = []
    
    for met_id, supply_rate in bootstrap_components.items():
        try:
            met = model.metabolites.get_by_id(met_id)
            
            # 기존 demand 반응 확인
            dm_id = f'DM_{{met_id}}'
            try:
                dm_rxn = model.reactions.get_by_id(dm_id)
                print(f"  [UPDATE] {{dm_id}}: lower_bound를 {{supply_rate}}로 설정")
                dm_rxn.lower_bound = supply_rate
            except KeyError:
                # 새 반응 추가
                dm_rxn = cobra.Reaction(dm_id)
                dm_rxn.name = f'{{met_id}} bootstrap demand'
                dm_rxn.lower_bound = supply_rate
                dm_rxn.upper_bound = 1000
                dm_rxn.add_metabolites({{met: -1}})
                model.add_reactions([dm_rxn])
                print(f"  [ADD] {{dm_id}}: {{supply_rate}} mmol/gDCW/h")
                added_reactions.append(dm_id)
        except KeyError:
            print(f"  [SKIP] {{met_id}} metabolite 없음")
    
    print(f"\\n총 {{len(added_reactions)}}개 반응 추가됨")
    return added_reactions

def main():
    model = load_model("BaseModel.xml")
    
    print("="*70)
    print("부트스트랩 반응 추가")
    print("="*70)
    
    added = add_bootstrap_reactions(model)
    
    # 모델 저장
    output_path = "BaseModel_with_bootstrap.xml"
    cobra.io.write_sbml_model(model, output_path)
    print(f"\\n[OK] 모델 저장: {{output_path}}")
    
    print("="*70)

if __name__ == "__main__":
    main()
"""
    
    with open('add_bootstrap_reactions.py', 'w', encoding='utf-8') as f:
        f.write(script_content)
    
    print(f"\n[OK] 부트스트랩 추가 스크립트 생성: add_bootstrap_reactions.py")
    print(f"\n최적 부트스트랩 전략:")
    print(f"  - ATP: {optimal_strategy['atp_c']} mmol/gDCW/h")
    print(f"  - CoA: {optimal_strategy['coa_c']} mmol/gDCW/h")

def main():
    print("="*70)
    print("최종 부트스트랩 전략 수립 및 검증")
    print("="*70)
    
    model = load_model("BaseModel.xml")
    
    # Biomass 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass: {biomass_rxn.id}")
    
    # 1. NADH 생산 경로 문제 분석
    analyze_nadh_production_issue(model)
    
    # 2. 최종 부트스트랩 전략 테스트
    bootstrap_results = test_final_bootstrap_strategy(model, biomass_rxn)
    
    # 3. 부트스트랩 반응 추가 스크립트 생성
    create_bootstrap_reaction_addition_script(model, biomass_rxn)
    
    # 결과 저장
    if bootstrap_results:
        df_results = pd.DataFrame(bootstrap_results)
        df_results.to_csv('final_bootstrap_strategy_results.csv', index=False)
        print(f"\n[OK] 최종 부트스트랩 전략 결과 저장: final_bootstrap_strategy_results.csv")
    
    # 최종 요약
    print("\n" + "="*70)
    print("최종 요약")
    print("="*70)
    
    working_strategies = [r for r in bootstrap_results if r['Can_Grow']]
    
    if working_strategies:
        print(f"\n[SUCCESS] {len(working_strategies)}개 부트스트랩 전략이 작동:")
        for r in working_strategies:
            print(f"\n  [{r['Strategy']}]")
            print(f"    구성: {r['Components']}")
            print(f"    Biomass flux: {r['Biomass_Flux']:.6f} 1/h")
            if 'Biomass_Yield' in r and r['Biomass_Yield'] > 0:
                print(f"    Biomass yield: {r['Biomass_Yield']:.6f} gDW/mmol acetate")
        
        # 최적 전략 추천
        best = max(working_strategies, key=lambda x: x.get('Biomass_Flux', 0))
        print(f"\n[추천] 최적 전략: {best['Strategy']}")
        print(f"  구성: {best['Components']}")
        print(f"  Biomass flux: {best['Biomass_Flux']:.6f} 1/h")
    else:
        print("\n[FAIL] 작동하는 부트스트랩 전략 없음")
        print("  → 추가 분석 필요")
    
    print("\n" + "="*70)

if __name__ == "__main__":
    main()
