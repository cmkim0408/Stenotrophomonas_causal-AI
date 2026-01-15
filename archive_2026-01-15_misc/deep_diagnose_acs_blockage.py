#!/usr/bin/env python
"""
ACS가 작동하지 않는 깊은 원인 분석
- 부트스트랩 제공했지만 여전히 작동 안 함
- Not blocked인데도 플럭스 0
- 왜 최적화가 ACS를 선택하지 않는가?
"""

import cobra
from pathlib import Path

def load_model(model_path):
    model = cobra.io.read_sbml_model(model_path)
    return model

def setup_acetate_medium(model):
    """Acetate 미디어 설정"""
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    ex_ac = model.reactions.get_by_id('EX_ac_e')
    ex_ac.upper_bound = 1000
    ex_ac.lower_bound = -1000
    
    essential = {
        'EX_nh4_e': (-1000, 1000),
        'EX_h2o_e': (-1000, 1000),
        'EX_h_e': (-1000, 1000),
        'EX_pi_e': (-1000, 1000),
        'EX_so4_e': (-1000, 1000),
        'EX_k_e': (-1000, 1000),
        'EX_na1_e': (-1000, 1000),
        'EX_mg2_e': (-1000, 1000),
        'EX_ca2_e': (-1000, 1000),
        'EX_fe2_e': (-1000, 1000),
        'EX_mn2_e': (-1000, 1000),
        'EX_zn2_e': (-1000, 1000),
        'EX_co2_e': (-1000, 1000),
        'EX_o2_e': (-1000, 1000),
    }
    
    for ex_id, (lb, ub) in essential.items():
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            ex_rxn.upper_bound = ub
            ex_rxn.lower_bound = lb
        except KeyError:
            pass
    
    return model

def add_bootstrap_exchanges(model, bootstrap_amount=0.001):
    """부트스트랩 Exchange 추가"""
    bootstrap_metabolites = {
        'oaa_c': 'OAA',
        'mal__L_c': 'Malate',
        'asp__L_c': 'Aspartate',
    }
    
    for met_id, met_name in bootstrap_metabolites.items():
        try:
            met = model.metabolites.get_by_id(met_id)
            ex_id = f'EX_{met_id}'
            
            if ex_id not in [r.id for r in model.exchanges]:
                ex_rxn = cobra.Reaction(ex_id)
                ex_rxn.name = f'{met_name} exchange'
                ex_rxn.lower_bound = -bootstrap_amount
                ex_rxn.upper_bound = 1000
                
                met_e_id = met_id.replace('_c', '_e')
                try:
                    met_e = model.metabolites.get_by_id(met_e_id)
                except KeyError:
                    met_e = cobra.Metabolite(met_e_id, name=met.name, compartment='e')
                    model.add_metabolites([met_e])
                
                ex_rxn.add_metabolites({met_e: -1.0})
                model.add_reactions([ex_rxn])
            else:
                ex_rxn = model.reactions.get_by_id(ex_id)
                ex_rxn.lower_bound = -bootstrap_amount
                ex_rxn.upper_bound = 1000
        except KeyError:
            pass

def check_coa_availability(model):
    """CoA 가용성 확인"""
    print("="*80)
    print("CoA 가용성 확인")
    print("="*80)
    
    model = setup_acetate_medium(model)
    add_bootstrap_exchanges(model)
    
    # ATPM=0 설정
    atpm_rxn = model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    try:
        coa_c = model.metabolites.get_by_id('coa_c')
        
        # CoA demand 반응 생성
        coa_demand = cobra.Reaction('DM_coa_c')
        coa_demand.name = 'CoA demand'
        coa_demand.lower_bound = 0
        coa_demand.upper_bound = 1000
        coa_demand.add_metabolites({coa_c: -1.0})
        model.add_reactions([coa_demand])
        
        model.objective = 'DM_coa_c'
        solution = model.optimize()
        
        print(f"\n[CoA 생성 가능 여부]")
        print(f"  상태: {solution.status}")
        print(f"  CoA 최대 생산량: {solution.objective_value:.6f}")
        
        if solution.objective_value > 1e-6:
            print(f"  -> CoA를 생성할 수 있음!")
            
            # CoA 생성 반응 확인
            coa_producing = []
            for rxn in coa_c.reactions:
                if 'DM_' in rxn.id:
                    continue
                flux = solution.fluxes.get(rxn.id, 0.0)
                if abs(flux) > 1e-6:
                    coeff = rxn.metabolites.get(coa_c, 0)
                    if coeff > 0:
                        coa_producing.append((rxn.id, flux * coeff))
            
            print(f"\n[CoA 생성 반응]")
            for rxn_id, net_flux in sorted(coa_producing, key=lambda x: x[1], reverse=True)[:5]:
                print(f"  {rxn_id}: {net_flux:.6f}")
        else:
            print(f"  -> CoA를 생성할 수 없음!")
            print(f"  -> 이것이 ACS가 작동하지 않는 이유!")
        
        model.remove_reactions([coa_demand])
        
    except KeyError:
        print(f"  coa_c 메타볼라이트 없음")

def check_atp_availability(model):
    """ATP 가용성 확인"""
    print("\n" + "="*80)
    print("ATP 가용성 확인")
    print("="*80)
    
    model = setup_acetate_medium(model)
    add_bootstrap_exchanges(model)
    
    # ATPM=0 설정
    atpm_rxn = model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    try:
        atp_c = model.metabolites.get_by_id('atp_c')
        
        # ATP demand 반응 생성
        atp_demand = cobra.Reaction('DM_atp_c')
        atp_demand.name = 'ATP demand'
        atp_demand.lower_bound = 0
        atp_demand.upper_bound = 1000
        atp_demand.add_metabolites({atp_c: -1.0})
        model.add_reactions([atp_demand])
        
        model.objective = 'DM_atp_c'
        solution = model.optimize()
        
        print(f"\n[ATP 생성 가능 여부]")
        print(f"  상태: {solution.status}")
        print(f"  ATP 최대 생산량: {solution.objective_value:.6f}")
        
        if solution.objective_value > 1e-6:
            print(f"  -> ATP를 생성할 수 있음!")
        else:
            print(f"  -> ATP를 생성할 수 없음!")
            print(f"  -> 이것이 ACS가 작동하지 않는 이유!")
        
        model.remove_reactions([atp_demand])
        
    except KeyError:
        print(f"  atp_c 메타볼라이트 없음")

def test_acs_with_coa_bootstrap(model):
    """CoA 부트스트랩 제공 후 ACS 테스트"""
    print("\n" + "="*80)
    print("CoA 부트스트랩 제공 후 ACS 테스트")
    print("="*80)
    
    model = setup_acetate_medium(model)
    add_bootstrap_exchanges(model)
    
    # CoA 부트스트랩 추가
    try:
        coa_c = model.metabolites.get_by_id('coa_c')
        ex_coa_id = 'EX_coa_c'
        
        if ex_coa_id not in [r.id for r in model.exchanges]:
            ex_coa = cobra.Reaction(ex_coa_id)
            ex_coa.name = 'CoA exchange'
            ex_coa.lower_bound = -0.001
            ex_coa.upper_bound = 1000
            
            coa_e_id = 'coa_e'
            try:
                coa_e = model.metabolites.get_by_id(coa_e_id)
            except KeyError:
                coa_e = cobra.Metabolite(coa_e_id, name='CoA', compartment='e')
                model.add_metabolites([coa_e])
            
            ex_coa.add_metabolites({coa_e: -1.0})
            model.add_reactions([ex_coa])
            print(f"  CoA 부트스트랩 추가: EX_coa_c")
        else:
            ex_coa = model.reactions.get_by_id(ex_coa_id)
            ex_coa.lower_bound = -0.001
            ex_coa.upper_bound = 1000
            print(f"  CoA 부트스트랩 설정: EX_coa_c")
    except KeyError:
        print(f"  coa_c 메타볼라이트 없음")
    
    # ATPM=0 설정
    atpm_rxn = model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    model.objective = 'Growth'
    solution = model.optimize()
    
    print(f"\n[FBA 결과]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    
    acs_flux = solution.fluxes.get('ACS', 0.0)
    cs_flux = solution.fluxes.get('CS', 0.0)
    adk1_flux = solution.fluxes.get('ADK1', 0.0)
    
    print(f"\n[주요 반응 플럭스]")
    print(f"  ACS: {acs_flux:.6f}")
    print(f"  CS: {cs_flux:.6f}")
    print(f"  ADK1: {adk1_flux:.6f}")
    
    if abs(acs_flux) > 1e-6:
        print(f"\n[OK] CoA 부트스트랩으로 ACS 작동!")
    else:
        print(f"\n[문제] CoA 부트스트랩으로도 ACS 작동 안 함")

def test_acs_with_atp_bootstrap(model):
    """ATP 부트스트랩 제공 후 ACS 테스트"""
    print("\n" + "="*80)
    print("ATP 부트스트랩 제공 후 ACS 테스트")
    print("="*80)
    
    model = setup_acetate_medium(model)
    add_bootstrap_exchanges(model)
    
    # ATP 부트스트랩 추가
    try:
        atp_c = model.metabolites.get_by_id('atp_c')
        ex_atp_id = 'EX_atp_c'
        
        if ex_atp_id not in [r.id for r in model.exchanges]:
            ex_atp = cobra.Reaction(ex_atp_id)
            ex_atp.name = 'ATP exchange'
            ex_atp.lower_bound = -0.001
            ex_atp.upper_bound = 1000
            
            atp_e_id = 'atp_e'
            try:
                atp_e = model.metabolites.get_by_id(atp_e_id)
            except KeyError:
                atp_e = cobra.Metabolite(atp_e_id, name='ATP', compartment='e')
                model.add_metabolites([atp_e])
            
            ex_atp.add_metabolites({atp_e: -1.0})
            model.add_reactions([ex_atp])
            print(f"  ATP 부트스트랩 추가: EX_atp_c")
        else:
            ex_atp = model.reactions.get_by_id(ex_atp_id)
            ex_atp.lower_bound = -0.001
            ex_atp.upper_bound = 1000
            print(f"  ATP 부트스트랩 설정: EX_atp_c")
    except KeyError:
        print(f"  atp_c 메타볼라이트 없음")
    
    # ATPM=0 설정
    atpm_rxn = model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    model.objective = 'Growth'
    solution = model.optimize()
    
    print(f"\n[FBA 결과]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    
    acs_flux = solution.fluxes.get('ACS', 0.0)
    cs_flux = solution.fluxes.get('CS', 0.0)
    adk1_flux = solution.fluxes.get('ADK1', 0.0)
    
    print(f"\n[주요 반응 플럭스]")
    print(f"  ACS: {acs_flux:.6f}")
    print(f"  CS: {cs_flux:.6f}")
    print(f"  ADK1: {adk1_flux:.6f}")
    
    if abs(acs_flux) > 1e-6:
        print(f"\n[OK] ATP 부트스트랩으로 ACS 작동!")
    else:
        print(f"\n[문제] ATP 부트스트랩으로도 ACS 작동 안 함")

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    
    model = load_model(str(model_path))
    
    print("="*80)
    print("ACS 작동 안 함 깊은 원인 분석")
    print("="*80)
    
    # CoA 가용성 확인
    check_coa_availability(model)
    
    # ATP 가용성 확인
    check_atp_availability(model)
    
    # CoA 부트스트랩 테스트
    test_acs_with_coa_bootstrap(model)
    
    # ATP 부트스트랩 테스트
    test_acs_with_atp_bootstrap(model)
    
    print("\n" + "="*80)
    print("결론")
    print("="*80)
    
    print("\n[확인 필요]")
    print("  1. CoA 초기화 문제인지")
    print("  2. ATP 초기화 문제인지")
    print("  3. 다른 제약 조건인지")

if __name__ == "__main__":
    main()
