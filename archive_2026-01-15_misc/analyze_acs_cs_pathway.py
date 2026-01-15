#!/usr/bin/env python
"""
ACS와 CS 경로 연결 문제 분석
- ACS만 강제 활성화: CS 작동 안 함
- CS도 강제 활성화: 둘 다 작동
- 왜 ACS만으로는 CS가 작동하지 않는가?
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

def analyze_acs_only(model):
    """ACS만 강제 활성화 시 분석"""
    print("="*80)
    print("ACS만 강제 활성화 시 분석")
    print("="*80)
    
    model = setup_acetate_medium(model)
    
    # ATPM=0 설정
    atpm_rxn = model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    # ACS만 강제 활성화
    acs = model.reactions.get_by_id('ACS')
    acs.lower_bound = 0.1
    
    model.objective = 'Growth'
    solution = model.optimize()
    
    print(f"\n[FBA 결과]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    
    # 주요 반응 플럭스
    acs_flux = solution.fluxes.get('ACS', 0.0)
    cs_flux = solution.fluxes.get('CS', 0.0)
    adk1_flux = solution.fluxes.get('ADK1', 0.0)
    
    print(f"\n[주요 반응 플럭스]")
    print(f"  ACS: {acs_flux:.6f}")
    print(f"  CS: {cs_flux:.6f}")
    print(f"  ADK1: {adk1_flux:.6f}")
    
    # OAA 관련 플럭스
    print(f"\n[OAA 관련 반응 플럭스]")
    oaa_reactions = ['PC', 'PEPCK', 'PEPCK_ATP', 'PPC', 'MDH']
    oaa_found = False
    for rxn_id in oaa_reactions:
        try:
            flux = solution.fluxes.get(rxn_id, 0.0)
            if abs(flux) > 1e-6:
                print(f"  {rxn_id}: {flux:.6f}")
                oaa_found = True
        except:
            pass
    
    if not oaa_found:
        print(f"  OAA 생성 반응 플럭스 없음!")
        print(f"  -> OAA가 없어서 CS가 작동하지 않음!")
    
    # Acetyl-CoA 플럭스 확인
    try:
        accoa_c = model.metabolites.get_by_id('accoa_c')
        accoa_producing = []
        accoa_consuming = []
        
        for rxn in accoa_c.reactions:
            flux = solution.fluxes.get(rxn.id, 0.0)
            if abs(flux) > 1e-6:
                coeff = rxn.metabolites.get(accoa_c, 0)
                net_flux = flux * coeff
                if net_flux > 0:
                    accoa_producing.append((rxn.id, net_flux))
                elif net_flux < 0:
                    accoa_consuming.append((rxn.id, abs(net_flux)))
        
        print(f"\n[Acetyl-CoA 생성/소모]")
        print(f"  생성 반응:")
        for rxn_id, net_flux in sorted(accoa_producing, key=lambda x: x[1], reverse=True)[:5]:
            print(f"    {rxn_id}: {net_flux:.6f}")
        print(f"  소모 반응:")
        for rxn_id, net_flux in sorted(accoa_consuming, key=lambda x: x[1], reverse=True)[:5]:
            print(f"    {rxn_id}: {net_flux:.6f}")
    except KeyError:
        pass

def analyze_acs_and_cs_together(model):
    """ACS와 CS 함께 강제 활성화 시 분석"""
    print("\n" + "="*80)
    print("ACS와 CS 함께 강제 활성화 시 분석")
    print("="*80)
    
    model = setup_acetate_medium(model)
    
    # ATPM=0 설정
    atpm_rxn = model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    # ACS와 CS 함께 강제 활성화
    acs = model.reactions.get_by_id('ACS')
    acs.lower_bound = 0.1
    
    cs = model.reactions.get_by_id('CS')
    cs.lower_bound = 0.1
    
    model.objective = 'Growth'
    solution = model.optimize()
    
    print(f"\n[FBA 결과]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    
    # 주요 반응 플럭스
    acs_flux = solution.fluxes.get('ACS', 0.0)
    cs_flux = solution.fluxes.get('CS', 0.0)
    adk1_flux = solution.fluxes.get('ADK1', 0.0)
    
    print(f"\n[주요 반응 플럭스]")
    print(f"  ACS: {acs_flux:.6f}")
    print(f"  CS: {cs_flux:.6f}")
    print(f"  ADK1: {adk1_flux:.6f}")
    
    # OAA 관련 플럭스
    print(f"\n[OAA 관련 반응 플럭스]")
    oaa_reactions = ['PC', 'PEPCK', 'PEPCK_ATP', 'PPC', 'MDH']
    for rxn_id in oaa_reactions:
        try:
            flux = solution.fluxes.get(rxn_id, 0.0)
            if abs(flux) > 1e-6:
                print(f"  {rxn_id}: {flux:.6f}")
        except:
            pass

def check_oaa_availability(model):
    """OAA 가용성 확인"""
    print("\n" + "="*80)
    print("OAA 가용성 확인")
    print("="*80)
    
    model = setup_acetate_medium(model)
    
    # ATPM=0 설정
    atpm_rxn = model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    # OAA demand 반응 추가 (OAA를 만들 수 있는지 테스트)
    try:
        oaa_c = model.metabolites.get_by_id('oaa_c')
        
        # OAA demand 반응 생성
        oaa_demand = cobra.Reaction('DM_oaa_c')
        oaa_demand.name = 'OAA demand'
        oaa_demand.lower_bound = 0
        oaa_demand.upper_bound = 1000
        oaa_demand.add_metabolites({oaa_c: -1.0})
        model.add_reactions([oaa_demand])
        
        model.objective = 'DM_oaa_c'
        solution = model.optimize()
        
        print(f"\n[OAA 생성 가능 여부]")
        print(f"  상태: {solution.status}")
        print(f"  OAA 최대 생산량: {solution.objective_value:.6f}")
        
        if solution.objective_value > 1e-6:
            print(f"  -> OAA를 생성할 수 있음!")
            
            # OAA 생성 반응 확인
            oaa_producing = []
            for rxn in oaa_c.reactions:
                if 'DM_' in rxn.id:
                    continue
                flux = solution.fluxes.get(rxn.id, 0.0)
                if abs(flux) > 1e-6:
                    coeff = rxn.metabolites.get(oaa_c, 0)
                    if coeff > 0:
                        oaa_producing.append((rxn.id, flux * coeff))
            
            print(f"\n[OAA 생성 반응]")
            for rxn_id, net_flux in sorted(oaa_producing, key=lambda x: x[1], reverse=True)[:5]:
                print(f"  {rxn_id}: {net_flux:.6f}")
        else:
            print(f"  -> OAA를 생성할 수 없음!")
            print(f"  -> 이것이 CS가 작동하지 않는 이유!")
        
        # Demand 반응 제거
        model.remove_reactions([oaa_demand])
        
    except KeyError:
        print(f"  oaa_c 메타볼라이트 없음")

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    
    model = load_model(str(model_path))
    
    print("="*80)
    print("ACS와 CS 경로 연결 문제 분석")
    print("="*80)
    
    print("\n[사용자 질문]")
    print("  왜 Acetyl-CoA가 Citrate로 안 가는가?")
    print("  ADK1은 발현이 되어야 한다")
    
    # ACS만 강제 활성화 시 분석
    analyze_acs_only(model)
    
    # ACS와 CS 함께 강제 활성화 시 분석
    analyze_acs_and_cs_together(model)
    
    # OAA 가용성 확인
    check_oaa_availability(model)
    
    print("\n" + "="*80)
    print("결론")
    print("="*80)
    
    print("\n[핵심 발견]")
    print("  1. ACS만 강제 활성화: CS 작동 안 함")
    print("  2. ACS + CS 함께 강제 활성화: 둘 다 작동!")
    print("  3. ADK1도 CS와 함께 작동함")
    print("\n[원인 추정]")
    print("  -> OAA가 없어서 CS가 작동하지 않음")
    print("  -> CS를 강제 활성화하면 OAA 생성 경로가 활성화됨")
    print("\n[해결]")
    print("  OAA 생성 경로 확인 및 확보 필요")

if __name__ == "__main__":
    main()
