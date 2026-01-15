#!/usr/bin/env python
"""
ACS가 자동으로 작동하지 않는 이유 분석
- 강제 활성화하면 작동하지만
- 자동으로는 작동하지 않음
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

def test_without_force(model):
    """강제 활성화 없이 테스트"""
    print("="*80)
    print("강제 활성화 없이 테스트 (자동 최적화)")
    print("="*80)
    
    model = setup_acetate_medium(model)
    
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
    icl_flux = solution.fluxes.get('ICL', 0.0)
    mals_flux = solution.fluxes.get('MALS', 0.0)
    
    print(f"\n[주요 반응 플럭스 (자동)]")
    print(f"  ACS: {acs_flux:.6f}")
    print(f"  CS: {cs_flux:.6f}")
    print(f"  ADK1: {adk1_flux:.6f}")
    print(f"  ICL: {icl_flux:.6f}")
    print(f"  MALS: {mals_flux:.6f}")
    
    if abs(acs_flux) < 1e-6:
        print(f"\n[문제] ACS가 자동으로 작동하지 않음")
        return False
    else:
        print(f"\n[OK] ACS가 자동으로 작동함")
        return True

def test_with_force(model):
    """강제 활성화 시 테스트"""
    print("\n" + "="*80)
    print("강제 활성화 시 테스트")
    print("="*80)
    
    model = setup_acetate_medium(model)
    
    # ATPM=0 설정
    atpm_rxn = model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    # ACS 강제 활성화
    acs = model.reactions.get_by_id('ACS')
    acs.lower_bound = 0.1
    
    model.objective = 'Growth'
    solution = model.optimize()
    
    print(f"\n[FBA 결과 (ACS 강제 활성화)]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    
    acs_flux = solution.fluxes.get('ACS', 0.0)
    cs_flux = solution.fluxes.get('CS', 0.0)
    adk1_flux = solution.fluxes.get('ADK1', 0.0)
    
    print(f"\n[주요 반응 플럭스 (강제)]")
    print(f"  ACS: {acs_flux:.6f}")
    print(f"  CS: {cs_flux:.6f}")
    print(f"  ADK1: {adk1_flux:.6f}")
    
    if abs(cs_flux) > 1e-6 and abs(adk1_flux) > 1e-6:
        print(f"\n[OK] ACS 강제 활성화 시 CS와 ADK1도 작동!")
        return True
    else:
        print(f"\n[문제] ACS 강제 활성화해도 CS나 ADK1이 작동 안 함")
        return False

def analyze_blocked_reactions(model):
    """Blocked reaction 분석"""
    print("\n" + "="*80)
    print("Blocked Reaction 분석")
    print("="*80)
    
    model = setup_acetate_medium(model)
    
    # ATPM=0 설정
    atpm_rxn = model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    from cobra.flux_analysis import find_blocked_reactions
    blocked = find_blocked_reactions(model)
    
    key_reactions = ['ACS', 'CS', 'ADK1', 'ICL', 'MALS', 'ICDHx']
    
    print(f"\n[주요 반응 Blocked 여부]")
    for rxn_id in key_reactions:
        if rxn_id in blocked:
            print(f"  {rxn_id}: Blocked!")
        else:
            print(f"  {rxn_id}: Not blocked")
    
    # ACS가 blocked인 이유 확인
    if 'ACS' in blocked:
        print(f"\n[ACS가 Blocked인 이유 확인 필요]")
    else:
        print(f"\n[ACS는 Not blocked - 작동 가능한 상태]")
        print(f"  -> 하지만 최적화에서 선택되지 않음")
        print(f"  -> 아마도 더 효율적인 경로가 없거나")
        print(f"  -> 경로가 완전하지 않아서 선택되지 않음")

def check_alternative_pathways(model):
    """대체 경로 확인"""
    print("\n" + "="*80)
    print("대체 경로 확인")
    print("="*80)
    
    model = setup_acetate_medium(model)
    
    # ATPM=0 설정
    atpm_rxn = model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    model.objective = 'Growth'
    solution = model.optimize()
    
    # Acetyl-CoA 생성 경로 확인
    print(f"\n[Acetyl-CoA 생성 경로]")
    try:
        accoa_c = model.metabolites.get_by_id('accoa_c')
        accoa_producing = []
        
        for rxn in accoa_c.reactions:
            flux = solution.fluxes.get(rxn.id, 0.0)
            if abs(flux) > 1e-6:
                coeff = rxn.metabolites.get(accoa_c, 0)
                if coeff > 0:  # 생성
                    accoa_producing.append((rxn.id, flux * coeff))
        
        if accoa_producing:
            print(f"  Acetyl-CoA 생성 반응 (플럭스 > 0):")
            for rxn_id, net_flux in sorted(accoa_producing, key=lambda x: x[1], reverse=True):
                print(f"    {rxn_id}: {net_flux:.6f}")
        else:
            print(f"  Acetyl-CoA 생성 반응 없음 (모든 플럭스 0)")
    except KeyError:
        print(f"  accoa_c 메타볼라이트 없음")
    
    # Acetate 소모 경로 확인
    print(f"\n[Acetate 소모 경로]")
    try:
        ac_c = model.metabolites.get_by_id('ac_c')
        ac_consuming = []
        
        for rxn in ac_c.reactions:
            flux = solution.fluxes.get(rxn.id, 0.0)
            if abs(flux) > 1e-6:
                coeff = rxn.metabolites.get(ac_c, 0)
                if coeff < 0:  # 소모
                    ac_consuming.append((rxn.id, abs(flux * coeff)))
        
        if ac_consuming:
            print(f"  Acetate 소모 반응 (플럭스 > 0):")
            for rxn_id, net_flux in sorted(ac_consuming, key=lambda x: x[1], reverse=True):
                print(f"    {rxn_id}: {net_flux:.6f}")
        else:
            print(f"  Acetate 소모 반응 없음 (모든 플럭스 0)")
    except KeyError:
        print(f"  ac_c 메타볼라이트 없음")

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    
    model = load_model(str(model_path))
    
    print("="*80)
    print("ACS가 자동으로 작동하지 않는 이유 분석")
    print("="*80)
    
    print("\n[사용자 질문]")
    print("  왜 Acetyl-CoA가 Citrate로 안 가는가?")
    print("  ADK1은 발현이 되어야 한다")
    
    # 강제 활성화 없이 테스트
    auto_works = test_without_force(model)
    
    # 강제 활성화 시 테스트
    force_works = test_with_force(model)
    
    # Blocked reaction 분석
    analyze_blocked_reactions(model)
    
    # 대체 경로 확인
    check_alternative_pathways(model)
    
    print("\n" + "="*80)
    print("결론")
    print("="*80)
    
    if not auto_works and force_works:
        print("\n[핵심 발견]")
        print("  ACS는 작동 가능하지만 자동으로 선택되지 않음")
        print("  강제 활성화하면 CS, ADK1 모두 작동!")
        print("\n[원인]")
        print("  FBA 최적화가 ACS를 선택하지 않음")
        print("  -> 아마도 경로가 완전하지 않아서")
        print("  -> 또는 초기화 문제 (CoA, OAA 등)")
        print("\n[해결]")
        print("  경로 완전성 확보 필요")
        print("  또는 초기화 문제 해결 필요")
    elif auto_works:
        print("\n[OK] ACS가 자동으로 작동합니다!")
    else:
        print("\n[문제] ACS가 작동하지 않습니다")
        print("  -> 다른 원인 조사 필요")

if __name__ == "__main__":
    main()
