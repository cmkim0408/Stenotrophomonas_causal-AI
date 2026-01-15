#!/usr/bin/env python
"""
ACS 반응이 작동하지 않는 실제 원인 디버깅
- ADK1은 있음
- ACS는 있는데 작동 안 함
- 실제 원인 찾기 (CoA, 제약 조건 등)
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

def check_acs_requirements(model):
    """ACS 반응에 필요한 메타볼라이트 확인"""
    print("="*80)
    print("ACS 반응 요구사항 확인")
    print("="*80)
    
    acs = model.reactions.get_by_id('ACS')
    print(f"\n[ACS 반응]")
    print(f"  반응식: {acs.reaction}")
    print(f"  bounds: [{acs.lower_bound}, {acs.upper_bound}]")
    
    # 필요한 메타볼라이트
    required_mets = ['ac_c', 'atp_c', 'coa_c']
    produced_mets = ['accoa_c', 'amp_c', 'ppi_c']
    
    print(f"\n[필요한 메타볼라이트 (소모)]")
    for met_id in required_mets:
        try:
            met = model.metabolites.get_by_id(met_id)
            coeff = acs.metabolites.get(met, 0)
            print(f"  {met_id}: 계수={coeff}, 존재={True}")
        except KeyError:
            print(f"  {met_id}: 없음!")
    
    print(f"\n[생성되는 메타볼라이트]")
    for met_id in produced_mets:
        try:
            met = model.metabolites.get_by_id(met_id)
            coeff = acs.metabolites.get(met, 0)
            print(f"  {met_id}: 계수={coeff}, 존재={True}")
        except KeyError:
            print(f"  {met_id}: 없음!")

def test_acs_directly(model):
    """ACS 반응을 직접 테스트 (CoA 공급 등)"""
    print("\n" + "="*80)
    print("ACS 반응 직접 테스트")
    print("="*80)
    
    model = setup_acetate_medium(model)
    
    # ATPM=0 설정
    atpm_rxn = model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    # CoA 초기화 테스트 (CoA를 공급해보기)
    print("\n[테스트 1] CoA 공급 테스트")
    try:
        ex_coa = model.reactions.get_by_id('EX_coa_e')
        ex_coa.upper_bound = 1000
        ex_coa.lower_bound = -0.1  # CoA 소량 공급
        print(f"  EX_coa_e bounds: [{ex_coa.lower_bound}, {ex_coa.upper_bound}]")
    except KeyError:
        print(f"  EX_coa_e 없음 (CoA는 내부에서 합성해야 함)")
    
    model.objective = 'Growth'
    solution = model.optimize()
    
    acs_flux = solution.fluxes.get('ACS', 0.0)
    print(f"  ACS 플럭스: {acs_flux:.6f}")
    
    if abs(acs_flux) > 1e-6:
        print(f"  -> CoA 공급 시 ACS 작동!")
    else:
        print(f"  -> CoA 공급해도 ACS 작동 안 함")
    
    # ATP 초기화 테스트
    print("\n[테스트 2] ATP 초기화 테스트")
    try:
        ex_atp = model.reactions.get_by_id('EX_atp_e')
        ex_atp.upper_bound = 1000
        ex_atp.lower_bound = -0.1  # ATP 소량 공급
        print(f"  EX_atp_e bounds: [{ex_atp.lower_bound}, {ex_atp.upper_bound}]")
    except KeyError:
        print(f"  EX_atp_e 없음 (ATP는 내부에서 생성해야 함)")
    
    solution2 = model.optimize()
    acs_flux2 = solution2.fluxes.get('ACS', 0.0)
    print(f"  ACS 플럭스: {acs_flux2:.6f}")
    
    if abs(acs_flux2) > 1e-6:
        print(f"  -> ATP 공급 시 ACS 작동!")
    else:
        print(f"  -> ATP 공급해도 ACS 작동 안 함")
    
    # ACS 반응을 강제로 활성화해보기
    print("\n[테스트 3] ACS 반응 강제 활성화 테스트")
    acs = model.reactions.get_by_id('ACS')
    acs.lower_bound = 0.1  # 최소 플럭스 강제
    
    model.objective = 'Growth'
    solution3 = model.optimize()
    
    print(f"  상태: {solution3.status}")
    acs_flux3 = solution3.fluxes.get('ACS', 0.0)
    print(f"  ACS 플럭스: {acs_flux3:.6f}")
    
    if solution3.status == 'optimal':
        print(f"  -> ACS 강제 활성화 시 작동 가능!")
        print(f"     성장률: {solution3.objective_value:.6f}")
    else:
        print(f"  -> ACS 강제 활성화 불가 (제약 조건 위반)")
        print(f"     상태: {solution3.status}")

def check_blocked_reactions(model):
    """Blocked reaction 확인"""
    print("\n" + "="*80)
    print("Blocked Reaction 확인")
    print("="*80)
    
    model = setup_acetate_medium(model)
    
    # ATPM=0 설정
    atpm_rxn = model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    # 주요 반응이 blocked인지 확인
    key_reactions = ['ACS', 'ADK1', 'ICL', 'MALS', 'CS', 'ICDHx']
    
    from cobra.flux_analysis import find_blocked_reactions
    blocked = find_blocked_reactions(model)
    
    print(f"\n[Blocked Reactions]")
    print(f"  총 blocked 반응: {len(blocked)}개")
    
    for rxn_id in key_reactions:
        if rxn_id in blocked:
            print(f"  {rxn_id}: Blocked!")
        else:
            print(f"  {rxn_id}: Not blocked")

def compare_with_reference_media():
    """레퍼런스 미디어와 비교"""
    print("\n" + "="*80)
    print("레퍼런스 미디어와 비교")
    print("="*80)
    
    base_path = Path(__file__).parent.parent
    media_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5" / "Acetate_YE0p5__nocmnt__normalized.tsv"
    
    if media_path.exists():
        print(f"  레퍼런스 미디어 파일: {media_path}")
        print(f"  -> 이 미디어에는 CoA 합성에 필요한 전구체가 포함되어 있을 수 있음")
        print(f"  -> 예: pnto__R_e (판토텐산), ncam_e (니코틴아마이드) 등")

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    
    model = load_model(str(model_path))
    
    print("="*80)
    print("ACS 반응 작동 안 함 원인 디버깅")
    print("="*80)
    
    # ACS 요구사항 확인
    check_acs_requirements(model)
    
    # ACS 직접 테스트
    test_acs_directly(model)
    
    # Blocked reaction 확인
    check_blocked_reactions(model)
    
    # 레퍼런스 미디어와 비교
    compare_with_reference_media()
    
    print("\n" + "="*80)
    print("결론")
    print("="*80)
    
    print("\n[사용자 지적]")
    print("  ACS + ADK1 조합으로 작동 가능해야 함 (에너지 효율이 낮아도 OK)")
    print("  -> 에너지 효율이 낮으면 TCA cycle을 통해 보상 (NADH 생성, 탄소 손실)")
    
    print("\n[확인 필요]")
    print("  1. CoA 초기화 문제?")
    print("  2. 미디어 설정 차이? (레퍼런스 미디어 사용)")
    print("  3. 다른 제약 조건?")
    print("  4. Blocked reaction 문제?")

if __name__ == "__main__":
    main()
