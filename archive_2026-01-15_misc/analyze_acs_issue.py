#!/usr/bin/env python
"""
ACS 반응이 작동하지 않는 이유 분석
- ACS: ATP → AMP + PPi
- AMP → ATP 재생성 경로 확인 (Adenylate kinase 등)
"""

import cobra
from pathlib import Path

def load_model(model_path):
    model = cobra.io.read_sbml_model(model_path)
    return model

def find_adenylate_kinase(model):
    """Adenylate kinase 반응 찾기 (AMP + ATP <=> 2ADP)"""
    adk_keywords = ['ADK', 'adenylate', 'AMP', 'ATP', 'ADP']
    
    adk_reactions = []
    for rxn in model.reactions:
        rxn_id = rxn.id.upper()
        rxn_reaction = rxn.reaction.upper()
        
        # Adenylate kinase 패턴: AMP + ATP <=> 2ADP 또는 2ADP <=> AMP + ATP
        if 'ADK' in rxn_id or ('AMP' in rxn_reaction and 'ATP' in rxn_reaction and 'ADP' in rxn_reaction):
            if '2.0 ADP' in rxn.reaction or '2 ADP' in rxn.reaction:
                adk_reactions.append({
                    'id': rxn.id,
                    'name': rxn.name if rxn.name else '',
                    'reaction': rxn.reaction,
                    'bounds': (rxn.lower_bound, rxn.upper_bound)
                })
    
    return adk_reactions

def analyze_acs_reaction(model):
    """ACS 반응 분석"""
    print("="*80)
    print("ACS 반응 분석")
    print("="*80)
    
    # ACS 반응 찾기
    acs_reactions = []
    for rxn in model.reactions:
        if 'ACS' in rxn.id and 'ACS_ADP' not in rxn.id:
            acs_reactions.append({
                'id': rxn.id,
                'name': rxn.name if rxn.name else '',
                'reaction': rxn.reaction,
                'bounds': (rxn.lower_bound, rxn.upper_bound)
            })
    
    print(f"\n[ACS 반응] {len(acs_reactions)}개 발견")
    for acs_info in acs_reactions:
        print(f"\n  {acs_info['id']}")
        if acs_info['name']:
            print(f"    이름: {acs_info['name']}")
        print(f"    반응식: {acs_info['reaction']}")
        print(f"    bounds: [{acs_info['bounds'][0]}, {acs_info['bounds'][1]}]")
    
    # Adenylate kinase 찾기
    adk_reactions = find_adenylate_kinase(model)
    print(f"\n[Adenylate kinase 반응] {len(adk_reactions)}개 발견")
    for adk_info in adk_reactions:
        print(f"\n  {adk_info['id']}")
        if adk_info['name']:
            print(f"    이름: {adk_info['name']}")
        print(f"    반응식: {adk_info['reaction']}")
        print(f"    bounds: [{adk_info['bounds'][0]}, {adk_info['bounds'][1]}]")
    
    return acs_reactions, adk_reactions

def test_acs_with_atpm0(model):
    """ATPM=0일 때 ACS 반응 테스트"""
    print("\n" + "="*80)
    print("ATPM=0일 때 ACS 반응 테스트")
    print("="*80)
    
    # 미디어 설정 (간단한 설정)
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    ex_ac = model.reactions.get_by_id('EX_ac_e')
    ex_ac.upper_bound = 1000
    ex_ac.lower_bound = -1000
    
    # 필수 무기염
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
    
    # ATPM=0 설정
    atpm_rxn = model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 0
    
    model.objective = 'Growth'
    
    # FBA 수행
    solution = model.optimize()
    print(f"\n[FBA 결과]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    
    # ACS 플럭스 확인
    acs_flux = solution.fluxes.get('ACS', 0.0)
    print(f"\n[ACS 반응 플럭스]")
    print(f"  ACS: {acs_flux:.6f}")
    
    # Adenylate kinase 플럭스 확인
    adk_reactions = find_adenylate_kinase(model)
    if adk_reactions:
        print(f"\n[Adenylate kinase 플럭스]")
        for adk_info in adk_reactions:
            adk_flux = solution.fluxes.get(adk_info['id'], 0.0)
            print(f"  {adk_info['id']}: {adk_flux:.6f}")
    
    # ATP, AMP, ADP 플럭스 확인
    print(f"\n[ATP 관련 반응 플럭스]")
    try:
        atp_c = model.metabolites.get_by_id('atp_c')
        amp_c = model.metabolites.get_by_id('amp_c')
        adp_c = model.metabolites.get_by_id('adp_c')
        
        print(f"\n  ATP 생성/소모:")
        atp_producing = []
        atp_consuming = []
        for rxn in atp_c.reactions:
            flux = solution.fluxes.get(rxn.id, 0.0)
            if abs(flux) > 1e-6:
                coeff = rxn.metabolites.get(atp_c, 0)
                if coeff > 0:
                    atp_producing.append((rxn.id, flux * coeff))
                elif coeff < 0:
                    atp_consuming.append((rxn.id, abs(flux * coeff)))
        
        for rxn_id, net_flux in sorted(atp_producing, key=lambda x: x[1], reverse=True)[:5]:
            print(f"    생성: {rxn_id}: {net_flux:.6f}")
        for rxn_id, net_flux in sorted(atp_consuming, key=lambda x: x[1], reverse=True)[:5]:
            print(f"    소모: {rxn_id}: {net_flux:.6f}")
            
    except KeyError as e:
        print(f"  메타볼라이트 없음: {e}")

def compare_acs_vs_acs_adp(model):
    """ACS vs ACS_ADP 비교"""
    print("\n" + "="*80)
    print("ACS vs ACS_ADP 비교")
    print("="*80)
    
    try:
        acs = model.reactions.get_by_id('ACS')
        print(f"\n[ACS]")
        print(f"  반응식: {acs.reaction}")
        print(f"  ATP → AMP + PPi (에너지 비용: 높음)")
    except KeyError:
        print("\n[ACS] 반응 없음")
    
    try:
        acs_adp = model.reactions.get_by_id('ACS_ADP')
        print(f"\n[ACS_ADP]")
        print(f"  반응식: {acs_adp.reaction}")
        print(f"  ATP → ADP + Pi (에너지 비용: 낮음)")
    except KeyError:
        print("\n[ACS_ADP] 반응 없음 (레퍼런스 모델에만 있음)")

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    
    model = load_model(str(model_path))
    
    print("="*80)
    print("ACS 반응 작동 분석")
    print("="*80)
    
    # ACS 반응 분석
    acs_reactions, adk_reactions = analyze_acs_reaction(model)
    
    # ACS vs ACS_ADP 비교
    compare_acs_vs_acs_adp(model)
    
    # ATPM=0에서 ACS 테스트
    test_acs_with_atpm0(model)
    
    print("\n" + "="*80)
    print("결론")
    print("="*80)
    
    if adk_reactions:
        print("\n[확인] Adenylate kinase 반응이 있습니다.")
        print("  -> AMP + ATP <=> 2ADP 경로로 AMP를 ATP로 재생성 가능")
        print("\n[문제] 그런데도 ACS가 작동하지 않는다면:")
        print("  1. 다른 제약 조건이 있을 수 있음")
        print("  2. ACS 반응 자체에 문제가 있을 수 있음")
        print("  3. CoA 초기화 문제일 수 있음")
    else:
        print("\n[문제] Adenylate kinase 반응이 없을 수 있습니다.")
        print("  -> AMP를 ATP로 재생성할 수 없음")
        print("  -> ACS 반응이 작동하지 않을 수 있음")

if __name__ == "__main__":
    main()
