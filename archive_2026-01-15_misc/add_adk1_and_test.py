#!/usr/bin/env python
"""
ADK1 (Adenylate kinase) 반응 추가 후 테스트
ADK1: AMP + ATP <=> 2ADP
"""

import cobra
from pathlib import Path

def load_model(model_path):
    model = cobra.io.read_sbml_model(model_path)
    return model

def add_adk1(model):
    """ADK1 반응 추가: AMP + ATP <=> 2ADP"""
    if 'ADK1' in [r.id for r in model.reactions]:
        print("[SKIP] ADK1 반응이 이미 존재합니다")
        return False
    
    # 메타볼라이트 확인
    amp_c = model.metabolites.get_by_id('amp_c')
    atp_c = model.metabolites.get_by_id('atp_c')
    adp_c = model.metabolites.get_by_id('adp_c')
    
    # 반응 생성
    adk1 = cobra.Reaction('ADK1')
    adk1.name = 'Adenylate kinase'
    adk1.lower_bound = -1000.0
    adk1.upper_bound = 1000.0
    adk1.add_metabolites({
        amp_c: -1.0,
        atp_c: -1.0,
        adp_c: 2.0
    })
    
    model.add_reactions([adk1])
    print(f"[ADDED] ADK1: {adk1.reaction}")
    return True

def setup_acetate_medium(model):
    """Acetate 미디어 설정"""
    # 모든 exchange 차단
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    # Acetate 허용
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
    
    return model

def test_with_adk1(model):
    """ADK1 추가 후 테스트"""
    print("="*80)
    print("ADK1 추가 후 ATPM=0 테스트")
    print("="*80)
    
    # ADK1 추가
    added = add_adk1(model)
    
    # 미디어 설정
    model = setup_acetate_medium(model)
    
    # ATPM=0 설정
    atpm_rxn = model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    model.objective = 'Growth'
    
    # FBA 수행
    solution = model.optimize()
    print(f"\n[FBA 결과]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    
    if solution.status == 'optimal':
        # 주요 반응 플럭스
        acs_flux = solution.fluxes.get('ACS', 0.0)
        adk1_flux = solution.fluxes.get('ADK1', 0.0)
        icl_flux = solution.fluxes.get('ICL', 0.0)
        mals_flux = solution.fluxes.get('MALS', 0.0)
        icdhx_flux = solution.fluxes.get('ICDHx', 0.0)
        cs_flux = solution.fluxes.get('CS', 0.0)
        
        print(f"\n[주요 반응 플럭스]")
        print(f"  ACS: {acs_flux:.6f}")
        print(f"  ADK1: {adk1_flux:.6f}")
        print(f"  ICL: {icl_flux:.6f}")
        print(f"  MALS: {mals_flux:.6f}")
        print(f"  ICDHx: {icdhx_flux:.6f}")
        print(f"  CS: {cs_flux:.6f}")
        
        if abs(acs_flux) > 1e-6:
            print("\n[OK] ACS 반응이 작동합니다!")
            if abs(adk1_flux) > 1e-6:
                print(f"  -> ADK1도 작동 중 (플럭스: {adk1_flux:.6f})")
                print(f"     AMP + ATP -> 2ADP로 AMP를 재활용 중")
        else:
            print("\n[문제] ADK1을 추가해도 ACS가 작동하지 않습니다.")
            print("  -> 다른 원인이 있을 수 있습니다.")
    
    return solution

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    
    model = load_model(str(model_path))
    
    print("="*80)
    print("ADK1 추가 및 테스트")
    print("="*80)
    
    print("\n[레퍼런스 모델 정보]")
    print("  ADK1 (Adenylate kinase): AMP + ATP <=> 2ADP")
    print("  레퍼런스 모델 FBA에서 ADK1 플럭스: 0.2096 (ATPM=0)")
    print("  -> AMP를 ADP로 전환하여 ATP 재생성 경로 제공")
    
    # ADK1 추가 후 테스트
    solution = test_with_adk1(model)
    
    print("\n" + "="*80)
    print("결론")
    print("="*80)
    
    if solution and solution.objective_value > 1e-6:
        print("\n[성공] ADK1 추가로 문제 해결!")
        print("  -> ACS 반응이 작동하여 AMP가 생성됨")
        print("  -> ADK1이 AMP를 ADP로 전환하여 ATP 재생성")
    elif solution:
        acs_flux = solution.fluxes.get('ACS', 0.0)
        if abs(acs_flux) > 1e-6:
            print("\n[부분 성공] ACS는 작동하지만 성장률이 0")
            print("  -> 추가 반응이 필요할 수 있습니다")
        else:
            print("\n[실패] ADK1 추가로도 해결되지 않음")
            print("  -> 다른 원인 조사 필요")

if __name__ == "__main__":
    main()
