#!/usr/bin/env python
"""
레퍼런스 모델을 미디어 없이 테스트
- 레퍼런스 모델의 기본 설정에서 FBA 실행
- 어떤 조건에서 작동하는지 확인
"""

import cobra
from pathlib import Path

def load_model(model_path):
    model = cobra.io.read_sbml_model(model_path)
    return model

def test_reference_without_media():
    """레퍼런스 모델을 미디어 없이 테스트"""
    base_path = Path(__file__).parent.parent
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    
    print("="*80)
    print("레퍼런스 모델 기본 설정 테스트 (미디어 없이)")
    print("="*80)
    
    # 모델 로드
    ref_model = load_model(str(ref_model_path))
    
    print(f"\n[모델 정보]")
    print(f"  반응 수: {len(ref_model.reactions)}")
    print(f"  메타볼라이트 수: {len(ref_model.metabolites)}")
    
    # 기본 ATPM 확인
    atpm_rxn = ref_model.reactions.get_by_id('ATPM')
    print(f"\n[기본 ATPM 설정]")
    print(f"  bounds: [{atpm_rxn.lower_bound}, {atpm_rxn.upper_bound}]")
    print(f"  반응식: {atpm_rxn.reaction}")
    
    # Exchange bounds 확인 (기본값)
    print(f"\n[주요 Exchange bounds (기본값)]")
    key_exchanges = ['EX_ac_e', 'EX_o2_e', 'EX_hco3_e', 'EX_nh4_e', 'EX_pi_e', 'EX_h2o_e', 'EX_h_e']
    for ex_id in key_exchanges:
        if ex_id in ref_model.reactions:
            rxn = ref_model.reactions.get_by_id(ex_id)
            print(f"  {ex_id:20s}: [{rxn.lower_bound:>8.1f}, {rxn.upper_bound:>8.1f}]")
    
    # ATPM=0으로 설정하고 테스트
    print(f"\n[ATPM=0으로 설정하고 테스트]")
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    ref_model.objective = 'Growth'
    solution = ref_model.optimize()
    
    print(f"\n[FBA 결과]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    
    key_reactions = ['ACS_ADP', 'ACS', 'CS', 'ICL', 'MALS', 'SUCDi', 'ATPS4rpp', 'NADH16pp']
    print(f"\n[주요 반응 플럭스]")
    active_count = 0
    for rxn_id in key_reactions:
        if rxn_id in ref_model.reactions:
            flux = solution.fluxes.get(rxn_id, 0.0)
            if abs(flux) > 1e-6:
                print(f"  {rxn_id:20s}: {flux:>12.6f}")
                active_count += 1
    
    if active_count == 0:
        print(f"  [문제] 모든 주요 반응 플럭스가 0")
    else:
        print(f"  [OK] {active_count}개 반응이 작동함")
    
    # Exchange 플럭스 확인
    print(f"\n[Exchange 플럭스]")
    for ex_id in key_exchanges:
        if ex_id in ref_model.reactions:
            flux = solution.fluxes.get(ex_id, 0.0)
            if abs(flux) > 1e-6:
                print(f"  {ex_id:20s}: {flux:>12.6f}")
    
    # ATPM=10으로 설정하고 비교
    print(f"\n[ATPM=10으로 설정하고 비교]")
    atpm_rxn.lower_bound = 10
    atpm_rxn.upper_bound = 1000
    
    solution_10 = ref_model.optimize()
    
    print(f"\n[FBA 결과 (ATPM=10)]")
    print(f"  상태: {solution_10.status}")
    print(f"  성장률: {solution_10.objective_value:.6f}")
    
    print(f"\n[주요 반응 플럭스 (ATPM=10)]")
    active_count_10 = 0
    for rxn_id in key_reactions:
        if rxn_id in ref_model.reactions:
            flux = solution_10.fluxes.get(rxn_id, 0.0)
            if abs(flux) > 1e-6:
                print(f"  {rxn_id:20s}: {flux:>12.6f}")
                active_count_10 += 1
    
    print(f"\n[비교]")
    print(f"  ATPM=0: {active_count}개 반응 작동")
    print(f"  ATPM=10: {active_count_10}개 반응 작동")

def main():
    test_reference_without_media()
    
    print("\n" + "="*80)
    print("결론")
    print("="*80)
    print(f"\n[핵심 질문]")
    print(f"  레퍼런스 모델에서는 FBA가 잘 돌아가는데, 신규 모델에서는 왜 안 될까?")
    print(f"  -> 레퍼런스 모델의 기본 설정에서 테스트 완료")

if __name__ == "__main__":
    main()
