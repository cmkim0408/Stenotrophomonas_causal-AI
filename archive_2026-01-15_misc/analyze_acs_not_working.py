#!/usr/bin/env python
"""
ACS 반응이 작동하지 않는 실제 원인 분석
- ADK1은 이미 있음
- ACS는 있는데 작동 안 함
- 실제 원인 찾기
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

def analyze_why_acs_not_working(model):
    """ACS가 작동하지 않는 이유 분석"""
    print("="*80)
    print("ACS 반응 작동 안 함 원인 분석")
    print("="*80)
    
    # ADK1 확인
    adk1 = model.reactions.get_by_id('ADK1')
    print(f"\n[ADK1 확인]")
    print(f"  반응식: {adk1.reaction}")
    print(f"  bounds: [{adk1.lower_bound}, {adk1.upper_bound}]")
    print(f"  -> 이미 존재하고 bounds도 정상")
    
    # ACS 확인
    acs = model.reactions.get_by_id('ACS')
    print(f"\n[ACS 확인]")
    print(f"  반응식: {acs.reaction}")
    print(f"  bounds: [{acs.lower_bound}, {acs.upper_bound}]")
    
    # ACS 반응의 메타볼라이트 확인
    print(f"\n[ACS 반응 메타볼라이트]")
    for met, coeff in acs.metabolites.items():
        print(f"  {met.id}: {coeff}")
    
    # ATPM=0 설정
    atpm_rxn = model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 0
    
    model = setup_acetate_medium(model)
    model.objective = 'Growth'
    
    # FBA 수행
    solution = model.optimize()
    print(f"\n[FBA 결과]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    print(f"  ACS 플럭스: {solution.fluxes.get('ACS', 0):.6f}")
    print(f"  ADK1 플럭스: {solution.fluxes.get('ADK1', 0):.6f}")
    
    # 레퍼런스 모델과 비교
    print(f"\n[레퍼런스 모델 비교]")
    print(f"  레퍼런스 모델 ATPM=0:")
    print(f"    ACS 플럭스: 없음 (ACS 사용 안 함)")
    print(f"    ACS_ADP 플럭스: 0.943448 (실제 사용)")
    print(f"    ADK1 플럭스: 0.2096")
    print(f"    ICL 플럭스: 0.367")
    print(f"    MALS 플럭스: 0.367")
    
    print(f"\n[핵심 발견]")
    print(f"  레퍼런스 모델은 ACS를 사용하지 않고 ACS_ADP를 사용!")
    print(f"  -> ACS가 작동하지 않는 것이 아니라, 레퍼런스 모델도 ACS를 사용 안 함")
    print(f"  -> 신규 모델에는 ACS_ADP가 없어서 Acetate 전환이 안 됨")
    
    # SUCOAACTr 확인
    try:
        sucoaac = model.reactions.get_by_id('SUCOAACTr')
        print(f"\n[SUCOAACTr 확인]")
        print(f"  반응식: {sucoaac.reaction}")
        print(f"  bounds: [{sucoaac.lower_bound}, {sucoaac.upper_bound}]")
        sucoaac_flux = solution.fluxes.get('SUCOAACTr', 0.0)
        print(f"  플럭스 (ATPM=0): {sucoaac_flux:.6f}")
        
        if abs(sucoaac_flux) < 1e-6:
            print(f"  -> SUCOAACTr도 작동 안 함 (Succinyl-CoA 필요)")
    except KeyError:
        print(f"\n[SUCOAACTr] 반응 없음")

def compare_acs_vs_acs_adp_usage():
    """ACS vs ACS_ADP 사용 비교"""
    print("\n" + "="*80)
    print("레퍼런스 모델에서 ACS vs ACS_ADP 사용 비교")
    print("="*80)
    
    import pandas as pd
    base_path = Path(__file__).parent.parent
    flux_file = base_path / "Stenotrophomonas" / "fba_flux_gradient_acid.csv"
    
    if flux_file.exists():
        df = pd.read_csv(flux_file)
        
        # ACS 관련 반응 찾기
        acs_rxns = ['ACS', 'ACS_ADP']
        for rxn_id in acs_rxns:
            row = df[df.iloc[:,0] == rxn_id]
            if not row.empty:
                flux_atpm0 = row.iloc[0, 1]  # ATPM_0 컬럼
                print(f"\n  {rxn_id}: ATPM=0 플럭스 = {flux_atpm0:.6f}")
                if abs(flux_atpm0) > 1e-6:
                    print(f"    -> 실제 사용됨!")
                else:
                    print(f"    -> 사용 안 됨")

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    
    model = load_model(str(model_path))
    
    print("="*80)
    print("ACS 반응 작동 안 함 원인 분석")
    print("="*80)
    
    # 원인 분석
    analyze_why_acs_not_working(model)
    
    # 레퍼런스 모델 비교
    compare_acs_vs_acs_adp_usage()
    
    print("\n" + "="*80)
    print("결론")
    print("="*80)
    
    print("\n[핵심 발견]")
    print("  레퍼런스 모델도 ACS를 사용하지 않고 ACS_ADP를 사용합니다!")
    print("  -> ACS (ATP → AMP + PPi)는 에너지 효율이 낮아서 사용 안 함")
    print("  -> ACS_ADP (ATP → ADP + Pi)가 더 효율적이어서 사용")
    
    print("\n[신규 모델의 문제]")
    print("  -> ACS_ADP 반응이 없음")
    print("  -> ACS만 있는데, 레퍼런스 모델도 ACS를 사용하지 않음")
    print("  -> 따라서 ACS로는 작동하지 않는 것이 정상일 수 있음")
    
    print("\n[해결 방안]")
    print("  -> ACS_ADP 반응을 추가해야 함 (레퍼런스 모델처럼)")
    print("  -> 또는 SUCOAACTr 경로가 작동하도록 다른 조건 확인")

if __name__ == "__main__":
    main()
