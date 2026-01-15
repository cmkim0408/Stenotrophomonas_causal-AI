#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ACS 반응이 사용되지 않는 이유 분석

ACS: ac_c + atp_c + coa_c --> accoa_c + amp_c + ppi_c
- ACS_ADP와 SUCOAACTr를 제거했는데도 ACS가 사용되지 않음
- 대신 AADb/AADa 경로가 사용됨
"""

import cobra
from pathlib import Path
import sys

def load_model(model_path):
    try:
        model = cobra.io.read_sbml_model(str(model_path))
        return model
    except Exception as e:
        print(f"[ERROR] 모델 로드 실패: {e}")
        sys.exit(1)

def setup_media_forced(model):
    """배지 조건을 강제로 고정"""
    if 'EX_ac_e' in model.reactions:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.lower_bound = -19.0
        ex_ac.upper_bound = -19.0
    
    if 'EX_o2_e' in model.reactions:
        ex_o2 = model.reactions.get_by_id('EX_o2_e')
        ex_o2.lower_bound = -100.0
        ex_o2.upper_bound = 1000.0
    
    essential_exchanges = {
        'EX_nh4_e': (-1000.0, 1000.0),
        'EX_pi_e': (-1000.0, 1000.0),
        'EX_so4_e': (-1000.0, 1000.0),
        'EX_mg2_e': (-1000.0, 1000.0),
        'EX_k_e': (-1000.0, 1000.0),
        'EX_na1_e': (-1000.0, 1000.0),
        'EX_fe2_e': (-1000.0, 1000.0),
        'EX_fe3_e': (-1000.0, 1000.0),
        'EX_h2o_e': (-1000.0, 1000.0),
        'EX_h_e': (-1000.0, 1000.0),
        'EX_co2_e': (-1000.0, 1000.0),
        'EX_hco3_e': (-1000.0, 1000.0),
        'EX_nac_e': (-1000.0, 1000.0),
        'EX_ncam_e': (-1000.0, 1000.0),
    }
    
    for ex_id, (lb, ub) in essential_exchanges.items():
        if ex_id in model.reactions:
            ex_rxn = model.reactions.get_by_id(ex_id)
            ex_rxn.lower_bound = lb
            ex_rxn.upper_bound = ub
    
    return model

def check_acs_reaction(model):
    """ACS 반응 확인"""
    print("\n" + "="*70)
    print("ACS 반응 확인")
    print("="*70)
    
    if 'ACS' not in model.reactions:
        print("[ERROR] ACS 반응이 모델에 없습니다.")
        return None
    
    acs = model.reactions.get_by_id('ACS')
    print(f"\n[ACS 반응]")
    print(f"  ID: {acs.id}")
    print(f"  이름: {acs.name if acs.name else ''}")
    print(f"  반응식: {acs.reaction}")
    print(f"  bounds: [{acs.lower_bound}, {acs.upper_bound}]")
    
    return acs

def check_aadb_aada_pathway(model, solution):
    """AADb/AADa 경로 확인"""
    print("\n" + "="*70)
    print("AADb/AADa 경로 확인")
    print("="*70)
    
    if 'AADb' in model.reactions:
        aadb = model.reactions.get_by_id('AADb')
        aadb_flux = solution.fluxes.get('AADb', 0.0)
        print(f"\n[AADb]")
        print(f"  반응식: {aadb.reaction}")
        print(f"  플럭스: {aadb_flux:.6f}")
        print(f"  bounds: [{aadb.lower_bound}, {aadb.upper_bound}]")
        print(f"  설명: Acetate + ATP → Acetyl-adenylate + PPi")
    
    if 'AADa' in model.reactions:
        aada = model.reactions.get_by_id('AADa')
        aada_flux = solution.fluxes.get('AADa', 0.0)
        print(f"\n[AADa]")
        print(f"  반응식: {aada.reaction}")
        print(f"  플럭스: {aada_flux:.6f}")
        print(f"  bounds: [{aada.lower_bound}, {aada.upper_bound}]")
        print(f"  설명: Acetyl-adenylate + CoA → Acetyl-CoA + AMP")
    
    print(f"\n[경로 비교]")
    print(f"  AADb/AADa 경로: Acetate + ATP + CoA → Acetyl-CoA + AMP + PPi (2단계)")
    print(f"  ACS 경로: Acetate + ATP + CoA → Acetyl-CoA + AMP + PPi (1단계)")
    print(f"  -> 두 경로는 동일한 결과이지만, ACS가 더 직접적")

def test_force_acs(model):
    """ACS를 강제로 사용하도록 테스트"""
    print("\n" + "="*70)
    print("ACS 강제 사용 테스트")
    print("="*70)
    
    # ACS의 하한을 설정하여 강제로 사용
    if 'ACS' in model.reactions:
        acs = model.reactions.get_by_id('ACS')
        original_lb = acs.lower_bound
        
        # ACS를 최소 10.0으로 설정
        acs.lower_bound = 10.0
        
        solution = model.optimize()
        
        print(f"\n[ACS 하한 = {acs.lower_bound}]")
        print(f"  상태: {solution.status}")
        print(f"  성장률: {solution.objective_value:.6f}")
        
        acs_flux = solution.fluxes.get('ACS', 0.0)
        print(f"  ACS 플럭스: {acs_flux:.6f}")
        
        # AADb/AADa 플럭스 확인
        if 'AADb' in model.reactions:
            aadb_flux = solution.fluxes.get('AADb', 0.0)
            print(f"  AADb 플럭스: {aadb_flux:.6f}")
        
        if 'AADa' in model.reactions:
            aada_flux = solution.fluxes.get('AADa', 0.0)
            print(f"  AADa 플럭스: {aada_flux:.6f}")
        
        # 원래 값 복원
        acs.lower_bound = original_lb
        
        return solution

def compare_pathways(model):
    """경로 효율성 비교"""
    print("\n" + "="*70)
    print("경로 효율성 비교")
    print("="*70)
    
    print(f"\n[1. AADb/AADa 경로]")
    print(f"  AADb: ac_c + atp_c <=> aad_c + ppi_c")
    print(f"  AADa: aad_c + coa_c <=> accoa_c + amp_c")
    print(f"  합계: ac_c + atp_c + coa_c <=> accoa_c + amp_c + ppi_c")
    print(f"  -> 2단계, 중간체 (aad_c) 생성")
    
    print(f"\n[2. ACS 경로]")
    print(f"  ACS: ac_c + atp_c + coa_c --> accoa_c + amp_c + ppi_c")
    print(f"  -> 1단계, 직접 전환")
    
    print(f"\n[3. ACS_ADP 경로 (제거됨)]")
    print(f"  ACS_ADP: ac_c + atp_c + coa_c <=> accoa_c + adp_c + pi_c")
    print(f"  -> 1단계, ATP → ADP (더 효율적)")

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_without_ACS_ADP_SUCOAACTr.xml"
    
    print("="*70)
    print("ACS 반응이 사용되지 않는 이유 분석")
    print("="*70)
    
    # 모델 로드
    model = load_model(model_path)
    model = setup_media_forced(model)
    
    # ATPM=0 설정
    if 'ATPM' in model.reactions:
        atpm = model.reactions.get_by_id('ATPM')
        atpm.lower_bound = 0.0
        atpm.upper_bound = 0.0
    
    # FBA 실행
    solution = model.optimize()
    print(f"\n[FBA 결과]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    
    # ACS 반응 확인
    acs = check_acs_reaction(model)
    
    # AADb/AADa 경로 확인
    check_aadb_aada_pathway(model, solution)
    
    # 경로 비교
    compare_pathways(model)
    
    # ACS 강제 사용 테스트
    test_force_acs(model)
    
    # 결론
    print("\n" + "="*70)
    print("결론")
    print("="*70)
    
    acs_flux = solution.fluxes.get('ACS', 0.0) if 'ACS' in model.reactions else 0.0
    aadb_flux = solution.fluxes.get('AADb', 0.0) if 'AADb' in model.reactions else 0.0
    
    print(f"\n[현재 상황]")
    print(f"  ACS 플럭스: {acs_flux:.6f}")
    print(f"  AADb 플럭스: {aadb_flux:.6f}")
    
    if abs(acs_flux) < 1e-6 and abs(aadb_flux) > 1e-6:
        print(f"\n[이유]")
        print(f"  모델이 AADb/AADa 경로를 선호합니다.")
        print(f"  ACS와 AADb/AADa는 동일한 결과를 생성하지만,")
        print(f"  FBA가 AADb/AADa를 선택하는 이유는:")
        print(f"  1. 중간체 (aad_c)가 다른 경로에서도 사용될 수 있음")
        print(f"  2. 또는 반응 bounds/제약 조건 때문일 수 있음")
        print(f"  -> ACS를 사용하려면 ACS_ADP를 추가하는 것이 더 나을 수 있습니다.")
    
    return model, solution

if __name__ == "__main__":
    model, solution = main()
