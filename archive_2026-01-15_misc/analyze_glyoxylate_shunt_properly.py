#!/usr/bin/env python
"""
Glyoxylate shunt 활성화 여부를 제대로 분석
- ATPM 값 확인
- TCA vs Glyoxylate 플럭스 비교
- 레퍼런스 모델과 비교
"""

import cobra
import pandas as pd
from pathlib import Path

def load_model(model_path):
    print(f"모델 로드 중: {model_path}")
    model = cobra.io.read_sbml_model(model_path)
    print(f"[OK] 모델 로드 완료: {model.id}")
    return model

def find_biomass_reaction(model):
    """Biomass 반응 찾기"""
    for name in ['Growth', 'BIOMASS', 'BIOMASS_Ecoli_core', 'R_BIOMASS']:
        try:
            return model.reactions.get_by_id(name)
        except KeyError:
            continue
    
    for rxn in model.reactions:
        if 'biomass' in rxn.id.lower() or 'growth' in rxn.id.lower():
            return rxn
    
    return None

def find_atpm_reaction(model):
    """ATPM 반응 찾기"""
    atpm_keywords = ['ATPM', 'ATPM', 'ATP_maintenance', 'atp_maintenance']
    for rxn in model.reactions:
        if any(kw.lower() in rxn.id.lower() for kw in atpm_keywords):
            return rxn
    return None

def setup_acetate_medium(model):
    """Acetate 미디어 설정"""
    # 모든 exchange 차단
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    # Acetate 허용
    try:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.upper_bound = 1000
        ex_ac.lower_bound = -1000
    except KeyError:
        for rxn in model.exchanges:
            if 'ac_e' in rxn.id.lower() and 'EX' in rxn.id:
                rxn.upper_bound = 1000
                rxn.lower_bound = -1000
                break
    
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

def analyze_key_reactions(model, solution):
    """주요 반응 플럭스 분석"""
    print("\n" + "="*80)
    print("주요 반응 플럭스 분석")
    print("="*80)
    
    key_reactions = {
        'Acetate 전환': ['ACS', 'ACS_ADP', 'SUCOAACTr'],
        'TCA Cycle': ['CS', 'ACONT', 'ICDHx', 'ICDHyr', 'AKGDH', 'SUCOAS', 'SUCD', 'SUCDi', 'FUM', 'MDH'],
        'Glyoxylate Shunt': ['ICL', 'MALS'],
        '에너지': ['ATPM']
    }
    
    for category, rxn_ids in key_reactions.items():
        print(f"\n[{category}]")
        for rxn_id in rxn_ids:
            try:
                rxn = model.reactions.get_by_id(rxn_id)
                flux = solution.fluxes.get(rxn_id, 0.0) if solution else 0.0
                print(f"  {rxn_id}: {flux:.6f}")
                if abs(flux) > 1e-6:
                    print(f"    반응식: {rxn.reaction}")
            except KeyError:
                print(f"  {rxn_id}: 반응 없음")

def analyze_glyoxylate_vs_tca(model, solution):
    """Glyoxylate shunt vs TCA cycle 비교"""
    print("\n" + "="*80)
    print("Glyoxylate Shunt vs TCA Cycle 비교")
    print("="*80)
    
    # ICL과 ICDHx 플럭스 비교
    try:
        icl = model.reactions.get_by_id('ICL')
        icl_flux = solution.fluxes.get('ICL', 0.0) if solution else 0.0
        print(f"\nICL (Isocitrate Lyase): {icl_flux:.6f}")
        print(f"  반응식: {icl.reaction}")
    except KeyError:
        print("\nICL: 반응 없음")
        icl_flux = 0.0
    
    try:
        icdhx = model.reactions.get_by_id('ICDHx')
        icdhx_flux = solution.fluxes.get('ICDHx', 0.0) if solution else 0.0
        print(f"\nICDHx (Isocitrate Dehydrogenase): {icdhx_flux:.6f}")
        print(f"  반응식: {icdhx.reaction}")
    except KeyError:
        print("\nICDHx: 반응 없음")
        icdhx_flux = 0.0
    
    try:
        mals = model.reactions.get_by_id('MALS')
        mals_flux = solution.fluxes.get('MALS', 0.0) if solution else 0.0
        print(f"\nMALS (Malate Synthase): {mals_flux:.6f}")
        print(f"  반응식: {mals.reaction}")
    except KeyError:
        print("\nMALS: 반응 없음")
        mals_flux = 0.0
    
    print("\n[분석]")
    if abs(icl_flux) > 1e-6 and abs(mals_flux) > 1e-6:
        print("  -> Glyoxylate shunt가 활성화되어 있습니다!")
        print(f"     ICL 플럭스: {icl_flux:.6f}")
        print(f"     MALS 플럭스: {mals_flux:.6f}")
    else:
        print("  -> Glyoxylate shunt가 비활성화되어 있습니다.")
        if abs(icdhx_flux) > 1e-6:
            print(f"     대신 ICDHx가 활성화됨 (플럭스: {icdhx_flux:.6f})")
            print("     -> TCA cycle 경로 사용 (탄소 손실 발생)")

def check_atpm_value(model):
    """ATPM 값 확인"""
    print("\n" + "="*80)
    print("ATPM (ATP Maintenance) 확인")
    print("="*80)
    
    atpm_rxn = find_atpm_reaction(model)
    if atpm_rxn:
        print(f"\nATPM 반응: {atpm_rxn.id}")
        print(f"  하한: {atpm_rxn.lower_bound}")
        print(f"  상한: {atpm_rxn.upper_bound}")
        print(f"  반응식: {atpm_rxn.reaction}")
        return atpm_rxn
    else:
        print("\nATPM 반응을 찾을 수 없습니다.")
        return None

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    
    model = load_model(str(model_path))
    
    # Biomass 반응 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("\n[ERROR] Biomass 반응을 찾을 수 없습니다.")
        return
    
    print(f"\n[Biomass 반응] {biomass_rxn.id}")
    model.objective = biomass_rxn.id
    
    # ATPM 확인
    atpm_rxn = check_atpm_value(model)
    
    # 미디어 설정
    model = setup_acetate_medium(model)
    
    # FBA 수행
    print("\n" + "="*80)
    print("FBA 수행")
    print("="*80)
    
    solution = model.optimize()
    print(f"\n[FBA 결과]")
    print(f"  상태: {solution.status}")
    print(f"  목적 함수 값: {solution.objective_value:.6f}")
    
    if solution.status == 'optimal':
        # ATPM 플럭스 확인
        if atpm_rxn:
            atpm_flux = solution.fluxes.get(atpm_rxn.id, 0.0)
            print(f"\n  ATPM 플럭스: {atpm_flux:.6f}")
            if abs(atpm_flux) < 10:
                print("    -> ATPM 값이 작습니다. Glyoxylate shunt가 메인이 되어야 합니다.")
        
        # 주요 반응 분석
        analyze_key_reactions(model, solution)
        
        # Glyoxylate vs TCA 비교
        analyze_glyoxylate_vs_tca(model, solution)
        
    else:
        print("\n[ERROR] FBA 최적화 실패")

if __name__ == "__main__":
    main()
