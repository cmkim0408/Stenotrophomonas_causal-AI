#!/usr/bin/env python
"""
신규 모델에서 탄소 플럭스 경로 분석
- TCA cycle
- Glyoxylate shunt
- Acetate 대사 경로
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
    # 이름으로 찾기
    for name in ['Growth', 'BIOMASS', 'BIOMASS_Ecoli_core', 'R_BIOMASS']:
        try:
            return model.reactions.get_by_id(name)
        except KeyError:
            continue
    
    # 패턴으로 찾기
    for rxn in model.reactions:
        if 'biomass' in rxn.id.lower() or 'growth' in rxn.id.lower():
            return rxn
    
    return None

def setup_acetate_medium(model):
    """Acetate 미디어 설정"""
    # 모든 exchange 차단 (upper_bound 먼저 설정)
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    # Acetate 허용
    try:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.upper_bound = 1000
        ex_ac.lower_bound = -1000
    except KeyError:
        # 다른 이름으로 찾기
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

def find_tca_reactions(model):
    """TCA cycle 관련 반응 찾기"""
    tca_keywords = {
        'CS': ['citrate synthase', 'CS'],
        'ACONT': ['aconitase', 'ACONT'],
        'ICDHxr': ['isocitrate dehydrogenase', 'ICDH'],
        'AKGD': ['alpha-ketoglutarate dehydrogenase', 'AKGD'],
        'SUCOAS': ['succinyl-CoA synthetase', 'SUCOAS'],
        'SUCDi': ['succinate dehydrogenase', 'SUCD'],
        'FUM': ['fumarase', 'FUM'],
        'MDH': ['malate dehydrogenase', 'MDH'],
        'ICL': ['isocitrate lyase', 'ICL'],
        'MALS': ['malate synthase', 'MALS', 'malate synthase']
    }
    
    tca_reactions = {}
    for rxn in model.reactions:
        rxn_name = (rxn.name or '').lower()
        rxn_id = rxn.id.upper()
        
        for key, keywords in tca_keywords.items():
            if any(kw.lower() in rxn_name or kw in rxn_id for kw in keywords):
                if key not in tca_reactions:
                    tca_reactions[key] = []
                tca_reactions[key].append({
                    'id': rxn.id,
                    'name': rxn.name if rxn.name else '',
                    'equation': rxn.reaction
                })
    
    return tca_reactions

def analyze_tca_cycle(model, solution):
    """TCA cycle 플럭스 분석"""
    print("\n" + "="*80)
    print("TCA Cycle 분석")
    print("="*80)
    
    tca_reactions = find_tca_reactions(model)
    
    # 핵심 TCA cycle 반응
    core_tca = ['CS', 'ACONT', 'ICDHxr', 'AKGD', 'SUCOAS', 'SUCDi', 'FUM', 'MDH']
    glyoxylate = ['ICL', 'MALS']
    
    print("\n[핵심 TCA Cycle 반응]")
    for rxn_type in core_tca:
        if rxn_type in tca_reactions:
            for rxn_info in tca_reactions[rxn_type]:
                rxn_id = rxn_info['id']
                flux = solution.fluxes.get(rxn_id, 0.0) if solution else 0.0
                print(f"\n  {rxn_id} ({rxn_type})")
                if rxn_info['name']:
                    print(f"    이름: {rxn_info['name']}")
                print(f"    반응식: {rxn_info['equation']}")
                if solution:
                    print(f"    플럭스: {flux:.6f}")
        else:
            print(f"\n  [{rxn_type}] 반응 없음")
    
    print("\n[Glyoxylate Shunt 반응]")
    for rxn_type in glyoxylate:
        if rxn_type in tca_reactions:
            for rxn_info in tca_reactions[rxn_type]:
                rxn_id = rxn_info['id']
                flux = solution.fluxes.get(rxn_id, 0.0) if solution else 0.0
                print(f"\n  {rxn_id} ({rxn_type})")
                if rxn_info['name']:
                    print(f"    이름: {rxn_info['name']}")
                print(f"    반응식: {rxn_info['equation']}")
                if solution:
                    print(f"    플럭스: {flux:.6f}")
        else:
            print(f"\n  [{rxn_type}] 반응 없음")
    
    return tca_reactions

def analyze_acetate_pathway(model, solution):
    """Acetate 대사 경로 분석"""
    print("\n" + "="*80)
    print("Acetate 대사 경로 분석")
    print("="*80)
    
    acetate_rxns = []
    for rxn in model.reactions:
        if 'ac_c' in [m.id for m in rxn.reactants + rxn.products]:
            flux = solution.fluxes.get(rxn.id, 0.0) if solution else 0.0
            if abs(flux) > 1e-6 or 'accoa' in str(rxn.reaction).lower():
                acetate_rxns.append({
                    'id': rxn.id,
                    'name': rxn.name if rxn.name else '',
                    'equation': rxn.reaction,
                    'flux': flux
                })
    
    print(f"\n[Acetate 관련 활성 반응] (플럭스 > 1e-6)")
    for rxn_info in sorted(acetate_rxns, key=lambda x: abs(x['flux']), reverse=True)[:10]:
        print(f"\n  {rxn_info['id']}: {rxn_info['flux']:.6f}")
        if rxn_info['name']:
            print(f"    이름: {rxn_info['name']}")
        print(f"    반응식: {rxn_info['equation']}")
    
    return acetate_rxns

def analyze_carbon_flow(model, solution):
    """탄소 플럭스 흐름 분석"""
    print("\n" + "="*80)
    print("탄소 플럭스 흐름 분석")
    print("="*80)
    
    # 주요 탄소 대사물질 추적
    key_metabolites = ['ac_c', 'accoa_c', 'cit_c', 'icit_c', 'akg_c', 'succ_c', 'fum_c', 'mal_c', 'oaa_c']
    
    print("\n[주요 탄소 대사물질 플럭스]")
    for met_id in key_metabolites:
        try:
            met = model.metabolites.get_by_id(met_id)
            producing_rxns = []
            consuming_rxns = []
            
            for rxn in met.reactions:
                flux = solution.fluxes.get(rxn.id, 0.0) if solution else 0.0
                if abs(flux) < 1e-6:
                    continue
                
                # 생성/소비 판단
                coeff = rxn.metabolites[met]
                if coeff > 0 and flux > 0:  # 생성
                    producing_rxns.append((rxn.id, flux * coeff))
                elif coeff < 0 and flux < 0:  # 생성 (역방향)
                    producing_rxns.append((rxn.id, abs(flux * coeff)))
                elif coeff > 0 and flux < 0:  # 소비 (역방향)
                    consuming_rxns.append((rxn.id, abs(flux * coeff)))
                elif coeff < 0 and flux > 0:  # 소비
                    consuming_rxns.append((rxn.id, abs(flux * coeff)))
            
            if producing_rxns or consuming_rxns:
                print(f"\n  {met_id} ({met.name if met.name else '이름 없음'})")
                if producing_rxns:
                    print(f"    생성 반응:")
                    for rxn_id, net_flux in sorted(producing_rxns, key=lambda x: x[1], reverse=True)[:3]:
                        print(f"      {rxn_id}: {net_flux:.6f}")
                if consuming_rxns:
                    print(f"    소비 반응:")
                    for rxn_id, net_flux in sorted(consuming_rxns, key=lambda x: x[1], reverse=True)[:3]:
                        print(f"      {rxn_id}: {net_flux:.6f}")
        
        except KeyError:
            print(f"\n  {met_id}: 메타볼라이트 없음")

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    
    model = load_model(str(model_path))
    
    # Biomass 반응 찾기
    biomass_rxn = find_biomass_reaction(model)
    if biomass_rxn:
        print(f"\n[Biomass 반응] {biomass_rxn.id}")
        model.objective = biomass_rxn.id
    else:
        print("\n[WARNING] Biomass 반응을 찾을 수 없습니다.")
        return
    
    # 미디어 설정
    model = setup_acetate_medium(model)
    print("\n[미디어 설정] Acetate minimal medium")
    
    # FBA 수행
    print("\n" + "="*80)
    print("FBA 수행")
    print("="*80)
    
    try:
        solution = model.optimize()
        print(f"\n[FBA 결과]")
        print(f"  상태: {solution.status}")
        print(f"  목적 함수 값 (성장률): {solution.objective_value:.6f}")
        
        if solution.status == 'optimal' and solution.objective_value > 0:
            print("\n[OK] 모델이 성장합니다!")
            
            # TCA cycle 분석
            tca_reactions = analyze_tca_cycle(model, solution)
            
            # Acetate 경로 분석
            acetate_rxns = analyze_acetate_pathway(model, solution)
            
            # 탄소 플럭스 분석
            analyze_carbon_flow(model, solution)
            
        else:
            print("\n[ERROR] 모델이 성장하지 않습니다.")
            print("  TCA cycle 반응만 확인합니다.")
            tca_reactions = analyze_tca_cycle(model, None)
    
    except Exception as e:
        print(f"\n[ERROR] FBA 수행 중 오류: {e}")
        print("  TCA cycle 반응만 확인합니다.")
        tca_reactions = analyze_tca_cycle(model, None)

if __name__ == "__main__":
    main()
