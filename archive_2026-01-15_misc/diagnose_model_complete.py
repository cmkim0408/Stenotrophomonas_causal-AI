#!/usr/bin/env python
"""
모델 완전 진단 - 3단계 통합
1. Biomass 반응 검증
2. 기본 대사 경로 연결성 점검
3. 모델 구조 검토
"""

import cobra
import pandas as pd
from cobra.flux_analysis import find_blocked_reactions

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    print(f"[OK] 모델 로드: {model.id}")
    return model

def find_biomass_reaction(model):
    """Biomass 반응 찾기"""
    for name in ['Growth', 'BIOMASS', 'BIOMASS_Ecoli_core']:
        try:
            return model.reactions.get_by_id(name)
        except KeyError:
            continue
    return None

def analyze_biomass_reaction(model, biomass_rxn):
    """1단계: Biomass 반응 분석"""
    print("\n" + "="*70)
    print("1단계: Biomass 반응 분석")
    print("="*70)
    
    print(f"\nBiomass 반응: {biomass_rxn.id}")
    print(f"이름: {biomass_rxn.name}")
    print(f"\n반응식:")
    print(f"  {biomass_rxn.reaction}")
    
    metabolites = biomass_rxn.metabolites
    
    print(f"\n총 구성 요소: {len(metabolites)}개")
    
    # 주요 구성 요소 분석 (계수 절대값 기준)
    sorted_mets = sorted(metabolites.items(), key=lambda x: abs(x[1]), reverse=True)
    
    print(f"\n주요 구성 요소 (계수 상위 20개):")
    print(f"{'대사물질 ID':<25} {'계수':<15} {'구분':<20}")
    print("-" * 65)
    
    component_analysis = []
    
    # 카테고리 분류
    amino_acids = ['ala', 'arg', 'asn', 'asp', 'cys', 'gln', 'glu', 'gly', 
                   'his', 'ile', 'leu', 'lys', 'met', 'phe', 'pro', 'ser', 
                   'thr', 'trp', 'tyr', 'val']
    
    nucleotides = ['atp', 'ctp', 'gtp', 'utp', 'datp', 'dctp', 'dgtp', 'dttp']
    
    cofactors = ['coa', 'nad', 'nadp', 'fad', 'thf', 'mlthf', '10fthf', 'pydx5p', 'ribflv']
    
    for met, coeff in sorted_mets[:20]:
        met_id_lower = met.id.lower()
        
        category = "기타"
        if any(aa in met_id_lower for aa in amino_acids):
            category = "아미노산"
        elif any(nt in met_id_lower for nt in nucleotides):
            category = "뉴클레오티드"
        elif any(cf in met_id_lower for cf in cofactors):
            category = "보조인자"
        elif 'h2o' in met_id_lower:
            category = "물"
        elif 'h_c' in met_id_lower or 'h_e' in met_id_lower:
            category = "양성자"
        elif 'pi' in met_id_lower or 'ppi' in met_id_lower:
            category = "인산"
        elif any(ion in met_id_lower for ion in ['ca2', 'k', 'mg2', 'fe2', 'fe3', 'cu2', 'mn2', 'zn2', 'co2']):
            category = "이온"
        
        component_analysis.append({
            'Metabolite_ID': met.id,
            'Coefficient': coeff,
            'Category': category
        })
        
        print(f"{met.id:<25} {coeff:<15.6f} {category:<20}")
    
    return component_analysis

def check_basic_pathways(model):
    """2단계: 기본 대사 경로 연결성 점검"""
    print("\n" + "="*70)
    print("2단계: 기본 대사 경로 연결성 점검")
    print("="*70)
    
    # 포도당 → Pyruvate 경로 (Glycolysis)
    print("\n[Glycolysis 경로]")
    glycolysis_reactions = {
        'HEX1': 'Hexokinase (Glucose → G6P)',
        'PGI': 'Phosphoglucose Isomerase (G6P → F6P)',
        'PFK': 'Phosphofructokinase (F6P → FBP)',
        'FBA': 'Fructose-bisphosphate Aldolase (FBP → G3P + DHAP)',
        'TPI': 'Triosephosphate Isomerase (DHAP ↔ G3P)',
        'GAPD': 'Glyceraldehyde-3-phosphate Dehydrogenase (G3P → 1,3BPG)',
        'PGK': 'Phosphoglycerate Kinase (1,3BPG → 3PG + ATP)',
        'PGM': 'Phosphoglycerate Mutase (3PG → 2PG)',
        'ENO': 'Enolase (2PG → PEP)',
        'PYK': 'Pyruvate Kinase (PEP → Pyruvate + ATP)'
    }
    
    found_glycolysis = []
    missing_glycolysis = []
    
    for rxn_id, rxn_name in glycolysis_reactions.items():
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            found_glycolysis.append(rxn_id)
            print(f"  [OK] {rxn_id}: {rxn_name}")
        except KeyError:
            missing_glycolysis.append(rxn_id)
            print(f"  [MISSING] {rxn_id}: {rxn_name}")
    
    print(f"\nGlycolysis: {len(found_glycolysis)}/{len(glycolysis_reactions)} 반응 존재")
    
    # Pyruvate → Acetyl-CoA 경로
    print("\n[Pyruvate → Acetyl-CoA 경로]")
    pyr_to_accoa = {
        'PDH': 'Pyruvate Dehydrogenase (Pyruvate → Acetyl-CoA)',
        'PTAr': 'Phosphotransacetylase (Acetyl-P → Acetyl-CoA)',
        'ACKr': 'Acetate Kinase (Acetyl-P + ADP ↔ Acetate + ATP)'
    }
    
    for rxn_id, rxn_name in pyr_to_accoa.items():
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            print(f"  [OK] {rxn_id}: {rxn_name}")
        except KeyError:
            print(f"  [MISSING] {rxn_id}: {rxn_name}")
    
    # TCA cycle 경로
    print("\n[TCA Cycle 경로]")
    tca_reactions = {
        'CS': 'Citrate Synthase',
        'ACONT': 'Aconitase',
        'ICDHx': 'Isocitrate Dehydrogenase (NAD+)',
        'AKGDH': 'α-Ketoglutarate Dehydrogenase',
        'SUCOAS': 'Succinyl-CoA Synthetase',
        'SUCD': 'Succinate Dehydrogenase',
        'FUM': 'Fumarase',
        'MDH': 'Malate Dehydrogenase'
    }
    
    found_tca = []
    for rxn_id, rxn_name in tca_reactions.items():
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            found_tca.append(rxn_id)
            print(f"  [OK] {rxn_id}: {rxn_name}")
        except KeyError:
            print(f"  [MISSING] {rxn_id}: {rxn_name}")
    
    print(f"\nTCA Cycle: {len(found_tca)}/{len(tca_reactions)} 반응 존재")
    
    # Glyoxylate shunt
    print("\n[Glyoxylate Shunt 경로]")
    gly_reactions = {
        'ICL': 'Isocitrate Lyase',
        'MALS': 'Malate Synthase'
    }
    
    found_gly = []
    for rxn_id, rxn_name in gly_reactions.items():
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            found_gly.append(rxn_id)
            print(f"  [OK] {rxn_id}: {rxn_name}")
        except KeyError:
            print(f"  [MISSING] {rxn_id}: {rxn_name}")
    
    print(f"\nGlyoxylate Shunt: {len(found_gly)}/{len(gly_reactions)} 반응 존재")
    
    return {
        'glycolysis': {'found': found_glycolysis, 'missing': missing_glycolysis},
        'tca': {'found': found_tca},
        'glyoxylate': {'found': found_gly}
    }

def check_pathway_connectivity(model):
    """경로 연결성 확인"""
    print("\n" + "="*70)
    print("경로 연결성 확인 (대사물질 기반)")
    print("="*70)
    
    # 주요 대사물질 확인
    key_metabolites = {
        'glc__D_c': 'Glucose',
        'pyr_c': 'Pyruvate',
        'accoa_c': 'Acetyl-CoA',
        'cit_c': 'Citrate',
        'icit_c': 'Isocitrate',
        'akg_c': 'α-Ketoglutarate',
        'succoa_c': 'Succinyl-CoA',
        'succ_c': 'Succinate',
        'fum_c': 'Fumarate',
        'mal__L_c': 'Malate',
        'oaa_c': 'Oxaloacetate',
        'glx_c': 'Glyoxylate'
    }
    
    print("\n주요 대사물질 존재 및 연결성:")
    metabolite_status = []
    
    for met_id, met_name in key_metabolites.items():
        try:
            met = model.metabolites.get_by_id(met_id)
            reactions = [r.id for r in met.reactions]
            metabolite_status.append({
                'Metabolite_ID': met_id,
                'Name': met_name,
                'Exists': True,
                'Reaction_Count': len(reactions),
                'Sample_Reactions': ', '.join(reactions[:5])
            })
            print(f"  [OK] {met_name} ({met_id}): {len(reactions)}개 반응 참여")
        except KeyError:
            metabolite_status.append({
                'Metabolite_ID': met_id,
                'Name': met_name,
                'Exists': False,
                'Reaction_Count': 0,
                'Sample_Reactions': ''
            })
            print(f"  [MISSING] {met_name} ({met_id})")
    
    return metabolite_status

def check_blocked_reactions_analysis(model):
    """3단계: 블록된 반응 분석"""
    print("\n" + "="*70)
    print("3단계: 블록된 반응 분석")
    print("="*70)
    
    # Exchange 설정 (포도당 허용)
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    # 포도당 허용
    try:
        ex_glc = model.reactions.get_by_id('EX_glc__D_e')
        ex_glc.lower_bound = -100
        ex_glc.upper_bound = 1000
    except KeyError:
        for rxn in model.exchanges:
            if 'glc' in rxn.id.lower() or 'glucose' in rxn.id.lower():
                rxn.lower_bound = -100
                rxn.upper_bound = 1000
                break
    
    # 필수 영양소
    essentials = ['EX_nh4_e', 'EX_h2o_e', 'EX_h_e', 'EX_pi_e', 'EX_so4_e',
                  'EX_k_e', 'EX_o2_e', 'EX_co2_e']
    
    for ex_id in essentials:
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            if ex_id == 'EX_co2_e':
                ex_rxn.lower_bound = -1000
                ex_rxn.upper_bound = 1000
            else:
                ex_rxn.lower_bound = -1000
        except KeyError:
            pass
    
    print("\n블록된 반응 찾기 중...")
    blocked = find_blocked_reactions(model)
    
    blocked_ids = [rxn.id if hasattr(rxn, 'id') else rxn for rxn in blocked]
    
    print(f"\n총 블록된 반응: {len(blocked_ids)}개")
    
    # 핵심 반응의 블록 상태 확인
    print("\n핵심 반응 블록 상태:")
    key_reactions = {
        'Glycolysis': ['HEX1', 'PGI', 'PFK', 'FBA', 'GAPD', 'PGK', 'PYK'],
        'TCA': ['CS', 'ACONT', 'ICDHx', 'AKGDH', 'SUCOAS', 'SUCD', 'FUM', 'MDH'],
        'Glyoxylate': ['ICL', 'MALS'],
        'Biomass': ['Growth'],
        'Exchange': ['EX_glc__D_e', 'EX_ac_e']
    }
    
    blocked_analysis = []
    
    for category, rxn_ids in key_reactions.items():
        print(f"\n[{category}]")
        for rxn_id in rxn_ids:
            try:
                rxn = model.reactions.get_by_id(rxn_id)
                is_blocked = rxn_id in blocked_ids
                status = "[BLOCKED]" if is_blocked else "[OK]"
                blocked_analysis.append({
                    'Category': category,
                    'Reaction_ID': rxn_id,
                    'Is_Blocked': is_blocked,
                    'Status': status
                })
                print(f"  {status} {rxn_id}")
            except KeyError:
                print(f"  [MISSING] {rxn_id}")
                blocked_analysis.append({
                    'Category': category,
                    'Reaction_ID': rxn_id,
                    'Is_Blocked': False,
                    'Status': 'MISSING'
                })
    
    return blocked_analysis, blocked_ids

def test_growth_with_unlimited(model, biomass_rxn):
    """무제한 영양소로 성장 테스트"""
    print("\n" + "="*70)
    print("성장 테스트 (무제한 영양소)")
    print("="*70)
    
    # 모든 exchange 허용
    for rxn in model.exchanges:
        rxn.lower_bound = -1000
        rxn.upper_bound = 1000
    
    model.objective = biomass_rxn.id
    
    solution = model.optimize()
    
    print(f"\n최적화 상태: {solution.status}")
    
    if solution.status == 'optimal':
        biomass_flux = solution.fluxes.get(biomass_rxn.id, 0)
        print(f"Biomass flux: {biomass_flux:.6f} 1/h")
        
        if biomass_flux > 1e-6:
            print("\n[SUCCESS] 무제한 영양소로 성장 가능!")
            print("  → 모델 자체는 정상 작동합니다")
            print("  → 문제는 영양소 제약 조건에 있습니다")
            return True
        else:
            print("\n[FAIL] 무제한 영양소로도 성장 불가")
            print("  → 모델 구조 자체에 문제가 있을 수 있습니다")
            return False
    else:
        print(f"\n[ERROR] 최적화 실패: {solution.status}")
        print("  → 모델 구조에 심각한 문제가 있을 수 있습니다")
        return False

def main():
    print("="*70)
    print("모델 완전 진단 - 3단계 통합")
    print("="*70)
    
    # 모델 로드
    model = load_model("BaseModel.xml")
    
    # Biomass 반응 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass 반응: {biomass_rxn.id}")
    
    # 1단계: Biomass 반응 분석
    component_analysis = analyze_biomass_reaction(model, biomass_rxn)
    
    # 2단계: 기본 대사 경로 연결성
    pathway_status = check_basic_pathways(model)
    
    # 경로 연결성 확인
    metabolite_status = check_pathway_connectivity(model)
    
    # 무제한 영양소로 성장 테스트 (먼저 확인)
    can_grow_unlimited = test_growth_with_unlimited(model, biomass_rxn)
    
    # 3단계: 블록된 반응 분석 (성장 가능한 경우에만)
    blocked_analysis = []
    blocked_ids = []
    if can_grow_unlimited:
        try:
            blocked_analysis, blocked_ids = check_blocked_reactions_analysis(model)
        except Exception as e:
            print(f"\n[WARNING] 블록된 반응 분석 실패: {e}")
    else:
        print("\n[SKIP] 블록된 반응 분석 건너뜀 (모델이 infeasible 상태)")
    
    # 결과 저장
    if component_analysis:
        df_components = pd.DataFrame(component_analysis)
        df_components.to_csv('biomass_component_analysis.csv', index=False)
        print(f"\n[OK] Biomass 구성 요소 분석 저장: biomass_component_analysis.csv")
    
    if blocked_analysis:
        df_blocked = pd.DataFrame(blocked_analysis)
        df_blocked.to_csv('blocked_reactions_analysis.csv', index=False)
        print(f"[OK] 블록된 반응 분석 저장: blocked_reactions_analysis.csv")
    
    if metabolite_status:
        df_metabolites = pd.DataFrame(metabolite_status)
        df_metabolites.to_csv('metabolite_status.csv', index=False)
        print(f"[OK] 대사물질 상태 저장: metabolite_status.csv")
    
    # 최종 요약
    print("\n" + "="*70)
    print("최종 진단 결과")
    print("="*70)
    
    print(f"\n[1단계] Biomass 반응: {len(component_analysis)}개 구성 요소 분석 완료")
    print(f"[2단계] 기본 경로:")
    print(f"  - Glycolysis: {len(pathway_status['glycolysis']['found'])}개 반응 존재")
    print(f"  - TCA Cycle: {len(pathway_status['tca']['found'])}개 반응 존재")
    print(f"  - Glyoxylate Shunt: {len(pathway_status['glyoxylate']['found'])}개 반응 존재")
    print(f"[3단계] 블록된 반응: {len(blocked_ids)}개")
    
    if can_grow_unlimited:
        print("\n[결론] 모델 자체는 정상 작동합니다")
        print("  → 문제는 영양소 제약 조건 또는 부트스트랩 문제입니다")
    else:
        print("\n[결론] 모델 구조 자체에 문제가 있을 수 있습니다")
        print("  → Biomass 반응 또는 기본 대사 경로를 점검해야 합니다")
    
    print("="*70)

if __name__ == "__main__":
    main()
