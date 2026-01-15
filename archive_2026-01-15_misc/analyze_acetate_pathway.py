#!/usr/bin/env python
"""
아세트산(acetate) 기반 대사 경로 분석
Acetate → Acetyl-CoA → TCA cycle/Glyoxylate shunt → Malate → Pyruvate/PEP
"""

import cobra
from cobra import Model
import pandas as pd

def load_model(model_path="BaseModel.xml"):
    """모델 로드"""
    print(f"모델 로드 중: {model_path}")
    model = cobra.io.read_sbml_model(model_path)
    print(f"[OK] 모델 로드 완료: {model.id}\n")
    return model

def find_reaction_by_id(model, reaction_id):
    """반응 ID로 반응 찾기"""
    try:
        return model.reactions.get_by_id(reaction_id)
    except KeyError:
        return None

def find_reactions_by_pattern(model, pattern):
    """패턴으로 반응 찾기"""
    results = []
    for rxn in model.reactions:
        if pattern.lower() in rxn.id.lower() or pattern.lower() in rxn.name.lower():
            results.append(rxn)
    return results

def print_reaction_info(reaction):
    """반응 정보 출력"""
    if reaction is None:
        print("  [NOT FOUND]")
        return
    
    print(f"  반응 ID: {reaction.id}")
    print(f"  이름: {reaction.name}")
    print(f"  반응식: {reaction.reaction}")
    print(f"  유전자: {[gene.id for gene in reaction.genes]}")
    print(f"  가역성: {reaction.reversibility}")
    print(f"  하한: {reaction.lower_bound}, 상한: {reaction.upper_bound}")
    print()

def analyze_acetate_to_acetylcoa(model):
    """Acetate → Acetyl-CoA 경로 분석"""
    print("="*70)
    print("1. Acetate → Acetyl-CoA 경로")
    print("="*70)
    
    # AcsA 반응 (ACS 또는 R_ACS)
    acs_rxn = find_reaction_by_id(model, 'ACS')
    if not acs_rxn:
        acs_rxn = find_reaction_by_id(model, 'R_ACS')
    if acs_rxn:
        print("[OK] ACS (Acetyl-CoA synthetase) 반응 발견:")
        print_reaction_info(acs_rxn)
    else:
        print("[ERROR] ACS 반응을 찾을 수 없습니다.")
    
    # Acetate exchange reaction
    ac_ex = find_reaction_by_id(model, 'EX_ac_e')
    if ac_ex:
        print("Acetate 교환 반응:")
        print_reaction_info(ac_ex)
    else:
        # 다른 이름으로 찾아보기
        ac_ex = find_reaction_by_id(model, 'EX_ac_c')
        if ac_ex:
            print("Acetate 관련 반응:")
            print_reaction_info(ac_ex)

def analyze_tca_cycle(model):
    """TCA Cycle 경로 분석"""
    print("="*70)
    print("2. TCA Cycle 경로")
    print("="*70)
    
    tca_reactions = {
        'CS': 'Citrate synthase',
        'ACONT': 'Aconitase',
        'ICDH': 'Isocitrate dehydrogenase',
        'AKGDH': 'Alpha-ketoglutarate dehydrogenase',
        'SUCOAS': 'Succinate-CoA ligase',
        'SUCD': 'Succinate dehydrogenase',
        'FUM': 'Fumarase',
        'MDH': 'Malate dehydrogenase'
    }
    
    found_reactions = {}
    for rxn_id, rxn_name in tca_reactions.items():
        rxn = find_reaction_by_id(model, rxn_id)
        if rxn:
            found_reactions[rxn_id] = rxn
    
    # 패턴 검색으로 추가 찾기
    for pattern in ['citrate', 'aconit', 'isocitrate', 'akg', 'succ', 'fumar', 'malate']:
        rxns = find_reactions_by_pattern(model, pattern)
        for rxn in rxns:
            if rxn.id not in found_reactions:
                # TCA 관련인지 확인
                if any(keyword in rxn.name.lower() for keyword in ['citrate', 'aconit', 'isocitrate', 
                                                                     'ketoglutarate', 'succinate', 
                                                                     'fumarate', 'malate', 'tca']):
                    found_reactions[rxn.id] = rxn
    
    print(f"발견된 TCA Cycle 관련 반응: {len(found_reactions)}개\n")
    
    # 주요 TCA 반응 출력
    key_reactions = ['CS', 'ICDH', 'AKGDH', 'SUCD', 'FUM', 'MDH']
    for key in key_reactions:
        if key in found_reactions:
            print(f"[{key}] {found_reactions[key].name}:")
            print_reaction_info(found_reactions[key])
    
    # Acetyl-CoA를 사용하는 반응 찾기
    accoa = model.metabolites.get_by_id('accoa_c')
    print("\nAcetyl-CoA를 사용하는 주요 반응:")
    accoa_reactions = [rxn for rxn in accoa.reactions if 'CS' in rxn.id or 'citrate' in rxn.name.lower()]
    for rxn in accoa_reactions[:5]:  # 상위 5개만
        print_reaction_info(rxn)

def analyze_glyoxylate_shunt(model):
    """Glyoxylate Shunt 경로 분석"""
    print("="*70)
    print("3. Glyoxylate Shunt 경로")
    print("="*70)
    
    # 주요 glyoxylate shunt 반응
    glyoxylate_reactions = {
        'ICL': 'Isocitrate lyase',
        'MALS': 'Malate synthase',
        'ACONT': 'Aconitase (공통)'
    }
    
    found_reactions = {}
    
    # Isocitrate lyase 찾기
    icl_rxns = find_reactions_by_pattern(model, 'isocitrate lyase')
    if not icl_rxns:
        icl_rxns = find_reactions_by_pattern(model, 'ICL')
    
    # Malate synthase 찾기
    mals_rxns = find_reactions_by_pattern(model, 'malate synthase')
    if not mals_rxns:
        mals_rxns = find_reactions_by_pattern(model, 'MALS')
    
    # Glyoxylate 관련 반응 찾기
    glyoxylate_rxns = find_reactions_by_pattern(model, 'glyoxylate')
    
    all_glyoxylate = list(set(icl_rxns + mals_rxns + glyoxylate_rxns))
    
    if all_glyoxylate:
        print(f"발견된 Glyoxylate Shunt 관련 반응: {len(all_glyoxylate)}개\n")
        for rxn in all_glyoxylate:
            print_reaction_info(rxn)
    else:
        print("[WARNING] Glyoxylate shunt 반응을 명확히 찾을 수 없습니다.")
        print("대안: Isocitrate와 Glyoxylate 관련 반응 검색 중...\n")
        
        # Isocitrate 관련 반응
        icit_rxns = find_reactions_by_pattern(model, 'isocitrate')
        for rxn in icit_rxns:
            if 'lyase' in rxn.name.lower() or 'ICL' in rxn.id:
                print_reaction_info(rxn)

def analyze_malate_to_pyr_pep(model):
    """Malate → Pyruvate/PEP 경로 분석"""
    print("="*70)
    print("4. Malate → Pyruvate/PEP 경로")
    print("="*70)
    
    # Malate 관련 반응 찾기
    malate = model.metabolites.get_by_id('mal__L_c')
    
    print("Malate에서 생성되는 반응들:")
    malate_producing = [rxn for rxn in malate.reactions if malate in rxn.products]
    
    # Pyruvate 관련 반응
    pyr_rxns = find_reactions_by_pattern(model, 'pyruvate')
    
    # PEP 관련 반응
    pep_rxns = find_reactions_by_pattern(model, 'phosphoenolpyruvate')
    if not pep_rxns:
        pep_rxns = find_reactions_by_pattern(model, 'PEP')
    
    # Malate → Pyruvate (Malic enzyme)
    malic_enzyme = find_reactions_by_pattern(model, 'malic enzyme')
    if not malic_enzyme:
        malic_enzyme = find_reactions_by_pattern(model, 'ME')
    
    # Malate → OAA → PEP (PEP carboxykinase)
    pepck_rxns = find_reactions_by_pattern(model, 'PEPCK')
    if not pepck_rxns:
        pepck_rxns = find_reactions_by_pattern(model, 'phosphoenolpyruvate carboxykinase')
    
    print("\n[1] Malic Enzyme (Malate → Pyruvate):")
    if malic_enzyme:
        for rxn in malic_enzyme:
            print_reaction_info(rxn)
    else:
        print("  [NOT FOUND] Malic enzyme 반응을 찾을 수 없습니다.\n")
    
    print("\n[2] PEP Carboxykinase (OAA → PEP):")
    if pepck_rxns:
        for rxn in pepck_rxns:
            print_reaction_info(rxn)
    else:
        print("  [NOT FOUND] PEPCK 반응을 찾을 수 없습니다.\n")
    
    # Malate dehydrogenase (Malate ↔ OAA)
    mdh_rxns = find_reactions_by_pattern(model, 'malate dehydrogenase')
    if not mdh_rxns:
        mdh_rxns = [rxn for rxn in model.reactions if rxn.id == 'MDH']
    
    print("\n[3] Malate Dehydrogenase (Malate ↔ OAA):")
    if mdh_rxns:
        for rxn in mdh_rxns:
            print_reaction_info(rxn)
    else:
        print("  [NOT FOUND] Malate dehydrogenase 반응을 찾을 수 없습니다.\n")

def find_pathway_connections(model):
    """경로 간 연결 확인"""
    print("="*70)
    print("5. 경로 연결 확인")
    print("="*70)
    
    try:
        # 주요 대사물질 확인
        metabolites_to_check = {
            'ac_c': 'Acetate',
            'accoa_c': 'Acetyl-CoA',
            'cit_c': 'Citrate',
            'icit_c': 'Isocitrate',
            'akg_c': 'Alpha-ketoglutarate',
            'succ_c': 'Succinate',
            'fum_c': 'Fumarate',
            'mal__L_c': 'Malate',
            'oaa_c': 'Oxaloacetate',
            'pyr_c': 'Pyruvate',
            'pep_c': 'PEP'
        }
        
        print("\n주요 대사물질 존재 확인:")
        for met_id, met_name in metabolites_to_check.items():
            try:
                met = model.metabolites.get_by_id(met_id)
                print(f"  [OK] {met_name} ({met_id}): {len(met.reactions)}개 반응에 참여")
            except KeyError:
                print(f"  [NOT FOUND] {met_name} ({met_id})")
        
        print()
        
    except Exception as e:
        print(f"[ERROR] 대사물질 확인 중 오류: {e}\n")

def main():
    """메인 함수"""
    # 모델 로드
    model = load_model("BaseModel.xml")
    
    # 1. Acetate → Acetyl-CoA
    analyze_acetate_to_acetylcoa(model)
    
    # 2. TCA Cycle
    analyze_tca_cycle(model)
    
    # 3. Glyoxylate Shunt
    analyze_glyoxylate_shunt(model)
    
    # 4. Malate → Pyruvate/PEP
    analyze_malate_to_pyr_pep(model)
    
    # 5. 경로 연결 확인
    find_pathway_connections(model)
    
    print("\n" + "="*70)
    print("경로 분석 완료!")
    print("="*70)

if __name__ == "__main__":
    main()

