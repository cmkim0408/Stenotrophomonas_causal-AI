#!/usr/bin/env python
"""
TCA cycle 및 Glyoxylate shunt의 모든 반응 분석
각 반응의 존재 여부, GPR, 방향성 확인
"""

import cobra

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    return model

def get_reaction_details(model, rxn_id):
    """반응 상세 정보 가져오기"""
    try:
        rxn = model.reactions.get_by_id(rxn_id)
        return {
            'id': rxn.id,
            'name': rxn.name,
            'reaction': rxn.reaction,
            'reversible': rxn.reversibility,
            'lower_bound': rxn.lower_bound,
            'upper_bound': rxn.upper_bound,
            'genes': [g.id for g in rxn.genes],
            'gpr': str(rxn.gene_reaction_rule) if hasattr(rxn, 'gene_reaction_rule') else 'N/A'
        }
    except KeyError:
        return None

def analyze_tca_cycle(model):
    """TCA cycle 모든 반응 분석"""
    print("="*70)
    print("TCA Cycle (Tricarboxylic Acid Cycle) - 모든 반응")
    print("="*70)
    
    # TCA cycle 표준 반응들
    tca_reactions = {
        'CS': {
            'name': 'Citrate Synthase',
            'reaction': 'Acetyl-CoA + Oxaloacetate + H2O → Citrate + CoA + H+',
            'description': 'TCA cycle 시작점'
        },
        'ACONT': {
            'name': 'Aconitase',
            'reaction': 'Citrate ↔ Isocitrate',
            'description': 'Citrate를 Isocitrate로 변환'
        },
        'ICDHx': {
            'name': 'Isocitrate Dehydrogenase (NAD+)',
            'reaction': 'Isocitrate + NAD+ ↔ α-Ketoglutarate + CO2 + NADH',
            'description': 'NAD+ 의존성 Isocitrate 탈수소화'
        },
        'ICDHyr': {
            'name': 'Isocitrate Dehydrogenase (NADP+)',
            'reaction': 'Isocitrate + NADP+ ↔ α-Ketoglutarate + CO2 + NADPH',
            'description': 'NADP+ 의존성 Isocitrate 탈수소화'
        },
        'AKGDH': {
            'name': 'α-Ketoglutarate Dehydrogenase Complex',
            'reaction': 'α-Ketoglutarate + CoA + NAD+ → Succinyl-CoA + CO2 + NADH',
            'description': 'α-KG를 Succinyl-CoA로 변환'
        },
        'SUCOAS': {
            'name': 'Succinyl-CoA Synthetase',
            'reaction': 'Succinyl-CoA + ADP + Pi ↔ Succinate + ATP + CoA',
            'description': 'Succinyl-CoA에서 GTP/ATP 생성'
        },
        'SUCD': {
            'name': 'Succinate Dehydrogenase',
            'reaction': 'Succinate + FAD ↔ Fumarate + FADH2',
            'description': 'Succinate를 Fumarate로 변환, FADH2 생성'
        },
        'FUM': {
            'name': 'Fumarase',
            'reaction': 'Fumarate + H2O ↔ Malate',
            'description': 'Fumarate를 Malate로 변환'
        },
        'MDH': {
            'name': 'Malate Dehydrogenase',
            'reaction': 'Malate + NAD+ ↔ Oxaloacetate + NADH + H+',
            'description': 'Malate를 Oxaloacetate로 변환, NADH 생성'
        }
    }
    
    print(f"\n{'반응 ID':<15} {'이름':<35} {'상태':<10} {'가역성':<8} {'GPR':<30}")
    print("-" * 110)
    
    found_reactions = []
    missing_reactions = []
    
    for rxn_id, info in tca_reactions.items():
        details = get_reaction_details(model, rxn_id)
        
        if details:
            found_reactions.append(rxn_id)
            gpr_str = details['gpr'][:28] + '...' if len(details['gpr']) > 28 else details['gpr']
            reversible_str = 'Yes' if details['reversible'] else 'No'
            print(f"{rxn_id:<15} {details['name']:<35} {'[OK]':<10} {reversible_str:<8} {gpr_str:<30}")
            print(f"               반응식: {details['reaction']}")
            print(f"               유전자: {len(details['genes'])}개 - {', '.join(details['genes'][:5])}")
            if len(details['genes']) > 5:
                print(f"                        ... 외 {len(details['genes'])-5}개")
        else:
            missing_reactions.append(rxn_id)
            print(f"{rxn_id:<15} {info['name']:<35} {'[MISSING]':<10} {'N/A':<8} {'N/A':<30}")
            print(f"               설명: {info['description']}")
        
        print()
    
    print(f"\n요약:")
    print(f"  존재하는 반응: {len(found_reactions)}/{len(tca_reactions)}")
    print(f"  존재: {', '.join(found_reactions)}")
    if missing_reactions:
        print(f"  누락: {', '.join(missing_reactions)}")
    
    return found_reactions, missing_reactions

def analyze_glyoxylate_shunt(model):
    """Glyoxylate shunt 모든 반응 분석"""
    print("\n" + "="*70)
    print("Glyoxylate Shunt (Glyoxylate Cycle) - 모든 반응")
    print("="*70)
    
    glyoxylate_reactions = {
        'ICL': {
            'name': 'Isocitrate Lyase',
            'reaction': 'Isocitrate → Glyoxylate + Succinate',
            'description': 'Isocitrate를 Glyoxylate와 Succinate로 분해'
        },
        'MALS': {
            'name': 'Malate Synthase',
            'reaction': 'Glyoxylate + Acetyl-CoA + H2O → Malate + CoA + H+',
            'description': 'Glyoxylate와 Acetyl-CoA로부터 Malate 합성'
        }
    }
    
    print(f"\n{'반응 ID':<15} {'이름':<35} {'상태':<10} {'가역성':<8} {'GPR':<30}")
    print("-" * 110)
    
    found_reactions = []
    missing_reactions = []
    
    for rxn_id, info in glyoxylate_reactions.items():
        details = get_reaction_details(model, rxn_id)
        
        if details:
            found_reactions.append(rxn_id)
            gpr_str = details['gpr'][:28] + '...' if len(details['gpr']) > 28 else details['gpr']
            reversible_str = 'Yes' if details['reversible'] else 'No'
            print(f"{rxn_id:<15} {details['name']:<35} {'[OK]':<10} {reversible_str:<8} {gpr_str:<30}")
            print(f"               반응식: {details['reaction']}")
            print(f"               유전자: {len(details['genes'])}개 - {', '.join(details['genes'][:5])}")
            if len(details['genes']) > 5:
                print(f"                        ... 외 {len(details['genes'])-5}개")
        else:
            missing_reactions.append(rxn_id)
            print(f"{rxn_id:<15} {info['name']:<35} {'[MISSING]':<10} {'N/A':<8} {'N/A':<30}")
            print(f"               설명: {info['description']}")
        
        print()
    
    print(f"\n요약:")
    print(f"  존재하는 반응: {len(found_reactions)}/{len(glyoxylate_reactions)}")
    print(f"  존재: {', '.join(found_reactions)}")
    if missing_reactions:
        print(f"  누락: {', '.join(missing_reactions)}")
    
    # Glyoxylate shunt의 생물학적 의미 설명
    print(f"\n[Glyoxylate Shunt의 중요성]")
    print(f"  - Acetate와 같은 2탄소 기질로부터 성장할 때 필수")
    print(f"  - TCA cycle만으로는 탄소 손실 (CO2 방출)")
    print(f"  - Glyoxylate shunt를 통해 탄소를 보존하면서 OAA 재생성")
    print(f"  - 순환: Acetyl-CoA → Isocitrate → Glyoxylate → Malate → OAA")
    
    return found_reactions, missing_reactions

def check_alternative_reactions(model):
    """대체 반응 또는 변형 반응 확인"""
    print("\n" + "="*70)
    print("TCA/Glyoxylate 관련 대체 반응 또는 변형 확인")
    print("="*70)
    
    # 비표준 이름으로 검색
    search_terms = ['citrate', 'aconit', 'isocitrate', 'ketoglutarate', 
                   'succinyl', 'succinate', 'fumarate', 'malate', 
                   'oxaloacetate', 'glyoxylate', 'lyase', 'synthase']
    
    alternative_rxns = []
    for rxn in model.reactions:
        rxn_lower = (rxn.id + ' ' + rxn.name).lower()
        if any(term in rxn_lower for term in search_terms):
            # TCA cycle의 핵심 대사물질 포함 여부 확인
            metabolites = [m.id.lower() for m in rxn.metabolites]
            tca_metabolites = ['cit_c', 'icit_c', 'akg_c', 'succoa_c', 'succ_c', 
                             'fum_c', 'mal__l_c', 'oaa_c', 'glx_c']
            if any(m in metabolites for m in tca_metabolites):
                if rxn.id not in ['CS', 'ACONT', 'ICDHx', 'ICDHyr', 'AKGDH', 
                                 'SUCOAS', 'SUCD', 'FUM', 'MDH', 'ICL', 'MALS']:
                    alternative_rxns.append(rxn)
    
    print(f"\nTCA/Glyoxylate 관련 추가 반응: {len(alternative_rxns)}개")
    for rxn in alternative_rxns[:20]:  # 처음 20개만 표시
        genes = [g.id for g in rxn.genes]
        gpr = str(rxn.gene_reaction_rule) if hasattr(rxn, 'gene_reaction_rule') else 'N/A'
        print(f"\n  {rxn.id}: {rxn.name}")
        print(f"    반응식: {rxn.reaction}")
        print(f"    유전자: {len(genes)}개 - {', '.join(genes[:3])}")
        if len(genes) > 3:
            print(f"            ... 외 {len(genes)-3}개")
        print(f"    GPR: {gpr[:60]}")

def summarize_gap_filling_status(model):
    """갭필링 상태 요약"""
    print("\n" + "="*70)
    print("갭필링 상태 요약")
    print("="*70)
    
    print("\n[지금까지 추가한 반응]")
    print("  1. R_ACS (기존 반응) - Smlt4623 유전자 추가 (GPR 수정)")
    print("  2. EX_ac_e (Exchange 반응) - 기본 구조 반응")
    print("  3. ACt (Transport 반응) - 기본 구조 반응")
    print("\n  → 갭필링 없음, 유전자가 있는 반응만 수정/추가")
    
    tca_found, tca_missing = analyze_tca_cycle(model)
    gly_found, gly_missing = analyze_glyoxylate_shunt(model)
    
    print("\n[전체 TCA/Glyoxylate 반응 상태]")
    total_tca = len(tca_found) + len(tca_missing)
    total_gly = len(gly_found) + len(gly_missing)
    print(f"  TCA cycle: {len(tca_found)}/{total_tca} 존재")
    print(f"  Glyoxylate shunt: {len(gly_found)}/{total_gly} 존재")
    
    if tca_missing or gly_missing:
        print(f"\n  [주의] 누락된 반응이 있습니다:")
        if tca_missing:
            print(f"    TCA cycle: {', '.join(tca_missing)}")
        if gly_missing:
            print(f"    Glyoxylate shunt: {', '.join(gly_missing)}")
        print(f"  → 갭필링이 필요할 수 있습니다")
    else:
        print(f"\n  [OK] 모든 필수 반응이 존재합니다")
        print(f"  → 갭필링 불필요")

def main():
    model = load_model("BaseModel.xml")
    
    print("="*70)
    print("TCA Cycle 및 Glyoxylate Shunt 완전 분석")
    print("="*70)
    
    tca_found, tca_missing = analyze_tca_cycle(model)
    gly_found, gly_missing = analyze_glyoxylate_shunt(model)
    check_alternative_reactions(model)
    summarize_gap_filling_status(model)
    
    print("\n" + "="*70)
    print("분석 완료")
    print("="*70)

if __name__ == "__main__":
    main()

