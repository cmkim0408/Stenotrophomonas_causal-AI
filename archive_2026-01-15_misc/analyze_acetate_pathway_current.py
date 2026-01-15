#!/usr/bin/env python
"""
현재 신규 모델에서 아세트산(acetate) 전환 경로 분석
"""

import cobra
from pathlib import Path

def load_model(model_path):
    print(f"모델 로드 중: {model_path}")
    model = cobra.io.read_sbml_model(model_path)
    print(f"[OK] 모델 로드 완료: {model.id}")
    return model

def find_acetate_reactions(model):
    """아세트산 관련 반응 찾기"""
    acetate_reactions = []
    
    # acetate 메타볼라이트 찾기
    acetate_mets = []
    for met in model.metabolites:
        if 'ac_c' in met.id or 'acetate' in met.id.lower():
            acetate_mets.append(met)
    
    print(f"\n[Acetate 메타볼라이트]")
    for met in acetate_mets:
        print(f"  {met.id} ({met.name if met.name else '이름 없음'})")
    
    # acetate를 포함하는 반응 찾기
    for rxn in model.reactions:
        reactants = [met.id for met in rxn.reactants]
        products = [met.id for met in rxn.products]
        all_mets = reactants + products
        
        if 'ac_c' in all_mets or any('acetate' in met.lower() for met in all_mets):
            acetate_reactions.append({
                'id': rxn.id,
                'name': rxn.name if rxn.name else '',
                'equation': rxn.reaction,
                'reactants': reactants,
                'products': products,
                'reversible': rxn.reversibility,
                'lower_bound': rxn.lower_bound,
                'upper_bound': rxn.upper_bound,
                'subsystem': rxn.subsystem if hasattr(rxn, 'subsystem') else ''
            })
    
    return acetate_reactions, acetate_mets

def analyze_acetate_to_accoa(model):
    """Acetate -> Acetyl-CoA 전환 경로 분석"""
    print("\n" + "="*80)
    print("Acetate -> Acetyl-CoA 전환 경로 분석")
    print("="*80)
    
    # ac_c -> accoa_c 직접 전환 반응 찾기
    direct_conversion = []
    for rxn in model.reactions:
        reactants = [met.id for met in rxn.reactants]
        products = [met.id for met in rxn.products]
        
        if 'ac_c' in reactants and 'accoa_c' in products:
            direct_conversion.append({
                'id': rxn.id,
                'name': rxn.name if rxn.name else '',
                'equation': rxn.reaction,
                'subsystem': rxn.subsystem if hasattr(rxn, 'subsystem') else ''
            })
    
    if direct_conversion:
        print(f"\n[직접 전환 반응] {len(direct_conversion)}개 발견")
        for rxn_info in direct_conversion:
            print(f"\n  {rxn_info['id']}")
            if rxn_info['name']:
                print(f"    이름: {rxn_info['name']}")
            print(f"    반응식: {rxn_info['equation']}")
            if rxn_info['subsystem']:
                print(f"    경로: {rxn_info['subsystem']}")
    else:
        print("\n[직접 전환 반응] 없음")
    
    # Acetyl-CoA 관련 반응도 찾기
    accoa_reactions = []
    for rxn in model.reactions:
        reactants = [met.id for met in rxn.reactants]
        products = [met.id for met in rxn.products]
        
        if 'accoa_c' in reactants or 'accoa_c' in products:
            accoa_reactions.append({
                'id': rxn.id,
                'name': rxn.name if rxn.name else '',
                'equation': rxn.reaction,
                'subsystem': rxn.subsystem if hasattr(rxn, 'subsystem') else ''
            })
    
    print(f"\n[Acetyl-CoA 관련 반응] 총 {len(accoa_reactions)}개")
    print("  주요 반응:")
    for rxn_info in accoa_reactions[:10]:  # 처음 10개만
        print(f"    {rxn_info['id']}: {rxn_info['equation'][:80]}")
        if rxn_info['subsystem']:
            print(f"      경로: {rxn_info['subsystem']}")
    
    return direct_conversion, accoa_reactions

def find_alternative_pathways(model):
    """대체 경로 찾기 (예: Acetyl-phosphate 경유)"""
    print("\n" + "="*80)
    print("대체 경로 분석 (Acetyl-phosphate 경유 등)")
    print("="*80)
    
    # Acetyl-phosphate 관련 반응
    actp_reactions = []
    for rxn in model.reactions:
        reactants = [met.id for met in rxn.reactants]
        products = [met.id for met in rxn.products]
        
        if 'actp_c' in reactants or 'actp_c' in products:
            actp_reactions.append({
                'id': rxn.id,
                'name': rxn.name if rxn.name else '',
                'equation': rxn.reaction
            })
    
    if actp_reactions:
        print(f"\n[Acetyl-phosphate 관련 반응] {len(actp_reactions)}개")
        for rxn_info in actp_reactions:
            print(f"\n  {rxn_info['id']}")
            if rxn_info['name']:
                print(f"    이름: {rxn_info['name']}")
            print(f"    반응식: {rxn_info['equation']}")
    else:
        print("\n[Acetyl-phosphate 관련 반응] 없음")

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    
    model = load_model(str(model_path))
    
    print("\n" + "="*80)
    print("신규 모델에서 아세트산 전환 경로 분석")
    print("="*80)
    
    # 1. Acetate 관련 반응 찾기
    acetate_rxns, acetate_mets = find_acetate_reactions(model)
    
    print(f"\n[Acetate 관련 반응] 총 {len(acetate_rxns)}개")
    print("\n주요 반응:")
    for i, rxn_info in enumerate(acetate_rxns[:20], 1):  # 처음 20개만
        print(f"\n{i}. {rxn_info['id']}")
        if rxn_info['name']:
            print(f"   이름: {rxn_info['name']}")
        print(f"   반응식: {rxn_info['equation']}")
        if rxn_info['subsystem']:
            print(f"   경로: {rxn_info['subsystem']}")
        print(f"   방향: {'Reversible' if rxn_info['reversible'] else 'Irreversible'}")
    
    if len(acetate_rxns) > 20:
        print(f"\n  ... 외 {len(acetate_rxns) - 20}개 반응")
    
    # 2. Acetate -> Acetyl-CoA 직접 전환 분석
    direct_conv, accoa_rxns = analyze_acetate_to_accoa(model)
    
    # 3. 대체 경로 찾기
    find_alternative_pathways(model)
    
    # 4. Exchange 반응 확인
    print("\n" + "="*80)
    print("Acetate Exchange 반응")
    print("="*80)
    exchange_rxns = [rxn_info for rxn_info in acetate_rxns if 'EX_' in rxn_info['id'] or 'exchange' in rxn_info['id'].lower()]
    if exchange_rxns:
        for rxn_info in exchange_rxns:
            print(f"\n  {rxn_info['id']}")
            print(f"    반응식: {rxn_info['equation']}")
    else:
        print("\n  Exchange 반응 없음")
    
    # 5. 결론
    print("\n" + "="*80)
    print("결론")
    print("="*80)
    
    if direct_conv:
        print(f"\n[OK] 신규 모델에는 Acetate -> Acetyl-CoA 직접 전환 반응이 {len(direct_conv)}개 있습니다:")
        for rxn_info in direct_conv:
            print(f"  - {rxn_info['id']}: {rxn_info['equation']}")
    else:
        print("\n[ERROR] 신규 모델에는 Acetate -> Acetyl-CoA 직접 전환 반응이 없습니다.")
        print("  -> ACS_ADP 반응이 누락되어 있습니다 (레퍼런스 모델에는 존재)")
        print("  -> 대체 경로가 있는지 확인 필요")

if __name__ == "__main__":
    main()
