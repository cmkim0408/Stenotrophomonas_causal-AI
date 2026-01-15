#!/usr/bin/env python
"""
SUCD (Succinate Dehydrogenase) 반응 찾기
다른 이름이나 대체 반응이 있는지 확인
"""

import cobra

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    return model

def find_succd_alternatives(model):
    """SUCD 대체 반응 찾기"""
    print("="*70)
    print("SUCD (Succinate Dehydrogenase) 반응 찾기")
    print("="*70)
    
    # Succinate와 Fumarate 관련 반응 찾기
    try:
        succ_c = model.metabolites.get_by_id('succ_c')
        fum_c = model.metabolites.get_by_id('fum_c')
        
        print(f"\nSuccinate (succ_c) 관련 반응: {len(succ_c.reactions)}개")
        print(f"Fumarate (fum_c) 관련 반응: {len(fum_c.reactions)}개")
        
        # Succinate를 소비하고 Fumarate를 생성하는 반응
        succ_to_fum = []
        for rxn in succ_c.reactions:
            if succ_c in rxn.reactants and fum_c in rxn.products:
                succ_to_fum.append(rxn)
            elif rxn.reversibility and succ_c in rxn.reactants and fum_c in rxn.reactants:
                # 가역 반응이면 양쪽 모두 가능
                succ_to_fum.append(rxn)
        
        print(f"\nSuccinate → Fumarate 변환 반응: {len(succ_to_fum)}개")
        for rxn in succ_to_fum:
            genes = [g.id for g in rxn.genes]
            gpr = str(rxn.gene_reaction_rule) if hasattr(rxn, 'gene_reaction_rule') else 'N/A'
            print(f"\n  {rxn.id}: {rxn.name}")
            print(f"    반응식: {rxn.reaction}")
            print(f"    가역성: {rxn.reversibility}")
            print(f"    유전자: {len(genes)}개")
            if genes:
                print(f"      {', '.join(genes)}")
            print(f"    GPR: {gpr}")
        
        # FADH2 관련 반응 확인
        try:
            fadh2_c = model.metabolites.get_by_id('fadh2_c')
            print(f"\nFADH2 생성 반응 (Succinate 관련):")
            fadh2_producing = [r for r in fadh2_c.reactions 
                             if fadh2_c in r.products and succ_c in r.reactants]
            for rxn in fadh2_producing:
                genes = [g.id for g in rxn.genes]
                print(f"  {rxn.id}: {rxn.reaction}")
                print(f"    유전자: {len(genes)}개 - {', '.join(genes[:5])}")
        except KeyError:
            print("\n[WARN] fadh2_c metabolite 없음")
        
        # 표준 SUCD 반응 패턴 검색
        print(f"\nSuccinate Dehydrogenase 관련 반응 검색:")
        sdh_rxns = [r for r in model.reactions 
                   if 'succinate dehydrogenase' in r.name.lower() or 
                      'SDH' in r.id or 'SUCD' in r.id]
        
        if sdh_rxns:
            print(f"  발견: {len(sdh_rxns)}개")
            for rxn in sdh_rxns:
                print(f"    {rxn.id}: {rxn.name}")
                print(f"      {rxn.reaction}")
        else:
            print(f"  [NOT FOUND] SUCD 관련 반응 없음")
        
        # Complex II (전자 전달계 복합체 II) 확인
        print(f"\n전자 전달계 복합체 II (Complex II) 관련 반응:")
        complex2_rxns = [r for r in model.reactions 
                        if 'complex' in r.name.lower() and 'ii' in r.name.lower()]
        if not complex2_rxns:
            complex2_rxns = [r for r in model.reactions 
                           if 'complex 2' in r.name.lower()]
        
        if complex2_rxns:
            for rxn in complex2_rxns:
                print(f"    {rxn.id}: {rxn.name}")
                print(f"      {rxn.reaction}")
        else:
            print(f"  [NOT FOUND] Complex II 반응 없음")
        
        return succ_to_fum
        
    except KeyError as e:
        print(f"[ERROR] {e}")
        return []

def check_tca_completeness(model):
    """TCA cycle 완성도 확인"""
    print("\n" + "="*70)
    print("TCA Cycle 완성도 확인")
    print("="*70)
    
    tca_steps = {
        '1. CS': 'Acetyl-CoA + OAA → Citrate',
        '2. ACONT': 'Citrate ↔ Isocitrate',
        '3. ICDH': 'Isocitrate → α-KG + NADH',
        '4. AKGDH': 'α-KG → Succinyl-CoA + NADH',
        '5. SUCOAS': 'Succinyl-CoA → Succinate + ATP',
        '6. SUCD': 'Succinate → Fumarate + FADH2',
        '7. FUM': 'Fumarate ↔ Malate',
        '8. MDH': 'Malate → OAA + NADH'
    }
    
    found_steps = []
    missing_steps = []
    
    for step_id, description in tca_steps.items():
        rxn_id = step_id.split('. ')[1]
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            found_steps.append(step_id)
            print(f"  {step_id}: {description} - [OK] ({rxn.id})")
        except KeyError:
            missing_steps.append(step_id)
            print(f"  {step_id}: {description} - [MISSING]")
    
    print(f"\n요약:")
    print(f"  완성된 단계: {len(found_steps)}/8")
    print(f"  누락된 단계: {len(missing_steps)}/8")
    if missing_steps:
        print(f"  누락: {', '.join(missing_steps)}")
    
    if '6. SUCD' in missing_steps:
        print(f"\n  [중요] SUCD가 누락되어 있습니다.")
        print(f"    → Succinate에서 Fumarate로의 직접 변환 경로가 없음")
        print(f"    → TCA cycle이 불완전할 수 있음")
        print(f"    → 하지만 다른 경로로 우회 가능할 수도 있음")

def main():
    model = load_model("BaseModel.xml")
    
    succ_to_fum = find_succd_alternatives(model)
    check_tca_completeness(model)
    
    print("\n" + "="*70)
    print("결론")
    print("="*70)
    if succ_to_fum:
        print(f"  [OK] Succinate → Fumarate 변환 반응 발견: {len(succ_to_fum)}개")
        print(f"    → SUCD가 다른 이름으로 존재하거나 대체 경로가 있을 수 있음")
    else:
        print(f"  [WARNING] SUCD 반응이 명확히 존재하지 않음")
        print(f"    → 갭필링이 필요할 수 있음")
        print(f"    → 또는 전자 전달계 복합체 II로 통합되어 있을 수 있음")
    print("="*70)

if __name__ == "__main__":
    main()

