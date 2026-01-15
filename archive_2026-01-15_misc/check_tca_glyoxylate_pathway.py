#!/usr/bin/env python
"""
TCA Cycle과 Glyoxylate Shunt 경로 점검 및 연결성 테스트
- 모든 반응 존재 여부 확인
- 경로 연결성 확인
- 대사물질 플럭스 분석
"""

import cobra
from cobra.flux_analysis import find_blocked_reactions

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    return model

def check_reaction_exists(model, rxn_id):
    """반응 존재 여부 확인"""
    try:
        rxn = model.reactions.get_by_id(rxn_id)
        return True, rxn
    except KeyError:
        return False, None

def check_succd_alternatives(model):
    """SUCD 또는 대체 반응 확인"""
    print("="*70)
    print("SUCD (Succinate Dehydrogenase) 확인")
    print("="*70)
    
    # 표준 SUCD 확인
    exists, rxn = check_reaction_exists(model, 'SUCD')
    if exists:
        print(f"\n[OK] SUCD 반응 존재: {rxn.id}")
        print(f"  반응식: {rxn.reaction}")
        print(f"  유전자: {len(rxn.genes)}개")
        return True, rxn
    
    # R_SUCD4 확인
    exists, rxn = check_reaction_exists(model, 'R_SUCD4')
    if exists:
        print(f"\n[OK] R_SUCD4 반응 발견 (SUCD 대체 가능): {rxn.id}")
        print(f"  이름: {rxn.name}")
        print(f"  반응식: {rxn.reaction}")
        print(f"  유전자: {len(rxn.genes)}개")
        if rxn.genes:
            print(f"  GPR: {rxn.gene_reaction_rule}")
        return True, rxn
    
    # Succinate → Fumarate 변환 반응 찾기
    try:
        succ_c = model.metabolites.get_by_id('succ_c')
        fum_c = model.metabolites.get_by_id('fum_c')
        
        succ_to_fum = []
        for rxn in succ_c.reactions:
            if succ_c in rxn.reactants and fum_c in rxn.products:
                succ_to_fum.append(rxn)
            elif rxn.reversibility and fum_c in rxn.metabolites:
                # 가역 반응인 경우
                if (succ_c in rxn.reactants and fum_c in rxn.reactants) or \
                   (succ_c in rxn.products and fum_c in rxn.products):
                    succ_to_fum.append(rxn)
        
        if succ_to_fum:
            print(f"\n[OK] Succinate → Fumarate 변환 반응 발견: {len(succ_to_fum)}개")
            for rxn in succ_to_fum:
                print(f"  {rxn.id}: {rxn.name}")
                print(f"    반응식: {rxn.reaction}")
        return len(succ_to_fum) > 0, succ_to_fum[0] if succ_to_fum else None
        
    except KeyError:
        print("\n[ERROR] succ_c 또는 fum_c metabolite를 찾을 수 없음")
        return False, None

def check_tca_reactions(model):
    """TCA cycle 모든 반응 확인"""
    print("\n" + "="*70)
    print("TCA Cycle 반응 확인")
    print("="*70)
    
    tca_reactions = {
        'CS': 'Citrate Synthase',
        'ACONT': 'Aconitase',
        'ICDHx': 'Isocitrate Dehydrogenase (NAD+)',
        'ICDHyr': 'Isocitrate Dehydrogenase (NADP+)',
        'AKGDH': 'α-Ketoglutarate Dehydrogenase',
        'SUCOAS': 'Succinyl-CoA Synthetase',
        'SUCD': 'Succinate Dehydrogenase',
        'FUM': 'Fumarase',
        'MDH': 'Malate Dehydrogenase'
    }
    
    found = {}
    missing = []
    
    for rxn_id, rxn_name in tca_reactions.items():
        exists, rxn = check_reaction_exists(model, rxn_id)
        if exists:
            found[rxn_id] = rxn
            genes = [g.id for g in rxn.genes]
            print(f"\n[{rxn_id}] {rxn_name}: [OK]")
            print(f"  반응식: {rxn.reaction}")
            print(f"  가역성: {'Yes' if rxn.reversibility else 'No'}")
            print(f"  유전자: {len(genes)}개")
            if genes:
                print(f"    {', '.join(genes[:3])}" + (f" ... 외 {len(genes)-3}개" if len(genes) > 3 else ""))
        else:
            missing.append(rxn_id)
            print(f"\n[{rxn_id}] {rxn_name}: [MISSING]")
    
    # SUCD 대체 확인
    if 'SUCD' in missing:
        print("\n[SUCD 대체 반응 확인 중...]")
        succd_exists, succd_rxn = check_succd_alternatives(model)
        if succd_exists:
            found['SUCD_alt'] = succd_rxn
            missing.remove('SUCD')
    
    print(f"\n요약: {len(found)}/{len(tca_reactions)} 반응 존재")
    if missing:
        print(f"누락: {', '.join(missing)}")
    
    return found, missing

def check_glyoxylate_reactions(model):
    """Glyoxylate shunt 반응 확인"""
    print("\n" + "="*70)
    print("Glyoxylate Shunt 반응 확인")
    print("="*70)
    
    gly_reactions = {
        'ICL': 'Isocitrate Lyase',
        'MALS': 'Malate Synthase'
    }
    
    found = {}
    missing = []
    
    for rxn_id, rxn_name in gly_reactions.items():
        exists, rxn = check_reaction_exists(model, rxn_id)
        if exists:
            found[rxn_id] = rxn
            genes = [g.id for g in rxn.genes]
            print(f"\n[{rxn_id}] {rxn_name}: [OK]")
            print(f"  반응식: {rxn.reaction}")
            print(f"  가역성: {'Yes' if rxn.reversibility else 'No'}")
            print(f"  유전자: {len(genes)}개")
            if genes:
                print(f"    {', '.join(genes)}")
        else:
            missing.append(rxn_id)
            print(f"\n[{rxn_id}] {rxn_name}: [MISSING]")
    
    print(f"\n요약: {len(found)}/{len(gly_reactions)} 반응 존재")
    if missing:
        print(f"누락: {', '.join(missing)}")
    
    return found, missing

def test_pathway_connectivity(model):
    """TCA/Glyoxylate 경로 연결성 테스트"""
    print("\n" + "="*70)
    print("경로 연결성 테스트")
    print("="*70)
    
    # 주요 대사물질 확인
    metabolites = {
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
    
    print("\n주요 대사물질 확인:")
    found_mets = {}
    for met_id, met_name in metabolites.items():
        try:
            met = model.metabolites.get_by_id(met_id)
            found_mets[met_id] = met
            reactions = [r.id for r in met.reactions]
            print(f"  [OK] {met_name} ({met_id}): {len(reactions)}개 반응 참여")
        except KeyError:
            print(f"  [MISSING] {met_name} ({met_id})")
    
    # TCA cycle 경로: Acetyl-CoA → Citrate → ... → OAA
    print("\n[TCA Cycle 경로 연결 확인]")
    tca_path = ['accoa_c', 'cit_c', 'icit_c', 'akg_c', 'succoa_c', 
                'succ_c', 'fum_c', 'mal__L_c', 'oaa_c']
    
    all_connected = True
    for i in range(len(tca_path) - 1):
        met1_id = tca_path[i]
        met2_id = tca_path[i+1]
        if met1_id in found_mets and met2_id in found_mets:
            met1 = found_mets[met1_id]
            met2 = found_mets[met2_id]
            # 공통 반응 찾기
            common_rxns = set(met1.reactions) & set(met2.reactions)
            if common_rxns:
                print(f"  [{met1_id}] → [{met2_id}]: 연결됨 ({len(common_rxns)}개 반응)")
            else:
                print(f"  [{met1_id}] → [{met2_id}]: [WARNING] 직접 연결 반응 없음")
                all_connected = False
        else:
            print(f"  [{met1_id}] → [{met2_id}]: [ERROR] 대사물질 없음")
            all_connected = False
    
    # Glyoxylate shunt 경로: Isocitrate → Glyoxylate → Malate
    print("\n[Glyoxylate Shunt 경로 연결 확인]")
    gly_path = ['icit_c', 'glx_c', 'mal__L_c']
    
    for i in range(len(gly_path) - 1):
        met1_id = gly_path[i]
        met2_id = gly_path[i+1]
        if met1_id in found_mets and met2_id in found_mets:
            met1 = found_mets[met1_id]
            met2 = found_mets[met2_id]
            common_rxns = set(met1.reactions) & set(met2.reactions)
            if common_rxns:
                print(f"  [{met1_id}] → [{met2_id}]: 연결됨 ({len(common_rxns)}개 반응)")
            else:
                print(f"  [{met1_id}] → [{met2_id}]: [WARNING] 직접 연결 반응 없음")
                all_connected = False
        else:
            print(f"  [{met1_id}] → [{met2_id}]: [ERROR] 대사물질 없음")
            all_connected = False
    
    return all_connected

def check_blocked_reactions(model):
    """블록된 반응 확인"""
    print("\n" + "="*70)
    print("블록된 반응 확인 (TCA/Glyoxylate 관련)")
    print("="*70)
    
    tca_gly_rxns = ['CS', 'ACONT', 'ICDHx', 'ICDHyr', 'AKGDH', 'SUCOAS', 
                    'SUCD', 'FUM', 'MDH', 'ICL', 'MALS']
    
    blocked = find_blocked_reactions(model)
    
    # blocked는 반응 객체 또는 반응 ID 리스트일 수 있음
    blocked_ids = [rxn.id if hasattr(rxn, 'id') else rxn for rxn in blocked]
    blocked_tca_gly_ids = [rxn_id for rxn_id in blocked_ids if rxn_id in tca_gly_rxns]
    
    if blocked_tca_gly_ids:
        print(f"\n[WARNING] 블록된 TCA/Glyoxylate 반응: {len(blocked_tca_gly_ids)}개")
        for rxn_id in blocked_tca_gly_ids:
            try:
                rxn = model.reactions.get_by_id(rxn_id)
                print(f"  {rxn.id}: {rxn.name}")
            except KeyError:
                print(f"  {rxn_id}")
    else:
        print(f"\n[OK] TCA/Glyoxylate 반응 블록 없음")
    
    return blocked_tca_gly_ids

def main():
    model = load_model("BaseModel.xml")
    
    print("="*70)
    print("TCA Cycle과 Glyoxylate Shunt 종합 점검")
    print("="*70)
    
    # 1. TCA cycle 반응 확인
    tca_found, tca_missing = check_tca_reactions(model)
    
    # 2. Glyoxylate shunt 반응 확인
    gly_found, gly_missing = check_glyoxylate_reactions(model)
    
    # 3. 경로 연결성 확인
    is_connected = test_pathway_connectivity(model)
    
    # 4. 블록된 반응 확인
    blocked = check_blocked_reactions(model)
    
    # 최종 요약
    print("\n" + "="*70)
    print("최종 요약")
    print("="*70)
    print(f"\nTCA Cycle: {len(tca_found)}/9 반응 존재")
    if tca_missing:
        print(f"  누락: {', '.join(tca_missing)}")
    
    print(f"\nGlyoxylate Shunt: {len(gly_found)}/2 반응 존재")
    if gly_missing:
        print(f"  누락: {', '.join(gly_missing)}")
    
    print(f"\n경로 연결성: {'OK' if is_connected else 'WARNING'}")
    print(f"블록된 반응: {len(blocked)}개")
    
    if not tca_missing and not gly_missing and is_connected and len(blocked) == 0:
        print("\n[결론] TCA Cycle과 Glyoxylate Shunt가 완전히 구성되어 있습니다.")
    else:
        print("\n[결론] 일부 문제가 발견되었습니다.")
        if tca_missing:
            print(f"  - TCA cycle 누락 반응: {', '.join(tca_missing)}")
        if gly_missing:
            print(f"  - Glyoxylate shunt 누락 반응: {', '.join(gly_missing)}")
        if not is_connected:
            print(f"  - 경로 연결성 문제 있음")
        if len(blocked) > 0:
            print(f"  - 블록된 반응: {len(blocked)}개")
    
    print("="*70)

if __name__ == "__main__":
    main()
