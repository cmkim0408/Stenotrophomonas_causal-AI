#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 4-3: 막힌 전구체 주변을 우선 복구

우선순위 1: BCAA (Leucine, Valine) 합성 경로
- 레퍼런스 모델에서 BCAA 관련 반응 확인
- 신규 모델에 없는 반응 추가
"""

import cobra
from pathlib import Path
import sys

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
    }
    
    for ex_id, (lb, ub) in essential_exchanges.items():
        if ex_id in model.reactions:
            ex_rxn = model.reactions.get_by_id(ex_id)
            ex_rxn.lower_bound = lb
            ex_rxn.upper_bound = ub
    
    return model

def check_bcaa_reactions_in_model(model):
    """모델에서 BCAA 관련 반응 확인"""
    print("\n" + "="*70)
    print("BCAA 관련 반응 확인 (신규 모델)")
    print("="*70)
    
    bcaa_reactions = {
        'ALS': 'Acetolactate synthase',
        'KARI': 'Ketol-acid reductoisomerase',
        'DHAD': 'Dihydroxyacid dehydratase',
        'IPMS': 'Isopropylmalate synthase',
        'IPMI': 'Isopropylmalate isomerase',
        'IPMDH': 'Isopropylmalate dehydrogenase',
        'IPMD': '3-Isopropylmalate dehydrogenase',
        'BCAT_VAL': 'Branched-chain aminotransferase (Val)',
        'BCAT_LEU': 'Branched-chain aminotransferase (Leu)',
    }
    
    found = []
    missing = []
    
    for rxn_id, rxn_name in bcaa_reactions.items():
        if rxn_id in model.reactions:
            rxn = model.reactions.get_by_id(rxn_id)
            found.append(rxn_id)
            print(f"[OK] {rxn_id:12s}: {rxn_name}")
            print(f"      반응식: {rxn.reaction}")
        else:
            missing.append(rxn_id)
            print(f"[MISSING] {rxn_id:12s}: {rxn_name}")
    
    print(f"\n[요약]")
    print(f"  발견: {len(found)}개")
    print(f"  누락: {len(missing)}개")
    
    return found, missing

def check_bcaa_reactions_reference(ref_model):
    """레퍼런스 모델에서 BCAA 관련 반응 확인"""
    print("\n" + "="*70)
    print("BCAA 관련 반응 확인 (레퍼런스 모델)")
    print("="*70)
    
    bcaa_reactions = {
        'ALS': 'Acetolactate synthase',
        'KARI': 'Ketol-acid reductoisomerase',
        'DHAD': 'Dihydroxyacid dehydratase',
        'IPMS': 'Isopropylmalate synthase',
        'IPMI': 'Isopropylmalate isomerase',
        'IPMDH': 'Isopropylmalate dehydrogenase',
        'IPMD': '3-Isopropylmalate dehydrogenase',
        'BCAT_VAL': 'Branched-chain aminotransferase (Val)',
        'BCAT_LEU': 'Branched-chain aminotransferase (Leu)',
    }
    
    found = []
    
    for rxn_id, rxn_name in bcaa_reactions.items():
        if rxn_id in ref_model.reactions:
            rxn = ref_model.reactions.get_by_id(rxn_id)
            found.append(rxn_id)
            print(f"[OK] {rxn_id:12s}: {rxn_name}")
            print(f"      반응식: {rxn.reaction}")
            print(f"      bounds: [{rxn.lower_bound}, {rxn.upper_bound}]")
        else:
            print(f"[NOT FOUND] {rxn_id:12s}: {rxn_name}")
    
    print(f"\n[요약]")
    print(f"  발견: {len(found)}개")
    
    return found

def get_reaction_from_reference(ref_model, rxn_id):
    """레퍼런스 모델에서 반응 정보 가져오기"""
    if rxn_id not in ref_model.reactions:
        return None
    
    ref_rxn = ref_model.reactions.get_by_id(rxn_id)
    
    reaction_info = {
        'id': rxn_id,
        'name': ref_rxn.name if hasattr(ref_rxn, 'name') else rxn_id,
        'reaction': ref_rxn.reaction,
        'lower_bound': ref_rxn.lower_bound,
        'upper_bound': ref_rxn.upper_bound,
        'metabolites': dict(ref_rxn.metabolites)
    }
    
    return reaction_info

def add_reaction_from_reference(model, ref_model, rxn_id):
    """레퍼런스 모델에서 반응을 신규 모델에 추가"""
    if rxn_id in model.reactions:
        return False, "already_exists"
    
    reaction_info = get_reaction_from_reference(ref_model, rxn_id)
    if reaction_info is None:
        return False, "not_in_reference"
    
    try:
        # 새 반응 생성
        new_rxn = cobra.Reaction(rxn_id)
        new_rxn.name = reaction_info['name']
        new_rxn.lower_bound = reaction_info['lower_bound']
        new_rxn.upper_bound = reaction_info['upper_bound']
        
        # 대사물질 추가
        metabolites_dict = {}
        for met_id, coeff in reaction_info['metabolites'].items():
            # 대사물질이 모델에 있는지 확인
            if met_id in model.metabolites:
                metabolites_dict[model.metabolites.get_by_id(met_id)] = coeff
            else:
                # 대사물질이 없으면 레퍼런스에서 복사
                ref_met = ref_model.metabolites.get_by_id(met_id)
                new_met = cobra.Metabolite(
                    id=ref_met.id,
                    formula=ref_met.formula if hasattr(ref_met, 'formula') else None,
                    name=ref_met.name if hasattr(ref_met, 'name') else met_id,
                    compartment=ref_met.compartment
                )
                model.add_metabolites([new_met])
                metabolites_dict[new_met] = coeff
        
        new_rxn.add_metabolites(metabolites_dict)
        model.add_reactions([new_rxn])
        
        return True, "added"
    except Exception as e:
        return False, str(e)

def main():
    base_path = Path(__file__).parent.parent
    new_model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_with_ACtexi.xml"
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    
    print("="*70)
    print("Step 4-3: 막힌 전구체 주변 복구 - BCAA 합성 경로")
    print("="*70)
    
    # 모델 로드
    print(f"\n[모델 로드]")
    print(f"  신규 모델: {new_model_path}")
    print(f"  레퍼런스 모델: {ref_model_path}")
    
    if not new_model_path.exists():
        print(f"[ERROR] 신규 모델 파일이 없습니다: {new_model_path}")
        return
    
    if not ref_model_path.exists():
        print(f"[ERROR] 레퍼런스 모델 파일이 없습니다: {ref_model_path}")
        return
    
    new_model = cobra.io.read_sbml_model(str(new_model_path))
    ref_model = cobra.io.read_sbml_model(str(ref_model_path))
    
    print(f"[OK] 모델 로드 완료")
    print(f"  신규 모델: {len(new_model.reactions)}개 반응")
    print(f"  레퍼런스 모델: {len(ref_model.reactions)}개 반응")
    
    # 배지 조건 설정
    new_model = setup_media_forced(new_model)
    
    # 신규 모델에서 BCAA 반응 확인
    new_found, new_missing = check_bcaa_reactions_in_model(new_model)
    
    # 레퍼런스 모델에서 BCAA 반응 확인
    ref_found = check_bcaa_reactions_reference(ref_model)
    
    # 레퍼런스에는 있지만 신규 모델에는 없는 반응 찾기
    print("\n" + "="*70)
    print("레퍼런스에만 있는 BCAA 반응 (추가 필요)")
    print("="*70)
    
    to_add = []
    for rxn_id in ref_found:
        if rxn_id not in new_found:
            to_add.append(rxn_id)
            reaction_info = get_reaction_from_reference(ref_model, rxn_id)
            print(f"\n[{rxn_id}]")
            print(f"  이름: {reaction_info['name']}")
            print(f"  반응식: {reaction_info['reaction']}")
            print(f"  bounds: [{reaction_info['lower_bound']}, {reaction_info['upper_bound']}]")
    
    if len(to_add) == 0:
        print("\n[OK] 모든 BCAA 반응이 이미 신규 모델에 있습니다!")
    else:
        print(f"\n[추가 필요] {len(to_add)}개 반응")
        print(f"  {', '.join(to_add)}")
    
    return new_model, ref_model, to_add

if __name__ == "__main__":
    new_model, ref_model, to_add = main()
