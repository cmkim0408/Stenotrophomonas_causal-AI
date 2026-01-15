#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NAD/NADP 생산 경로 자세히 확인

문제: NADK, NADS1, NADS2 같은 반응이 있는데도 NAD/NADP가 생성되지 않음
원인 분석:
1. NAD/NADP 생산 경로의 실제 반응 확인
2. 필요한 전구체 확인
3. 반응식/대사물질 불일치 확인
"""

import cobra
from pathlib import Path
import pandas as pd

def check_nad_nadp_reactions(model, ref_model):
    """NAD/NADP 관련 반응 확인"""
    print("\n" + "="*70)
    print("NAD/NADP 관련 반응 확인")
    print("="*70)
    
    nad_reactions = {
        'NADK': 'NAD kinase (NAD -> NADP)',
        'NADPPPS': 'NADP phosphatase (NADP -> NAD)',
        'NADS1': 'NAD synthase (ammonia)',
        'NADS2': 'NAD synthase (glutamine)',
        'THD2pp': 'NADP transhydrogenase',
        'NADTRHD': 'NAD transhydrogenase',
        'NAPRT': 'Nicotinate phosphoribosyltransferase',
        'NAMPT': 'Nicotinamide phosphoribosyltransferase',
        'NMNAT': 'NMN adenylyltransferase',
        'NADS_LUMP': 'NAD synthetase lumped',
        'NMNAT_LUMP': 'NMN adenylyltransferase lumped',
    }
    
    print("\n[NAD/NADP 관련 반응]")
    found_new = {}
    found_ref = {}
    
    for rxn_id, rxn_name in nad_reactions.items():
        in_new = rxn_id in model.reactions
        in_ref = rxn_id in ref_model.reactions
        
        found_new[rxn_id] = in_new
        found_ref[rxn_id] = in_ref
        
        status_new = "[OK]" if in_new else "[MISSING]"
        status_ref = "[OK]" if in_ref else "[NOT FOUND]"
        print(f"  {rxn_id:12s}: 신규={status_new:10s} 레퍼런스={status_ref}")
        
        if in_new:
            rxn = model.reactions.get_by_id(rxn_id)
            print(f"             반응식: {rxn.reaction}")
            print(f"             bounds: [{rxn.lower_bound}, {rxn.upper_bound}]")
        elif in_ref:
            rxn = ref_model.reactions.get_by_id(rxn_id)
            print(f"             반응식: {rxn.reaction}")
            print(f"             bounds: [{rxn.lower_bound}, {rxn.upper_bound}]")
    
    return found_new, found_ref

def test_nad_pathway(model):
    """NAD 생산 경로 테스트"""
    print("\n" + "="*70)
    print("NAD 생산 경로 테스트")
    print("="*70)
    
    # NAD 전구체 확인
    nad_precursors = {
        'nicotinate_c': 'Nicotinate',
        'nicotinamide_c': 'Nicotinamide',
        'namn_c': 'NaMN (Nicotinate mononucleotide)',
        'nmn_c': 'NMN (Nicotinamide mononucleotide)',
    }
    
    print("\n[NAD 전구체 생산 가능성]")
    for met_id, met_name in nad_precursors.items():
        if met_id in model.metabolites:
            with model:
                demand_id = f'DM_{met_id}'
                if demand_id in model.reactions:
                    model.remove_reactions([demand_id])
                
                demand_rxn = cobra.Reaction(demand_id)
                demand_rxn.add_metabolites({model.metabolites.get_by_id(met_id): -1})
                model.add_reactions([demand_rxn])
                
                model.objective = demand_id
                model.objective_direction = 'max'
                
                solution = model.optimize()
                max_prod = solution.objective_value if solution.status == 'optimal' else 0.0
                status_str = "[OK]" if max_prod > 1e-6 else "[NO]"
                print(f"  {met_id:20s}: {status_str} max_production = {max_prod:.6f}")
        else:
            print(f"  {met_id:20s}: [NOT IN MODEL]")

def check_nad_metabolites(model, ref_model):
    """NAD 관련 대사물질 확인"""
    print("\n" + "="*70)
    print("NAD 관련 대사물질 확인")
    print("="*70)
    
    nad_metabolites = [
        'nad_c', 'nadp_c', 'nadh_c', 'nadph_c',
        'nicotinate_c', 'nicotinamide_c',
        'namn_c', 'nmn_c',
        'nac_c', 'ncam_c',  # 외부 형태
    ]
    
    print("\n[대사물질 존재 여부]")
    for met_id in nad_metabolites:
        in_new = met_id in model.metabolites
        in_ref = met_id in ref_model.metabolites
        status_new = "[OK]" if in_new else "[MISSING]"
        status_ref = "[OK]" if in_ref else "[NOT FOUND]"
        print(f"  {met_id:20s}: 신규={status_new:10s} 레퍼런스={status_ref}")

def check_exchange_for_nad(model):
    """NAD 전구체 Exchange 확인"""
    print("\n" + "="*70)
    print("NAD 전구체 Exchange 확인")
    print("="*70)
    
    exchanges = ['EX_nac_e', 'EX_ncam_e', 'EX_nad_e', 'EX_nadp_e']
    
    print("\n[Exchange 반응]")
    for ex_id in exchanges:
        if ex_id in model.reactions:
            ex_rxn = model.reactions.get_by_id(ex_id)
            print(f"  {ex_id:15s}: bounds = [{ex_rxn.lower_bound}, {ex_rxn.upper_bound}]")
        else:
            print(f"  {ex_id:15s}: [NOT IN MODEL]")

def main():
    base_path = Path(__file__).parent.parent
    new_model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_with_BCAA_cofactors_ions.xml"
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    
    print("="*70)
    print("NAD/NADP 생산 경로 자세히 확인")
    print("="*70)
    
    # 모델 로드
    print(f"\n[모델 로드]")
    new_model = cobra.io.read_sbml_model(str(new_model_path))
    ref_model = cobra.io.read_sbml_model(str(ref_model_path))
    print(f"[OK] 모델 로드 완료")
    
    # NAD/NADP 반응 확인
    found_new, found_ref = check_nad_nadp_reactions(new_model, ref_model)
    
    # NAD 관련 대사물질 확인
    check_nad_metabolites(new_model, ref_model)
    
    # NAD 전구체 Exchange 확인
    check_exchange_for_nad(new_model)
    
    # NAD 생산 경로 테스트
    test_nad_pathway(new_model)
    
    # 레퍼런스에는 있지만 신규 모델에는 없는 반응 찾기
    print("\n" + "="*70)
    print("레퍼런스에만 있는 NAD/NADP 반응")
    print("="*70)
    
    missing = []
    for rxn_id, in_ref in found_ref.items():
        if in_ref and not found_new.get(rxn_id, False):
            missing.append(rxn_id)
            ref_rxn = ref_model.reactions.get_by_id(rxn_id)
            print(f"\n[{rxn_id}]")
            print(f"  이름: {ref_rxn.name if hasattr(ref_rxn, 'name') else rxn_id}")
            print(f"  반응식: {ref_rxn.reaction}")
            print(f"  bounds: [{ref_rxn.lower_bound}, {ref_rxn.upper_bound}]")
    
    if len(missing) == 0:
        print("\n[OK] 모든 NAD/NADP 반응이 이미 신규 모델에 있습니다!")
    else:
        print(f"\n[추가 필요] {len(missing)}개 반응")
        print(f"  {', '.join(missing)}")
    
    return new_model, ref_model, missing

if __name__ == "__main__":
    new_model, ref_model, missing = main()
