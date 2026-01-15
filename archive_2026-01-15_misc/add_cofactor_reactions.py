#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
조효소 생산 경로 복구: NAD, NADP, CoA

Step 4-2에서 발견된 막힌 조효소:
- nad_c (NAD)
- nadp_c (NADP)
- coa_c (CoA)
"""

import cobra
from pathlib import Path
import sys

def check_cofactor_reactions(model, ref_model):
    """조효소 관련 반응 확인"""
    print("\n" + "="*70)
    print("조효소 관련 반응 확인")
    print("="*70)
    
    # NAD/NADP 관련 반응
    nad_reactions = {
        'NADK': 'NAD kinase (NAD -> NADP)',
        'NADPPPS': 'NADP phosphatase (NADP -> NAD)',
        'NADS1': 'NAD synthase (ammonia)',
        'NADS2': 'NAD synthase (glutamine)',
        'THD2pp': 'NADP transhydrogenase',
        'NADTRHD': 'NAD transhydrogenase',
    }
    
    # CoA 관련 반응
    coa_reactions = {
        'PPAT': 'Pantothenate phosphoribosyltransferase',
        'DPCK': 'Dephospho-CoA kinase',
        'COASY': 'CoA synthase',
    }
    
    print("\n[NAD/NADP 관련 반응]")
    nad_found_new = []
    nad_found_ref = []
    for rxn_id, rxn_name in nad_reactions.items():
        in_new = rxn_id in model.reactions
        in_ref = rxn_id in ref_model.reactions
        nad_found_new.append(rxn_id) if in_new else None
        nad_found_ref.append(rxn_id) if in_ref else None
        status_new = "[OK]" if in_new else "[MISSING]"
        status_ref = "[OK]" if in_ref else "[NOT FOUND]"
        print(f"  {rxn_id:12s}: 신규={status_new:10s} 레퍼런스={status_ref}")
    
    print("\n[CoA 관련 반응]")
    coa_found_new = []
    coa_found_ref = []
    for rxn_id, rxn_name in coa_reactions.items():
        in_new = rxn_id in model.reactions
        in_ref = rxn_id in ref_model.reactions
        coa_found_new.append(rxn_id) if in_new else None
        coa_found_ref.append(rxn_id) if in_ref else None
        status_new = "[OK]" if in_new else "[MISSING]"
        status_ref = "[OK]" if in_ref else "[NOT FOUND]"
        print(f"  {rxn_id:12s}: 신규={status_new:10s} 레퍼런스={status_ref}")
    
    return nad_found_new, nad_found_ref, coa_found_new, coa_found_ref

def add_reaction_from_reference(model, ref_model, rxn_id):
    """레퍼런스 모델에서 반응을 신규 모델에 추가"""
    if rxn_id in model.reactions:
        return False, "already_exists"
    
    if rxn_id not in ref_model.reactions:
        return False, "not_in_reference"
    
    try:
        ref_rxn = ref_model.reactions.get_by_id(rxn_id)
        
        # 새 반응 생성
        new_rxn = cobra.Reaction(rxn_id)
        new_rxn.name = ref_rxn.name if hasattr(ref_rxn, 'name') and ref_rxn.name else rxn_id
        new_rxn.lower_bound = ref_rxn.lower_bound
        new_rxn.upper_bound = ref_rxn.upper_bound
        
        # 대사물질 추가
        metabolites_dict = {}
        for met, coeff in ref_rxn.metabolites.items():
            met_id = met.id
            
            if met_id in model.metabolites:
                metabolites_dict[model.metabolites.get_by_id(met_id)] = coeff
            else:
                # 대사물질이 없으면 레퍼런스에서 복사
                ref_met = ref_model.metabolites.get_by_id(met_id)
                new_met = cobra.Metabolite(
                    id=ref_met.id,
                    formula=getattr(ref_met, 'formula', None),
                    name=getattr(ref_met, 'name', met_id),
                    compartment=ref_met.compartment
                )
                model.add_metabolites([new_met])
                metabolites_dict[new_met] = coeff
        
        new_rxn.add_metabolites(metabolites_dict)
        model.add_reactions([new_rxn])
        
        return True, "added"
    except Exception as e:
        return False, str(e)

def test_metabolite_production(model, metabolite_id):
    """특정 대사물질의 최대 생산량 테스트"""
    if metabolite_id not in model.metabolites:
        return None, "metabolite_not_in_model"
    
    with model:
        demand_id = f'DM_{metabolite_id}'
        if demand_id in model.reactions:
            model.remove_reactions([demand_id])
        
        demand_rxn = cobra.Reaction(demand_id)
        demand_rxn.add_metabolites({model.metabolites.get_by_id(metabolite_id): -1})
        model.add_reactions([demand_rxn])
        
        model.objective = demand_id
        model.objective_direction = 'max'
        
        solution = model.optimize()
        
        if solution.status == 'optimal':
            return solution.objective_value, "optimal"
        else:
            return 0.0, solution.status

def main():
    base_path = Path(__file__).parent.parent
    new_model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_with_BCAA.xml"
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    output_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_with_BCAA_cofactors.xml"
    
    print("="*70)
    print("조효소 생산 경로 복구: NAD, NADP, CoA")
    print("="*70)
    
    # 모델 로드
    print(f"\n[모델 로드]")
    new_model = cobra.io.read_sbml_model(str(new_model_path))
    ref_model = cobra.io.read_sbml_model(str(ref_model_path))
    print(f"[OK] 모델 로드 완료")
    
    # 조효소 반응 확인
    nad_found_new, nad_found_ref, coa_found_new, coa_found_ref = check_cofactor_reactions(new_model, ref_model)
    
    # 레퍼런스에는 있지만 신규 모델에는 없는 반응 찾기
    nad_to_add = [rxn_id for rxn_id in nad_found_ref if rxn_id not in nad_found_new]
    coa_to_add = [rxn_id for rxn_id in coa_found_ref if rxn_id not in coa_found_new]
    
    all_to_add = nad_to_add + coa_to_add
    
    if len(all_to_add) == 0:
        print("\n[OK] 모든 조효소 반응이 이미 신규 모델에 있습니다!")
    else:
        print(f"\n[추가할 반응] {len(all_to_add)}개")
        if nad_to_add:
            print(f"  NAD/NADP: {', '.join(nad_to_add)}")
        if coa_to_add:
            print(f"  CoA: {', '.join(coa_to_add)}")
        
        # 반응 추가
        print(f"\n[반응 추가 중...]")
        added_count = 0
        for rxn_id in all_to_add:
            success, msg = add_reaction_from_reference(new_model, ref_model, rxn_id)
            if success:
                print(f"[OK] {rxn_id} 추가")
                added_count += 1
            else:
                print(f"[WARNING] {rxn_id} 추가 실패: {msg}")
        
        print(f"\n[결과]")
        print(f"  추가된 반응: {added_count}/{len(all_to_add)}개")
        
        # 모델 저장
        cobra.io.write_sbml_model(new_model, str(output_path))
        print(f"\n[모델 저장] {output_path}")
    
    # 조효소 생산 테스트
    print(f"\n[조효소 생산 테스트]")
    cofactors_to_test = ['nad_c', 'nadp_c', 'coa_c']
    for met_id in cofactors_to_test:
        max_prod, status = test_metabolite_production(new_model, met_id)
        if max_prod is not None:
            status_str = "[OK]" if max_prod > 1e-6 else "[NO]"
            print(f"  {met_id:12s}: {status_str} max_production = {max_prod:.6f}")
        else:
            print(f"  {met_id:12s}: [ERROR] {status}")
    
    return new_model

if __name__ == "__main__":
    model = main()
