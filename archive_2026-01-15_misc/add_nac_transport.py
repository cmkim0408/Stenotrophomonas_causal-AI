#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nac transport 반응 추가

문제: nac_c가 생성되지 않음
원인: T_nac_e_to_nac_c transport 반응이 없음
해결: 레퍼런스 모델에서 T_nac_e_to_nac_c 반응 추가
"""

import cobra
from pathlib import Path
import sys

def load_model(model_path):
    try:
        model = cobra.io.read_sbml_model(str(model_path))
        print(f"[OK] 모델 로드 완료 (반응 수: {len(model.reactions)})")
        return model
    except Exception as e:
        print(f"[ERROR] 모델 로드 실패: {e}")
        sys.exit(1)

def check_nac_metabolites_and_reactions(model):
    """nac 관련 metabolite와 반응 확인"""
    print("\n" + "="*70)
    print("nac 관련 metabolite 및 반응 확인")
    print("="*70)
    
    # Metabolites 확인
    nac_mets = [m.id for m in model.metabolites if 'nac' in m.id.lower()]
    print(f"\n[nac 관련 metabolite]")
    for met_id in nac_mets:
        met = model.metabolites.get_by_id(met_id)
        print(f"  {met_id}: {met.name} (compartment: {met.compartment})")
    
    # Reactions 확인
    nac_rxns = [r.id for r in model.reactions if 'nac' in r.id.lower()]
    print(f"\n[nac 관련 반응]")
    for rxn_id in nac_rxns:
        rxn = model.reactions.get_by_id(rxn_id)
        print(f"  {rxn_id}: {rxn.reaction}")
        print(f"    bounds: [{rxn.lower_bound}, {rxn.upper_bound}]")
    
    return nac_mets, nac_rxns

def get_reference_nac_transport(ref_model):
    """레퍼런스 모델에서 nac transport 반응 정보 가져오기"""
    print("\n" + "="*70)
    print("레퍼런스 모델에서 nac transport 반응 확인")
    print("="*70)
    
    # T_nac_e_to_nac_c 반응 찾기
    nac_transport_ids = [
        'T_nac_e_to_nac_c',
        'R_T_nac_e_to_nac_c',
        'T_nac_e_to_c',
    ]
    
    for rxn_id in nac_transport_ids:
        if rxn_id in ref_model.reactions:
            rxn = ref_model.reactions.get_by_id(rxn_id)
            print(f"\n[찾은 반응] {rxn_id}")
            print(f"  반응식: {rxn.reaction}")
            print(f"  bounds: [{rxn.lower_bound}, {rxn.upper_bound}]")
            return rxn_id, rxn
    else:
        # 직접 찾기
        for rxn in ref_model.reactions:
            if 'nac' in rxn.id.lower() and 'T_' in rxn.id:
                mets = [m.id for m in rxn.metabolites]
                if 'nac_e' in mets and 'nac_c' in mets:
                    print(f"\n[찾은 반응] {rxn.id}")
                    print(f"  반응식: {rxn.reaction}")
                    print(f"  bounds: [{rxn.lower_bound}, {rxn.upper_bound}]")
                    return rxn.id, rxn
    
    return None, None

def add_nac_transport(model, ref_rxn):
    """nac transport 반응 추가"""
    print("\n" + "="*70)
    print("nac transport 반응 추가")
    print("="*70)
    
    # 필요한 metabolite 확인
    nac_e = None
    nac_c = None
    
    if 'nac_e' in model.metabolites:
        nac_e = model.metabolites.get_by_id('nac_e')
        print(f"[OK] nac_e 존재: {nac_e.id}")
    else:
        print("[ERROR] nac_e가 모델에 없습니다.")
        return False
    
    if 'nac_c' in model.metabolites:
        nac_c = model.metabolites.get_by_id('nac_c')
        print(f"[OK] nac_c 존재: {nac_c.id}")
    else:
        print("[ERROR] nac_c가 모델에 없습니다.")
        return False
    
    # Transport 반응 ID
    transport_id = 'T_nac_e_to_nac_c'
    
    if transport_id in model.reactions:
        print(f"[INFO] {transport_id} 반응이 이미 모델에 있습니다.")
        rxn = model.reactions.get_by_id(transport_id)
        print(f"  반응식: {rxn.reaction}")
        return True
    
    # Transport 반응 생성
    transport_rxn = cobra.Reaction(transport_id)
    transport_rxn.name = "Nicotinate transport (e<->c)"
    
    # 반응식: nac_e <=> nac_c
    transport_rxn.add_metabolites({
        nac_e: -1,
        nac_c: 1,
    })
    
    transport_rxn.lower_bound = -1000.0
    transport_rxn.upper_bound = 1000.0
    
    model.add_reactions([transport_rxn])
    
    print(f"[OK] {transport_id} 반응 추가 완료")
    print(f"  반응식: {transport_rxn.reaction}")
    print(f"  bounds: [{transport_rxn.lower_bound}, {transport_rxn.upper_bound}]")
    
    return True

def test_nac_c_production(model):
    """nac_c 생산 가능성 테스트"""
    print("\n" + "="*70)
    print("nac_c 생산 가능성 테스트")
    print("="*70)
    
    if 'nac_c' not in model.metabolites:
        print("[ERROR] nac_c가 모델에 없습니다.")
        return 0.0
    
    with model:
        # EX_nac_e 열기
        if 'EX_nac_e' in model.reactions:
            ex_nac = model.reactions.get_by_id('EX_nac_e')
            ex_nac.lower_bound = -1000.0
            ex_nac.upper_bound = 1000.0
        
        # Demand 반응 추가
        demand_id = 'DM_nac_c'
        if demand_id in model.reactions:
            model.remove_reactions([demand_id])
        
        demand_rxn = cobra.Reaction(demand_id)
        demand_rxn.add_metabolites({model.metabolites.get_by_id('nac_c'): -1})
        model.add_reactions([demand_rxn])
        
        model.objective = demand_id
        model.objective_direction = 'max'
        
        solution = model.optimize()
        
        max_prod = solution.objective_value if solution.status == 'optimal' else 0.0
        print(f"\n[결과]")
        print(f"  상태: {solution.status}")
        print(f"  nac_c 최대 생산량: {max_prod:.6f}")
        
        if max_prod > 1e-6:
            print(f"  [OK] nac_c 생성 가능!")
        else:
            print(f"  [NO] nac_c 생성 불가")
        
        return max_prod

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_with_BCAA_cofactors_ions_nad.xml"
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    output_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_with_BCAA_cofactors_ions_nad_transport.xml"
    
    print("="*70)
    print("nac transport 반응 추가")
    print("="*70)
    
    # 모델 로드
    model = load_model(model_path)
    ref_model = load_model(ref_model_path)
    
    # nac 관련 metabolite 및 반응 확인
    nac_mets, nac_rxns = check_nac_metabolites_and_reactions(model)
    
    # 레퍼런스 모델에서 transport 반응 확인
    ref_rxn_id, ref_rxn = get_reference_nac_transport(ref_model)
    
    if ref_rxn is None:
        print("[ERROR] 레퍼런스 모델에서 nac transport 반응을 찾을 수 없습니다.")
        return
    
    # Transport 반응 추가
    success = add_nac_transport(model, ref_rxn)
    
    if not success:
        print("[ERROR] nac transport 반응 추가 실패")
        return
    
    # nac_c 생산 가능성 테스트
    max_prod = test_nac_c_production(model)
    
    # 모델 저장
    cobra.io.write_sbml_model(model, str(output_path))
    print(f"\n[모델 저장] {output_path}")
    
    print("\n" + "="*70)
    print("완료")
    print("="*70)

if __name__ == "__main__":
    main()
