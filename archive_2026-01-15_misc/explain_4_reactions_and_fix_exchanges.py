#!/usr/bin/env python
"""
4개 반응 설명 및 Exchange bounds 수정
"""

import cobra
from pathlib import Path
import pandas as pd

def load_model(model_path):
    model = cobra.io.read_sbml_model(model_path)
    return model

def explain_4_reactions(ref_model):
    """4개 반응 설명"""
    print("="*80)
    print("4개 반응 설명")
    print("="*80)
    
    reactions = {
        'ACS_ADP': 'Acetyl-CoA synthetase (ADP-forming)',
        'SUCDi': 'Succinate dehydrogenase (irreversible)',
        'PEPCK_ATP': 'Phosphoenolpyruvate carboxykinase (ATP)',
        'ACtexi': 'Acetate transport (exchange)'
    }
    
    for rxn_id, description in reactions.items():
        if rxn_id in ref_model.reactions:
            rxn = ref_model.reactions.get_by_id(rxn_id)
            print(f"\n[{rxn_id}] {description}")
            print(f"  반응식: {rxn.reaction}")
            print(f"  bounds: [{rxn.lower_bound}, {rxn.upper_bound}]")
            
            # 메타볼라이트 확인
            print(f"  메타볼라이트:")
            for met, coeff in rxn.metabolites.items():
                if abs(coeff) > 0.1:  # 주요 메타볼라이트만
                    print(f"    {met.id}: {coeff}")
        else:
            print(f"\n[{rxn_id}] {description}")
            print(f"  [없음] 레퍼런스 모델에 없음")

def add_4_reactions(new_model, ref_model):
    """4개 반응 추가"""
    print("\n" + "="*80)
    print("4개 반응 추가")
    print("="*80)
    
    reactions_to_add = ['ACS_ADP', 'SUCDi', 'PEPCK_ATP', 'ACtexi']
    added = []
    
    for rxn_id in reactions_to_add:
        if rxn_id in [r.id for r in new_model.reactions]:
            print(f"\n  {rxn_id}: 이미 존재")
            continue
        
        if rxn_id not in ref_model.reactions:
            print(f"\n  {rxn_id}: 레퍼런스 모델에 없음")
            continue
        
        try:
            ref_rxn = ref_model.reactions.get_by_id(rxn_id)
            
            new_rxn = cobra.Reaction(rxn_id)
            new_rxn.name = ref_rxn.name
            new_rxn.lower_bound = ref_rxn.lower_bound
            new_rxn.upper_bound = ref_rxn.upper_bound
            
            metabolites_dict = {}
            for met, coeff in ref_rxn.metabolites.items():
                if met.id in [m.id for m in new_model.metabolites]:
                    metabolites_dict[new_model.metabolites.get_by_id(met.id)] = coeff
                else:
                    new_met = cobra.Metabolite(
                        met.id,
                        name=met.name,
                        compartment=met.compartment,
                        formula=met.formula if hasattr(met, 'formula') else None,
                        charge=met.charge if hasattr(met, 'charge') else None
                    )
                    new_model.add_metabolites([new_met])
                    metabolites_dict[new_met] = coeff
            
            new_rxn.add_metabolites(metabolites_dict)
            new_model.add_reactions([new_rxn])
            added.append(rxn_id)
            print(f"\n  {rxn_id}: 추가됨")
            print(f"    반응식: {new_rxn.reaction}")
        except Exception as e:
            print(f"\n  {rxn_id}: 추가 실패 - {e}")
    
    return added

def fix_exchange_bounds(new_model, ref_model):
    """Exchange bounds를 레퍼런스 모델과 동일하게 설정"""
    print("\n" + "="*80)
    print("Exchange bounds 수정 (레퍼런스 모델과 동일하게)")
    print("="*80)
    
    key_exchanges = [
        'EX_ac_e', 'EX_o2_e', 'EX_hco3_e', 'EX_nh4_e', 
        'EX_pi_e', 'EX_h2o_e', 'EX_h_e', 'EX_so4_e',
        'EX_k_e', 'EX_na1_e', 'EX_mg2_e', 'EX_ca2_e',
        'EX_fe2_e', 'EX_fe3_e', 'EX_co2_e'
    ]
    
    fixed = []
    
    for ex_id in key_exchanges:
        if ex_id in ref_model.reactions and ex_id in new_model.reactions:
            ref_rxn = ref_model.reactions.get_by_id(ex_id)
            new_rxn = new_model.reactions.get_by_id(ex_id)
            
            old_bounds = f"[{new_rxn.lower_bound}, {new_rxn.upper_bound}]"
            new_rxn.lower_bound = ref_rxn.lower_bound
            new_rxn.upper_bound = ref_rxn.upper_bound
            new_bounds = f"[{new_rxn.lower_bound}, {new_rxn.upper_bound}]"
            
            if old_bounds != new_bounds:
                fixed.append(ex_id)
                print(f"  {ex_id:20s}: {old_bounds} -> {new_bounds}")
        elif ex_id in ref_model.reactions and ex_id not in new_model.reactions:
            # 레퍼런스에만 있는 Exchange 추가
            ref_rxn = ref_model.reactions.get_by_id(ex_id)
            
            # Exchange 반응 생성
            new_ex = cobra.Reaction(ex_id)
            new_ex.name = ref_rxn.name
            new_ex.lower_bound = ref_rxn.lower_bound
            new_ex.upper_bound = ref_rxn.upper_bound
            
            # 메타볼라이트 확인
            for met, coeff in ref_rxn.metabolites.items():
                if met.id in [m.id for m in new_model.metabolites]:
                    new_ex.add_metabolites({new_model.metabolites.get_by_id(met.id): coeff})
                else:
                    # 메타볼라이트 추가
                    new_met = cobra.Metabolite(
                        met.id,
                        name=met.name,
                        compartment=met.compartment,
                        formula=met.formula if hasattr(met, 'formula') else None,
                        charge=met.charge if hasattr(met, 'charge') else None
                    )
                    new_model.add_metabolites([new_met])
                    new_ex.add_metabolites({new_met: coeff})
            
            new_model.add_reactions([new_ex])
            fixed.append(ex_id)
            print(f"  {ex_id:20s}: 추가됨 [{new_ex.lower_bound}, {new_ex.upper_bound}]")
    
    return fixed

def fix_atpm_bounds(new_model, ref_model):
    """ATPM bounds를 레퍼런스 모델과 동일하게 설정"""
    print("\n" + "="*80)
    print("ATPM bounds 수정 (레퍼런스 모델과 동일하게)")
    print("="*80)
    
    ref_atpm = ref_model.reactions.get_by_id('ATPM')
    new_atpm = new_model.reactions.get_by_id('ATPM')
    
    old_bounds = f"[{new_atpm.lower_bound}, {new_atpm.upper_bound}]"
    new_atpm.lower_bound = ref_atpm.lower_bound
    new_atpm.upper_bound = ref_atpm.upper_bound
    new_bounds = f"[{new_atpm.lower_bound}, {new_atpm.upper_bound}]"
    
    print(f"  ATPM: {old_bounds} -> {new_bounds}")

def test_after_fixes(new_model):
    """수정 후 테스트"""
    print("\n" + "="*80)
    print("수정 후 FBA 테스트")
    print("="*80)
    
    # ATPM=0 설정
    atpm_rxn = new_model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    new_model.objective = 'Growth'
    solution = new_model.optimize()
    
    print(f"\n[FBA 결과]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    
    key_reactions = ['ACS_ADP', 'SUCDi', 'PEPCK_ATP', 'ACtexi', 'CS', 'ICL', 'MALS', 'ATPS4rpp']
    print(f"\n[주요 반응 플럭스]")
    active_count = 0
    for rxn_id in key_reactions:
        if rxn_id in new_model.reactions:
            flux = solution.fluxes.get(rxn_id, 0.0)
            if abs(flux) > 1e-6:
                print(f"  {rxn_id:20s}: {flux:>12.6f}")
                active_count += 1
    
    if active_count > 0:
        print(f"\n[성공] {active_count}개 반응이 작동함!")
        return True
    else:
        print(f"\n[실패] 모든 반응 플럭스가 0")
        return False

def main():
    base_path = Path(__file__).parent.parent
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    new_model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    
    print("="*80)
    print("4개 반응 설명 및 Exchange bounds 수정")
    print("="*80)
    
    # 모델 로드
    ref_model = load_model(str(ref_model_path))
    new_model = load_model(str(new_model_path))
    
    # 4개 반응 설명
    explain_4_reactions(ref_model)
    
    # 4개 반응 추가
    added_reactions = add_4_reactions(new_model, ref_model)
    
    # Exchange bounds 수정
    fixed_exchanges = fix_exchange_bounds(new_model, ref_model)
    
    # ATPM bounds 수정
    fix_atpm_bounds(new_model, ref_model)
    
    # 수정 후 테스트
    success = test_after_fixes(new_model)
    
    print("\n" + "="*80)
    print("최종 결과")
    print("="*80)
    
    print(f"\n[추가된 반응] {len(added_reactions)}개: {', '.join(added_reactions)}")
    print(f"[수정된 Exchange] {len(fixed_exchanges)}개")
    print(f"[ATPM bounds 수정] 완료")
    
    if success:
        print(f"\n[성공] 경로가 작동함!")
    else:
        print(f"\n[실패] 추가 조사 필요")

if __name__ == "__main__":
    main()
