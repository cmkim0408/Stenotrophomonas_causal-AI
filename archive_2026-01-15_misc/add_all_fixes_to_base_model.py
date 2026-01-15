#!/usr/bin/env python
"""
BaseModel_with_ACtexi.xml에 모든 수정사항 적용
- 4개 반응 추가 (ACS_ADP, SUCDi, PEPCK_ATP, ACtexi)
- Exchange bounds 수정
- ATPM bounds 수정
"""

import cobra
from pathlib import Path

def load_model(model_path):
    model = cobra.io.read_sbml_model(model_path)
    return model

def add_reaction(new_model, ref_model, rxn_id):
    """반응 추가"""
    if rxn_id in [r.id for r in new_model.reactions]:
        return False
    
    if rxn_id not in ref_model.reactions:
        return False
    
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
        return True
    except:
        return False

def fix_exchange_bounds(new_model, ref_model):
    """Exchange bounds를 레퍼런스 모델과 동일하게 설정"""
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
            
            if new_rxn.lower_bound != ref_rxn.lower_bound or new_rxn.upper_bound != ref_rxn.upper_bound:
                new_rxn.lower_bound = ref_rxn.lower_bound
                new_rxn.upper_bound = ref_rxn.upper_bound
                fixed.append(ex_id)
        elif ex_id in ref_model.reactions and ex_id not in new_model.reactions:
            # Exchange 반응 추가
            ref_rxn = ref_model.reactions.get_by_id(ex_id)
            new_ex = cobra.Reaction(ex_id)
            new_ex.name = ref_rxn.name
            new_ex.lower_bound = ref_rxn.lower_bound
            new_ex.upper_bound = ref_rxn.upper_bound
            
            for met, coeff in ref_rxn.metabolites.items():
                if met.id in [m.id for m in new_model.metabolites]:
                    new_ex.add_metabolites({new_model.metabolites.get_by_id(met.id): coeff})
                else:
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
    
    return fixed

def fix_atpm_bounds(new_model, ref_model):
    """ATPM bounds를 레퍼런스 모델과 동일하게 설정"""
    ref_atpm = ref_model.reactions.get_by_id('ATPM')
    new_atpm = new_model.reactions.get_by_id('ATPM')
    new_atpm.lower_bound = ref_atpm.lower_bound
    new_atpm.upper_bound = ref_atpm.upper_bound

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_with_ACtexi.xml"
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    
    print("="*80)
    print("BaseModel_with_ACtexi.xml에 모든 수정사항 적용")
    print("="*80)
    
    # 모델 로드
    new_model = load_model(str(model_path))
    ref_model = load_model(str(ref_model_path))
    
    # 4개 반응 추가
    reactions_to_add = ['ACS_ADP', 'SUCDi', 'PEPCK_ATP', 'ACtexi']
    
    print(f"\n[반응 추가]")
    added = []
    for rxn_id in reactions_to_add:
        if add_reaction(new_model, ref_model, rxn_id):
            added.append(rxn_id)
            print(f"  {rxn_id}: 추가됨")
        else:
            if rxn_id in [r.id for r in new_model.reactions]:
                print(f"  {rxn_id}: 이미 존재")
            else:
                print(f"  {rxn_id}: 추가 실패 또는 레퍼런스에 없음")
    
    # Exchange bounds 수정
    print(f"\n[Exchange bounds 수정]")
    fixed_exchanges = fix_exchange_bounds(new_model, ref_model)
    print(f"  {len(fixed_exchanges)}개 Exchange 수정/추가됨")
    
    # ATPM bounds 수정
    print(f"\n[ATPM bounds 수정]")
    fix_atpm_bounds(new_model, ref_model)
    atpm_rxn = new_model.reactions.get_by_id('ATPM')
    print(f"  ATPM bounds: [{atpm_rxn.lower_bound}, {atpm_rxn.upper_bound}]")
    
    # 모델 저장
    output_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_with_ACtexi.xml"
    cobra.io.write_sbml_model(new_model, str(output_path))
    print(f"\n[모델 저장 완료] {output_path}")
    
    print("\n" + "="*80)
    print("최종 결과")
    print("="*80)
    print(f"\n[추가된 반응] {len(added)}개: {', '.join(added)}")
    print(f"[수정된 Exchange] {len(fixed_exchanges)}개")
    print(f"[ATPM bounds 수정] 완료")
    print(f"[모델 저장] BaseModel_with_ACtexi.xml")

if __name__ == "__main__":
    main()
