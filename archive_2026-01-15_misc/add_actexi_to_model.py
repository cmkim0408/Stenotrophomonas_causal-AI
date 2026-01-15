#!/usr/bin/env python
"""
ACtexi 반응을 모델에 추가
"""

import cobra
from pathlib import Path

def load_model(model_path):
    model = cobra.io.read_sbml_model(model_path)
    return model

def add_actexi(new_model, ref_model):
    """ACtexi 반응 추가"""
    rxn_id = 'ACtexi'
    
    if rxn_id in [r.id for r in new_model.reactions]:
        print(f"[ACtexi] 이미 모델에 존재")
        rxn = new_model.reactions.get_by_id(rxn_id)
        print(f"  반응식: {rxn.reaction}")
        print(f"  bounds: [{rxn.lower_bound}, {rxn.upper_bound}]")
        return False
    
    if rxn_id not in ref_model.reactions:
        print(f"[ACtexi] 레퍼런스 모델에 없음")
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
        
        print(f"[ACtexi] 추가 완료")
        print(f"  반응식: {new_rxn.reaction}")
        print(f"  bounds: [{new_rxn.lower_bound}, {new_rxn.upper_bound}]")
        return True
        
    except Exception as e:
        print(f"[ACtexi] 추가 실패: {e}")
        return False

def main():
    base_path = Path(__file__).parent.parent
    new_model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    
    print("="*80)
    print("ACtexi 반응 추가")
    print("="*80)
    
    # 모델 로드
    new_model = load_model(str(new_model_path))
    ref_model = load_model(str(ref_model_path))
    
    # ACtexi 추가
    added = add_actexi(new_model, ref_model)
    
    if added:
        # 모델 저장
        output_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_with_ACtexi.xml"
        cobra.io.write_sbml_model(new_model, str(output_path))
        print(f"\n[모델 저장] {output_path}")
        print(f"  -> BaseModel_with_ACtexi.xml로 저장됨")
    else:
        print(f"\n[모델 저장 생략] ACtexi가 이미 존재하거나 추가 실패")

if __name__ == "__main__":
    main()
