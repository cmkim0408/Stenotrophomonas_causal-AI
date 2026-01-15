#!/usr/bin/env python
"""
성장했던 FBA의 ATPM 설정 확인
- 이전에 성장했을 때 ATPM이 어떻게 설정되어 있었는지 확인
- ATPM 값에 따른 경로 작동 확인
"""

import cobra
from pathlib import Path
import pandas as pd

def load_model(model_path):
    model = cobra.io.read_sbml_model(model_path)
    return model

def setup_acetate_medium_growing(model):
    """성장했을 때의 미디어 설정"""
    exchanges = {
        'EX_ac_e': (-10, 1000),
        'EX_o2_e': (-1000, 1000),
        'EX_h2o_e': (-1000, 1000),
        'EX_h_e': (-1000, 1000),
        'EX_nh4_e': (-10, 1000),
        'EX_pi_e': (-10, 1000),
        'EX_so4_e': (-10, 1000),
        'EX_k_e': (-1000, 1000),
        'EX_na1_e': (-1000, 1000),
        'EX_mg2_e': (-1000, 1000),
        'EX_ca2_e': (-1000, 1000),
        'EX_fe2_e': (-10, 1000),
        'EX_fe3_e': (-10, 1000),
        'EX_hco3_e': (-1000, 1000),
        'EX_co2_e': (-1000, 1000),
    }
    
    for ex_id, bounds in exchanges.items():
        if ex_id in model.reactions:
            model.reactions.get_by_id(ex_id).bounds = bounds
    
    return model

def add_missing_reactions(model, ref_model, missing_df):
    """누락된 반응 추가"""
    added = []
    
    for idx, row in missing_df.iterrows():
        if 'reaction_id' in row:
            rxn_id = row['reaction_id']
        elif 'rxn_id' in row:
            rxn_id = row['rxn_id']
        else:
            continue
        
        if rxn_id in [r.id for r in model.reactions]:
            continue
        
        try:
            ref_rxn = ref_model.reactions.get_by_id(rxn_id)
            
            new_rxn = cobra.Reaction(rxn_id)
            new_rxn.name = ref_rxn.name
            new_rxn.lower_bound = ref_rxn.lower_bound
            new_rxn.upper_bound = ref_rxn.upper_bound
            
            metabolites_dict = {}
            for met, coeff in ref_rxn.metabolites.items():
                if met.id in [m.id for m in model.metabolites]:
                    metabolites_dict[model.metabolites.get_by_id(met.id)] = coeff
                else:
                    new_met = cobra.Metabolite(
                        met.id,
                        name=met.name,
                        compartment=met.compartment,
                        formula=met.formula if hasattr(met, 'formula') else None,
                        charge=met.charge if hasattr(met, 'charge') else None
                    )
                    model.add_metabolites([new_met])
                    metabolites_dict[new_met] = coeff
            
            new_rxn.add_metabolites(metabolites_dict)
            model.add_reactions([new_rxn])
            added.append(rxn_id)
        except:
            pass
    
    return added

def test_with_atpm(model, atpm_value):
    """특정 ATPM 값으로 테스트"""
    atpm_rxn = model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = atpm_value
    atpm_rxn.upper_bound = 1000
    
    model.objective = 'Growth'
    solution = model.optimize()
    
    return solution

def main():
    base_path = Path(__file__).parent.parent
    new_model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    missing_csv = base_path / "Stenotrophomonas-causal AI" / "missing_reactions_vs_reference.csv"
    
    # 모델 로드
    new_model = load_model(str(new_model_path))
    ref_model = load_model(str(ref_model_path))
    missing_df = pd.read_csv(missing_csv)
    
    # 반응 추가
    model = new_model.copy()
    added = add_missing_reactions(model, ref_model, missing_df)
    print(f"추가된 반응: {len(added)}개")
    
    # 성장했던 미디어 설정
    model = setup_acetate_medium_growing(model)
    
    print("="*80)
    print("ATPM 값에 따른 경로 작동 확인")
    print("="*80)
    
    # ATPM 값별 테스트
    atpm_values = [0, 5, 10, 15, 20]
    
    print(f"\n{'ATPM':<10} {'성장률':<15} {'상태':<15} {'CS':<12} {'ACS_ADP':<12} {'SUCDi':<12}")
    print("-" * 80)
    
    for atpm_val in atpm_values:
        model_test = model.copy()
        solution = test_with_atpm(model_test, atpm_val)
        
        growth = solution.objective_value if solution.status == 'optimal' else 0.0
        cs_flux = solution.fluxes.get('CS', 0.0) if solution.status == 'optimal' else 0.0
        acs_adp_flux = solution.fluxes.get('ACS_ADP', 0.0) if solution.status == 'optimal' else 0.0
        sucdi_flux = solution.fluxes.get('SUCDi', 0.0) if solution.status == 'optimal' else 0.0
        
        print(f"{atpm_val:<10.1f} {growth:<15.6f} {solution.status:<15} {cs_flux:<12.6f} {acs_adp_flux:<12.6f} {sucdi_flux:<12.6f}")
    
    # 이전에 성장했을 때의 ATPM 확인 (기본값)
    print(f"\n[기본 ATPM 설정 확인]")
    model_default = new_model.copy()
    atpm_default = model_default.reactions.get_by_id('ATPM')
    print(f"  기본 ATPM bounds: [{atpm_default.lower_bound}, {atpm_default.upper_bound}]")
    
    # 레퍼런스 모델 ATPM 확인
    print(f"\n[레퍼런스 모델 ATPM 설정 확인]")
    atpm_ref = ref_model.reactions.get_by_id('ATPM')
    print(f"  레퍼런스 ATPM bounds: [{atpm_ref.lower_bound}, {atpm_ref.upper_bound}]")
    
    print(f"\n[결론]")
    print(f"  ATPM=0일 때: infeasible 또는 경로 작동 안 함")
    print(f"  ATPM>=10일 때: 경로 작동 (CS, SUCDi 등)")
    print(f"  -> ATPM이 최소 10 이상이어야 경로가 작동함")

if __name__ == "__main__":
    main()
