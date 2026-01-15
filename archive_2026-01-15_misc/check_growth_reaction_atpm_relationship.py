#!/usr/bin/env python
"""
Growth 반응과 ATPM의 관계 확인
- ATPM=0일 때 Growth를 목적함수로 하면 경로가 작동하지 않는 이유
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

def check_growth_reaction(model):
    """Growth 반응 확인"""
    print("="*80)
    print("Growth 반응 확인")
    print("="*80)
    
    try:
        growth_rxn = model.reactions.get_by_id('Growth')
        print(f"\n[Growth 반응]")
        print(f"  ID: {growth_rxn.id}")
        print(f"  이름: {growth_rxn.name}")
        print(f"  bounds: [{growth_rxn.lower_bound}, {growth_rxn.upper_bound}]")
        print(f"  반응식: {growth_rxn.reaction}")
        
        # ATP 관련 확인
        try:
            atp_c = model.metabolites.get_by_id('atp_c')
            atp_coeff = growth_rxn.metabolites.get(atp_c, 0)
            if atp_coeff != 0:
                print(f"\n  ATP 계수: {atp_coeff}")
                if atp_coeff < 0:
                    print(f"    -> Growth는 ATP를 소모함")
                elif atp_coeff > 0:
                    print(f"    -> Growth는 ATP를 생성함")
        except KeyError:
            pass
        
        # ATPM 관련 확인
        try:
            atpm_rxn = model.reactions.get_by_id('ATPM')
            atpm_coeff = growth_rxn.metabolites.get(atpm_rxn, 0)
            if atpm_coeff != 0:
                print(f"  ATPM 계수: {atpm_coeff}")
        except:
            pass
        
    except KeyError:
        print(f"  [오류] Growth 반응 없음")

def test_growth_with_different_atpm(model):
    """다양한 ATPM 값에서 Growth 테스트"""
    print("\n" + "="*80)
    print("다양한 ATPM 값에서 Growth 테스트")
    print("="*80)
    
    atpm_values = [0, 1, 5, 10, 20]
    
    print(f"\n{'ATPM':<10} {'Growth':<15} {'CS':<12} {'SUCDi':<12} {'ATPS4rpp':<12}")
    print("-" * 60)
    
    for atpm_val in atpm_values:
        model_test = model.copy()
        atpm_rxn = model_test.reactions.get_by_id('ATPM')
        atpm_rxn.lower_bound = atpm_val
        atpm_rxn.upper_bound = 1000
        
        model_test.objective = 'Growth'
        solution = model_test.optimize()
        
        growth = solution.objective_value if solution.status == 'optimal' else 0.0
        cs_flux = solution.fluxes.get('CS', 0.0) if solution.status == 'optimal' else 0.0
        sucdi_flux = solution.fluxes.get('SUCDi', 0.0) if solution.status == 'optimal' else 0.0
        atps_flux = solution.fluxes.get('ATPS4rpp', 0.0) if solution.status == 'optimal' else 0.0
        
        print(f"{atpm_val:<10.1f} {growth:<15.6f} {cs_flux:<12.6f} {sucdi_flux:<12.6f} {atps_flux:<12.6f}")

def check_atp_balance_with_growth(model):
    """Growth를 목적함수로 할 때 ATP 균형 확인"""
    print("\n" + "="*80)
    print("Growth를 목적함수로 할 때 ATP 균형 확인")
    print("="*80)
    
    # ATPM=0
    model_0 = model.copy()
    atpm_0 = model_0.reactions.get_by_id('ATPM')
    atpm_0.lower_bound = 0
    atpm_0.upper_bound = 1000
    
    model_0.objective = 'Growth'
    sol_0 = model_0.optimize()
    
    # ATPM=10
    model_10 = model.copy()
    atpm_10 = model_10.reactions.get_by_id('ATPM')
    atpm_10.lower_bound = 10
    atpm_10.upper_bound = 1000
    
    model_10.objective = 'Growth'
    sol_10 = model_10.optimize()
    
    print(f"\n[ATPM=0일 때 ATP 균형]")
    try:
        atp_c = model_0.metabolites.get_by_id('atp_c')
        total_prod = 0
        total_cons = 0
        
        for rxn in atp_c.reactions:
            flux = sol_0.fluxes.get(rxn.id, 0.0)
            coeff = rxn.metabolites.get(atp_c, 0)
            if coeff > 0:
                total_prod += flux * coeff
            elif coeff < 0:
                total_cons += flux * abs(coeff)
        
        print(f"  ATP 총 생성량: {total_prod:.6f}")
        print(f"  ATP 총 소모량: {total_cons:.6f}")
        print(f"  ATP 균형: {total_prod - total_cons:.6f}")
    except:
        pass
    
    print(f"\n[ATPM=10일 때 ATP 균형]")
    try:
        atp_c = model_10.metabolites.get_by_id('atp_c')
        total_prod = 0
        total_cons = 0
        
        for rxn in atp_c.reactions:
            flux = sol_10.fluxes.get(rxn.id, 0.0)
            coeff = rxn.metabolites.get(atp_c, 0)
            if coeff > 0:
                total_prod += flux * coeff
            elif coeff < 0:
                total_cons += flux * abs(coeff)
        
        print(f"  ATP 총 생성량: {total_prod:.6f}")
        print(f"  ATP 총 소모량: {total_cons:.6f}")
        print(f"  ATP 균형: {total_prod - total_cons:.6f}")
    except:
        pass

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
    
    # Growth 반응 확인
    check_growth_reaction(model)
    
    # 다양한 ATPM 값에서 테스트
    test_growth_with_different_atpm(model)
    
    # ATP 균형 확인
    check_atp_balance_with_growth(model)
    
    print("\n" + "="*80)
    print("결론")
    print("="*80)
    print(f"\n[핵심 발견]")
    print(f"  ATPM=0일 때 ATP는 생성 가능하지만, Growth를 목적함수로 하면 경로가 작동하지 않음")
    print(f"  ATPM>=5일 때부터 경로가 작동하기 시작함")
    print(f"  -> ATPM이 단순히 ATP 소모량이 아니라, 모델의 제약 조건으로 작동하는 것 같음")

if __name__ == "__main__":
    main()
