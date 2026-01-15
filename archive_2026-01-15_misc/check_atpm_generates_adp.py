#!/usr/bin/env python
"""
ATPM 반응이 ADP를 생성하는지 확인
- ATPM: atp_c + h2o_c --> adp_c + h_c + pi_c
- ATPM이 ADP를 생성하면, ATPM=0일 때 ADP가 부족할 수 있음
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

def check_atpm_reaction(model):
    """ATPM 반응 확인"""
    print("="*80)
    print("ATPM 반응 확인")
    print("="*80)
    
    try:
        atpm_rxn = model.reactions.get_by_id('ATPM')
        print(f"\n[ATPM 반응]")
        print(f"  ID: {atpm_rxn.id}")
        print(f"  이름: {atpm_rxn.name}")
        print(f"  bounds: [{atpm_rxn.lower_bound}, {atpm_rxn.upper_bound}]")
        print(f"  반응식: {atpm_rxn.reaction}")
        
        # ATP, ADP 계수 확인
        try:
            atp_c = model.metabolites.get_by_id('atp_c')
            adp_c = model.metabolites.get_by_id('adp_c')
            
            atp_coeff = atpm_rxn.metabolites.get(atp_c, 0)
            adp_coeff = atpm_rxn.metabolites.get(adp_c, 0)
            
            print(f"\n  ATP 계수: {atp_coeff}")
            print(f"  ADP 계수: {adp_coeff}")
            
            if atp_coeff < 0:
                print(f"    -> ATPM은 ATP를 소모함")
            if adp_coeff > 0:
                print(f"    -> ATPM은 ADP를 생성함!")
                print(f"    -> 이것이 핵심! ATPM=0이면 ADP가 생성되지 않음")
        except KeyError:
            pass
        
    except KeyError:
        print(f"  [오류] ATPM 반응 없음")

def check_adp_balance_with_growth(model, atpm_value):
    """Growth를 목적함수로 할 때 ADP 균형 확인"""
    model_test = model.copy()
    atpm_rxn = model_test.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = atpm_value
    atpm_rxn.upper_bound = 1000
    
    model_test.objective = 'Growth'
    solution = model_test.optimize()
    
    try:
        adp_c = model_test.metabolites.get_by_id('adp_c')
        
        total_prod = 0
        total_cons = 0
        
        for rxn in adp_c.reactions:
            flux = solution.fluxes.get(rxn.id, 0.0)
            coeff = rxn.metabolites.get(adp_c, 0)
            if coeff > 0:
                total_prod += flux * coeff
            elif coeff < 0:
                total_cons += flux * abs(coeff)
        
        return total_prod, total_cons, solution
        
    except KeyError:
        return None, None, solution

def compare_atpm0_vs_atpm10_adp_balance(model):
    """ATPM=0 vs ATPM=10에서 ADP 균형 비교"""
    print("\n" + "="*80)
    print("ATPM=0 vs ATPM=10에서 ADP 균형 비교 (Growth 목적함수)")
    print("="*80)
    
    # ATPM=0
    prod_0, cons_0, sol_0 = check_adp_balance_with_growth(model, 0)
    
    print(f"\n[ATPM=0일 때]")
    if prod_0 is not None:
        print(f"  ADP 총 생성량: {prod_0:.6f}")
        print(f"  ADP 총 소모량: {cons_0:.6f}")
        print(f"  ADP 균형: {prod_0 - cons_0:.6f}")
        print(f"  성장률: {sol_0.objective_value:.6f}")
        
        # 주요 반응 플럭스
        cs_flux = sol_0.fluxes.get('CS', 0.0)
        sucdi_flux = sol_0.fluxes.get('SUCDi', 0.0)
        atps_flux = sol_0.fluxes.get('ATPS4rpp', 0.0)
        print(f"  CS: {cs_flux:.6f}")
        print(f"  SUCDi: {sucdi_flux:.6f}")
        print(f"  ATPS4rpp: {atps_flux:.6f}")
    
    # ATPM=10
    prod_10, cons_10, sol_10 = check_adp_balance_with_growth(model, 10)
    
    print(f"\n[ATPM=10일 때]")
    if prod_10 is not None:
        print(f"  ADP 총 생성량: {prod_10:.6f}")
        print(f"  ADP 총 소모량: {cons_10:.6f}")
        print(f"  ADP 균형: {prod_10 - cons_10:.6f}")
        print(f"  성장률: {sol_10.objective_value:.6f}")
        
        # 주요 반응 플럭스
        cs_flux = sol_10.fluxes.get('CS', 0.0)
        sucdi_flux = sol_10.fluxes.get('SUCDi', 0.0)
        atps_flux = sol_10.fluxes.get('ATPS4rpp', 0.0)
        print(f"  CS: {cs_flux:.6f}")
        print(f"  SUCDi: {sucdi_flux:.6f}")
        print(f"  ATPS4rpp: {atps_flux:.6f}")
    
    # 비교
    print(f"\n[비교]")
    if prod_0 is not None and prod_10 is not None:
        print(f"  ADP 생성량 차이: {prod_10 - prod_0:.6f}")
        print(f"  ADP 소모량 차이: {cons_10 - cons_0:.6f}")
        
        if abs(prod_10 - prod_0) > 1e-6:
            print(f"  [핵심 발견] ATPM=10일 때 ADP 생성량이 더 많음!")
            print(f"  -> ATPM이 ADP를 생성하므로, ATPM=0이면 ADP가 부족할 수 있음")

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
    
    # ATPM 반응 확인
    check_atpm_reaction(model)
    
    # ADP 균형 비교
    compare_atpm0_vs_atpm10_adp_balance(model)
    
    print("\n" + "="*80)
    print("결론")
    print("="*80)
    print(f"\n[핵심 발견]")
    print(f"  ATPM 반응식: atp_c + h2o_c --> adp_c + h_c + pi_c")
    print(f"  -> ATPM은 ADP를 생성함!")
    print(f"  -> ATPM=0이면 ADP가 생성되지 않아서 ATP 생성 경로가 작동하지 않을 수 있음")
    print(f"  -> ATP 생성 경로(ATPS4rpp, PYK 등)는 ADP를 필요로 함")

if __name__ == "__main__":
    main()
