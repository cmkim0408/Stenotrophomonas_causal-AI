#!/usr/bin/env python
"""
ATPM=0일 때 ATP가 생성되지 않는 실제 원인 찾기
- ATPM은 ATP 소모량인데, 소모가 없는데도 왜 ATP가 생성되지 않는지
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

def check_atp_generation_without_atpm(model):
    """ATPM 없이 ATP 생성 확인"""
    print("="*80)
    print("ATPM 없이 ATP 생성 확인")
    print("="*80)
    
    # ATPM 반응 확인 (비활성화하지 않음, lower_bound만 0)
    atpm_rxn = model.reactions.get_by_id('ATPM')
    original_lb = atpm_rxn.lower_bound
    original_ub = atpm_rxn.upper_bound
    atpm_rxn.lower_bound = 0
    # upper_bound는 그대로 유지
    
    # ATP demand로 최대 생산량 확인
    try:
        atp_c = model.metabolites.get_by_id('atp_c')
        atp_demand = cobra.Reaction('DM_atp_c')
        atp_demand.lower_bound = 0
        atp_demand.upper_bound = 1000
        atp_demand.add_metabolites({atp_c: -1.0})
        model.add_reactions([atp_demand])
        
        model.objective = 'DM_atp_c'
        solution = model.optimize()
        
        print(f"\n[ATP 최대 생산량 (ATPM=0)]")
        print(f"  상태: {solution.status}")
        print(f"  ATP 최대 생산량: {solution.objective_value:.6f}")
        
        if solution.objective_value < 1e-6:
            print(f"  [문제] ATP 생성 불가")
            
            # ATP 생성 반응 확인
            print(f"\n[ATP 생성 반응 확인]")
            atp_producing_rxns = ['ATPS4rpp', 'PGK', 'PYK', 'ACKr', 'PTAr']
            for rxn_id in atp_producing_rxns:
                if rxn_id in model.reactions:
                    rxn = model.reactions.get_by_id(rxn_id)
                    flux = solution.fluxes.get(rxn_id, 0.0)
                    bounds = f"[{rxn.lower_bound}, {rxn.upper_bound}]"
                    print(f"  {rxn_id:20s}: 플럭스={flux:>12.6f}, bounds={bounds}")
                    
                    # ATP 계수 확인
                    atp_coeff = rxn.metabolites.get(atp_c, 0)
                    if atp_coeff > 0:
                        print(f"    -> ATP 생성 (계수: {atp_coeff})")
                        print(f"    -> 반응식: {rxn.reaction[:100]}")
            
            # ATP 생성에 필요한 것들 확인
            print(f"\n[ATP 생성에 필요한 것들 확인]")
            
            # NADH 확인
            try:
                nadh_c = model.metabolites.get_by_id('nadh_c')
                nadh_demand = cobra.Reaction('DM_nadh_c')
                nadh_demand.lower_bound = 0
                nadh_demand.upper_bound = 1000
                nadh_demand.add_metabolites({nadh_c: -1.0})
                model_test = model.copy()
                model_test.add_reactions([nadh_demand])
                model_test.objective = nadh_demand.id
                sol_nadh = model_test.optimize()
                nadh_max = sol_nadh.objective_value
                print(f"  NADH 최대 생산량: {nadh_max:.6f}")
            except:
                pass
            
            # 전자전달사슬 확인
            etc_rxns = ['NADH16pp', 'CYTBO3_4pp', 'ATPS4rpp']
            print(f"\n[전자전달사슬 반응 확인]")
            for rxn_id in etc_rxns:
                if rxn_id in model.reactions:
                    rxn = model.reactions.get_by_id(rxn_id)
                    flux = solution.fluxes.get(rxn_id, 0.0)
                    bounds = f"[{rxn.lower_bound}, {rxn.upper_bound}]"
                    print(f"  {rxn_id:20s}: 플럭스={flux:>12.6f}, bounds={bounds}")
                    if abs(flux) < 1e-6:
                        print(f"    -> [문제] 플럭스 0")
                        print(f"    -> 반응식: {rxn.reaction[:100]}")
        else:
            print(f"  [OK] ATP 생성 가능")
        
        model.remove_reactions([atp_demand])
        # ATPM bounds 복원
        atpm_rxn.lower_bound = original_lb
        
    except KeyError:
        print(f"  [오류] atp_c 메타볼라이트 없음")

def compare_atpm0_vs_atpm10(model):
    """ATPM=0 vs ATPM=10 비교"""
    print("\n" + "="*80)
    print("ATPM=0 vs ATPM=10 비교")
    print("="*80)
    
    # ATPM=0
    model_0 = model.copy()
    atpm_0 = model_0.reactions.get_by_id('ATPM')
    atpm_0.lower_bound = 0
    atpm_0.upper_bound = 0
    
    model_0.objective = 'Growth'
    sol_0 = model_0.optimize()
    
    # ATPM=10
    model_10 = model.copy()
    atpm_10 = model_10.reactions.get_by_id('ATPM')
    atpm_10.lower_bound = 10
    atpm_10.upper_bound = 1000
    # 이전에 upper_bound를 0으로 설정했을 수 있으므로 다시 설정
    
    model_10.objective = 'Growth'
    sol_10 = model_10.optimize()
    
    print(f"\n[성장률 비교]")
    print(f"  ATPM=0: {sol_0.objective_value:.6f} (상태: {sol_0.status})")
    print(f"  ATPM=10: {sol_10.objective_value:.6f} (상태: {sol_10.status})")
    
    # 주요 반응 플럭스 비교
    key_reactions = ['CS', 'SUCDi', 'ATPS4rpp', 'NADH16pp', 'CYTBO3_4pp', 'ACS_ADP']
    print(f"\n[주요 반응 플럭스 비교]")
    print(f"{'반응':<20} {'ATPM=0':<15} {'ATPM=10':<15}")
    print("-" * 50)
    for rxn_id in key_reactions:
        if rxn_id in model.reactions:
            flux_0 = sol_0.fluxes.get(rxn_id, 0.0) if sol_0.status == 'optimal' else 0.0
            flux_10 = sol_10.fluxes.get(rxn_id, 0.0) if sol_10.status == 'optimal' else 0.0
            print(f"{rxn_id:<20} {flux_0:<15.6f} {flux_10:<15.6f}")

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
    
    # ATP 생성 확인
    check_atp_generation_without_atpm(model)
    
    # 비교
    compare_atpm0_vs_atpm10(model)
    
    print("\n" + "="*80)
    print("결론")
    print("="*80)
    print(f"\n[핵심 문제]")
    print(f"  ATPM=0일 때 ATP 생성이 0")
    print(f"  -> ATPM이 단순히 ATP 소모량이 아니라, 모델의 제약 조건으로 작동하는 것 같음")
    print(f"  -> 또는 ATP 생성 경로 자체가 ATPM이 있어야 작동하도록 설정되어 있을 수 있음")

if __name__ == "__main__":
    main()
