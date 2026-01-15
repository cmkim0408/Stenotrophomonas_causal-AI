#!/usr/bin/env python
"""
ATPM=0일 때 경로가 작동하지 않는 실제 원인 찾기
- ATPM은 ATP 소모량인데, 소모가 없는데도 왜 작동 안 하는지
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

def diagnose_atpm0_blockage(model):
    """ATPM=0일 때 경로가 막히는 실제 원인 진단"""
    print("="*80)
    print("ATPM=0일 때 경로 막힘 원인 진단")
    print("="*80)
    
    # ATPM=0 설정
    atpm_rxn = model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    model.objective = 'Growth'
    solution = model.optimize()
    
    print(f"\n[FBA 결과]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    
    # ATPM 반응 확인
    print(f"\n[ATPM 반응 확인]")
    print(f"  ATPM bounds: [{atpm_rxn.lower_bound}, {atpm_rxn.upper_bound}]")
    print(f"  ATPM 반응식: {atpm_rxn.reaction}")
    atpm_flux = solution.fluxes.get('ATPM', 0.0)
    print(f"  ATPM 플럭스: {atpm_flux:.6f}")
    
    # ATP 생성/소모 확인
    print(f"\n[ATP 생성/소모 확인]")
    try:
        atp_c = model.metabolites.get_by_id('atp_c')
        
        # ATP를 생성하는 반응
        atp_producing = []
        for rxn in atp_c.reactions:
            if 'ATPM' in rxn.id:
                continue
            flux = solution.fluxes.get(rxn.id, 0.0)
            if abs(flux) > 1e-6:
                coeff = rxn.metabolites.get(atp_c, 0)
                if coeff > 0:
                    atp_producing.append((rxn.id, flux * coeff))
        
        # ATP를 소모하는 반응
        atp_consuming = []
        for rxn in atp_c.reactions:
            if 'ATPM' in rxn.id:
                continue
            flux = solution.fluxes.get(rxn.id, 0.0)
            if abs(flux) > 1e-6:
                coeff = rxn.metabolites.get(atp_c, 0)
                if coeff < 0:
                    atp_consuming.append((rxn.id, flux * abs(coeff)))
        
        print(f"  ATP 생성 반응 (플럭스 > 0): {len(atp_producing)}개")
        for rxn_id, net_flux in sorted(atp_producing, key=lambda x: x[1], reverse=True)[:5]:
            print(f"    {rxn_id}: {net_flux:.6f}")
        
        print(f"  ATP 소모 반응 (플럭스 > 0): {len(atp_consuming)}개")
        for rxn_id, net_flux in sorted(atp_consuming, key=lambda x: x[1], reverse=True)[:5]:
            print(f"    {rxn_id}: {net_flux:.6f}")
        
        # ATP 균형 확인
        total_production = sum([f for _, f in atp_producing])
        total_consumption = sum([f for _, f in atp_consuming])
        print(f"\n  ATP 총 생성량: {total_production:.6f}")
        print(f"  ATP 총 소모량: {total_consumption:.6f}")
        print(f"  ATP 균형: {total_production - total_consumption:.6f}")
        
    except KeyError:
        print(f"  [오류] atp_c 메타볼라이트 없음")
    
    # 주요 경로 반응 확인
    print(f"\n[주요 경로 반응 확인]")
    key_reactions = ['ACS_ADP', 'ACS', 'CS', 'ICL', 'MALS', 'SUCDi', 'PEPCK_ATP']
    for rxn_id in key_reactions:
        if rxn_id in model.reactions:
            rxn = model.reactions.get_by_id(rxn_id)
            flux = solution.fluxes.get(rxn_id, 0.0)
            bounds = f"[{rxn.lower_bound}, {rxn.upper_bound}]"
            print(f"  {rxn_id:20s}: 플럭스={flux:>12.6f}, bounds={bounds}")
            
            # ATP 관련 확인
            if 'atp_c' in [m.id for m in rxn.metabolites]:
                atp_coeff = rxn.metabolites.get(model.metabolites.get_by_id('atp_c'), 0)
                if atp_coeff != 0:
                    print(f"    -> ATP 계수: {atp_coeff}")
    
    # Blocked reactions 확인
    print(f"\n[Blocked reactions 확인]")
    try:
        blocked = cobra.flux_analysis.find_blocked_reactions(model, open_exchanges=True)
        key_blocked = [r for r in blocked if r in key_reactions]
        if key_blocked:
            print(f"  막힌 주요 반응: {key_blocked}")
        else:
            print(f"  주요 반응은 blocked 아님")
    except:
        pass
    
    # ATPM=10으로 설정했을 때 비교
    print(f"\n[ATPM=10으로 설정했을 때 비교]")
    model_test = model.copy()
    atpm_test = model_test.reactions.get_by_id('ATPM')
    atpm_test.lower_bound = 10
    atpm_test.upper_bound = 1000
    
    model_test.objective = 'Growth'
    solution_test = model_test.optimize()
    
    print(f"  상태: {solution_test.status}")
    print(f"  성장률: {solution_test.objective_value:.6f}")
    
    for rxn_id in key_reactions:
        if rxn_id in model_test.reactions:
            flux = solution_test.fluxes.get(rxn_id, 0.0)
            if abs(flux) > 1e-6:
                print(f"  {rxn_id:20s}: {flux:>12.6f} (ATPM=10일 때 작동)")

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
    
    # 진단
    diagnose_atpm0_blockage(model)

if __name__ == "__main__":
    main()
