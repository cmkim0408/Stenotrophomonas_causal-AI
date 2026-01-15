#!/usr/bin/env python
"""
ATPM=0일 때 ADP 생성 문제 확인
- ATPM이 없으면 ATP가 소모되지 않아서 ADP가 생성되지 않을 수 있음
- ATP 생성 경로는 ADP를 필요로 함
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

def check_adp_generation(model, atpm_value):
    """특정 ATPM 값에서 ADP 생성 확인"""
    model_test = model.copy()
    atpm_rxn = model_test.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = atpm_value
    atpm_rxn.upper_bound = 1000
    
    try:
        adp_c = model_test.metabolites.get_by_id('adp_c')
        
        # ADP demand로 최대 생산량 확인
        adp_demand = cobra.Reaction('DM_adp_c')
        adp_demand.lower_bound = 0
        adp_demand.upper_bound = 1000
        adp_demand.add_metabolites({adp_c: -1.0})
        model_test.add_reactions([adp_demand])
        
        model_test.objective = 'DM_adp_c'
        solution = model_test.optimize()
        
        adp_max = solution.objective_value if solution.status == 'optimal' else 0.0
        
        # ADP를 생성하는 반응 확인
        adp_producing = []
        for rxn in adp_c.reactions:
            if 'DM_' in rxn.id:
                continue
            flux = solution.fluxes.get(rxn.id, 0.0)
            if abs(flux) > 1e-6:
                coeff = rxn.metabolites.get(adp_c, 0)
                if coeff > 0:
                    adp_producing.append((rxn.id, flux * coeff))
        
        # ADP를 소모하는 반응 확인
        adp_consuming = []
        for rxn in adp_c.reactions:
            if 'DM_' in rxn.id:
                continue
            flux = solution.fluxes.get(rxn.id, 0.0)
            if abs(flux) > 1e-6:
                coeff = rxn.metabolites.get(adp_c, 0)
                if coeff < 0:
                    adp_consuming.append((rxn.id, flux * abs(coeff)))
        
        model_test.remove_reactions([adp_demand])
        
        return adp_max, adp_producing, adp_consuming
        
    except KeyError:
        return None, [], []

def check_atp_generation_needs_adp(model):
    """ATP 생성 경로가 ADP를 필요로 하는지 확인"""
    print("="*80)
    print("ATP 생성 경로가 ADP를 필요로 하는지 확인")
    print("="*80)
    
    try:
        atp_c = model.metabolites.get_by_id('atp_c')
        adp_c = model.metabolites.get_by_id('adp_c')
        
        # 주요 ATP 생성 반응
        atp_producing_rxns = ['ATPS4rpp', 'PGK', 'PYK', 'ACKr', 'PTAr']
        
        print(f"\n[ATP 생성 반응과 ADP 관계]")
        for rxn_id in atp_producing_rxns:
            if rxn_id in model.reactions:
                rxn = model.reactions.get_by_id(rxn_id)
                atp_coeff = rxn.metabolites.get(atp_c, 0)
                adp_coeff = rxn.metabolites.get(adp_c, 0)
                
                if atp_coeff > 0:  # ATP 생성
                    print(f"\n  {rxn_id}:")
                    print(f"    반응식: {rxn.reaction[:100]}")
                    print(f"    ATP 계수: {atp_coeff}")
                    print(f"    ADP 계수: {adp_coeff}")
                    if adp_coeff < 0:
                        print(f"    -> [필수] ADP를 소모함 (계수: {adp_coeff})")
                    elif adp_coeff > 0:
                        print(f"    -> ADP를 생성함 (계수: {adp_coeff})")
                    else:
                        print(f"    -> ADP 관련 없음")
    except KeyError:
        print(f"  [오류] atp_c 또는 adp_c 메타볼라이트 없음")

def compare_atpm0_vs_atpm10_adp(model):
    """ATPM=0 vs ATPM=10에서 ADP 생성 비교"""
    print("\n" + "="*80)
    print("ATPM=0 vs ATPM=10에서 ADP 생성 비교")
    print("="*80)
    
    # ATPM=0
    adp_max_0, adp_prod_0, adp_cons_0 = check_adp_generation(model, 0)
    
    print(f"\n[ATPM=0일 때]")
    if adp_max_0 is not None:
        print(f"  ADP 최대 생산량: {adp_max_0:.6f}")
        print(f"  ADP 생성 반응 (플럭스 > 0): {len(adp_prod_0)}개")
        for rxn_id, net_flux in sorted(adp_prod_0, key=lambda x: x[1], reverse=True)[:5]:
            print(f"    {rxn_id}: {net_flux:.6f}")
        print(f"  ADP 소모 반응 (플럭스 > 0): {len(adp_cons_0)}개")
        for rxn_id, net_flux in sorted(adp_cons_0, key=lambda x: x[1], reverse=True)[:5]:
            print(f"    {rxn_id}: {net_flux:.6f}")
    else:
        print(f"  [오류] ADP 확인 실패")
    
    # ATPM=10
    adp_max_10, adp_prod_10, adp_cons_10 = check_adp_generation(model, 10)
    
    print(f"\n[ATPM=10일 때]")
    if adp_max_10 is not None:
        print(f"  ADP 최대 생산량: {adp_max_10:.6f}")
        print(f"  ADP 생성 반응 (플럭스 > 0): {len(adp_prod_10)}개")
        for rxn_id, net_flux in sorted(adp_prod_10, key=lambda x: x[1], reverse=True)[:5]:
            print(f"    {rxn_id}: {net_flux:.6f}")
        print(f"  ADP 소모 반응 (플럭스 > 0): {len(adp_cons_10)}개")
        for rxn_id, net_flux in sorted(adp_cons_10, key=lambda x: x[1], reverse=True)[:5]:
            print(f"    {rxn_id}: {net_flux:.6f}")
    else:
        print(f"  [오류] ADP 확인 실패")
    
    # 비교
    print(f"\n[비교]")
    if adp_max_0 is not None and adp_max_10 is not None:
        if adp_max_0 < 1e-6 and adp_max_10 > 1e-6:
            print(f"  [핵심 발견] ATPM=0일 때 ADP 생성 불가!")
            print(f"  -> ATPM=10일 때 ADP 생성 가능")
            print(f"  -> 이것이 경로가 작동하지 않는 원인!")

def test_with_adp_bootstrap(model):
    """ADP 부트스트랩 제공 시 경로 작동 확인"""
    print("\n" + "="*80)
    print("ADP 부트스트랩 제공 시 경로 작동 확인")
    print("="*80)
    
    model_test = model.copy()
    atpm_rxn = model_test.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    # ADP 부트스트랩 제공
    try:
        adp_c = model_test.metabolites.get_by_id('adp_c')
        ex_adp_id = 'EX_adp_c'
        
        if ex_adp_id not in [r.id for r in model_test.exchanges]:
            ex_adp = cobra.Reaction(ex_adp_id)
            ex_adp.name = 'ADP exchange'
            ex_adp.lower_bound = -0.001
            ex_adp.upper_bound = 1000
            
            adp_e_id = 'adp_e'
            try:
                adp_e = model_test.metabolites.get_by_id(adp_e_id)
            except KeyError:
                adp_e = cobra.Metabolite(adp_e_id, name='ADP', compartment='e')
                model_test.add_metabolites([adp_e])
            
            ex_adp.add_metabolites({adp_e: -1.0})
            model_test.add_reactions([ex_adp])
            
            # Transport 추가
            trans_id = 'T_adp_e_to_adp_c'
            if trans_id not in [r.id for r in model_test.reactions]:
                trans = cobra.Reaction(trans_id)
                trans.name = 'ADP transport'
                trans.lower_bound = -1000
                trans.upper_bound = 1000
                trans.add_metabolites({adp_e: -1.0, adp_c: 1.0})
                model_test.add_reactions([trans])
        
        print(f"  ADP 부트스트랩 제공: EX_adp_c 하한=-0.001")
        
        model_test.objective = 'Growth'
        solution = model_test.optimize()
        
        print(f"\n[FBA 결과]")
        print(f"  상태: {solution.status}")
        print(f"  성장률: {solution.objective_value:.6f}")
        
        cs_flux = solution.fluxes.get('CS', 0.0)
        sucdi_flux = solution.fluxes.get('SUCDi', 0.0)
        atps_flux = solution.fluxes.get('ATPS4rpp', 0.0)
        acs_adp_flux = solution.fluxes.get('ACS_ADP', 0.0)
        
        print(f"\n[주요 반응 플럭스]")
        print(f"  CS: {cs_flux:.6f}")
        print(f"  SUCDi: {sucdi_flux:.6f}")
        print(f"  ATPS4rpp: {atps_flux:.6f}")
        print(f"  ACS_ADP: {acs_adp_flux:.6f}")
        
        if abs(cs_flux) > 1e-6 or abs(sucdi_flux) > 1e-6:
            print(f"\n[성공] ADP 부트스트랩으로 경로 작동!")
            return True
        else:
            print(f"\n[실패] ADP 부트스트랩으로도 경로 작동 안 함")
            return False
            
    except KeyError:
        print(f"  [오류] adp_c 메타볼라이트 없음")
        return False

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
    
    # ATP 생성 경로가 ADP를 필요로 하는지 확인
    check_atp_generation_needs_adp(model)
    
    # ATPM=0 vs ATPM=10에서 ADP 생성 비교
    compare_atpm0_vs_atpm10_adp(model)
    
    # ADP 부트스트랩 테스트
    success = test_with_adp_bootstrap(model)
    
    print("\n" + "="*80)
    print("결론")
    print("="*80)
    
    if success:
        print(f"\n[핵심 발견]")
        print(f"  ATPM=0일 때 ADP 생성 불가")
        print(f"  -> ATP 생성 경로는 ADP를 필요로 함")
        print(f"  -> ADP 부트스트랩 제공 시 경로 작동!")
        print(f"  -> 이것이 ATPM=0에서 경로가 작동하지 않는 원인!")
    else:
        print(f"\n[추가 조사 필요]")
        print(f"  ADP 부트스트랩으로도 경로가 작동하지 않음")
        print(f"  -> 다른 문제가 있을 수 있음")

if __name__ == "__main__":
    main()
