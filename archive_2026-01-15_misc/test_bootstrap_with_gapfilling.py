#!/usr/bin/env python
"""
부트스트랩 문제 해결 테스트 (Gap-filling 후)
논문 방법론 적용: ATPM을 설정하고 FBA 수행
"""

import cobra
import pandas as pd

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    print(f"[OK] 모델 로드: {model.id}")
    return model

def find_biomass_reaction(model):
    for name in ['Growth', 'BIOMASS']:
        try:
            return model.reactions.get_by_id(name)
        except KeyError:
            continue
    return None

def find_atp_maintenance_reaction(model):
    """ATP Maintenance 반응 찾기"""
    atpm_candidates = ['ATPM', 'ATPM_c', 'ATPS4rpp', 'ATPS4r']
    
    for atpm_id in atpm_candidates:
        try:
            return model.reactions.get_by_id(atpm_id)
        except KeyError:
            continue
    
    # ATP maintenance 반응 생성
    try:
        atp_c = model.metabolites.get_by_id('atp_c')
        adp_c = model.metabolites.get_by_id('adp_c')
        pi_c = model.metabolites.get_by_id('pi_c')
        h_c = model.metabolites.get_by_id('h_c')
        
        atpm_rxn = cobra.Reaction('ATPM')
        atpm_rxn.name = 'ATP maintenance'
        atpm_rxn.lower_bound = 0
        atpm_rxn.upper_bound = 1000
        atpm_rxn.add_metabolites({
            atp_c: -1,
            adp_c: 1,
            pi_c: 1,
            h_c: 1
        })
        model.add_reactions([atpm_rxn])
        print(f"[OK] ATPM 반응 생성")
        return atpm_rxn
    except KeyError as e:
        print(f"[ERROR] ATPM 생성 실패: {e}")
        return None

def setup_acetate_medium(model, acetate_uptake=-19.0, oxygen_uptake=-100.0):
    """Acetate minimal medium 설정 (논문 방법)"""
    print("\n" + "="*70)
    print("Acetate Minimal Medium 설정 (논문 방법)")
    print("="*70)
    
    # 모든 exchange 차단
    for rxn in model.exchanges:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    
    # Acetate 설정 (논문 값)
    try:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.lower_bound = acetate_uptake
        ex_ac.upper_bound = acetate_uptake
        print(f"[OK] Acetate uptake: {acetate_uptake} mmol/gDCW/h")
    except KeyError:
        print("[WARNING] EX_ac_e 없음")
    
    # Oxygen 설정
    try:
        ex_o2 = model.reactions.get_by_id('EX_o2_e')
        ex_o2.lower_bound = oxygen_uptake
        ex_o2.upper_bound = 1000
        print(f"[OK] Oxygen uptake: {oxygen_uptake} ~ 1000 mmol/gDCW/h")
    except KeyError:
        print("[WARNING] EX_o2_e 없음")
    
    # 필수 영양소
    essentials = ['EX_nh4_e', 'EX_h2o_e', 'EX_h_e', 'EX_pi_e', 'EX_so4_e',
                  'EX_k_e', 'EX_na1_e', 'EX_mg2_e', 'EX_ca2_e', 'EX_fe2_e',
                  'EX_mn2_e', 'EX_zn2_e', 'EX_co2_e']
    
    for ex_id in essentials:
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            if ex_id == 'EX_co2_e':
                ex_rxn.lower_bound = -1000
                ex_rxn.upper_bound = 1000
            else:
                ex_rxn.lower_bound = -1000
        except KeyError:
            pass
    
    print(f"[OK] 필수 영양소 {len(essentials)}개 설정")
    return model

def test_growth_with_atpm(model, biomass_rxn, atpm_rxn, atpm_value=0):
    """ATPM 설정 후 성장 테스트"""
    print("\n" + "="*70)
    print(f"성장 테스트 (ATPM = {atpm_value} mmol ATP/gDCW/h)")
    print("="*70)
    
    # ATPM 설정
    atpm_rxn.lower_bound = atpm_value
    atpm_rxn.upper_bound = max(atpm_value, 1000)
    
    model.objective = biomass_rxn.id
    
    solution = model.optimize()
    
    print(f"\n최적화 상태: {solution.status}")
    
    if solution.status == 'optimal':
        biomass_flux = solution.fluxes.get(biomass_rxn.id, 0)
        print(f"Biomass flux: {biomass_flux:.6f} 1/h")
        
        if biomass_flux > 1e-6:
            print("\n[SUCCESS] 성장 가능!")
            
            # Exchange 플럭스
            print("\nExchange 플럭스 (절대값 > 0.001):")
            exchange_fluxes = []
            for rxn in model.exchanges:
                flux = solution.fluxes.get(rxn.id, 0)
                if abs(flux) > 0.001:
                    exchange_fluxes.append((rxn.id, flux))
            
            if exchange_fluxes:
                for rxn_id, flux in sorted(exchange_fluxes, key=lambda x: abs(x[1]), reverse=True):
                    direction = "Uptake" if flux < 0 else "Secretion"
                    print(f"  {rxn_id}: {flux:.6f} ({direction})")
            
            # 주요 경로 플럭스
            print("\n주요 경로 반응 플럭스:")
            key_reactions = {
                'EX_ac_e': 'Acetate Exchange',
                'ACt': 'Acetate Transport',
                'SUCOAACTr': 'Acetate-CoA Transferase',
                'ACS': 'Acetyl-CoA Synthetase',
                'CS': 'Citrate Synthase',
                'ICL': 'Isocitrate Lyase',
                'MALS': 'Malate Synthase',
                'SUCD': 'Succinate Dehydrogenase',
                'FUM': 'Fumarase',
                'MDH': 'Malate Dehydrogenase',
                'Growth': 'Biomass'
            }
            
            for rxn_id, rxn_name in key_reactions.items():
                try:
                    flux = solution.fluxes.get(rxn_id, 0)
                    if abs(flux) > 1e-8:
                        print(f"  {rxn_id} ({rxn_name}): {flux:.6f}")
                except KeyError:
                    pass
            
            # Yield 계산
            acetate_uptake = abs(solution.fluxes.get('EX_ac_e', 0))
            if acetate_uptake > 0:
                yield_biomass = biomass_flux / acetate_uptake
                print(f"\n성장 속도: {biomass_flux:.6f} 1/h")
                print(f"Acetate 섭취: {acetate_uptake:.6f} mmol/gDCW/h")
                print(f"Biomass yield: {yield_biomass:.6f} gDW/mmol acetate")
            
            return solution, True
        else:
            print("\n[FAIL] 성장 불가 (Biomass flux = 0)")
            return solution, False
    else:
        print(f"\n[ERROR] 최적화 실패: {solution.status}")
        return solution, False

def test_bootstrap_with_demand(model, biomass_rxn):
    """Demand reaction으로 부트스트랩 테스트"""
    print("\n" + "="*70)
    print("부트스트랩 테스트 (Demand Reaction 방법)")
    print("="*70)
    
    # CoA demand reaction 추가
    try:
        coa_c = model.metabolites.get_by_id('coa_c')
        
        try:
            dm_coa = model.reactions.get_by_id('DM_coa_c')
            dm_coa.lower_bound = -0.1
            print(f"[OK] DM_coa_c 존재, LB=-0.1로 설정")
        except KeyError:
            dm_coa = cobra.Reaction('DM_coa_c')
            dm_coa.name = 'CoA demand (bootstrap)'
            dm_coa.lower_bound = -0.1
            dm_coa.upper_bound = 1000
            dm_coa.add_metabolites({coa_c: -1})
            model.add_reactions([dm_coa])
            print(f"[OK] DM_coa_c 추가 (LB=-0.1)")
    except KeyError:
        print("[WARNING] coa_c metabolite 없음")
    
    # FBA 수행
    model.objective = biomass_rxn.id
    solution = model.optimize()
    
    if solution.status == 'optimal':
        biomass_flux = solution.fluxes.get(biomass_rxn.id, 0)
        print(f"Biomass flux: {biomass_flux:.6f} 1/h")
        
        if biomass_flux > 1e-6:
            print("[SUCCESS] Demand reaction으로 성장 가능!")
            return solution, True
    
    return solution, False

def main():
    print("="*70)
    print("부트스트랩 문제 해결 테스트 (Gap-filling 후)")
    print("="*70)
    
    # 모델 로드
    model = load_model("BaseModel.xml")
    
    # Biomass 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass: {biomass_rxn.id}")
    
    # ATP Maintenance 찾기/생성
    atpm_rxn = find_atp_maintenance_reaction(model)
    if not atpm_rxn:
        print("[ERROR] ATPM 반응 없음")
        return
    
    print(f"[OK] ATP Maintenance: {atpm_rxn.id}")
    
    # Medium 설정 (논문 방법)
    model = setup_acetate_medium(model, acetate_uptake=-19.0, oxygen_uptake=-100.0)
    
    # 방법 1: ATPM=0으로 테스트
    print("\n[방법 1] ATPM=0으로 테스트")
    solution1, success1 = test_growth_with_atpm(model, biomass_rxn, atpm_rxn, atpm_value=0)
    
    if not success1:
        # 방법 2: ATPM=10으로 테스트
        print("\n[방법 2] ATPM=10으로 테스트")
        solution2, success2 = test_growth_with_atpm(model, biomass_rxn, atpm_rxn, atpm_value=10)
        
        if not success2:
            # 방법 3: Demand reaction으로 부트스트랩
            print("\n[방법 3] Demand reaction으로 부트스트랩")
            solution3, success3 = test_bootstrap_with_demand(model, biomass_rxn)
    
    print("\n" + "="*70)
    print("부트스트랩 테스트 완료")
    print("="*70)
    
    # 최종 모델 저장
    output_path = "BaseModel.xml"
    cobra.io.write_sbml_model(model, output_path)
    print(f"\n[OK] 모델 저장: {output_path}")

if __name__ == "__main__":
    main()
