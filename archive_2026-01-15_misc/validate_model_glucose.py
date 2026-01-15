#!/usr/bin/env python
"""
모델 검증: 포도당 기반 성장 테스트
모델 자체가 작동하는지 확인
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
    
    return None

def setup_glucose_medium(model, glucose_uptake=-10.0, oxygen_uptake=-100.0):
    """포도당 minimal medium 설정"""
    print("\n" + "="*70)
    print("포도당 Minimal Medium 설정")
    print("="*70)
    
    # 모든 exchange 차단 (상한을 먼저 설정)
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    # 포도당 설정
    glucose_exchanges = ['EX_glc_e', 'EX_glc__D_e', 'EX_glucose_e']
    glucose_set = False
    
    for ex_id in glucose_exchanges:
        try:
            ex_glc = model.reactions.get_by_id(ex_id)
            ex_glc.lower_bound = glucose_uptake
            ex_glc.upper_bound = 1000
            print(f"[OK] Glucose uptake: {ex_id} ({glucose_uptake} ~ 1000 mmol/gDCW/h)")
            glucose_set = True
            break
        except KeyError:
            continue
    
    if not glucose_set:
        print("[WARNING] 포도당 exchange 반응을 찾을 수 없습니다")
        # 패턴으로 찾기
        for rxn in model.exchanges:
            if 'glc' in rxn.id.lower() or 'glucose' in rxn.id.lower():
                rxn.lower_bound = glucose_uptake
                rxn.upper_bound = 1000
                print(f"[OK] Glucose uptake: {rxn.id} ({glucose_uptake} ~ 1000 mmol/gDCW/h)")
                glucose_set = True
                break
    
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
    
    print("\n필수 영양소 설정:")
    for ex_id in essentials:
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            if ex_id == 'EX_co2_e':
                ex_rxn.lower_bound = -1000
                ex_rxn.upper_bound = 1000
            else:
                ex_rxn.lower_bound = -1000
            print(f"  [OK] {ex_id}")
        except KeyError:
            pass
    
    return model

def test_growth_on_glucose(model, biomass_rxn, atpm_rxn=None, atpm_value=0):
    """포도당 기반 성장 테스트"""
    print("\n" + "="*70)
    print("포도당 기반 성장 테스트")
    print("="*70)
    
    # ATPM 설정 (있는 경우)
    if atpm_rxn:
        atpm_rxn.lower_bound = atpm_value
        atpm_rxn.upper_bound = max(atpm_value, 1000)
        print(f"ATPM: {atpm_value} mmol ATP/gDCW/h")
    
    model.objective = biomass_rxn.id
    
    solution = model.optimize()
    
    print(f"\n최적화 상태: {solution.status}")
    
    if solution.status == 'optimal':
        biomass_flux = solution.fluxes.get(biomass_rxn.id, 0)
        print(f"Biomass flux: {biomass_flux:.6f} 1/h")
        
        if biomass_flux > 1e-6:
            print("\n[SUCCESS] 포도당 기반 성장 가능!")
            
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
            
            # Yield 계산
            for rxn in model.exchanges:
                if 'glc' in rxn.id.lower() or 'glucose' in rxn.id.lower():
                    glucose_uptake = abs(solution.fluxes.get(rxn.id, 0))
                    if glucose_uptake > 0:
                        yield_biomass = biomass_flux / glucose_uptake
                        print(f"\n성장 속도: {biomass_flux:.6f} 1/h")
                        print(f"Glucose 섭취: {glucose_uptake:.6f} mmol/gDCW/h")
                        print(f"Biomass yield: {yield_biomass:.6f} gDW/mmol glucose")
                    break
            
            return solution, True
        else:
            print("\n[FAIL] 포도당 기반 성장 불가 (Biomass flux = 0)")
            return solution, False
    else:
        print(f"\n[ERROR] 최적화 실패: {solution.status}")
        return solution, False

def compare_substrates(model, biomass_rxn, atpm_rxn=None):
    """다양한 탄소원 비교 테스트"""
    print("\n" + "="*70)
    print("다양한 탄소원 비교 테스트")
    print("="*70)
    
    substrates = {
        'Glucose': ['EX_glc_e', 'EX_glc__D_e', 'EX_glucose_e'],
        'Acetate': ['EX_ac_e'],
        'Pyruvate': ['EX_pyr_e', 'EX_pyruvate_e'],
        'Succinate': ['EX_succ_e', 'EX_succinate_e']
    }
    
    results = []
    
    for sub_name, ex_ids in substrates.items():
        print(f"\n[{sub_name}] 테스트 중...")
        
        # 모든 exchange 차단
        for rxn in model.exchanges:
            rxn.lower_bound = 0
            rxn.upper_bound = 0
        
        # 해당 탄소원 설정
        found = False
        for ex_id in ex_ids:
            try:
                ex_rxn = model.reactions.get_by_id(ex_id)
                ex_rxn.lower_bound = -10
                ex_rxn.upper_bound = 1000
                found = True
                print(f"  [OK] {ex_id} 설정")
                break
            except KeyError:
                continue
        
        if not found:
            print(f"  [SKIP] {sub_name} exchange 반응 없음")
            results.append({
                'Substrate': sub_name,
                'Status': 'No exchange',
                'Biomass_flux': 0
            })
            continue
        
        # 필수 영양소 설정
        for ex_id in ['EX_nh4_e', 'EX_h2o_e', 'EX_h_e', 'EX_pi_e', 'EX_so4_e',
                      'EX_k_e', 'EX_o2_e', 'EX_co2_e']:
            try:
                ex_rxn = model.reactions.get_by_id(ex_id)
                if ex_id == 'EX_co2_e':
                    ex_rxn.lower_bound = -1000
                    ex_rxn.upper_bound = 1000
                else:
                    ex_rxn.lower_bound = -1000
            except KeyError:
                pass
        
        # ATPM 설정
        if atpm_rxn:
            atpm_rxn.lower_bound = 0
            atpm_rxn.upper_bound = 1000
        
        # FBA 수행
        model.objective = biomass_rxn.id
        solution = model.optimize()
        
        if solution.status == 'optimal':
            biomass_flux = solution.fluxes.get(biomass_rxn.id, 0)
            status = "Growth" if biomass_flux > 1e-6 else "No growth"
            results.append({
                'Substrate': sub_name,
                'Status': status,
                'Biomass_flux': biomass_flux
            })
            print(f"  결과: {status} (Biomass flux: {biomass_flux:.6f})")
        else:
            results.append({
                'Substrate': sub_name,
                'Status': solution.status,
                'Biomass_flux': 0
            })
            print(f"  결과: {solution.status}")
    
    df_results = pd.DataFrame(results)
    
    print("\n" + "="*70)
    print("탄소원 비교 결과")
    print("="*70)
    print(df_results.to_string(index=False))
    
    return df_results

def main():
    print("="*70)
    print("모델 검증: 포도당 및 다양한 탄소원 테스트")
    print("="*70)
    
    # 모델 로드
    model = load_model("BaseModel.xml")
    
    # Biomass 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass: {biomass_rxn.id}")
    
    # ATP Maintenance 찾기
    atpm_rxn = find_atp_maintenance_reaction(model)
    if atpm_rxn:
        print(f"[OK] ATP Maintenance: {atpm_rxn.id}")
    
    # 포도당 medium 설정
    model = setup_glucose_medium(model, glucose_uptake=-10.0, oxygen_uptake=-100.0)
    
    # 포도당 기반 성장 테스트
    solution, success = test_growth_on_glucose(model, biomass_rxn, atpm_rxn, atpm_value=0)
    
    # 다양한 탄소원 비교
    df_comparison = compare_substrates(model, biomass_rxn, atpm_rxn)
    
    # 결과 저장
    df_comparison.to_csv('substrate_comparison.csv', index=False)
    print(f"\n[OK] 결과 저장: substrate_comparison.csv")
    
    print("\n" + "="*70)
    print("모델 검증 완료")
    print("="*70)
    
    if success:
        print("\n[결론] 모델 자체는 정상 작동합니다")
        print("  → 포도당 기반 성장 가능")
        print("  → Acetate 경로에 추가 문제가 있을 수 있습니다")
    else:
        print("\n[결론] 모델에 구조적 문제가 있을 수 있습니다")
        print("  → Biomass 반응이나 기본 대사 경로 점검 필요")

if __name__ == "__main__":
    main()
