#!/usr/bin/env python
"""
무제한 영양소 상태에서 사용되는 구성 요소 분석
실제로 작동하는 조건에서 무엇이 사용되는지 확인
"""

import cobra
import pandas as pd

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    return model

def find_biomass_reaction(model):
    for name in ['Growth', 'BIOMASS']:
        try:
            return model.reactions.get_by_id(name)
        except KeyError:
            continue
    return None

def analyze_unlimited_growth(model, biomass_rxn):
    """무제한 영양소에서 생장 분석"""
    print("="*70)
    print("무제한 영양소 상태 분석")
    print("="*70)
    
    # 모든 exchange 무제한
    for rxn in model.exchanges:
        rxn.lower_bound = -1000
        rxn.upper_bound = 1000
    
    model.objective = biomass_rxn.id
    solution = model.optimize()
    
    if solution.status == 'optimal':
        biomass_flux = solution.fluxes.get(biomass_rxn.id, 0)
        print(f"\n[OK] 무제한 영양소에서 생장 가능")
        print(f"  Biomass flux: {biomass_flux:.6f} 1/h")
        
        # 사용된 exchange 확인
        print("\n[사용된 Exchange (Uptake)]:")
        uptake_exchanges = []
        for rxn in model.exchanges:
            flux = solution.fluxes.get(rxn.id, 0)
            if flux < -0.001:  # Uptake
                uptake_exchanges.append({
                    'Exchange': rxn.id,
                    'Flux': flux,
                    'Uptake_Rate': abs(flux)
                })
                print(f"  {rxn.id}: {flux:.6f} mmol/gDCW/h")
        
        # 생성된 exchange 확인
        print("\n[생성된 Exchange (Secretion)]:")
        secretion_exchanges = []
        for rxn in model.exchanges:
            flux = solution.fluxes.get(rxn.id, 0)
            if flux > 0.001:  # Secretion
                secretion_exchanges.append({
                    'Exchange': rxn.id,
                    'Flux': flux,
                    'Secretion_Rate': flux
                })
                print(f"  {rxn.id}: {flux:.6f} mmol/gDCW/h")
        
        # 결과 저장
        if uptake_exchanges:
            df_uptake = pd.DataFrame(uptake_exchanges)
            df_uptake = df_uptake.sort_values('Uptake_Rate', ascending=False)
            df_uptake.to_csv('unlimited_uptake_exchanges.csv', index=False)
            print(f"\n[OK] Uptake exchange 저장: unlimited_uptake_exchanges.csv")
        
        if secretion_exchanges:
            df_secretion = pd.DataFrame(secretion_exchanges)
            df_secretion = df_secretion.sort_values('Secretion_Rate', ascending=False)
            df_secretion.to_csv('unlimited_secretion_exchanges.csv', index=False)
            print(f"[OK] Secretion exchange 저장: unlimited_secretion_exchanges.csv")
        
        return solution, uptake_exchanges, secretion_exchanges
    else:
        print(f"\n[FAIL] 최적화 실패: {solution.status}")
        return solution, [], []

def test_minimal_nutrients(model, biomass_rxn, uptake_list):
    """최소 영양소만으로 테스트"""
    print("\n" + "="*70)
    print("최소 영양소만으로 테스트")
    print("="*70)
    
    # 모든 exchange 초기화
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    # 무제한 상태에서 사용된 영양소만 허용
    print("\n[허용된 영양소]:")
    for exchange_info in uptake_list[:20]:  # 상위 20개
        ex_id = exchange_info['Exchange']
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            if ex_id in ['EX_co2_e', 'EX_o2_e']:
                ex_rxn.lower_bound = -1000
                ex_rxn.upper_bound = 1000
            else:
                ex_rxn.lower_bound = -1000
            print(f"  {ex_id}: {exchange_info['Uptake_Rate']:.6f} mmol/gDCW/h")
        except KeyError:
            pass
    
    model.objective = biomass_rxn.id
    solution = model.optimize()
    
    print(f"\n최적화 상태: {solution.status}")
    
    if solution.status == 'optimal':
        biomass_flux = solution.fluxes.get(biomass_rxn.id, 0)
        if biomass_flux > 1e-6:
            print(f"[SUCCESS] 최소 영양소로 생장 가능!")
            print(f"  Biomass flux: {biomass_flux:.6f} 1/h")
            return True, solution
        else:
            print(f"[FAIL] Biomass flux: {biomass_flux:.6f} 1/h")
            return False, solution
    else:
        print(f"[FAIL] 최적화 실패: {solution.status}")
        return False, solution

def main():
    model = load_model("BaseModel.xml")
    
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass: {biomass_rxn.id}")
    
    # 무제한 영양소 분석
    solution, uptake_list, secretion_list = analyze_unlimited_growth(model, biomass_rxn)
    
    if solution.status == 'optimal' and solution.objective_value > 1e-6:
        # 최소 영양소 테스트
        can_grow, min_solution = test_minimal_nutrients(model, biomass_rxn, uptake_list)
        
        print("\n" + "="*70)
        print("결과 요약")
        print("="*70)
        
        if can_grow:
            print("\n[SUCCESS] 최소 영양소로 생장 가능!")
            print("  → 무제한 영양소에서 사용된 구성 요소만으로 생장 가능")
            print("\n[최소 영양소 목록]:")
            for exchange_info in uptake_list[:20]:
                print(f"  - {exchange_info['Exchange']}: {exchange_info['Uptake_Rate']:.6f} mmol/gDCW/h")
        else:
            print("\n[FAIL] 최소 영양소로 생장 불가능")
            print("  → 추가 구성 요소 필요")
    
    print("\n" + "="*70)

if __name__ == "__main__":
    main()
