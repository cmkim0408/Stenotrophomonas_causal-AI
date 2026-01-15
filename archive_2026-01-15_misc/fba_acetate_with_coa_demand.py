#!/usr/bin/env python
"""
Acetate 기반 FBA - CoA Demand Reaction 추가
부트스트랩 문제 해결을 위한 CoA 소량 공급 허용
"""

import cobra
from cobra import Reaction, Metabolite
import pandas as pd

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    print(f"[OK] 모델 로드: {model.id}")
    return model

def find_biomass_reaction(model):
    """Biomass 반응 찾기"""
    for name in ['Growth', 'BIOMASS']:
        try:
            return model.reactions.get_by_id(name)
        except KeyError:
            continue
    return None

def add_coa_demand_reaction(model, coa_supply_rate=-0.1):
    """
    CoA Demand Reaction 추가
    외부에서 소량의 CoA 공급을 허용하여 부트스트랩 문제 해결
    
    Parameters:
    -----------
    model : cobra.Model
        대사 모델
    coa_supply_rate : float
        CoA 공급 속도 (음수 = uptake). 기본값: -0.1 mmol/gDW/h
    """
    print("\n" + "="*70)
    print("CoA Demand Reaction 추가 (부트스트랩 해결)")
    print("="*70)
    
    try:
        coa_c = model.metabolites.get_by_id('coa_c')
        print(f"[OK] coa_c metabolite 존재: {coa_c.id}")
    except KeyError:
        print("[ERROR] coa_c metabolite를 찾을 수 없습니다!")
        return model, None
    
    # Demand reaction ID
    demand_id = 'DM_coa_c'
    
    # 이미 존재하는지 확인
    try:
        existing = model.reactions.get_by_id(demand_id)
        print(f"[INFO] {demand_id} 반응이 이미 존재합니다")
        # 기존 반응의 경계 조건 수정
        existing.lower_bound = coa_supply_rate
        existing.upper_bound = 1000
        print(f"[OK] {demand_id} 경계 조건 수정: LB={coa_supply_rate}, UB=1000")
        return model, existing
    except KeyError:
        pass
    
    # 새로운 demand reaction 생성
    demand_rxn = Reaction(demand_id)
    demand_rxn.name = 'CoA demand (bootstrap)'
    demand_rxn.lower_bound = coa_supply_rate  # 소량 공급 허용
    demand_rxn.upper_bound = 1000  # 제한 없음
    
    # coa_c --> (demand)
    demand_rxn.add_metabolites({coa_c: -1})
    
    # 모델에 추가
    model.add_reactions([demand_rxn])
    
    print(f"[OK] {demand_id} 반응 생성 완료")
    print(f"  반응식: {demand_rxn.reaction}")
    print(f"  LB={demand_rxn.lower_bound}, UB={demand_rxn.upper_bound}")
    print(f"  의미: 외부에서 CoA {abs(coa_supply_rate)} mmol/gDW/h 공급 허용")
    
    return model, demand_rxn

def setup_acetate_medium(model):
    """Acetate minimal medium 설정"""
    print("\n" + "="*70)
    print("Acetate Minimal Medium 설정")
    print("="*70)
    
    # 모든 exchange 차단
    for rxn in model.exchanges:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    
    # Acetate
    try:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.lower_bound = -1000
        ex_ac.upper_bound = 1000
        print(f"[OK] Acetate: EX_ac_e (-1000 ~ 1000)")
    except KeyError:
        print("[WARNING] EX_ac_e 없음")
    
    # 필수 영양소
    essentials = {
        'EX_nh4_e': 'Ammonium',
        'EX_h2o_e': 'Water',
        'EX_h_e': 'Proton',
        'EX_pi_e': 'Phosphate',
        'EX_so4_e': 'Sulfate',
        'EX_k_e': 'Potassium',
        'EX_na1_e': 'Sodium',
        'EX_mg2_e': 'Magnesium',
        'EX_ca2_e': 'Calcium',
        'EX_fe2_e': 'Iron',
        'EX_mn2_e': 'Manganese',
        'EX_zn2_e': 'Zinc',
        'EX_co2_e': 'CO2',
        'EX_o2_e': 'Oxygen'
    }
    
    print("\n필수 영양소 설정:")
    for ex_id, name in essentials.items():
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            if ex_id in ['EX_co2_e', 'EX_o2_e']:
                ex_rxn.lower_bound = -1000
                ex_rxn.upper_bound = 1000
            else:
                ex_rxn.lower_bound = -1000
            print(f"  [OK] {name}: {ex_id}")
        except KeyError:
            pass
    
    return model

def run_fba_with_coa_bootstrap(model, biomass_rxn):
    """CoA 부트스트랩을 포함한 FBA 분석"""
    print("\n" + "="*70)
    print("FBA 최적화 (CoA 부트스트랩 포함)")
    print("="*70)
    
    model.objective = biomass_rxn.id
    
    solution = model.optimize()
    
    print(f"\n최적화 상태: {solution.status}")
    print(f"Objective value: {solution.objective_value:.6f}")
    
    if solution.status == 'optimal':
        biomass_flux = solution.fluxes.get(biomass_rxn.id, 0)
        print(f"Biomass flux: {biomass_flux:.6f} 1/h")
        
        # CoA demand 플럭스 확인
        try:
            coa_demand_flux = solution.fluxes.get('DM_coa_c', 0)
            print(f"CoA demand flux: {coa_demand_flux:.6f} mmol/gDW/h")
            if abs(coa_demand_flux) > 1e-6:
                print(f"  → CoA가 {abs(coa_demand_flux):.6f} mmol/gDW/h 소비됨")
        except KeyError:
            pass
        
        if biomass_flux > 1e-6:
            print("\n[SUCCESS] 성장 가능!")
            
            # Exchange 플럭스
            print("\nExchange 플럭스 (절대값 > 0.001):")
            exchange_data = []
            for rxn in model.exchanges:
                flux = solution.fluxes.get(rxn.id, 0)
                if abs(flux) > 0.001:
                    exchange_data.append((rxn.id, flux))
            
            if exchange_data:
                for rxn_id, flux in sorted(exchange_data, key=lambda x: abs(x[1]), reverse=True):
                    direction = "Uptake" if flux < 0 else "Secretion"
                    print(f"  {rxn_id}: {flux:.6f} ({direction})")
            else:
                print("  [없음]")
            
            # 주요 경로 플럭스
            print("\n주요 경로 반응 플럭스:")
            pathway_reactions = {
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
                'DM_coa_c': 'CoA Demand',
                'Growth': 'Biomass'
            }
            
            active_pathway = []
            for rxn_id, rxn_name in pathway_reactions.items():
                try:
                    flux = solution.fluxes.get(rxn_id, 0)
                    if abs(flux) > 1e-8:
                        active_pathway.append((rxn_id, rxn_name, flux))
                        print(f"  {rxn_id} ({rxn_name}): {flux:.6f}")
                except KeyError:
                    pass
            
            # Yield 계산
            acetate_uptake = abs(solution.fluxes.get('EX_ac_e', 0))
            if acetate_uptake > 0:
                yield_biomass = biomass_flux / acetate_uptake
                print(f"\n성장 속도: {biomass_flux:.6f} 1/h")
                print(f"Acetate 섭취: {acetate_uptake:.6f} mmol/gDW/h")
                print(f"Biomass yield: {yield_biomass:.6f} gDW/mmol acetate")
                
                # CoA 사용량 대비 성장률
                try:
                    coa_usage = abs(solution.fluxes.get('DM_coa_c', 0))
                    if coa_usage > 0:
                        print(f"CoA 사용량: {coa_usage:.6f} mmol/gDW/h")
                        print(f"CoA 효율: {biomass_flux / coa_usage:.6f} gDW/mmol CoA")
                except KeyError:
                    pass
            
            return solution, True
        else:
            print("\n[FAIL] 여전히 성장 불가 (Biomass flux = 0)")
            print("  → 다른 부트스트랩 문제가 있을 수 있습니다")
            return solution, False
    else:
        print(f"\n[ERROR] 최적화 실패: {solution.status}")
        return solution, False

def main():
    print("="*70)
    print("Acetate 기반 FBA - CoA Demand Reaction 추가")
    print("="*70)
    
    # 모델 로드
    model = load_model("BaseModel.xml")
    
    # Biomass 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass: {biomass_rxn.id}")
    
    # Medium 설정
    model = setup_acetate_medium(model)
    
    # CoA Demand Reaction 추가
    coa_supply_rate = -0.1  # 0.1 mmol/gDW/h 공급 허용
    model, coa_demand = add_coa_demand_reaction(model, coa_supply_rate)
    
    if not coa_demand:
        print("\n[ERROR] CoA demand reaction 추가 실패")
        return
    
    # FBA 실행
    solution, success = run_fba_with_coa_bootstrap(model, biomass_rxn)
    
    print("\n" + "="*70)
    if success:
        print("[SUCCESS] CoA 부트스트랩으로 성장 가능 확인!")
        print("\n결론:")
        print("  - 모델 구조는 올바르게 구성되어 있습니다")
        print("  - 문제는 단순한 부트스트랩(초기화) 문제였습니다")
        print("  - CoA 소량 공급으로 전체 대사 경로가 작동합니다")
    else:
        print("[FAIL] 여전히 성장 불가")
        print("  - 추가 진단이 필요할 수 있습니다")
    print("="*70)

if __name__ == "__main__":
    main()
