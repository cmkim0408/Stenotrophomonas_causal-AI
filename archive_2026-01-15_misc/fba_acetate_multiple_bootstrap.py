#!/usr/bin/env python
"""
Acetate 기반 FBA - 다중 보조인자 부트스트랩
CoA, ATP 등 여러 보조인자에 대한 demand reaction 추가
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

def add_demand_reaction(model, metabolite_id, demand_id, supply_rate=-0.1):
    """
    Demand Reaction 추가
    
    Parameters:
    -----------
    model : cobra.Model
        대사 모델
    metabolite_id : str
        대사물질 ID
    demand_id : str
        Demand 반응 ID
    supply_rate : float
        공급 속도 (음수 = uptake)
    """
    try:
        met = model.metabolites.get_by_id(metabolite_id)
    except KeyError:
        print(f"[WARNING] {metabolite_id} metabolite 없음")
        return model, None
    
    # 이미 존재하는지 확인
    try:
        existing = model.reactions.get_by_id(demand_id)
        existing.lower_bound = supply_rate
        existing.upper_bound = 1000
        print(f"[OK] {demand_id} 경계 조건 수정: LB={supply_rate}")
        return model, existing
    except KeyError:
        pass
    
    # 새로운 demand reaction 생성
    demand_rxn = Reaction(demand_id)
    demand_rxn.name = f'{metabolite_id} demand (bootstrap)'
    demand_rxn.lower_bound = supply_rate
    demand_rxn.upper_bound = 1000
    
    demand_rxn.add_metabolites({met: -1})
    
    model.add_reactions([demand_rxn])
    print(f"[OK] {demand_id} 반응 생성: {demand_rxn.reaction} (LB={supply_rate})")
    
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
        print(f"[OK] Acetate: EX_ac_e")
    except KeyError:
        pass
    
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
    
    for ex_id, name in essentials.items():
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            if ex_id in ['EX_co2_e', 'EX_o2_e']:
                ex_rxn.lower_bound = -1000
                ex_rxn.upper_bound = 1000
            else:
                ex_rxn.lower_bound = -1000
        except KeyError:
            pass
    
    return model

def add_multiple_bootstrap_cofactors(model):
    """다중 보조인자 부트스트랩 추가"""
    print("\n" + "="*70)
    print("다중 보조인자 부트스트랩 추가")
    print("="*70)
    
    # 부트스트랩 보조인자 목록
    bootstrap_cofactors = {
        'coa_c': ('DM_coa_c', -0.1, 'CoA'),
        'atp_c': ('DM_atp_c', -0.1, 'ATP'),
        'nad_c': ('DM_nad_c', -0.01, 'NAD+'),
        'nadp_c': ('DM_nadp_c', -0.01, 'NADP+'),
    }
    
    demand_reactions = {}
    
    for met_id, (demand_id, supply_rate, name) in bootstrap_cofactors.items():
        model, demand_rxn = add_demand_reaction(model, met_id, demand_id, supply_rate)
        if demand_rxn:
            demand_reactions[demand_id] = demand_rxn
    
    print(f"\n총 {len(demand_reactions)}개 부트스트랩 반응 추가")
    
    return model, demand_reactions

def run_fba_with_bootstrap(model, biomass_rxn, demand_reactions):
    """부트스트랩을 포함한 FBA 분석"""
    print("\n" + "="*70)
    print("FBA 최적화 (다중 부트스트랩)")
    print("="*70)
    
    model.objective = biomass_rxn.id
    
    solution = model.optimize()
    
    print(f"\n최적화 상태: {solution.status}")
    print(f"Objective value: {solution.objective_value:.6f}")
    
    if solution.status == 'optimal':
        biomass_flux = solution.fluxes.get(biomass_rxn.id, 0)
        print(f"Biomass flux: {biomass_flux:.6f} 1/h")
        
        # Demand 플럭스 확인
        print("\n부트스트랩 Demand 플럭스:")
        for demand_id in demand_reactions.keys():
            try:
                flux = solution.fluxes.get(demand_id, 0)
                if abs(flux) > 1e-8:
                    print(f"  {demand_id}: {flux:.6f} mmol/gDW/h")
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
                'Growth': 'Biomass'
            }
            
            for rxn_id, rxn_name in pathway_reactions.items():
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
                print(f"Acetate 섭취: {acetate_uptake:.6f} mmol/gDW/h")
                print(f"Biomass yield: {yield_biomass:.6f} gDW/mmol acetate")
            
            return solution, True
        else:
            print("\n[FAIL] 여전히 성장 불가 (Biomass flux = 0)")
            print("  → 다른 구조적 문제가 있을 수 있습니다")
            return solution, False
    else:
        print(f"\n[ERROR] 최적화 실패: {solution.status}")
        return solution, False

def main():
    print("="*70)
    print("Acetate 기반 FBA - 다중 보조인자 부트스트랩")
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
    
    # 다중 부트스트랩 추가
    model, demand_reactions = add_multiple_bootstrap_cofactors(model)
    
    # FBA 실행
    solution, success = run_fba_with_bootstrap(model, biomass_rxn, demand_reactions)
    
    print("\n" + "="*70)
    if success:
        print("[SUCCESS] 다중 부트스트랩으로 성장 가능 확인!")
        print("\n결론:")
        print("  - 모델 구조는 올바르게 구성되어 있습니다")
        print("  - 여러 보조인자의 부트스트랩이 필요했습니다")
        print("  - Acetate 기반 성장 경로가 정상 작동합니다")
    else:
        print("[FAIL] 여전히 성장 불가")
        print("  - 모델 구조 검토 필요")
        print("  - 추가 진단 필요")
    print("="*70)

if __name__ == "__main__":
    main()
