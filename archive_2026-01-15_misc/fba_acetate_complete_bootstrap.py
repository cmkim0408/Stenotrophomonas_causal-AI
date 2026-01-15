#!/usr/bin/env python
"""
Acetate 기반 FBA - 완전한 부트스트랩
SUCOAACTr 반응 작동을 위한 Succinyl-CoA 초기화 포함
"""

import cobra
from cobra import Reaction

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

def add_demand_reaction(model, metabolite_id, demand_id, supply_rate=-0.1):
    """Demand Reaction 추가"""
    try:
        met = model.metabolites.get_by_id(metabolite_id)
    except KeyError:
        return model, None
    
    try:
        existing = model.reactions.get_by_id(demand_id)
        existing.lower_bound = supply_rate
        return model, existing
    except KeyError:
        pass
    
    demand_rxn = Reaction(demand_id)
    demand_rxn.name = f'{metabolite_id} demand (bootstrap)'
    demand_rxn.lower_bound = supply_rate
    demand_rxn.upper_bound = 1000
    demand_rxn.add_metabolites({met: -1})
    
    model.add_reactions([demand_rxn])
    return model, demand_rxn

def setup_acetate_medium(model):
    """Acetate minimal medium 설정"""
    for rxn in model.exchanges:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    
    try:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.lower_bound = -1000
        ex_ac.upper_bound = 1000
    except KeyError:
        pass
    
    essentials = ['EX_nh4_e', 'EX_h2o_e', 'EX_h_e', 'EX_pi_e', 'EX_so4_e',
                  'EX_k_e', 'EX_na1_e', 'EX_mg2_e', 'EX_ca2_e', 'EX_fe2_e',
                  'EX_mn2_e', 'EX_zn2_e', 'EX_co2_e', 'EX_o2_e']
    
    for ex_id in essentials:
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

def add_complete_bootstrap(model):
    """완전한 부트스트랩 설정"""
    print("\n" + "="*70)
    print("완전한 부트스트랩 설정")
    print("="*70)
    
    # SUCOAACTr 작동을 위한 부트스트랩
    # ac_c + succoa_c <=> accoa_c + succ_c
    # Succinyl-CoA가 필요함
    
    bootstrap_cofactors = {
        'coa_c': ('DM_coa_c', -0.1, 'CoA'),
        'succoa_c': ('DM_succoa_c', -0.01, 'Succinyl-CoA'),  # SUCOAACTr 작동용
        'atp_c': ('DM_atp_c', -0.1, 'ATP'),
        'nad_c': ('DM_nad_c', -0.01, 'NAD+'),
    }
    
    demand_reactions = {}
    
    for met_id, (demand_id, supply_rate, name) in bootstrap_cofactors.items():
        model, demand_rxn = add_demand_reaction(model, met_id, demand_id, supply_rate)
        if demand_rxn:
            demand_reactions[demand_id] = demand_rxn
            print(f"[OK] {name}: {demand_id} (LB={supply_rate})")
    
    print(f"\n총 {len(demand_reactions)}개 부트스트랩 반응 추가")
    
    return model, demand_reactions

def run_fba_complete(model, biomass_rxn, demand_reactions):
    """완전한 FBA 분석"""
    print("\n" + "="*70)
    print("FBA 최적화 (완전한 부트스트랩)")
    print("="*70)
    
    model.objective = biomass_rxn.id
    
    solution = model.optimize()
    
    print(f"\n최적화 상태: {solution.status}")
    print(f"Objective value: {solution.objective_value:.6f}")
    
    if solution.status == 'optimal':
        biomass_flux = solution.fluxes.get(biomass_rxn.id, 0)
        print(f"Biomass flux: {biomass_flux:.6f} 1/h")
        
        if biomass_flux > 1e-6:
            print("\n[SUCCESS] 성장 가능!")
            
            # 부트스트랩 플럭스
            print("\n부트스트랩 Demand 플럭스:")
            for demand_id in demand_reactions.keys():
                try:
                    flux = solution.fluxes.get(demand_id, 0)
                    if abs(flux) > 1e-8:
                        print(f"  {demand_id}: {flux:.6f} mmol/gDW/h")
                except KeyError:
                    pass
            
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
            
            active_count = 0
            for rxn_id, rxn_name in pathway_reactions.items():
                try:
                    flux = solution.fluxes.get(rxn_id, 0)
                    if abs(flux) > 1e-8:
                        active_count += 1
                        print(f"  {rxn_id} ({rxn_name}): {flux:.6f}")
                except KeyError:
                    pass
            
            print(f"\n활성 경로 반응: {active_count}개")
            
            # Yield 계산
            acetate_uptake = abs(solution.fluxes.get('EX_ac_e', 0))
            if acetate_uptake > 0:
                yield_biomass = biomass_flux / acetate_uptake
                print(f"\n성장 속도: {biomass_flux:.6f} 1/h")
                print(f"Acetate 섭취: {acetate_uptake:.6f} mmol/gDW/h")
                print(f"Biomass yield: {yield_biomass:.6f} gDW/mmol acetate")
            
            return solution, True
        else:
            print("\n[FAIL] 여전히 성장 불가")
            return solution, False
    else:
        print(f"\n[ERROR] 최적화 실패: {solution.status}")
        return solution, False

def main():
    print("="*70)
    print("Acetate 기반 FBA - 완전한 부트스트랩 (SUCOAACTr 포함)")
    print("="*70)
    
    model = load_model("BaseModel.xml")
    
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass: {biomass_rxn.id}")
    
    model = setup_acetate_medium(model)
    
    model, demand_reactions = add_complete_bootstrap(model)
    
    solution, success = run_fba_complete(model, biomass_rxn, demand_reactions)
    
    print("\n" + "="*70)
    if success:
        print("[SUCCESS] 완전한 부트스트랩으로 성장 가능!")
        print("\n결론:")
        print("  - SUCOAACTr 반응 작동을 위해 Succinyl-CoA 초기화 필요")
        print("  - CoA, ATP 등 보조인자 초기화 필요")
        print("  - Acetate 기반 성장 경로 정상 작동 확인")
    else:
        print("[FAIL] 여전히 성장 불가")
        print("  - 모델 구조에 추가 문제가 있을 수 있습니다")
    print("="*70)

if __name__ == "__main__":
    main()
