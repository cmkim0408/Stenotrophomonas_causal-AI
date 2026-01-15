#!/usr/bin/env python
"""
G6P 직접 공급 + 부트스트랩 테스트
"""

import cobra

def main():
    model = cobra.io.read_sbml_model("BaseModel.xml")
    
    # 모든 exchange 초기화
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    # 필수 영양소
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
    
    # G6P 직접 공급
    g6p = model.metabolites.get_by_id('g6p_c')
    dm_g6p = cobra.Reaction('DM_g6p')
    dm_g6p.add_metabolites({g6p: -1})
    dm_g6p.lower_bound = -100
    model.add_reactions([dm_g6p])
    
    # 부트스트랩 추가
    atp = model.metabolites.get_by_id('atp_c')
    dm_atp = cobra.Reaction('DM_atp')
    dm_atp.add_metabolites({atp: -1})
    dm_atp.lower_bound = -0.1
    model.add_reactions([dm_atp])
    
    nad = model.metabolites.get_by_id('nad_c')
    dm_nad = cobra.Reaction('DM_nad')
    dm_nad.add_metabolites({nad: -1})
    dm_nad.lower_bound = -0.1
    model.add_reactions([dm_nad])
    
    coa = model.metabolites.get_by_id('coa_c')
    dm_coa = cobra.Reaction('DM_coa')
    dm_coa.add_metabolites({coa: -1})
    dm_coa.lower_bound = -0.01
    model.add_reactions([dm_coa])
    
    print("="*70)
    print("G6P 직접 공급 + ATP/NAD+/CoA 부트스트랩")
    print("="*70)
    print("\n설정:")
    print("  G6P: -100 mmol/gDCW/h")
    print("  ATP: -0.1 mmol/gDCW/h")
    print("  NAD+: -0.1 mmol/gDCW/h")
    print("  CoA: -0.01 mmol/gDCW/h")
    
    model.objective = 'Growth'
    solution = model.optimize()
    
    print(f"\n최적화 상태: {solution.status}")
    
    if solution.status == 'optimal':
        biomass = solution.objective_value
        print(f"Biomass flux: {biomass:.6f} 1/h")
        
        if biomass > 1e-6:
            print("\n[SUCCESS] 생장 가능!")
            
            # ETC 경로 확인
            etc_rxns = ['NADH16pp', 'ATPS4rpp']
            print("\nETC 경로:")
            for rxn_id in etc_rxns:
                try:
                    flux = solution.fluxes.get(rxn_id, 0)
                    if abs(flux) > 1e-8:
                        print(f"  {rxn_id}: {flux:.6f}")
                except:
                    pass
        else:
            print("\n[FAIL] 생장 불가능 (flux = 0)")
    else:
        print("\n[FAIL] 최적화 실패")

if __name__ == "__main__":
    main()
