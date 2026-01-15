#!/usr/bin/env python
"""
최종 Glyoxylate shunt 분석 및 정리
"""

import cobra
import pandas as pd
from pathlib import Path

def load_model(model_path):
    model = cobra.io.read_sbml_model(model_path)
    return model

def setup_acetate_medium(model):
    """Acetate 미디어 설정"""
    # 모든 exchange 차단
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    # Acetate 허용
    ex_ac = model.reactions.get_by_id('EX_ac_e')
    ex_ac.upper_bound = 1000
    ex_ac.lower_bound = -1000
    
    # 필수 무기염
    essential = {
        'EX_nh4_e': (-1000, 1000),
        'EX_h2o_e': (-1000, 1000),
        'EX_h_e': (-1000, 1000),
        'EX_pi_e': (-1000, 1000),
        'EX_so4_e': (-1000, 1000),
        'EX_k_e': (-1000, 1000),
        'EX_na1_e': (-1000, 1000),
        'EX_mg2_e': (-1000, 1000),
        'EX_ca2_e': (-1000, 1000),
        'EX_fe2_e': (-1000, 1000),
        'EX_mn2_e': (-1000, 1000),
        'EX_zn2_e': (-1000, 1000),
        'EX_co2_e': (-1000, 1000),
        'EX_o2_e': (-1000, 1000),
    }
    
    for ex_id, (lb, ub) in essential.items():
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            ex_rxn.upper_bound = ub
            ex_rxn.lower_bound = lb
        except KeyError:
            pass
    
    return model

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    
    model = load_model(str(model_path))
    model.objective = 'Growth'
    
    atpm_rxn = model.reactions.get_by_id('ATPM')
    biomass_rxn = model.reactions.get_by_id('Growth')
    
    model = setup_acetate_medium(model)
    
    print("="*80)
    print("신규 모델 Glyoxylate Shunt 분석")
    print("="*80)
    
    print("\n[레퍼런스 모델 결과 (Stenotrophomonas/fba_flux_gradient_acid.csv)]")
    print("  ATPM=0:   ICL=0.367, MALS=0.367, ICDHx=0.0   -> Glyoxylate shunt만")
    print("  ATPM=5:   ICL=0.074, MALS=0.074, ICDHx=0.793 -> 혼합")
    print("  ATPM=10:  ICL=0.0,   MALS=0.0,   ICDHx=1.0   -> TCA cycle만")
    
    print("\n[신규 모델 결과]")
    results = []
    
    for atpm_val in [0, 5, 10]:
        atpm_rxn.lower_bound = atpm_val
        atpm_rxn.upper_bound = 1000.0
        
        solution = model.optimize()
        
        growth = solution.objective_value
        icl_flux = solution.fluxes.get('ICL', 0.0)
        mals_flux = solution.fluxes.get('MALS', 0.0)
        icdhx_flux = solution.fluxes.get('ICDHx', 0.0)
        cs_flux = solution.fluxes.get('CS', 0.0)
        
        print(f"\n  ATPM={atpm_val}:")
        print(f"    성장률: {growth:.6f}")
        print(f"    ICL: {icl_flux:.6f}")
        print(f"    MALS: {mals_flux:.6f}")
        print(f"    ICDHx: {icdhx_flux:.6f}")
        print(f"    CS: {cs_flux:.6f}")
        
        if abs(icl_flux) > 1e-6 and abs(mals_flux) > 1e-6:
            if abs(icdhx_flux) < 1e-6:
                print(f"    -> Glyoxylate shunt만 사용")
            else:
                print(f"    -> Glyoxylate shunt + TCA cycle 혼합")
        elif abs(icdhx_flux) > 1e-6:
            print(f"    -> TCA cycle만 사용")
        
        results.append({
            'ATPM': atpm_val,
            'Growth': growth,
            'ICL': icl_flux,
            'MALS': mals_flux,
            'ICDHx': icdhx_flux,
            'CS': cs_flux
        })
    
    print("\n" + "="*80)
    print("결론")
    print("="*80)
    
    df = pd.DataFrame(results)
    print("\n요약:")
    print(df.to_string(index=False))
    
    # 레퍼런스와 비교
    print("\n[레퍼런스 모델과 비교]")
    
    if results[0]['Growth'] > 1e-6:
        print("\n  ATPM=0:")
        if abs(results[0]['ICL']) > 1e-6 and abs(results[0]['ICDHx']) < 1e-6:
            print("    [OK] 레퍼런스와 유사: Glyoxylate shunt만 사용")
        else:
            print("    [차이] 레퍼런스와 다름")
    
    if results[2]['Growth'] > 1e-6:
        print("\n  ATPM=10:")
        if abs(results[2]['ICL']) < 1e-6 and abs(results[2]['ICDHx']) > 1e-6:
            print("    [OK] 레퍼런스와 유사: TCA cycle만 사용")
        elif abs(results[2]['ICL']) > 1e-6:
            print("    [차이] 레퍼런스는 ICL=0이지만 신규 모델은 ICL>0 (혼합 사용)")

if __name__ == "__main__":
    main()
