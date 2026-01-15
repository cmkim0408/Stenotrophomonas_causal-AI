#!/usr/bin/env python
"""
ATPM 값을 낮춰서 Glyoxylate shunt 활성화 확인
레퍼런스 모델에서는 ATPM=0일 때 ICL/MALS가 활성화됨
"""

import cobra
import pandas as pd
from pathlib import Path

def load_model(model_path):
    print(f"모델 로드 중: {model_path}")
    model = cobra.io.read_sbml_model(model_path)
    print(f"[OK] 모델 로드 완료: {model.id}")
    return model

def find_biomass_reaction(model):
    """Biomass 반응 찾기"""
    for name in ['Growth', 'BIOMASS', 'BIOMASS_Ecoli_core', 'R_BIOMASS']:
        try:
            return model.reactions.get_by_id(name)
        except KeyError:
            continue
    
    for rxn in model.reactions:
        if 'biomass' in rxn.id.lower() or 'growth' in rxn.id.lower():
            return rxn
    
    return None

def find_atpm_reaction(model):
    """ATPM 반응 찾기"""
    atpm_keywords = ['ATPM', 'ATPM', 'ATP_maintenance', 'atp_maintenance']
    for rxn in model.reactions:
        if any(kw.lower() in rxn.id.lower() for kw in atpm_keywords):
            return rxn
    return None

def setup_acetate_medium(model):
    """Acetate 미디어 설정"""
    # 모든 exchange 차단
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    # Acetate 허용
    try:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.upper_bound = 1000
        ex_ac.lower_bound = -1000
    except KeyError:
        for rxn in model.exchanges:
            if 'ac_e' in rxn.id.lower() and 'EX' in rxn.id:
                rxn.upper_bound = 1000
                rxn.lower_bound = -1000
                break
    
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

def test_different_atpm_values(model, biomass_rxn, atpm_rxn, atpm_values=[0, 5, 10]):
    """다양한 ATPM 값에서 테스트"""
    print("\n" + "="*80)
    print("다양한 ATPM 값에서 테스트 (레퍼런스 모델과 비교)")
    print("="*80)
    
    results = []
    
    for atpm_val in atpm_values:
        print(f"\n[ATPM = {atpm_val}]")
        
        # ATPM 설정
        if atpm_rxn:
            atpm_rxn.lower_bound = atpm_val
            atpm_rxn.upper_bound = 1000.0
        
        # FBA 수행
        solution = model.optimize()
        
        if solution.status == 'optimal':
            growth = solution.objective_value
            
            # 주요 반응 플럭스
            icl_flux = solution.fluxes.get('ICL', 0.0)
            mals_flux = solution.fluxes.get('MALS', 0.0)
            icdhx_flux = solution.fluxes.get('ICDHx', 0.0)
            cs_flux = solution.fluxes.get('CS', 0.0)
            atpm_actual = solution.fluxes.get(atpm_rxn.id, 0.0) if atpm_rxn else 0.0
            
            print(f"  성장률: {growth:.6f}")
            print(f"  ATPM 플럭스: {atpm_actual:.6f}")
            print(f"  ICL: {icl_flux:.6f}")
            print(f"  MALS: {mals_flux:.6f}")
            print(f"  ICDHx: {icdhx_flux:.6f}")
            print(f"  CS: {cs_flux:.6f}")
            
            if abs(icl_flux) > 1e-6 and abs(mals_flux) > 1e-6:
                print(f"  -> Glyoxylate shunt 활성!")
            elif abs(icdhx_flux) > 1e-6:
                print(f"  -> TCA cycle (ICDHx) 활성")
            
            results.append({
                'ATPM': atpm_val,
                'Growth': growth,
                'ICL': icl_flux,
                'MALS': mals_flux,
                'ICDHx': icdhx_flux,
                'CS': cs_flux
            })
        else:
            print(f"  [ERROR] FBA 실패: {solution.status}")
    
    # 레퍼런스 모델 결과와 비교
    print("\n" + "="*80)
    print("레퍼런스 모델 결과 비교")
    print("="*80)
    print("\n레퍼런스 모델 (Stenotrophomonas/fba_flux_gradient_acid.csv):")
    print("  ATPM=0:  ICL=0.367, MALS=0.367, ICDHx=0.0   -> Glyoxylate shunt 활성")
    print("  ATPM=5:  ICL=0.074, MALS=0.074, ICDHx=0.793 -> 혼합")
    print("  ATPM=10: ICL=0.0,   MALS=0.0,   ICDHx=1.0   -> TCA cycle 활성")
    
    print("\n신규 모델 결과:")
    df_results = pd.DataFrame(results)
    print(df_results.to_string(index=False))
    
    return results

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    
    model = load_model(str(model_path))
    
    # Biomass 반응 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("\n[ERROR] Biomass 반응을 찾을 수 없습니다.")
        return
    
    print(f"\n[Biomass 반응] {biomass_rxn.id}")
    model.objective = biomass_rxn.id
    
    # ATPM 찾기
    atpm_rxn = find_atpm_reaction(model)
    if atpm_rxn:
        print(f"\n[ATPM 반응] {atpm_rxn.id}")
        print(f"  현재 하한: {atpm_rxn.lower_bound}")
    else:
        print("\n[WARNING] ATPM 반응을 찾을 수 없습니다.")
    
    # 미디어 설정
    model = setup_acetate_medium(model)
    
    # 다양한 ATPM 값에서 테스트
    results = test_different_atpm_values(model, biomass_rxn, atpm_rxn, atpm_values=[0, 5, 10])
    
    print("\n" + "="*80)
    print("결론")
    print("="*80)
    if results:
        atpm_0_result = [r for r in results if r['ATPM'] == 0]
        if atpm_0_result and (abs(atpm_0_result[0]['ICL']) > 1e-6 or abs(atpm_0_result[0]['MALS']) > 1e-6):
            print("\n[OK] ATPM=0일 때 Glyoxylate shunt가 활성화됩니다!")
        else:
            print("\n[문제] ATPM=0일 때도 Glyoxylate shunt가 활성화되지 않습니다.")
            print("  -> 모델에 문제가 있을 수 있습니다.")

if __name__ == "__main__":
    main()
