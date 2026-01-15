#!/usr/bin/env python
"""
Acetate 기반 FBA - ATP Maintenance (ATPM) Sweep 분석
논문 방법론에 따른 ATPM 파라미터 변화에 따른 대사 플럭스 분석
"""

import cobra
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def load_model(model_path="BaseModel.xml"):
    """모델 로드"""
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

def find_atp_maintenance_reaction(model):
    """ATP Maintenance 반응 찾기"""
    # 일반적인 ATP maintenance 반응 이름들
    atpm_candidates = ['ATPM', 'ATPM_c', 'ATPS4rpp', 'ATPS4r', 'ATP_maintenance']
    
    for atpm_id in atpm_candidates:
        try:
            atpm_rxn = model.reactions.get_by_id(atpm_id)
            return atpm_rxn
        except KeyError:
            continue
    
    # ATP 소비 반응 찾기 (atp_c --> adp_c + pi_c + h_c)
    try:
        atp_c = model.metabolites.get_by_id('atp_c')
        adp_c = model.metabolites.get_by_id('adp_c')
        pi_c = model.metabolites.get_by_id('pi_c')
        
        # ATP를 소비하고 ADP를 생성하는 반응 찾기
        for rxn in atp_c.reactions:
            if (atp_c in rxn.reactants and 
                adp_c in rxn.products and 
                'maintenance' in rxn.name.lower()):
                return rxn
    except KeyError:
        pass
    
    # 없으면 ATP maintenance 반응 생성
    print("[INFO] ATP Maintenance 반응을 찾을 수 없습니다. 생성 중...")
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
        print(f"[OK] ATPM 반응 생성: {atpm_rxn.reaction}")
        return atpm_rxn
    except KeyError as e:
        print(f"[ERROR] ATPM 반응 생성 실패: {e}")
        return None

def setup_acetate_medium(model, acetate_uptake=-19.0, oxygen_uptake=-100.0):
    """
    Acetate minimal medium 설정
    
    Parameters:
    -----------
    acetate_uptake : float
        Acetate 섭취 속도 (mmol/gDCW/h). 기본값: -19 (논문 값)
    oxygen_uptake : float
        Oxygen 섭취 속도 (mmol/gDCW/h). 기본값: -100 (논문 값)
    """
    print("\n" + "="*70)
    print("Acetate Minimal Medium 설정")
    print("="*70)
    
    # 모든 exchange 차단
    for rxn in model.exchanges:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    
    # Acetate 설정
    try:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.lower_bound = acetate_uptake
        ex_ac.upper_bound = acetate_uptake  # 고정
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
        'EX_co2_e': 'CO2'
    }
    
    for ex_id, name in essentials.items():
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            if ex_id == 'EX_co2_e':
                ex_rxn.lower_bound = -1000  # 생성 허용
                ex_rxn.upper_bound = 1000
            else:
                ex_rxn.lower_bound = -1000  # Uptake 허용
            print(f"  [OK] {name}: {ex_id}")
        except KeyError:
            pass
    
    return model

def run_atpm_sweep(model, biomass_rxn, atpm_rxn, atpm_range=None):
    """
    ATPM sweep 수행
    
    Parameters:
    -----------
    model : cobra.Model
        대사 모델
    biomass_rxn : cobra.Reaction
        Biomass 반응
    atpm_rxn : cobra.Reaction
        ATP Maintenance 반응
    atpm_range : list or None
        ATPM 값 리스트. None이면 기본 범위 사용
    """
    print("\n" + "="*70)
    print("ATPM Sweep 분석")
    print("="*70)
    
    if atpm_range is None:
        # 논문과 유사한 범위: 0, 20, 40, 60, 80, 100
        atpm_range = [0, 20, 40, 60, 80, 100, 110, 120]
    
    model.objective = biomass_rxn.id
    
    results = []
    
    print(f"\nATPM 범위: {min(atpm_range)} ~ {max(atpm_range)} mmol ATP/gDCW/h")
    print(f"\n{'ATPM':<10} {'Biomass flux':<15} {'CO2 production':<18} {'Status':<10}")
    print("-" * 60)
    
    for atpm_value in atpm_range:
        # ATPM 설정
        atpm_rxn.lower_bound = atpm_value
        atpm_rxn.upper_bound = max(atpm_value, 1000)  # 최소한 현재 값 이상
        
        # FBA 수행
        solution = model.optimize()
        
        if solution.status == 'optimal':
            biomass_flux = solution.fluxes.get(biomass_rxn.id, 0)
            
            # CO2 생성량
            try:
                co2_exchange = model.reactions.get_by_id('EX_co2_e')
                co2_production = solution.fluxes.get('EX_co2_e', 0)
            except KeyError:
                co2_production = 0
            
            status = "OK" if biomass_flux > 1e-6 else "No growth"
            
            results.append({
                'ATPM': atpm_value,
                'Biomass_flux': biomass_flux,
                'CO2_production': co2_production,
                'Status': status
            })
            
            print(f"{atpm_value:<10.1f} {biomass_flux:<15.6f} {co2_production:<18.6f} {status:<10}")
        else:
            results.append({
                'ATPM': atpm_value,
                'Biomass_flux': 0,
                'CO2_production': 0,
                'Status': solution.status
            })
            print(f"{atpm_value:<10.1f} {'N/A':<15} {'N/A':<18} {solution.status:<10}")
    
    df_results = pd.DataFrame(results)
    
    return df_results, solution

def analyze_pathway_fluxes(model, biomass_rxn, atpm_rxn, atpm_values=[0, 40, 80, 100]):
    """
    특정 ATPM 값에서 경로 플럭스 분석
    """
    print("\n" + "="*70)
    print("경로 플럭스 분석 (다양한 ATPM 값)")
    print("="*70)
    
    model.objective = biomass_rxn.id
    
    # 분석할 반응들
    pathway_reactions = {
        'CS': 'Citrate Synthase',
        'ICDHx': 'Isocitrate Dehydrogenase (NAD+)',
        'AKGDH': 'α-Ketoglutarate Dehydrogenase',
        'SUCD': 'Succinate Dehydrogenase',
        'FUM': 'Fumarase',
        'ICL': 'Isocitrate Lyase',
        'MALS': 'Malate Synthase',
        'SUCOAACTr': 'Acetate-CoA Transferase',
        'ACS': 'Acetyl-CoA Synthetase'
    }
    
    flux_data = []
    
    for atpm_value in atpm_values:
        atpm_rxn.lower_bound = atpm_value
        atpm_rxn.upper_bound = max(atpm_value, 1000)
        
        solution = model.optimize()
        
        if solution.status == 'optimal':
            for rxn_id, rxn_name in pathway_reactions.items():
                try:
                    flux = solution.fluxes.get(rxn_id, 0)
                    flux_data.append({
                        'ATPM': atpm_value,
                        'Reaction_ID': rxn_id,
                        'Reaction_Name': rxn_name,
                        'Flux': flux
                    })
                except KeyError:
                    pass
    
    df_fluxes = pd.DataFrame(flux_data)
    
    # 피벗 테이블 생성
    if len(df_fluxes) > 0:
        pivot = df_fluxes.pivot(index='Reaction_ID', columns='ATPM', values='Flux')
        print("\n경로 플럭스 변화 (ATPM별):")
        print(pivot.to_string())
    
    return df_fluxes

def main():
    print("="*70)
    print("Acetate 기반 FBA - ATP Maintenance Sweep 분석")
    print("논문 방법론 적용")
    print("="*70)
    
    # 모델 로드
    model = load_model("BaseModel.xml")
    
    # Biomass 반응 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass: {biomass_rxn.id}")
    
    # ATP Maintenance 반응 찾기/생성
    atpm_rxn = find_atp_maintenance_reaction(model)
    if not atpm_rxn:
        print("[ERROR] ATP Maintenance 반응을 찾거나 생성할 수 없습니다")
        return
    
    print(f"[OK] ATP Maintenance: {atpm_rxn.id}")
    
    # Medium 설정 (논문 값 사용)
    model = setup_acetate_medium(model, acetate_uptake=-19.0, oxygen_uptake=-100.0)
    
    # ATPM Sweep 수행
    df_results, solution = run_atpm_sweep(model, biomass_rxn, atpm_rxn)
    
    # 결과 저장
    df_results.to_csv('atpm_sweep_results.csv', index=False)
    print(f"\n[OK] 결과 저장: atpm_sweep_results.csv")
    
    # 경로 플럭스 분석
    df_fluxes = analyze_pathway_fluxes(model, biomass_rxn, atpm_rxn, atpm_values=[0, 40, 80, 100])
    
    if len(df_fluxes) > 0:
        df_fluxes.to_csv('pathway_fluxes_atpm.csv', index=False)
        print(f"[OK] 경로 플럭스 저장: pathway_fluxes_atpm.csv")
    
    # 최종 요약
    print("\n" + "="*70)
    print("최종 요약")
    print("="*70)
    
    if len(df_results) > 0:
        max_growth = df_results['Biomass_flux'].max()
        max_growth_atpm = df_results.loc[df_results['Biomass_flux'].idxmax(), 'ATPM']
        
        print(f"\n최대 성장률: {max_growth:.6f} 1/h (ATPM = {max_growth_atpm:.1f})")
        
        # 성장 가능한 ATPM 범위
        growth_positive = df_results[df_results['Biomass_flux'] > 1e-6]
        if len(growth_positive) > 0:
            print(f"성장 가능한 ATPM 범위: {growth_positive['ATPM'].min():.1f} ~ {growth_positive['ATPM'].max():.1f} mmol ATP/gDCW/h")
        
        # CO2 production 증가 확인
        if 'CO2_production' in df_results.columns:
            co2_increase = df_results['CO2_production'].max() - df_results['CO2_production'].min()
            print(f"CO2 생성량 증가: {co2_increase:.6f} mmol/gDCW/h")
    
    print("\n[OK] ATPM Sweep 분석 완료!")
    print("="*70)

if __name__ == "__main__":
    main()
