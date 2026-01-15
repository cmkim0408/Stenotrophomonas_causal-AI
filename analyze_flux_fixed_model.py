#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ACt2rpp 제거 후 모델의 Flux 분석 및 Export

분석 내용:
1. 주요 대사 경로의 flux 확인
2. 탄소 전환 경로 분석 (Acetate → Acetyl-CoA → TCA/Glyoxylate)
3. 주요 반응의 flux를 CSV로 export
"""

import cobra
from pathlib import Path
import pandas as pd
import sys

def load_model(model_path):
    try:
        model = cobra.io.read_sbml_model(str(model_path))
        print(f"[OK] 모델 로드 완료 (반응 수: {len(model.reactions)})")
        return model
    except Exception as e:
        print(f"[ERROR] 모델 로드 실패: {e}")
        sys.exit(1)

def setup_media_forced(model):
    """배지 조건을 강제로 고정"""
    if 'EX_ac_e' in model.reactions:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.lower_bound = -19.0
        ex_ac.upper_bound = -19.0
    
    if 'EX_o2_e' in model.reactions:
        ex_o2 = model.reactions.get_by_id('EX_o2_e')
        ex_o2.lower_bound = -100.0
        ex_o2.upper_bound = 1000.0
    
    essential_exchanges = {
        'EX_nh4_e': (-1000.0, 1000.0),
        'EX_pi_e': (-1000.0, 1000.0),
        'EX_so4_e': (-1000.0, 1000.0),
        'EX_mg2_e': (-1000.0, 1000.0),
        'EX_k_e': (-1000.0, 1000.0),
        'EX_na1_e': (-1000.0, 1000.0),
        'EX_fe2_e': (-1000.0, 1000.0),
        'EX_fe3_e': (-1000.0, 1000.0),
        'EX_h2o_e': (-1000.0, 1000.0),
        'EX_h_e': (-1000.0, 1000.0),
        'EX_co2_e': (-1000.0, 1000.0),
        'EX_hco3_e': (-1000.0, 1000.0),
        'EX_nac_e': (-1000.0, 1000.0),
        'EX_ncam_e': (-1000.0, 1000.0),
    }
    
    for ex_id, (lb, ub) in essential_exchanges.items():
        if ex_id in model.reactions:
            ex_rxn = model.reactions.get_by_id(ex_id)
            ex_rxn.lower_bound = lb
            ex_rxn.upper_bound = ub
    
    return model

def analyze_carbon_pathway(model, solution):
    """탄소 전환 경로 분석"""
    print("\n" + "="*70)
    print("탄소 전환 경로 분석")
    print("="*70)
    
    # 1. Acetate uptake
    print("\n[1. Acetate Uptake]")
    if 'EX_ac_e' in model.reactions:
        ex_ac_flux = solution.fluxes.get('EX_ac_e', 0.0)
        print(f"  EX_ac_e: {ex_ac_flux:.6f}")
    
    # 2. Acetate → Acetyl-CoA
    print("\n[2. Acetate → Acetyl-CoA]")
    ac_to_accoa_rxns = {
        'ACS': 'Acetyl-CoA synthetase (AMP-forming)',
        'ACS_ADP': 'Acetyl-CoA synthetase (ADP-forming)',
        'AADb': 'Acetyl-adenylate synthetase',
        'AADa': 'Acetyl-CoA synthetase (adenylate intermediate)',
    }
    
    for rxn_id, desc in ac_to_accoa_rxns.items():
        if rxn_id in model.reactions:
            flux = solution.fluxes.get(rxn_id, 0.0)
            if abs(flux) > 1e-6:
                rxn = model.reactions.get_by_id(rxn_id)
                print(f"  {rxn_id:15s}: {flux:10.6f} ({desc})")
                print(f"                 반응식: {rxn.reaction}")
    
    # 3. Acetyl-CoA → TCA vs Glyoxylate
    print("\n[3. Acetyl-CoA → TCA Cycle vs Glyoxylate Shunt]")
    
    # TCA cycle (CS)
    if 'CS' in model.reactions:
        cs_flux = solution.fluxes.get('CS', 0.0)
        if abs(cs_flux) > 1e-6:
            cs_rxn = model.reactions.get_by_id('CS')
            print(f"  CS (Citrate Synthase): {cs_flux:.6f}")
            print(f"    반응식: {cs_rxn.reaction}")
            print(f"    경로: Acetyl-CoA + OAA → Citrate (TCA cycle)")
    
    # Glyoxylate shunt
    if 'ICL' in model.reactions:
        icl_flux = solution.fluxes.get('ICL', 0.0)
        if abs(icl_flux) > 1e-6:
            icl_rxn = model.reactions.get_by_id('ICL')
            print(f"  ICL (Isocitrate Lyase): {icl_flux:.6f}")
            print(f"    반응식: {icl_rxn.reaction}")
            print(f"    경로: Isocitrate → Glyoxylate + Succinate (Glyoxylate shunt)")
    
    if 'MALS' in model.reactions:
        mals_flux = solution.fluxes.get('MALS', 0.0)
        if abs(mals_flux) > 1e-6:
            mals_rxn = model.reactions.get_by_id('MALS')
            print(f"  MALS (Malate Synthase): {mals_flux:.6f}")
            print(f"    반응식: {mals_rxn.reaction}")
            print(f"    경로: Glyoxylate + Acetyl-CoA → Malate (Glyoxylate shunt)")
    
    # 4. TCA cycle 주요 반응
    print("\n[4. TCA Cycle 주요 반응]")
    tca_rxns = {
        'ACONT': 'Aconitase (Citrate → Isocitrate)',
        'ICDHx': 'Isocitrate Dehydrogenase (Isocitrate → α-KG)',
        'AKGDH': 'α-Ketoglutarate Dehydrogenase (α-KG → Succinyl-CoA)',
        'SUCDi': 'Succinate Dehydrogenase (Succinate → Fumarate)',
        'FUM': 'Fumarase (Fumarate → Malate)',
        'MDH': 'Malate Dehydrogenase (Malate → OAA)',
    }
    
    for rxn_id, desc in tca_rxns.items():
        if rxn_id in model.reactions:
            flux = solution.fluxes.get(rxn_id, 0.0)
            if abs(flux) > 1e-6:
                print(f"  {rxn_id:10s}: {flux:10.6f} ({desc})")
    
    # 5. Gluconeogenesis
    print("\n[5. Gluconeogenesis 경로]")
    glc_rxns = {
        'PEPCK_ATP': 'PEP Carboxykinase (ATP) (OAA → PEP)',
        'PCK': 'PEP Carboxykinase',
        'PYCK': 'PEP Carboxykinase',
        'PYK': 'Pyruvate Kinase (PEP → Pyruvate)',
        'PPC': 'PEP Carboxylase (PEP → OAA)',
    }
    
    for rxn_id, desc in glc_rxns.items():
        if rxn_id in model.reactions:
            flux = solution.fluxes.get(rxn_id, 0.0)
            if abs(flux) > 1e-6:
                rxn = model.reactions.get_by_id(rxn_id)
                print(f"  {rxn_id:15s}: {flux:10.6f} ({desc})")
                print(f"                 반응식: {rxn.reaction}")

def analyze_energy_pathway(model, solution):
    """에너지 생성 경로 분석"""
    print("\n" + "="*70)
    print("에너지 생성 경로 분석")
    print("="*70)
    
    # ATP 생성
    print("\n[ATP 생성 경로]")
    atp_gen_rxns = {
        'ATPS4rpp': 'ATP Synthase (4H+)',
        'ATPS': 'ATP Synthase',
        'PYK': 'Pyruvate Kinase',
        'PGK': 'Phosphoglycerate Kinase',
    }
    
    for rxn_id, desc in atp_gen_rxns.items():
        if rxn_id in model.reactions:
            flux = solution.fluxes.get(rxn_id, 0.0)
            if abs(flux) > 1e-6:
                rxn = model.reactions.get_by_id(rxn_id)
                # ATP 생성 여부 확인
                if 'atp_c' in [m.id for m in rxn.metabolites]:
                    atp_coeff = rxn.metabolites.get(model.metabolites.get_by_id('atp_c'), 0)
                    if atp_coeff > 0:  # ATP 생성
                        atp_gen = flux * atp_coeff
                        print(f"  {rxn_id:15s}: 플럭스 {flux:10.6f}, ATP 생성 {atp_gen:.6f} ({desc})")
    
    # ATP 소비
    print("\n[ATP 소비 경로]")
    atp_cons_rxns = {
        'ATPM': 'ATP Maintenance',
        'Growth': 'Biomass Growth',
        'ACS': 'Acetyl-CoA Synthetase',
        'PEPCK_ATP': 'PEP Carboxykinase',
    }
    
    for rxn_id, desc in atp_cons_rxns.items():
        if rxn_id in model.reactions:
            flux = solution.fluxes.get(rxn_id, 0.0)
            if abs(flux) > 1e-6:
                rxn = model.reactions.get_by_id(rxn_id)
                if 'atp_c' in [m.id for m in rxn.metabolites]:
                    atp_coeff = rxn.metabolites.get(model.metabolites.get_by_id('atp_c'), 0)
                    if atp_coeff < 0:  # ATP 소비
                        atp_cons = abs(flux * atp_coeff)
                        print(f"  {rxn_id:15s}: 플럭스 {flux:10.6f}, ATP 소비 {atp_cons:.6f} ({desc})")

def export_fluxes(model, solution, output_file):
    """주요 반응의 flux를 CSV로 export"""
    print("\n" + "="*70)
    print("Flux Export")
    print("="*70)
    
    # 주요 반응 카테고리
    categories = {
        'Acetate_Uptake': ['EX_ac_e', 'ACtexi', 'ACt'],
        'Acetyl_CoA_Synthesis': ['ACS', 'ACS_ADP', 'AADb', 'AADa'],
        'TCA_Cycle': ['CS', 'ACONT', 'ICDHx', 'AKGDH', 'SUCDi', 'FUM', 'MDH'],
        'Glyoxylate_Shunt': ['ICL', 'MALS'],
        'Gluconeogenesis': ['MAEB', 'PYK', 'PPC'],
        'Energy_Generation': ['ATPS4rpp', 'ATPS', 'ATPM'],
        'Biomass': ['Growth'],
        'Exchange': ['EX_o2_e', 'EX_co2_e', 'EX_nh4_e', 'EX_pi_e'],
    }
    
    flux_data = []
    
    for category, rxn_ids in categories.items():
        for rxn_id in rxn_ids:
            if rxn_id in model.reactions:
                flux = solution.fluxes.get(rxn_id, 0.0)
                rxn = model.reactions.get_by_id(rxn_id)
                flux_data.append({
                    'category': category,
                    'reaction_id': rxn_id,
                    'reaction_name': rxn.name if rxn.name else '',
                    'reaction_equation': rxn.reaction,
                    'flux': flux,
                    'lower_bound': rxn.lower_bound,
                    'upper_bound': rxn.upper_bound,
                })
    
    # 모든 non-zero flux 반응도 포함
    all_nonzero = []
    for rxn in model.reactions:
        flux = solution.fluxes.get(rxn.id, 0.0)
        if abs(flux) > 1e-6:
            all_nonzero.append({
                'reaction_id': rxn.id,
                'reaction_name': rxn.name if rxn.name else '',
                'reaction_equation': rxn.reaction,
                'flux': flux,
                'lower_bound': rxn.lower_bound,
                'upper_bound': rxn.upper_bound,
            })
    
    # DataFrame 생성
    df_main = pd.DataFrame(flux_data)
    df_all = pd.DataFrame(all_nonzero)
    
    # CSV 저장
    output_path = Path(output_file)
    df_main.to_csv(output_path, index=False, encoding='utf-8-sig')
    print(f"\n[주요 반응 Flux 저장] {output_path}")
    print(f"  총 {len(df_main)}개 반응")
    
    # 모든 non-zero flux 저장
    output_all_path = output_path.parent / (output_path.stem + '_all_nonzero.csv')
    df_all.to_csv(output_all_path, index=False, encoding='utf-8-sig')
    print(f"[모든 Non-zero Flux 저장] {output_all_path}")
    print(f"  총 {len(df_all)}개 반응")
    
    return df_main, df_all

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_final_cleaned.xml"
    output_dir = base_path / "Stenotrophomonas-causal AI"
    
    print("="*70)
    print("최종 모델 Flux 분석 (모든 아티팩트 제거)")
    print("="*70)
    
    # 모델 로드
    model = load_model(model_path)
    
    # 배지 조건 설정
    model = setup_media_forced(model)
    
    # ATPM=0 설정
    if 'ATPM' in model.reactions:
        atpm = model.reactions.get_by_id('ATPM')
        atpm.lower_bound = 0.0
        atpm.upper_bound = 0.0
        print(f"\n[ATPM 설정] bounds: [{atpm.lower_bound}, {atpm.upper_bound}]")
    
    # FBA 실행
    solution = model.optimize()
    print(f"\n[FBA 결과]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    
    # 탄소 전환 경로 분석
    analyze_carbon_pathway(model, solution)
    
    # 에너지 생성 경로 분석
    analyze_energy_pathway(model, solution)
    
    # Flux export
    output_file = output_dir / "flux_analysis_final_complete.csv"
    df_main, df_all = export_fluxes(model, solution, output_file)
    
    print("\n" + "="*70)
    print("완료")
    print("="*70)
    
    return model, solution, df_main, df_all

if __name__ == "__main__":
    model, solution, df_main, df_all = main()
