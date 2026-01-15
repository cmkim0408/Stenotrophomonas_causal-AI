#!/usr/bin/env python
"""
성장했던 FBA와 현재 FBA 비교
- 이전에 성장했던 설정 재현
- 현재 설정과 비교
- 레퍼런스 모델과 비교
"""

import cobra
from pathlib import Path
import pandas as pd
import csv

def load_model(model_path):
    model = cobra.io.read_sbml_model(model_path)
    return model

def setup_acetate_medium_growing(model):
    """성장했을 때의 미디어 설정 (test_fba_with_missing_reactions.py 기준)"""
    exchanges = {
        'EX_ac_e': (-10, 1000),      # Acetate uptake
        'EX_o2_e': (-1000, 1000),    # Oxygen
        'EX_h2o_e': (-1000, 1000),   # Water
        'EX_h_e': (-1000, 1000),     # H+
        'EX_nh4_e': (-10, 1000),     # Ammonium
        'EX_pi_e': (-10, 1000),      # Phosphate
        'EX_so4_e': (-10, 1000),     # Sulfate
        'EX_k_e': (-1000, 1000),     # Potassium
        'EX_na1_e': (-1000, 1000),   # Sodium
        'EX_mg2_e': (-1000, 1000),   # Magnesium
        'EX_ca2_e': (-1000, 1000),   # Calcium
        'EX_fe2_e': (-10, 1000),     # Fe2+
        'EX_fe3_e': (-10, 1000),     # Fe3+
        'EX_hco3_e': (-1000, 1000),  # Bicarbonate
        'EX_co2_e': (-1000, 1000),   # CO2
    }
    
    for ex_id, bounds in exchanges.items():
        if ex_id in model.reactions:
            model.reactions.get_by_id(ex_id).bounds = bounds
    
    return model

def apply_media_from_tsv(model, tsv_path):
    """TSV 파일에서 미디어 적용 (현재 사용 중)"""
    with open(tsv_path, "r", encoding="utf-8-sig", newline="") as f:
        raw_lines = f.readlines()
    
    lines = []
    for ln in raw_lines:
        ln_stripped = ln.strip()
        if ln_stripped and not ln_stripped.startswith("#"):
            lines.append(ln)
    
    rdr = csv.DictReader(lines, delimiter="\t")
    fieldnames = list(rdr.fieldnames)
    
    if len(fieldnames) >= 3 and fieldnames[1].replace(".", "").replace("-", "").isdigit():
        col_ex = fieldnames[0]
        col_lb = fieldnames[1]
        col_ub = fieldnames[2]
    else:
        col_ex = fieldnames[0]
        col_lb = fieldnames[1] if len(fieldnames) > 1 else fieldnames[0]
        col_ub = fieldnames[2] if len(fieldnames) > 2 else fieldnames[1]
    
    applied = 0
    for row in rdr:
        ex_id = row[col_ex].strip()
        if not ex_id or ex_id.lower() == "exchange_id":
            continue
        
        try:
            lb = float(row[col_lb])
            ub = float(row[col_ub])
            
            if ex_id in model.reactions:
                rxn = model.reactions.get_by_id(ex_id)
                rxn.lower_bound = lb
                rxn.upper_bound = ub
                applied += 1
        except (KeyError, ValueError, TypeError):
            continue
    
    return applied

def add_missing_reactions(model, ref_model, missing_df):
    """누락된 반응 추가"""
    added = []
    
    for idx, row in missing_df.iterrows():
        # 컬럼명 확인
        if 'rxn_id' in row:
            rxn_id = row['rxn_id']
        elif 'reaction_id' in row:
            rxn_id = row['reaction_id']
        else:
            continue
        
        if rxn_id in [r.id for r in model.reactions]:
            continue
        
        try:
            ref_rxn = ref_model.reactions.get_by_id(rxn_id)
            
            new_rxn = cobra.Reaction(rxn_id)
            new_rxn.name = ref_rxn.name
            new_rxn.lower_bound = ref_rxn.lower_bound
            new_rxn.upper_bound = ref_rxn.upper_bound
            
            metabolites_dict = {}
            for met, coeff in ref_rxn.metabolites.items():
                if met.id in [m.id for m in model.metabolites]:
                    metabolites_dict[model.metabolites.get_by_id(met.id)] = coeff
                else:
                    new_met = cobra.Metabolite(
                        met.id,
                        name=met.name,
                        compartment=met.compartment,
                        formula=met.formula if hasattr(met, 'formula') else None,
                        charge=met.charge if hasattr(met, 'charge') else None
                    )
                    model.add_metabolites([new_met])
                    metabolites_dict[new_met] = coeff
            
            new_rxn.add_metabolites(metabolites_dict)
            model.add_reactions([new_rxn])
            added.append(rxn_id)
        except:
            pass
    
    return added

def run_fba_and_analyze(model, label):
    """FBA 실행 및 분석"""
    print(f"\n{'='*80}")
    print(f"{label} - FBA 실행")
    print(f"{'='*80}")
    
    model.objective = 'Growth'
    solution = model.optimize()
    
    print(f"\n[FBA 결과]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    
    # 주요 반응 플럭스
    key_reactions = [
        'ACS_ADP', 'ACS', 'CS', 'ICL', 'MALS', 'SUCDi', 'PEPCK_ATP',
        'EX_ac_e', 'EX_o2_e', 'EX_hco3_e',
        'ATPS4rpp', 'NADH16pp', 'CYTBO3_4pp',
        'ACtexi', 'ADK1'
    ]
    
    print(f"\n[주요 반응 플럭스]")
    active_fluxes = {}
    for rxn_id in key_reactions:
        if rxn_id in model.reactions:
            flux = solution.fluxes.get(rxn_id, 0.0)
            if abs(flux) > 1e-6:
                print(f"  {rxn_id:20s}: {flux:>12.6f}")
                active_fluxes[rxn_id] = flux
    
    return solution, active_fluxes

def compare_scenarios():
    """성장했던 시나리오와 현재 시나리오 비교"""
    base_path = Path(__file__).parent.parent
    new_model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    missing_csv = base_path / "Stenotrophomonas-causal AI" / "missing_reactions_vs_reference.csv"
    media_tsv = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5" / "Acetate_YE0p5__nocmnt__normalized.tsv"
    
    # 모델 로드
    new_model = load_model(str(new_model_path))
    ref_model = load_model(str(ref_model_path))
    
    # 누락된 반응 로드
    missing_df = pd.read_csv(missing_csv)
    
    print("="*80)
    print("성장했던 FBA vs 현재 FBA 비교")
    print("="*80)
    
    # 시나리오 1: 성장했던 설정 (99개 또는 114개 반응 추가 + 성장 미디어)
    print(f"\n[시나리오 1] 성장했던 설정 재현")
    model_growing = new_model.copy()
    
    # 99개 반응 추가 (9 active + 90 medium)
    # 또는 114개 모두 추가
    added = add_missing_reactions(model_growing, ref_model, missing_df)
    print(f"  추가된 반응: {len(added)}개")
    
    # 성장했던 미디어 설정
    model_growing = setup_acetate_medium_growing(model_growing)
    
    # ATPM 설정 확인
    atpm_rxn = model_growing.reactions.get_by_id('ATPM')
    print(f"  ATPM bounds: [{atpm_rxn.lower_bound}, {atpm_rxn.upper_bound}]")
    
    sol_growing, fluxes_growing = run_fba_and_analyze(model_growing, "시나리오 1: 성장했던 설정")
    
    # 시나리오 2: 현재 설정 (TSV 미디어)
    print(f"\n[시나리오 2] 현재 설정 (TSV 미디어)")
    model_current = new_model.copy()
    
    # 같은 반응 추가
    added_current = add_missing_reactions(model_current, ref_model, missing_df)
    print(f"  추가된 반응: {len(added_current)}개")
    
    # 현재 미디어 설정 (TSV)
    for rxn in model_current.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    if media_tsv.exists():
        applied = apply_media_from_tsv(model_current, media_tsv)
        print(f"  TSV 미디어 적용: {applied}개 반응")
    
    # ATPM=0 설정
    atpm_rxn = model_current.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    print(f"  ATPM bounds: [{atpm_rxn.lower_bound}, {atpm_rxn.upper_bound}]")
    
    sol_current, fluxes_current = run_fba_and_analyze(model_current, "시나리오 2: 현재 설정")
    
    # 비교
    print(f"\n{'='*80}")
    print("비교 결과")
    print(f"{'='*80}")
    
    print(f"\n[성장률 비교]")
    print(f"  성장했던 설정: {sol_growing.objective_value:.6f}")
    print(f"  현재 설정: {sol_current.objective_value:.6f}")
    
    print(f"\n[주요 반응 플럭스 비교]")
    all_reactions = set(list(fluxes_growing.keys()) + list(fluxes_current.keys()))
    
    for rxn_id in sorted(all_reactions):
        flux_growing = fluxes_growing.get(rxn_id, 0.0)
        flux_current = fluxes_current.get(rxn_id, 0.0)
        
        if abs(flux_growing) > 1e-6 or abs(flux_current) > 1e-6:
            status = ""
            if abs(flux_growing) > 1e-6 and abs(flux_current) < 1e-6:
                status = " [성장했을 때만 작동]"
            elif abs(flux_growing) < 1e-6 and abs(flux_current) > 1e-6:
                status = " [현재만 작동]"
            elif abs(flux_growing) > 1e-6 and abs(flux_current) > 1e-6:
                status = " [둘 다 작동]"
            
            print(f"  {rxn_id:20s}: 성장={flux_growing:>12.6f}, 현재={flux_current:>12.6f}{status}")
    
    # 차이점 분석
    print(f"\n[차이점 분석]")
    
    if sol_growing.objective_value > 1e-6 and sol_current.objective_value < 1e-6:
        print(f"  [핵심 차이] 성장했던 설정은 성장하지만, 현재 설정은 성장하지 않음")
        print(f"  -> 미디어 설정 차이 또는 ATPM 설정 차이")
        
        # Exchange bounds 비교
        print(f"\n[Exchange bounds 비교]")
        key_exchanges = ['EX_ac_e', 'EX_o2_e', 'EX_hco3_e', 'EX_nh4_e', 'EX_pi_e']
        for ex_id in key_exchanges:
            if ex_id in model_growing.reactions and ex_id in model_current.reactions:
                rxn_g = model_growing.reactions.get_by_id(ex_id)
                rxn_c = model_current.reactions.get_by_id(ex_id)
                print(f"  {ex_id:20s}: 성장=[{rxn_g.lower_bound:>8.1f}, {rxn_g.upper_bound:>8.1f}], 현재=[{rxn_c.lower_bound:>8.1f}, {rxn_c.upper_bound:>8.1f}]")

def main():
    compare_scenarios()

if __name__ == "__main__":
    main()
