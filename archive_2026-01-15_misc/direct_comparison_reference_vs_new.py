#!/usr/bin/env python
"""
레퍼런스 모델 vs 신규 모델 직접 비교
- 레퍼런스 모델의 실제 FBA 설정 확인
- 신규 모델과 정확히 동일한 조건으로 비교
- 차이점 명확히 파악
"""

import cobra
from pathlib import Path
import pandas as pd
import csv

def load_model(model_path):
    model = cobra.io.read_sbml_model(model_path)
    return model

def apply_media_from_tsv(model, tsv_path):
    """TSV 파일에서 미디어 적용"""
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

def run_fba_simple(model, atpm_value=0):
    """간단한 FBA 실행"""
    atpm_rxn = model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = atpm_value
    atpm_rxn.upper_bound = 1000
    
    model.objective = 'Growth'
    solution = model.optimize()
    
    return solution

def compare_reference_vs_new():
    """레퍼런스 모델 vs 신규 모델 직접 비교"""
    base_path = Path(__file__).parent.parent
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    new_model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    media_tsv = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5" / "Acetate_YE0p5__nocmnt__normalized.tsv"
    
    # 모델 로드
    print("="*80)
    print("레퍼런스 모델 vs 신규 모델 직접 비교")
    print("="*80)
    
    ref_model = load_model(str(ref_model_path))
    new_model = load_model(str(new_model_path))
    
    print(f"\n[모델 정보]")
    print(f"  레퍼런스 모델: {len(ref_model.reactions)}개 반응, {len(ref_model.metabolites)}개 메타볼라이트")
    print(f"  신규 모델: {len(new_model.reactions)}개 반응, {len(new_model.metabolites)}개 메타볼라이트")
    
    # 미디어 설정
    print(f"\n[미디어 설정]")
    for rxn in ref_model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    for rxn in new_model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    if media_tsv.exists():
        ref_applied = apply_media_from_tsv(ref_model, media_tsv)
        new_applied = apply_media_from_tsv(new_model, media_tsv)
        print(f"  레퍼런스 모델: {ref_applied}개 exchange 설정")
        print(f"  신규 모델: {new_applied}개 exchange 설정")
    
    # ATPM 설정 확인
    print(f"\n[ATPM 설정 확인]")
    ref_atpm = ref_model.reactions.get_by_id('ATPM')
    new_atpm = new_model.reactions.get_by_id('ATPM')
    print(f"  레퍼런스 ATPM bounds: [{ref_atpm.lower_bound}, {ref_atpm.upper_bound}]")
    print(f"  신규 ATPM bounds: [{new_atpm.lower_bound}, {new_atpm.upper_bound}]")
    print(f"  레퍼런스 ATPM 반응식: {ref_atpm.reaction}")
    print(f"  신규 ATPM 반응식: {new_atpm.reaction}")
    
    # ATPM=0으로 설정하고 테스트
    print(f"\n[ATPM=0으로 설정하고 테스트]")
    
    ref_sol = run_fba_simple(ref_model, atpm_value=0)
    new_sol = run_fba_simple(new_model, atpm_value=0)
    
    print(f"\n[레퍼런스 모델 FBA 결과]")
    print(f"  상태: {ref_sol.status}")
    print(f"  성장률: {ref_sol.objective_value:.6f}")
    
    key_reactions = ['ACS_ADP', 'ACS', 'CS', 'ICL', 'MALS', 'SUCDi', 'ATPS4rpp', 'NADH16pp']
    print(f"\n  주요 반응 플럭스:")
    for rxn_id in key_reactions:
        if rxn_id in ref_model.reactions:
            flux = ref_sol.fluxes.get(rxn_id, 0.0)
            if abs(flux) > 1e-6:
                print(f"    {rxn_id:20s}: {flux:>12.6f}")
    
    print(f"\n[신규 모델 FBA 결과]")
    print(f"  상태: {new_sol.status}")
    print(f"  성장률: {new_sol.objective_value:.6f}")
    
    print(f"\n  주요 반응 플럭스:")
    for rxn_id in key_reactions:
        if rxn_id in new_model.reactions:
            flux = new_sol.fluxes.get(rxn_id, 0.0)
            if abs(flux) > 1e-6:
                print(f"    {rxn_id:20s}: {flux:>12.6f}")
    
    # 차이점 분석
    print(f"\n[차이점 분석]")
    
    # ACS_ADP 확인
    ref_acs_adp = ref_sol.fluxes.get('ACS_ADP', 0.0) if 'ACS_ADP' in ref_model.reactions else None
    new_acs_adp = new_sol.fluxes.get('ACS_ADP', 0.0) if 'ACS_ADP' in new_model.reactions else None
    
    print(f"\n  ACS_ADP:")
    print(f"    레퍼런스: {ref_acs_adp}")
    print(f"    신규: {new_acs_adp}")
    if ref_acs_adp is None:
        print(f"    -> 레퍼런스 모델에 ACS_ADP 없음")
    elif new_acs_adp is None:
        print(f"    -> 신규 모델에 ACS_ADP 없음 (추가 필요)")
    elif abs(ref_acs_adp) > 1e-6 and abs(new_acs_adp) < 1e-6:
        print(f"    -> 레퍼런스는 작동하지만 신규는 작동 안 함")
    
    # CS 확인
    ref_cs = ref_sol.fluxes.get('CS', 0.0)
    new_cs = new_sol.fluxes.get('CS', 0.0)
    
    print(f"\n  CS:")
    print(f"    레퍼런스: {ref_cs:.6f}")
    print(f"    신규: {new_cs:.6f}")
    if abs(ref_cs) > 1e-6 and abs(new_cs) < 1e-6:
        print(f"    -> 레퍼런스는 작동하지만 신규는 작동 안 함")
    
    # ICL 확인
    ref_icl = ref_sol.fluxes.get('ICL', 0.0)
    new_icl = new_sol.fluxes.get('ICL', 0.0)
    
    print(f"\n  ICL:")
    print(f"    레퍼런스: {ref_icl:.6f}")
    print(f"    신규: {new_icl:.6f}")
    if abs(ref_icl) > 1e-6 and abs(new_icl) < 1e-6:
        print(f"    -> 레퍼런스는 작동하지만 신규는 작동 안 함")
    
    # Exchange 플럭스 비교
    print(f"\n[Exchange 플럭스 비교]")
    key_exchanges = ['EX_ac_e', 'EX_o2_e', 'EX_hco3_e', 'EX_nh4_e']
    for ex_id in key_exchanges:
        if ex_id in ref_model.reactions and ex_id in new_model.reactions:
            ref_flux = ref_sol.fluxes.get(ex_id, 0.0)
            new_flux = new_sol.fluxes.get(ex_id, 0.0)
            ref_bounds = f"[{ref_model.reactions.get_by_id(ex_id).lower_bound}, {ref_model.reactions.get_by_id(ex_id).upper_bound}]"
            new_bounds = f"[{new_model.reactions.get_by_id(ex_id).lower_bound}, {new_model.reactions.get_by_id(ex_id).upper_bound}]"
            
            if abs(ref_flux - new_flux) > 1e-6:
                print(f"  {ex_id:20s}: 레퍼런스={ref_flux:>12.6f} {ref_bounds}, 신규={new_flux:>12.6f} {new_bounds}")

def main():
    compare_reference_vs_new()
    
    print("\n" + "="*80)
    print("결론")
    print("="*80)
    print(f"\n[핵심 질문]")
    print(f"  레퍼런스 모델에서는 FBA가 잘 돌아가는데, 신규 모델에서는 왜 안 될까?")
    print(f"  -> 위의 비교 결과를 참고하여 차이점 확인")

if __name__ == "__main__":
    main()
