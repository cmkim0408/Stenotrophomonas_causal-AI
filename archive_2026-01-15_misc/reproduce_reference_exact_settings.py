#!/usr/bin/env python
"""
레퍼런스 FBA 결과의 정확한 설정 재현
- 레퍼런스 FBA 결과 파일의 실제 플럭스 값 확인
- 그 설정을 재현하여 신규 모델과 비교
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

def apply_reference_exchange_fluxes(model, ref_fluxes):
    """레퍼런스 FBA 결과의 Exchange 플럭스를 bounds로 설정"""
    # Exchange 플럭스를 bounds로 변환
    # 음수 플럭스 = uptake -> lower_bound 설정
    # 양수 플럭스 = secretion -> upper_bound 설정
    
    for ex_id, flux in ref_fluxes.items():
        if ex_id in model.reactions:
            rxn = model.reactions.get_by_id(ex_id)
            if flux < 0:  # Uptake
                rxn.lower_bound = flux
                rxn.upper_bound = 1000
            elif flux > 0:  # Secretion
                rxn.lower_bound = -1000
                rxn.upper_bound = flux
            else:  # flux = 0
                rxn.lower_bound = 0
                rxn.upper_bound = 0

def reproduce_reference_settings():
    """레퍼런스 설정 재현"""
    base_path = Path(__file__).parent.parent
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    new_model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    media_tsv = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5" / "Acetate_YE0p5__nocmnt__normalized.tsv"
    fba_file = base_path / "Stenotrophomonas" / "fba_flux_gradient_acid.csv"
    
    # 레퍼런스 FBA 결과 로드
    df_fba = pd.read_csv(fba_file)
    
    print("="*80)
    print("레퍼런스 FBA 결과의 정확한 설정 재현")
    print("="*80)
    
    # 레퍼런스 FBA 결과에서 ATPM=0일 때 Exchange 플럭스 추출
    ref_exchanges = {}
    ex_list = ['EX_ac_e', 'EX_o2_e', 'EX_hco3_e', 'EX_nh4_e', 'EX_pi_e']
    
    for ex_id in ex_list:
        row = df_fba[df_fba.iloc[:,0] == ex_id]
        if not row.empty:
            flux = row['ATPM_0'].iloc[0]
            ref_exchanges[ex_id] = flux
    
    print(f"\n[레퍼런스 FBA 결과 (ATPM=0) Exchange 플럭스]")
    for ex_id, flux in ref_exchanges.items():
        print(f"  {ex_id:20s}: {flux:>12.6f}")
    
    # 모델 로드
    ref_model = load_model(str(ref_model_path))
    new_model = load_model(str(new_model_path))
    
    # 미디어 설정
    for rxn in ref_model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    for rxn in new_model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    if media_tsv.exists():
        apply_media_from_tsv(ref_model, media_tsv)
        apply_media_from_tsv(new_model, media_tsv)
    
    # 레퍼런스 Exchange 플럭스를 bounds로 설정
    apply_reference_exchange_fluxes(ref_model, ref_exchanges)
    apply_reference_exchange_fluxes(new_model, ref_exchanges)
    
    # ATPM=0 설정
    ref_atpm = ref_model.reactions.get_by_id('ATPM')
    ref_atpm.lower_bound = 0
    ref_atpm.upper_bound = 1000
    
    new_atpm = new_model.reactions.get_by_id('ATPM')
    new_atpm.lower_bound = 0
    new_atpm.upper_bound = 1000
    
    # FBA 실행
    print(f"\n[레퍼런스 모델 FBA (레퍼런스 Exchange 플럭스 사용)]")
    ref_model.objective = 'Growth'
    ref_sol = ref_model.optimize()
    
    print(f"  상태: {ref_sol.status}")
    print(f"  성장률: {ref_sol.objective_value:.6f}")
    
    key_reactions = ['ACS_ADP', 'CS', 'ICL', 'MALS', 'SUCDi', 'ATPS4rpp']
    print(f"\n  주요 반응 플럭스:")
    for rxn_id in key_reactions:
        if rxn_id in ref_model.reactions:
            flux = ref_sol.fluxes.get(rxn_id, 0.0)
            if abs(flux) > 1e-6:
                print(f"    {rxn_id:20s}: {flux:>12.6f}")
    
    print(f"\n[신규 모델 FBA (레퍼런스 Exchange 플럭스 사용)]")
    new_model.objective = 'Growth'
    new_sol = new_model.optimize()
    
    print(f"  상태: {new_sol.status}")
    print(f"  성장률: {new_sol.objective_value:.6f}")
    
    print(f"\n  주요 반응 플럭스:")
    for rxn_id in key_reactions:
        if rxn_id in new_model.reactions:
            flux = new_sol.fluxes.get(rxn_id, 0.0)
            if abs(flux) > 1e-6:
                print(f"    {rxn_id:20s}: {flux:>12.6f}")
    
    # 비교
    print(f"\n[비교]")
    print(f"{'반응':<20} {'레퍼런스':<15} {'신규':<15} {'차이':<15}")
    print("-" * 65)
    
    for rxn_id in key_reactions:
        ref_flux = ref_sol.fluxes.get(rxn_id, 0.0) if rxn_id in ref_model.reactions else None
        new_flux = new_sol.fluxes.get(rxn_id, 0.0) if rxn_id in new_model.reactions else None
        
        if ref_flux is not None and new_flux is not None:
            diff = abs(ref_flux - new_flux)
            if abs(ref_flux) > 1e-6 or abs(new_flux) > 1e-6:
                print(f"{rxn_id:<20} {ref_flux:>15.6f} {new_flux:>15.6f} {diff:>15.6f}")
        elif ref_flux is None:
            print(f"{rxn_id:<20} {'없음':<15} {new_flux:>15.6f} {'N/A':<15}")
        elif new_flux is None:
            print(f"{rxn_id:<20} {ref_flux:>15.6f} {'없음':<15} {'N/A':<15}")

def main():
    reproduce_reference_settings()
    
    print("\n" + "="*80)
    print("결론")
    print("="*80)
    print(f"\n[핵심 질문]")
    print(f"  레퍼런스 모델에서는 FBA가 잘 돌아가는데, 신규 모델에서는 왜 안 될까?")
    print(f"  -> 레퍼런스 FBA 결과의 정확한 설정을 재현하여 비교 완료")

if __name__ == "__main__":
    main()
