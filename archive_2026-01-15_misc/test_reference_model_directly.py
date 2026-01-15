#!/usr/bin/env python
"""
레퍼런스 모델을 직접 테스트
- 레퍼런스 모델을 그대로 사용하여 FBA 실행
- 어떤 설정에서 작동하는지 확인
"""

import cobra
from pathlib import Path
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

def test_reference_model():
    """레퍼런스 모델 직접 테스트"""
    base_path = Path(__file__).parent.parent
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    media_tsv = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5" / "Acetate_YE0p5__nocmnt__normalized.tsv"
    
    print("="*80)
    print("레퍼런스 모델 직접 테스트")
    print("="*80)
    
    # 모델 로드
    ref_model = load_model(str(ref_model_path))
    
    print(f"\n[모델 정보]")
    print(f"  반응 수: {len(ref_model.reactions)}")
    print(f"  메타볼라이트 수: {len(ref_model.metabolites)}")
    
    # 기본 ATPM 확인
    atpm_rxn = ref_model.reactions.get_by_id('ATPM')
    print(f"\n[기본 ATPM 설정]")
    print(f"  bounds: [{atpm_rxn.lower_bound}, {atpm_rxn.upper_bound}]")
    print(f"  반응식: {atpm_rxn.reaction}")
    
    # 미디어 설정
    print(f"\n[미디어 설정]")
    for rxn in ref_model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    if media_tsv.exists():
        applied = apply_media_from_tsv(ref_model, media_tsv)
        print(f"  TSV 미디어 적용: {applied}개 exchange")
    
    # Exchange bounds 확인
    key_exchanges = ['EX_ac_e', 'EX_o2_e', 'EX_hco3_e', 'EX_nh4_e', 'EX_pi_e']
    print(f"\n[주요 Exchange bounds]")
    for ex_id in key_exchanges:
        if ex_id in ref_model.reactions:
            rxn = ref_model.reactions.get_by_id(ex_id)
            print(f"  {ex_id:20s}: [{rxn.lower_bound:>8.1f}, {rxn.upper_bound:>8.1f}]")
    
    # ATPM=0으로 설정하고 테스트
    print(f"\n[ATPM=0으로 설정하고 테스트]")
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    ref_model.objective = 'Growth'
    solution = ref_model.optimize()
    
    print(f"\n[FBA 결과]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    
    key_reactions = ['ACS_ADP', 'ACS', 'CS', 'ICL', 'MALS', 'SUCDi', 'ATPS4rpp', 'NADH16pp']
    print(f"\n[주요 반응 플럭스]")
    active_count = 0
    for rxn_id in key_reactions:
        if rxn_id in ref_model.reactions:
            flux = solution.fluxes.get(rxn_id, 0.0)
            if abs(flux) > 1e-6:
                print(f"  {rxn_id:20s}: {flux:>12.6f}")
                active_count += 1
    
    if active_count == 0:
        print(f"  [문제] 모든 주요 반응 플럭스가 0")
        print(f"  -> 레퍼런스 모델도 이 설정으로는 작동하지 않음")
    else:
        print(f"  [OK] {active_count}개 반응이 작동함")
    
    # Exchange 플럭스 확인
    print(f"\n[Exchange 플럭스]")
    for ex_id in key_exchanges:
        if ex_id in ref_model.reactions:
            flux = solution.fluxes.get(ex_id, 0.0)
            if abs(flux) > 1e-6:
                print(f"  {ex_id:20s}: {flux:>12.6f}")

def main():
    test_reference_model()
    
    print("\n" + "="*80)
    print("결론")
    print("="*80)
    print(f"\n[핵심 질문]")
    print(f"  레퍼런스 모델에서는 FBA가 잘 돌아가는데, 신규 모델에서는 왜 안 될까?")
    print(f"  -> 레퍼런스 모델 자체를 테스트하여 확인 완료")

if __name__ == "__main__":
    main()
