#!/usr/bin/env python
"""
레퍼런스 미디어로 ACS 테스트
- CoA 합성 전구체 포함 (pnto__R_e, ncam_e 등)
- ACS + ADK1 조합이 작동하는지 확인
"""

import cobra
from pathlib import Path
import csv

def load_model(model_path):
    model = cobra.io.read_sbml_model(model_path)
    return model

def apply_media_from_tsv(model, tsv_path):
    """TSV 파일에서 미디어 적용 (레퍼런스 스크립트 방식)"""
    with open(tsv_path, "r", encoding="utf-8-sig", newline="") as f:
        raw_lines = f.readlines()
    
    # 주석 및 빈 줄 제거
    lines = []
    for ln in raw_lines:
        ln_stripped = ln.strip()
        if ln_stripped and not ln_stripped.startswith("#"):
            lines.append(ln)
    
    rdr = csv.DictReader(lines, delimiter="\t")
    
    # 헤더 확인
    fieldnames = list(rdr.fieldnames)
    if len(fieldnames) >= 3 and fieldnames[1].replace(".", "").replace("-", "").isdigit():
        # normalized 파일 형식
        col_ex = fieldnames[0]
        col_lb = fieldnames[1]
        col_ub = fieldnames[2]
    else:
        # 일반 형식
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

def test_with_reference_media(model, media_path):
    """레퍼런스 미디어로 테스트"""
    print("="*80)
    print("레퍼런스 미디어로 ACS 테스트")
    print("="*80)
    
    # 모든 exchange 초기화
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    # 레퍼런스 미디어 적용
    try:
        applied = apply_media_from_tsv(model, media_path)
        print(f"\n[미디어 적용] {applied}개 반응 설정")
    except Exception as e:
        print(f"\n[미디어 적용 실패] {e}")
        return None
    
    # ATPM=0 설정
    atpm_rxn = model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    model.objective = 'Growth'
    
    # FBA 수행
    solution = model.optimize()
    print(f"\n[FBA 결과]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    
    if solution.status == 'optimal':
        # 주요 반응 플럭스
        acs_flux = solution.fluxes.get('ACS', 0.0)
        adk1_flux = solution.fluxes.get('ADK1', 0.0)
        icl_flux = solution.fluxes.get('ICL', 0.0)
        mals_flux = solution.fluxes.get('MALS', 0.0)
        icdhx_flux = solution.fluxes.get('ICDHx', 0.0)
        cs_flux = solution.fluxes.get('CS', 0.0)
        
        print(f"\n[주요 반응 플럭스]")
        print(f"  ACS: {acs_flux:.6f}")
        print(f"  ADK1: {adk1_flux:.6f}")
        print(f"  ICL: {icl_flux:.6f}")
        print(f"  MALS: {mals_flux:.6f}")
        print(f"  ICDHx: {icdhx_flux:.6f}")
        print(f"  CS: {cs_flux:.6f}")
        
        # CoA 관련 플럭스 확인
        print(f"\n[CoA 관련 Exchange 플럭스]")
        coa_exchanges = ['EX_pnto__R_e', 'EX_ncam_e', 'EX_coa_e']
        for ex_id in coa_exchanges:
            try:
                ex_rxn = model.reactions.get_by_id(ex_id)
                flux = solution.fluxes.get(ex_id, 0.0)
                bounds = f"[{ex_rxn.lower_bound}, {ex_rxn.upper_bound}]"
                print(f"  {ex_id}: 플럭스={flux:.6f}, bounds={bounds}")
            except KeyError:
                print(f"  {ex_id}: 반응 없음")
        
        if abs(acs_flux) > 1e-6:
            print("\n[OK] ACS 반응이 작동합니다!")
            if abs(adk1_flux) > 1e-6:
                print(f"  -> ADK1도 작동 중 (플럭스: {adk1_flux:.6f})")
                print(f"     ACS + ADK1 조합으로 작동!")
            
            if solution.objective_value > 1e-6:
                print(f"\n[성공] 성장 가능! (성장률: {solution.objective_value:.6f})")
            else:
                print(f"\n[부분 성공] ACS는 작동하지만 성장률이 0")
                print(f"  -> 다른 반응이 필요할 수 있음")
        else:
            print("\n[문제] ACS 반응이 여전히 작동하지 않음")
            print(f"  -> 레퍼런스 미디어로도 해결되지 않음")
    
    return solution

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    media_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5" / "Acetate_YE0p5__nocmnt__normalized.tsv"
    
    print("="*80)
    print("레퍼런스 미디어로 ACS 테스트")
    print("="*80)
    
    print("\n[사용자 지적]")
    print("  ACS + ADK1 조합으로 작동 가능해야 함 (에너지 효율이 낮아도 OK)")
    print("  -> 에너지 효율이 낮으면 TCA cycle을 통해 보상 (NADH 생성, 탄소 손실)")
    
    model = load_model(str(model_path))
    print(f"\n[모델 로드 완료]")
    
    if media_path.exists():
        print(f"[레퍼런스 미디어] {media_path}")
        print(f"  -> CoA 합성 전구체 포함 (pnto__R_e, ncam_e 등)")
        
        solution = test_with_reference_media(model, media_path)
        
        print("\n" + "="*80)
        print("결론")
        print("="*80)
        
        if solution:
            acs_flux = solution.fluxes.get('ACS', 0.0)
            growth = solution.objective_value
            
            if abs(acs_flux) > 1e-6:
                print("\n[성공] ACS + ADK1 조합으로 작동!")
                print(f"  ACS 플럭스: {acs_flux:.6f}")
                print(f"  사용자 지적이 맞습니다!")
                
                if growth > 1e-6:
                    print(f"\n[완전 성공] 성장도 가능! (성장률: {growth:.6f})")
                else:
                    print(f"\n[부분 성공] ACS는 작동하지만 성장률이 0")
                    print(f"  -> 추가 반응이 필요할 수 있음")
            else:
                print("\n[실패] 레퍼런스 미디어로도 ACS가 작동하지 않음")
                print(f"  -> 다른 원인이 있을 수 있음")
    else:
        print(f"\n[레퍼런스 미디어 파일 없음] {media_path}")

if __name__ == "__main__":
    main()
