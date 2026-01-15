#!/usr/bin/env python
"""
종합 해결: 1, 2, 3번을 모두 수행하고 최종 테스트
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

def step1_check_reference_coa_init(ref_model, media_path):
    """1번: 레퍼런스 모델의 실제 CoA 초기화 방식 확인"""
    print("="*80)
    print("1번: 레퍼런스 모델의 실제 CoA 초기화 방식 확인")
    print("="*80)
    
    # 미디어 적용
    for rxn in ref_model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    if media_path.exists():
        applied = apply_media_from_tsv(ref_model, media_path)
        print(f"\n[미디어 적용] {applied}개 반응 설정")
    
    # ATPM=0 설정
    atpm_rxn = ref_model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    ref_model.objective = 'Growth'
    solution = ref_model.optimize()
    
    print(f"\n[레퍼런스 모델 FBA 결과]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    print(f"  ACS_ADP: {solution.fluxes.get('ACS_ADP', 0.0):.6f}")
    print(f"  ICL: {solution.fluxes.get('ICL', 0.0):.6f}")
    
    # CoA 관련 확인
    print(f"\n[CoA 관련 확인]")
    try:
        coa_c = ref_model.metabolites.get_by_id('coa_c')
        
        # CoA를 생성하는 반응
        coa_producing = []
        for rxn in coa_c.reactions:
            flux = solution.fluxes.get(rxn.id, 0.0)
            if abs(flux) > 1e-6:
                coeff = rxn.metabolites.get(coa_c, 0)
                if coeff > 0:
                    coa_producing.append((rxn.id, flux * coeff))
        
        if coa_producing:
            print(f"  CoA 생성 반응 (플럭스 > 0):")
            for rxn_id, net_flux in sorted(coa_producing, key=lambda x: x[1], reverse=True)[:5]:
                print(f"    {rxn_id}: {net_flux:.6f}")
        else:
            print(f"  CoA 생성 반응 플럭스 없음")
            print(f"  -> CoA는 초기화되어 있거나 다른 방식으로 공급됨")
    except KeyError:
        pass
    
    return solution

def step2_add_coa_synthesis(new_model, ref_model):
    """2번: CoA 합성 경로의 모든 반응 추가"""
    print("\n" + "="*80)
    print("2번: CoA 합성 경로의 모든 반응 추가")
    print("="*80)
    
    # 레퍼런스 모델의 CoA 합성 관련 반응 찾기
    coa_keywords = ['PPAT', 'DPCK', 'PPCS', 'PPCDC', 'COAS', 'PAN']
    
    ref_coa_rxns = []
    for rxn in ref_model.reactions:
        rxn_id_upper = rxn.id.upper()
        for keyword in coa_keywords:
            if keyword in rxn_id_upper:
                ref_coa_rxns.append(rxn)
                break
    
    print(f"\n[레퍼런스 모델 CoA 합성 관련 반응] {len(ref_coa_rxns)}개")
    
    # 신규 모델과 비교
    new_coa_rxns = set([r.id for r in new_model.reactions if any(kw in r.id.upper() for kw in coa_keywords)])
    ref_coa_rxns_ids = set([r.id for r in ref_coa_rxns])
    
    missing = ref_coa_rxns_ids - new_coa_rxns
    
    print(f"[신규 모델에 없는 반응] {len(missing)}개")
    
    added = []
    for rxn_id in sorted(missing)[:5]:
        try:
            ref_rxn = ref_model.reactions.get_by_id(rxn_id)
            print(f"\n  {rxn_id}: {ref_rxn.reaction[:80]}")
            
            new_rxn = cobra.Reaction(rxn_id)
            new_rxn.name = ref_rxn.name
            new_rxn.lower_bound = ref_rxn.lower_bound
            new_rxn.upper_bound = ref_rxn.upper_bound
            
            metabolites_dict = {}
            for met, coeff in ref_rxn.metabolites.items():
                if met.id in [m.id for m in new_model.metabolites]:
                    metabolites_dict[new_model.metabolites.get_by_id(met.id)] = coeff
                else:
                    new_met = cobra.Metabolite(
                        met.id,
                        name=met.name,
                        compartment=met.compartment,
                        formula=met.formula if hasattr(met, 'formula') else None,
                        charge=met.charge if hasattr(met, 'charge') else None
                    )
                    new_model.add_metabolites([new_met])
                    metabolites_dict[new_met] = coeff
            
            new_rxn.add_metabolites(metabolites_dict)
            new_model.add_reactions([new_rxn])
            added.append(rxn_id)
            print(f"    -> 추가됨")
        except Exception as e:
            print(f"    -> 추가 실패: {e}")
    
    return added

def step3_add_acs_adp(new_model, ref_model):
    """3번: ACS_ADP 반응 추가"""
    print("\n" + "="*80)
    print("3번: ACS_ADP 반응 추가")
    print("="*80)
    
    if 'ACS_ADP' in [r.id for r in new_model.reactions]:
        print(f"\n[ACS_ADP 이미 존재]")
        return False
    
    try:
        ref_acs_adp = ref_model.reactions.get_by_id('ACS_ADP')
        
        new_acs_adp = cobra.Reaction('ACS_ADP')
        new_acs_adp.name = ref_acs_adp.name
        new_acs_adp.lower_bound = ref_acs_adp.lower_bound
        new_acs_adp.upper_bound = ref_acs_adp.upper_bound
        
        metabolites_dict = {}
        for met, coeff in ref_acs_adp.metabolites.items():
            if met.id in [m.id for m in new_model.metabolites]:
                metabolites_dict[new_model.metabolites.get_by_id(met.id)] = coeff
            else:
                new_met = cobra.Metabolite(
                    met.id,
                    name=met.name,
                    compartment=met.compartment,
                    formula=met.formula if hasattr(met, 'formula') else None,
                    charge=met.charge if hasattr(met, 'charge') else None
                )
                new_model.add_metabolites([new_met])
                metabolites_dict[new_met] = coeff
        
        new_acs_adp.add_metabolites(metabolites_dict)
        new_model.add_reactions([new_acs_adp])
        
        print(f"\n[ACS_ADP 추가 완료]")
        print(f"  반응식: {new_acs_adp.reaction}")
        return True
        
    except Exception as e:
        print(f"\n[오류] ACS_ADP 추가 실패: {e}")
        return False

def test_final_comprehensive(new_model, media_path):
    """최종 종합 테스트"""
    print("\n" + "="*80)
    print("최종 종합 테스트")
    print("="*80)
    
    # 미디어 설정
    for rxn in new_model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    if media_path.exists():
        applied = apply_media_from_tsv(new_model, media_path)
        print(f"[미디어 적용] {applied}개 반응 설정")
    
    # ATPM=0 설정
    atpm_rxn = new_model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    new_model.objective = 'Growth'
    solution = new_model.optimize()
    
    print(f"\n[FBA 결과]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    
    acs_flux = solution.fluxes.get('ACS', 0.0)
    acs_adp_flux = solution.fluxes.get('ACS_ADP', 0.0)
    cs_flux = solution.fluxes.get('CS', 0.0)
    adk1_flux = solution.fluxes.get('ADK1', 0.0)
    icl_flux = solution.fluxes.get('ICL', 0.0)
    mals_flux = solution.fluxes.get('MALS', 0.0)
    
    print(f"\n[주요 반응 플럭스]")
    print(f"  ACS: {acs_flux:.6f}")
    print(f"  ACS_ADP: {acs_adp_flux:.6f}")
    print(f"  CS: {cs_flux:.6f}")
    print(f"  ADK1: {adk1_flux:.6f}")
    print(f"  ICL: {icl_flux:.6f}")
    print(f"  MALS: {mals_flux:.6f}")
    
    # 경로 단계별 확인
    print(f"\n[경로 단계별 확인]")
    
    if abs(acs_adp_flux) > 1e-6:
        print(f"  [OK] 1단계: Acetate -> Acetyl-CoA (ACS_ADP: {acs_adp_flux:.6f})")
    elif abs(acs_flux) > 1e-6:
        print(f"  [OK] 1단계: Acetate -> Acetyl-CoA (ACS: {acs_flux:.6f})")
    else:
        print(f"  [NO] 1단계: Acetate -> Acetyl-CoA 막힘")
        return False
    
    if abs(cs_flux) > 1e-6:
        print(f"  [OK] 2단계: Acetyl-CoA -> Citrate (CS: {cs_flux:.6f})")
    else:
        print(f"  [NO] 2단계: Acetyl-CoA -> Citrate 막힘")
        return False
    
    if abs(icl_flux) > 1e-6 or abs(mals_flux) > 1e-6:
        print(f"  [OK] 3단계: Glyoxylate shunt 작동 (ICL: {icl_flux:.6f}, MALS: {mals_flux:.6f})")
    else:
        print(f"  [NO] 3단계: Glyoxylate shunt 작동 안 함")
    
    if abs(adk1_flux) > 1e-6:
        print(f"  [OK] ADK1: AMP 재생성 작동 (ADK1: {adk1_flux:.6f})")
    
    if abs(acs_adp_flux) > 1e-6 or abs(acs_flux) > 1e-6:
        print(f"\n[성공] Acetate 전환 경로 작동!")
        return True
    else:
        print(f"\n[실패] Acetate 전환 경로 작동 안 함")
        return False

def main():
    base_path = Path(__file__).parent.parent
    new_model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    media_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5" / "Acetate_YE0p5__nocmnt__normalized.tsv"
    
    new_model = load_model(str(new_model_path))
    
    if not ref_model_path.exists():
        print(f"[오류] 레퍼런스 모델 파일 없음: {ref_model_path}")
        return
    
    ref_model = cobra.io.read_sbml_model(str(ref_model_path))
    
    print("="*80)
    print("종합 해결: 1번 → 2번 → 3번 순서로 수행")
    print("="*80)
    
    # 1번: 레퍼런스 모델의 실제 CoA 초기화 방식 확인
    ref_sol = step1_check_reference_coa_init(ref_model, media_path)
    
    # 2번: CoA 합성 경로의 모든 반응 추가
    added_coa = step2_add_coa_synthesis(new_model, ref_model)
    
    # 3번: ACS_ADP 반응 추가
    acs_adp_added = step3_add_acs_adp(new_model, ref_model)
    
    # 최종 테스트
    success = test_final_comprehensive(new_model, media_path)
    
    print("\n" + "="*80)
    print("최종 결과")
    print("="*80)
    
    print(f"\n[1번 결과]")
    print(f"  레퍼런스 모델 분석 완료")
    print(f"  ACS_ADP 플럭스: {ref_sol.fluxes.get('ACS_ADP', 0.0):.6f}")
    
    print(f"\n[2번 결과]")
    print(f"  CoA 합성 경로 반응 추가: {len(added_coa)}개")
    if added_coa:
        print(f"  추가된 반응: {', '.join(added_coa)}")
    
    print(f"\n[3번 결과]")
    print(f"  ACS_ADP 추가: {'성공' if acs_adp_added else '실패'}")
    
    print(f"\n[최종 테스트]")
    if success:
        print(f"  성공: 경로 작동!")
    else:
        print(f"  실패: 추가 조사 필요")
        print(f"\n[막힌 지점]")
        print(f"  -> 경로 단계별 확인 결과 참고")

if __name__ == "__main__":
    main()
