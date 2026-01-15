#!/usr/bin/env python
"""
레퍼런스 모델의 실제 FBA 설정 재현
- fba_flux_gradient_acid.csv와 동일한 설정으로 재현
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

def test_with_reference_settings(model, media_path, atpm_value=0):
    """레퍼런스 설정으로 테스트"""
    # 모든 exchange 초기화
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    # 레퍼런스 미디어 적용
    applied = apply_media_from_tsv(model, media_path)
    
    # ATPM 설정
    atpm_rxn = model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = atpm_value
    atpm_rxn.upper_bound = 1000
    
    model.objective = 'Growth'
    solution = model.optimize()
    
    return solution

def step1_reproduce_reference(new_model, ref_model, media_path):
    """1번: 레퍼런스 모델의 실제 CoA 초기화 방식 확인 및 재현"""
    print("="*80)
    print("1번: 레퍼런스 모델의 실제 FBA 설정 재현")
    print("="*80)
    
    # 레퍼런스 모델로 테스트
    print(f"\n[레퍼런스 모델 테스트]")
    ref_solution = test_with_reference_settings(ref_model.copy(), media_path, atpm_value=0)
    
    print(f"  상태: {ref_solution.status}")
    print(f"  성장률: {ref_solution.objective_value:.6f}")
    print(f"  ACS_ADP: {ref_solution.fluxes.get('ACS_ADP', 0.0):.6f}")
    print(f"  ICL: {ref_solution.fluxes.get('ICL', 0.0):.6f}")
    print(f"  MALS: {ref_solution.fluxes.get('MALS', 0.0):.6f}")
    
    # 신규 모델로 테스트
    print(f"\n[신규 모델 테스트 (레퍼런스 설정)]")
    new_solution = test_with_reference_settings(new_model.copy(), media_path, atpm_value=0)
    
    print(f"  상태: {new_solution.status}")
    print(f"  성장률: {new_solution.objective_value:.6f}")
    print(f"  ACS_ADP: {new_solution.fluxes.get('ACS_ADP', 0.0):.6f}")
    print(f"  ACS: {new_solution.fluxes.get('ACS', 0.0):.6f}")
    print(f"  ICL: {new_solution.fluxes.get('ICL', 0.0):.6f}")
    print(f"  MALS: {new_solution.fluxes.get('MALS', 0.0):.6f}")
    
    # 차이점 확인
    print(f"\n[차이점]")
    if abs(ref_solution.fluxes.get('ACS_ADP', 0.0)) > 1e-6:
        if abs(new_solution.fluxes.get('ACS_ADP', 0.0)) < 1e-6:
            print(f"  레퍼런스: ACS_ADP 작동 (플럭스: {ref_solution.fluxes.get('ACS_ADP', 0.0):.6f})")
            print(f"  신규: ACS_ADP 작동 안 함")
            print(f"  -> ACS_ADP가 신규 모델에 없거나 작동하지 않음")
    
    return ref_solution, new_solution

def step2_add_coa_synthesis_pathway(new_model, ref_model):
    """2번: CoA 합성 경로의 모든 반응 추가"""
    print("\n" + "="*80)
    print("2번: CoA 합성 경로의 모든 반응 추가")
    print("="*80)
    
    # CoA 합성 경로 전체 확인
    print(f"\n[CoA 합성 경로 전체 확인]")
    
    # 표준 CoA 합성 경로
    # 1. Pantothenate → 4'-Phosphopantothenate
    # 2. 4'-Phosphopantothenate + Cysteine → 4'-Phosphopantothenoylcysteine
    # 3. 4'-Phosphopantothenoylcysteine → 4'-Phosphopantetheine
    # 4. 4'-Phosphopantetheine + ATP → Dephospho-CoA
    # 5. Dephospho-CoA + ATP → CoA
    
    coa_pathway_keywords = [
        'PPAT',  # Phosphopantetheine adenylyltransferase
        'DPCK',  # Dephospho-CoA kinase
        'PPCS',  # Phosphopantothenate--cysteine ligase
        'PPCDC', # Phosphopantothenoylcysteine decarboxylase
        'COAS',  # CoA synthase
        'PAN',   # Pantothenate 관련
    ]
    
    # 레퍼런스 모델에서 CoA 합성 관련 반응 찾기
    ref_coa_rxns = []
    for rxn in ref_model.reactions:
        rxn_id_upper = rxn.id.upper()
        for keyword in coa_pathway_keywords:
            if keyword in rxn_id_upper:
                ref_coa_rxns.append(rxn)
                break
    
    print(f"  레퍼런스 모델 CoA 합성 관련 반응: {len(ref_coa_rxns)}개")
    for rxn in ref_coa_rxns[:10]:
        print(f"    {rxn.id}: {rxn.reaction[:80]}")
    
    # 신규 모델과 비교
    new_coa_rxns = []
    for rxn in new_model.reactions:
        rxn_id_upper = rxn.id.upper()
        for keyword in coa_pathway_keywords:
            if keyword in rxn_id_upper:
                new_coa_rxns.append(rxn.id)
                break
    
    missing_rxns = [r for r in ref_coa_rxns if r.id not in new_coa_rxns]
    
    print(f"\n[신규 모델에 없는 CoA 합성 반응] {len(missing_rxns)}개")
    
    added = []
    for ref_rxn in missing_rxns[:5]:  # 상위 5개만
        print(f"\n  {ref_rxn.id}: {ref_rxn.reaction[:80]}")
        
        try:
            new_rxn = cobra.Reaction(ref_rxn.id)
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
            added.append(ref_rxn.id)
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
        print(f"\n[레퍼런스 모델 ACS_ADP]")
        print(f"  반응식: {ref_acs_adp.reaction}")
        print(f"  bounds: [{ref_acs_adp.lower_bound}, {ref_acs_adp.upper_bound}]")
        
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
        return True
        
    except Exception as e:
        print(f"\n[오류] ACS_ADP 추가 실패: {e}")
        return False

def test_final(new_model, media_path):
    """최종 테스트"""
    print("\n" + "="*80)
    print("최종 테스트 (레퍼런스 설정)")
    print("="*80)
    
    solution = test_with_reference_settings(new_model.copy(), media_path, atpm_value=0)
    
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
    
    if abs(acs_adp_flux) > 1e-6:
        print(f"\n[성공] ACS_ADP 작동!")
        if abs(icl_flux) > 1e-6 or abs(mals_flux) > 1e-6:
            print(f"  Glyoxylate shunt도 작동!")
        return True
    elif abs(acs_flux) > 1e-6:
        print(f"\n[부분 성공] ACS 작동!")
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
    print("CoA 문제 해결: 1번 → 2번 → 3번 순서로 수행")
    print("="*80)
    
    # 1번: 레퍼런스 모델의 실제 FBA 설정 재현
    ref_sol, new_sol = step1_reproduce_reference(new_model, ref_model, media_path)
    
    # 2번: CoA 합성 경로의 모든 반응 추가
    added_coa = step2_add_coa_synthesis_pathway(new_model, ref_model)
    
    # 3번: ACS_ADP 반응 추가
    acs_adp_added = step3_add_acs_adp(new_model, ref_model)
    
    # 최종 테스트
    success = test_final(new_model, media_path)
    
    print("\n" + "="*80)
    print("최종 결과")
    print("="*80)
    
    print(f"\n[1번 결과]")
    print(f"  레퍼런스 ACS_ADP 플럭스: {ref_sol.fluxes.get('ACS_ADP', 0.0):.6f}")
    print(f"  신규 ACS_ADP 플럭스: {new_sol.fluxes.get('ACS_ADP', 0.0):.6f}")
    
    print(f"\n[2번 결과]")
    print(f"  CoA 합성 경로 반응 추가: {len(added_coa)}개")
    
    print(f"\n[3번 결과]")
    print(f"  ACS_ADP 추가: {'성공' if acs_adp_added else '실패'}")
    
    print(f"\n[최종 테스트]")
    if success:
        print(f"  성공: 경로 작동!")
    else:
        print(f"  실패: 추가 조사 필요")

if __name__ == "__main__":
    main()
