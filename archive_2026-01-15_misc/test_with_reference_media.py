#!/usr/bin/env python
"""
레퍼런스 모델과 동일한 미디어 설정으로 FBA 테스트
레퍼런스 모델이 사용한 정확한 미디어 설정을 적용
"""

import cobra
from pathlib import Path
import pandas as pd
import csv
import sys

def load_model(model_path):
    """모델 로드"""
    print(f"모델 로드 중: {model_path}")
    model = cobra.io.read_sbml_model(str(model_path))
    print(f"[OK] 모델 로드 완료: {model.id}")
    return model

def add_reaction_to_model(model, ref_model, rxn_id):
    """레퍼런스 모델에서 반응을 신규 모델에 추가"""
    if rxn_id in model.reactions:
        return False
    
    if rxn_id not in ref_model.reactions:
        return False
    
    ref_rxn = ref_model.reactions.get_by_id(rxn_id)
    
    # 대사물질 먼저 확인/추가
    metabolites_dict = {}
    for met in ref_rxn.metabolites:
        if met.id not in model.metabolites:
            new_met = cobra.Metabolite(met.id, formula=met.formula, name=met.name, 
                                      compartment=met.compartment, charge=met.charge)
            model.add_metabolites([new_met])
        metabolites_dict[model.metabolites.get_by_id(met.id)] = ref_rxn.metabolites[met]
    
    # 반응 생성 및 추가
    new_rxn = cobra.Reaction(rxn_id)
    new_rxn.name = ref_rxn.name
    new_rxn.lower_bound = ref_rxn.lower_bound
    new_rxn.upper_bound = ref_rxn.upper_bound
    new_rxn.add_metabolites(metabolites_dict)
    
    model.add_reactions([new_rxn])
    
    # 유전자 연결
    for gene in ref_rxn.genes:
        if gene.id in model.genes:
            model.genes.get_by_id(gene.id).reactions.add(new_rxn)
    
    return True

def apply_media_from_tsv(model, tsv_path):
    """레퍼런스 모델의 미디어 TSV 파일을 적용"""
    print(f"\n미디어 적용 중: {tsv_path}")
    
    p = Path(tsv_path)
    if not p.exists():
        raise FileNotFoundError(f"Media TSV not found: {p}")
    
    def _clean_text(s):
        if s is None:
            return ""
        return str(s).lstrip("\ufeff").strip()
    
    def _clean_num(s):
        if s is None or s == "":
            return None
        try:
            return float(str(s).strip())
        except:
            return None
    
    def _resolve_cols(fieldnames):
        hdr = [_clean_text(h).lower() for h in fieldnames]
        def pick(key, synonyms):
            for syn in synonyms:
                if syn in hdr:
                    return hdr.index(syn)
            return None
        
        i_id = pick("exchange_id", ["exchange_id", "exchange", "id", "rxn", "rxn_id"])
        i_lb = pick("minflux", ["minflux", "lb", "lower", "lower_bound", "0.0"])
        i_ub = pick("maxflux", ["maxflux", "ub", "upper", "upper_bound", "0.0"])
        
        if i_id is None:
            raise ValueError("exchange_id column not found")
        if i_lb is None or i_ub is None:
            raise ValueError("minflux/maxflux columns not found")
        
        return {"exchange_id": i_id, "minflux": i_lb, "maxflux": i_ub}
    
    with open(p, "r", encoding="utf-8-sig", newline="") as f:
        raw = f.readlines()
    
    lines = []
    for ln in raw:
        if not ln or ln.strip() == "":
            continue
        if ln.lstrip().startswith("#"):
            continue
        lines.append(ln)
    
    rdr = csv.DictReader(lines, delimiter="\t")
    if not rdr.fieldnames:
        raise ValueError(f"Media TSV has no header: {p}")
    
    cols = _resolve_cols(rdr.fieldnames)
    n = 0
    missing = []
    
    for i, row in enumerate(rdr, start=2):
        ex_id = _clean_text(row.get(list(rdr.fieldnames)[cols["exchange_id"]], ""))
        if ex_id == "" or ex_id.lower() == "exchange_id":
            continue
        
        # 컬럼 인덱스로 접근
        fieldnames_list = list(rdr.fieldnames)
        lb_val = row.get(fieldnames_list[cols["minflux"]], None)
        ub_val = row.get(fieldnames_list[cols["maxflux"]], None)
        
        lb = _clean_num(lb_val)
        ub = _clean_num(ub_val)
        
        if lb is None and ub is None:
            continue
        if lb is None:
            lb = 0.0
        if ub is None:
            ub = 0.0
        
        if ex_id not in model.reactions:
            missing.append(ex_id)
            continue
        
        rxn = model.reactions.get_by_id(ex_id)
        rxn.lower_bound = lb
        rxn.upper_bound = ub
        n += 1
    
    print(f"  {n}개 exchange 반응 설정 완료")
    if missing:
        print(f"  [WARN] 모델에 없는 exchange: {len(missing)}개 (처음 10개: {missing[:10]})")
    
    return model

def test_fba(model):
    """FBA 테스트"""
    # Biomass 반응 찾기
    biomass_rxns = [r for r in model.reactions if 'growth' in r.id.lower() or 'biomass' in r.id.lower()]
    if not biomass_rxns:
        print("[ERROR] Biomass 반응을 찾을 수 없습니다")
        return None, "NO_BIOMASS"
    
    biomass_rxn = biomass_rxns[0]
    
    # 목적함수 설정 (레퍼런스와 동일하게)
    for r in model.reactions:
        try:
            r.objective_coefficient = 0.0
        except:
            pass
    
    model.reactions.get_by_id(biomass_rxn.id).objective_coefficient = 1.0
    model.objective = biomass_rxn.id
    try:
        model.objective_direction = "max"
        model.solver.objective.direction = "max"
    except:
        pass
    
    print(f"목적함수: {biomass_rxn.id} (최대화)")
    
    # FBA 실행
    try:
        solution = model.optimize()
        print(f"\nFBA 상태: {solution.status}")
        print(f"목적함수 값 (성장률): {solution.objective_value}")
        
        if solution.status == 'optimal' and solution.objective_value and solution.objective_value > 1e-6:
            return solution, "SUCCESS"
        elif solution.status == 'optimal':
            return solution, "OPTIMAL_BUT_ZERO"
        else:
            return solution, solution.status
    except Exception as e:
        print(f"[ERROR] FBA 실행 중 오류: {e}")
        return None, str(e)

def main():
    base_path = Path(__file__).parent.parent
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    new_model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    media_tsv = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5" / "Acetate_YE0p5__nocmnt__normalized.tsv"
    
    print("="*70)
    print("레퍼런스 모델과 동일한 미디어 설정으로 FBA 테스트")
    print("="*70)
    
    # 모델 로드
    ref_model = load_model(ref_model_path)
    new_model = load_model(new_model_path)
    
    # 누락된 9개 반응 추가 (실제 사용된 반응)
    missing_reactions = [
        'ACS_ADP', 'PEPCK_ATP', 'SUCDi', 'ACtexi',
        'EX_hco3_e', 'T_hco3_e_to_c', 'T_o2_e_to_o2_c',
        'T_nh4_e_to_nh4_c', 'T_fe3_e_to_fe3_c'
    ]
    
    print("\n" + "="*70)
    print("누락된 9개 반응 추가")
    print("="*70)
    
    added_count = 0
    for rxn_id in missing_reactions:
        if add_reaction_to_model(new_model, ref_model, rxn_id):
            added_count += 1
            print(f"  [ADDED] {rxn_id}")
    
    print(f"\n총 {added_count}/{len(missing_reactions)}개 반응 추가 완료")
    
    # 레퍼런스 모델의 미디어 설정 적용
    try:
        new_model = apply_media_from_tsv(new_model, media_tsv)
    except Exception as e:
        print(f"[ERROR] 미디어 적용 실패: {e}")
        print("기본 미디어 설정으로 대체...")
        # 기본 설정
        if "EX_ac_e" in new_model.reactions:
            new_model.reactions.get_by_id("EX_ac_e").bounds = (-100.0, 0.0)
        if "EX_o2_e" in new_model.reactions:
            new_model.reactions.get_by_id("EX_o2_e").bounds = (-200.0, 0.0)
    
    # FBA 테스트
    print("\n" + "="*70)
    print("FBA 테스트")
    print("="*70)
    
    solution, status = test_fba(new_model)
    
    # 결과 요약
    print("\n" + "="*70)
    print("결과 요약")
    print("="*70)
    print(f"추가된 반응 수: {added_count}/{len(missing_reactions)}")
    print(f"FBA 상태: {status}")
    
    if solution:
        print(f"성장률: {solution.objective_value:.6f} 1/h")
        
        if status == "SUCCESS":
            print("\n[SUCCESS] 레퍼런스 미디어 설정으로 FBA가 성공했습니다!")
            
            # 주요 반응 플럭스 확인
            print("\n주요 반응 플럭스:")
            key_reactions = ['ACS_ADP', 'PEPCK_ATP', 'SUCDi', 'ACtexi', 'EX_ac_e', 'EX_o2_e']
            for rxn_id in key_reactions:
                if rxn_id in new_model.reactions:
                    flux = solution.fluxes.get(rxn_id, 0)
                    if abs(flux) > 1e-6:
                        print(f"  {rxn_id:20s}: {flux:>12.6f}")
        else:
            print("\n[FAIL] 레퍼런스 미디어 설정으로도 FBA가 실패합니다.")
    else:
        print("성장률: N/A")
        print("\n[FAIL] FBA 실행 실패")
    
    print("\n" + "="*70)

if __name__ == "__main__":
    main()
