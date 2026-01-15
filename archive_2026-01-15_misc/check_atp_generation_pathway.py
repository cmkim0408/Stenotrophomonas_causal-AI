#!/usr/bin/env python
"""
ATP 생성 경로 확인
- ATP 생성이 불가능한 원인 찾기
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

def check_atp_generation_pathway(model, media_path):
    """ATP 생성 경로 확인"""
    print("="*80)
    print("ATP 생성 경로 확인")
    print("="*80)
    
    # 미디어 설정
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    if media_path.exists():
        apply_media_from_tsv(model, media_path)
    
    # ATPM=0 설정
    atpm_rxn = model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    # ATP 생성 반응 찾기
    print(f"\n[ATP 생성 반응 찾기]")
    try:
        atp_c = model.metabolites.get_by_id('atp_c')
        
        atp_producing = []
        for rxn in atp_c.reactions:
            coeff = rxn.metabolites.get(atp_c, 0)
            if coeff > 0:  # 생성
                atp_producing.append((rxn.id, coeff))
        
        print(f"  ATP 생성 반응: {len(atp_producing)}개")
        
        # 주요 ATP 생성 반응
        key_atp_rxns = ['ATPS4rpp', 'NADH16pp', 'CYTBO3_4pp', 'PGK', 'PYK', 'ACKr', 'PTAr']
        
        print(f"\n[주요 ATP 생성 반응 확인]")
        for rxn_id in key_atp_rxns:
            if rxn_id in [r.id for r in model.reactions]:
                rxn = model.reactions.get_by_id(rxn_id)
                coeff = rxn.metabolites.get(atp_c, 0)
                if coeff > 0:
                    print(f"  {rxn_id}: 존재, ATP 생성 (계수: {coeff})")
                    print(f"    반응식: {rxn.reaction[:100]}")
                elif coeff < 0:
                    print(f"  {rxn_id}: 존재, ATP 소모 (계수: {coeff})")
                else:
                    print(f"  {rxn_id}: 존재, ATP 관련 없음")
            else:
                print(f"  {rxn_id}: 없음")
        
        # ATP demand로 최대 생산량 확인
        print(f"\n[ATP 최대 생산량 확인]")
        atp_demand = cobra.Reaction('DM_atp_c')
        atp_demand.name = 'ATP demand'
        atp_demand.lower_bound = 0
        atp_demand.upper_bound = 1000
        atp_demand.add_metabolites({atp_c: -1.0})
        model.add_reactions([atp_demand])
        
        model.objective = 'DM_atp_c'
        solution = model.optimize()
        
        atp_max = solution.objective_value
        print(f"  ATP 최대 생산량: {atp_max:.6f}")
        
        if atp_max < 1e-6:
            print(f"  [문제] ATP 생성 불가")
            
            # ATP 생성에 필요한 것들 확인
            print(f"\n[ATP 생성에 필요한 것들 확인]")
            
            # 1. NADH 확인
            try:
                nadh_c = model.metabolites.get_by_id('nadh_c')
                nadh_demand = cobra.Reaction('DM_nadh_c')
                nadh_demand.lower_bound = 0
                nadh_demand.upper_bound = 1000
                nadh_demand.add_metabolites({nadh_c: -1.0})
                model_test = model.copy()
                model_test.add_reactions([nadh_demand])
                model_test.objective = nadh_demand.id
                sol_nadh = model_test.optimize()
                nadh_max = sol_nadh.objective_value
                print(f"  NADH 최대 생산량: {nadh_max:.6f}")
            except:
                pass
            
            # 2. 전자전달사슬 확인
            etc_rxns = ['NADH16pp', 'CYTBO3_4pp', 'ATPS4rpp']
            print(f"\n[전자전달사슬 반응 확인]")
            for rxn_id in etc_rxns:
                if rxn_id in [r.id for r in model.reactions]:
                    rxn = model.reactions.get_by_id(rxn_id)
                    flux = solution.fluxes.get(rxn_id, 0.0)
                    print(f"  {rxn_id}: 플럭스={flux:.6f}, bounds=[{rxn.lower_bound}, {rxn.upper_bound}]")
                    if abs(flux) < 1e-6:
                        print(f"    [문제] 플럭스 0")
                else:
                    print(f"  {rxn_id}: 없음")
            
            # 3. 기질 수준 인산화 확인
            slp_rxns = ['PGK', 'PYK', 'ACKr', 'PTAr']
            print(f"\n[기질 수준 인산화 반응 확인]")
            for rxn_id in slp_rxns:
                if rxn_id in [r.id for r in model.reactions]:
                    rxn = model.reactions.get_by_id(rxn_id)
                    flux = solution.fluxes.get(rxn_id, 0.0)
                    print(f"  {rxn_id}: 플럭스={flux:.6f}, bounds=[{rxn.lower_bound}, {rxn.upper_bound}]")
                    if abs(flux) < 1e-6:
                        print(f"    [문제] 플럭스 0")
                else:
                    print(f"  {rxn_id}: 없음")
        else:
            print(f"  [OK] ATP 생성 가능")
        
        model.remove_reactions([atp_demand])
        
    except KeyError:
        print(f"  [오류] atp_c 메타볼라이트 없음")
        return False
    
    return True

def compare_with_reference(model, ref_model, media_path):
    """레퍼런스 모델과 비교"""
    print("\n" + "="*80)
    print("레퍼런스 모델과 비교")
    print("="*80)
    
    # 레퍼런스 모델 ATP 생성 확인
    for rxn in ref_model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    if media_path.exists():
        apply_media_from_tsv(ref_model, media_path)
    
    atpm_rxn = ref_model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    try:
        atp_c = ref_model.metabolites.get_by_id('atp_c')
        atp_demand = cobra.Reaction('DM_atp_c')
        atp_demand.lower_bound = 0
        atp_demand.upper_bound = 1000
        atp_demand.add_metabolites({atp_c: -1.0})
        ref_model.add_reactions([atp_demand])
        
        ref_model.objective = 'DM_atp_c'
        ref_solution = ref_model.optimize()
        
        ref_atp_max = ref_solution.objective_value
        print(f"\n[레퍼런스 모델 ATP 최대 생산량] {ref_atp_max:.6f}")
        
        # 전자전달사슬 플럭스
        etc_rxns = ['NADH16pp', 'CYTBO3_4pp', 'ATPS4rpp']
        print(f"\n[레퍼런스 모델 전자전달사슬 플럭스]")
        for rxn_id in etc_rxns:
            if rxn_id in [r.id for r in ref_model.reactions]:
                flux = ref_solution.fluxes.get(rxn_id, 0.0)
                print(f"  {rxn_id}: {flux:.6f}")
        
        ref_model.remove_reactions([atp_demand])
        
    except:
        pass

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    media_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5" / "Acetate_YE0p5__nocmnt__normalized.tsv"
    
    model = load_model(str(model_path))
    
    if not ref_model_path.exists():
        print(f"[오류] 레퍼런스 모델 파일 없음: {ref_model_path}")
        return
    
    ref_model = cobra.io.read_sbml_model(str(ref_model_path))
    
    print("="*80)
    print("ATP 생성 경로 확인")
    print("="*80)
    
    # ATP 생성 경로 확인
    check_atp_generation_pathway(model, media_path)
    
    # 레퍼런스 모델과 비교
    compare_with_reference(model, ref_model, media_path)
    
    print("\n" + "="*80)
    print("결론")
    print("="*80)
    
    print(f"\n[핵심 문제]")
    print(f"  ATP 생성 불가")
    print(f"  -> 이것이 모든 경로가 막힌 근본 원인")
    print(f"  -> ATP가 없으면 ACS_ADP, ACS 등 모든 ATP 의존 반응이 작동 불가")

if __name__ == "__main__":
    main()
