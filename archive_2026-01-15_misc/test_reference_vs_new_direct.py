#!/usr/bin/env python
"""
레퍼런스 모델과 신규 모델을 동일한 조건에서 직접 비교
레퍼런스 모델이 실제로 성장하는지 먼저 확인하고,
같은 조건으로 신규 모델 테스트
"""

import cobra
from pathlib import Path
import pandas as pd
import sys

def load_model(model_path):
    """모델 로드"""
    model = cobra.io.read_sbml_model(str(model_path))
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
    """미디어 TSV 파일 적용"""
    df = pd.read_csv(tsv_path, sep='\t', comment='#', skip_blank_lines=True)
    
    # 컬럼 찾기
    id_col = 'exchange_id' if 'exchange_id' in df.columns else df.columns[0]
    
    # minflux, maxflux 컬럼 찾기
    lb_col = None
    ub_col = None
    for col in df.columns:
        col_lower = str(col).lower()
        if 'min' in col_lower or 'lower' in col_lower or 'lb' in col_lower:
            lb_col = col
        elif 'max' in col_lower or 'upper' in col_lower or 'ub' in col_lower:
            ub_col = col
    
    if lb_col is None and len(df.columns) > 1:
        lb_col = df.columns[1]
    if ub_col is None and len(df.columns) > 2:
        ub_col = df.columns[2]
    
    n = 0
    for _, row in df.iterrows():
        ex_id = str(row[id_col]).strip()
        if ex_id == "" or ex_id.lower() == "exchange_id" or ex_id == "nan":
            continue
        
        try:
            lb = float(row[lb_col])
            ub = float(row[ub_col])
        except:
            continue
        
        if ex_id in model.reactions:
            model.reactions.get_by_id(ex_id).bounds = (lb, ub)
            n += 1
    
    return model, n

def test_fba_simple(model):
    """간단한 FBA 테스트"""
    biomass_rxns = [r for r in model.reactions if 'growth' in r.id.lower() or 'biomass' in r.id.lower()]
    if not biomass_rxns:
        return None, "NO_BIOMASS"
    
    biomass_rxn = biomass_rxns[0]
    model.objective = biomass_rxn.id
    
    try:
        solution = model.optimize()
        if solution.status == 'optimal' and solution.objective_value and solution.objective_value > 1e-6:
            return solution, "SUCCESS"
        elif solution.status == 'optimal':
            return solution, "OPTIMAL_BUT_ZERO"
        else:
            return solution, solution.status
    except Exception as e:
        return None, str(e)

def main():
    base_path = Path(__file__).parent.parent
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    new_model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    media_tsv = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5" / "Acetate_YE0p5__nocmnt__normalized.tsv"
    
    print("="*70)
    print("레퍼런스 모델 vs 신규 모델 직접 비교")
    print("="*70)
    
    # 모델 로드
    print("\n모델 로드 중...")
    ref_model = load_model(ref_model_path)
    new_model = load_model(new_model_path)
    
    # Step 1: 레퍼런스 모델 테스트 (미디어 없이)
    print("\n" + "="*70)
    print("Step 1: 레퍼런스 모델 기본 테스트 (미디어 없이)")
    print("="*70)
    
    ref_test = cobra.io.read_sbml_model(str(ref_model_path))
    ref_test.objective = "Growth"
    
    # 기본 exchange 열기
    for ex in ref_test.exchanges:
        ex.bounds = (-1000, 1000)
    
    ref_sol, ref_status = test_fba_simple(ref_test)
    if ref_sol:
        print(f"레퍼런스 모델 (무제한 미디어): {ref_status}, 성장률 = {ref_sol.objective_value:.6f}")
    
    # Step 2: 레퍼런스 모델 + 미디어
    print("\n" + "="*70)
    print("Step 2: 레퍼런스 모델 + 미디어")
    print("="*70)
    
    ref_test2 = cobra.io.read_sbml_model(str(ref_model_path))
    ref_test2, n_applied = apply_media_from_tsv(ref_test2, media_tsv)
    print(f"미디어 적용: {n_applied}개 exchange")
    
    ref_test2.objective = "Growth"
    ref_sol2, ref_status2 = test_fba_simple(ref_test2)
    if ref_sol2:
        print(f"레퍼런스 모델 (미디어 적용): {ref_status2}, 성장률 = {ref_sol2.objective_value:.6f}")
    
    # Step 3: 신규 모델 + 9개 반응 추가 + 미디어
    print("\n" + "="*70)
    print("Step 3: 신규 모델 + 9개 반응 추가 + 미디어")
    print("="*70)
    
    # 9개 반응 추가
    missing_reactions = [
        'ACS_ADP', 'PEPCK_ATP', 'SUCDi', 'ACtexi',
        'EX_hco3_e', 'T_hco3_e_to_c', 'T_o2_e_to_o2_c',
        'T_nh4_e_to_nh4_c', 'T_fe3_e_to_fe3_c'
    ]
    
    added_count = 0
    for rxn_id in missing_reactions:
        if add_reaction_to_model(new_model, ref_model, rxn_id):
            added_count += 1
    
    print(f"추가된 반응: {added_count}/{len(missing_reactions)}개")
    
    # 미디어 적용
    new_model, n_applied = apply_media_from_tsv(new_model, media_tsv)
    print(f"미디어 적용: {n_applied}개 exchange")
    
    new_model.objective = "Growth"
    new_sol, new_status = test_fba_simple(new_model)
    if new_sol:
        print(f"신규 모델 (9개 반응 + 미디어): {new_status}, 성장률 = {new_sol.objective_value:.6f}")
    
    # Step 4: 신규 모델 + 모든 누락 반응 추가
    print("\n" + "="*70)
    print("Step 4: 신규 모델 + 모든 누락 반응 추가 + 미디어")
    print("="*70)
    
    # 모든 누락 반응 추가
    missing_df = pd.read_csv(base_path / "Stenotrophomonas-causal AI" / "comprehensive_missing_reactions.csv")
    
    all_added = 0
    for _, row in missing_df.iterrows():
        rxn_id = row['reaction_id']
        if add_reaction_to_model(new_model, ref_model, rxn_id):
            all_added += 1
    
    print(f"추가된 반응: {all_added}개 (전체 누락 반응)")
    
    # 미디어 재적용
    new_model2 = cobra.io.read_sbml_model(str(new_model_path))
    for rxn_id in missing_df['reaction_id']:
        add_reaction_to_model(new_model2, ref_model, rxn_id)
    
    new_model2, n_applied = apply_media_from_tsv(new_model2, media_tsv)
    new_model2.objective = "Growth"
    new_sol2, new_status2 = test_fba_simple(new_model2)
    if new_sol2:
        print(f"신규 모델 (모든 반응 + 미디어): {new_status2}, 성장률 = {new_sol2.objective_value:.6f}")
    
    # 최종 비교
    print("\n" + "="*70)
    print("최종 비교")
    print("="*70)
    print(f"{'조건':<40} {'상태':<20} {'성장률':>15}")
    print("-" * 75)
    
    if ref_sol:
        print(f"{'레퍼런스 (무제한 미디어)':<40} {ref_status:<20} {ref_sol.objective_value:>15.6f}")
    if ref_sol2:
        print(f"{'레퍼런스 (미디어 적용)':<40} {ref_status2:<20} {ref_sol2.objective_value:>15.6f}")
    if new_sol:
        print(f"{'신규 (9개 반응 + 미디어)':<40} {new_status:<20} {new_sol.objective_value:>15.6f}")
    if new_sol2:
        print(f"{'신규 (모든 반응 + 미디어)':<40} {new_status2:<20} {new_sol2.objective_value:>15.6f}")
    
    print("\n" + "="*70)

if __name__ == "__main__":
    main()
