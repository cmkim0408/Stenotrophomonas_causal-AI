#!/usr/bin/env python
"""
최소한의 반응 집합 찾기
모든 반응을 추가하면 성공하므로, 최소한 몇 개의 반응이 필요한지 찾기
"""

import cobra
from pathlib import Path
import pandas as pd
import copy

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
    
    id_col = 'exchange_id' if 'exchange_id' in df.columns else df.columns[0]
    
    lb_col = None
    ub_col = None
    for col in df.columns:
        col_lower = str(col).lower()
        if 'min' in col_lower or 'lower' in col_lower:
            lb_col = col
        elif 'max' in col_lower or 'upper' in col_lower:
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
    
    return model

def test_fba(model):
    """FBA 테스트"""
    biomass_rxns = [r for r in model.reactions if 'growth' in r.id.lower() or 'biomass' in r.id.lower()]
    if not biomass_rxns:
        return None, 0.0
    
    model.objective = biomass_rxns[0].id
    
    try:
        solution = model.optimize()
        growth = solution.objective_value if solution.objective_value else 0.0
        if solution.status == 'optimal' and growth > 1e-6:
            return solution, growth
        else:
            return solution, 0.0
    except:
        return None, 0.0

def find_minimal_set(ref_model, new_model, missing_df, media_tsv):
    """최소한의 반응 집합 찾기 (이진 탐색)"""
    print("\n" + "="*70)
    print("최소한의 반응 집합 찾기")
    print("="*70)
    
    # 실제 사용된 반응 먼저 추가
    active_missing = missing_df[missing_df['is_active'] == True].copy()
    print(f"\n1단계: 실제 사용된 반응 {len(active_missing)}개 추가")
    
    base_model = copy.deepcopy(new_model)
    for _, row in active_missing.iterrows():
        add_reaction_to_model(base_model, ref_model, row['reaction_id'])
    
    base_model = apply_media_from_tsv(base_model, media_tsv)
    sol, growth = test_fba(base_model)
    print(f"  성장률: {growth:.6f}")
    
    if growth > 1e-6:
        print("[SUCCESS] 실제 사용된 반응만으로도 성공!")
        return active_missing['reaction_id'].tolist()
    
    # 실제 사용 안 된 반응들을 우선순위별로 추가
    inactive_missing = missing_df[missing_df['is_active'] == False].copy()
    
    # HIGH -> MEDIUM -> LOW 순서
    for priority in ['HIGH', 'MEDIUM', 'LOW']:
        priority_rxns = inactive_missing[inactive_missing['priority'] == priority]
        if len(priority_rxns) == 0:
            continue
        
        print(f"\n{priority} 우선순위 반응 {len(priority_rxns)}개 추가 중...")
        
        for _, row in priority_rxns.iterrows():
            add_reaction_to_model(base_model, ref_model, row['reaction_id'])
        
        test_model = apply_media_from_tsv(copy.deepcopy(base_model), media_tsv)
        sol, growth = test_fba(test_model)
        print(f"  성장률: {growth:.6f}")
        
        if growth > 1e-6:
            print(f"[SUCCESS] {priority} 우선순위까지 추가하면 성공!")
            return base_model.reactions.list_attr('id')
    
    # 모든 반응 추가
    print("\n모든 반응 추가...")
    for _, row in inactive_missing.iterrows():
        add_reaction_to_model(base_model, ref_model, row['reaction_id'])
    
    test_model = apply_media_from_tsv(copy.deepcopy(base_model), media_tsv)
    sol, growth = test_fba(test_model)
    print(f"  성장률: {growth:.6f}")
    
    if growth > 1e-6:
        print("[SUCCESS] 모든 반응 추가로 성공!")
        return base_model.reactions.list_attr('id')
    
    return None

def main():
    base_path = Path(__file__).parent.parent
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    new_model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    media_tsv = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5" / "Acetate_YE0p5__nocmnt__normalized.tsv"
    missing_file = base_path / "Stenotrophomonas-causal AI" / "comprehensive_missing_reactions.csv"
    
    print("="*70)
    print("최소한의 반응 집합 찾기")
    print("="*70)
    
    # 모델 로드
    ref_model = load_model(ref_model_path)
    new_model = load_model(new_model_path)
    missing_df = pd.read_csv(missing_file)
    
    # 최소 집합 찾기
    minimal_set = find_minimal_set(ref_model, new_model, missing_df, media_tsv)
    
    if minimal_set:
        print("\n" + "="*70)
        print("최소 반응 집합")
        print("="*70)
        print(f"총 {len(minimal_set)}개 반응이 필요합니다")
        
        # 추가된 반응만 필터링
        original_rxns = set(new_model.reactions.list_attr('id'))
        added_rxns = [r for r in minimal_set if r not in original_rxns]
        
        print(f"\n추가된 반응: {len(added_rxns)}개")
        
        # 결과 저장
        result_df = missing_df[missing_df['reaction_id'].isin(added_rxns)].copy()
        output_file = base_path / "Stenotrophomonas-causal AI" / "minimal_required_reactions.csv"
        result_df.to_csv(output_file, index=False, encoding='utf-8-sig')
        print(f"\n[OK] 최소 필요 반응 리스트 저장: {output_file}")
        
        # 경로별 요약
        print("\n경로별 최소 필요 반응:")
        pathway_summary = result_df.groupby('pathway').size().sort_values(ascending=False)
        for pathway, count in pathway_summary.items():
            print(f"  {pathway}: {count}개")

if __name__ == "__main__":
    main()
