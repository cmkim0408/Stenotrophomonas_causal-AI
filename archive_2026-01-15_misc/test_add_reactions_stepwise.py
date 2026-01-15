#!/usr/bin/env python
"""
레퍼런스 모델에 있는 누락 반응들을 단계적으로 추가하며 FBA 테스트
우선순위에 따라 반응을 추가하고 각 단계에서 성장 가능 여부 확인
"""

import cobra
from pathlib import Path
import pandas as pd
import copy

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

def setup_acetate_medium(model):
    """Acetate 기반 미디어 설정"""
    exchanges = {
        'EX_ac_e': (-10, 1000),
        'EX_o2_e': (-1000, 1000),
        'EX_h2o_e': (-1000, 1000),
        'EX_h_e': (-1000, 1000),
        'EX_nh4_e': (-10, 1000),
        'EX_pi_e': (-10, 1000),
        'EX_so4_e': (-10, 1000),
        'EX_k_e': (-1000, 1000),
        'EX_na1_e': (-1000, 1000),
        'EX_mg2_e': (-1000, 1000),
        'EX_ca2_e': (-1000, 1000),
        'EX_fe2_e': (-10, 1000),
        'EX_fe3_e': (-10, 1000),
        'EX_hco3_e': (-1000, 1000),
        'EX_co2_e': (-1000, 1000),
    }
    
    for ex_id, bounds in exchanges.items():
        if ex_id in model.reactions:
            model.reactions.get_by_id(ex_id).bounds = bounds
    
    return model

def test_fba(model):
    """FBA 테스트"""
    biomass_rxns = [r for r in model.reactions if 'growth' in r.id.lower() or 'biomass' in r.id.lower()]
    if not biomass_rxns:
        return None, "Biomass reaction not found"
    
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

def add_reactions_by_priority(model, ref_model, missing_df, priority_list):
    """우선순위에 따라 반응 추가"""
    results = []
    
    for priority in priority_list:
        # 해당 우선순위의 반응들 필터링
        priority_rxns = missing_df[missing_df['priority'] == priority]
        
        print(f"\n{priority} 우선순위 반응 추가 중... ({len(priority_rxns)}개)")
        
        added_count = 0
        for _, row in priority_rxns.iterrows():
            rxn_id = row['reaction_id']
            if add_reaction_to_model(model, ref_model, rxn_id):
                added_count += 1
                if added_count % 10 == 0:
                    print(f"  {added_count}/{len(priority_rxns)}개 추가 완료...")
        
        print(f"  {added_count}개 반응 추가 완료")
        
        # FBA 테스트
        test_model = setup_acetate_medium(copy.deepcopy(model))
        solution, status = test_fba(test_model)
        
        growth = solution.objective_value if solution and solution.objective_value else 0.0
        
        results.append({
            'priority': priority,
            'added_reactions': added_count,
            'total_reactions': len(priority_rxns),
            'fba_status': status,
            'growth_rate': growth
        })
        
        print(f"  FBA 상태: {status}, 성장률: {growth:.6f}")
        
        # 성공하면 중단
        if status == "SUCCESS":
            print(f"\n[SUCCESS] {priority} 우선순위까지 추가하면 FBA가 성공합니다!")
            break
    
    return results

def main():
    base_path = Path(__file__).parent.parent
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    new_model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    missing_file = base_path / "Stenotrophomonas-causal AI" / "comprehensive_missing_reactions.csv"
    
    print("="*70)
    print("누락된 반응들을 우선순위별로 추가하며 FBA 테스트")
    print("="*70)
    
    # 모델 로드
    ref_model = load_model(ref_model_path)
    new_model = load_model(new_model_path)
    
    # 누락된 반응 리스트 로드
    missing_df = pd.read_csv(missing_file)
    print(f"\n누락된 반응 총 {len(missing_df)}개")
    
    # 우선순위별 개수
    print("\n우선순위별 분류:")
    priority_counts = missing_df.groupby('priority').size()
    for priority, count in priority_counts.items():
        print(f"  {priority}: {count}개")
    
    # 실제 사용된 반응 먼저 추가 (HIGH)
    print("\n" + "="*70)
    print("Step 1: 실제 사용된 반응 추가 (HIGH - 실제 사용됨)")
    print("="*70)
    
    active_missing = missing_df[missing_df['is_active'] == True].copy()
    active_missing['priority'] = 'HIGH_ACTIVE'  # 실제 사용된 것 먼저
    
    # 실제 사용된 반응 추가
    print(f"\n실제 사용된 반응 {len(active_missing)}개 추가 중...")
    for _, row in active_missing.iterrows():
        rxn_id = row['reaction_id']
        add_reaction_to_model(new_model, ref_model, rxn_id)
    
    # 실제 사용된 반응 추가 후 테스트
    test_model = setup_acetate_medium(copy.deepcopy(new_model))
    solution, status = test_fba(test_model)
    growth = solution.objective_value if solution and solution.objective_value else 0.0
    print(f"FBA 상태: {status}, 성장률: {growth:.6f}")
    
    if status == "SUCCESS":
        print("\n[SUCCESS] 실제 사용된 반응만으로도 FBA가 성공합니다!")
        return
    
    # 우선순위별로 추가
    print("\n" + "="*70)
    print("Step 2: 우선순위별로 반응 추가")
    print("="*70)
    
    # 실제 사용 안 된 반응들만
    inactive_missing = missing_df[missing_df['is_active'] == False].copy()
    
    # HIGH -> MEDIUM -> LOW 순서로 추가
    results = add_reactions_by_priority(new_model, ref_model, inactive_missing, ['HIGH', 'MEDIUM', 'LOW'])
    
    # 최종 요약
    print("\n" + "="*70)
    print("최종 결과 요약")
    print("="*70)
    
    print("\n단계별 추가 결과:")
    print(f"{'단계':<20} {'추가된 반응':>15} {'FBA 상태':<20} {'성장률':>15}")
    print("-" * 70)
    
    print(f"{'실제 사용된 반응':<20} {len(active_missing):>15} {status:<20} {growth:>15.6f}")
    
    for result in results:
        print(f"{result['priority']:<20} {result['added_reactions']:>15} {result['fba_status']:<20} {result['growth_rate']:>15.6f}")
    
    # 최종 상태
    final_model = setup_acetate_medium(copy.deepcopy(new_model))
    final_solution, final_status = test_fba(final_model)
    final_growth = final_solution.objective_value if final_solution and final_solution.objective_value else 0.0
    
    print("\n" + "="*70)
    print("최종 상태")
    print("="*70)
    print(f"추가된 총 반응 수: {len(missing_df)}개")
    print(f"최종 FBA 상태: {final_status}")
    print(f"최종 성장률: {final_growth:.6f}")
    
    if final_status == "SUCCESS":
        print("\n[SUCCESS] 모든 반응을 추가하면 FBA가 성공합니다!")
    else:
        print("\n[FAIL] 모든 반응을 추가해도 FBA가 실패합니다.")
        print("  -> 추가 진단 필요")

if __name__ == "__main__":
    main()
