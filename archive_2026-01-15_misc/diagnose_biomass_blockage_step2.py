#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 4-2: "Biomass가 요구하는 성분 중 뭐가 안 만들어지는지" 찾기

가장 빠른 방법:
- Biomass 반응이 소비(–)하는 각 metabolite m_i에 대해 DM_mi (demand)를 만들고
- max DM_mi가 0이면 그 전구체/보조인자 생산 경로가 끊긴 것

이 작업을 하면 "114개 중 뭘 넣어야 하나?"가 바로 핵심 후보 몇 개로 줄어듭니다.
"""

import cobra
from pathlib import Path
import pandas as pd

def find_biomass_reaction(model):
    """Biomass 반응 찾기"""
    biomass_keywords = ['biomass', 'growth', 'BIOMASS', 'Growth']
    for rxn in model.reactions:
        if any(keyword in rxn.id for keyword in biomass_keywords):
            return rxn
    return None

def setup_media_forced(model):
    """배지 조건을 강제로 고정 (Step 4-1과 동일)"""
    # Acetate uptake 고정
    if 'EX_ac_e' in model.reactions:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.lower_bound = -19.0
        ex_ac.upper_bound = -19.0
    
    # O₂ uptake 제한/고정
    if 'EX_o2_e' in model.reactions:
        ex_o2 = model.reactions.get_by_id('EX_o2_e')
        ex_o2.lower_bound = -100.0
        ex_o2.upper_bound = 1000.0
    
    # 필수 무기물 열어두기
    essential_exchanges = {
        'EX_nh4_e': (-1000.0, 1000.0),
        'EX_pi_e': (-1000.0, 1000.0),
        'EX_so4_e': (-1000.0, 1000.0),
        'EX_mg2_e': (-1000.0, 1000.0),
        'EX_k_e': (-1000.0, 1000.0),
        'EX_na1_e': (-1000.0, 1000.0),
        'EX_fe2_e': (-1000.0, 1000.0),
        'EX_fe3_e': (-1000.0, 1000.0),
        'EX_h2o_e': (-1000.0, 1000.0),
        'EX_h_e': (-1000.0, 1000.0),
        'EX_co2_e': (-1000.0, 1000.0),
        'EX_hco3_e': (-1000.0, 1000.0),
    }
    
    for ex_id, (lb, ub) in essential_exchanges.items():
        if ex_id in model.reactions:
            ex_rxn = model.reactions.get_by_id(ex_id)
            ex_rxn.lower_bound = lb
            ex_rxn.upper_bound = ub
    
    return model

def test_metabolite_production(model, metabolite_id):
    """특정 대사물질의 최대 생산량 테스트"""
    if metabolite_id not in model.metabolites:
        return None, "metabolite_not_in_model"
    
    with model:
        # Demand 반응 추가
        demand_id = f'DM_{metabolite_id}'
        
        # 기존 demand 반응이 있으면 제거
        if demand_id in model.reactions:
            model.remove_reactions([demand_id])
        
        demand_rxn = cobra.Reaction(demand_id)
        demand_rxn.add_metabolites({model.metabolites.get_by_id(metabolite_id): -1})
        model.add_reactions([demand_rxn])
        
        # Objective를 demand로 설정
        model.objective = demand_id
        model.objective_direction = 'max'
        
        # FBA 실행
        solution = model.optimize()
        
        if solution.status == 'optimal':
            max_prod = solution.objective_value
            return max_prod, "optimal"
        else:
            return 0.0, solution.status

def test_biomass_components(model):
    """Biomass 반응이 소비하는 각 metabolite의 생산 가능성 테스트"""
    print("\n" + "="*70)
    print("Step 4-2: Biomass 성분 생산 가능성 테스트")
    print("="*70)
    
    # Biomass 반응 찾기
    biomass_rxn = find_biomass_reaction(model)
    if biomass_rxn is None:
        print("[ERROR] Biomass 반응을 찾을 수 없습니다.")
        return None
    
    print(f"\n[Biomass 반응] {biomass_rxn.id}")
    
    # Biomass 반응이 소비하는 metabolite 추출
    consumed_metabolites = []
    for met, coeff in biomass_rxn.metabolites.items():
        if coeff < 0:  # 소비되는 metabolite (음수 계수)
            consumed_metabolites.append((met.id, abs(coeff), met.name if hasattr(met, 'name') else met.id))
    
    print(f"\n[Biomass가 소비하는 metabolite 수] {len(consumed_metabolites)}개")
    
    # 각 metabolite의 생산 가능성 테스트
    results = []
    blocked_count = 0
    
    print(f"\n[테스트 진행 중...] (총 {len(consumed_metabolites)}개)")
    
    for i, (met_id, coeff, met_name) in enumerate(consumed_metabolites):
        if (i + 1) % 10 == 0:
            print(f"  진행: {i+1}/{len(consumed_metabolites)}")
        
        max_prod, status = test_metabolite_production(model, met_id)
        
        if max_prod is None:
            results.append({
                'metabolite_id': met_id,
                'metabolite_name': met_name,
                'biomass_coeff': coeff,
                'max_production': None,
                'status': status,
                'is_blocked': True
            })
            blocked_count += 1
        elif max_prod < 1e-6:
            results.append({
                'metabolite_id': met_id,
                'metabolite_name': met_name,
                'biomass_coeff': coeff,
                'max_production': max_prod,
                'status': status,
                'is_blocked': True
            })
            blocked_count += 1
        else:
            results.append({
                'metabolite_id': met_id,
                'metabolite_name': met_name,
                'biomass_coeff': coeff,
                'max_production': max_prod,
                'status': status,
                'is_blocked': False
            })
    
    # 결과 정리
    df = pd.DataFrame(results)
    blocked_df = df[df['is_blocked'] == True].copy()
    
    print(f"\n[결과 요약]")
    print(f"  총 테스트: {len(consumed_metabolites)}개")
    print(f"  막힌 metabolite: {blocked_count}개")
    print(f"  생성 가능: {len(consumed_metabolites) - blocked_count}개")
    
    if len(blocked_df) > 0:
        print(f"\n[막힌 metabolite 목록]")
        print(f"{'Metabolite ID':<30} {'Name':<40} {'Biomass Coeff':>15} {'Max Production':>15}")
        print("-" * 100)
        
        # 계수 순으로 정렬 (큰 것부터)
        blocked_df_sorted = blocked_df.sort_values('biomass_coeff', ascending=False)
        
        for idx, row in blocked_df_sorted.iterrows():
            met_id = row['metabolite_id']
            met_name = row['metabolite_name'] if pd.notna(row['metabolite_name']) else ''
            coeff = row['biomass_coeff']
            max_prod = row['max_production'] if pd.notna(row['max_production']) else 0.0
            
            # 이름이 너무 길면 자르기
            if len(met_name) > 37:
                met_name = met_name[:34] + "..."
            
            max_prod_str = f"{max_prod:.6f}" if max_prod is not None else "N/A"
            print(f"{met_id:<30} {met_name:<40} {coeff:>15.6f} {max_prod_str:>15}")
        
        # 계수 합계로 우선순위 추정
        print(f"\n[계수 합계]")
        total_blocked_coeff = blocked_df['biomass_coeff'].sum()
        print(f"  막힌 metabolite의 총 계수 합: {total_blocked_coeff:.6f}")
        print(f"  (계수가 클수록 biomass에 더 많이 필요)")
    else:
        print("\n[OK] 모든 metabolite가 생성 가능합니다!")
    
    return df, blocked_df

def categorize_blocked_metabolites(blocked_df):
    """막힌 metabolite를 카테고리별로 분류"""
    print("\n" + "="*70)
    print("막힌 metabolite 카테고리별 분류")
    print("="*70)
    
    categories = {
        'Amino acids': [],
        'Nucleotides': [],
        'Cofactors': [],
        'Ions': [],
        'Others': []
    }
    
    aa_keywords = ['__L_c', '__L_e', 'ala', 'arg', 'asn', 'asp', 'cys', 'gln', 'glu', 'gly', 
                   'his', 'ile', 'leu', 'lys', 'met', 'phe', 'pro', 'ser', 'thr', 'trp', 'tyr', 'val']
    nuc_keywords = ['atp', 'ctp', 'gtp', 'utp', 'dttp', 'datp', 'dctp', 'dgtp']
    cofactor_keywords = ['nad', 'nadp', 'fad', 'coa', 'thf', 'pydx', 'ribflv', 'thm']
    ion_keywords = ['ca2', 'cl', 'fe2', 'fe3', 'k', 'mg2', 'mn2', 'na1', 'so4', 'zn2', 'cobalt2', 'cu2', 'ni2']
    
    for idx, row in blocked_df.iterrows():
        met_id = row['metabolite_id'].lower()
        classified = False
        
        if any(keyword in met_id for keyword in aa_keywords):
            categories['Amino acids'].append((row['metabolite_id'], row['biomass_coeff']))
            classified = True
        elif any(keyword in met_id for keyword in nuc_keywords):
            categories['Nucleotides'].append((row['metabolite_id'], row['biomass_coeff']))
            classified = True
        elif any(keyword in met_id for keyword in cofactor_keywords):
            categories['Cofactors'].append((row['metabolite_id'], row['biomass_coeff']))
            classified = True
        elif any(keyword in met_id for keyword in ion_keywords):
            categories['Ions'].append((row['metabolite_id'], row['biomass_coeff']))
            classified = True
        
        if not classified:
            categories['Others'].append((row['metabolite_id'], row['biomass_coeff']))
    
    for cat_name, items in categories.items():
        if len(items) > 0:
            print(f"\n[{cat_name}] ({len(items)}개)")
            total_coeff = sum(coeff for _, coeff in items)
            print(f"  총 계수 합: {total_coeff:.6f}")
            # 계수 순으로 정렬
            items_sorted = sorted(items, key=lambda x: x[1], reverse=True)
            for met_id, coeff in items_sorted[:10]:  # 상위 10개만 표시
                print(f"    {met_id:<30} (coeff: {coeff:.6f})")
            if len(items) > 10:
                print(f"    ... 외 {len(items)-10}개")
    
    return categories

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_with_BCAA_cofactors_ions_nad_transport.xml"
    output_dir = base_path / "Stenotrophomonas-causal AI"
    
    print("="*70)
    print("Step 4-2: Biomass 성분 생산 가능성 테스트")
    print("="*70)
    print(f"\n모델: {model_path}")
    
    # 모델 로드
    if not model_path.exists():
        print(f"[ERROR] 모델 파일이 없습니다: {model_path}")
        return
    
    model = cobra.io.read_sbml_model(str(model_path))
    print(f"[OK] 모델 로드 완료 (반응 수: {len(model.reactions)})")
    
    # 배지 조건 강제 고정
    model = setup_media_forced(model)
    
    # Biomass 성분 생산 가능성 테스트
    df, blocked_df = test_biomass_components(model)
    
    if blocked_df is not None and len(blocked_df) > 0:
        # 카테고리별 분류
        categories = categorize_blocked_metabolites(blocked_df)
        
        # 결과 저장
        output_file = output_dir / "biomass_blocked_metabolites.csv"
        blocked_df.to_csv(output_file, index=False, encoding='utf-8-sig')
        print(f"\n[결과 저장] {output_file}")
        
        print("\n" + "="*70)
        print("다음 단계: Step 4-3 - 막힌 전구체 주변을 우선 복구")
        print("="*70)
    else:
        print("\n[OK] 모든 metabolite가 생성 가능합니다!")
    
    return df, blocked_df

if __name__ == "__main__":
    df, blocked_df = main()
