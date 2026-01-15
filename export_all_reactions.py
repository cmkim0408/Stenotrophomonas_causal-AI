#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
모든 반응의 Flux Export

최종 모델의 모든 반응의 flux를 CSV로 export
"""

import cobra
from pathlib import Path
import pandas as pd
import sys

def load_model(model_path):
    try:
        model = cobra.io.read_sbml_model(str(model_path))
        print(f"[OK] 모델 로드 완료 (반응 수: {len(model.reactions)})")
        return model
    except Exception as e:
        print(f"[ERROR] 모델 로드 실패: {e}")
        sys.exit(1)

def setup_media_forced(model):
    """배지 조건을 강제로 고정"""
    if 'EX_ac_e' in model.reactions:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.lower_bound = -19.0
        ex_ac.upper_bound = -19.0
    
    if 'EX_o2_e' in model.reactions:
        ex_o2 = model.reactions.get_by_id('EX_o2_e')
        ex_o2.lower_bound = -100.0
        ex_o2.upper_bound = 1000.0
    
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
        'EX_nac_e': (-1000.0, 1000.0),
        'EX_ncam_e': (-1000.0, 1000.0),
    }
    
    for ex_id, (lb, ub) in essential_exchanges.items():
        if ex_id in model.reactions:
            ex_rxn = model.reactions.get_by_id(ex_id)
            ex_rxn.lower_bound = lb
            ex_rxn.upper_bound = ub
    
    return model

def export_all_reactions(model, solution, output_file):
    """모든 반응의 flux를 CSV로 export"""
    print("\n" + "="*70)
    print("모든 반응 Flux Export")
    print("="*70)
    
    all_reactions = []
    
    for rxn in model.reactions:
        flux = solution.fluxes.get(rxn.id, 0.0)
        
        all_reactions.append({
            'reaction_id': rxn.id,
            'reaction_name': rxn.name if rxn.name else '',
            'reaction_equation': rxn.reaction,
            'flux': flux,
            'lower_bound': rxn.lower_bound,
            'upper_bound': rxn.upper_bound,
            'reversible': rxn.reversibility,
        })
    
    # DataFrame 생성
    df = pd.DataFrame(all_reactions)
    df = df.sort_values('reaction_id')
    
    # CSV 저장
    output_path = Path(output_file)
    df.to_csv(output_path, index=False, encoding='utf-8-sig')
    print(f"\n[모든 반응 Flux 저장] {output_path}")
    print(f"  총 {len(df)}개 반응")
    
    # 통계
    nonzero_count = len(df[df['flux'].abs() > 1e-6])
    zero_count = len(df) - nonzero_count
    
    print(f"\n[통계]")
    print(f"  Non-zero flux: {nonzero_count}개")
    print(f"  Zero flux: {zero_count}개")
    print(f"  총: {len(df)}개")
    
    return df

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_final_cleaned.xml"
    output_dir = base_path / "Stenotrophomonas-causal AI"
    
    print("="*70)
    print("모든 반응 Flux Export")
    print("="*70)
    
    # 모델 로드
    model = load_model(model_path)
    model = setup_media_forced(model)
    
    # ATPM=0 설정
    if 'ATPM' in model.reactions:
        atpm = model.reactions.get_by_id('ATPM')
        atpm.lower_bound = 0.0
        atpm.upper_bound = 0.0
        print(f"\n[ATPM 설정] bounds: [{atpm.lower_bound}, {atpm.upper_bound}]")
    
    # FBA 실행
    solution = model.optimize()
    print(f"\n[FBA 결과]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    
    # 모든 반응 export
    output_file = output_dir / "all_reactions_flux.csv"
    df = export_all_reactions(model, solution, output_file)
    
    print("\n" + "="*70)
    print("완료")
    print("="*70)
    
    return model, solution, df

if __name__ == "__main__":
    model, solution, df = main()
