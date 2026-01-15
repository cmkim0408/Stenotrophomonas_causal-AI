#!/usr/bin/env python
"""
ATPM=0일 때 성장 불가 문제 해결 시도
- 레퍼런스 모델의 미디어 설정 사용
- 누락된 반응 확인
"""

import cobra
from pathlib import Path
import pandas as pd

def load_model(model_path):
    model = cobra.io.read_sbml_model(model_path)
    return model

def apply_media_from_tsv(model, tsv_path):
    """TSV 파일에서 미디어 적용"""
    df = pd.read_csv(tsv_path, sep='\t', comment='#')
    
    # 헤더 확인
    if 'rxn_id' in df.columns:
        rxn_col = 'rxn_id'
    elif 'exchange_id' in df.columns:
        rxn_col = 'exchange_id'
    else:
        rxn_col = df.columns[0]
    
    # bounds 컬럼 찾기
    if 'lower' in df.columns and 'upper' in df.columns:
        lb_col = 'lower'
        ub_col = 'upper'
    elif 'lower_bound' in df.columns and 'upper_bound' in df.columns:
        lb_col = 'lower_bound'
        ub_col = 'upper_bound'
    else:
        # 첫 번째 데이터 행이 헤더인 경우
        df.columns = [rxn_col, 'lower', 'upper']
        lb_col = 'lower'
        ub_col = 'upper'
    
    applied = 0
    for _, row in df.iterrows():
        rxn_id = str(row[rxn_col]).strip()
        try:
            lb = float(row[lb_col])
            ub = float(row[ub_col])
            if rxn_id in model.reactions:
                rxn = model.reactions.get_by_id(rxn_id)
                rxn.lower_bound = lb
                rxn.upper_bound = ub
                applied += 1
        except (KeyError, ValueError, TypeError):
            continue
    
    return applied

def test_with_reference_media(model, media_path):
    """레퍼런스 모델의 미디어로 테스트"""
    print("="*80)
    print("레퍼런스 모델 미디어로 테스트")
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
        print(f"  ICL: {solution.fluxes.get('ICL', 0):.6f}")
        print(f"  MALS: {solution.fluxes.get('MALS', 0):.6f}")
        print(f"  ICDHx: {solution.fluxes.get('ICDHx', 0):.6f}")
        print(f"  CS: {solution.fluxes.get('CS', 0):.6f}")
        
        if solution.objective_value > 1e-6:
            print("\n[OK] 레퍼런스 미디어로 성장 가능!")
        else:
            print("\n[문제] 레퍼런스 미디어로도 성장 불가")
    
    return solution

def compare_media_settings(new_model, ref_model, media_path):
    """미디어 설정 비교"""
    print("\n" + "="*80)
    print("미디어 설정 비교")
    print("="*80)
    
    # 레퍼런스 모델에 미디어 적용
    for rxn in ref_model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    try:
        apply_media_from_tsv(ref_model, media_path)
        ref_atpm = ref_model.reactions.get_by_id('ATPM')
        ref_atpm.lower_bound = 0
        ref_atpm.upper_bound = 1000
        
        ref_model.objective = 'Growth'
        ref_sol = ref_model.optimize()
        
        print(f"\n[레퍼런스 모델] 같은 미디어 사용")
        print(f"  성장률: {ref_sol.objective_value:.6f}")
        print(f"  ICL: {ref_sol.fluxes.get('ICL', 0):.6f}")
        print(f"  MALS: {ref_sol.fluxes.get('MALS', 0):.6f}")
        
        # 신규 모델 테스트
        new_sol = test_with_reference_media(new_model, media_path)
        
        if new_sol and new_sol.objective_value > 1e-6:
            print("\n[결론] 레퍼런스 미디어를 사용하면 성장 가능!")
            return True
        else:
            print("\n[결론] 미디어 차이가 원인일 수 있음")
            return False
            
    except Exception as e:
        print(f"\n[비교 실패] {e}")
        return False

def main():
    base_path = Path(__file__).parent.parent
    
    # 모델 로드
    new_model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    media_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5" / "Acetate_YE0p5__nocmnt__normalized.tsv"
    
    print("="*80)
    print("ATPM=0 성장 불가 문제 해결 시도")
    print("="*80)
    
    new_model = load_model(str(new_model_path))
    print(f"\n[신규 모델 로드 완료]")
    
    if ref_model_path.exists() and media_path.exists():
        ref_model = load_model(str(ref_model_path))
        print(f"[레퍼런스 모델 로드 완료]")
        
        # 미디어 비교
        compare_media_settings(new_model, ref_model, media_path)
    else:
        print(f"\n[레퍼런스 파일 없음]")
        # 레퍼런스 미디어만 테스트
        if media_path.exists():
            test_with_reference_media(new_model, media_path)

if __name__ == "__main__":
    main()
