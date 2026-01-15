#!/usr/bin/env python
"""
레퍼런스 모델 vs 신규 모델 구조적 차이 찾기
- 핵심 반응의 존재 여부 확인
- 반응식 차이 확인
- bounds 차이 확인
"""

import cobra
from pathlib import Path

def load_model(model_path):
    model = cobra.io.read_sbml_model(model_path)
    return model

def compare_key_reactions(ref_model, new_model):
    """핵심 반응 비교"""
    print("="*80)
    print("핵심 반응 비교")
    print("="*80)
    
    key_reactions = [
        'ACS_ADP', 'ACS', 'CS', 'ICL', 'MALS', 'SUCDi', 
        'ATPS4rpp', 'NADH16pp', 'CYTBO3_4pp',
        'PEPCK_ATP', 'ACtexi', 'ADK1'
    ]
    
    print(f"\n{'반응':<20} {'레퍼런스':<15} {'신규':<15} {'차이':<20}")
    print("-" * 70)
    
    differences = []
    
    for rxn_id in key_reactions:
        ref_exists = rxn_id in ref_model.reactions
        new_exists = rxn_id in new_model.reactions
        
        ref_info = ""
        new_info = ""
        diff_info = ""
        
        if ref_exists and new_exists:
            ref_rxn = ref_model.reactions.get_by_id(rxn_id)
            new_rxn = new_model.reactions.get_by_id(rxn_id)
            
            ref_bounds = f"[{ref_rxn.lower_bound}, {ref_rxn.upper_bound}]"
            new_bounds = f"[{new_rxn.lower_bound}, {new_rxn.upper_bound}]"
            
            ref_info = f"존재 {ref_bounds}"
            new_info = f"존재 {new_bounds}"
            
            # 반응식 비교
            if ref_rxn.reaction != new_rxn.reaction:
                diff_info = "반응식 다름"
                differences.append(f"{rxn_id}: 반응식 다름")
            elif ref_bounds != new_bounds:
                diff_info = "bounds 다름"
                differences.append(f"{rxn_id}: bounds 다름")
            else:
                diff_info = "동일"
        elif ref_exists and not new_exists:
            ref_info = "존재"
            new_info = "없음"
            diff_info = "신규 모델에 없음"
            differences.append(f"{rxn_id}: 신규 모델에 없음")
        elif not ref_exists and new_exists:
            ref_info = "없음"
            new_info = "존재"
            diff_info = "레퍼런스에 없음"
        else:
            ref_info = "없음"
            new_info = "없음"
            diff_info = "둘 다 없음"
        
        if diff_info != "동일" and diff_info != "":
            print(f"{rxn_id:<20} {ref_info:<15} {new_info:<15} {diff_info:<20}")
    
    if differences:
        print(f"\n[차이점 요약]")
        for diff in differences:
            print(f"  - {diff}")
    else:
        print(f"\n[차이점 없음]")

def compare_exchange_bounds(ref_model, new_model):
    """Exchange bounds 비교"""
    print("\n" + "="*80)
    print("Exchange bounds 비교")
    print("="*80)
    
    key_exchanges = ['EX_ac_e', 'EX_o2_e', 'EX_hco3_e', 'EX_nh4_e', 'EX_pi_e', 'EX_h2o_e', 'EX_h_e']
    
    print(f"\n{'Exchange':<20} {'레퍼런스 bounds':<25} {'신규 bounds':<25} {'차이':<15}")
    print("-" * 85)
    
    differences = []
    
    for ex_id in key_exchanges:
        ref_exists = ex_id in ref_model.reactions
        new_exists = ex_id in new_model.reactions
        
        if ref_exists and new_exists:
            ref_rxn = ref_model.reactions.get_by_id(ex_id)
            new_rxn = new_model.reactions.get_by_id(ex_id)
            
            ref_bounds = f"[{ref_rxn.lower_bound}, {ref_rxn.upper_bound}]"
            new_bounds = f"[{new_rxn.lower_bound}, {new_rxn.upper_bound}]"
            
            if ref_bounds != new_bounds:
                print(f"{ex_id:<20} {ref_bounds:<25} {new_bounds:<25} {'다름':<15}")
                differences.append(f"{ex_id}: bounds 다름")
        elif ref_exists and not new_exists:
            print(f"{ex_id:<20} {'존재':<25} {'없음':<25} {'신규에 없음':<15}")
            differences.append(f"{ex_id}: 신규 모델에 없음")
        elif not ref_exists and new_exists:
            print(f"{ex_id:<20} {'없음':<25} {'존재':<25} {'레퍼런스에 없음':<15}")
    
    if differences:
        print(f"\n[차이점 요약]")
        for diff in differences:
            print(f"  - {diff}")

def compare_atpm_setting(ref_model, new_model):
    """ATPM 설정 비교"""
    print("\n" + "="*80)
    print("ATPM 설정 비교")
    print("="*80)
    
    ref_atpm = ref_model.reactions.get_by_id('ATPM')
    new_atpm = new_model.reactions.get_by_id('ATPM')
    
    print(f"\n[ATPM 반응]")
    print(f"  레퍼런스 bounds: [{ref_atpm.lower_bound}, {ref_atpm.upper_bound}]")
    print(f"  신규 bounds: [{new_atpm.lower_bound}, {new_atpm.upper_bound}]")
    print(f"  레퍼런스 반응식: {ref_atpm.reaction}")
    print(f"  신규 반응식: {new_atpm.reaction}")
    
    if ref_atpm.reaction != new_atpm.reaction:
        print(f"  [차이] 반응식이 다름!")
    if ref_atpm.lower_bound != new_atpm.lower_bound:
        print(f"  [차이] lower_bound가 다름! ({ref_atpm.lower_bound} vs {new_atpm.lower_bound})")

def main():
    base_path = Path(__file__).parent.parent
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    new_model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    
    print("="*80)
    print("레퍼런스 모델 vs 신규 모델 구조적 차이 찾기")
    print("="*80)
    
    # 모델 로드
    ref_model = load_model(str(ref_model_path))
    new_model = load_model(str(new_model_path))
    
    print(f"\n[모델 정보]")
    print(f"  레퍼런스: {len(ref_model.reactions)}개 반응, {len(ref_model.metabolites)}개 메타볼라이트")
    print(f"  신규: {len(new_model.reactions)}개 반응, {len(new_model.metabolites)}개 메타볼라이트")
    
    # 핵심 반응 비교
    compare_key_reactions(ref_model, new_model)
    
    # Exchange bounds 비교
    compare_exchange_bounds(ref_model, new_model)
    
    # ATPM 설정 비교
    compare_atpm_setting(ref_model, new_model)
    
    print("\n" + "="*80)
    print("결론")
    print("="*80)
    print(f"\n[핵심 질문]")
    print(f"  레퍼런스 모델에서는 FBA가 잘 돌아가는데, 신규 모델에서는 왜 안 될까?")
    print(f"  -> 구조적 차이 확인 완료")

if __name__ == "__main__":
    main()
