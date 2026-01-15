#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
레퍼런스 모델에서 BCAA 관련 반응을 신규 모델에 추가

추가할 반응 (6개):
- KARI (Ketol-acid reductoisomerase)
- DHAD (Dihydroxyacid dehydratase)
- IPMI (Isopropylmalate isomerase)
- IPMDH (Isopropylmalate dehydrogenase)
- BCAT_VAL (Branched-chain aminotransferase, Val)
- BCAT_LEU (Branched-chain aminotransferase, Leu)
"""

import cobra
from pathlib import Path
import sys

def get_metabolite(model, met_id):
    """대사물질 가져오기 (없으면 생성)"""
    if met_id in model.metabolites:
        return model.metabolites.get_by_id(met_id)
    else:
        # 레퍼런스 모델에서 복사 필요 (일단 기본 정보만)
        new_met = cobra.Metabolite(id=met_id, compartment='c')
        model.add_metabolites([new_met])
        return new_met

def add_reaction_from_reference(model, ref_model, rxn_id):
    """레퍼런스 모델에서 반응을 신규 모델에 추가"""
    if rxn_id in model.reactions:
        print(f"[SKIP] {rxn_id}: 이미 모델에 있음")
        return False
    
    if rxn_id not in ref_model.reactions:
        print(f"[ERROR] {rxn_id}: 레퍼런스 모델에 없음")
        return False
    
    try:
        ref_rxn = ref_model.reactions.get_by_id(rxn_id)
        
        # 새 반응 생성
        new_rxn = cobra.Reaction(rxn_id)
        new_rxn.name = ref_rxn.name if hasattr(ref_rxn, 'name') and ref_rxn.name else rxn_id
        new_rxn.lower_bound = ref_rxn.lower_bound
        new_rxn.upper_bound = ref_rxn.upper_bound
        
        # 대사물질 추가
        metabolites_dict = {}
        for met, coeff in ref_rxn.metabolites.items():
            met_id = met.id
            
            # 대사물질이 모델에 있는지 확인
            if met_id in model.metabolites:
                metabolites_dict[model.metabolites.get_by_id(met_id)] = coeff
            else:
                # 대사물질이 없으면 레퍼런스에서 복사
                ref_met = ref_model.metabolites.get_by_id(met_id)
                new_met = cobra.Metabolite(
                    id=ref_met.id,
                    formula=getattr(ref_met, 'formula', None),
                    name=getattr(ref_met, 'name', met_id),
                    compartment=ref_met.compartment
                )
                model.add_metabolites([new_met])
                metabolites_dict[new_met] = coeff
        
        new_rxn.add_metabolites(metabolites_dict)
        model.add_reactions([new_rxn])
        
        print(f"[OK] {rxn_id} 추가: {new_rxn.reaction}")
        print(f"     bounds: [{new_rxn.lower_bound}, {new_rxn.upper_bound}]")
        return True
    except Exception as e:
        print(f"[ERROR] {rxn_id} 추가 실패: {e}")
        import traceback
        traceback.print_exc()
        return False

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
    }
    
    for ex_id, (lb, ub) in essential_exchanges.items():
        if ex_id in model.reactions:
            ex_rxn = model.reactions.get_by_id(ex_id)
            ex_rxn.lower_bound = lb
            ex_rxn.upper_bound = ub
    
    return model

def test_biomass_growth(model):
    """Biomass 성장 테스트"""
    biomass_keywords = ['biomass', 'growth', 'BIOMASS', 'Growth']
    for rxn in model.reactions:
        if any(keyword in rxn.id for keyword in biomass_keywords):
            biomass_rxn = rxn
            break
    else:
        return None
    
    model.objective = biomass_rxn.id
    model.objective_direction = 'max'
    solution = model.optimize()
    
    return solution.status, solution.objective_value

def main():
    base_path = Path(__file__).parent.parent
    new_model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_with_ACtexi.xml"
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    output_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_with_BCAA.xml"
    
    print("="*70)
    print("BCAA 반응 추가: 레퍼런스 모델 -> 신규 모델")
    print("="*70)
    
    # 모델 로드
    print(f"\n[모델 로드]")
    new_model = cobra.io.read_sbml_model(str(new_model_path))
    ref_model = cobra.io.read_sbml_model(str(ref_model_path))
    print(f"[OK] 모델 로드 완료")
    
    # 추가할 반응 목록
    reactions_to_add = ['KARI', 'DHAD', 'IPMI', 'IPMDH', 'BCAT_VAL', 'BCAT_LEU']
    
    print(f"\n[추가할 반응] {len(reactions_to_add)}개")
    print(f"  {', '.join(reactions_to_add)}")
    
    # 배지 조건 설정
    new_model = setup_media_forced(new_model)
    
    # 반응 추가
    print(f"\n[반응 추가 중...]")
    added_count = 0
    for rxn_id in reactions_to_add:
        if add_reaction_from_reference(new_model, ref_model, rxn_id):
            added_count += 1
    
    print(f"\n[결과]")
    print(f"  추가된 반응: {added_count}/{len(reactions_to_add)}개")
    
    # 모델 저장
    cobra.io.write_sbml_model(new_model, str(output_path))
    print(f"\n[모델 저장] {output_path}")
    
    # 간단한 성장 테스트
    print(f"\n[성장 테스트]")
    status, growth = test_biomass_growth(new_model)
    print(f"  상태: {status}")
    print(f"  성장률: {growth:.6f}")
    
    if growth > 1e-6:
        print(f"\n[OK] 성장 가능! (성장률: {growth:.6f})")
    else:
        print(f"\n[NO] 아직 성장 불가 (성장률: {growth:.6f})")
        print("  -> 추가 반응이 더 필요할 수 있습니다")
    
    return new_model

if __name__ == "__main__":
    model = main()
