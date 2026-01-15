#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NNATr 반응 추가

NNATr (Nicotinate-nucleotide adenylyltransferase):
- 반응식: atp_c + h_c + nicrnt_c <=> dnad_c + ppi_c
- 이것이 dnad_c를 생성하는 핵심 반응!
- dnad_c → NADS1/NADS2 → nad_c → NADK → nadp_c
"""

import cobra
from pathlib import Path
import sys

def add_reaction_from_reference(model, ref_model, rxn_id):
    """레퍼런스 모델에서 반응을 신규 모델에 추가"""
    if rxn_id in model.reactions:
        return False, "already_exists"
    
    if rxn_id not in ref_model.reactions:
        return False, "not_in_reference"
    
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
        
        return True, "added"
    except Exception as e:
        return False, str(e)

def test_metabolite_production(model, metabolite_id):
    """특정 대사물질의 최대 생산량 테스트"""
    if metabolite_id not in model.metabolites:
        return None, "metabolite_not_in_model"
    
    with model:
        demand_id = f'DM_{metabolite_id}'
        if demand_id in model.reactions:
            model.remove_reactions([demand_id])
        
        demand_rxn = cobra.Reaction(demand_id)
        demand_rxn.add_metabolites({model.metabolites.get_by_id(metabolite_id): -1})
        model.add_reactions([demand_rxn])
        
        model.objective = demand_id
        model.objective_direction = 'max'
        
        solution = model.optimize()
        
        if solution.status == 'optimal':
            return solution.objective_value, "optimal"
        else:
            return 0.0, solution.status

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
    
    # 필수 무기물 + NAD 전구체 Exchange 열기
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

def main():
    base_path = Path(__file__).parent.parent
    new_model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_with_BCAA_cofactors_ions_nad.xml"
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    output_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_with_BCAA_cofactors_ions_nad_NNATr.xml"
    
    print("="*70)
    print("NNATr 반응 추가")
    print("="*70)
    print("\nNNATr (Nicotinate-nucleotide adenylyltransferase):")
    print("  반응식: atp_c + h_c + nicrnt_c <=> dnad_c + ppi_c")
    print("  이것이 dnad_c를 생성하는 핵심 반응!")
    print("  경로: nicrnt_c → NNATr → dnad_c → NADS1/NADS2 → nad_c → NADK → nadp_c")
    
    # 모델 로드
    print(f"\n[모델 로드]")
    new_model = cobra.io.read_sbml_model(str(new_model_path))
    ref_model = cobra.io.read_sbml_model(str(ref_model_path))
    print(f"[OK] 모델 로드 완료")
    
    # 배지 조건 설정
    new_model = setup_media_forced(new_model)
    
    # NNATr 확인
    if 'NNATr' in new_model.reactions:
        print(f"\n[OK] NNATr가 이미 모델에 있습니다!")
        ref_rxn = ref_model.reactions.get_by_id('NNATr')
        new_rxn = new_model.reactions.get_by_id('NNATr')
        print(f"  레퍼런스 반응식: {ref_rxn.reaction}")
        print(f"  신규 반응식: {new_rxn.reaction}")
        print(f"  레퍼런스 bounds: [{ref_rxn.lower_bound}, {ref_rxn.upper_bound}]")
        print(f"  신규 bounds: [{new_rxn.lower_bound}, {new_rxn.upper_bound}]")
    else:
        print(f"\n[MISSING] NNATr가 모델에 없습니다!")
        print(f"  레퍼런스 반응식: {ref_model.reactions.get_by_id('NNATr').reaction}")
        
        # NNATr 추가
        print(f"\n[NNATr 추가 중...]")
        success, msg = add_reaction_from_reference(new_model, ref_model, 'NNATr')
        if success:
            print(f"[OK] NNATr 추가 성공")
            new_rxn = new_model.reactions.get_by_id('NNATr')
            print(f"  반응식: {new_rxn.reaction}")
            print(f"  bounds: [{new_rxn.lower_bound}, {new_rxn.upper_bound}]")
            
            # 모델 저장
            cobra.io.write_sbml_model(new_model, str(output_path))
            print(f"\n[모델 저장] {output_path}")
        else:
            print(f"[ERROR] NNATr 추가 실패: {msg}")
            return new_model
    
    # NAD/NADP 생산 테스트
    print(f"\n[NAD/NADP 생산 테스트]")
    cofactors_to_test = ['dnad_c', 'nad_c', 'nadp_c']
    for met_id in cofactors_to_test:
        max_prod, status = test_metabolite_production(new_model, met_id)
        if max_prod is not None:
            status_str = "[OK]" if max_prod > 1e-6 else "[NO]"
            print(f"  {met_id:12s}: {status_str} max_production = {max_prod:.6f}")
        else:
            print(f"  {met_id:12s}: [ERROR] {status}")
    
    return new_model

if __name__ == "__main__":
    model = main()
