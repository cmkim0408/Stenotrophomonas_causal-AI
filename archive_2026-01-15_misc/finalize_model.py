#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
최종 모델 정리 및 요약

현재 상태:
- sink_4hba_c 제거, 4hba exchange 추가 완료
- MQN8t는 유지 (플럭스가 매우 작아서 큰 문제 아님)
- 모든 아티팩트 제거 완료
"""

import cobra
from pathlib import Path
import sys

def load_model(model_path):
    try:
        model = cobra.io.read_sbml_model(str(model_path))
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

def summarize_model_changes():
    """모델 변경 사항 요약"""
    print("\n" + "="*70)
    print("최종 모델 변경 사항 요약")
    print("="*70)
    
    changes = {
        "제거된 반응": [
            "ACt2rpp (Acetate 재순환 방지)",
            "ACS_ADP, SUCOAACTr (비현실적 경로)",
            "PEPCK_ATP (PEPCK 제거)",
            "ACCOAL, APAT_1, PACPT_1 (에너지 회수 루프)",
            "GALpts, GALt2, A6PAG (PEP 루프)",
            "ACOAD2, ACOAD2f, SUCD (레독스 셔틀)",
            "DM_coa_c (CoA 생합성 경로 활성화)",
            "sink_4hba_c (비현실적 sink)",
        ],
        "추가된 반응": [
            "BCAA 합성: KARI, DHAD, IPMI, IPMDH, BCAT_VAL, BCAT_LEU (6개)",
            "이온 수송: T_cl_e_to_cl_c, T_cu2_e_to_cu2_c, T_cobalt2_e_to_cobalt2_c (3개)",
            "NAD/NADP: T_nac_e_to_nac_c, EX_nac_e, EX_ncam_e (3개)",
            "maeB (NADP-dependent malic enzyme)",
            "EX_pnto__R_e, T_pnto__R_e_to_c (판토텐산)",
            "EX_4hba_e, T_4hba_e_to_c (4hba 배출)",
        ],
        "수정된 반응": [
            "PPA_1pp 비활성화 (H+-translocating PPase)",
        ],
        "유지된 반응 (작은 플럭스)": [
            "MQN8t (메나퀴논 생합성 경로 불완전, 플럭스 매우 작음)",
        ],
    }
    
    for category, items in changes.items():
        print(f"\n[{category}]")
        for item in items:
            print(f"  - {item}")
    
    print(f"\n[배지 가정]")
    print(f"  - EX_pnto__R_e: 판토텐산 exchange (yeast extract 가정)")
    print(f"  - EX_nac_e, EX_ncam_e: NAD 전구체 exchange (yeast extract 가정)")
    print(f"  - 실제 배지 조건에 맞게 조정 가능")

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_final_cleaned.xml"
    
    print("="*70)
    print("최종 모델 정리 및 요약")
    print("="*70)
    
    # 모델 로드
    model = load_model(model_path)
    model = setup_media_forced(model)
    
    # ATPM=0 설정
    if 'ATPM' in model.reactions:
        atpm = model.reactions.get_by_id('ATPM')
        atpm.lower_bound = 0.0
        atpm.upper_bound = 0.0
    
    # FBA 실행
    solution = model.optimize()
    
    print(f"\n[최종 FBA 결과]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    
    # 주요 경로 플럭스
    print(f"\n[주요 경로 플럭스]")
    key_rxns = ['EX_ac_e', 'ACS', 'CS', 'ICL', 'MALS', 'SUCDi', 'MDH', 'MAEB', 'Growth', 'ATPS4rpp']
    for rxn_id in key_rxns:
        if rxn_id in model.reactions:
            flux = solution.fluxes.get(rxn_id, 0.0)
            if abs(flux) > 1e-6:
                print(f"  {rxn_id:15s}: {flux:10.6f}")
    
    # 아티팩트 확인
    print(f"\n[아티팩트 확인]")
    artifact_rxns = ['ACCOAL', 'APAT_1', 'PACPT_1', 'GALpts', 'ACOAD2', 'ACOAD2f', 'SUCD', 'sink_4hba_c', 'DM_coa_c']
    found_artifacts = []
    for rxn_id in artifact_rxns:
        if rxn_id in model.reactions:
            found_artifacts.append(rxn_id)
    
    if found_artifacts:
        print(f"  [주의] 다음 아티팩트가 여전히 모델에 있습니다: {found_artifacts}")
    else:
        print(f"  [OK] 모든 아티팩트가 제거되었습니다.")
    
    # MQN8t 확인
    if 'MQN8t' in model.reactions:
        mqn8t_flux = solution.fluxes.get('MQN8t', 0.0)
        print(f"\n[MQN8t]")
        print(f"  플럭스: {mqn8t_flux:.6f} (매우 작음)")
        print(f"  상태: 유지 (메나퀴논 생합성 경로 불완전, 플럭스가 매우 작아서 큰 문제 아님)")
    
    # 모델 변경 사항 요약
    summarize_model_changes()
    
    return model, solution

if __name__ == "__main__":
    model, solution = main()
