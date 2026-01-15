#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
판토텐산 Exchange 추가 및 DM_coa_c 제거

1. EX_pnto__R_e 추가 (판토텐산 exchange)
2. DM_coa_c 제거
3. CoA 생합성 경로가 작동하는지 확인
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

def add_pantothenate_exchange(model):
    """판토텐산 Exchange 추가"""
    print("\n" + "="*70)
    print("판토텐산 Exchange 추가")
    print("="*70)
    
    # pnto__R_e metabolite 확인/추가
    if 'pnto__R_e' not in model.metabolites:
        pnto_e = cobra.Metabolite('pnto__R_e', name='D-Pantothenate (extracellular)', compartment='e')
        model.add_metabolites([pnto_e])
        print(f"[추가] pnto__R_e metabolite 추가")
    
    # EX_pnto__R_e 확인/추가
    if 'EX_pnto__R_e' not in model.reactions:
        ex_pnto = cobra.Reaction('EX_pnto__R_e')
        ex_pnto.name = "Pantothenate exchange"
        ex_pnto.add_metabolites({model.metabolites.get_by_id('pnto__R_e'): -1})
        ex_pnto.lower_bound = -1000.0  # Uptake 허용
        ex_pnto.upper_bound = 1000.0
        model.add_reactions([ex_pnto])
        print(f"[추가] EX_pnto__R_e 반응 추가")
        print(f"  반응식: {ex_pnto.reaction}")
        print(f"  bounds: [{ex_pnto.lower_bound}, {ex_pnto.upper_bound}]")
    else:
        ex_pnto = model.reactions.get_by_id('EX_pnto__R_e')
        ex_pnto.lower_bound = -1000.0
        ex_pnto.upper_bound = 1000.0
        print(f"[수정] EX_pnto__R_e bounds: [{ex_pnto.lower_bound}, {ex_pnto.upper_bound}]")
    
    # 판토텐산 transport 확인/추가
    if 'pnto__R_c' not in model.metabolites:
        pnto_c = cobra.Metabolite('pnto__R_c', name='D-Pantothenate (cytosol)', compartment='c')
        model.add_metabolites([pnto_c])
        print(f"[추가] pnto__R_c metabolite 추가")
    
    # Transport 반응 확인
    transport_id = 'T_pnto__R_e_to_c'
    if transport_id not in model.reactions:
        if 'pnto__R_e' in model.metabolites and 'pnto__R_c' in model.metabolites:
            transport = cobra.Reaction(transport_id)
            transport.name = "Pantothenate transport (e<->c)"
            transport.add_metabolites({
                model.metabolites.get_by_id('pnto__R_e'): -1,
                model.metabolites.get_by_id('pnto__R_c'): 1,
            })
            transport.lower_bound = -1000.0
            transport.upper_bound = 1000.0
            model.add_reactions([transport])
            print(f"[추가] {transport_id} 반응 추가")
            print(f"  반응식: {transport.reaction}")
    
    return True

def remove_dm_coa(model):
    """DM_coa_c 제거"""
    print("\n" + "="*70)
    print("DM_coa_c 제거")
    print("="*70)
    
    if 'DM_coa_c' not in model.reactions:
        print("[INFO] DM_coa_c가 모델에 없습니다.")
        return False
    
    dm_coa = model.reactions.get_by_id('DM_coa_c')
    print(f"[제거] DM_coa_c")
    print(f"  반응식: {dm_coa.reaction}")
    print(f"  이유: CoA 생합성 경로가 작동하도록 하기 위해")
    
    model.remove_reactions(['DM_coa_c'])
    print(f"[OK] DM_coa_c 제거 완료")
    
    return True

def test_model(model):
    """모델 테스트"""
    print("\n" + "="*70)
    print("모델 테스트")
    print("="*70)
    
    # ATPM=0 설정
    if 'ATPM' in model.reactions:
        atpm = model.reactions.get_by_id('ATPM')
        atpm.lower_bound = 0.0
        atpm.upper_bound = 0.0
    
    solution = model.optimize()
    
    print(f"\n[FBA 결과]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    
    # DM_coa_c 확인
    if 'DM_coa_c' in model.reactions:
        dm_coa_flux = solution.fluxes.get('DM_coa_c', 0.0)
        print(f"  DM_coa_c 플럭스: {dm_coa_flux:.6f} (여전히 존재)")
    else:
        print(f"  DM_coa_c: 제거됨")
    
    # 판토텐산 관련 플럭스 확인
    if 'EX_pnto__R_e' in model.reactions:
        ex_pnto_flux = solution.fluxes.get('EX_pnto__R_e', 0.0)
        print(f"  EX_pnto__R_e 플럭스: {ex_pnto_flux:.6f}")
    
    # CoA 생산 경로 확인
    coa_producing_rxns = []
    for rxn in model.reactions:
        if 'coa_c' in [m.id for m in rxn.metabolites]:
            coa_coeff = rxn.metabolites.get(model.metabolites.get_by_id('coa_c'), 0)
            if coa_coeff > 0:  # CoA 생성
                flux = solution.fluxes.get(rxn.id, 0.0)
                if abs(flux) > 1e-6:
                    coa_producing_rxns.append((rxn.id, flux, coa_coeff))
    
    if coa_producing_rxns:
        print(f"\n[CoA 생산 반응] (플럭스 > 1e-6)")
        for rxn_id, flux, coeff in sorted(coa_producing_rxns, key=lambda x: abs(x[1]), reverse=True)[:5]:
            coa_prod = flux * coeff
            print(f"  {rxn_id:20s}: 플럭스 {flux:10.6f}, CoA 생산 {coa_prod:.6f}")
    
    return solution

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_no_redox_shuttle.xml"
    output_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel_coa_synthesis_fixed.xml"
    
    print("="*70)
    print("판토텐산 Exchange 추가 및 DM_coa_c 제거")
    print("="*70)
    
    # 모델 로드
    model = load_model(model_path)
    model = setup_media_forced(model)
    
    # ATPM=0 설정
    if 'ATPM' in model.reactions:
        atpm = model.reactions.get_by_id('ATPM')
        atpm.lower_bound = 0.0
        atpm.upper_bound = 0.0
    
    # 제거 전 FBA
    solution_before = model.optimize()
    print(f"\n[수정 전 FBA]")
    print(f"  상태: {solution_before.status}")
    print(f"  성장률: {solution_before.objective_value:.6f}")
    
    if 'DM_coa_c' in model.reactions:
        dm_coa_flux_before = solution_before.fluxes.get('DM_coa_c', 0.0)
        print(f"  DM_coa_c 플럭스: {dm_coa_flux_before:.6f}")
    
    # 판토텐산 Exchange 추가
    add_pantothenate_exchange(model)
    
    # DM_coa_c 제거
    remove_dm_coa(model)
    
    # 수정 후 FBA
    solution_after = test_model(model)
    
    # 결과 요약
    print("\n" + "="*70)
    print("결과 요약")
    print("="*70)
    
    growth_before = solution_before.objective_value
    growth_after = solution_after.objective_value
    
    print(f"\n[성장률 비교]")
    print(f"  수정 전: {growth_before:.6f}")
    print(f"  수정 후: {growth_after:.6f}")
    print(f"  차이: {growth_after - growth_before:.6f}")
    
    if growth_after > 1e-6:
        print(f"\n[결론] 판토텐산 Exchange 추가 및 DM_coa_c 제거 후에도 성장 가능")
        print(f"       (성장률: {growth_after:.6f})")
        
        # 모델 저장
        cobra.io.write_sbml_model(model, str(output_path))
        print(f"\n[모델 저장] {output_path}")
    else:
        print(f"\n[주의] 수정 후 성장 불가 (성장률: {growth_after:.6f})")
        print(f"       -> CoA 생합성 경로가 완전하지 않을 수 있음")
    
    return model, solution_before, solution_after

if __name__ == "__main__":
    model, sol_before, sol_after = main()
