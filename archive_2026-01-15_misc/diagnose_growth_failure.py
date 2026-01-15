#!/usr/bin/env python
"""
성장 실패 원인 분석
단계별 상세 진단
"""

import cobra
import pandas as pd
import numpy as np

def load_model(model_path="BaseModel.xml"):
    """모델 로드"""
    print(f"모델 로드 중: {model_path}")
    model = cobra.io.read_sbml_model(model_path)
    print(f"[OK] 모델 로드 완료\n")
    return model

def setup_acetate_medium(model):
    """Acetate medium 설정"""
    for rxn in model.exchanges:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    
    essentials = {
        'EX_ac_e': (-1000, 1000),
        'EX_nh4_e': (-1000, 0),
        'EX_h2o_e': (-1000, 0),
        'EX_h_e': (-1000, 0),
        'EX_pi_e': (-1000, 0),
        'EX_so4_e': (-1000, 0),
        'EX_k_e': (-1000, 0),
        'EX_mg2_e': (-1000, 0),
        'EX_fe2_e': (-1000, 0),
        'EX_mn2_e': (-1000, 0),
        'EX_zn2_e': (-1000, 0),
        'EX_co2_e': (-1000, 1000),
        'EX_o2_e': (-1000, 1000)
    }
    
    for ex_id, bounds in essentials.items():
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            ex_rxn.lower_bound = bounds[0]
            ex_rxn.upper_bound = bounds[1]
        except KeyError:
            pass
    
    model.objective = 'Growth'
    return model

def check_acetate_pathway(model):
    """1. Acetate 경로 확인"""
    print("="*70)
    print("1. Acetate 경로 확인")
    print("="*70)
    
    # Exchange 및 Transport 확인
    try:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        print(f"[OK] EX_ac_e: {ex_ac.reaction}")
        print(f"  하한: {ex_ac.lower_bound}, 상한: {ex_ac.upper_bound}")
    except KeyError:
        print("[FAIL] EX_ac_e 없음")
    
    try:
        ac_transport = model.reactions.get_by_id('ACt')
        print(f"[OK] ACt: {ac_transport.reaction}")
    except KeyError:
        print("[FAIL] ACt transport 없음")
    
    # ACS 반응 확인
    try:
        acs = model.reactions.get_by_id('ACS')
        print(f"[OK] ACS: {acs.reaction}")
        print(f"  유전자: {[g.id for g in acs.genes]}")
        print(f"  하한: {acs.lower_bound}, 상한: {acs.upper_bound}")
        
        # 필요한 대사물질 확인
        print(f"  Reactants: {[(m.id, coeff) for m, coeff in acs.metabolites.items() if coeff < 0]}")
        print(f"  Products: {[(m.id, coeff) for m, coeff in acs.metabolites.items() if coeff > 0]}")
    except KeyError:
        print("[FAIL] ACS 반응 없음")
    
    print()

def check_tca_glyoxylate(model):
    """2. TCA/Glyoxylate 경로 확인"""
    print("="*70)
    print("2. TCA Cycle 및 Glyoxylate Shunt 경로 확인")
    print("="*70)
    
    key_reactions = {
        'CS': 'Citrate synthase',
        'ICL': 'Isocitrate lyase',
        'MALS': 'Malate synthase',
        'MDH': 'Malate dehydrogenase',
        'ME1': 'Malic enzyme (NAD)',
        'ME2': 'Malic enzyme (NADP)'
    }
    
    for rxn_id, name in key_reactions.items():
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            print(f"[OK] {rxn_id} ({name}): {rxn.reaction}")
            print(f"  가역성: {rxn.reversibility}, 하한: {rxn.lower_bound}, 상한: {rxn.upper_bound}")
        except KeyError:
            print(f"[FAIL] {rxn_id} 없음")
    
    print()

def check_biomass_reaction(model):
    """3. Biomass reaction 분석"""
    print("="*70)
    print("3. Biomass Reaction 분석")
    print("="*70)
    
    try:
        biomass = model.reactions.get_by_id('Growth')
        print(f"[OK] Growth 반응: {biomass.id}")
        print(f"  반응식: {biomass.reaction}")
        
        # 필요한 대사물질 목록 (음수 = 소비)
        required = {m.id: -coeff for m, coeff in biomass.metabolites.items() if coeff < 0}
        
        print(f"\n  필요한 대사물질 (소비): {len(required)}개")
        print("  주요 구성 요소 (상위 20개):")
        sorted_req = sorted(required.items(), key=lambda x: abs(x[1]), reverse=True)
        for met_id, coeff in sorted_req[:20]:
            print(f"    {met_id}: {coeff:.4f}")
        
        # 생성되는 대사물질 (양수)
        produced = {m.id: coeff for m, coeff in biomass.metabolites.items() if coeff > 0}
        if produced:
            print(f"\n  생성되는 대사물질: {len(produced)}개")
            for met_id, coeff in produced.items():
                print(f"    {met_id}: {coeff:.4f}")
        
    except KeyError:
        print("[FAIL] Growth 반응 없음")
    
    print()

def perform_gap_analysis(model):
    """4. Gap Analysis - 필수 대사물질 생성 가능 여부"""
    print("="*70)
    print("4. Gap Analysis: 필수 대사물질 생성 가능 여부")
    print("="*70)
    
    try:
        biomass = model.reactions.get_by_id('Growth')
        
        # Biomass에 필요한 대사물질 중 소비되는 것들 (음수)
        required_mets = {m.id: -coeff for m, coeff in biomass.metabolites.items() if coeff < 0}
        
        print(f"Biomass에 필요한 대사물질: {len(required_mets)}개\n")
        
        # 각 대사물질의 생산 가능 여부 확인
        gap_mets = []
        
        for met_id, req_coeff in list(required_mets.items())[:30]:  # 상위 30개만 확인
            try:
                met = model.metabolites.get_by_id(met_id)
                
                # 해당 대사물질을 생성하는 반응 찾기
                producing_rxns = [rxn for rxn in met.reactions 
                                 if met in rxn.products or (rxn.reversibility and met in rxn.reactants)]
                
                if len(producing_rxns) == 0:
                    gap_mets.append(met_id)
                    print(f"  [GAP] {met_id}: 생산 반응 없음")
                else:
                    print(f"  [OK] {met_id}: 생산 반응 {len(producing_rxns)}개 존재")
                    
            except KeyError:
                gap_mets.append(met_id)
                print(f"  [GAP] {met_id}: 모델에 없음")
        
        if gap_mets:
            print(f"\n[WARNING] Gap 대사물질 발견: {len(gap_mets)}개")
        else:
            print("\n[OK] 주요 대사물질 생산 경로 존재")
        
    except KeyError:
        print("[ERROR] Biomass 반응을 찾을 수 없음")
    
    print()

def check_reaction_bounds(model):
    """5. 반응 경계 조건 확인"""
    print("="*70)
    print("5. 주요 반응 경계 조건 확인")
    print("="*70)
    
    key_rxns = ['EX_ac_e', 'ACt', 'ACS', 'CS', 'ICL', 'MALS', 'MDH', 'ME1', 'ME2']
    
    for rxn_id in key_rxns:
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            print(f"{rxn_id}: 하한={rxn.lower_bound}, 상한={rxn.upper_bound}, 가역={rxn.reversibility}")
        except KeyError:
            print(f"{rxn_id}: [NOT FOUND]")
    
    print()

def perform_fva(model):
    """6. Flux Variability Analysis"""
    print("="*70)
    print("6. Flux Variability Analysis (FVA)")
    print("="*70)
    print("각 반응의 가능한 플럭스 범위 확인 중...\n")
    
    try:
        # 주요 반응에 대해 FVA 수행
        reactions_of_interest = ['EX_ac_e', 'ACt', 'ACS', 'CS', 'ICL', 'MALS', 'Growth']
        
        fva_result = cobra.flux_analysis.flux_variability_analysis(
            model, 
            reactions_of_interest,
            fraction_of_optimum=0.0  # 최적값의 0% (현재 최적값이 0이므로)
        )
        
        print("FVA 결과:")
        print(fva_result)
        
        # 각 반응이 활성화 가능한지 확인
        print("\n활성화 가능 여부:")
        for rxn_id in reactions_of_interest:
            try:
                row = fva_result.loc[rxn_id]
                min_flux = row['minimum']
                max_flux = row['maximum']
                
                if abs(min_flux) > 1e-6 or abs(max_flux) > 1e-6:
                    print(f"  {rxn_id}: 활성화 가능 (min={min_flux:.6f}, max={max_flux:.6f})")
                else:
                    print(f"  {rxn_id}: 활성화 불가능 (min={min_flux:.6f}, max={max_flux:.6f})")
            except KeyError:
                print(f"  {rxn_id}: [NOT FOUND]")
                
    except Exception as e:
        print(f"[ERROR] FVA 수행 중 오류: {e}")
        import traceback
        traceback.print_exc()
    
    print()

def check_blocked_reactions_in_pathway(model):
    """7. Acetate 경로에서 blocked reactions 확인"""
    print("="*70)
    print("7. Acetate 경로에서 Blocked Reactions 확인")
    print("="*70)
    
    # Acetate 경로의 주요 반응들
    pathway_reactions = ['EX_ac_e', 'ACt', 'ACS', 'CS', 'ACONT', 'ICL', 'MALS', 
                        'SUCD', 'FUM', 'MDH', 'ME1', 'ME2']
    
    print("경로 반응들의 blocked 여부 확인:")
    for rxn_id in pathway_reactions:
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            # 간단한 테스트: 양수 플럭스 가능 여부
            rxn.lower_bound = 0
            rxn.upper_bound = 1000
            
            # 음수 플럭스 가능 여부
            if rxn.reversibility:
                rxn.lower_bound = -1000
            
            # FBA로 테스트
            test_solution = model.optimize()
            if test_solution.status == 'optimal' and test_solution.objective_value > 0:
                flux = test_solution.fluxes.get(rxn_id, 0)
                if abs(flux) > 1e-6:
                    print(f"  {rxn_id}: 활성화 가능 (flux={flux:.6f})")
                else:
                    print(f"  {rxn_id}: 활성화 가능하지만 플럭스 0")
            else:
                print(f"  {rxn_id}: 활성화 어려움")
        except KeyError:
            print(f"  {rxn_id}: [NOT FOUND]")
        except Exception as e:
            print(f"  {rxn_id}: [ERROR] {e}")
    
    print()

def main():
    """메인 함수"""
    model = load_model("BaseModel.xml")
    model = setup_acetate_medium(model)
    
    # 단계별 진단
    check_acetate_pathway(model)
    check_tca_glyoxylate(model)
    check_biomass_reaction(model)
    perform_gap_analysis(model)
    check_reaction_bounds(model)
    perform_fva(model)
    check_blocked_reactions_in_pathway(model)
    
    print("="*70)
    print("진단 완료")
    print("="*70)

if __name__ == "__main__":
    main()

