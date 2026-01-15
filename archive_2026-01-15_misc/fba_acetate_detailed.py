#!/usr/bin/env python
"""
Acetate 기반 FBA 상세 분석
문제 해결을 위한 상세 분석
"""

import cobra
import pandas as pd

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    return model

def find_biomass_reaction(model):
    for name in ['Growth', 'BIOMASS']:
        try:
            return model.reactions.get_by_id(name)
        except KeyError:
            continue
    return None

def setup_acetate_medium(model):
    """Acetate medium 설정"""
    for rxn in model.exchanges:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    
    # Acetate
    try:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.lower_bound = -1000
        ex_ac.upper_bound = 1000
    except KeyError:
        pass
    
    # 필수 영양소
    essentials = ['EX_nh4_e', 'EX_h2o_e', 'EX_h_e', 'EX_pi_e', 'EX_so4_e',
                  'EX_k_e', 'EX_na1_e', 'EX_mg2_e', 'EX_ca2_e', 'EX_fe2_e',
                  'EX_mn2_e', 'EX_zn2_e', 'EX_co2_e', 'EX_o2_e']
    
    for ex_id in essentials:
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            if ex_id in ['EX_co2_e', 'EX_o2_e']:
                ex_rxn.lower_bound = -1000
                ex_rxn.upper_bound = 1000
            else:
                ex_rxn.lower_bound = -1000
        except KeyError:
            pass
    
    return model

def force_acetate_uptake(model):
    """Acetate uptake 강제 설정"""
    print("\n" + "="*70)
    print("Acetate Uptake 강제 테스트")
    print("="*70)
    
    try:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        # Uptake 강제
        ex_ac.lower_bound = -10  # 10 mmol/gDW/h로 제한
        ex_ac.upper_bound = 1000
        
        print(f"[OK] EX_ac_e: LB={ex_ac.lower_bound}, UB={ex_ac.upper_bound}")
        
        # Transport 확인
        try:
            ac_transport = model.reactions.get_by_id('ACt')
            ac_transport.lower_bound = -1000
            ac_transport.upper_bound = 1000
            print(f"[OK] ACt: LB={ac_transport.lower_bound}, UB={ac_transport.upper_bound}")
        except KeyError:
            print("[WARNING] ACt transport 없음")
        
    except KeyError:
        print("[ERROR] EX_ac_e 없음")

def check_cofactors(model):
    """보조인자 확인"""
    print("\n" + "="*70)
    print("보조인자 확인")
    print("="*70)
    
    cofactors = {
        'atp_c': 'ATP',
        'adp_c': 'ADP',
        'amp_c': 'AMP',
        'coa_c': 'CoA',
        'accoa_c': 'Acetyl-CoA',
        'nad_c': 'NAD+',
        'nadh_c': 'NADH',
        'nadp_c': 'NADP+',
        'nadph_c': 'NADPH'
    }
    
    for met_id, met_name in cofactors.items():
        try:
            met = model.metabolites.get_by_id(met_id)
            # Exchange 확인
            exchange_rxns = [r for r in met.reactions if 'EX' in r.id]
            if exchange_rxns:
                for rxn in exchange_rxns:
                    print(f"  {met_name} ({met_id}): {rxn.id} (LB={rxn.lower_bound}, UB={rxn.upper_bound})")
        except KeyError:
            pass

def run_fba_with_analysis(model, biomass_rxn):
    """FBA 실행 및 분석"""
    print("\n" + "="*70)
    print("FBA 최적화")
    print("="*70)
    
    model.objective = biomass_rxn.id
    
    solution = model.optimize()
    
    print(f"\n최적화 상태: {solution.status}")
    print(f"Objective value: {solution.objective_value:.6f}")
    
    if solution.status == 'optimal':
        biomass_flux = solution.fluxes.get(biomass_rxn.id, 0)
        print(f"Biomass flux: {biomass_flux:.6f} 1/h")
        
        # 주요 반응 플럭스
        print("\n주요 반응 플럭스:")
        key_reactions = ['EX_ac_e', 'ACt', 'ACS', 'R_ACS', 'CS', 'ICL', 'MALS', 'ICDHx', 
                         'AKGDH', 'SUCD', 'FUM', 'MDH', 'SUCOAS', 'Growth']
        
        for rxn_id in key_reactions:
            try:
                flux = solution.fluxes.get(rxn_id, 0)
                if abs(flux) > 1e-8:
                    rxn = model.reactions.get_by_id(rxn_id)
                    print(f"  {rxn_id}: {flux:.6f} | {rxn.reaction}")
            except KeyError:
                pass
        
        # Exchange 플럭스
        print("\nExchange 플럭스 (절대값 > 1e-6):")
        exchange_fluxes = []
        for rxn in model.exchanges:
            flux = solution.fluxes.get(rxn.id, 0)
            if abs(flux) > 1e-6:
                exchange_fluxes.append((rxn.id, flux))
        
        if exchange_fluxes:
            for rxn_id, flux in sorted(exchange_fluxes, key=lambda x: abs(x[1]), reverse=True):
                print(f"  {rxn_id}: {flux:.6f}")
        else:
            print("  [없음]")
    
    return solution

def test_bootstrap(model):
    """부트스트랩 테스트 - 초기 CoA 생성 가능 여부"""
    print("\n" + "="*70)
    print("부트스트랩 테스트")
    print("="*70)
    
    try:
        coa_c = model.metabolites.get_by_id('coa_c')
        
        # CoA 생산 가능 여부 확인
        test_rxn = cobra.Reaction('TEST_COA')
        test_rxn.add_metabolites({coa_c: 1})
        test_rxn.lower_bound = 0
        test_rxn.upper_bound = 1000
        
        model.add_reactions([test_rxn])
        model.objective = test_rxn.id
        
        solution = model.optimize()
        
        model.remove_reactions([test_rxn])
        
        if solution.status == 'optimal' and solution.objective_value > 1e-6:
            print("[OK] CoA 생산 가능")
            return True
        else:
            print("[FAIL] CoA 생산 불가")
            return False
            
    except KeyError:
        print("[ERROR] coa_c metabolite 없음")
        return False

def main():
    print("="*70)
    print("Acetate 기반 FBA 상세 분석")
    print("="*70)
    
    # 모델 로드
    model = load_model("BaseModel.xml")
    print(f"[OK] 모델 로드: {model.id}")
    
    # Biomass 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass: {biomass_rxn.id}")
    
    # Medium 설정
    model = setup_acetate_medium(model)
    
    # Acetate uptake 강제
    force_acetate_uptake(model)
    
    # 보조인자 확인
    check_cofactors(model)
    
    # 부트스트랩 테스트
    can_bootstrap = test_bootstrap(model)
    
    # FBA 실행
    solution = run_fba_with_analysis(model, biomass_rxn)
    
    print("\n" + "="*70)
    print("분석 완료")
    print("="*70)
    
    if solution and solution.status == 'optimal':
        if solution.objective_value > 1e-6:
            print("\n[SUCCESS] 성장 가능!")
        else:
            print("\n[FAIL] 성장 불가 - 추가 진단 필요")

if __name__ == "__main__":
    main()
