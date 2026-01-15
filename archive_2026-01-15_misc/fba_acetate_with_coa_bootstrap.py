#!/usr/bin/env python
"""
Acetate 기반 FBA - CoA 부트스트랩 문제 해결
초기 보조인자 허용 또는 CoA exchange 허용
"""

import cobra
import pandas as pd

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    print(f"[OK] 모델 로드: {model.id}")
    return model

def find_biomass_reaction(model):
    for name in ['Growth', 'BIOMASS']:
        try:
            return model.reactions.get_by_id(name)
        except KeyError:
            continue
    return None

def setup_acetate_medium_with_bootstrap(model):
    """Acetate medium 설정 + 부트스트랩 보조인자 허용"""
    print("\n" + "="*70)
    print("Acetate Medium 설정 (부트스트랩 포함)")
    print("="*70)
    
    # 모든 exchange 차단
    for rxn in model.exchanges:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    
    # Acetate
    try:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.lower_bound = -1000
        ex_ac.upper_bound = 1000
        print(f"[OK] Acetate: EX_ac_e")
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
    
    # 부트스트랩 보조인자 허용 (선택적)
    bootstrap_cofactors = {
        'EX_coa_e': 'CoA',
        'EX_atp_e': 'ATP',
        'EX_nad_e': 'NAD+',
    }
    
    print("\n부트스트랩 보조인자 (선택적):")
    for ex_id, name in bootstrap_cofactors.items():
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            ex_rxn.lower_bound = -10  # 소량만 허용
            ex_rxn.upper_bound = 1000
            print(f"  [OK] {name}: {ex_id} (LB=-10)")
        except KeyError:
            print(f"  [SKIP] {name}: {ex_id} (없음)")
    
    return model

def find_coa_production_reactions(model):
    """CoA 생산 반응 찾기"""
    print("\n" + "="*70)
    print("CoA 생산 반응 검색")
    print("="*70)
    
    try:
        coa_c = model.metabolites.get_by_id('coa_c')
        
        # CoA를 생성하는 반응 찾기
        coa_producing = []
        for rxn in coa_c.reactions:
            if coa_c in rxn.products:
                coa_producing.append(rxn)
        
        print(f"\nCoA 생산 반응: {len(coa_producing)}개")
        for rxn in coa_producing[:10]:
            print(f"  {rxn.id}: {rxn.reaction}")
            print(f"    LB={rxn.lower_bound}, UB={rxn.upper_bound}")
        
        # Exchange 반응 확인
        coa_exchange = []
        for rxn in model.exchanges:
            if 'coa' in rxn.id.lower():
                coa_exchange.append(rxn)
        
        if coa_exchange:
            print(f"\nCoA Exchange 반응: {len(coa_exchange)}개")
            for rxn in coa_exchange:
                print(f"  {rxn.id}: LB={rxn.lower_bound}, UB={rxn.upper_bound}")
        
        return coa_producing, coa_exchange
        
    except KeyError:
        print("[ERROR] coa_c metabolite 없음")
        return [], []

def run_fba_detailed(model, biomass_rxn):
    """상세 FBA 분석"""
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
        
        if biomass_flux > 1e-6:
            print("\n[SUCCESS] 성장 가능!")
            
            # 주요 플럭스
            print("\n주요 반응 플럭스:")
            key_reactions = ['EX_ac_e', 'ACt', 'ACS', 'CS', 'ICL', 'MALS', 
                             'SUCD', 'FUM', 'MDH', 'Growth']
            
            for rxn_id in key_reactions:
                try:
                    flux = solution.fluxes.get(rxn_id, 0)
                    if abs(flux) > 1e-8:
                        print(f"  {rxn_id}: {flux:.6f}")
                except KeyError:
                    pass
            
            # Exchange 플럭스
            print("\nExchange 플럭스 (절대값 > 0.001):")
            exchange_fluxes = []
            for rxn in model.exchanges:
                flux = solution.fluxes.get(rxn.id, 0)
                if abs(flux) > 0.001:
                    exchange_fluxes.append((rxn.id, flux))
            
            if exchange_fluxes:
                for rxn_id, flux in sorted(exchange_fluxes, key=lambda x: abs(x[1]), reverse=True):
                    print(f"  {rxn_id}: {flux:.6f}")
            
            # Yield 계산
            acetate_uptake = abs(solution.fluxes.get('EX_ac_e', 0))
            if acetate_uptake > 0:
                yield_biomass = biomass_flux / acetate_uptake
                print(f"\nBiomass yield: {yield_biomass:.6f} gDW/mmol acetate")
        else:
            print("\n[FAIL] 성장 불가 (Biomass flux = 0)")
    
    return solution

def main():
    print("="*70)
    print("Acetate 기반 FBA (CoA 부트스트랩 해결)")
    print("="*70)
    
    # 모델 로드
    model = load_model("BaseModel.xml")
    
    # Biomass 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass: {biomass_rxn.id}")
    
    # CoA 생산 반응 확인
    coa_prod, coa_ex = find_coa_production_reactions(model)
    
    # Medium 설정
    model = setup_acetate_medium_with_bootstrap(model)
    
    # FBA 실행
    solution = run_fba_detailed(model, biomass_rxn)
    
    print("\n" + "="*70)
    print("분석 완료")
    print("="*70)

if __name__ == "__main__":
    main()
