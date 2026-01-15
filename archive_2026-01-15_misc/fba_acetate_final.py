#!/usr/bin/env python
"""
Acetate 기반 FBA - 최종 버전
SUCOAACTr 반응을 활용한 부트스트랩 해결
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

def setup_acetate_medium(model):
    """Acetate minimal medium 설정"""
    print("\n" + "="*70)
    print("Acetate Minimal Medium 설정")
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
    essentials = {
        'EX_nh4_e': 'Ammonium',
        'EX_h2o_e': 'Water',
        'EX_h_e': 'Proton',
        'EX_pi_e': 'Phosphate',
        'EX_so4_e': 'Sulfate',
        'EX_k_e': 'Potassium',
        'EX_na1_e': 'Sodium',
        'EX_mg2_e': 'Magnesium',
        'EX_ca2_e': 'Calcium',
        'EX_fe2_e': 'Iron',
        'EX_mn2_e': 'Manganese',
        'EX_zn2_e': 'Zinc',
        'EX_co2_e': 'CO2',
        'EX_o2_e': 'Oxygen'
    }
    
    for ex_id, name in essentials.items():
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            if ex_id in ['EX_co2_e', 'EX_o2_e']:
                ex_rxn.lower_bound = -1000
                ex_rxn.upper_bound = 1000
            else:
                ex_rxn.lower_bound = -1000
            print(f"  [OK] {name}: {ex_id}")
        except KeyError:
            pass
    
    return model

def check_key_reactions(model):
    """핵심 반응 확인"""
    print("\n" + "="*70)
    print("핵심 반응 확인")
    print("="*70)
    
    # Acetate → Acetyl-CoA 경로
    print("\n[Acetate → Acetyl-CoA 경로]")
    ac_to_accoa = []
    
    try:
        ac_c = model.metabolites.get_by_id('ac_c')
        accoa_c = model.metabolites.get_by_id('accoa_c')
        
        for rxn in ac_c.reactions:
            if accoa_c in rxn.metabolites:
                ac_to_accoa.append(rxn)
                print(f"  {rxn.id}: {rxn.reaction}")
                print(f"    LB={rxn.lower_bound}, UB={rxn.upper_bound}")
                print(f"    유전자: {[g.id for g in rxn.genes]}")
    except KeyError:
        pass
    
    # SUCOAACTr 확인 (부트스트랩 가능)
    try:
        sucoaacr = model.reactions.get_by_id('SUCOAACTr')
        print(f"\n[부트스트랩 반응] SUCOAACTr:")
        print(f"  반응식: {sucoaacr.reaction}")
        print(f"  LB={sucoaacr.lower_bound}, UB={sucoaacr.upper_bound}")
        print(f"  이것은 CoA 없이 Acetyl-CoA를 생성할 수 있습니다!")
    except KeyError:
        print("\n[WARNING] SUCOAACTr 반응 없음")

def run_fba_complete(model, biomass_rxn):
    """완전한 FBA 분석"""
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
            
            # Exchange 플럭스
            print("\nExchange 플럭스 (절대값 > 0.001):")
            exchange_data = []
            for rxn in model.exchanges:
                flux = solution.fluxes.get(rxn.id, 0)
                if abs(flux) > 0.001:
                    exchange_data.append((rxn.id, flux))
            
            if exchange_data:
                for rxn_id, flux in sorted(exchange_data, key=lambda x: abs(x[1]), reverse=True):
                    direction = "Uptake" if flux < 0 else "Secretion"
                    print(f"  {rxn_id}: {flux:.6f} ({direction})")
            
            # 주요 경로 플럭스
            print("\n주요 경로 반응 플럭스:")
            pathway_reactions = {
                'EX_ac_e': 'Acetate Exchange',
                'ACt': 'Acetate Transport',
                'SUCOAACTr': 'Acetate-CoA Transferase',
                'ACS': 'Acetyl-CoA Synthetase',
                'CS': 'Citrate Synthase',
                'ICL': 'Isocitrate Lyase',
                'MALS': 'Malate Synthase',
                'SUCD': 'Succinate Dehydrogenase',
                'FUM': 'Fumarase',
                'MDH': 'Malate Dehydrogenase',
                'Growth': 'Biomass'
            }
            
            for rxn_id, rxn_name in pathway_reactions.items():
                try:
                    flux = solution.fluxes.get(rxn_id, 0)
                    if abs(flux) > 1e-8:
                        print(f"  {rxn_id} ({rxn_name}): {flux:.6f}")
                except KeyError:
                    pass
            
            # Yield 계산
            acetate_uptake = abs(solution.fluxes.get('EX_ac_e', 0))
            if acetate_uptake > 0:
                yield_biomass = biomass_flux / acetate_uptake
                print(f"\n성장 속도: {biomass_flux:.6f} 1/h")
                print(f"Acetate 섭취: {acetate_uptake:.6f} mmol/gDW/h")
                print(f"Biomass yield: {yield_biomass:.6f} gDW/mmol acetate")
        else:
            print("\n[FAIL] 성장 불가 (Biomass flux = 0)")
    
    return solution

def main():
    print("="*70)
    print("Acetate 기반 FBA 계산 (최종)")
    print("="*70)
    
    # 모델 로드
    model = load_model("BaseModel.xml")
    
    # Biomass 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass: {biomass_rxn.id}")
    
    # Medium 설정
    model = setup_acetate_medium(model)
    
    # 핵심 반응 확인
    check_key_reactions(model)
    
    # FBA 실행
    solution = run_fba_complete(model, biomass_rxn)
    
    print("\n" + "="*70)
    print("FBA 계산 완료")
    print("="*70)

if __name__ == "__main__":
    main()
