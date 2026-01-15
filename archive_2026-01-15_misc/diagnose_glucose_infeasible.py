#!/usr/bin/env python
"""
포도당 infeasible 원인 진단
왜 포도당만으로 infeasible인지 상세 분석
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

def diagnose_glucose_infeasible(model, biomass_rxn):
    """포도당 infeasible 원인 진단"""
    print("="*70)
    print("포도당 Infeasible 원인 진단")
    print("="*70)
    
    # 1. 포도당 transport 경로 확인
    print("\n[1] 포도당 Transport 경로 확인:")
    print("-" * 70)
    
    glucose_transports = ['GLCabc', 'GLCabcpp', 'GLCpts', 'GLCt2rpp', 'GLCtex']
    
    for rxn_id in glucose_transports:
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            print(f"  [OK] {rxn_id}: {rxn.name}")
            print(f"    {rxn.reaction}")
        except KeyError:
            # 패턴으로 찾기
            for rxn in model.reactions:
                if 'glc' in rxn.id.lower() and ('abc' in rxn.id.lower() or 'pts' in rxn.id.lower() or 't' in rxn.id.lower()):
                    if rxn not in [model.reactions.get_by_id(id) for id in glucose_transports if id in model.reactions]:
                        print(f"  [FOUND] {rxn.id}: {rxn.name}")
                        print(f"    {rxn.reaction}")
                        break
    
    # 2. 포도당이 실제로 세포 내로 들어가는지 확인
    print("\n[2] 포도당 세포 내 유입 확인:")
    print("-" * 70)
    
    # 포도당 exchange 설정
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    try:
        ex_glc = model.reactions.get_by_id('EX_glc__D_e')
        ex_glc.lower_bound = -100
        ex_glc.upper_bound = 1000
    except KeyError:
        pass
    
    # 필수 영양소 (최소한만)
    minimal_nutrients = ['EX_nh4_e', 'EX_h2o_e', 'EX_pi_e', 'EX_o2_e']
    
    for ex_id in minimal_nutrients:
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            if ex_id == 'EX_o2_e':
                ex_rxn.lower_bound = -1000
                ex_rxn.upper_bound = 1000
            else:
                ex_rxn.lower_bound = -1000
        except KeyError:
            pass
    
    # 포도당 세포 내 농도 테스트
    try:
        glc__D_c = model.metabolites.get_by_id('glc__D_c')
        
        dm_glc = cobra.Reaction('DM_glc__D_c')
        dm_glc.name = 'Glucose demand'
        dm_glc.lower_bound = 0
        dm_glc.upper_bound = 1000
        dm_glc.add_metabolites({glc__D_c: -1})
        model.add_reactions([dm_glc])
        
        model.objective = dm_glc.id
        solution = model.optimize()
        
        model.remove_reactions([dm_glc])
        
        if solution.status == 'optimal' and solution.objective_value > 1e-6:
            print(f"  [OK] 포도당 세포 내 유입 가능: {solution.objective_value:.6f} mmol/gDCW/h")
        else:
            print(f"  [FAIL] 포도당 세포 내 유입 불가능 ({solution.status})")
            print("    → Transport 경로 문제 가능성")
    except KeyError:
        print("  [ERROR] glc__D_c metabolite 없음")
    
    # 3. 단계별 영양소 추가 테스트
    print("\n[3] 단계별 영양소 추가 테스트:")
    print("-" * 70)
    
    nutrient_groups = [
        {
            'name': '기본 영양소',
            'nutrients': ['EX_nh4_e', 'EX_h2o_e', 'EX_h_e', 'EX_pi_e', 'EX_o2_e']
        },
        {
            'name': '기본 영양소 + CO2',
            'nutrients': ['EX_nh4_e', 'EX_h2o_e', 'EX_h_e', 'EX_pi_e', 'EX_o2_e', 'EX_co2_e']
        },
        {
            'name': '기본 영양소 + CO2 + Sulfate',
            'nutrients': ['EX_nh4_e', 'EX_h2o_e', 'EX_h_e', 'EX_pi_e', 'EX_o2_e', 'EX_co2_e', 'EX_so4_e']
        },
        {
            'name': '기본 영양소 + 무기염',
            'nutrients': ['EX_nh4_e', 'EX_h2o_e', 'EX_h_e', 'EX_pi_e', 'EX_o2_e', 'EX_co2_e', 'EX_so4_e',
                         'EX_k_e', 'EX_na1_e', 'EX_mg2_e', 'EX_ca2_e']
        },
        {
            'name': '기본 영양소 + 무기염 + 미량원소',
            'nutrients': ['EX_nh4_e', 'EX_h2o_e', 'EX_h_e', 'EX_pi_e', 'EX_o2_e', 'EX_co2_e', 'EX_so4_e',
                         'EX_k_e', 'EX_na1_e', 'EX_mg2_e', 'EX_ca2_e', 'EX_fe2_e', 'EX_mn2_e', 'EX_zn2_e']
        }
    ]
    
    for group in nutrient_groups:
        # 모든 exchange 초기화
        for rxn in model.exchanges:
            rxn.upper_bound = 0
            rxn.lower_bound = 0
        
        # 포도당 설정
        try:
            ex_glc = model.reactions.get_by_id('EX_glc__D_e')
            ex_glc.lower_bound = -100
            ex_glc.upper_bound = 1000
        except KeyError:
            pass
        
        # 영양소 추가
        for ex_id in group['nutrients']:
            try:
                ex_rxn = model.reactions.get_by_id(ex_id)
                if ex_id in ['EX_co2_e', 'EX_o2_e']:
                    ex_rxn.lower_bound = -1000
                    ex_rxn.upper_bound = 1000
                else:
                    ex_rxn.lower_bound = -1000
            except KeyError:
                pass
        
        # 최적화
        model.objective = biomass_rxn.id
        solution = model.optimize()
        
        print(f"\n  [{group['name']}]:")
        print(f"    상태: {solution.status}")
        
        if solution.status == 'optimal':
            biomass_flux = solution.fluxes.get(biomass_rxn.id, 0)
            if biomass_flux > 1e-6:
                print(f"    [SUCCESS] 생장 가능: {biomass_flux:.6f} 1/h")
                print(f"    → 이 조합으로 생장 가능!")
                break
            else:
                print(f"    Biomass flux: {biomass_flux:.6f} 1/h (생장 불가)")
        else:
            print(f"    [FAIL] 최적화 실패: {solution.status}")
    
    # 4. 무제한 영양소와 비교
    print("\n[4] 무제한 영양소 상태:")
    print("-" * 70)
    
    for rxn in model.exchanges:
        rxn.lower_bound = -1000
        rxn.upper_bound = 1000
    
    model.objective = biomass_rxn.id
    solution_unlimited = model.optimize()
    
    if solution_unlimited.status == 'optimal':
        biomass_unlimited = solution_unlimited.objective_value
        print(f"  [SUCCESS] Biomass flux: {biomass_unlimited:.6f} 1/h")
        
        # 포도당 사용량 확인
        try:
            glc_flux = solution_unlimited.fluxes.get('EX_glc__D_e', 0)
            print(f"  Glucose flux: {glc_flux:.6f} mmol/gDCW/h")
            if glc_flux < 0:
                print(f"  → 포도당을 섭취하고 있음 (절대값: {abs(glc_flux):.6f})")
        except:
            pass
    
    # 5. Biomass 구성 요소 생산 가능 여부 확인
    print("\n[5] Biomass 구성 요소 생산 가능 여부 (포도당 기반):")
    print("-" * 70)
    
    # 포도당 + 필수 영양소 설정
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    try:
        ex_glc = model.reactions.get_by_id('EX_glc__D_e')
        ex_glc.lower_bound = -100
        ex_glc.upper_bound = 1000
    except KeyError:
        pass
    
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
    
    # 주요 구성 요소 테스트
    key_components = ['atp_c', 'nad_c', 'coa_c', 'ala__L_c', 'gly_c', 'g6p_c', 'pep_c']
    
    print("\n주요 구성 요소 생산 가능 여부:")
    component_status = []
    
    for met_id in key_components:
        try:
            met = model.metabolites.get_by_id(met_id)
            
            test_rxn = cobra.Reaction(f'TEST_{met_id}')
            test_rxn.add_metabolites({met: 1})
            test_rxn.lower_bound = 0
            test_rxn.upper_bound = 1000
            
            model.add_reactions([test_rxn])
            model.objective = test_rxn.id
            
            solution = model.optimize()
            model.remove_reactions([test_rxn])
            
            can_produce = solution.status == 'optimal' and solution.objective_value > 1e-6
            
            component_status.append({
                'Component': met_id,
                'Can_Produce': can_produce,
                'Status': solution.status,
                'Max_Flux': solution.objective_value if solution.status == 'optimal' else 0
            })
            
            status_icon = "[OK]" if can_produce else "[FAIL]"
            print(f"  {status_icon} {met_id}")
            
        except KeyError:
            component_status.append({
                'Component': met_id,
                'Can_Produce': False,
                'Status': 'Metabolite missing',
                'Max_Flux': 0
            })
            print(f"  [MISSING] {met_id}")
    
    # 결과 저장
    if component_status:
        df_components = pd.DataFrame(component_status)
        df_components.to_csv('glucose_component_production_status.csv', index=False)
        print(f"\n[OK] 구성 요소 생산 상태 저장: glucose_component_production_status.csv")
    
    return component_status

def main():
    print("="*70)
    print("포도당 Infeasible 원인 진단")
    print("="*70)
    
    model = load_model("BaseModel.xml")
    
    # Biomass 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass: {biomass_rxn.id}")
    
    # 진단 수행
    component_status = diagnose_glucose_infeasible(model, biomass_rxn)
    
    # 최종 요약
    print("\n" + "="*70)
    print("최종 요약")
    print("="*70)
    
    if component_status:
        can_produce_count = sum(1 for s in component_status if s['Can_Produce'])
        total_count = len(component_status)
        
        print(f"\n포도당 기반 구성 요소 생산:")
        print(f"  생산 가능: {can_produce_count}/{total_count}")
        
        if can_produce_count < total_count:
            print("\n생산 불가능한 구성 요소:")
            for s in component_status:
                if not s['Can_Produce']:
                    print(f"  - {s['Component']}: {s['Status']}")
    
    print("\n" + "="*70)

if __name__ == "__main__":
    main()
