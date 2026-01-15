#!/usr/bin/env python
"""
포도당 완전 산화 경로 테스트
포도당 → Glycolysis → TCA → NADH → ETC → ATP 경로 확인
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

def test_glucose_to_nadh_pathway(model):
    """포도당 → NADH 경로 테스트"""
    print("="*70)
    print("포도당 완전 산화 경로 테스트")
    print("="*70)
    
    # 모든 exchange 초기화
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    # 포도당 설정 (단순 확산 경로 사용)
    # GLCtex는 ATP 필요 없음
    try:
        ex_glc = model.reactions.get_by_id('EX_glc__D_e')
        ex_glc.lower_bound = -100
        ex_glc.upper_bound = 1000
        print(f"[OK] 포도당 설정: EX_glc__D_e")
    except KeyError:
        pass
    
    # 필수 영양소 (최소한만)
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
    
    # 1. 포도당 → Pyruvate 경로 확인
    print("\n[1] 포도당 → Pyruvate 경로 확인:")
    print("-" * 70)
    
    glycolysis_pathway = [
        ('EX_glc__D_e', 'Glucose Exchange'),
        ('GLCtex', 'Glucose transport (diffusion)'),
        ('GLCabc', 'Glucose ABC transport'),
        ('HEX1', 'Hexokinase'),
        ('PGI', 'Phosphoglucose Isomerase'),
        ('PFK', 'Phosphofructokinase'),
        ('FBA', 'Fructose-bisphosphate Aldolase'),
        ('TPI', 'Triosephosphate Isomerase'),
        ('GAPD', 'Glyceraldehyde-3-phosphate Dehydrogenase'),
        ('PGK', 'Phosphoglycerate Kinase'),
        ('PGM', 'Phosphoglycerate Mutase'),
        ('ENO', 'Enolase'),
        ('PYK', 'Pyruvate Kinase'),
        ('PDH', 'Pyruvate Dehydrogenase')
    ]
    
    print("\nGlycolysis 경로 반응 확인:")
    for rxn_id, rxn_name in glycolysis_pathway:
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            print(f"  [OK] {rxn_id}: {rxn.name}")
        except KeyError:
            print(f"  [MISSING] {rxn_id}: {rxn_name}")
    
    # Pyruvate 생산 테스트
    try:
        pyr_c = model.metabolites.get_by_id('pyr_c')
        
        dm_pyr = cobra.Reaction('DM_pyr_c')
        dm_pyr.name = 'Pyruvate demand'
        dm_pyr.lower_bound = 0
        dm_pyr.upper_bound = 1000
        dm_pyr.add_metabolites({pyr_c: -1})
        model.add_reactions([dm_pyr])
        
        model.objective = dm_pyr.id
        solution = model.optimize()
        
        model.remove_reactions([dm_pyr])
        
        if solution.status == 'optimal' and solution.objective_value > 1e-6:
            print(f"\n  [SUCCESS] Pyruvate 생산 가능: {solution.objective_value:.6f} mmol/gDCW/h")
        else:
            print(f"\n  [FAIL] Pyruvate 생산 불가능 ({solution.status})")
    except KeyError:
        print("\n  [ERROR] pyr_c metabolite 없음")
    
    # 2. TCA Cycle → NADH 경로 확인
    print("\n[2] TCA Cycle → NADH 생산 경로 확인:")
    print("-" * 70)
    
    tca_nadh_producers = {
        'ICDHx': 'Isocitrate dehydrogenase (NAD)',
        'AKGDH': '2-Oxoglutarate dehydrogenase',
        'MDH': 'Malate dehydrogenase'
    }
    
    print("\nTCA Cycle NADH 생산 반응:")
    for rxn_id, rxn_name in tca_nadh_producers.items():
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            print(f"  [OK] {rxn_id}: {rxn.name}")
            print(f"    {rxn.reaction}")
        except KeyError:
            print(f"  [MISSING] {rxn_id}: {rxn_name}")
    
    # NADH 생산 테스트
    try:
        nadh_c = model.metabolites.get_by_id('nadh_c')
        
        dm_nadh = cobra.Reaction('DM_nadh_c')
        dm_nadh.name = 'NADH demand'
        dm_nadh.lower_bound = 0
        dm_nadh.upper_bound = 1000
        dm_nadh.add_metabolites({nadh_c: -1})
        model.add_reactions([dm_nadh])
        
        model.objective = dm_nadh.id
        solution = model.optimize()
        
        model.remove_reactions([dm_nadh])
        
        if solution.status == 'optimal':
            nadh_flux = solution.objective_value
            
            if nadh_flux > 1e-6:
                print(f"\n  [SUCCESS] NADH 생산 가능: {nadh_flux:.6f} mmol/gDCW/h")
                
                # TCA 경로 플럭스 확인
                print("\n  TCA Cycle 플럭스:")
                tca_rxns = ['CS', 'ACONT', 'ICDHx', 'AKGDH', 'SUCOAS', 'SUCD', 'FUM', 'MDH']
                for rxn_id in tca_rxns:
                    try:
                        flux = solution.fluxes.get(rxn_id, 0)
                        if abs(flux) > 1e-8:
                            print(f"    {rxn_id}: {flux:.6f}")
                    except KeyError:
                        pass
            else:
                print(f"\n  [FAIL] NADH 생산 불가능 (flux: {nadh_flux:.6f})")
        else:
            print(f"\n  [FAIL] NADH 생산 불가능 ({solution.status})")
    except KeyError:
        print("\n  [ERROR] nadh_c metabolite 없음")
    
    # 3. NADH → ATP (ETC 경로) 확인
    print("\n[3] NADH → ATP (ETC 경로) 확인:")
    print("-" * 70)
    
    etc_pathway = {
        'NADH16': 'NADH dehydrogenase (Complex I)',
        'SUCD': 'Succinate dehydrogenase (Complex II)',
        'QCR': 'Cytochrome bc1 complex (Complex III)',
        'CYO3': 'Cytochrome c oxidase (Complex IV)',
        'ATPS4rpp': 'ATP synthase (Complex V)'
    }
    
    print("\nETC Complex 반응:")
    etc_found = []
    for rxn_id, rxn_name in etc_pathway.items():
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            etc_found.append(rxn_id)
            print(f"  [OK] {rxn_id}: {rxn.name}")
            print(f"    {rxn.reaction}")
            
            # h_p 생성/소비 확인
            if 'h_p' in rxn.reaction:
                try:
                    h_p = model.metabolites.get_by_id('h_p')
                    h_p_coeff = rxn.metabolites.get(h_p, 0)
                    if h_p_coeff > 0:
                        print(f"    [OK] h_p 생성: {h_p_coeff}")
                    elif h_p_coeff < 0:
                        print(f"    [OK] h_p 소비: {abs(h_p_coeff)}")
                except:
                    pass
        except KeyError:
            print(f"  [MISSING] {rxn_id}: {rxn_name}")
    
    # ATP 생산 테스트 (NADH 생산 경로 포함)
    try:
        atp_c = model.metabolites.get_by_id('atp_c')
        
        dm_atp = cobra.Reaction('DM_atp_c')
        dm_atp.name = 'ATP demand'
        dm_atp.lower_bound = 0
        dm_atp.upper_bound = 1000
        dm_atp.add_metabolites({atp_c: -1})
        model.add_reactions([dm_atp])
        
        model.objective = dm_atp.id
        solution = model.optimize()
        
        model.remove_reactions([dm_atp])
        
        if solution.status == 'optimal':
            atp_flux = solution.objective_value
            
            if atp_flux > 1e-6:
                print(f"\n  [SUCCESS] ATP 생산 가능: {atp_flux:.6f} mmol/gDCW/h")
                
                # 전체 경로 플럭스 확인
                print("\n  포도당 완전 산화 경로 플럭스:")
                
                # Exchange
                glucose_uptake = abs(solution.fluxes.get('EX_glc__D_e', 0))
                if glucose_uptake > 1e-8:
                    print(f"    EX_glc__D_e: {solution.fluxes.get('EX_glc__D_e', 0):.6f}")
                
                # Glycolysis
                print("\n    [Glycolysis]:")
                gly_rxns = ['HEX1', 'PGI', 'PFK', 'GAPD', 'PGK', 'PYK']
                for rxn_id in gly_rxns:
                    try:
                        flux = solution.fluxes.get(rxn_id, 0)
                        if abs(flux) > 1e-8:
                            print(f"      {rxn_id}: {flux:.6f}")
                    except KeyError:
                        pass
                
                # PDH
                try:
                    pdh_flux = solution.fluxes.get('PDH', 0)
                    if abs(pdh_flux) > 1e-8:
                        print(f"\n    [Pyruvate → Acetyl-CoA]:")
                        print(f"      PDH: {pdh_flux:.6f}")
                except:
                    pass
                
                # TCA
                print("\n    [TCA Cycle]:")
                tca_rxns = ['CS', 'ACONT', 'ICDHx', 'AKGDH', 'SUCOAS', 'SUCD', 'FUM', 'MDH']
                for rxn_id in tca_rxns:
                    try:
                        flux = solution.fluxes.get(rxn_id, 0)
                        if abs(flux) > 1e-8:
                            print(f"      {rxn_id}: {flux:.6f}")
                    except KeyError:
                        pass
                
                # ETC
                print("\n    [ETC]:")
                for rxn_id in etc_found:
                    try:
                        flux = solution.fluxes.get(rxn_id, 0)
                        if abs(flux) > 1e-8:
                            print(f"      {rxn_id}: {flux:.6f}")
                    except KeyError:
                        pass
                
                # ATP 생산 경로
                print("\n    [ATP 생산]:")
                atp_producers = ['ATPS4rpp', 'PYK', 'PGK']
                for rxn_id in atp_producers:
                    try:
                        flux = solution.fluxes.get(rxn_id, 0)
                        if abs(flux) > 1e-8:
                            print(f"      {rxn_id}: {flux:.6f}")
                    except KeyError:
                        pass
                
                return True, solution
            else:
                print(f"\n  [FAIL] ATP 생산 불가능 (flux: {atp_flux:.6f})")
                
                # 차단 원인 분석
                print("\n  차단 원인 분석:")
                
                # NADH 생산 가능 여부
                try:
                    nadh_test = cobra.Reaction('TEST_nadh')
                    nadh_test.add_metabolites({nadh_c: -1})
                    nadh_test.lower_bound = 0
                    nadh_test.upper_bound = 1000
                    model.add_reactions([nadh_test])
                    model.objective = nadh_test.id
                    nadh_sol = model.optimize()
                    model.remove_reactions([nadh_test])
                    
                    if nadh_sol.status == 'optimal' and nadh_sol.objective_value > 1e-6:
                        print("    [OK] NADH 생산 가능")
                    else:
                        print("    [FAIL] NADH 생산 불가능 -> TCA 경로 차단")
                except:
                    pass
                
                # h_p 생산 가능 여부
                try:
                    h_p = model.metabolites.get_by_id('h_p')
                    h_p_test = cobra.Reaction('TEST_h_p')
                    h_p_test.add_metabolites({h_p: -1})
                    h_p_test.lower_bound = 0
                    h_p_test.upper_bound = 1000
                    model.add_reactions([h_p_test])
                    model.objective = h_p_test.id
                    h_p_sol = model.optimize()
                    model.remove_reactions([h_p_test])
                    
                    if h_p_sol.status == 'optimal' and h_p_sol.objective_value > 1e-6:
                        print("    [OK] h_p 생산 가능")
                    else:
                        print("    [FAIL] h_p 생산 불가능 -> ETC 차단")
                except:
                    pass
                
                return False, solution
        else:
            print(f"\n  [FAIL] ATP 생산 불가능 ({solution.status})")
            return False, solution
    
    except KeyError:
        print("\n  [ERROR] atp_c metabolite 없음")
        return False, None

def test_complete_glucose_pathway(model, biomass_rxn):
    """포도당 완전 산화를 통한 Biomass 생산 테스트"""
    print("\n" + "="*70)
    print("포도당 완전 산화 → ATP → Biomass 경로 테스트")
    print("="*70)
    
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
    
    # Biomass 최적화
    model.objective = biomass_rxn.id
    
    print("\n최적화 수행 중...")
    solution = model.optimize()
    
    print(f"\n최적화 상태: {solution.status}")
    
    if solution.status == 'optimal':
        biomass_flux = solution.fluxes.get(biomass_rxn.id, 0)
        
        if biomass_flux > 1e-6:
            print(f"[SUCCESS] 포도당으로 생장 가능!")
            print(f"  Biomass flux: {biomass_flux:.6f} 1/h")
            
            # 전체 경로 플럭스
            print("\n전체 경로 플럭스:")
            
            # 포도당 섭취
            glucose_uptake = abs(solution.fluxes.get('EX_glc__D_e', 0))
            if glucose_uptake > 0:
                print(f"\n  Glucose 섭취: {glucose_uptake:.6f} mmol/gDCW/h")
            
            # ATP 생산
            try:
                atp_producing_rxns = []
                atp_c = model.metabolites.get_by_id('atp_c')
                for rxn in atp_c.reactions:
                    if atp_c in rxn.products:
                        flux = solution.fluxes.get(rxn.id, 0)
                        if abs(flux) > 1e-8:
                            atp_producing_rxns.append((rxn.id, flux))
                
                if atp_producing_rxns:
                    print(f"\n  ATP 생산 반응:")
                    for rxn_id, flux in sorted(atp_producing_rxns, key=lambda x: abs(x[1]), reverse=True)[:5]:
                        print(f"    {rxn_id}: {flux:.6f}")
            except:
                pass
            
            # NADH 생산
            try:
                nadh_c = model.metabolites.get_by_id('nadh_c')
                nadh_producing_rxns = []
                for rxn in nadh_c.reactions:
                    if nadh_c in rxn.products:
                        flux = solution.fluxes.get(rxn.id, 0)
                        if abs(flux) > 1e-8:
                            nadh_producing_rxns.append((rxn.id, flux))
                
                if nadh_producing_rxns:
                    print(f"\n  NADH 생산 반응:")
                    for rxn_id, flux in sorted(nadh_producing_rxns, key=lambda x: abs(x[1]), reverse=True)[:5]:
                        print(f"    {rxn_id}: {flux:.6f}")
            except:
                pass
            
            # ETC 경로
            print(f"\n  ETC 경로:")
            etc_rxns = ['NADH16', 'SUCD', 'QCR', 'CYO3', 'ATPS4rpp']
            etc_active = False
            for rxn_id in etc_rxns:
                try:
                    flux = solution.fluxes.get(rxn_id, 0)
                    if abs(flux) > 1e-8:
                        print(f"    {rxn_id}: {flux:.6f}")
                        etc_active = True
                except KeyError:
                    pass
            
            if etc_active:
                print("    [OK] ETC 경로 활성화!")
            
            return True, solution
        else:
            print(f"[FAIL] Biomass flux: {biomass_flux:.6f} 1/h")
            return False, solution
    
    elif solution.status == 'infeasible':
        print(f"[FAIL] 최적화 실패: infeasible")
        
        # ATP 생산은 되는지 확인
        can_produce_atp, atp_solution = test_glucose_to_nadh_pathway(model)
        
        if can_produce_atp:
            print("\n  [발견] ATP 생산은 가능하지만 Biomass 생산은 불가능")
            print("    → Biomass 구성 요소 생산 경로 문제")
        
        return False, solution
    
    else:
        print(f"[FAIL] 최적화 실패: {solution.status}")
        return False, solution

def main():
    print("="*70)
    print("포도당 완전 산화 경로 테스트")
    print("="*70)
    
    model = load_model("BaseModel.xml")
    
    # Biomass 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass: {biomass_rxn.id}")
    
    # 1. 포도당 → NADH → ATP 경로 테스트
    can_produce_atp, atp_solution = test_glucose_to_nadh_pathway(model)
    
    # 2. 포도당 완전 산화 → Biomass 테스트
    can_grow, biomass_solution = test_complete_glucose_pathway(model, biomass_rxn)
    
    # 최종 요약
    print("\n" + "="*70)
    print("최종 요약")
    print("="*70)
    
    if can_produce_atp:
        print("\n[SUCCESS] 포도당 → NADH → ATP 경로 작동 가능!")
        print("  → 포도당 완전 산화를 통해 ATP 생산 가능")
    else:
        print("\n[FAIL] 포도당 → NADH → ATP 경로 작동 불가능")
        print("  → 경로 차단 지점 확인 필요")
    
    if can_grow:
        print("\n[SUCCESS] 포도당으로 생장 가능!")
        print("  → 포도당 완전 산화 경로가 정상 작동")
    else:
        print("\n[FAIL] 포도당으로 생장 불가능")
        if can_produce_atp:
            print("  → ATP 생산은 가능하지만 Biomass 구성 요소 생산 경로 문제")
        else:
            print("  → ATP 생산 경로 문제")
    
    print("\n" + "="*70)

if __name__ == "__main__":
    main()
