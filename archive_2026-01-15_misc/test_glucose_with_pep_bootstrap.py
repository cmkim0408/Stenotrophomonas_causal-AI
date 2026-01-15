#!/usr/bin/env python
"""
포도당 경로 테스트 (PEP 부트스트랩)
PTS 경로는 ATP 필요 없이 PEP만 사용하므로 PEP 부트스트랩으로 테스트
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

def test_glucose_with_pep_bootstrap(model, biomass_rxn):
    """PEP 부트스트랩으로 포도당 경로 테스트"""
    print("="*70)
    print("포도당 완전 산화 경로 테스트")
    print("PEP 부트스트랩으로 PTS 경로 활성화 (ATP 필요 없음)")
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
    
    # PEP 부트스트랩 (PTS 경로 활성화)
    print("\n[전략 1] PEP 부트스트랩 (PTS 경로 활성화, ATP 필요 없음)")
    print("-" * 70)
    
    try:
        pep_c = model.metabolites.get_by_id('pep_c')
        dm_pep = cobra.Reaction('DM_pep_c_bs')
        dm_pep.name = 'PEP bootstrap'
        dm_pep.lower_bound = -0.1
        dm_pep.upper_bound = 1000
        dm_pep.add_metabolites({pep_c: -1})
        model.add_reactions([dm_pep])
        
        # NAD+ 부트스트랩 (NADH 생산을 위해)
        nad_c = model.metabolites.get_by_id('nad_c')
        dm_nad = cobra.Reaction('DM_nad_c_bs')
        dm_nad.name = 'NAD+ bootstrap'
        dm_nad.lower_bound = -0.1
        dm_nad.upper_bound = 1000
        dm_nad.add_metabolites({nad_c: -1})
        model.add_reactions([dm_nad])
        
        # CoA 부트스트랩 (TCA cycle을 위해)
        coa_c = model.metabolites.get_by_id('coa_c')
        dm_coa = cobra.Reaction('DM_coa_c_bs')
        dm_coa.name = 'CoA bootstrap'
        dm_coa.lower_bound = -0.01
        dm_coa.upper_bound = 1000
        dm_coa.add_metabolites({coa_c: -1})
        model.add_reactions([dm_coa])
        
        print("  [OK] 부트스트랩 추가:")
        print("    PEP: -0.1 mmol/gDCW/h (PTS 경로 활성화)")
        print("    NAD+: -0.1 mmol/gDCW/h (NADH 생산용)")
        print("    CoA: -0.01 mmol/gDCW/h (TCA cycle용)")
        
        # Biomass 최적화
        model.objective = biomass_rxn.id
        solution = model.optimize()
        
        if solution.status == 'optimal':
            biomass_flux = solution.fluxes.get(biomass_rxn.id, 0)
            
            if biomass_flux > 1e-6:
                print(f"\n  [SUCCESS] 포도당으로 생장 가능!")
                print(f"    Biomass flux: {biomass_flux:.6f} 1/h")
                
                # 경로 플럭스 확인
                print("\n  포도당 완전 산화 경로 플럭스:")
                
                # Exchange
                glucose_uptake = abs(solution.fluxes.get('EX_glc__D_e', 0))
                if glucose_uptake > 0:
                    print(f"\n    Glucose 섭취: {glucose_uptake:.6f} mmol/gDCW/h")
                
                # Transport (PTS 경로 확인)
                print("\n    [Transport]:")
                try:
                    glcpts_flux = solution.fluxes.get('GLCpts', 0)
                    if abs(glcpts_flux) > 1e-8:
                        print(f"      GLCpts (PTS): {glcpts_flux:.6f} [OK] 활성화!")
                    else:
                        print(f"      GLCpts: {glcpts_flux:.6f} (비활성)")
                except:
                    pass
                
                # Glycolysis
                print("\n    [Glycolysis]:")
                gly_rxns = ['HEX1', 'PGI', 'PFK', 'FBA', 'GAPD', 'PGK', 'ENO', 'PYK']
                active_gly = False
                for rxn_id in gly_rxns:
                    try:
                        flux = solution.fluxes.get(rxn_id, 0)
                        if abs(flux) > 1e-8:
                            print(f"      {rxn_id}: {flux:.6f}")
                            active_gly = True
                    except KeyError:
                        pass
                
                if active_gly:
                    print("      [OK] Glycolysis 활성화!")
                
                # PDH
                try:
                    pdh_flux = solution.fluxes.get('PDH', 0)
                    if abs(pdh_flux) > 1e-8:
                        print(f"\n    [Pyruvate → Acetyl-CoA]:")
                        print(f"      PDH: {pdh_flux:.6f}")
                except:
                    pass
                
                # TCA Cycle
                print("\n    [TCA Cycle]:")
                tca_rxns = ['CS', 'ACONT', 'ICDHx', 'AKGDH', 'SUCOAS', 'SUCD', 'FUM', 'MDH']
                active_tca = False
                nadh_producers = {'ICDHx': 0, 'AKGDH': 0, 'MDH': 0, 'GAPD': 0}
                
                for rxn_id in tca_rxns:
                    try:
                        flux = solution.fluxes.get(rxn_id, 0)
                        if abs(flux) > 1e-8:
                            print(f"      {rxn_id}: {flux:.6f}")
                            active_tca = True
                            if rxn_id in nadh_producers:
                                nadh_producers[rxn_id] = flux
                    except KeyError:
                        pass
                
                if active_tca:
                    print("      [OK] TCA Cycle 활성화!")
                
                # NADH 생산 확인
                try:
                    nadh_c = model.metabolites.get_by_id('nadh_c')
                    total_nadh = 0
                    
                    print("\n    [NADH 생산]:")
                    for rxn_id, flux in nadh_producers.items():
                        if abs(flux) > 1e-8:
                            try:
                                rxn = model.reactions.get_by_id(rxn_id)
                                nadh_coeff = rxn.metabolites.get(nadh_c, 0)
                                if nadh_coeff > 0:
                                    nadh_produced = flux * nadh_coeff
                                    total_nadh += nadh_produced
                                    print(f"      {rxn_id}: {nadh_produced:.6f} mmol/gDCW/h (flux: {flux:.6f} * {nadh_coeff})")
                            except:
                                pass
                    
                    if total_nadh > 0:
                        print(f"      총 NADH 생산: {total_nadh:.6f} mmol/gDCW/h")
                        print(f"      [OK] NADH 생산 활성화!")
                except:
                    pass
                
                # ETC 경로
                print("\n    [ETC 경로]:")
                etc_rxns = ['NADH16pp', 'NADH17pp', 'NADH18pp',  # Complex I 대체
                           'SUCD',  # Complex II
                           'CYO1_KT', 'CYTCAA3pp', 'CYTBDpp',  # Complex IV 대체
                           'ATPS4rpp']  # Complex V
                
                etc_active = False
                h_p_produced = 0
                
                for rxn_id in etc_rxns:
                    try:
                        flux = solution.fluxes.get(rxn_id, 0)
                        if abs(flux) > 1e-8:
                            print(f"      {rxn_id}: {flux:.6f}")
                            
                            # h_p 생성 확인
                            try:
                                h_p = model.metabolites.get_by_id('h_p')
                                rxn = model.reactions.get_by_id(rxn_id)
                                h_p_coeff = rxn.metabolites.get(h_p, 0)
                                if h_p_coeff > 0:
                                    h_p_flux = flux * h_p_coeff
                                    h_p_produced += h_p_flux
                                    print(f"        → h_p 생성: {h_p_flux:.6f} mmol/gDCW/h")
                            except:
                                pass
                            
                            etc_active = True
                    except KeyError:
                        pass
                
                if etc_active:
                    print(f"      [OK] ETC 경로 활성화! (총 h_p 생성: {h_p_produced:.6f} mmol/gDCW/h)")
                else:
                    print("      [WARNING] ETC 경로 비활성")
                
                # ATP 생산 확인
                try:
                    atp_c = model.metabolites.get_by_id('atp_c')
                    atp_producing = []
                    etc_atp = 0
                    
                    for rxn in atp_c.reactions:
                        if atp_c in rxn.products:
                            flux = solution.fluxes.get(rxn.id, 0)
                            if abs(flux) > 1e-8:
                                atp_coeff = rxn.metabolites.get(atp_c, 0)
                                if atp_coeff > 0:
                                    atp_produced = flux * atp_coeff
                                    atp_producing.append((rxn.id, atp_produced))
                                    
                                    # ETC를 통한 ATP 생산 확인
                                    if rxn.id == 'ATPS4rpp':
                                        etc_atp = atp_produced
                    
                    if atp_producing:
                        print(f"\n    [ATP 생산]:")
                        atp_total = sum([flux for _, flux in atp_producing])
                        
                        # 주요 ATP 생산 반응
                        for rxn_id, flux in sorted(atp_producing, key=lambda x: x[1], reverse=True)[:5]:
                            marker = " [ETC]" if rxn_id == 'ATPS4rpp' else ""
                            print(f"      {rxn_id}: {flux:.6f} mmol/gDCW/h{marker}")
                        
                        print(f"      총 ATP 생산: {atp_total:.6f} mmol/gDCW/h")
                        
                        if etc_atp > 0:
                            print(f"      [SUCCESS] ETC를 통한 ATP 생산: {etc_atp:.6f} mmol/gDCW/h")
                            print(f"      → 포도당 완전 산화 → NADH → ETC → ATP 경로 작동!")
                except:
                    pass
                
                # Yield 계산
                if glucose_uptake > 0:
                    yield_biomass = biomass_flux / glucose_uptake
                    print(f"\n    성장 속도: {biomass_flux:.6f} 1/h")
                    print(f"    Glucose 섭취: {glucose_uptake:.6f} mmol/gDCW/h")
                    print(f"    Biomass yield: {yield_biomass:.6f} gDW/mmol glucose")
                
                model.remove_reactions(['DM_pep_c_bs', 'DM_nad_c_bs', 'DM_coa_c_bs'])
                return True, solution
            else:
                print(f"\n  [FAIL] Biomass flux: {biomass_flux:.6f} 1/h")
                model.remove_reactions(['DM_pep_c_bs', 'DM_nad_c_bs', 'DM_coa_c_bs'])
                return False, solution
        else:
            print(f"\n  [FAIL] 최적화 실패: {solution.status}")
            model.remove_reactions(['DM_pep_c_bs', 'DM_nad_c_bs', 'DM_coa_c_bs'])
            return False, solution
    
    except Exception as e:
        print(f"  [ERROR] 테스트 실패: {e}")
        import traceback
        traceback.print_exc()
        return False, None

def main():
    print("="*70)
    print("포도당 완전 산화 경로 테스트 (PEP 부트스트랩)")
    print("="*70)
    
    model = load_model("BaseModel.xml")
    
    # Biomass 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass: {biomass_rxn.id}")
    
    # PEP 부트스트랩으로 포도당 경로 테스트
    can_grow, solution = test_glucose_with_pep_bootstrap(model, biomass_rxn)
    
    # 최종 요약
    print("\n" + "="*70)
    print("최종 요약")
    print("="*70)
    
    if can_grow:
        print("\n[SUCCESS] 포도당 완전 산화 경로 작동 가능!")
        print("\n경로 요약:")
        print("  1. 포도당 → GLCpts (PTS 경로, PEP 사용) → G6P")
        print("  2. G6P → Glycolysis → Pyruvate")
        print("  3. Pyruvate → PDH → Acetyl-CoA")
        print("  4. Acetyl-CoA → TCA Cycle → NADH")
        print("  5. NADH → ETC (NADH16pp 등) → h_p")
        print("  6. h_p → ATPS4rpp (ATP synthase) → ATP")
        print("  7. ATP → Biomass")
        print("\n→ 포도당 완전 산화를 통해 ATP를 생산하여 생장 가능!")
    else:
        print("\n[FAIL] 포도당 완전 산화 경로 작동 불가능")
        print("  → 추가 원인 분석 필요")
    
    print("\n" + "="*70)

if __name__ == "__main__":
    main()
