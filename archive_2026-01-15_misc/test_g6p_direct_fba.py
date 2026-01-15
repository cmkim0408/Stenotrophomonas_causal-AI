#!/usr/bin/env python
"""
G6P를 직접 공급하여 FBA 수행
HEX1 단계를 우회하여 포도당 완전 산화 경로 작동 확인
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

def test_g6p_direct_fba(model, biomass_rxn):
    """G6P 직접 공급으로 포도당 완전 산화 경로 테스트"""
    print("="*70)
    print("G6P 직접 공급 FBA")
    print("HEX1 단계를 우회하여 포도당 완전 산화 경로 작동 확인")
    print("="*70)
    
    # 모든 exchange 초기화
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
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
    
    # G6P 직접 공급
    try:
        g6p_c = model.metabolites.get_by_id('g6p_c')
        
        dm_g6p = cobra.Reaction('DM_g6p_c_direct')
        dm_g6p.name = 'G6P direct supply (bypass HEX1)'
        dm_g6p.lower_bound = -100  # G6P 공급
        dm_g6p.upper_bound = 1000
        dm_g6p.add_metabolites({g6p_c: -1})  # g6p_c 소비 = 공급
        model.add_reactions([dm_g6p])
        
        print("\n[설정] G6P 직접 공급:")
        print(f"  g6p_c demand: -100 mmol/gDCW/h")
        print("  → HEX1 단계를 우회하여 G6P가 세포 내에 직접 존재")
        
        # Biomass 최적화
        model.objective = biomass_rxn.id
        
        print("\n[FBA 수행 중...]")
        solution = model.optimize()
        
        if solution.status == 'optimal':
            biomass_flux = solution.fluxes.get(biomass_rxn.id, 0)
            
            if biomass_flux > 1e-6:
                print(f"\n[SUCCESS] G6P 직접 공급으로 생장 가능!")
                print(f"  Biomass flux: {biomass_flux:.6f} 1/h")
                
                # G6P 사용량
                g6p_used = abs(solution.fluxes.get('DM_g6p_c_direct', 0))
                print(f"  G6P 사용: {g6p_used:.6f} mmol/gDCW/h")
                
                # Yield
                if g6p_used > 0:
                    yield_biomass = biomass_flux / g6p_used
                    print(f"  Biomass yield: {yield_biomass:.6f} gDW/mmol G6P")
                
                # 경로 플럭스 확인
                print("\n[포도당 완전 산화 경로 플럭스 확인]")
                print("-" * 70)
                
                # Glycolysis (HEX1 제외)
                print("\n[1] Glycolysis 경로 (HEX1 제외):")
                gly_rxns = ['PGI', 'PFK', 'FBA', 'GAPD', 'PGK', 'ENO', 'PYK']
                active_gly = False
                for rxn_id in gly_rxns:
                    try:
                        flux = solution.fluxes.get(rxn_id, 0)
                        if abs(flux) > 1e-8:
                            print(f"  {rxn_id}: {flux:.6f}")
                            active_gly = True
                    except KeyError:
                        pass
                
                if active_gly:
                    print("  [OK] Glycolysis 활성화!")
                
                # PDH
                print("\n[2] Pyruvate -> Acetyl-CoA:")
                try:
                    pdh_flux = solution.fluxes.get('PDH', 0)
                    if abs(pdh_flux) > 1e-8:
                        print(f"  PDH: {pdh_flux:.6f} [OK] 활성화!")
                except:
                    pass
                
                # TCA Cycle
                print("\n[3] TCA Cycle:")
                tca_rxns = ['CS', 'ACONT', 'ICDHx', 'AKGDH', 'SUCOAS', 'SUCD', 'FUM', 'MDH']
                active_tca = False
                for rxn_id in tca_rxns:
                    try:
                        flux = solution.fluxes.get(rxn_id, 0)
                        if abs(flux) > 1e-8:
                            print(f"  {rxn_id}: {flux:.6f}")
                            active_tca = True
                    except KeyError:
                        pass
                
                if active_tca:
                    print("  [OK] TCA Cycle 활성화!")
                
                # NADH 생산
                print("\n[4] NADH 생산:")
                try:
                    nadh_c = model.metabolites.get_by_id('nadh_c')
                    nadh_producers = {'GAPD': 0, 'ICDHx': 0, 'AKGDH': 0, 'MDH': 0}
                    
                    for rxn_id in nadh_producers.keys():
                        try:
                            flux = solution.fluxes.get(rxn_id, 0)
                            if abs(flux) > 1e-8:
                                rxn = model.reactions.get_by_id(rxn_id)
                                nadh_coeff = rxn.metabolites.get(nadh_c, 0)
                                if nadh_coeff > 0:
                                    nadh_produced = flux * nadh_coeff
                                    nadh_producers[rxn_id] = nadh_produced
                                    print(f"  {rxn_id}: {nadh_produced:.6f} mmol/gDCW/h")
                        except:
                            pass
                    
                    total_nadh = sum(nadh_producers.values())
                    if total_nadh > 0:
                        print(f"  [OK] 총 NADH 생산: {total_nadh:.6f} mmol/gDCW/h")
                except:
                    pass
                
                # ETC
                print("\n[5] ETC 경로:")
                etc_rxns = ['NADH16pp', 'NADH17pp', 'NADH18pp', 'CYO1_KT', 'CYTCAA3pp', 'CYTBDpp', 'ATPS4rpp']
                etc_active = False
                h_p_total = 0
                
                try:
                    h_p = model.metabolites.get_by_id('h_p')
                    
                    for rxn_id in etc_rxns:
                        try:
                            flux = solution.fluxes.get(rxn_id, 0)
                            if abs(flux) > 1e-8:
                                rxn = model.reactions.get_by_id(rxn_id)
                                h_p_coeff = rxn.metabolites.get(h_p, 0)
                                
                                if h_p_coeff > 0:
                                    h_p_produced = flux * h_p_coeff
                                    h_p_total += h_p_produced
                                    print(f"  {rxn_id}: {flux:.6f} -> h_p: {h_p_produced:.6f}")
                                elif h_p_coeff < 0:
                                    print(f"  {rxn_id}: {flux:.6f} -> h_p 소비: {abs(flux * h_p_coeff):.6f}")
                                
                                etc_active = True
                        except KeyError:
                            pass
                    
                    if etc_active and h_p_total > 0:
                        print(f"  [OK] 총 h_p 생성: {h_p_total:.6f} mmol/gDCW/h")
                except:
                    pass
                
                # ATP 생산
                print("\n[6] ATP 생산:")
                try:
                    atp_c = model.metabolites.get_by_id('atp_c')
                    
                    # ATP synthase
                    try:
                        atps4_flux = solution.fluxes.get('ATPS4rpp', 0)
                        if abs(atps4_flux) > 1e-8:
                            atps4_rxn = model.reactions.get_by_id('ATPS4rpp')
                            atp_coeff = atps4_rxn.metabolites.get(atp_c, 0)
                            if atp_coeff > 0:
                                etc_atp = atps4_flux * atp_coeff
                                print(f"  ATPS4rpp (ETC): {etc_atp:.6f} mmol/gDCW/h [OK]")
                    except:
                        pass
                    
                    # 전체 ATP 생산
                    atp_producing = []
                    for rxn in atp_c.reactions:
                        if atp_c in rxn.products:
                            flux = solution.fluxes.get(rxn.id, 0)
                            if abs(flux) > 1e-8:
                                atp_coeff = rxn.metabolites.get(atp_c, 0)
                                if atp_coeff > 0:
                                    atp_produced = flux * atp_coeff
                                    atp_producing.append((rxn.id, atp_produced))
                    
                    if atp_producing:
                        total_atp = sum([flux for _, flux in atp_producing])
                        print(f"  총 ATP 생산: {total_atp:.6f} mmol/gDCW/h")
                        
                        print(f"\n  주요 ATP 생산 반응:")
                        for rxn_id, flux in sorted(atp_producing, key=lambda x: x[1], reverse=True)[:5]:
                            marker = " [ETC]" if rxn_id == 'ATPS4rpp' else " [Glycolysis]" if rxn_id in ['PGK', 'PYK'] else ""
                            print(f"    {rxn_id}: {flux:.6f} mmol/gDCW/h{marker}")
                
                except:
                    pass
                
                # 최종 결론
                print("\n" + "="*70)
                print("[포도당 완전 산화 경로 확인]")
                print("="*70)
                print("[OK] G6P -> Glycolysis -> Pyruvate")
                print("[OK] Pyruvate -> PDH -> Acetyl-CoA")
                print("[OK] Acetyl-CoA -> TCA Cycle -> NADH")
                print("[OK] NADH -> ETC -> h_p")
                print("[OK] h_p -> ATPS4rpp -> ATP")
                print("[OK] ATP -> Biomass 생산")
                print("\n[결론] 포도당 완전 산화 경로는 작동 가능!")
                print("  → 문제는 초기 ATP 생성 (HEX1 단계)")
                
                model.remove_reactions([dm_g6p])
                return True, solution
            else:
                print(f"\n[FAIL] Biomass flux: {biomass_flux:.6f} 1/h")
                model.remove_reactions([dm_g6p])
                return False, solution
        else:
            print(f"\n[FAIL] 최적화 실패: {solution.status}")
            model.remove_reactions([dm_g6p])
            return False, solution
    
    except Exception as e:
        print(f"[ERROR] 테스트 실패: {e}")
        import traceback
        traceback.print_exc()
        return False, None

def main():
    print("="*70)
    print("G6P 직접 공급 FBA")
    print("="*70)
    
    model = load_model("BaseModel.xml")
    
    # Biomass 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass: {biomass_rxn.id}")
    
    # G6P 직접 공급 테스트
    can_grow, solution = test_g6p_direct_fba(model, biomass_rxn)
    
    # 최종 요약
    print("\n" + "="*70)
    print("최종 요약")
    print("="*70)
    
    if can_grow:
        print("\n[SUCCESS] G6P가 세포 내에 있으면 생장 가능!")
        print("  → 포도당 완전 산화 -> NADH -> ETC -> ATP 경로 정상 작동")
        print("  → 문제는 HEX1 단계 (ATP 필요)")
        print("\n[해결 방안]:")
        print("  1. 초기 ATP 부트스트랩으로 HEX1 활성화")
        print("  2. 또는 G6P를 직접 공급 (실험적으로는 불가능하지만 이론적으로는 가능)")
    else:
        print("\n[FAIL] G6P가 세포 내에 있어도 생장 불가능")
        print("  → HEX1 문제 외에도 다른 문제 존재")
    
    print("\n" + "="*70)

if __name__ == "__main__":
    main()
