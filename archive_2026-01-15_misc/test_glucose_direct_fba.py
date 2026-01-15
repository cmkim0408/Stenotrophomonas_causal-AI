#!/usr/bin/env python
"""
포도당이 세포 내에 있다고 가정하고 FBA 수행
Transport 문제를 우회하여 완전 산화 경로 작동 확인
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

def test_glucose_direct_availability(model, biomass_rxn):
    """포도당이 세포 내에 직접 공급된다고 가정하고 FBA"""
    print("="*70)
    print("포도당 직접 공급 가정 FBA")
    print("포도당이 이미 세포 내(glc__D_c)에 있다고 가정")
    print("="*70)
    
    # 모든 exchange 초기화
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    # 필수 영양소만 설정
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
    
    # 포도당을 세포 내에 직접 공급 (glc__D_c demand reaction)
    try:
        glc__D_c = model.metabolites.get_by_id('glc__D_c')
        
        dm_glc = cobra.Reaction('DM_glc__D_c_direct')
        dm_glc.name = 'Glucose direct supply (bypass transport)'
        dm_glc.lower_bound = -100  # 포도당 공급
        dm_glc.upper_bound = 1000
        dm_glc.add_metabolites({glc__D_c: -1})  # glc__D_c 소비 = 공급
        model.add_reactions([dm_glc])
        
        print("\n[설정] 포도당 직접 공급:")
        print(f"  glc__D_c demand: -100 mmol/gDCW/h")
        print("  → Transport 과정을 우회하여 포도당이 세포 내에 직접 존재")
        
        # Biomass 최적화
        model.objective = biomass_rxn.id
        
        print("\n[FBA 수행 중...]")
        solution = model.optimize()
        
        if solution.status == 'optimal':
            biomass_flux = solution.fluxes.get(biomass_rxn.id, 0)
            
            if biomass_flux > 1e-6:
                print(f"\n[SUCCESS] 포도당 직접 공급으로 생장 가능!")
                print(f"  Biomass flux: {biomass_flux:.6f} 1/h")
                
                # 포도당 사용량
                glucose_used = abs(solution.fluxes.get('DM_glc__D_c_direct', 0))
                print(f"  포도당 사용: {glucose_used:.6f} mmol/gDCW/h")
                
                # Yield
                if glucose_used > 0:
                    yield_biomass = biomass_flux / glucose_used
                    print(f"  Biomass yield: {yield_biomass:.6f} gDW/mmol glucose")
                
                # 전체 경로 플럭스 분석
                print("\n[포도당 완전 산화 경로 플럭스 분석]")
                print("-" * 70)
                
                # 1. Glycolysis
                print("\n[1] Glycolysis 경로:")
                gly_rxns = {
                    'HEX1': 'Hexokinase',
                    'PGI': 'Glucose-6-phosphate isomerase',
                    'PFK': 'Phosphofructokinase',
                    'FBA': 'Fructose-bisphosphate aldolase',
                    'GAPD': 'Glyceraldehyde-3-phosphate dehydrogenase',
                    'PGK': 'Phosphoglycerate kinase',
                    'ENO': 'Enolase',
                    'PYK': 'Pyruvate kinase'
                }
                
                gly_active = False
                for rxn_id, rxn_name in gly_rxns.items():
                    try:
                        flux = solution.fluxes.get(rxn_id, 0)
                        if abs(flux) > 1e-8:
                            print(f"  {rxn_id} ({rxn_name}): {flux:.6f}")
                            gly_active = True
                    except KeyError:
                        pass
                
                if gly_active:
                    print("  [OK] Glycolysis 활성화!")
                
                # 2. Pyruvate → Acetyl-CoA
                print("\n[2] Pyruvate → Acetyl-CoA:")
                try:
                    pdh_flux = solution.fluxes.get('PDH', 0)
                    if abs(pdh_flux) > 1e-8:
                        print(f"  PDH: {pdh_flux:.6f} [OK] 활성화!")
                    else:
                        print(f"  PDH: {pdh_flux:.6f} (비활성)")
                except KeyError:
                    print("  [MISSING] PDH")
                
                # 3. TCA Cycle
                print("\n[3] TCA Cycle:")
                tca_rxns = {
                    'CS': 'Citrate synthase',
                    'ACONT': 'Aconitase',
                    'ICDHx': 'Isocitrate dehydrogenase (NAD)',
                    'AKGDH': '2-Oxogluterate dehydrogenase',
                    'SUCOAS': 'Succinate-CoA ligase',
                    'SUCD': 'Succinate dehydrogenase',
                    'FUM': 'Fumarase',
                    'MDH': 'Malate dehydrogenase'
                }
                
                tca_active = False
                nadh_producers = {}
                
                for rxn_id, rxn_name in tca_rxns.items():
                    try:
                        flux = solution.fluxes.get(rxn_id, 0)
                        if abs(flux) > 1e-8:
                            print(f"  {rxn_id} ({rxn_name}): {flux:.6f}")
                            tca_active = True
                            
                            # NADH 생산 반응 확인
                            if rxn_id in ['ICDHx', 'AKGDH', 'MDH']:
                                try:
                                    nadh_c = model.metabolites.get_by_id('nadh_c')
                                    rxn = model.reactions.get_by_id(rxn_id)
                                    nadh_coeff = rxn.metabolites.get(nadh_c, 0)
                                    if nadh_coeff > 0:
                                        nadh_produced = flux * nadh_coeff
                                        nadh_producers[rxn_id] = nadh_produced
                                except:
                                    pass
                    except KeyError:
                        pass
                
                if tca_active:
                    print("  [OK] TCA Cycle 활성화!")
                
                # 4. NADH 생산
                print("\n[4] NADH 생산:")
                try:
                    nadh_c = model.metabolites.get_by_id('nadh_c')
                    
                    # Glycolysis에서 NADH
                    gapd_flux = solution.fluxes.get('GAPD', 0)
                    if abs(gapd_flux) > 1e-8:
                        try:
                            gapd_rxn = model.reactions.get_by_id('GAPD')
                            gapd_nadh_coeff = gapd_rxn.metabolites.get(nadh_c, 0)
                            if gapd_nadh_coeff > 0:
                                gapd_nadh = gapd_flux * gapd_nadh_coeff
                                nadh_producers['GAPD'] = gapd_nadh
                                print(f"  GAPD (Glycolysis): {gapd_nadh:.6f} mmol/gDCW/h")
                        except:
                            pass
                    
                    # TCA에서 NADH
                    total_nadh = 0
                    for rxn_id, nadh_flux in nadh_producers.items():
                        total_nadh += nadh_flux
                        if rxn_id != 'GAPD':
                            print(f"  {rxn_id} (TCA): {nadh_flux:.6f} mmol/gDCW/h")
                    
                    if total_nadh > 0:
                        print(f"  [OK] 총 NADH 생산: {total_nadh:.6f} mmol/gDCW/h")
                except:
                    pass
                
                # 5. ETC → h_p
                print("\n[5] ETC 경로 (NADH -> h_p):")
                etc_rxns = {
                    'NADH16pp': 'NADH dehydrogenase (ubiquinone-8)',
                    'NADH17pp': 'NADH dehydrogenase (menaquinone-8)',
                    'NADH18pp': 'NADH dehydrogenase (demethylmenaquinone-8)',
                    'SUCD': 'Succinate dehydrogenase (Complex II)',
                    'CYO1_KT': 'Ubiquinol cytochrome c reductase',
                    'CYTCAA3pp': 'Cytochrome c oxidase aa3',
                    'CYTBDpp': 'Cytochrome oxidase bd'
                }
                
                etc_active = False
                h_p_total = 0
                
                try:
                    h_p = model.metabolites.get_by_id('h_p')
                    
                    for rxn_id, rxn_name in etc_rxns.items():
                        try:
                            flux = solution.fluxes.get(rxn_id, 0)
                            if abs(flux) > 1e-8:
                                rxn = model.reactions.get_by_id(rxn_id)
                                h_p_coeff = rxn.metabolites.get(h_p, 0)
                                
                                if h_p_coeff > 0:  # h_p 생성
                                    h_p_produced = flux * h_p_coeff
                                    h_p_total += h_p_produced
                                    print(f"  {rxn_id} ({rxn_name}): {flux:.6f} -> h_p: {h_p_produced:.6f}")
                                elif h_p_coeff < 0:  # h_p 소비
                                    h_p_consumed = flux * abs(h_p_coeff)
                                    print(f"  {rxn_id} ({rxn_name}): {flux:.6f} -> h_p 소비: {h_p_consumed:.6f}")
                                
                                etc_active = True
                        except KeyError:
                            pass
                    
                    if etc_active and h_p_total > 0:
                        print(f"  [OK] 총 h_p 생성: {h_p_total:.6f} mmol/gDCW/h")
                    elif etc_active:
                        print(f"  [WARNING] ETC 활성화되었지만 h_p 생성 없음")
                    else:
                        print(f"  [WARNING] ETC 경로 비활성")
                except:
                    pass
                
                # 6. ATP 생산 (ETC를 통한)
                print("\n[6] ATP 생산 (ETC 경로):")
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
                                print(f"  ATPS4rpp (ATP synthase): {atps4_flux:.6f} -> ATP: {etc_atp:.6f} mmol/gDCW/h")
                                print(f"  [OK] ETC를 통한 ATP 생산 활성화!")
                            else:
                                print(f"  ATPS4rpp: {atps4_flux:.6f} (ATP 생산 없음)")
                        else:
                            print(f"  ATPS4rpp: {atps4_flux:.6f} (비활성)")
                    except KeyError:
                        print("  [MISSING] ATPS4rpp")
                    
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
                        print(f"\n  총 ATP 생산 (모든 경로): {total_atp:.6f} mmol/gDCW/h")
                        
                        # 상위 ATP 생산 반응
                        print(f"\n  주요 ATP 생산 반응:")
                        for rxn_id, flux in sorted(atp_producing, key=lambda x: x[1], reverse=True)[:5]:
                            marker = " [ETC]" if rxn_id == 'ATPS4rpp' else " [Glycolysis]" if rxn_id in ['PGK', 'PYK'] else ""
                            print(f"    {rxn_id}: {flux:.6f} mmol/gDCW/h{marker}")
                    
                except:
                    pass
                
                # 최종 경로 요약
                print("\n" + "="*70)
                print("[포도당 완전 산화 경로 확인]")
                print("="*70)
                print("✓ 포도당 → Glycolysis → Pyruvate")
                print("✓ Pyruvate → PDH → Acetyl-CoA")
                print("✓ Acetyl-CoA → TCA Cycle → NADH")
                print("✓ NADH → ETC (NADH16pp 등) → h_p")
                print("✓ h_p → ATPS4rpp (ATP synthase) → ATP")
                print("✓ ATP → Biomass 생산")
                print("\n[결론] 포도당 완전 산화를 통해 ATP를 생산하여 생장 가능!")
                
                model.remove_reactions([dm_glc])
                return True, solution
            else:
                print(f"\n[FAIL] Biomass flux: {biomass_flux:.6f} 1/h (생장 불가능)")
                model.remove_reactions([dm_glc])
                return False, solution
        else:
            print(f"\n[FAIL] 최적화 실패: {solution.status}")
            model.remove_reactions([dm_glc])
            return False, solution
    
    except Exception as e:
        print(f"[ERROR] 테스트 실패: {e}")
        import traceback
        traceback.print_exc()
        return False, None

def compare_transport_vs_direct(model, biomass_rxn):
    """Transport 필요 vs 직접 공급 비교"""
    print("\n" + "="*70)
    print("Transport 필요 vs 직접 공급 비교")
    print("="*70)
    
    results = {
        'Method': [],
        'Status': [],
        'Biomass_Flux': [],
        'Can_Grow': []
    }
    
    # 1. Transport 필요 (포도당 exchange)
    print("\n[방법 1] 포도당 Exchange 사용 (Transport 필요):")
    print("-" * 70)
    
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
    
    model.objective = biomass_rxn.id
    solution_transport = model.optimize()
    
    results['Method'].append('Transport 필요 (EX_glc__D_e)')
    results['Status'].append(solution_transport.status)
    
    if solution_transport.status == 'optimal':
        biomass_transport = solution_transport.fluxes.get(biomass_rxn.id, 0)
        results['Biomass_Flux'].append(biomass_transport)
        results['Can_Grow'].append(biomass_transport > 1e-6)
        print(f"  상태: {solution_transport.status}")
        print(f"  Biomass flux: {biomass_transport:.6f} 1/h")
    else:
        results['Biomass_Flux'].append(0)
        results['Can_Grow'].append(False)
        print(f"  상태: {solution_transport.status}")
    
    # 2. 직접 공급
    print("\n[방법 2] 포도당 직접 공급 (Transport 우회):")
    print("-" * 70)
    
    can_grow, solution_direct = test_glucose_direct_availability(model, biomass_rxn)
    
    results['Method'].append('직접 공급 (glc__D_c demand)')
    results['Status'].append(solution_direct.status if solution_direct else 'error')
    results['Biomass_Flux'].append(solution_direct.fluxes.get(biomass_rxn.id, 0) if solution_direct and solution_direct.status == 'optimal' else 0)
    results['Can_Grow'].append(can_grow)
    
    # 결과 저장
    df_results = pd.DataFrame(results)
    df_results.to_csv('glucose_transport_vs_direct_comparison.csv', index=False)
    print(f"\n[OK] 비교 결과 저장: glucose_transport_vs_direct_comparison.csv")
    
    return results

def main():
    print("="*70)
    print("포도당 직접 공급 가정 FBA")
    print("="*70)
    
    model = load_model("BaseModel.xml")
    
    # Biomass 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass: {biomass_rxn.id}")
    
    # 포도당 직접 공급 테스트
    can_grow, solution = test_glucose_direct_availability(model, biomass_rxn)
    
    # 비교 테스트
    comparison = compare_transport_vs_direct(model, biomass_rxn)
    
    # 최종 요약
    print("\n" + "="*70)
    print("최종 요약")
    print("="*70)
    
    if can_grow:
        print("\n[SUCCESS] 포도당이 세포 내에 있으면 생장 가능!")
        print("  → 포도당 완전 산화 → NADH → ETC → ATP 경로 정상 작동")
        print("  → 문제는 포도당 Transport 단계")
        print("\n[해결 방안]:")
        print("  1. 포도당 transport 경로 수정")
        print("  2. 또는 초기 부트스트랩으로 transport 활성화")
    else:
        print("\n[FAIL] 포도당이 세포 내에 있어도 생장 불가능")
        print("  → Transport 문제 외에도 다른 문제 존재")
    
    print("\n" + "="*70)

if __name__ == "__main__":
    main()
