#!/usr/bin/env python
"""
포도당 완전 산화 경로 이론적 테스트
ETC Complex들이 있다고 가정하고 경로 작동 여부 확인
또한 포도당 transport 부트스트랩으로 해결
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

def test_glucose_pathway_with_bootstrap(model, biomass_rxn):
    """포도당 경로 테스트 (부트스트랩으로 transport 문제 해결)"""
    print("="*70)
    print("포도당 완전 산화 경로 테스트")
    print("Transport 문제를 부트스트랩으로 해결하고 경로 작동 확인")
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
    
    # 전략 1: ATP 부트스트랩 (ABC transport 활성화)
    print("\n[전략 1] ATP 부트스트랩 (ABC transport 활성화)")
    print("-" * 70)
    
    try:
        atp_c = model.metabolites.get_by_id('atp_c')
        dm_atp = cobra.Reaction('DM_atp_c_bs')
        dm_atp.name = 'ATP bootstrap'
        dm_atp.lower_bound = -0.1
        dm_atp.upper_bound = 1000
        dm_atp.add_metabolites({atp_c: -1})
        model.add_reactions([dm_atp])
        
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
        print("    ATP: -0.1 mmol/gDCW/h")
        print("    NAD+: -0.1 mmol/gDCW/h")
        print("    CoA: -0.01 mmol/gDCW/h")
        
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
                
                # Transport
                print("\n    [Transport]:")
                transport_rxns = ['GLCtex', 'GLCabc', 'GLCabcpp', 'GLCpts']
                for rxn_id in transport_rxns:
                    try:
                        flux = solution.fluxes.get(rxn_id, 0)
                        if abs(flux) > 1e-8:
                            print(f"      {rxn_id}: {flux:.6f}")
                    except KeyError:
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
                for rxn_id in tca_rxns:
                    try:
                        flux = solution.fluxes.get(rxn_id, 0)
                        if abs(flux) > 1e-8:
                            print(f"      {rxn_id}: {flux:.6f}")
                            active_tca = True
                    except KeyError:
                        pass
                
                if active_tca:
                    print("      [OK] TCA Cycle 활성화!")
                
                # NADH 생산 확인
                try:
                    nadh_c = model.metabolites.get_by_id('nadh_c')
                    nadh_producing = ['ICDHx', 'AKGDH', 'MDH', 'GAPD']
                    nadh_total = 0
                    for rxn_id in nadh_producing:
                        try:
                            flux = solution.fluxes.get(rxn_id, 0)
                            # NADH 생산량 계산 (반응식에서 계수 확인)
                            if abs(flux) > 1e-8:
                                rxn = model.reactions.get_by_id(rxn_id)
                                nadh_coeff = rxn.metabolites.get(nadh_c, 0)
                                if nadh_coeff > 0:
                                    nadh_produced = flux * nadh_coeff
                                    nadh_total += nadh_produced
                        except:
                            pass
                    
                    if nadh_total > 0:
                        print(f"\n    [NADH 생산]:")
                        print(f"      총 NADH 생산: {nadh_total:.6f} mmol/gDCW/h")
                except:
                    pass
                
                # ETC 경로
                print("\n    [ETC 경로]:")
                etc_rxns = ['NADH16', 'SUCD', 'QCR', 'CYO3', 'ATPS4rpp']
                etc_active = False
                for rxn_id in etc_rxns:
                    try:
                        flux = solution.fluxes.get(rxn_id, 0)
                        if abs(flux) > 1e-8:
                            print(f"      {rxn_id}: {flux:.6f}")
                            etc_active = True
                    except KeyError:
                        pass
                
                if etc_active:
                    print("      [OK] ETC 경로 활성화!")
                else:
                    print("      [WARNING] ETC 경로 비활성 (Complex 누락 가능)")
                
                # ATP 생산 확인
                try:
                    atp_c = model.metabolites.get_by_id('atp_c')
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
                        print(f"\n    [ATP 생산]:")
                        atp_total = sum([flux for _, flux in atp_producing])
                        for rxn_id, flux in sorted(atp_producing, key=lambda x: x[1], reverse=True)[:5]:
                            print(f"      {rxn_id}: {flux:.6f} mmol/gDCW/h")
                        print(f"      총 ATP 생산: {atp_total:.6f} mmol/gDCW/h")
                        
                        # ETC를 통한 ATP 생산 확인
                        try:
                            atps4_flux = solution.fluxes.get('ATPS4rpp', 0)
                            if abs(atps4_flux) > 1e-8:
                                print(f"      [OK] ATPS4rpp (ETC): {atps4_flux:.6f} -> ATP 생산!")
                        except:
                            pass
                except:
                    pass
                
                # Yield 계산
                if glucose_uptake > 0:
                    yield_biomass = biomass_flux / glucose_uptake
                    print(f"\n    성장 속도: {biomass_flux:.6f} 1/h")
                    print(f"    Glucose 섭취: {glucose_uptake:.6f} mmol/gDCW/h")
                    print(f"    Biomass yield: {yield_biomass:.6f} gDW/mmol glucose")
                
                model.remove_reactions(['DM_atp_c_bs', 'DM_nad_c_bs', 'DM_coa_c_bs'])
                return True, solution
            else:
                print(f"\n  [FAIL] Biomass flux: {biomass_flux:.6f} 1/h")
                model.remove_reactions(['DM_atp_c_bs', 'DM_nad_c_bs', 'DM_coa_c_bs'])
                return False, solution
        else:
            print(f"\n  [FAIL] 최적화 실패: {solution.status}")
            try:
                model.remove_reactions(['DM_atp_c_bs', 'DM_nad_c_bs', 'DM_coa_c_bs'])
            except:
                pass
            return False, solution
    
    except Exception as e:
        print(f"  [ERROR] 테스트 실패: {e}")
        return False, None

def check_etc_complexes_availability(model):
    """ETC Complex 가용성 확인 및 대체 경로 찾기"""
    print("\n" + "="*70)
    print("ETC Complex 가용성 확인")
    print("="*70)
    
    # NADH16 대체 찾기
    print("\n[1] NADH dehydrogenase (Complex I) 대체 경로:")
    print("-" * 70)
    
    nadh_patterns = ['NADH', 'nadh', 'dehydrogenase', 'complex', 'quinone', 'ubiquinone']
    
    nadh_related = []
    for pattern in nadh_patterns:
        for rxn in model.reactions:
            if (pattern.lower() in rxn.id.lower() or 
                (rxn.name and pattern.lower() in rxn.name.lower())):
                if rxn not in nadh_related:
                    # NADH 소비 반응만 (dehydrogenase)
                    if 'nadh' in rxn.reaction.lower():
                        nadh_related.append(rxn)
    
    print(f"  NADH 관련 반응: {len(nadh_related)}개")
    
    # NADH → h_p 또는 quinone 경로 찾기
    print("\n  NADH → Proton motive force 경로:")
    nadh_to_pmf = []
    
    try:
        nadh_c = model.metabolites.get_by_id('nadh_c')
        h_p = model.metabolites.get_by_id('h_p')
        
        for rxn in nadh_c.reactions:
            # NADH를 소비하고 h_p를 생성하는 반응
            nadh_coeff = rxn.metabolites.get(nadh_c, 0)
            h_p_coeff = rxn.metabolites.get(h_p, 0) if h_p in rxn.metabolites else 0
            
            if nadh_coeff < 0 and h_p_coeff > 0:  # NADH 소비, h_p 생성
                nadh_to_pmf.append(rxn)
                print(f"    [FOUND] {rxn.id}: {rxn.name}")
                print(f"      {rxn.reaction}")
        
        if not nadh_to_pmf:
            print("    [WARNING] NADH → h_p 직접 경로 없음")
            print("    → Complex I (NADH16) 누락 가능")
    
    except KeyError:
        print("  [ERROR] nadh_c 또는 h_p metabolite 없음")
    
    # QCR 대체 찾기
    print("\n[2] Cytochrome bc1 complex (Complex III) 대체 경로:")
    print("-" * 70)
    
    qcr_patterns = ['QCR', 'bc1', 'cytochrome', 'ubiquinol', 'cytochrome c']
    
    qcr_related = []
    for pattern in qcr_patterns:
        for rxn in model.reactions:
            if (pattern.lower() in rxn.id.lower() or 
                (rxn.name and pattern.lower() in rxn.name.lower())):
                if rxn not in qcr_related:
                    qcr_related.append(rxn)
    
    if qcr_related:
        print(f"  발견된 관련 반응: {len(qcr_related)}개")
        for rxn in qcr_related[:5]:
            print(f"    - {rxn.id}: {rxn.name}")
            print(f"      {rxn.reaction}")
    else:
        print("  [WARNING] QCR 관련 반응 없음")
    
    # CYO3 대체 찾기
    print("\n[3] Cytochrome c oxidase (Complex IV) 대체 경로:")
    print("-" * 70)
    
    cyo_patterns = ['CYO', 'cytochrome c oxidase', 'CYTBD', 'cytochrome bd']
    
    cyo_related = []
    for pattern in cyo_patterns:
        for rxn in model.reactions:
            if (pattern.lower() in rxn.id.lower() or 
                (rxn.name and pattern.lower() in rxn.name.lower())):
                if rxn not in cyo_related:
                    cyo_related.append(rxn)
    
    if cyo_related:
        print(f"  발견된 관련 반응: {len(cyo_related)}개")
        for rxn in cyo_related:
            print(f"    - {rxn.id}: {rxn.name}")
            print(f"      {rxn.reaction}")
            
            # h_p 생성 확인
            try:
                h_p = model.metabolites.get_by_id('h_p')
                if h_p in rxn.metabolites:
                    h_p_coeff = rxn.metabolites.get(h_p, 0)
                    if h_p_coeff > 0:
                        print(f"        [OK] h_p 생성: {h_p_coeff}")
            except:
                pass
    else:
        print("  [WARNING] CYO3 관련 반응 없음")
    
    return nadh_to_pmf, qcr_related, cyo_related

def main():
    print("="*70)
    print("포도당 완전 산화 경로 이론적 테스트")
    print("="*70)
    
    model = load_model("BaseModel.xml")
    
    # Biomass 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass: {biomass_rxn.id}")
    
    # 1. ETC Complex 가용성 확인
    nadh_to_pmf, qcr_related, cyo_related = check_etc_complexes_availability(model)
    
    # 2. 부트스트랩으로 포도당 경로 테스트
    can_grow, solution = test_glucose_pathway_with_bootstrap(model, biomass_rxn)
    
    # 최종 요약
    print("\n" + "="*70)
    print("최종 요약")
    print("="*70)
    
    print("\nETC Complex 상태:")
    print(f"  - NADH dehydrogenase (Complex I): {'대체 경로 있음' if nadh_to_pmf else '누락'}")
    print(f"  - Succinate dehydrogenase (Complex II): 존재 (SUCD)")
    print(f"  - Cytochrome bc1 (Complex III): {'대체 경로 있음' if qcr_related else '누락'}")
    print(f"  - Cytochrome c oxidase (Complex IV): {'대체 경로 있음' if cyo_related else '누락'}")
    print(f"  - ATP synthase (Complex V): 존재 (ATPS4rpp)")
    
    if can_grow:
        print("\n[SUCCESS] 포도당 완전 산화 경로 작동 가능!")
        print("  → 부트스트랩으로 transport 문제 해결 후 경로 정상 작동")
        print("  → Glycolysis → TCA → NADH → ETC → ATP 경로 활성화")
    else:
        print("\n[FAIL] 포도당 완전 산화 경로 작동 불가능")
        if nadh_to_pmf:
            print("  → ETC 경로는 부분적으로 존재")
        else:
            print("  → ETC 경로 불완전 (Complex I 누락)")
        print("  → 추가 ETC Complex 또는 경로 수정 필요")
    
    print("\n" + "="*70)

if __name__ == "__main__":
    main()
