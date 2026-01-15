#!/usr/bin/env python
"""
포도당 생장 문제 해결 방안 구현
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

def fix_glucose_transport(model):
    """방안 1: 포도당 transport 경로 수정 (ATP 필요 없게)"""
    print("="*70)
    print("방안 1: 포도당 Transport 경로 수정")
    print("="*70)
    
    # 방법 A: GLCabc를 가역적으로 만들기 (불가능 - 이미 단방향)
    # 방법 B: ATP 필요 없는 새로운 transport 경로 추가
    
    print("\n[방법 B] ATP 필요 없는 포도당 transport 경로 추가:")
    
    try:
        glc__D_e = model.metabolites.get_by_id('glc__D_e')
        glc__D_c = model.metabolites.get_by_id('glc__D_c')
        
        # 간단한 확산 경로: glc__D_e <=> glc__D_c
        try:
            glc_transport = model.reactions.get_by_id('GLCt')
            print(f"  [OK] GLCt 이미 존재: {glc_transport.reaction}")
        except KeyError:
            glc_transport = cobra.Reaction('GLCt')
            glc_transport.name = 'Glucose transport via diffusion (direct)'
            glc_transport.lower_bound = -1000  # 가역적
            glc_transport.upper_bound = 1000
            glc_transport.add_metabolites({
                glc__D_e: -1,
                glc__D_c: 1
            })
            model.add_reactions([glc_transport])
            print(f"  [OK] GLCt 추가: {glc_transport.reaction}")
        
        return True
    except Exception as e:
        print(f"  [ERROR] {e}")
        return False

def fix_hex1_issue(model):
    """방안 2: HEX1 문제 해결 - 가역적으로 만들기"""
    print("\n" + "="*70)
    print("방안 2: HEX1 가역성 문제 해결")
    print("="*70)
    
    try:
        hex1 = model.reactions.get_by_id('HEX1')
        
        # HEX1이 이미 가역적인지 확인
        if hex1.reversibility:
            print(f"  [OK] HEX1 이미 가역적")
            return True
        
        # 가역적으로 만들기 (G6P -> Glc + ATP)
        hex1.lower_bound = -1000
        print(f"  [OK] HEX1 하한 확장: {hex1.lower_bound}")
        print(f"  반응식: {hex1.reaction}")
        
        return True
    except Exception as e:
        print(f"  [ERROR] {e}")
        return False

def add_minimal_bootstrap_reactions(model):
    """방안 3: 최소 부트스트랩 반응 추가"""
    print("\n" + "="*70)
    print("방안 3: 최소 부트스트랩 반응 추가")
    print("="*70)
    
    print("\n초기 ATP, NAD+, CoA 생산을 위한 demand reaction 추가")
    print("(실제로는 모델 외부에서 공급되는 것으로 간주)")
    
    bootstrap_reactions = []
    
    # ATP demand (소량)
    try:
        atp_c = model.metabolites.get_by_id('atp_c')
        try:
            dm_atp = model.reactions.get_by_id('DM_atp_c_bootstrap')
            print(f"  [OK] DM_atp_c_bootstrap 이미 존재")
        except KeyError:
            dm_atp = cobra.Reaction('DM_atp_c_bootstrap')
            dm_atp.name = 'ATP bootstrap (minimal)'
            dm_atp.lower_bound = -0.01  # 매우 소량
            dm_atp.upper_bound = 0.01
            dm_atp.add_metabolites({atp_c: -1})
            model.add_reactions([dm_atp])
            bootstrap_reactions.append(dm_atp)
            print(f"  [OK] DM_atp_c_bootstrap 추가")
    except KeyError:
        print("  [ERROR] atp_c metabolite 없음")
    
    # NAD+ demand (소량)
    try:
        nad_c = model.metabolites.get_by_id('nad_c')
        try:
            dm_nad = model.reactions.get_by_id('DM_nad_c_bootstrap')
            print(f"  [OK] DM_nad_c_bootstrap 이미 존재")
        except KeyError:
            dm_nad = cobra.Reaction('DM_nad_c_bootstrap')
            dm_nad.name = 'NAD+ bootstrap (minimal)'
            dm_nad.lower_bound = -0.01
            dm_nad.upper_bound = 0.01
            dm_nad.add_metabolites({nad_c: -1})
            model.add_reactions([dm_nad])
            bootstrap_reactions.append(dm_nad)
            print(f"  [OK] DM_nad_c_bootstrap 추가")
    except KeyError:
        print("  [ERROR] nad_c metabolite 없음")
    
    # CoA demand (소량)
    try:
        coa_c = model.metabolites.get_by_id('coa_c')
        try:
            dm_coa = model.reactions.get_by_id('DM_coa_c_bootstrap')
            print(f"  [OK] DM_coa_c_bootstrap 이미 존재")
        except KeyError:
            dm_coa = cobra.Reaction('DM_coa_c_bootstrap')
            dm_coa.name = 'CoA bootstrap (minimal)'
            dm_coa.lower_bound = -0.001  # 매우 소량
            dm_coa.upper_bound = 0.001
            dm_coa.add_metabolites({coa_c: -1})
            model.add_reactions([dm_coa])
            bootstrap_reactions.append(dm_coa)
            print(f"  [OK] DM_coa_c_bootstrap 추가")
    except KeyError:
        print("  [ERROR] coa_c metabolite 없음")
    
    return len(bootstrap_reactions) > 0

def test_solutions(model, biomass_rxn):
    """해결 방안 테스트"""
    print("\n" + "="*70)
    print("해결 방안 테스트")
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
        print("\n[설정] 포도당: -100 mmol/gDCW/h")
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
    
    print("\n[FBA 수행 중...]")
    solution = model.optimize()
    
    print(f"\n최적화 상태: {solution.status}")
    
    if solution.status == 'optimal':
        biomass_flux = solution.fluxes.get(biomass_rxn.id, 0)
        
        if biomass_flux > 1e-6:
            print(f"\n[SUCCESS] 포도당으로 생장 가능!")
            print(f"  Biomass flux: {biomass_flux:.6f} 1/h")
            
            # 포도당 사용량
            try:
                glucose_used = abs(solution.fluxes.get('EX_glc__D_e', 0))
                if glucose_used > 0:
                    print(f"  포도당 사용: {glucose_used:.6f} mmol/gDCW/h")
                    yield_biomass = biomass_flux / glucose_used
                    print(f"  Biomass yield: {yield_biomass:.6f} gDW/mmol glucose")
            except:
                pass
            
            # Transport 경로 확인
            print("\n[Transport 경로 확인]:")
            transport_rxns = ['GLCt', 'GLCtex', 'GLCabc', 'GLCpts']
            for rxn_id in transport_rxns:
                try:
                    flux = solution.fluxes.get(rxn_id, 0)
                    if abs(flux) > 1e-8:
                        print(f"  {rxn_id}: {flux:.6f}")
                except KeyError:
                    pass
            
            # ATP 생산 확인
            try:
                atp_c = model.metabolites.get_by_id('atp_c')
                atp_total = 0
                for rxn in atp_c.reactions:
                    if atp_c in rxn.products:
                        flux = solution.fluxes.get(rxn.id, 0)
                        if abs(flux) > 1e-8:
                            atp_coeff = rxn.metabolites.get(atp_c, 0)
                            if atp_coeff > 0:
                                atp_produced = flux * atp_coeff
                                atp_total += atp_produced
                
                if atp_total > 0:
                    print(f"\n[ATP 생산]: {atp_total:.6f} mmol/gDCW/h")
                    
                    # ETC를 통한 ATP 생산
                    try:
                        atps4_flux = solution.fluxes.get('ATPS4rpp', 0)
                        if abs(atps4_flux) > 1e-8:
                            print(f"  ATPS4rpp (ETC): 활성화!")
                    except:
                        pass
            except:
                pass
            
            return True, solution
        else:
            print(f"\n[FAIL] Biomass flux: {biomass_flux:.6f} 1/h (생장 불가능)")
            return False, solution
    else:
        print(f"\n[FAIL] 최적화 실패: {solution.status}")
        return False, solution

def main():
    print("="*70)
    print("포도당 생장 문제 해결 방안")
    print("="*70)
    
    model = load_model("BaseModel.xml")
    
    # Biomass 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass: {biomass_rxn.id}")
    
    # 해결 방안 적용
    fix1 = fix_glucose_transport(model)
    fix2 = fix_hex1_issue(model)
    fix3 = add_minimal_bootstrap_reactions(model)
    
    # 테스트
    can_grow, solution = test_solutions(model, biomass_rxn)
    
    # 결과 요약
    print("\n" + "="*70)
    print("해결 방안 요약")
    print("="*70)
    
    print("\n[적용된 해결 방안]:")
    print(f"  1. 포도당 Transport 경로 수정: {'[OK]' if fix1 else '[FAIL]'}")
    print(f"  2. HEX1 가역성 문제 해결: {'[OK]' if fix2 else '[FAIL]'}")
    print(f"  3. 최소 부트스트랩 반응 추가: {'[OK]' if fix3 else '[FAIL]'}")
    
    if can_grow:
        print("\n[SUCCESS] 모든 해결 방안이 적용되었고 생장 가능!")
        
        # 모델 저장
        output_path = "BaseModel_fixed_glucose.xml"
        cobra.io.write_sbml_model(model, output_path)
        print(f"\n[OK] 수정된 모델 저장: {output_path}")
    else:
        print("\n[WARNING] 일부 해결 방안이 적용되었지만 여전히 생장 불가능")
        print("  → 추가 gap-filling 또는 다른 문제 존재")
    
    # 추가 권장 사항
    print("\n[추가 권장 사항]:")
    print("  1. 누락된 반응 추가 (gap-filling):")
    print("     - NAD+ 생산 경로 완성")
    print("     - CoA 생산 경로 완성")
    print("     - 뉴클레오티드 생합성 경로 완성")
    print("     - 아미노산 생합성 경로 완성")
    print("\n  2. 부트스트랩 전략 (실험적으로 검증):")
    print("     - 초기 소량의 ATP/NAD+/CoA 공급")
    print("     - 또는 부유한 배지에서 사전 배양")
    print("\n  3. Transport 경로 실험적 검증:")
    print("     - 실제로 어떤 transport 경로가 작동하는지 확인")
    print("     - ATP 필요 없는 transport 경로 존재 여부 확인")
    
    print("\n" + "="*70)

if __name__ == "__main__":
    main()
