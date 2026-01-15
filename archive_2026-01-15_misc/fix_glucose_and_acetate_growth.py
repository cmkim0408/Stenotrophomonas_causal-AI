#!/usr/bin/env python
"""
포도당 및 아세트산 생장 문제 해결
포도당과 아세트산 모두에서 작동하도록 모델 수정
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

def apply_common_fixes(model):
    """공통 수정 사항 적용"""
    print("="*70)
    print("공통 수정 사항 적용")
    print("="*70)
    
    fixes_applied = []
    
    # 1. 포도당 Transport 경로 추가
    try:
        glc__D_e = model.metabolites.get_by_id('glc__D_e')
        glc__D_c = model.metabolites.get_by_id('glc__D_c')
        try:
            model.reactions.get_by_id('GLCt')
            print("[OK] GLCt 이미 존재")
        except KeyError:
            glc_transport = cobra.Reaction('GLCt')
            glc_transport.name = 'Glucose transport via diffusion (direct)'
            glc_transport.lower_bound = -1000
            glc_transport.upper_bound = 1000
            glc_transport.add_metabolites({glc__D_e: -1, glc__D_c: 1})
            model.add_reactions([glc_transport])
            fixes_applied.append('GLCt')
            print("[OK] GLCt 추가")
    except Exception as e:
        print(f"[ERROR] GLCt 추가 실패: {e}")
    
    # 2. HEX1 가역성 확보
    try:
        hex1 = model.reactions.get_by_id('HEX1')
        if hex1.lower_bound == 0:
            hex1.lower_bound = -1000
            fixes_applied.append('HEX1')
            print("[OK] HEX1 가역성 확보")
        else:
            print("[OK] HEX1 이미 가역적")
    except KeyError:
        print("[WARNING] HEX1 없음")
    
    # 3. 아세트산 Transport 확인
    try:
        ac_e = model.metabolites.get_by_id('ac_e')
        ac_c = model.metabolites.get_by_id('ac_c')
        try:
            ac_transport = model.reactions.get_by_id('ACt')
            print(f"[OK] ACt 이미 존재: {ac_transport.reaction}")
        except KeyError:
            # 아세트산 확산 transport 추가
            ac_transport = cobra.Reaction('ACt')
            ac_transport.name = 'Acetate transport via diffusion'
            ac_transport.lower_bound = -1000
            ac_transport.upper_bound = 1000
            ac_transport.add_metabolites({ac_e: -1, ac_c: 1})
            model.add_reactions([ac_transport])
            fixes_applied.append('ACt')
            print("[OK] ACt 추가")
    except Exception as e:
        print(f"[WARNING] 아세트산 transport 확인 실패: {e}")
    
    return fixes_applied

def add_comprehensive_bootstrap(model):
    """포괄적 부트스트랩 추가"""
    print("\n" + "="*70)
    print("포괄적 부트스트랩 추가")
    print("="*70)
    
    bootstrap_components = {
        'atp_c': -0.01,
        'nad_c': -0.01,
        'coa_c': -0.001,
        'gtp_c': -0.001,
        'utp_c': -0.001,
        'ctp_c': -0.001,
        'pep_c': -0.01,
    }
    
    added_reactions = []
    
    for met_id, rate in bootstrap_components.items():
        try:
            met = model.metabolites.get_by_id(met_id)
            dm_rxn_id = f'DM_{met_id}_bootstrap'
            
            try:
                model.reactions.get_by_id(dm_rxn_id)
                print(f"  [SKIP] {dm_rxn_id} 이미 존재")
            except KeyError:
                dm_rxn = cobra.Reaction(dm_rxn_id)
                dm_rxn.name = f'{met_id} bootstrap (minimal)'
                dm_rxn.lower_bound = rate
                dm_rxn.upper_bound = abs(rate)  # 소량만 허용
                dm_rxn.add_metabolites({met: -1})
                model.add_reactions([dm_rxn])
                added_reactions.append(dm_rxn_id)
                print(f"  [OK] {dm_rxn_id} 추가: {rate} mmol/gDCW/h")
        except KeyError:
            print(f"  [SKIP] {met_id} metabolite 없음")
    
    return added_reactions

def test_glucose_growth(model, biomass_rxn):
    """포도당 생장 테스트"""
    print("\n" + "="*70)
    print("포도당 생장 테스트")
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
        print("[ERROR] EX_glc__D_e 없음")
        return False, None
    
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
    
    model.objective = biomass_rxn.id
    solution = model.optimize()
    
    print(f"\n최적화 상태: {solution.status}")
    
    if solution.status == 'optimal':
        biomass_flux = solution.fluxes.get(biomass_rxn.id, 0)
        if biomass_flux > 1e-6:
            print(f"[SUCCESS] 포도당으로 생장 가능!")
            print(f"  Biomass flux: {biomass_flux:.6f} 1/h")
            
            glucose_used = abs(solution.fluxes.get('EX_glc__D_e', 0))
            if glucose_used > 0:
                yield_biomass = biomass_flux / glucose_used
                print(f"  포도당 사용: {glucose_used:.6f} mmol/gDCW/h")
                print(f"  Biomass yield: {yield_biomass:.6f} gDW/mmol glucose")
            
            return True, solution
        else:
            print(f"[FAIL] Biomass flux: {biomass_flux:.6f} 1/h")
            return False, solution
    else:
        print(f"[FAIL] 최적화 실패: {solution.status}")
        return False, solution

def test_acetate_growth(model, biomass_rxn):
    """아세트산 생장 테스트"""
    print("\n" + "="*70)
    print("아세트산 생장 테스트")
    print("="*70)
    
    # 모든 exchange 초기화
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    # 아세트산 설정
    try:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.lower_bound = -100
        ex_ac.upper_bound = 1000
    except KeyError:
        print("[ERROR] EX_ac_e 없음")
        return False, None
    
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
    
    model.objective = biomass_rxn.id
    solution = model.optimize()
    
    print(f"\n최적화 상태: {solution.status}")
    
    if solution.status == 'optimal':
        biomass_flux = solution.fluxes.get(biomass_rxn.id, 0)
        if biomass_flux > 1e-6:
            print(f"[SUCCESS] 아세트산으로 생장 가능!")
            print(f"  Biomass flux: {biomass_flux:.6f} 1/h")
            
            acetate_used = abs(solution.fluxes.get('EX_ac_e', 0))
            if acetate_used > 0:
                yield_biomass = biomass_flux / acetate_used
                print(f"  아세트산 사용: {acetate_used:.6f} mmol/gDCW/h")
                print(f"  Biomass yield: {yield_biomass:.6f} gDW/mmol acetate")
            
            return True, solution
        else:
            print(f"[FAIL] Biomass flux: {biomass_flux:.6f} 1/h")
            return False, solution
    else:
        print(f"[FAIL] 최적화 실패: {solution.status}")
        return False, solution

def main():
    print("="*70)
    print("포도당 및 아세트산 생장 문제 해결")
    print("="*70)
    
    model = load_model("BaseModel.xml")
    
    # Biomass 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass: {biomass_rxn.id}")
    
    # 공통 수정 사항 적용
    fixes_applied = apply_common_fixes(model)
    
    # 부트스트랩 추가
    bootstrap_added = add_comprehensive_bootstrap(model)
    
    # 포도당 생장 테스트
    glucose_success, glucose_solution = test_glucose_growth(model, biomass_rxn)
    
    # 아세트산 생장 테스트
    acetate_success, acetate_solution = test_acetate_growth(model, biomass_rxn)
    
    # 결과 요약
    print("\n" + "="*70)
    print("결과 요약")
    print("="*70)
    
    print(f"\n[적용된 수정 사항]:")
    print(f"  Transport 경로: {', '.join(fixes_applied) if fixes_applied else '없음'}")
    print(f"  부트스트랩 반응: {len(bootstrap_added)}개 추가")
    
    print(f"\n[생장 테스트 결과]:")
    print(f"  포도당: {'[SUCCESS] 생장 가능' if glucose_success else '[FAIL] 생장 불가능'}")
    print(f"  아세트산: {'[SUCCESS] 생장 가능' if acetate_success else '[FAIL] 생장 불가능'}")
    
    if glucose_success or acetate_success:
        # 모델 저장
        output_path = "BaseModel_fixed_glucose_acetate.xml"
        cobra.io.write_sbml_model(model, output_path)
        print(f"\n[OK] 수정된 모델 저장: {output_path}")
    
    # 추가 권장 사항
    if not (glucose_success and acetate_success):
        print("\n[추가 권장 사항]:")
        if not glucose_success:
            print("  - 포도당 생장: 추가 gap-filling 필요")
        if not acetate_success:
            print("  - 아세트산 생장: 추가 gap-filling 필요")
        print("  - 누락된 반응 추가 (NAD+, CoA, 뉴클레오티드, 아미노산 생합성)")
    
    print("\n" + "="*70)

if __name__ == "__main__":
    main()
