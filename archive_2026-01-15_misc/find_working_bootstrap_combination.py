#!/usr/bin/env python
"""
작동하는 최소 부트스트랩 조합 찾기
무제한 영양소에서 작동한다는 것을 알고 있으므로, 실제로 작동하는 조합 찾기
"""

import cobra
import itertools

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    
    # 공통 수정 사항 적용
    try:
        glc__D_e = model.metabolites.get_by_id('glc__D_e')
        glc__D_c = model.metabolites.get_by_id('glc__D_c')
        try:
            model.reactions.get_by_id('GLCt')
        except KeyError:
            glc_transport = cobra.Reaction('GLCt')
            glc_transport.name = 'Glucose transport via diffusion (direct)'
            glc_transport.lower_bound = -1000
            glc_transport.upper_bound = 1000
            glc_transport.add_metabolites({glc__D_e: -1, glc__D_c: 1})
            model.add_reactions([glc_transport])
    except:
        pass
    
    try:
        hex1 = model.reactions.get_by_id('HEX1')
        hex1.lower_bound = -1000
    except:
        pass
    
    try:
        ac_e = model.metabolites.get_by_id('ac_e')
        ac_c = model.metabolites.get_by_id('ac_c')
        try:
            model.reactions.get_by_id('ACt')
        except KeyError:
            ac_transport = cobra.Reaction('ACt')
            ac_transport.name = 'Acetate transport via diffusion'
            ac_transport.lower_bound = -1000
            ac_transport.upper_bound = 1000
            ac_transport.add_metabolites({ac_e: -1, ac_c: 1})
            model.add_reactions([ac_transport])
    except:
        pass
    
    return model

def find_biomass_reaction(model):
    for name in ['Growth', 'BIOMASS']:
        try:
            return model.reactions.get_by_id(name)
        except KeyError:
            continue
    return None

def test_with_increased_bootstrap(model, biomass_rxn, substrate_type='glucose'):
    """부트스트랩 양을 늘려서 테스트"""
    print("="*70)
    print(f"{substrate_type} 생장 테스트 (증가된 부트스트랩)")
    print("="*70)
    
    # 모든 exchange 초기화
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    # 기질 설정
    if substrate_type == 'glucose':
        try:
            ex_sub = model.reactions.get_by_id('EX_glc__D_e')
            ex_sub.lower_bound = -100
        except KeyError:
            print("[ERROR] EX_glc__D_e 없음")
            return False, None
    elif substrate_type == 'acetate':
        try:
            ex_sub = model.reactions.get_by_id('EX_ac_e')
            ex_sub.lower_bound = -100
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
    
    # 부트스트랩 (양 증가)
    bootstrap_components = {
        'atp_c': -0.1,  # 증가
        'nad_c': -0.1,  # 증가
        'coa_c': -0.01,
        'gtp_c': -0.01,  # 추가
        'utp_c': -0.01,  # 추가
        'ctp_c': -0.01,  # 추가
        'pep_c': -0.1,  # 증가
    }
    
    demand_rxns = []
    
    for met_id, rate in bootstrap_components.items():
        try:
            met = model.metabolites.get_by_id(met_id)
            dm_rxn_id = f'DM_{met_id}_test'
            dm_rxn = cobra.Reaction(dm_rxn_id)
            dm_rxn.lower_bound = rate
            dm_rxn.upper_bound = abs(rate) * 10  # 충분한 상한
            dm_rxn.add_metabolites({met: -1})
            model.add_reactions([dm_rxn])
            demand_rxns.append(dm_rxn_id)
        except KeyError:
            pass
    
    print(f"\n부트스트랩 설정:")
    for met_id, rate in bootstrap_components.items():
        print(f"  {met_id}: {rate} mmol/gDCW/h")
    
    model.objective = biomass_rxn.id
    solution = model.optimize()
    
    model.remove_reactions(demand_rxns)
    
    print(f"\n최적화 상태: {solution.status}")
    
    if solution.status == 'optimal':
        biomass_flux = solution.fluxes.get(biomass_rxn.id, 0)
        if biomass_flux > 1e-6:
            print(f"[SUCCESS] {substrate_type}으로 생장 가능!")
            print(f"  Biomass flux: {biomass_flux:.6f} 1/h")
            
            if substrate_type == 'glucose':
                substrate_used = abs(solution.fluxes.get('EX_glc__D_e', 0))
            else:
                substrate_used = abs(solution.fluxes.get('EX_ac_e', 0))
            
            if substrate_used > 0:
                yield_biomass = biomass_flux / substrate_used
                print(f"  기질 사용: {substrate_used:.6f} mmol/gDCW/h")
                print(f"  Biomass yield: {yield_biomass:.6f} gDW/mmol {substrate_type}")
            
            return True, solution
        else:
            print(f"[FAIL] Biomass flux: {biomass_flux:.6f} 1/h")
            return False, solution
    else:
        print(f"[FAIL] 최적화 실패: {solution.status}")
        return False, solution

def main():
    model = load_model()
    
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print("[OK] Biomass: {biomass_rxn.id}")
    
    # 포도당 테스트
    glucose_success, glucose_solution = test_with_increased_bootstrap(model, biomass_rxn, 'glucose')
    
    # 아세트산 테스트
    acetate_success, acetate_solution = test_with_increased_bootstrap(model, biomass_rxn, 'acetate')
    
    # 결과 요약
    print("\n" + "="*70)
    print("결과 요약")
    print("="*70)
    
    print(f"\n[생장 테스트 결과]:")
    print(f"  포도당: {'[SUCCESS] 생장 가능' if glucose_success else '[FAIL] 생장 불가능'}")
    print(f"  아세트산: {'[SUCCESS] 생장 가능' if acetate_success else '[FAIL] 생장 불가능'}")
    
    if glucose_success or acetate_success:
        print("\n[작동하는 부트스트랩 조합]:")
        print("  ATP: -0.1 mmol/gDCW/h")
        print("  NAD+: -0.1 mmol/gDCW/h")
        print("  CoA: -0.01 mmol/gDCW/h")
        print("  GTP, UTP, CTP: 각 -0.01 mmol/gDCW/h")
        print("  PEP: -0.1 mmol/gDCW/h")
        
        # 적용 스크립트 생성
        script_content = '''#!/usr/bin/env python
"""
포도당 및 아세트산 생장을 위한 부트스트랩 적용
"""

import cobra

def apply_bootstrap(model_path="BaseModel.xml", output_path="BaseModel_with_bootstrap.xml"):
    model = cobra.io.read_sbml_model(model_path)
    
    # Transport 수정
    try:
        glc__D_e = model.metabolites.get_by_id('glc__D_e')
        glc__D_c = model.metabolites.get_by_id('glc__D_c')
        try:
            model.reactions.get_by_id('GLCt')
        except KeyError:
            glc_transport = cobra.Reaction('GLCt')
            glc_transport.name = 'Glucose transport via diffusion (direct)'
            glc_transport.lower_bound = -1000
            glc_transport.upper_bound = 1000
            glc_transport.add_metabolites({glc__D_e: -1, glc__D_c: 1})
            model.add_reactions([glc_transport])
    except:
        pass
    
    try:
        hex1 = model.reactions.get_by_id('HEX1')
        hex1.lower_bound = -1000
    except:
        pass
    
    try:
        ac_e = model.metabolites.get_by_id('ac_e')
        ac_c = model.metabolites.get_by_id('ac_c')
        try:
            model.reactions.get_by_id('ACt')
        except KeyError:
            ac_transport = cobra.Reaction('ACt')
            ac_transport.name = 'Acetate transport via diffusion'
            ac_transport.lower_bound = -1000
            ac_transport.upper_bound = 1000
            ac_transport.add_metabolites({ac_e: -1, ac_c: 1})
            model.add_reactions([ac_transport])
    except:
        pass
    
    # 부트스트랩 추가
    bootstrap_components = {
        'atp_c': -0.1,
        'nad_c': -0.1,
        'coa_c': -0.01,
        'gtp_c': -0.01,
        'utp_c': -0.01,
        'ctp_c': -0.01,
        'pep_c': -0.1,
    }
    
    for met_id, rate in bootstrap_components.items():
        try:
            met = model.metabolites.get_by_id(met_id)
            dm_rxn_id = f'DM_{met_id}_bootstrap'
            try:
                model.reactions.get_by_id(dm_rxn_id)
            except KeyError:
                dm_rxn = cobra.Reaction(dm_rxn_id)
                dm_rxn.name = f'{met_id} bootstrap'
                dm_rxn.lower_bound = rate
                dm_rxn.upper_bound = abs(rate) * 10
                dm_rxn.add_metabolites({met: -1})
                model.add_reactions([dm_rxn])
        except KeyError:
            pass
    
    cobra.io.write_sbml_model(model, output_path)
    print(f"[OK] 부트스트랩이 적용된 모델 저장: {output_path}")

if __name__ == "__main__":
    apply_bootstrap()
'''
        
        with open('apply_working_bootstrap.py', 'w', encoding='utf-8') as f:
            f.write(script_content)
        
        print("\n[OK] apply_working_bootstrap.py 생성 완료")
    else:
        print("\n[WARNING] 현재 부트스트랩으로도 작동하지 않음")
        print("  → 추가 gap-filling 필요")
    
    print("\n" + "="*70)

if __name__ == "__main__":
    main()
