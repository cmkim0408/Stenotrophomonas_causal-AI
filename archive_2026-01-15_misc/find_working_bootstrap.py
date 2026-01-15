#!/usr/bin/env python
"""
작동하는 최소 부트스트랩 조합 찾기
"""

import cobra
import itertools

def load_model(model_path="BaseModel_fixed_glucose.xml"):
    try:
        model = cobra.io.read_sbml_model(model_path)
    except:
        model = cobra.io.read_sbml_model("BaseModel.xml")
        
        # GLCt 추가
        try:
            glc__D_e = model.metabolites.get_by_id('glc__D_e')
            glc__D_c = model.metabolites.get_by_id('glc__D_c')
            glc_transport = cobra.Reaction('GLCt')
            glc_transport.name = 'Glucose transport via diffusion (direct)'
            glc_transport.lower_bound = -1000
            glc_transport.upper_bound = 1000
            glc_transport.add_metabolites({glc__D_e: -1, glc__D_c: 1})
            model.add_reactions([glc_transport])
        except:
            pass
        
        # HEX1 가역성
        try:
            hex1 = model.reactions.get_by_id('HEX1')
            hex1.lower_bound = -1000
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

def test_bootstrap_combinations(model, biomass_rxn):
    """다양한 부트스트랩 조합 테스트"""
    print("="*70)
    print("작동하는 최소 부트스트랩 조합 찾기")
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
    
    # 부트스트랩 가능한 구성 요소들
    bootstrap_candidates = {
        'atp_c': -0.1,
        'nad_c': -0.1,
        'coa_c': -0.01,
        'pep_c': -0.1,
        'gtp_c': -0.01,
        'utp_c': -0.01,
        'ctp_c': -0.01,
    }
    
    # 다양한 조합 테스트
    print("\n[1] 단일 구성 요소 테스트:")
    print("-" * 70)
    
    single_working = []
    for met_id, rate in bootstrap_candidates.items():
        try:
            met = model.metabolites.get_by_id(met_id)
            dm_rxn = cobra.Reaction(f'DM_{met_id}_test')
            dm_rxn.lower_bound = rate
            dm_rxn.upper_bound = 1000
            dm_rxn.add_metabolites({met: -1})
            model.add_reactions([dm_rxn])
            
            model.objective = biomass_rxn.id
            solution = model.optimize()
            
            model.remove_reactions([dm_rxn])
            
            if solution.status == 'optimal' and solution.objective_value > 1e-6:
                single_working.append(met_id)
                print(f"  [OK] {met_id}: Biomass = {solution.objective_value:.6f} 1/h")
            else:
                print(f"  [FAIL] {met_id}: {solution.status}")
        except:
            pass
    
    # 2개 조합 테스트
    print("\n[2] 2개 구성 요소 조합 테스트:")
    print("-" * 70)
    
    if not single_working:
        print("  단일 구성 요소로 작동하는 것 없음. 2개 조합 테스트...")
        
        pairs = list(itertools.combinations(bootstrap_candidates.keys(), 2))
        working_pairs = []
        
        for met1, met2 in pairs[:10]:  # 처음 10개만 테스트
            try:
                # 부트스트랩 추가
                demand_rxns = []
                for met_id in [met1, met2]:
                    met = model.metabolites.get_by_id(met_id)
                    rate = bootstrap_candidates[met_id]
                    dm_rxn = cobra.Reaction(f'DM_{met_id}_test')
                    dm_rxn.lower_bound = rate
                    dm_rxn.upper_bound = 1000
                    dm_rxn.add_metabolites({met: -1})
                    model.add_reactions([dm_rxn])
                    demand_rxns.append(dm_rxn.id)
                
                model.objective = biomass_rxn.id
                solution = model.optimize()
                
                model.remove_reactions(demand_rxns)
                
                if solution.status == 'optimal' and solution.objective_value > 1e-6:
                    working_pairs.append((met1, met2))
                    print(f"  [OK] {met1} + {met2}: Biomass = {solution.objective_value:.6f} 1/h")
                else:
                    print(f"  [FAIL] {met1} + {met2}: {solution.status}")
            except:
                pass
    
    # 필수 조합 테스트 (ATP + NAD+ + CoA)
    print("\n[3] 필수 조합 테스트 (ATP + NAD+ + CoA + 추가):")
    print("-" * 70)
    
    essential_combinations = [
        ['atp_c', 'nad_c', 'coa_c'],
        ['atp_c', 'nad_c', 'coa_c', 'pep_c'],
        ['atp_c', 'nad_c', 'coa_c', 'gtp_c'],
        ['atp_c', 'nad_c', 'coa_c', 'utp_c'],
        ['atp_c', 'nad_c', 'coa_c', 'gtp_c', 'utp_c', 'ctp_c'],
    ]
    
    working_combinations = []
    
    for combo in essential_combinations:
        try:
            demand_rxns = []
            for met_id in combo:
                met = model.metabolites.get_by_id(met_id)
                rate = bootstrap_candidates[met_id]
                dm_rxn = cobra.Reaction(f'DM_{met_id}_test')
                dm_rxn.lower_bound = rate
                dm_rxn.upper_bound = 1000
                dm_rxn.add_metabolites({met: -1})
                model.add_reactions([dm_rxn])
                demand_rxns.append(dm_rxn.id)
            
            model.objective = biomass_rxn.id
            solution = model.optimize()
            
            model.remove_reactions(demand_rxns)
            
            combo_str = ' + '.join(combo)
            if solution.status == 'optimal' and solution.objective_value > 1e-6:
                working_combinations.append((combo, solution.objective_value))
                print(f"  [OK] {combo_str}: Biomass = {solution.objective_value:.6f} 1/h")
            else:
                print(f"  [FAIL] {combo_str}: {solution.status}")
        except Exception as e:
            print(f"  [ERROR] {combo}: {e}")
    
    # 최적 조합 선택 및 스크립트 생성
    if working_combinations:
        # 가장 적은 구성 요소로 작동하는 조합 선택
        best_combo = min(working_combinations, key=lambda x: len(x[0]))
        
        print("\n" + "="*70)
        print(f"[최적 조합 발견]")
        print("="*70)
        print(f"구성 요소: {', '.join(best_combo[0])}")
        print(f"Biomass flux: {best_combo[1]:.6f} 1/h")
        
        # 적용 스크립트 생성
        print("\n[적용 스크립트 생성 중...]")
        script_content = f"""#!/usr/bin/env python
\"\"\"
포도당 생장을 위한 부트스트랩 적용 스크립트
\"\"\"

import cobra

def apply_glucose_bootstrap(model_path="BaseModel.xml", output_path="BaseModel_with_glucose_bootstrap.xml"):
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
            glc_transport.add_metabolites({{glc__D_e: -1, glc__D_c: 1}})
            model.add_reactions([glc_transport])
    except:
        pass
    
    # HEX1 가역성
    try:
        hex1 = model.reactions.get_by_id('HEX1')
        hex1.lower_bound = -1000
    except:
        pass
    
    # 부트스트랩 반응 추가
    bootstrap_components = {best_combo[0]}
    bootstrap_rates = {bootstrap_candidates}
    
    for met_id in bootstrap_components:
        try:
            met = model.metabolites.get_by_id(met_id)
            dm_rxn_id = f'DM_{{met_id}}_bootstrap'
            try:
                model.reactions.get_by_id(dm_rxn_id)
            except KeyError:
                dm_rxn = cobra.Reaction(dm_rxn_id)
                dm_rxn.name = f'{{met_id}} bootstrap'
                dm_rxn.lower_bound = bootstrap_rates[met_id]
                dm_rxn.upper_bound = 1000
                dm_rxn.add_metabolites({{met: -1}})
                model.add_reactions([dm_rxn])
        except KeyError:
            pass
    
    cobra.io.write_sbml_model(model, output_path)
    print(f"[OK] 부트스트랩이 적용된 모델 저장: {{output_path}}")

if __name__ == "__main__":
    apply_glucose_bootstrap()
"""
        
        with open('apply_glucose_bootstrap.py', 'w', encoding='utf-8') as f:
            f.write(script_content)
        
        print("[OK] apply_glucose_bootstrap.py 생성 완료")
        
        return best_combo
    else:
        print("\n[FAIL] 작동하는 부트스트랩 조합을 찾지 못했습니다")
        print("  → 추가 gap-filling 또는 모델 구조 수정 필요")
        return None

def main():
    model = load_model()
    
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    best_combo = test_bootstrap_combinations(model, biomass_rxn)
    
    print("\n" + "="*70)
    print("최종 해결 방안")
    print("="*70)
    
    if best_combo:
        print("\n[해결 방법]:")
        print(f"  1. 포도당 transport 경로 추가 (GLCt)")
        print(f"  2. HEX1 가역성 확보")
        print(f"  3. 부트스트랩 적용: {', '.join(best_combo[0])}")
        print(f"\n적용 스크립트: apply_glucose_bootstrap.py")
        print(f"\n사용법:")
        print(f"  python apply_glucose_bootstrap.py")
        print(f"  또는 코드에서 apply_glucose_bootstrap() 함수 호출")
    else:
        print("\n[추가 해결 방안 필요]:")
        print("  1. 누락된 반응 추가 (gap-filling)")
        print("  2. 모델 구조 검토 및 수정")
        print("  3. 실험적 검증 필요")
    
    print("\n" + "="*70)

if __name__ == "__main__":
    main()
