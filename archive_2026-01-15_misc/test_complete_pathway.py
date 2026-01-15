#!/usr/bin/env python
"""
완전한 Acetate → Biomass 경로 테스트
단계별로 경로가 연결되는지 확인
"""

import cobra

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    return model

def setup_medium(model):
    for rxn in model.exchanges:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    
    essentials = {
        'EX_ac_e': (-1000, 1000),
        'EX_nh4_e': (-1000, 0),
        'EX_h2o_e': (-1000, 0),
        'EX_h_e': (-1000, 0),
        'EX_pi_e': (-1000, 0),
        'EX_so4_e': (-1000, 0),
        'EX_k_e': (-1000, 0),
        'EX_mg2_e': (-1000, 0),
        'EX_fe2_e': (-1000, 0),
        'EX_mn2_e': (-1000, 0),
        'EX_zn2_e': (-1000, 0),
        'EX_co2_e': (-1000, 1000),
        'EX_o2_e': (-1000, 1000)
    }
    
    for ex_id, bounds in essentials.items():
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            ex_rxn.lower_bound = bounds[0]
            ex_rxn.upper_bound = bounds[1]
        except KeyError:
            pass
    
    return model

def test_step_by_step(model):
    """단계별 경로 테스트"""
    print("="*70)
    print("단계별 경로 연결성 테스트")
    print("="*70)
    
    model = setup_medium(model)
    
    # Step 1: Acetate → Acetyl-CoA
    print("\n[Step 1] Acetate → Acetyl-CoA")
    try:
        from cobra import Reaction
        ac_c = model.metabolites.get_by_id('ac_c')
        accoa_c = model.metabolites.get_by_id('accoa_c')
        
        # 임시 반응: accoa_c 생산
        test_accoa = Reaction('TEST_accoa')
        test_accoa.lower_bound = 0
        test_accoa.upper_bound = 1000
        test_accoa.add_metabolites({accoa_c: -1})
        model.add_reactions([test_accoa])
        model.objective = 'TEST_accoa'
        
        solution = model.optimize()
        if solution.status == 'optimal' and solution.objective_value > 1e-6:
            print(f"  [OK] Acetyl-CoA 생산 가능: {solution.objective_value:.6f}")
            print(f"    EX_ac_e: {solution.fluxes.get('EX_ac_e', 0):.4f}")
            print(f"    ACt: {solution.fluxes.get('ACt', 0):.4f}")
            print(f"    ACS: {solution.fluxes.get('ACS', 0):.4f}")
        else:
            print(f"  [FAIL] Acetyl-CoA 생산 불가: {solution.objective_value}")
        
        model.remove_reactions([test_accoa])
        
    except Exception as e:
        print(f"  [ERROR] {e}")
    
    # Step 2: Acetyl-CoA → Malate (via Glyoxylate)
    print("\n[Step 2] Acetyl-CoA → Malate (Glyoxylate Shunt)")
    try:
        mal_c = model.metabolites.get_by_id('mal__L_c')
        test_mal = Reaction('TEST_mal')
        test_mal.lower_bound = 0
        test_mal.upper_bound = 1000
        test_mal.add_metabolites({mal_c: -1})
        model.add_reactions([test_mal])
        model.objective = 'TEST_mal'
        
        solution = model.optimize()
        if solution.status == 'optimal' and solution.objective_value > 1e-6:
            print(f"  [OK] Malate 생산 가능: {solution.objective_value:.6f}")
            print(f"    ICL: {solution.fluxes.get('ICL', 0):.4f}")
            print(f"    MALS: {solution.fluxes.get('MALS', 0):.4f}")
            print(f"    CS: {solution.fluxes.get('CS', 0):.4f}")
        else:
            print(f"  [FAIL] Malate 생산 불가: {solution.objective_value}")
        
        model.remove_reactions([test_mal])
        
    except Exception as e:
        print(f"  [ERROR] {e}")
    
    # Step 3: Malate → Pyruvate
    print("\n[Step 3] Malate → Pyruvate")
    try:
        pyr_c = model.metabolites.get_by_id('pyr_c')
        test_pyr = Reaction('TEST_pyr')
        test_pyr.lower_bound = 0
        test_pyr.upper_bound = 1000
        test_pyr.add_metabolites({pyr_c: -1})
        model.add_reactions([test_pyr])
        model.objective = 'TEST_pyr'
        
        solution = model.optimize()
        if solution.status == 'optimal' and solution.objective_value > 1e-6:
            print(f"  [OK] Pyruvate 생산 가능: {solution.objective_value:.6f}")
            print(f"    ME1: {solution.fluxes.get('ME1', 0):.4f}")
            print(f"    ME2: {solution.fluxes.get('ME2', 0):.4f}")
        else:
            print(f"  [FAIL] Pyruvate 생산 불가: {solution.objective_value}")
        
        model.remove_reactions([test_pyr])
        
    except Exception as e:
        print(f"  [ERROR] {e}")

def check_biomass_requirements(model):
    """Biomass requirement 확인"""
    print("\n" + "="*70)
    print("Biomass Reaction 필수 구성 요소 생성 가능 여부")
    print("="*70)
    
    model = setup_medium(model)
    model.objective = 'Growth'
    
    try:
        biomass = model.reactions.get_by_id('Growth')
        required_mets = {m.id: -coeff for m, coeff in biomass.metabolites.items() if coeff < 0}
        
        print(f"\nBiomass에 필요한 주요 구성 요소 확인:")
        
        # 각 대사물질이 생성 가능한지 테스트
        critical_mets = ['atp_c', 'gln__L_c', 'glu__L_c', 'gly_c', 'ala__L_c', 
                        'leu__L_c', 'arg__L_c', 'asp__L_c', 'asn__L_c', 'thr__L_c']
        
        for met_id in critical_mets:
            if met_id not in required_mets:
                continue
                
            try:
                met = model.metabolites.get_by_id(met_id)
                from cobra import Reaction
                
                test_rxn = Reaction(f'TEST_{met_id}')
                test_rxn.lower_bound = 0
                test_rxn.upper_bound = 1000
                test_rxn.add_metabolites({met: -1})
                model.add_reactions([test_rxn])
                model.objective = test_rxn.id
                
                solution = model.optimize()
                if solution.status == 'optimal' and solution.objective_value > 1e-6:
                    print(f"  [OK] {met_id}: 생산 가능 ({solution.objective_value:.6f})")
                else:
                    print(f"  [FAIL] {met_id}: 생산 불가")
                
                model.remove_reactions([test_rxn])
                
            except KeyError:
                print(f"  [NOT FOUND] {met_id}")
            except Exception as e:
                print(f"  [ERROR] {met_id}: {e}")
        
        # 원래 objective 복원
        model.objective = 'Growth'
        
    except KeyError:
        print("[ERROR] Growth 반응 없음")

def find_missing_reactions(model):
    """누락된 반응 찾기"""
    print("\n" + "="*70)
    print("누락 가능한 반응 확인")
    print("="*70)
    
    model = setup_medium(model)
    
    # OAA 생성 경로 확인 (CS에 필요)
    print("\nOAA (Oxaloacetate) 생성 경로:")
    try:
        oaa_c = model.metabolites.get_by_id('oaa_c')
        oaa_producing = [rxn for rxn in oaa_c.reactions if oaa_c in rxn.products or (rxn.reversibility and oaa_c in rxn.reactants)]
        print(f"  OAA 생산 반응: {len(oaa_producing)}개")
        for rxn in oaa_producing[:5]:
            print(f"    {rxn.id}: {rxn.reaction}")
        
        # OAA 생산 테스트
        from cobra import Reaction
        test_oaa = Reaction('TEST_oaa')
        test_oaa.lower_bound = 0
        test_oaa.upper_bound = 1000
        test_oaa.add_metabolites({oaa_c: -1})
        model.add_reactions([test_oaa])
        model.objective = 'TEST_oaa'
        
        solution = model.optimize()
        if solution.status == 'optimal' and solution.objective_value > 1e-6:
            print(f"  [OK] OAA 생산 가능: {solution.objective_value:.6f}")
        else:
            print(f"  [FAIL] OAA 생산 불가 - 이게 CS를 막고 있을 수 있음!")
        
        model.remove_reactions([test_oaa])
        
    except KeyError:
        print("  [ERROR] oaa_c 없음")
    
    # Isocitrate 생성 경로 확인 (ICL에 필요)
    print("\nIsocitrate 생성 경로:")
    try:
        icit_c = model.metabolites.get_by_id('icit_c')
        icit_producing = [rxn for rxn in icit_c.reactions if icit_c in rxn.products or (rxn.reversibility and icit_c in rxn.reactants)]
        print(f"  Isocitrate 생산 반응: {len(icit_producing)}개")
        for rxn in icit_producing[:5]:
            print(f"    {rxn.id}: {rxn.reaction}")
        
        # Aconitase 확인
        try:
            aconit = model.reactions.get_by_id('ACONT')
            print(f"  [OK] ACONT (Aconitase): {aconit.reaction}")
        except KeyError:
            # 다른 이름으로 찾기
            aconit_rxns = [r for r in model.reactions if 'aconit' in r.name.lower() or 'ACONT' in r.id]
            if aconit_rxns:
                print(f"  [OK] Aconitase 관련 반응 발견: {len(aconit_rxns)}개")
            else:
                print(f"  [WARNING] Aconitase 반응을 찾을 수 없음!")
        
    except KeyError:
        print("  [ERROR] icit_c 없음")

def main():
    model = load_model("BaseModel.xml")
    
    test_step_by_step(model)
    check_biomass_requirements(model)
    find_missing_reactions(model)
    
    print("\n" + "="*70)
    print("경로 분석 완료")
    print("="*70)

if __name__ == "__main__":
    main()

