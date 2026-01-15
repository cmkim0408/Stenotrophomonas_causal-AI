#!/usr/bin/env python
"""
Gluconeogenesis 및 PPP 경로 점검
1. Acetate로부터 PEP 생성 확인
2. Acetate로부터 Glucose-6P 생성 확인
3. E4P 생성 확인
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

def test_gluconeogenesis_pathway(model):
    """Gluconeogenesis 경로 테스트"""
    print("="*70)
    print("1. Gluconeogenesis 경로 점검")
    print("="*70)
    
    # 모든 exchange 초기화
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    # Acetate 설정
    try:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.lower_bound = -100
        ex_ac.upper_bound = 1000
        print(f"[OK] Acetate 설정: EX_ac_e")
    except KeyError:
        print("[ERROR] EX_ac_e 없음")
        return None
    
    # Glucose 차단
    try:
        ex_glc = model.reactions.get_by_id('EX_glc__D_e')
        ex_glc.lower_bound = 0
        ex_glc.upper_bound = 0
        print(f"[OK] Glucose 차단: EX_glc__D_e")
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
    
    results = {}
    
    # 1. PEP 생산 테스트
    print("\n[테스트 1] PEP 생산 가능 여부 (Acetate 기반)")
    print("-" * 70)
    
    try:
        pep_c = model.metabolites.get_by_id('pep_c')
        
        # PEP sink 추가
        try:
            dm_pep = model.reactions.get_by_id('DM_pep_c')
            dm_pep.lower_bound = 0
            dm_pep.upper_bound = 1000
        except KeyError:
            dm_pep = cobra.Reaction('DM_pep_c')
            dm_pep.name = 'PEP demand (sink)'
            dm_pep.lower_bound = 0
            dm_pep.upper_bound = 1000
            dm_pep.add_metabolites({pep_c: -1})
            model.add_reactions([dm_pep])
        
        model.objective = dm_pep.id
        solution = model.optimize()
        
        if solution.status == 'optimal':
            pep_flux = solution.objective_value
            results['PEP_production'] = {
                'status': 'optimal',
                'max_flux': pep_flux,
                'can_produce': pep_flux > 1e-6
            }
            
            if pep_flux > 1e-6:
                print(f"  [SUCCESS] PEP 생산 가능: {pep_flux:.6f} mmol/gDCW/h")
                
                # 주요 경로 플럭스 확인
                key_reactions = ['EX_ac_e', 'ACt', 'ACt2rpp', 'SUCOAACTr', 'ACS', 
                               'CS', 'MDH', 'PPS', 'ENO', 'PGM', 'PEPCK']
                
                print("\n  주요 경로 플럭스:")
                for rxn_id in key_reactions:
                    try:
                        flux = solution.fluxes.get(rxn_id, 0)
                        if abs(flux) > 1e-8:
                            print(f"    {rxn_id}: {flux:.6f}")
                    except KeyError:
                        pass
            else:
                print(f"  [FAIL] PEP 생산 불가능 (flux: {pep_flux:.6f})")
        else:
            results['PEP_production'] = {
                'status': solution.status,
                'max_flux': 0,
                'can_produce': False
            }
            print(f"  [FAIL] 최적화 실패: {solution.status}")
        
    except KeyError:
        print("  [ERROR] pep_c metabolite 없음")
        results['PEP_production'] = {'status': 'metabolite_missing', 'can_produce': False}
    
    # 2. Glucose-6P 생산 테스트
    print("\n[테스트 2] Glucose-6P 생산 가능 여부 (Acetate 기반)")
    print("-" * 70)
    
    try:
        g6p_c = model.metabolites.get_by_id('g6p_c')
        
        # G6P sink 추가
        try:
            dm_g6p = model.reactions.get_by_id('DM_g6p_c')
            dm_g6p.lower_bound = 0
            dm_g6p.upper_bound = 1000
        except KeyError:
            dm_g6p = cobra.Reaction('DM_g6p_c')
            dm_g6p.name = 'Glucose-6P demand (sink)'
            dm_g6p.lower_bound = 0
            dm_g6p.upper_bound = 1000
            dm_g6p.add_metabolites({g6p_c: -1})
            model.add_reactions([dm_g6p])
        
        model.objective = dm_g6p.id
        solution = model.optimize()
        
        if solution.status == 'optimal':
            g6p_flux = solution.objective_value
            results['G6P_production'] = {
                'status': 'optimal',
                'max_flux': g6p_flux,
                'can_produce': g6p_flux > 1e-6
            }
            
            if g6p_flux > 1e-6:
                print(f"  [SUCCESS] G6P 생산 가능: {g6p_flux:.6f} mmol/gDCW/h")
                
                # 주요 경로 플럭스 확인
                key_reactions = ['EX_ac_e', 'ACt', 'SUCOAACTr', 'ACS', 'CS', 'MDH',
                               'PPS', 'ENO', 'PGM', 'PEPCK', 'PC', 'PGI']
                
                print("\n  주요 경로 플럭스:")
                for rxn_id in key_reactions:
                    try:
                        flux = solution.fluxes.get(rxn_id, 0)
                        if abs(flux) > 1e-8:
                            print(f"    {rxn_id}: {flux:.6f}")
                    except KeyError:
                        pass
            else:
                print(f"  [FAIL] G6P 생산 불가능 (flux: {g6p_flux:.6f})")
        else:
            results['G6P_production'] = {
                'status': solution.status,
                'max_flux': 0,
                'can_produce': False
            }
            print(f"  [FAIL] 최적화 실패: {solution.status}")
        
    except KeyError:
        print("  [ERROR] g6p_c metabolite 없음")
        results['G6P_production'] = {'status': 'metabolite_missing', 'can_produce': False}
    
    # 3. E4P 생산 테스트 (PPP 경로)
    print("\n[테스트 3] E4P 생산 가능 여부 (PPP 경로)")
    print("-" * 70)
    
    try:
        e4p_c = model.metabolites.get_by_id('e4p_c')
        
        # E4P sink 추가
        try:
            dm_e4p = model.reactions.get_by_id('DM_e4p_c')
            dm_e4p.lower_bound = 0
            dm_e4p.upper_bound = 1000
        except KeyError:
            dm_e4p = cobra.Reaction('DM_e4p_c')
            dm_e4p.name = 'Erythrose-4P demand (sink)'
            dm_e4p.lower_bound = 0
            dm_e4p.upper_bound = 1000
            dm_e4p.add_metabolites({e4p_c: -1})
            model.add_reactions([dm_e4p])
        
        model.objective = dm_e4p.id
        solution = model.optimize()
        
        if solution.status == 'optimal':
            e4p_flux = solution.objective_value
            results['E4P_production'] = {
                'status': 'optimal',
                'max_flux': e4p_flux,
                'can_produce': e4p_flux > 1e-6
            }
            
            if e4p_flux > 1e-6:
                print(f"  [SUCCESS] E4P 생산 가능: {e4p_flux:.6f} mmol/gDCW/h")
                
                # PPP 경로 플럭스 확인
                ppp_reactions = ['G6PDH2r', 'PGL', 'GND', 'RPE', 'RPI', 'TKT1', 'TKT2', 'TALA']
                
                print("\n  PPP 경로 플럭스:")
                for rxn_id in ppp_reactions:
                    try:
                        flux = solution.fluxes.get(rxn_id, 0)
                        if abs(flux) > 1e-8:
                            print(f"    {rxn_id}: {flux:.6f}")
                    except KeyError:
                        pass
            else:
                print(f"  [FAIL] E4P 생산 불가능 (flux: {e4p_flux:.6f})")
        else:
            results['E4P_production'] = {
                'status': solution.status,
                'max_flux': 0,
                'can_produce': False
            }
            print(f"  [FAIL] 최적화 실패: {solution.status}")
        
    except KeyError:
        print("  [ERROR] e4p_c metabolite 없음")
        results['E4P_production'] = {'status': 'metabolite_missing', 'can_produce': False}
    
    return results

def check_gluconeogenesis_reactions(model):
    """Gluconeogenesis 주요 반응 확인"""
    print("\n" + "="*70)
    print("2. Gluconeogenesis 주요 반응 존재 여부 확인")
    print("="*70)
    
    gluconeogenesis_reactions = {
        'PPS': 'Phosphoenolpyruvate synthase',
        'PEPCK': 'Phosphoenolpyruvate carboxykinase',
        'PC': 'Pyruvate carboxylase',
        'ENO': 'Enolase',
        'PGM': 'Phosphoglycerate mutase',
        'PGK': 'Phosphoglycerate kinase',
        'GAPD': 'Glyceraldehyde-3-phosphate dehydrogenase',
        'TPI': 'Triosephosphate isomerase',
        'FBA': 'Fructose-bisphosphate aldolase',
        'FBP': 'Fructose-1,6-bisphosphatase',
        'PGI': 'Phosphoglucose isomerase',
        'HEX1': 'Hexokinase (역반응)'
    }
    
    print("\nGluconeogenesis 반응:")
    reaction_status = []
    
    for rxn_id, rxn_name in gluconeogenesis_reactions.items():
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            reaction_status.append({
                'Reaction_ID': rxn_id,
                'Name': rxn_name,
                'Exists': True,
                'Reversible': rxn.reversibility,
                'Equation': rxn.reaction
            })
            status = "[OK]"
            print(f"  {status} {rxn_id}: {rxn_name}")
            print(f"    반응식: {rxn.reaction}")
            print(f"    가역성: {rxn.reversibility}")
        except KeyError:
            reaction_status.append({
                'Reaction_ID': rxn_id,
                'Name': rxn_name,
                'Exists': False,
                'Reversible': False,
                'Equation': ''
            })
            print(f"  [MISSING] {rxn_id}: {rxn_name}")
    
    return reaction_status

def check_ppp_reactions(model):
    """PPP (Pentose Phosphate Pathway) 반응 확인"""
    print("\n" + "="*70)
    print("3. PPP 경로 반응 존재 여부 확인")
    print("="*70)
    
    ppp_reactions = {
        'G6PDH2r': 'Glucose-6-phosphate dehydrogenase',
        'PGL': '6-phosphogluconolactonase',
        'GND': '6-phosphogluconate dehydrogenase',
        'RPE': 'Ribulose-5-phosphate epimerase',
        'RPI': 'Ribose-5-phosphate isomerase',
        'TKT1': 'Transketolase (G6P + G3P -> E4P + F6P)',
        'TKT2': 'Transketolase (F6P + G3P -> E4P + Xu5P)',
        'TALA': 'Transaldolase'
    }
    
    print("\nPPP 반응:")
    reaction_status = []
    
    for rxn_id, rxn_name in ppp_reactions.items():
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            reaction_status.append({
                'Reaction_ID': rxn_id,
                'Name': rxn_name,
                'Exists': True,
                'Reversible': rxn.reversibility,
                'Equation': rxn.reaction
            })
            status = "[OK]"
            print(f"  {status} {rxn_id}: {rxn_name}")
            print(f"    반응식: {rxn.reaction}")
        except KeyError:
            reaction_status.append({
                'Reaction_ID': rxn_id,
                'Name': rxn_name,
                'Exists': False,
                'Reversible': False,
                'Equation': ''
            })
            print(f"  [MISSING] {rxn_id}: {rxn_name}")
    
    return reaction_status

def test_folate_pathway(model):
    """Folate 경로 테스트"""
    print("\n" + "="*70)
    print("4. Folate 경로 점검 (THF 생합성)")
    print("="*70)
    
    # 모든 exchange 초기화
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    # Acetate 설정
    try:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.lower_bound = -100
        ex_ac.upper_bound = 1000
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
    
    print("\n[테스트] THF 생산 가능 여부 (Acetate 기반)")
    print("-" * 70)
    
    try:
        thf_c = model.metabolites.get_by_id('thf_c')
        
        # THF sink 추가
        try:
            dm_thf = model.reactions.get_by_id('DM_thf_c')
            dm_thf.lower_bound = 0
            dm_thf.upper_bound = 1000
        except KeyError:
            dm_thf = cobra.Reaction('DM_thf_c')
            dm_thf.name = 'THF demand (sink)'
            dm_thf.lower_bound = 0
            dm_thf.upper_bound = 1000
            dm_thf.add_metabolites({thf_c: -1})
            model.add_reactions([dm_thf])
        
        model.objective = dm_thf.id
        solution = model.optimize()
        
        if solution.status == 'optimal':
            thf_flux = solution.objective_value
            
            if thf_flux > 1e-6:
                print(f"  [SUCCESS] THF 생산 가능: {thf_flux:.6f} mmol/gDCW/h")
                
                # Folate 경로 플럭스 확인
                folate_reactions = ['FOLR', 'FOLR2', 'DHFR', 'FOLM', 'FOLA']
                
                print("\n  Folate 경로 플럭스:")
                for rxn_id in folate_reactions:
                    try:
                        flux = solution.fluxes.get(rxn_id, 0)
                        if abs(flux) > 1e-8:
                            print(f"    {rxn_id}: {flux:.6f}")
                    except KeyError:
                        pass
                
                return {'status': 'optimal', 'max_flux': thf_flux, 'can_produce': True}
            else:
                print(f"  [FAIL] THF 생산 불가능 (flux: {thf_flux:.6f})")
                return {'status': 'optimal', 'max_flux': thf_flux, 'can_produce': False}
        else:
            print(f"  [FAIL] 최적화 실패: {solution.status}")
            return {'status': solution.status, 'max_flux': 0, 'can_produce': False}
        
    except KeyError:
        print("  [ERROR] thf_c metabolite 없음")
        return {'status': 'metabolite_missing', 'can_produce': False}

def main():
    print("="*70)
    print("Gluconeogenesis 및 PPP 경로 점검")
    print("="*70)
    
    # 모델 로드
    model = load_model("BaseModel.xml")
    
    # Biomass 찾기
    biomass_rxn = find_biomass_reaction(model)
    if biomass_rxn:
        print(f"[OK] Biomass: {biomass_rxn.id}")
    
    # 1. Gluconeogenesis 경로 테스트
    gluconeogenesis_results = test_gluconeogenesis_pathway(model)
    
    # 2. Gluconeogenesis 반응 존재 확인
    gluconeogenesis_reactions = check_gluconeogenesis_reactions(model)
    
    # 3. PPP 반응 존재 확인
    ppp_reactions = check_ppp_reactions(model)
    
    # 4. Folate 경로 테스트
    folate_results = test_folate_pathway(model)
    
    # 결과 저장
    if gluconeogenesis_results:
        df_gluc = pd.DataFrame([gluconeogenesis_results])
        df_gluc.to_csv('gluconeogenesis_test_results.csv', index=False)
        print(f"\n[OK] Gluconeogenesis 테스트 결과 저장: gluconeogenesis_test_results.csv")
    
    if gluconeogenesis_reactions:
        df_gluc_rxns = pd.DataFrame(gluconeogenesis_reactions)
        df_gluc_rxns.to_csv('gluconeogenesis_reactions_status.csv', index=False)
    
    if ppp_reactions:
        df_ppp = pd.DataFrame(ppp_reactions)
        df_ppp.to_csv('ppp_reactions_status.csv', index=False)
        print(f"[OK] PPP 반응 상태 저장: ppp_reactions_status.csv")
    
    # 요약
    print("\n" + "="*70)
    print("요약")
    print("="*70)
    
    if gluconeogenesis_results:
        print("\nGluconeogenesis 테스트 결과:")
        if gluconeogenesis_results.get('PEP_production', {}).get('can_produce'):
            print("  [OK] PEP 생산 가능")
        else:
            print("  [FAIL] PEP 생산 불가능")
        
        if gluconeogenesis_results.get('G6P_production', {}).get('can_produce'):
            print("  [OK] G6P 생산 가능")
        else:
            print("  [FAIL] G6P 생산 불가능")
        
        if gluconeogenesis_results.get('E4P_production', {}).get('can_produce'):
            print("  [OK] E4P 생산 가능")
        else:
            print("  [FAIL] E4P 생산 불가능")
    
    if folate_results:
        if folate_results.get('can_produce'):
            print("\nFolate 경로: [OK] THF 생산 가능")
        else:
            print("\nFolate 경로: [FAIL] THF 생산 불가능")
    
    print("="*70)

if __name__ == "__main__":
    main()
