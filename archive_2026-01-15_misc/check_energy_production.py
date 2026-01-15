#!/usr/bin/env python
"""
아세트산 → ATP 생산 경로 확인
이론적으로는 가능해야 함: Acetate → Acetyl-CoA → TCA → ETC → ATP
"""

import cobra
from cobra import Reaction

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    return model

def setup_acetate_medium(model):
    """Acetate medium 설정"""
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

def check_atp_synthesis_pathway(model):
    """ATP 합성 경로 확인"""
    print("="*70)
    print("ATP 합성 경로 분석")
    print("="*70)
    
    model = setup_acetate_medium(model)
    
    # ATP 생산을 objective로 설정
    atp_c = model.metabolites.get_by_id('atp_c')
    
    # ATP 생산 테스트 반응
    test_atp = Reaction('TEST_atp_production')
    test_atp.lower_bound = 0
    test_atp.upper_bound = 1000
    test_atp.add_metabolites({atp_c: -1})
    model.add_reactions([test_atp])
    model.objective = 'TEST_atp_production'
    
    solution = model.optimize()
    
    print(f"\nATP 생산 테스트:")
    print(f"  상태: {solution.status}")
    print(f"  Objective (ATP 생산): {solution.objective_value:.6f}")
    
    if solution.status == 'optimal' and solution.objective_value > 1e-6:
        print(f"  [OK] ATP 생산 가능!")
        
        print(f"\n  주요 경로 플럭스:")
        key_rxns = ['EX_ac_e', 'ACt', 'ACS', 'CS', 'ACONT', 'ICDHx', 'ICDHyr', 
                   'AKGDH', 'SUCD', 'FUM', 'MDH', 'ICL', 'MALS']
        
        for rxn_id in key_rxns:
            try:
                flux = solution.fluxes.get(rxn_id, 0)
                if abs(flux) > 1e-6:
                    rxn = model.reactions.get_by_id(rxn_id)
                    print(f"    {rxn_id}: {flux:.4f} | {rxn.reaction}")
            except KeyError:
                pass
        
        # 전자 전달계 확인
        print(f"\n  전자 전달계 관련 플럭스:")
        etc_rxns = [r for r in model.reactions if 'nad' in r.id.lower() or 'fad' in r.id.lower() 
                   or 'q8' in r.id.lower() or 'cytochrome' in r.name.lower()]
        for rxn in etc_rxns[:10]:
            flux = solution.fluxes.get(rxn.id, 0)
            if abs(flux) > 1e-6:
                print(f"    {rxn.id}: {flux:.4f}")
        
    else:
        print(f"  [FAIL] ATP 생산 불가")
        print(f"\n  원인 분석:")
        
        # 각 단계별 확인
        print(f"  [1] Acetate uptake:")
        ex_ac_flux = solution.fluxes.get('EX_ac_e', 0)
        print(f"    EX_ac_e: {ex_ac_flux:.6f}")
        
        print(f"  [2] Transport:")
        ac_transport_flux = solution.fluxes.get('ACt', 0)
        print(f"    ACt: {ac_transport_flux:.6f}")
        
        print(f"  [3] ACS (Acetyl-CoA 생성):")
        acs_flux = solution.fluxes.get('ACS', 0)
        print(f"    ACS: {acs_flux:.6f}")
        if abs(acs_flux) < 1e-6:
            print(f"    [문제] ACS가 작동하지 않음")
            # ACS에 필요한 대사물질 확인
            try:
                acs = model.reactions.get_by_id('ACS')
                print(f"    필요한 대사물질:")
                for met, coeff in acs.metabolites.items():
                    if coeff < 0:
                        met_flux = solution.fluxes.get(met.id, 0)
                        print(f"      {met.id}: 필요량={-coeff}, 현재={met_flux:.6f}")
            except:
                pass
        
        print(f"  [4] TCA Cycle:")
        tca_rxns = ['CS', 'ACONT', 'ICDHx', 'AKGDH', 'SUCD', 'FUM', 'MDH']
        for rxn_id in tca_rxns:
            try:
                flux = solution.fluxes.get(rxn_id, 0)
                if abs(flux) > 1e-6:
                    print(f"    {rxn_id}: {flux:.4f} [OK]")
                else:
                    print(f"    {rxn_id}: {flux:.6f} [INACTIVE]")
            except KeyError:
                print(f"    {rxn_id}: [NOT FOUND]")
    
    model.remove_reactions([test_atp])
    return solution

def check_electron_transport_chain(model):
    """전자 전달계 확인"""
    print("\n" + "="*70)
    print("전자 전달계 (ETC) 확인")
    print("="*70)
    
    # NADH/FADH2를 ATP로 변환하는 경로 찾기
    try:
        nadh_c = model.metabolites.get_by_id('nadh_c')
        fadh2_c = model.metabolites.get_by_id('fadh2_c')
        atp_c = model.metabolites.get_by_id('atp_c')
        
        print(f"NADH 관련 반응: {len(nadh_c.reactions)}개")
        print(f"FADH2 관련 반응: {len(fadh2_c.reactions)}개")
        
        # NADH/FADH2를 소비하는 반응 찾기 (전자 전달계)
        nadh_consuming = [r for r in nadh_c.reactions if nadh_c in r.reactants]
        fadh2_consuming = [r for r in fadh2_c.reactions if fadh2_c in r.reactants]
        
        print(f"\nNADH 소비 반응 (전자 전달계 후보): {len(nadh_consuming)}개")
        for rxn in nadh_consuming[:10]:
            print(f"  {rxn.id}: {rxn.reaction}")
            # ATP와 관련이 있는지 확인
            if atp_c in rxn.metabolites:
                print(f"    [ATP 관련] 계수: {rxn.metabolites[atp_c]}")
        
        print(f"\nFADH2 소비 반응: {len(fadh2_consuming)}개")
        for rxn in fadh2_consuming[:10]:
            print(f"  {rxn.id}: {rxn.reaction}")
            if atp_c in rxn.metabolites:
                print(f"    [ATP 관련] 계수: {rxn.metabolites[atp_c]}")
        
        # ATP synthase 찾기
        atp_synthase = [r for r in model.reactions 
                       if 'synthase' in r.name.lower() and atp_c in r.products]
        if atp_synthase:
            print(f"\n[OK] ATP synthase 반응 발견: {len(atp_synthase)}개")
            for rxn in atp_synthase:
                print(f"  {rxn.id}: {rxn.reaction}")
        else:
            print(f"\n[WARNING] ATP synthase 반응을 명확히 찾을 수 없음")
            # ATP를 생성하는 반응 찾기
            atp_producing = [r for r in atp_c.reactions if atp_c in r.products]
            print(f"  ATP 생산 반응: {len(atp_producing)}개")
            # ADP를 사용하는 반응 찾기
            try:
                adp_c = model.metabolites.get_by_id('adp_c')
                adp_atp_rxns = [r for r in model.reactions 
                               if adp_c in r.reactants and atp_c in r.products]
                print(f"  ADP → ATP 반응: {len(adp_atp_rxns)}개")
                for rxn in adp_atp_rxns[:5]:
                    print(f"    {rxn.id}: {rxn.reaction}")
            except KeyError:
                pass
                
    except KeyError as e:
        print(f"[ERROR] {e} metabolite 없음")

def check_initial_atp_requirement(model):
    """초기 ATP 요구량 확인"""
    print("\n" + "="*70)
    print("초기 ATP 요구량 분석")
    print("="*70)
    
    model = setup_acetate_medium(model)
    
    try:
        acs = model.reactions.get_by_id('ACS')
        print(f"ACS 반응: {acs.reaction}")
        print(f"  필요한 ATP: 1 mol per Acetyl-CoA")
        
        # ACS가 1회 작동하면 생성되는 것
        print(f"  생성되는 것:")
        print(f"    - Acetyl-CoA: 1 mol")
        print(f"    - AMP: 1 mol")
        print(f"    - PPi: 1 mol")
        
        # TCA cycle에서 생성되는 ATP/NADH 확인
        print(f"\nTCA cycle (Acetyl-CoA 1 mol당):")
        print(f"  이론적 생성량:")
        print(f"    - NADH: ~3 mol (ICDH, AKGDH, MDH)")
        print(f"    - FADH2: ~1 mol (SUCD)")
        print(f"    - GTP: ~1 mol (SUCOAS)")
        print(f"    - CO2: ~2 mol")
        
        # 실제 모델에서 확인
        print(f"\n모델에서 TCA cycle 반응 확인:")
        tca_steps = {
            'CS': 'Acetyl-CoA + OAA → Citrate',
            'ACONT': 'Citrate → Isocitrate',
            'ICDHx': 'Isocitrate → α-KG + NADH',
            'AKGDH': 'α-KG → Succinyl-CoA + NADH',
            'SUCOAS': 'Succinyl-CoA → Succinate + GTP',
            'SUCD': 'Succinate → Fumarate + FADH2',
            'FUM': 'Fumarate → Malate',
            'MDH': 'Malate → OAA + NADH'
        }
        
        for rxn_id, desc in tca_steps.items():
            try:
                rxn = model.reactions.get_by_id(rxn_id)
                print(f"  {rxn_id}: {desc}")
                print(f"    반응식: {rxn.reaction}")
                
                # NADH/FADH2 생성 확인
                try:
                    nadh_c = model.metabolites.get_by_id('nadh_c')
                    if nadh_c in rxn.metabolites and rxn.metabolites[nadh_c] > 0:
                        print(f"    → NADH 생성: {rxn.metabolites[nadh_c]}")
                except:
                    pass
                
                try:
                    fadh2_c = model.metabolites.get_by_id('fadh2_c')
                    if fadh2_c in rxn.metabolites and fadh2_c.metabolites[fadh2_c] > 0:
                        print(f"    → FADH2 생성: {rxn.metabolites[fadh2_c]}")
                except:
                    pass
                    
            except KeyError:
                print(f"  {rxn_id}: [NOT FOUND]")
        
    except KeyError:
        print("[ERROR] ACS 반응 없음")

def test_minimal_bootstrap(model):
    """최소 bootstrap 테스트: 매우 작은 ATP 제공"""
    print("\n" + "="*70)
    print("최소 Bootstrap 테스트")
    print("="*70)
    
    model = setup_acetate_medium(model)
    
    # 매우 작은 ATP 제공 (0.001)
    atp_c = model.metabolites.get_by_id('atp_c')
    atp_boot = Reaction('EX_atp_minimal')
    atp_boot.lower_bound = -0.001  # 매우 작은 양
    atp_boot.upper_bound = 0
    atp_boot.add_metabolites({atp_c: -1})
    model.add_reactions([atp_boot])
    
    # ATP 생산을 objective로
    test_atp = Reaction('TEST_atp')
    test_atp.lower_bound = 0
    test_atp.upper_bound = 1000
    test_atp.add_metabolites({atp_c: -1})
    model.add_reactions([test_atp])
    model.objective = 'TEST_atp'
    
    solution = model.optimize()
    
    print(f"Bootstrap (ATP 0.001 제공) 후 ATP 생산:")
    print(f"  상태: {solution.status}")
    print(f"  ATP 생산: {solution.objective_value:.6f}")
    
    if solution.objective_value > 0.001:
        print(f"  [SUCCESS] Bootstrap으로 ATP 생산 가능!")
        print(f"    → 순환 가능: 초기 ATP → Acetyl-CoA → TCA → ATP")
        
        print(f"\n  경로 플럭스:")
        for rxn_id in ['EX_ac_e', 'ACt', 'ACS', 'CS', 'ICDHx', 'AKGDH', 'SUCOAS', 'MDH']:
            flux = solution.fluxes.get(rxn_id, 0)
            if abs(flux) > 1e-6:
                print(f"    {rxn_id}: {flux:.4f}")
    else:
        print(f"  [FAIL] Bootstrap으로도 ATP 생산 불가")
        print(f"    → 다른 문제가 있음 (OAA, CoA 등)")
    
    model.remove_reactions([atp_boot, test_atp])

def main():
    model = load_model("BaseModel.xml")
    
    check_atp_synthesis_pathway(model)
    check_electron_transport_chain(model)
    check_initial_atp_requirement(model)
    test_minimal_bootstrap(model)
    
    print("\n" + "="*70)
    print("에너지 생산 경로 분석 완료")
    print("="*70)

if __name__ == "__main__":
    main()

