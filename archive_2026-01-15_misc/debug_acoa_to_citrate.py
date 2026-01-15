#!/usr/bin/env python
"""
Acetyl-CoA가 Citrate로 안 가는 이유 디버깅
- CS (Citrate synthase) 반응 확인
- OAA (oxaloacetate) 확인
- ADK1 작동 여부 확인
"""

import cobra
from pathlib import Path

def load_model(model_path):
    model = cobra.io.read_sbml_model(model_path)
    return model

def setup_acetate_medium(model):
    """Acetate 미디어 설정"""
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    ex_ac = model.reactions.get_by_id('EX_ac_e')
    ex_ac.upper_bound = 1000
    ex_ac.lower_bound = -1000
    
    essential = {
        'EX_nh4_e': (-1000, 1000),
        'EX_h2o_e': (-1000, 1000),
        'EX_h_e': (-1000, 1000),
        'EX_pi_e': (-1000, 1000),
        'EX_so4_e': (-1000, 1000),
        'EX_k_e': (-1000, 1000),
        'EX_na1_e': (-1000, 1000),
        'EX_mg2_e': (-1000, 1000),
        'EX_ca2_e': (-1000, 1000),
        'EX_fe2_e': (-1000, 1000),
        'EX_mn2_e': (-1000, 1000),
        'EX_zn2_e': (-1000, 1000),
        'EX_co2_e': (-1000, 1000),
        'EX_o2_e': (-1000, 1000),
    }
    
    for ex_id, (lb, ub) in essential.items():
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            ex_rxn.upper_bound = ub
            ex_rxn.lower_bound = lb
        except KeyError:
            pass
    
    return model

def check_cs_reaction(model):
    """CS (Citrate synthase) 반응 확인"""
    print("="*80)
    print("CS (Citrate synthase) 반응 확인")
    print("="*80)
    
    try:
        cs = model.reactions.get_by_id('CS')
        print(f"\n[CS 반응]")
        print(f"  ID: {cs.id}")
        print(f"  이름: {cs.name}")
        print(f"  반응식: {cs.reaction}")
        print(f"  bounds: [{cs.lower_bound}, {cs.upper_bound}]")
        
        # 메타볼라이트 확인
        print(f"\n[CS 반응 메타볼라이트]")
        for met, coeff in cs.metabolites.items():
            print(f"  {met.id}: {coeff}")
        
        # 필요한 메타볼라이트: accoa_c, oaa_c, h2o_c
        required = ['accoa_c', 'oaa_c', 'h2o_c']
        produced = ['cit_c', 'coa_c', 'h_c']
        
        print(f"\n[필요한 메타볼라이트]")
        for met_id in required:
            try:
                met = model.metabolites.get_by_id(met_id)
                coeff = cs.metabolites.get(met, 0)
                print(f"  {met_id}: 계수={coeff}, 존재={True}")
            except KeyError:
                print(f"  {met_id}: 없음!")
        
        print(f"\n[생성되는 메타볼라이트]")
        for met_id in produced:
            try:
                met = model.metabolites.get_by_id(met_id)
                coeff = cs.metabolites.get(met, 0)
                print(f"  {met_id}: 계수={coeff}, 존재={True}")
            except KeyError:
                print(f"  {met_id}: 없음!")
        
    except KeyError:
        print(f"\n[CS 반응 없음]")

def check_oaa_generation(model):
    """OAA 생성 경로 확인"""
    print("\n" + "="*80)
    print("OAA (oxaloacetate) 생성 경로 확인")
    print("="*80)
    
    try:
        oaa_c = model.metabolites.get_by_id('oaa_c')
        
        # OAA를 생성하는 반응 찾기
        print(f"\n[OAA를 생성하는 반응]")
        oaa_producing = []
        for rxn in oaa_c.reactions:
            coeff = rxn.metabolites.get(oaa_c, 0)
            if coeff > 0:  # 생성
                oaa_producing.append((rxn.id, coeff))
        
        print(f"  총 {len(oaa_producing)}개 반응 발견")
        for rxn_id, coeff in oaa_producing[:10]:  # 상위 10개
            try:
                rxn = model.reactions.get_by_id(rxn_id)
                print(f"  {rxn_id}: 계수={coeff}, 반응식={rxn.reaction[:80]}")
            except:
                pass
        
        # 주요 OAA 생성 반응 확인
        key_oaa_reactions = ['PC', 'PEPCK', 'PEPCK_ATP', 'MDH']
        print(f"\n[주요 OAA 생성 반응 존재 여부]")
        for rxn_id in key_oaa_reactions:
            if rxn_id in [r.id for r in model.reactions]:
                rxn = model.reactions.get_by_id(rxn_id)
                print(f"  {rxn_id}: 있음 - {rxn.reaction[:80]}")
            else:
                print(f"  {rxn_id}: 없음")
        
    except KeyError:
        print(f"\n[oaa_c 메타볼라이트 없음]")

def test_acs_and_cs_together(model):
    """ACS와 CS를 함께 작동시키기"""
    print("\n" + "="*80)
    print("ACS와 CS 함께 작동 테스트")
    print("="*80)
    
    model = setup_acetate_medium(model)
    
    # ATPM=0 설정
    atpm_rxn = model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    # ACS 강제 활성화
    acs = model.reactions.get_by_id('ACS')
    acs.lower_bound = 0.1
    
    # CS 강제 활성화 (테스트)
    try:
        cs = model.reactions.get_by_id('CS')
        cs.lower_bound = 0.1
        print(f"  CS 하한 설정: 0.1")
    except KeyError:
        print(f"  CS 반응 없음")
        return None
    
    model.objective = 'Growth'
    solution = model.optimize()
    
    print(f"\n[FBA 결과]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    
    if solution.status == 'optimal':
        acs_flux = solution.fluxes.get('ACS', 0.0)
        cs_flux = solution.fluxes.get('CS', 0.0)
        adk1_flux = solution.fluxes.get('ADK1', 0.0)
        icl_flux = solution.fluxes.get('ICL', 0.0)
        mals_flux = solution.fluxes.get('MALS', 0.0)
        icdhx_flux = solution.fluxes.get('ICDHx', 0.0)
        
        print(f"\n[주요 반응 플럭스]")
        print(f"  ACS: {acs_flux:.6f}")
        print(f"  CS: {cs_flux:.6f}")
        print(f"  ADK1: {adk1_flux:.6f}")
        print(f"  ICL: {icl_flux:.6f}")
        print(f"  MALS: {mals_flux:.6f}")
        print(f"  ICDHx: {icdhx_flux:.6f}")
        
        # OAA 관련 플럭스 확인
        print(f"\n[OAA 관련 반응 플럭스]")
        oaa_reactions = ['PC', 'PEPCK', 'PEPCK_ATP', 'MDH']
        for rxn_id in oaa_reactions:
            try:
                flux = solution.fluxes.get(rxn_id, 0.0)
                if abs(flux) > 1e-6:
                    print(f"  {rxn_id}: {flux:.6f}")
            except:
                pass
        
        if abs(cs_flux) > 1e-6:
            print(f"\n[OK] CS가 작동합니다!")
        else:
            print(f"\n[문제] CS가 작동하지 않습니다")
            if abs(acs_flux) > 1e-6:
                print(f"  -> ACS는 작동하지만 CS는 작동 안 함")
                print(f"  -> OAA가 부족할 수 있음")
    
    return solution

def check_adk1_activation(model):
    """ADK1이 실제로 작동하는지 확인"""
    print("\n" + "="*80)
    print("ADK1 작동 여부 확인")
    print("="*80)
    
    model = setup_acetate_medium(model)
    
    # ATPM=0 설정
    atpm_rxn = model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    # ACS 강제 활성화
    acs = model.reactions.get_by_id('ACS')
    acs.lower_bound = 0.1
    
    # ADK1 강제 활성화 (테스트)
    adk1 = model.reactions.get_by_id('ADK1')
    adk1.lower_bound = 0.1  # ADK1도 강제 활성화
    
    model.objective = 'Growth'
    solution = model.optimize()
    
    print(f"\n[FBA 결과 (ADK1 강제 활성화)]")
    print(f"  상태: {solution.status}")
    
    if solution.status == 'optimal':
        acs_flux = solution.fluxes.get('ACS', 0.0)
        adk1_flux = solution.fluxes.get('ADK1', 0.0)
        cs_flux = solution.fluxes.get('CS', 0.0)
        
        print(f"  ACS: {acs_flux:.6f}")
        print(f"  ADK1: {adk1_flux:.6f}")
        print(f"  CS: {cs_flux:.6f}")
        
        if abs(adk1_flux) > 1e-6:
            print(f"\n[OK] ADK1이 작동합니다!")
        else:
            print(f"\n[문제] ADK1이 작동하지 않습니다")
            print(f"  -> ATP 부족일 수 있음")

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    
    model = load_model(str(model_path))
    
    print("="*80)
    print("Acetyl-CoA → Citrate 경로 디버깅")
    print("="*80)
    
    print("\n[사용자 질문]")
    print("  왜 Acetyl-CoA가 Citrate로 안 가는가?")
    print("  ADK1은 발현이 되어야 한다")
    
    # CS 반응 확인
    check_cs_reaction(model)
    
    # OAA 생성 경로 확인
    check_oaa_generation(model)
    
    # ACS와 CS 함께 테스트
    test_acs_and_cs_together(model)
    
    # ADK1 작동 여부 확인
    check_adk1_activation(model)
    
    print("\n" + "="*80)
    print("결론")
    print("="*80)
    
    print("\n[확인 필요]")
    print("  1. CS 반응이 존재하고 정상인지")
    print("  2. OAA를 생성할 수 있는지")
    print("  3. ADK1이 실제로 작동하는지")
    print("  4. ATP가 충분한지")

if __name__ == "__main__":
    main()
