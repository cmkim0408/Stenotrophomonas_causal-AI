#!/usr/bin/env python
"""
OAA 부트스트랩 제공 후 경로 막힘 지점 확인
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

def add_bootstrap_exchanges(model, bootstrap_amount=0.001):
    """OAA 초기화를 위한 부트스트랩 Exchange 추가"""
    print(f"\n[부트스트랩 추가] 양: {bootstrap_amount}")
    
    bootstrap_metabolites = {
        'oaa_c': 'OAA (oxaloacetate)',
        'mal__L_c': 'Malate',
        'asp__L_c': 'Aspartate',
    }
    
    added = []
    for met_id, met_name in bootstrap_metabolites.items():
        try:
            met = model.metabolites.get_by_id(met_id)
            ex_id = f'EX_{met_id}'
            
            # Exchange 반응이 없으면 생성
            if ex_id not in [r.id for r in model.exchanges]:
                ex_rxn = cobra.Reaction(ex_id)
                ex_rxn.name = f'{met_name} exchange'
                ex_rxn.lower_bound = -bootstrap_amount
                ex_rxn.upper_bound = 1000
                
                # 메타볼라이트의 extracellular 버전 찾기 또는 생성
                met_e_id = met_id.replace('_c', '_e')
                try:
                    met_e = model.metabolites.get_by_id(met_e_id)
                except KeyError:
                    met_e = cobra.Metabolite(met_e_id, name=met.name, compartment='e')
                    model.add_metabolites([met_e])
                
                ex_rxn.add_metabolites({met_e: -1.0})
                model.add_reactions([ex_rxn])
                added.append(ex_id)
                print(f"  추가: {ex_id} (하한: -{bootstrap_amount})")
            else:
                # Exchange가 이미 있으면 하한만 설정
                ex_rxn = model.reactions.get_by_id(ex_id)
                ex_rxn.lower_bound = -bootstrap_amount
                ex_rxn.upper_bound = 1000
                print(f"  설정: {ex_id} (하한: -{bootstrap_amount})")
        except KeyError:
            print(f"  건너뜀: {met_id} (메타볼라이트 없음)")
    
    return added

def diagnose_pathway_blockage(model):
    """경로 막힘 지점 진단"""
    print("\n" + "="*80)
    print("경로 막힘 지점 진단")
    print("="*80)
    
    # ATPM=0 설정
    atpm_rxn = model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    model.objective = 'Growth'
    solution = model.optimize()
    
    print(f"\n[FBA 결과]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    
    # 주요 반응 플럭스 확인
    key_reactions = {
        'ACS': 'Acetate → Acetyl-CoA',
        'CS': 'Acetyl-CoA + OAA → Citrate',
        'ADK1': 'AMP + ATP → 2ADP',
        'ICL': 'Isocitrate → Glyoxylate + Succinate',
        'MALS': 'Glyoxylate + Acetyl-CoA → Malate',
        'ICDHx': 'Isocitrate → α-KG (TCA)',
        'ACONT': 'Citrate → Isocitrate',
        'AKGDH': 'α-KG → Succinyl-CoA',
        'SUCDi': 'Succinate → Fumarate',
        'FUM': 'Fumarate → Malate',
        'MDH': 'Malate → OAA',
    }
    
    print(f"\n[주요 반응 플럭스]")
    active_reactions = []
    inactive_reactions = []
    
    for rxn_id, desc in key_reactions.items():
        try:
            flux = solution.fluxes.get(rxn_id, 0.0)
            if abs(flux) > 1e-6:
                print(f"  [OK] {rxn_id:8s}: {flux:10.6f}  ({desc})")
                active_reactions.append(rxn_id)
            else:
                print(f"  [NO] {rxn_id:8s}: {flux:10.6f}  ({desc})")
                inactive_reactions.append(rxn_id)
        except:
            print(f"  [??] {rxn_id:8s}: 반응 없음")
            inactive_reactions.append(rxn_id)
    
    # 경로 단계별 확인
    print(f"\n[경로 단계별 확인]")
    
    # 1단계: Acetate → Acetyl-CoA
    acs_flux = solution.fluxes.get('ACS', 0.0)
    if abs(acs_flux) > 1e-6:
        print(f"  [OK] 1단계: Acetate -> Acetyl-CoA (ACS: {acs_flux:.6f})")
    else:
        print(f"  [NO] 1단계: Acetate -> Acetyl-CoA 막힘 (ACS 플럭스: {acs_flux:.6f})")
        return "1단계: ACS"
    
    # 2단계: Acetyl-CoA → Citrate
    cs_flux = solution.fluxes.get('CS', 0.0)
    if abs(cs_flux) > 1e-6:
        print(f"  [OK] 2단계: Acetyl-CoA -> Citrate (CS: {cs_flux:.6f})")
    else:
        print(f"  [NO] 2단계: Acetyl-CoA -> Citrate 막힘 (CS 플럭스: {cs_flux:.6f})")
        # OAA 확인
        try:
            oaa_c = model.metabolites.get_by_id('oaa_c')
            oaa_producing = []
            for rxn in oaa_c.reactions:
                flux = solution.fluxes.get(rxn.id, 0.0)
                if abs(flux) > 1e-6:
                    coeff = rxn.metabolites.get(oaa_c, 0)
                    if coeff > 0:
                        oaa_producing.append((rxn.id, flux * coeff))
            
            if oaa_producing:
                print(f"    -> OAA 생성 반응 있음: {[r[0] for r in oaa_producing]}")
            else:
                print(f"    -> OAA 생성 반응 없음!")
        except:
            pass
        return "2단계: CS (OAA 부족 가능)"
    
    # 3단계: Citrate → Isocitrate
    aconit_flux = solution.fluxes.get('ACONT', 0.0)
    if abs(aconit_flux) > 1e-6:
        print(f"  [OK] 3단계: Citrate -> Isocitrate (ACONT: {aconit_flux:.6f})")
    else:
        print(f"  [NO] 3단계: Citrate -> Isocitrate 막힘")
        return "3단계: ACONT"
    
    # 4단계: 분기점 (ICL vs ICDHx)
    icl_flux = solution.fluxes.get('ICL', 0.0)
    icdhx_flux = solution.fluxes.get('ICDHx', 0.0)
    
    if abs(icl_flux) > 1e-6:
        print(f"  [OK] 4단계: Glyoxylate shunt 활성 (ICL: {icl_flux:.6f})")
    elif abs(icdhx_flux) > 1e-6:
        print(f"  [OK] 4단계: TCA cycle 활성 (ICDHx: {icdhx_flux:.6f})")
    else:
        print(f"  [NO] 4단계: Isocitrate 분기 막힘 (ICL: {icl_flux:.6f}, ICDHx: {icdhx_flux:.6f})")
        return "4단계: ICL/ICDHx"
    
    # 5단계: ADK1 확인
    adk1_flux = solution.fluxes.get('ADK1', 0.0)
    if abs(adk1_flux) > 1e-6:
        print(f"  [OK] ADK1: AMP 재생성 작동 (ADK1: {adk1_flux:.6f})")
    else:
        print(f"  [NO] ADK1: 작동 안 함 (ADK1 플럭스: {adk1_flux:.6f})")
        return "ADK1"
    
    # Exchange 플럭스 확인
    print(f"\n[부트스트랩 Exchange 플럭스]")
    bootstrap_exchanges = ['EX_oaa_c', 'EX_mal__L_c', 'EX_asp__L_c']
    for ex_id in bootstrap_exchanges:
        try:
            flux = solution.fluxes.get(ex_id, 0.0)
            if abs(flux) > 1e-6:
                print(f"  {ex_id}: {flux:.6f} (사용 중)")
        except:
            pass
    
    return None

def check_blocked_reactions(model):
    """Blocked reaction 확인"""
    print("\n" + "="*80)
    print("Blocked Reaction 확인")
    print("="*80)
    
    from cobra.flux_analysis import find_blocked_reactions
    
    model_copy = model.copy()
    atpm_rxn = model_copy.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    blocked = find_blocked_reactions(model_copy)
    
    key_reactions = ['ACS', 'CS', 'ADK1', 'ICL', 'MALS', 'ICDHx', 'ACONT', 'AKGDH', 'SUCDi', 'FUM', 'MDH']
    
    print(f"\n[주요 반응 Blocked 여부]")
    blocked_key = []
    for rxn_id in key_reactions:
        if rxn_id in blocked:
            print(f"  [NO] {rxn_id}: Blocked!")
            blocked_key.append(rxn_id)
        else:
            print(f"  [OK] {rxn_id}: Not blocked")
    
    if blocked_key:
        print(f"\n[Blocked 반응 발견] {len(blocked_key)}개")
        print(f"  -> {', '.join(blocked_key)}")
    else:
        print(f"\n[모든 주요 반응이 Not blocked]")

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    
    model = load_model(str(model_path))
    
    print("="*80)
    print("OAA 부트스트랩 제공 후 경로 막힘 지점 확인")
    print("="*80)
    
    # 미디어 설정
    model = setup_acetate_medium(model)
    
    # 부트스트랩 추가
    bootstrap_amount = 0.001  # 아주 약간의 부트스트랩
    added = add_bootstrap_exchanges(model, bootstrap_amount)
    
    # 경로 막힘 지점 진단
    blockage_point = diagnose_pathway_blockage(model)
    
    # Blocked reaction 확인
    check_blocked_reactions(model)
    
    print("\n" + "="*80)
    print("결론")
    print("="*80)
    
    if blockage_point:
        print(f"\n[막힘 지점] {blockage_point}")
        print(f"  -> 이 지점에서 경로가 막혀있음")
        print(f"  -> 추가 조사 필요")
    else:
        print(f"\n[경로 작동] 모든 단계가 작동함!")
        print(f"  -> 부트스트랩으로 경로가 활성화됨")
    
    print(f"\n[부트스트랩 정보]")
    print(f"  양: {bootstrap_amount}")
    print(f"  추가된 Exchange: {len(added)}개")

if __name__ == "__main__":
    main()
