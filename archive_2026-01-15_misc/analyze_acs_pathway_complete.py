#!/usr/bin/env python
"""
ACS가 작동하지 않는 경로 완전성 분석
- ACS를 강제 활성화했을 때 경로가 어디서 막히는지 확인
- CoA 합성, ATP 생성 경로 등 확인
"""

import cobra
from pathlib import Path

def load_model(model_path):
    model = cobra.io.read_sbml_model(model_path)
    return model

def setup_acetate_medium(model):
    """Acetate 미디어 설정 (기본)"""
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

def analyze_with_acs_forced(model):
    """ACS를 강제 활성화했을 때 경로 분석"""
    print("="*80)
    print("ACS 강제 활성화 시 경로 분석")
    print("="*80)
    
    model = setup_acetate_medium(model)
    
    # ATPM=0 설정
    atpm_rxn = model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    # ACS 강제 활성화
    acs = model.reactions.get_by_id('ACS')
    acs.lower_bound = 0.1  # 최소 플럭스 강제
    
    model.objective = 'Growth'
    solution = model.optimize()
    
    print(f"\n[FBA 결과]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    print(f"  ACS 플럭스: {solution.fluxes.get('ACS', 0.0):.6f}")
    
    if solution.status == 'optimal':
        # 주요 반응 플럭스
        print(f"\n[주요 대사 경로 플럭스]")
        key_reactions = {
            'ACS': 'Acetate → Acetyl-CoA',
            'ADK1': 'AMP + ATP → 2ADP',
            'CS': 'Citrate synthase',
            'ICL': 'Isocitrate lyase (Glyoxylate)',
            'MALS': 'Malate synthase (Glyoxylate)',
            'ICDHx': 'Isocitrate dehydrogenase (TCA)',
            'AKGDH': 'Alpha-ketoglutarate dehydrogenase',
            'SUCDi': 'Succinate dehydrogenase',
            'MDH': 'Malate dehydrogenase',
        }
        
        for rxn_id, desc in key_reactions.items():
            try:
                flux = solution.fluxes.get(rxn_id, 0.0)
                if abs(flux) > 1e-6:
                    print(f"  {rxn_id:8s}: {flux:10.6f}  ({desc})")
            except:
                pass
        
        # ATP 생성/소모 분석
        print(f"\n[ATP 수지]")
        try:
            atp_c = model.metabolites.get_by_id('atp_c')
            atp_producing = []
            atp_consuming = []
            
            for rxn in atp_c.reactions:
                if 'ATPM' in rxn.id:
                    continue
                flux = solution.fluxes.get(rxn.id, 0.0)
                if abs(flux) > 1e-6:
                    coeff = rxn.metabolites.get(atp_c, 0)
                    net_flux = flux * coeff
                    if net_flux > 0:
                        atp_producing.append((rxn.id, net_flux))
                    elif net_flux < 0:
                        atp_consuming.append((rxn.id, abs(net_flux)))
            
            print(f"  ATP 생성 반응 (상위 5개):")
            for rxn_id, net_flux in sorted(atp_producing, key=lambda x: x[1], reverse=True)[:5]:
                print(f"    {rxn_id}: {net_flux:.6f}")
            print(f"  ATP 소모 반응 (상위 5개):")
            for rxn_id, net_flux in sorted(atp_consuming, key=lambda x: x[1], reverse=True)[:5]:
                print(f"    {rxn_id}: {net_flux:.6f}")
        except KeyError:
            print(f"  atp_c 메타볼라이트 없음")
        
        # CoA 관련 반응 확인
        print(f"\n[CoA 관련 반응]")
        try:
            coa_c = model.metabolites.get_by_id('coa_c')
            coa_producing = []
            coa_consuming = []
            
            for rxn in coa_c.reactions:
                flux = solution.fluxes.get(rxn.id, 0.0)
                if abs(flux) > 1e-6:
                    coeff = rxn.metabolites.get(coa_c, 0)
                    net_flux = flux * coeff
                    if net_flux > 0:
                        coa_producing.append((rxn.id, net_flux))
                    elif net_flux < 0:
                        coa_consuming.append((rxn.id, abs(net_flux)))
            
            print(f"  CoA 생성 반응:")
            for rxn_id, net_flux in sorted(coa_producing, key=lambda x: x[1], reverse=True)[:5]:
                print(f"    {rxn_id}: {net_flux:.6f}")
            print(f"  CoA 소모 반응:")
            for rxn_id, net_flux in sorted(coa_consuming, key=lambda x: x[1], reverse=True)[:5]:
                print(f"    {rxn_id}: {net_flux:.6f}")
        except KeyError:
            print(f"  coa_c 메타볼라이트 없음")
        
        # Exchange 플럭스 확인
        print(f"\n[주요 Exchange 플럭스]")
        key_exchanges = ['EX_ac_e', 'EX_o2_e', 'EX_co2_e', 'EX_nh4_e', 'EX_h2o_e']
        for ex_id in key_exchanges:
            try:
                flux = solution.fluxes.get(ex_id, 0.0)
                if abs(flux) > 1e-6:
                    print(f"  {ex_id}: {flux:.6f}")
            except:
                pass

def check_coa_synthesis(model):
    """CoA 합성 경로 확인"""
    print("\n" + "="*80)
    print("CoA 합성 경로 확인")
    print("="*80)
    
    # CoA 합성 관련 반응 찾기
    coa_synthesis_keywords = ['COA', 'PNTO', 'panto', 'CoA.*synth']
    
    coa_reactions = []
    for rxn in model.reactions:
        rxn_id_upper = rxn.id.upper()
        rxn_name_upper = (rxn.name or "").upper()
        
        for keyword in coa_synthesis_keywords:
            if keyword in rxn_id_upper or keyword in rxn_name_upper:
                try:
                    if 'coa_c' in [met.id for met in rxn.metabolites]:
                        coa_reactions.append(rxn)
                        break
                except:
                    pass
    
    print(f"\n[CoA 합성 관련 반응] {len(coa_reactions)}개 발견")
    for rxn in coa_reactions[:10]:  # 상위 10개만
        print(f"  {rxn.id}: {rxn.reaction[:80]}")

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    
    model = load_model(str(model_path))
    
    print("="*80)
    print("ACS 경로 완전성 분석")
    print("="*80)
    
    print("\n[사용자 지적]")
    print("  ACS + ADK1 조합으로 작동 가능해야 함")
    print("  에너지 효율이 낮으면 TCA cycle로 보상 (NADH 생성, 탄소 손실)")
    
    # ACS 강제 활성화 시 분석
    analyze_with_acs_forced(model)
    
    # CoA 합성 경로 확인
    check_coa_synthesis(model)
    
    print("\n" + "="*80)
    print("결론")
    print("="*80)
    
    print("\n[확인]")
    print("  ACS를 강제 활성화하면 플럭스가 발생하지만 성장률이 0")
    print("  -> 경로가 완전하지 않을 수 있음")
    print("  -> CoA 합성, ATP 생성, 또는 다른 필수 반응이 필요할 수 있음")

if __name__ == "__main__":
    main()
