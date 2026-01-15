#!/usr/bin/env python
"""
ATP 생산 경로 문제 분석
Acetate만으로 ATP 생산 가능 여부 확인
"""

import cobra
import pandas as pd

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    return model

def analyze_atp_production(model):
    """ATP 생산 경로 분석"""
    print("="*70)
    print("ATP 생산 경로 문제 분석")
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
    
    # 필수 영양소 (H2O 포함)
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
    
    # ATP 생산 테스트
    print("\n[테스트 1] ATP 생산 가능 여부 (Acetate 기반)")
    print("-" * 70)
    
    try:
        atp_c = model.metabolites.get_by_id('atp_c')
        
        # ATP sink 추가
        try:
            dm_atp = model.reactions.get_by_id('DM_atp_c')
            dm_atp.lower_bound = 0
            dm_atp.upper_bound = 1000
        except KeyError:
            dm_atp = cobra.Reaction('DM_atp_c')
            dm_atp.name = 'ATP demand'
            dm_atp.lower_bound = 0
            dm_atp.upper_bound = 1000
            dm_atp.add_metabolites({atp_c: -1})
            model.add_reactions([dm_atp])
        
        model.objective = dm_atp.id
        solution = model.optimize()
        
        if solution.status == 'optimal':
            atp_flux = solution.objective_value
            
            if atp_flux > 1e-6:
                print(f"  [SUCCESS] ATP 생산 가능: {atp_flux:.6f} mmol/gDCW/h")
                
                # ATP 생산 경로 확인
                atp_producing_rxns = []
                for rxn in atp_c.reactions:
                    if atp_c in rxn.products:
                        flux = solution.fluxes.get(rxn.id, 0)
                        if abs(flux) > 1e-8:
                            atp_producing_rxns.append((rxn.id, flux))
                
                print(f"\n  ATP 생산 반응 ({len(atp_producing_rxns)}개):")
                for rxn_id, flux in sorted(atp_producing_rxns, key=lambda x: abs(x[1]), reverse=True):
                    try:
                        rxn = model.reactions.get_by_id(rxn_id)
                        print(f"    {rxn_id}: {flux:.6f} ({rxn.name})")
                        print(f"      {rxn.reaction}")
                    except KeyError:
                        pass
                
                # Exchange 플럭스
                print(f"\n  Exchange 플럭스:")
                for rxn in model.exchanges:
                    flux = solution.fluxes.get(rxn.id, 0)
                    if abs(flux) > 0.001:
                        print(f"    {rxn.id}: {flux:.6f}")
                
                return True, solution
            else:
                print(f"  [FAIL] ATP 생산 불가능 (flux: {atp_flux:.6f})")
                
                # 왜 생산 못하는지 분석
                print(f"\n  ATP 생산 차단 원인 분석:")
                
                # ATP 생산 반응들 확인
                atp_producing_rxns_list = [rxn for rxn in atp_c.reactions if atp_c in rxn.products]
                print(f"    ATP 생산 반응 수: {len(atp_producing_rxns_list)}개")
                
                for rxn in atp_producing_rxns_list[:5]:
                    print(f"\n    - {rxn.id}: {rxn.name}")
                    print(f"      {rxn.reaction}")
                    
                    # 반응물 확인
                    reactants = {m: abs(c) for m, c in rxn.metabolites.items() if c < 0}
                    print(f"      필요 반응물: {[str(m) for m in reactants.keys()]}")
                    
                    # 각 반응물 생산 가능 여부 (간단 확인)
                    for met in reactants.keys():
                        flux = solution.fluxes.get(met.id, None)
                        if flux is None:
                            # 대사물질이 참여하는 반응 플럭스 확인
                            try:
                                met_producing = [r for r in met.reactions if met in r.products]
                                if not met_producing:
                                    print(f"        [BLOCKED] {met.id}: 생산 반응 없음")
                            except:
                                pass
                
                return False, solution
        else:
            print(f"  [FAIL] 최적화 실패: {solution.status}")
            return False, solution
    
    except KeyError:
        print("  [ERROR] atp_c metabolite 없음")
        return False, None

def check_atp_production_reactions(model):
    """ATP 생산 반응들 확인"""
    print("\n" + "="*70)
    print("ATP 생산 반응 확인")
    print("="*70)
    
    try:
        atp_c = model.metabolites.get_by_id('atp_c')
        atp_producing = [rxn for rxn in atp_c.reactions if atp_c in rxn.products]
        
        print(f"\nATP 생산 반응: {len(atp_producing)}개")
        print("-" * 70)
        
        reaction_info = []
        
        for rxn in atp_producing:
            reactants = [str(m) for m, c in rxn.metabolites.items() if c < 0]
            products = [str(m) for m, c in rxn.metabolites.items() if c > 0]
            
            reaction_info.append({
                'Reaction_ID': rxn.id,
                'Name': rxn.name,
                'Reactants': ', '.join(reactants),
                'Products': ', '.join(products),
                'Reversible': rxn.reversibility,
                'GPR': rxn.gene_reaction_rule
            })
            
            print(f"\n{rxn.id}: {rxn.name}")
            print(f"  {rxn.reaction}")
            print(f"  반응물: {reactants}")
            print(f"  생성물: {products}")
            print(f"  가역성: {rxn.reversibility}")
            print(f"  GPR: {rxn.gene_reaction_rule}")
        
        # 결과 저장
        if reaction_info:
            df = pd.DataFrame(reaction_info)
            df.to_csv('atp_production_reactions.csv', index=False)
            print(f"\n[OK] ATP 생산 반응 리스트 저장: atp_production_reactions.csv")
        
        return reaction_info
    
    except KeyError:
        print("[ERROR] atp_c metabolite 없음")
        return []

def test_oxidative_phosphorylation(model):
    """Oxidative phosphorylation 경로 확인"""
    print("\n" + "="*70)
    print("Oxidative Phosphorylation 경로 확인")
    print("="*70)
    
    # ATP synthase 확인
    atp_synthase_patterns = ['ATPS', 'ATPS4', 'ATP synthase', 'F0F1']
    
    print("\nATP synthase 반응:")
    atp_synthase_found = []
    
    for pattern in atp_synthase_patterns:
        for rxn in model.reactions:
            if pattern.lower() in rxn.id.lower() or (rxn.name and pattern.lower() in rxn.name.lower()):
                if rxn not in atp_synthase_found:
                    atp_synthase_found.append(rxn)
    
    if atp_synthase_found:
        for rxn in atp_synthase_found[:5]:
            print(f"  - {rxn.id}: {rxn.name}")
            print(f"    {rxn.reaction}")
    else:
        print("  [WARNING] ATP synthase 반응 없음")
    
    # NADH -> ATP 경로 확인
    print("\nNADH 생산 경로:")
    try:
        nadh_c = model.metabolites.get_by_id('nadh_c')
        nadh_producing = [rxn for rxn in nadh_c.reactions if nadh_c in rxn.products]
        print(f"  NADH 생산 반응: {len(nadh_producing)}개")
        
        # TCA cycle에서 NADH 생산 확인
        tca_nadh = ['ICDHx', 'AKGDH', 'MDH']
        for rxn_id in tca_nadh:
            try:
                rxn = model.reactions.get_by_id(rxn_id)
                if 'nadh' in rxn.reaction.lower():
                    print(f"    - {rxn_id}: {rxn.reaction}")
            except KeyError:
                pass
    
    except KeyError:
        print("  [ERROR] nadh_c metabolite 없음")

def main():
    print("="*70)
    print("ATP 생산 경로 문제 분석")
    print("="*70)
    
    model = load_model("BaseModel.xml")
    
    # 1. ATP 생산 가능 여부 테스트
    can_produce_atp, solution = analyze_atp_production(model)
    
    # 2. ATP 생산 반응 리스트
    atp_reactions = check_atp_production_reactions(model)
    
    # 3. Oxidative phosphorylation 확인
    test_oxidative_phosphorylation(model)
    
    # 최종 결론
    print("\n" + "="*70)
    print("결론")
    print("="*70)
    
    if can_produce_atp:
        print("\n[OK] Acetate만으로 ATP 생산 가능")
        print("  → ATP 생산 경로는 정상")
        print("  → 문제는 다른 부분 (예: 뉴클레오티드 생합성 경로)")
    else:
        print("\n[FAIL] Acetate만으로 ATP 생산 불가능")
        print("  → 이것이 근본 원인!")
        print("  → ATP 생산 경로가 불완전하거나 부트스트랩 문제")
        print("\n해결 방안:")
        print("  1. ATP 생산 경로 확인 및 수정")
        print("  2. 부트스트랩: 소량의 ATP 공급")
        print("  3. Oxidative phosphorylation 경로 확인")
    
    print("="*70)

if __name__ == "__main__":
    main()
