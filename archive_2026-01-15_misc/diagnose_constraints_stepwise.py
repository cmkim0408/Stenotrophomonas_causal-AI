#!/usr/bin/env python
"""
제약 조건 단계별 진단
무제한 영양소에서 점진적으로 제약 조건을 추가하며 문제 찾기
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

def test_growth_stepwise(model, biomass_rxn):
    """단계별 제약 조건 테스트"""
    print("="*70)
    print("제약 조건 단계별 진단")
    print("="*70)
    
    results = []
    
    # 단계 1: 무제한 영양소
    print("\n[단계 1] 무제한 영양소")
    for rxn in model.exchanges:
        rxn.lower_bound = -1000
        rxn.upper_bound = 1000
    
    model.objective = biomass_rxn.id
    solution1 = model.optimize()
    
    if solution1.status == 'optimal':
        biomass_flux = solution1.fluxes.get(biomass_rxn.id, 0)
        results.append({
            'Step': '1. 무제한 영양소',
            'Status': solution1.status,
            'Biomass_flux': biomass_flux,
            'Can_Grow': biomass_flux > 1e-6
        })
        print(f"  상태: {solution1.status}")
        print(f"  Biomass flux: {biomass_flux:.6f} 1/h")
        if biomass_flux > 1e-6:
            print("  [SUCCESS] 성장 가능!")
        else:
            print("  [FAIL] 성장 불가")
    else:
        results.append({
            'Step': '1. 무제한 영양소',
            'Status': solution1.status,
            'Biomass_flux': 0,
            'Can_Grow': False
        })
        print(f"  상태: {solution1.status}")
    
    # 단계 2: 포도당만 허용
    print("\n[단계 2] 포도당만 허용 (필수 영양소 포함)")
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    # 포도당 허용
    try:
        ex_glc = model.reactions.get_by_id('EX_glc__D_e')
        ex_glc.lower_bound = -100
        ex_glc.upper_bound = 1000
    except KeyError:
        for rxn in model.exchanges:
            if 'glc' in rxn.id.lower() or 'glucose' in rxn.id.lower():
                rxn.lower_bound = -100
                rxn.upper_bound = 1000
                break
    
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
    
    solution2 = model.optimize()
    
    if solution2.status == 'optimal':
        biomass_flux = solution2.fluxes.get(biomass_rxn.id, 0)
        results.append({
            'Step': '2. 포도당만 허용',
            'Status': solution2.status,
            'Biomass_flux': biomass_flux,
            'Can_Grow': biomass_flux > 1e-6
        })
        print(f"  상태: {solution2.status}")
        print(f"  Biomass flux: {biomass_flux:.6f} 1/h")
        if biomass_flux > 1e-6:
            print("  [SUCCESS] 포도당으로 성장 가능!")
        else:
            print("  [FAIL] 포도당으로도 성장 불가")
    else:
        results.append({
            'Step': '2. 포도당만 허용',
            'Status': solution2.status,
            'Biomass_flux': 0,
            'Can_Grow': False
        })
        print(f"  상태: {solution2.status}")
        print("  [FAIL] 포도당으로도 성장 불가")
    
    # 단계 3: Acetate만 허용
    print("\n[단계 3] Acetate만 허용 (필수 영양소 포함)")
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    # Acetate 허용
    try:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.lower_bound = -100
        ex_ac.upper_bound = 1000
    except KeyError:
        for rxn in model.exchanges:
            if 'ac' in rxn.id.lower() and '_e' in rxn.id:
                rxn.lower_bound = -100
                rxn.upper_bound = 1000
                break
    
    # 필수 영양소
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
    
    solution3 = model.optimize()
    
    if solution3.status == 'optimal':
        biomass_flux = solution3.fluxes.get(biomass_rxn.id, 0)
        results.append({
            'Step': '3. Acetate만 허용',
            'Status': solution3.status,
            'Biomass_flux': biomass_flux,
            'Can_Grow': biomass_flux > 1e-6
        })
        print(f"  상태: {solution3.status}")
        print(f"  Biomass flux: {biomass_flux:.6f} 1/h")
        if biomass_flux > 1e-6:
            print("  [SUCCESS] Acetate로 성장 가능!")
        else:
            print("  [FAIL] Acetate로는 성장 불가 (부트스트랩 문제 가능)")
    else:
        results.append({
            'Step': '3. Acetate만 허용',
            'Status': solution3.status,
            'Biomass_flux': 0,
            'Can_Grow': False
        })
        print(f"  상태: {solution3.status}")
        print("  [FAIL] Acetate로는 성장 불가")
    
    # 단계 4: Acetate + CoA demand (부트스트랩)
    print("\n[단계 4] Acetate + CoA demand (부트스트랩)")
    
    # CoA demand 추가
    try:
        coa_c = model.metabolites.get_by_id('coa_c')
        try:
            dm_coa = model.reactions.get_by_id('DM_coa_c')
            dm_coa.lower_bound = -0.1
        except KeyError:
            dm_coa = cobra.Reaction('DM_coa_c')
            dm_coa.name = 'CoA demand (bootstrap)'
            dm_coa.lower_bound = -0.1
            dm_coa.upper_bound = 1000
            dm_coa.add_metabolites({coa_c: -1})
            model.add_reactions([dm_coa])
    except KeyError:
        pass
    
    solution4 = model.optimize()
    
    if solution4.status == 'optimal':
        biomass_flux = solution4.fluxes.get(biomass_rxn.id, 0)
        results.append({
            'Step': '4. Acetate + CoA demand',
            'Status': solution4.status,
            'Biomass_flux': biomass_flux,
            'Can_Grow': biomass_flux > 1e-6
        })
        print(f"  상태: {solution4.status}")
        print(f"  Biomass flux: {biomass_flux:.6f} 1/h")
        if biomass_flux > 1e-6:
            print("  [SUCCESS] CoA 부트스트랩으로 성장 가능!")
        else:
            print("  [FAIL] CoA 부트스트랩으로도 성장 불가")
    else:
        results.append({
            'Step': '4. Acetate + CoA demand',
            'Status': solution4.status,
            'Biomass_flux': 0,
            'Can_Grow': False
        })
        print(f"  상태: {solution4.status}")
    
    df_results = pd.DataFrame(results)
    
    return df_results, solution3, solution4

def analyze_acetate_fluxes(model, solution):
    """Acetate 경로 플럭스 분석"""
    if not solution or solution.status != 'optimal':
        return
    
    print("\n" + "="*70)
    print("Acetate 경로 플럭스 분석")
    print("="*70)
    
    # Acetate 관련 반응
    acetate_pathway = {
        'EX_ac_e': 'Acetate Exchange',
        'ACt': 'Acetate Transport',
        'ACt2rpp': 'Acetate Transport (periplasm)',
        'SUCOAACTr': 'Acetate-CoA Transferase',
        'ACS': 'Acetyl-CoA Synthetase',
        'CS': 'Citrate Synthase',
        'ICL': 'Isocitrate Lyase',
        'MALS': 'Malate Synthase',
        'MDH': 'Malate Dehydrogenase'
    }
    
    print("\nAcetate 경로 반응 플럭스:")
    print(f"{'반응 ID':<20} {'이름':<30} {'플럭스':<15} {'상태':<10}")
    print("-" * 80)
    
    for rxn_id, rxn_name in acetate_pathway.items():
        try:
            flux = solution.fluxes.get(rxn_id, 0)
            status = "활성" if abs(flux) > 1e-8 else "비활성"
            print(f"{rxn_id:<20} {rxn_name:<30} {flux:<15.6f} {status:<10}")
        except KeyError:
            print(f"{rxn_id:<20} {rxn_name:<30} {'[NOT FOUND]':<15} {'':<10}")

def check_blocked_critical_reactions(model, solution):
    """핵심 반응 블록 상태 확인"""
    if not solution or solution.status != 'optimal':
        return
    
    print("\n" + "="*70)
    print("핵심 반응 플럭스 확인")
    print("="*70)
    
    critical_reactions = {
        'EX_ac_e': 'Acetate Exchange',
        'ACt': 'Acetate Transport',
        'ACt2rpp': 'Acetate Transport (periplasm)',
        'SUCOAACTr': 'Acetate-CoA Transferase',
        'ACS': 'Acetyl-CoA Synthetase',
        'CS': 'Citrate Synthase',
        'ICL': 'Isocitrate Lyase',
        'MALS': 'Malate Synthase',
        'SUCD': 'Succinate Dehydrogenase',
        'FUM': 'Fumarase',
        'MDH': 'Malate Dehydrogenase',
        'Growth': 'Biomass'
    }
    
    print(f"\n{'반응 ID':<20} {'이름':<30} {'플럭스':<15} {'블록됨':<10}")
    print("-" * 80)
    
    blocked_critical = []
    
    for rxn_id, rxn_name in critical_reactions.items():
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            flux = solution.fluxes.get(rxn_id, 0)
            is_blocked = abs(flux) < 1e-8
            blocked_str = "Yes" if is_blocked else "No"
            
            if is_blocked and rxn_id != 'Growth':
                blocked_critical.append(rxn_id)
            
            print(f"{rxn_id:<20} {rxn_name:<30} {flux:<15.6f} {blocked_str:<10}")
        except KeyError:
            print(f"{rxn_id:<20} {rxn_name:<30} {'[NOT FOUND]':<15} {'N/A':<10}")
    
    if blocked_critical:
        print(f"\n블록된 핵심 반응: {len(blocked_critical)}개")
        print(f"  {', '.join(blocked_critical)}")
    else:
        print("\n[OK] 모든 핵심 반응이 활성입니다")
    
    return blocked_critical

def main():
    print("="*70)
    print("제약 조건 단계별 진단")
    print("="*70)
    
    # 모델 로드
    model = load_model("BaseModel.xml")
    
    # Biomass 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass: {biomass_rxn.id}")
    
    # 단계별 테스트
    df_results, solution3, solution4 = test_growth_stepwise(model, biomass_rxn)
    
    # 결과 저장
    df_results.to_csv('constraint_stepwise_test.csv', index=False)
    print(f"\n[OK] 결과 저장: constraint_stepwise_test.csv")
    
    # Acetate 경로 플럭스 분석
    if solution3 and solution3.status == 'optimal':
        analyze_acetate_fluxes(model, solution3)
        blocked = check_blocked_critical_reactions(model, solution3)
    
    # 최종 요약
    print("\n" + "="*70)
    print("최종 요약")
    print("="*70)
    
    print("\n단계별 성장 가능 여부:")
    for _, row in df_results.iterrows():
        status_icon = "[OK]" if row['Can_Grow'] else "[FAIL]"
        print(f"  {status_icon} {row['Step']}: {row['Status']} (Biomass: {row['Biomass_flux']:.6f})")
    
    # 결론
    if df_results.iloc[0]['Can_Grow']:
        print("\n[결론] 모델 자체는 정상 작동합니다")
        
        if df_results.iloc[2]['Can_Grow']:
            print("  → Acetate만으로도 성장 가능합니다!")
        else:
            print("  → Acetate만으로는 성장 불가 (부트스트랩 문제)")
            
            if len(df_results) > 3 and df_results.iloc[3]['Can_Grow']:
                print("  → CoA 부트스트랩으로 해결 가능합니다")
            else:
                print("  → 추가 부트스트랩 또는 gap-filling 필요합니다")
    
    print("="*70)

if __name__ == "__main__":
    main()
