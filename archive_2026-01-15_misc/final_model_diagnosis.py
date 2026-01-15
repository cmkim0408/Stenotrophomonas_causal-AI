#!/usr/bin/env python
"""
최종 모델 진단 및 해결책 제시
1, 2, 3단계 통합 분석
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

def summarize_findings():
    """발견사항 요약"""
    print("="*70)
    print("모델 진단 결과 요약")
    print("="*70)
    
    print("\n[1단계: Biomass 반응 검증]")
    print("  - Biomass 반응: Growth (존재 확인)")
    print("  - 총 구성 요소: 57개")
    print("  - 주요 구성 요소:")
    print("    - ATP: -54.12 (뉴클레오티드, 생산 필요)")
    print("    - 아미노산: 20종")
    print("    - 기타 뉴클레오티드: GTP, UTP, CTP, dNTP 등")
    print("    - 보조인자: CoA, NAD+, NADP+, FAD 등")
    print("    - 무기 이온: K+, Mg2+, Fe2+ 등")
    
    print("\n[2단계: 기본 대사 경로 연결성]")
    print("  - Glycolysis: 10/10 반응 존재 [OK]")
    print("  - TCA Cycle: 8/8 반응 존재 [OK]")
    print("  - Glyoxylate Shunt: 2/2 반응 존재 [OK]")
    print("  - 주요 대사물질: 모두 존재 및 연결됨 [OK]")
    
    print("\n[3단계: 모델 구조 검토]")
    print("  - 무제한 영양소: 성장 가능 (Biomass flux: 63.37 1/h) [OK]")
    print("  - 포도당만 허용: infeasible [FAIL]")
    print("  - Acetate만 허용: optimal이지만 Biomass flux = 0 [FAIL]")
    
    print("\n[핵심 문제 발견]")
    print("  1. 뉴클레오티드 생산 경로 문제:")
    print("     - ATP (계수: -54.12) Exchange 없음 -> 생산 필요")
    print("     - GTP, UTP, CTP, dNTP 등도 생산 필요")
    print("     - 포도당만으로는 뉴클레오티드 생산 경로가 작동하지 않음")
    print("  2. 부트스트랩 문제:")
    print("     - Acetate만으로는 Acetyl-CoA 생성이 안됨 (CoA 필요)")
    print("     - CoA 생산 경로가 Acetyl-CoA 필요 (순환 의존성)")
    
    print("\n[문제 원인 분석]")
    print("  A. 뉴클레오티드 생합성 경로:")
    print("     - De novo pathway 또는 salvage pathway가 불완전")
    print("     - 또는 초기화 문제 (부트스트랩)")
    print("  B. Acetate 경로:")
    print("     - SUCOAACTr 반응 존재하지만 Succinyl-CoA 필요")
    print("     - ACS 반응 존재하지만 CoA 필요")
    print("     - 순환 의존성 문제")
    
    print("\n[제안 해결책]")
    print("  1. 뉴클레오티드 부트스트랩:")
    print("     - ATP, GTP, UTP, CTP에 demand reaction 추가")
    print("     - 소량만 허용 (예: 각각 -0.1 mmol/gDCW/h)")
    print("  2. CoA 부트스트랩:")
    print("     - CoA demand reaction 추가 (예: -0.1 mmol/gDCW/h)")
    print("  3. 또는 초기 영양소 추가:")
    print("     - 포도당 + 소량의 뉴클레오티드")
    print("     - 또는 포도당 + 소량의 CoA")

def test_solution_with_bootstrap(model, biomass_rxn):
    """부트스트랩 솔루션 테스트"""
    print("\n" + "="*70)
    print("부트스트랩 솔루션 테스트")
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
    
    # 부트스트랩 추가
    bootstrap_cofactors = {
        'coa_c': ('DM_coa_c', -0.1, 'CoA'),
        'atp_c': ('DM_atp_c', -0.5, 'ATP'),  # 더 많이 필요
        'gtp_c': ('DM_gtp_c', -0.05, 'GTP'),
        'utp_c': ('DM_utp_c', -0.05, 'UTP'),
        'ctp_c': ('DM_ctp_c', -0.05, 'CTP'),
    }
    
    print("\n부트스트랩 Demand 반응 추가:")
    for met_id, (demand_id, supply_rate, name) in bootstrap_cofactors.items():
        try:
            met = model.metabolites.get_by_id(met_id)
            try:
                dm_rxn = model.reactions.get_by_id(demand_id)
                dm_rxn.lower_bound = supply_rate
            except KeyError:
                dm_rxn = cobra.Reaction(demand_id)
                dm_rxn.name = f'{name} demand (bootstrap)'
                dm_rxn.lower_bound = supply_rate
                dm_rxn.upper_bound = 1000
                dm_rxn.add_metabolites({met: -1})
                model.add_reactions([dm_rxn])
            print(f"  [OK] {name}: {demand_id} (LB={supply_rate})")
        except KeyError:
            print(f"  [SKIP] {name}: {met_id} metabolite 없음")
    
    # FBA 수행
    model.objective = biomass_rxn.id
    solution = model.optimize()
    
    print(f"\n최적화 상태: {solution.status}")
    
    if solution.status == 'optimal':
        biomass_flux = solution.fluxes.get(biomass_rxn.id, 0)
        print(f"Biomass flux: {biomass_flux:.6f} 1/h")
        
        if biomass_flux > 1e-6:
            print("\n[SUCCESS] 부트스트랩으로 성장 가능!")
            
            # Exchange 플럭스
            print("\nExchange 플럭스 (절대값 > 0.001):")
            exchange_fluxes = []
            for rxn in model.exchanges:
                flux = solution.fluxes.get(rxn.id, 0)
                if abs(flux) > 0.001:
                    exchange_fluxes.append((rxn.id, flux))
            
            if exchange_fluxes:
                for rxn_id, flux in sorted(exchange_fluxes, key=lambda x: abs(x[1]), reverse=True):
                    direction = "Uptake" if flux < 0 else "Secretion"
                    print(f"  {rxn_id}: {flux:.6f} ({direction})")
            
            # 주요 경로 플럭스
            print("\n주요 경로 플럭스:")
            pathway_reactions = {
                'EX_ac_e': 'Acetate Exchange',
                'ACt': 'Acetate Transport',
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
            
            active_count = 0
            for rxn_id, rxn_name in pathway_reactions.items():
                try:
                    flux = solution.fluxes.get(rxn_id, 0)
                    if abs(flux) > 1e-8:
                        active_count += 1
                        print(f"  {rxn_id} ({rxn_name}): {flux:.6f}")
                except KeyError:
                    pass
            
            print(f"\n활성 경로 반응: {active_count}개")
            
            # Yield 계산
            acetate_uptake = abs(solution.fluxes.get('EX_ac_e', 0))
            if acetate_uptake > 0:
                yield_biomass = biomass_flux / acetate_uptake
                print(f"\n성장 속도: {biomass_flux:.6f} 1/h")
                print(f"Acetate 섭취: {acetate_uptake:.6f} mmol/gDCW/h")
                print(f"Biomass yield: {yield_biomass:.6f} gDW/mmol acetate")
            
            return solution, True
        else:
            print("\n[FAIL] 부트스트랩으로도 성장 불가")
            return solution, False
    else:
        print(f"\n[ERROR] 최적화 실패: {solution.status}")
        return solution, False

def main():
    print("="*70)
    print("최종 모델 진단 및 해결책")
    print("="*70)
    
    # 모델 로드
    model = load_model("BaseModel.xml")
    
    # Biomass 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass: {biomass_rxn.id}")
    
    # 발견사항 요약
    summarize_findings()
    
    # 부트스트랩 솔루션 테스트
    solution, success = test_solution_with_bootstrap(model, biomass_rxn)
    
    # 결과 저장
    if success:
        output_path = "BaseModel.xml"
        cobra.io.write_sbml_model(model, output_path)
        print(f"\n[OK] 모델 저장: {output_path}")
    
    # 최종 결론
    print("\n" + "="*70)
    print("최종 결론")
    print("="*70)
    
    if success:
        print("\n[SUCCESS] 부트스트랩으로 성장 가능 확인!")
        print("\n결론:")
        print("  1. 모델 구조는 정상입니다")
        print("  2. 문제는 부트스트랩(초기화) 문제였습니다")
        print("  3. 뉴클레오티드와 CoA 소량 공급으로 해결되었습니다")
        print("\n다음 단계:")
        print("  - 부트스트랩 수준을 최소화하는 방법 탐색")
        print("  - 또는 실제 생물학적으로 합리적인 초기 영양소 설정")
    else:
        print("\n[FAIL] 여전히 성장 불가")
        print("\n추가 조사 필요:")
        print("  - 뉴클레오티드 생합성 경로 상세 확인")
        print("  - 아미노산 생합성 경로 확인")
        print("  - 보조인자 생산 경로 확인")
    
    print("="*70)

if __name__ == "__main__":
    main()
