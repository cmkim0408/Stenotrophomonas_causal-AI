#!/usr/bin/env python
"""
Biomass 반응 검증
Biomass 반응 구성 요소 생산 가능 여부 확인
"""

import cobra
import pandas as pd

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    print(f"[OK] 모델 로드: {model.id}")
    return model

def find_biomass_reaction(model):
    """Biomass 반응 찾기"""
    for name in ['Growth', 'BIOMASS', 'BIOMASS_Ecoli_core']:
        try:
            return model.reactions.get_by_id(name)
        except KeyError:
            continue
    return None

def check_biomass_components(model, biomass_rxn):
    """Biomass 구성 요소 확인"""
    print("="*70)
    print("Biomass 반응 구성 요소 분석")
    print("="*70)
    
    print(f"\nBiomass 반응: {biomass_rxn.id}")
    print(f"이름: {biomass_rxn.name}")
    print(f"\n반응식:")
    print(f"  {biomass_rxn.reaction}")
    
    # 구성 요소 분석
    metabolites = biomass_rxn.metabolites
    
    print(f"\n총 구성 요소: {len(metabolites)}개")
    
    # 절대값 기준으로 정렬
    sorted_mets = sorted(metabolites.items(), key=lambda x: abs(x[1]), reverse=True)
    
    print(f"\n주요 구성 요소 (계수 상위 15개):")
    print(f"{'대사물질 ID':<25} {'계수':<15} {'상태':<15}")
    print("-" * 60)
    
    component_status = []
    
    for met, coeff in sorted_mets[:15]:
        # 생산 가능 여부 테스트
        can_produce = test_metabolite_production(model, met.id)
        status = "[생산 가능]" if can_produce else "[생산 불가]"
        
        component_status.append({
            'Metabolite_ID': met.id,
            'Coefficient': coeff,
            'Can_Produce': can_produce,
            'Status': status
        })
        
        status_display = "OK" if can_produce else "FAIL"
        print(f"{met.id:<25} {coeff:<15.6f} {status_display:<15}")
    
    return component_status

def test_metabolite_production(model, metabolite_id):
    """특정 대사물질 생산 가능 여부 테스트"""
    met = None
    try:
        met = model.metabolites.get_by_id(metabolite_id)
    except KeyError:
        # 다른 구획 확인
        for comp in ['_c', '_e', '_p']:
            try:
                met_id_alt = metabolite_id.replace('_c', comp).replace('_e', comp).replace('_p', comp)
                if met_id_alt != metabolite_id:
                    try:
                        met = model.metabolites.get_by_id(met_id_alt)
                        metabolite_id = met_id_alt
                        break
                    except KeyError:
                        continue
            except:
                continue
    
    if met is None:
        return False
    
    # 모든 exchange 차단 (upper_bound 먼저 설정)
    original_bounds = {}
    for rxn in model.exchanges:
        original_bounds[rxn.id] = (rxn.lower_bound, rxn.upper_bound)
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    # 필수 영양소만 허용
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
    
    # 포도당 허용 (일반적인 테스트 탄소원)
    try:
        ex_glc = model.reactions.get_by_id('EX_glc__D_e')
        ex_glc.lower_bound = -100
        ex_glc.upper_bound = 1000
    except KeyError:
        try:
            for rxn in model.exchanges:
                if 'glc' in rxn.id.lower() or 'glucose' in rxn.id.lower():
                    rxn.lower_bound = -100
                    rxn.upper_bound = 1000
                    break
        except:
            pass
    
    # 임시 생산 반응 생성
    try:
        test_rxn = cobra.Reaction(f'TEST_{met.id}')
        test_rxn.add_metabolites({met: 1})
        test_rxn.lower_bound = 0
        test_rxn.upper_bound = 1000
        
        model.add_reactions([test_rxn])
        model.objective = test_rxn.id
        
        solution = model.optimize()
        
        # 테스트 반응 제거
        model.remove_reactions([test_rxn])
        
        # 경계 조건 복원
        for rxn_id, (lb, ub) in original_bounds.items():
            try:
                rxn = model.reactions.get_by_id(rxn_id)
                rxn.lower_bound = lb
                rxn.upper_bound = ub
            except KeyError:
                pass
        
        return solution.status == 'optimal' and solution.objective_value > 1e-6
        
    except Exception as e:
        # 경계 조건 복원
        for rxn_id, (lb, ub) in original_bounds.items():
            try:
                rxn = model.reactions.get_by_id(rxn_id)
                rxn.lower_bound = lb
                rxn.upper_bound = ub
            except KeyError:
                pass
        return False

def analyze_missing_components(model, biomass_rxn):
    """생산 불가능한 구성 요소 분석"""
    print("\n" + "="*70)
    print("생산 불가능한 구성 요소 분석")
    print("="*70)
    
    metabolites = biomass_rxn.metabolites
    missing_components = []
    
    print("\n생산 불가능한 구성 요소 검색 중...")
    
    for met, coeff in metabolites.items():
        can_produce = test_metabolite_production(model, met.id)
        if not can_produce:
            missing_components.append({
                'Metabolite_ID': met.id,
                'Coefficient': coeff,
                'Metabolite_Name': getattr(met, 'name', 'N/A')
            })
    
    if missing_components:
        print(f"\n총 {len(missing_components)}개 구성 요소가 생산 불가능:")
        print(f"{'대사물질 ID':<25} {'계수':<15} {'이름':<30}")
        print("-" * 75)
        
        for comp in sorted(missing_components, key=lambda x: abs(x['Coefficient']), reverse=True):
            print(f"{comp['Metabolite_ID']:<25} {comp['Coefficient']:<15.6f} {comp['Metabolite_Name']:<30}")
        
        # 가장 중요한 구성 요소 확인 (절대값 기준)
        most_important = sorted(missing_components, key=lambda x: abs(x['Coefficient']), reverse=True)[:10]
        print(f"\n중요도 상위 10개 생산 불가능 구성 요소:")
        for comp in most_important:
            print(f"  - {comp['Metabolite_ID']} (계수: {comp['Coefficient']:.6f})")
    else:
        print("\n[OK] 모든 구성 요소가 생산 가능합니다!")
    
    return missing_components

def main():
    print("="*70)
    print("1단계: Biomass 반응 검증")
    print("="*70)
    
    # 모델 로드
    model = load_model("BaseModel.xml")
    
    # Biomass 반응 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass 반응: {biomass_rxn.id}")
    
    # 구성 요소 확인
    component_status = check_biomass_components(model, biomass_rxn)
    
    # 생산 불가능한 구성 요소 분석
    missing_components = analyze_missing_components(model, biomass_rxn)
    
    # 결과 저장
    if component_status:
        df_status = pd.DataFrame(component_status)
        df_status.to_csv('biomass_component_status.csv', index=False)
        print(f"\n[OK] 결과 저장: biomass_component_status.csv")
    
    if missing_components:
        df_missing = pd.DataFrame(missing_components)
        df_missing.to_csv('biomass_missing_components.csv', index=False)
        print(f"[OK] 결과 저장: biomass_missing_components.csv")
    
    print("\n" + "="*70)
    print("Biomass 반응 검증 완료")
    print("="*70)
    
    if missing_components:
        print(f"\n[결론] {len(missing_components)}개 구성 요소가 생산 불가능")
        print("  → 이들이 Biomass 생성을 막고 있을 수 있습니다")
    else:
        print("\n[결론] 모든 구성 요소가 생산 가능")
        print("  → 다른 원인을 찾아야 합니다")

if __name__ == "__main__":
    main()
