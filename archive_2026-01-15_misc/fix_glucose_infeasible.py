#!/usr/bin/env python
"""
포도당 infeasible 문제 해결
포도당 경로와 필수 대사물질 생산 경로 확인
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

def test_glucose_medium_detailed(model, biomass_rxn):
    """포도당 medium 상세 테스트"""
    print("="*70)
    print("포도당 Medium 상세 진단")
    print("="*70)
    
    # 모든 exchange 초기화
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    # 포도당 찾기 및 설정
    glucose_found = False
    glucose_exchanges = ['EX_glc__D_e', 'EX_glc_e', 'EX_glucose_e']
    
    for ex_id in glucose_exchanges:
        try:
            ex_glc = model.reactions.get_by_id(ex_id)
            ex_glc.lower_bound = -100
            ex_glc.upper_bound = 1000
            print(f"[OK] 포도당 Exchange: {ex_id}")
            glucose_found = True
            break
        except KeyError:
            continue
    
    if not glucose_found:
        # 패턴으로 찾기
        for rxn in model.exchanges:
            if 'glc' in rxn.id.lower() or 'glucose' in rxn.id.lower():
                rxn.lower_bound = -100
                rxn.upper_bound = 1000
                print(f"[OK] 포도당 Exchange: {rxn.id}")
                glucose_found = True
                break
    
    if not glucose_found:
        print("[ERROR] 포도당 exchange 반응을 찾을 수 없습니다!")
        return None
    
    # 필수 영양소 단계별 추가
    print("\n필수 영양소 단계별 추가 테스트:")
    
    # 기본 영양소
    basic_nutrients = {
        'EX_nh4_e': 'Ammonium (질소원)',
        'EX_h2o_e': 'Water',
        'EX_h_e': 'Proton',
        'EX_pi_e': 'Phosphate',
        'EX_o2_e': 'Oxygen',
        'EX_co2_e': 'CO2'
    }
    
    print("\n[단계 1] 기본 영양소만 추가")
    for ex_id, name in basic_nutrients.items():
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            if ex_id in ['EX_co2_e', 'EX_o2_e']:
                ex_rxn.lower_bound = -1000
                ex_rxn.upper_bound = 1000
            else:
                ex_rxn.lower_bound = -1000
            print(f"  [OK] {name}: {ex_id}")
        except KeyError:
            print(f"  [MISSING] {name}: {ex_id}")
    
    model.objective = biomass_rxn.id
    solution1 = model.optimize()
    
    print(f"\n  최적화 상태: {solution1.status}")
    if solution1.status == 'optimal':
        biomass_flux = solution1.fluxes.get(biomass_rxn.id, 0)
        print(f"  Biomass flux: {biomass_flux:.6f} 1/h")
        if biomass_flux > 1e-6:
            print("  [SUCCESS] 기본 영양소만으로도 성장 가능!")
            return solution1
    
    # 무기염 추가
    print("\n[단계 2] 무기염 추가")
    minerals = {
        'EX_so4_e': 'Sulfate',
        'EX_k_e': 'Potassium',
        'EX_na1_e': 'Sodium',
        'EX_mg2_e': 'Magnesium',
        'EX_ca2_e': 'Calcium',
        'EX_fe2_e': 'Iron',
        'EX_mn2_e': 'Manganese',
        'EX_zn2_e': 'Zinc'
    }
    
    for ex_id, name in minerals.items():
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            ex_rxn.lower_bound = -1000
            print(f"  [OK] {name}: {ex_id}")
        except KeyError:
            print(f"  [MISSING] {name}: {ex_id}")
    
    solution2 = model.optimize()
    
    print(f"\n  최적화 상태: {solution2.status}")
    if solution2.status == 'optimal':
        biomass_flux = solution2.fluxes.get(biomass_rxn.id, 0)
        print(f"  Biomass flux: {biomass_flux:.6f} 1/h")
        if biomass_flux > 1e-6:
            print("  [SUCCESS] 무기염 추가 후 성장 가능!")
            return solution2
    else:
        print(f"  [FAIL] 여전히 {solution2.status}")
        return solution2
    
    return solution2

def check_biomass_requirements(model, biomass_rxn):
    """Biomass 반응 요구사항 확인"""
    print("\n" + "="*70)
    print("Biomass 반응 요구사항 확인")
    print("="*70)
    
    metabolites = biomass_rxn.metabolites
    
    # 생성물 확인 (양수 계수 = 생성물)
    products = [(met, coeff) for met, coeff in metabolites.items() if coeff > 0]
    # 반응물 확인 (음수 계수 = 반응물)
    reactants = [(met, abs(coeff)) for met, coeff in metabolites.items() if coeff < 0]
    
    print(f"\n반응물: {len(reactants)}개")
    print(f"생성물: {len(products)}개")
    
    # 반응물 중에서 exchange가 필요한 것 확인
    print("\n필수 exchange 반응 확인 (절대값 상위 20개):")
    sorted_reactants = sorted(reactants, key=lambda x: x[1], reverse=True)
    
    for met, coeff in sorted_reactants[:20]:
        met_id_base = met.id.replace('_c', '').replace('_e', '').replace('_p', '')
        
        # Exchange 반응 찾기
        exchange_found = False
        exchange_ids = [f'EX_{met.id}', f'EX_{met_id_base}_e', f'EX_{met.id.replace("_c", "_e")}']
        
        for ex_id in exchange_ids:
            try:
                ex_rxn = model.reactions.get_by_id(ex_id)
                exchange_found = True
                print(f"  [OK] {met.id} (계수: {coeff:.6f}): {ex_id} 존재")
                break
            except KeyError:
                continue
        
        if not exchange_found:
            # 패턴으로 찾기
            for rxn in model.exchanges:
                if met_id_base in rxn.id.lower() or met.id.split('_')[0] in rxn.id.lower():
                    print(f"  [OK] {met.id} (계수: {coeff:.6f}): {rxn.id} (패턴 매칭)")
                    exchange_found = True
                    break
            
            if not exchange_found:
                print(f"  [WARNING] {met.id} (계수: {coeff:.6f}): Exchange 없음 (생산 필요)")

def find_missing_exchanges(model, biomass_rxn):
    """Biomass에 필요한 누락된 exchange 확인"""
    print("\n" + "="*70)
    print("누락된 Exchange 반응 확인")
    print("="*70)
    
    metabolites = biomass_rxn.metabolites
    reactants = [(met, abs(coeff)) for met, coeff in metabolites.items() if coeff < 0]
    sorted_reactants = sorted(reactants, key=lambda x: x[1], reverse=True)
    
    missing_exchanges = []
    
    # 중요한 대사물질 중 exchange가 없는 것
    important_missing = []
    
    for met, coeff in sorted_reactants[:30]:  # 상위 30개 확인
        # 이미 존재하는 exchange 확인
        has_exchange = False
        
        # 직접 확인
        ex_ids = [f'EX_{met.id}', f'EX_{met.id.replace("_c", "_e")}']
        for ex_id in ex_ids:
            try:
                model.reactions.get_by_id(ex_id)
                has_exchange = True
                break
            except KeyError:
                continue
        
        # 패턴으로 확인
        if not has_exchange:
            met_base = met.id.split('_')[0]
            for rxn in model.exchanges:
                if met_base in rxn.id.lower():
                    has_exchange = True
                    break
        
        if not has_exchange and abs(coeff) > 0.001:  # 의미있는 계수만
            missing_exchanges.append({
                'Metabolite_ID': met.id,
                'Coefficient': coeff,
                'Name': getattr(met, 'name', 'N/A')
            })
            
            if abs(coeff) > 0.01:  # 중요도가 높은 것
                important_missing.append(met.id)
    
    if missing_exchanges:
        print(f"\n중요한 누락 Exchange (계수 > 0.01): {len(important_missing)}개")
        for met_id in important_missing[:10]:
            met_info = next((m for m in missing_exchanges if m['Metabolite_ID'] == met_id), None)
            if met_info:
                print(f"  - {met_id}: {met_info['Name']} (계수: {met_info['Coefficient']:.6f})")
    
    return missing_exchanges, important_missing

def main():
    print("="*70)
    print("포도당 Infeasible 문제 해결")
    print("="*70)
    
    # 모델 로드
    model = load_model("BaseModel.xml")
    
    # Biomass 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass: {biomass_rxn.id}")
    
    # 포도당 medium 상세 테스트
    solution = test_glucose_medium_detailed(model, biomass_rxn)
    
    # Biomass 요구사항 확인
    check_biomass_requirements(model, biomass_rxn)
    
    # 누락된 exchange 확인
    missing_exchanges, important_missing = find_missing_exchanges(model, biomass_rxn)
    
    # 결과 저장
    if missing_exchanges:
        df_missing = pd.DataFrame(missing_exchanges)
        df_missing.to_csv('missing_exchanges.csv', index=False)
        print(f"\n[OK] 결과 저장: missing_exchanges.csv")
    
    print("\n" + "="*70)
    print("진단 완료")
    print("="*70)
    
    if solution and solution.status == 'optimal' and solution.objective_value > 1e-6:
        print("\n[결론] 포도당으로 성장 가능합니다!")
    else:
        print("\n[결론] 포도당으로도 성장 불가")
        if important_missing:
            print(f"  → {len(important_missing)}개 중요한 대사물질의 Exchange 누락 가능")
        print("  → 추가 영양소 또는 경로가 필요할 수 있습니다")
    
    print("="*70)

if __name__ == "__main__":
    main()
