#!/usr/bin/env python
"""
GSM Quality Control 및 성장 테스트
Step 2: 성장 테스트(FBA) + 기본 QC
"""

import cobra
from cobra import Model
import pandas as pd
import numpy as np
from collections import defaultdict

def load_model(model_path="BaseModel.xml"):
    """모델 로드"""
    print(f"모델 로드 중: {model_path}")
    model = cobra.io.read_sbml_model(model_path)
    print(f"[OK] 모델 로드 완료: {model.id}\n")
    return model

def find_biomass_reaction(model):
    """Biomass 반응 찾기"""
    biomass_candidates = []
    
    # Objective로 설정된 반응 찾기
    if model.objective:
        # Objective의 변수들 확인
        try:
            objective_vars = list(model.objective.variables)
            for var in objective_vars:
                try:
                    rxn_id = var.name
                    rxn = model.reactions.get_by_id(rxn_id)
                    if rxn not in biomass_candidates:
                        biomass_candidates.append(rxn)
                except (KeyError, AttributeError):
                    continue
        except AttributeError:
            pass
    
    # 이름으로 찾기
    for rxn in model.reactions:
        if 'biomass' in rxn.id.lower() or 'growth' in rxn.id.lower():
            if rxn not in biomass_candidates:
                biomass_candidates.append(rxn)
    
    if biomass_candidates:
        return biomass_candidates[0]
    else:
        # 'BIOMASS' 또는 'Growth'로 찾기
        for name in ['BIOMASS', 'Growth', 'BIOMASS_Ecoli_core']:
            try:
                return model.reactions.get_by_id(name)
            except KeyError:
                continue
    
    return None

def setup_acetate_minimal_medium(model, biomass_rxn):
    """
    Acetate minimal medium 설정
    
    Acetate만 탄소원으로 허용
    질소원, 무기염, 산소 허용
    """
    print("="*70)
    print("2-1. Acetate Minimal Medium 성장 테스트")
    print("="*70)
    
    # 모든 exchange reaction을 0으로 설정 (기본적으로 차단)
    for rxn in model.exchanges:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    
    # Acetate 허용 (탄소원)
    try:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.lower_bound = -1000  # 무제한 uptake
        print(f"[OK] Acetate 교환 반응 설정: {ex_ac.id}")
    except KeyError:
        # 다른 이름으로 찾기
        for rxn in model.exchanges:
            if 'ac' in rxn.id.lower() and ('e' in rxn.id or '_e' in rxn.id):
                rxn.lower_bound = -1000
                print(f"[OK] Acetate 교환 반응 설정: {rxn.id}")
                break
    
    # 필수 무기염 및 질소원 허용
    essential_exchanges = {
        'EX_nh4_e': 'Ammonium (질소원)',
        'EX_h2o_e': 'Water',
        'EX_h_e': 'Proton',
        'EX_pi_e': 'Phosphate',
        'EX_so4_e': 'Sulfate',
        'EX_k_e': 'Potassium',
        'EX_na1_e': 'Sodium',
        'EX_mg2_e': 'Magnesium',
        'EX_ca2_e': 'Calcium',
        'EX_fe2_e': 'Iron',
        'EX_mn2_e': 'Manganese',
        'EX_zn2_e': 'Zinc',
        'EX_co2_e': 'CO2',
        'EX_o2_e': 'Oxygen'
    }
    
    print("\n필수 무기염/영양소 설정:")
    for ex_id, name in essential_exchanges.items():
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            if ex_id == 'EX_co2_e':
                ex_rxn.lower_bound = -1000  # CO2는 생성/소비 모두 가능
                ex_rxn.upper_bound = 1000
            elif ex_id == 'EX_o2_e':
                ex_rxn.lower_bound = -1000  # O2 uptake
                ex_rxn.upper_bound = 1000
            else:
                ex_rxn.lower_bound = -1000  # Uptake 허용
            print(f"  [OK] {name}: {ex_id}")
        except KeyError:
            # 다른 compartment 확인
            for rxn in model.exchanges:
                if ex_id.replace('_e', '') in rxn.id and ('_e' in rxn.id or '_p' in rxn.id):
                    rxn.lower_bound = -1000
                    print(f"  [OK] {name}: {rxn.id} (대체)")
                    break
    
    # Biomass objective 설정
    if biomass_rxn:
        model.objective = biomass_rxn.id
        print(f"\n[OK] Biomass objective 설정: {biomass_rxn.id}")
    else:
        print("\n[ERROR] Biomass 반응을 찾을 수 없습니다!")
        return None
    
    return model

def test_growth_on_acetate(model, biomass_rxn):
    """Acetate에서 성장 테스트"""
    print("\nFBA 최적화 수행 중...")
    
    try:
        solution = model.optimize()
        
        if solution.status == 'optimal':
            biomass_flux = solution.fluxes.get(biomass_rxn.id, 0)
            
            print(f"\n[결과]")
            print(f"  최적화 상태: {solution.status}")
            print(f"  Biomass flux: {biomass_flux:.6f} 1/h")
            print(f"  Objective value: {solution.objective_value:.6f}")
            
            if biomass_flux > 0:
                print(f"\n[OK] 성장 가능: Biomass flux > 0")
                
                # 성장률 평가 (일반적으로 0.01 ~ 2.0 1/h 범위)
                if 0.001 <= biomass_flux <= 2.0:
                    print(f"  [OK] 성장률이 합리적인 범위입니다 (0.001 ~ 2.0 1/h)")
                elif biomass_flux < 0.001:
                    print(f"  [WARNING] 성장률이 매우 낮습니다 (< 0.001 1/h)")
                else:
                    print(f"  [WARNING] 성장률이 비정상적으로 높습니다 (> 2.0 1/h)")
                
                # 주요 플럭스 확인
                print(f"\n주요 교환 플럭스:")
                for rxn in model.exchanges:
                    flux = solution.fluxes.get(rxn.id, 0)
                    if abs(flux) > 1e-6:
                        print(f"  {rxn.id}: {flux:.4f}")
                
                return True, biomass_flux
            else:
                print(f"\n[FAIL] 성장 불가: Biomass flux = 0")
                print(f"  → 모델에 갭이 있거나 경로가 완성되지 않았습니다.")
                return False, 0
                
        else:
            print(f"\n[FAIL] 최적화 실패: {solution.status}")
            return False, 0
            
    except Exception as e:
        print(f"\n[ERROR] FBA 최적화 중 오류: {e}")
        return False, 0

def test_energy_loop(model, biomass_rxn):
    """
    에너지 루프 체크: 무탄소 ATP 생성 테스트
    
    모든 탄소원 uptake = 0
    산소는 허용
    ATP 생성 또는 biomass 생성 확인
    """
    print("\n" + "="*70)
    print("2-2. 에너지 루프 체크 (무탄소 ATP 생성 테스트)")
    print("="*70)
    
    # 모든 exchange reaction 차단
    for rxn in model.exchanges:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    
    # 산소만 허용
    try:
        ex_o2 = model.reactions.get_by_id('EX_o2_e')
        ex_o2.lower_bound = -1000
        ex_o2.upper_bound = 1000
        print("[OK] 산소만 허용 설정")
    except KeyError:
        print("[WARNING] EX_o2_e를 찾을 수 없습니다")
    
    # Biomass objective
    if biomass_rxn:
        model.objective = biomass_rxn.id
    
    print("\nFBA 최적화 수행 (탄소원 없이)...")
    
    try:
        solution = model.optimize()
        
        if solution.status == 'optimal':
            biomass_flux = solution.fluxes.get(biomass_rxn.id, 0)
            
            # ATP 생성 확인
            atp_rxns = [rxn for rxn in model.reactions if 'atp' in rxn.id.lower() and 'EX' not in rxn.id]
            atp_production = False
            
            for rxn in atp_rxns:
                flux = solution.fluxes.get(rxn.id, 0)
                if abs(flux) > 1e-6:
                    atp_production = True
                    break
            
            print(f"\n[결과]")
            print(f"  최적화 상태: {solution.status}")
            print(f"  Biomass flux: {biomass_flux:.6f}")
            print(f"  ATP 생성 반응 플럭스: {'발견됨' if atp_production else '없음'}")
            
            if biomass_flux > 1e-6 or atp_production:
                print(f"\n[FAIL] 무탄소에서 ATP/Biomass 생성 발생!")
                print(f"  → 모델에 에너지 루프가 있습니다 (물리적으로 불가능)")
                print(f"  → 이 상태에서 AI/XAI/causal 분석은 신뢰할 수 없습니다.")
                return False
            else:
                print(f"\n[OK] 무탄소에서 ATP/Biomass 생성 없음")
                print(f"  → 에너지 루프 없음 (정상)")
                return True
                
        else:
            print(f"\n[OK] 최적화 불가능 (탄소원 없이 성장 불가)")
            return True
            
    except Exception as e:
        print(f"\n[ERROR] 테스트 중 오류: {e}")
        return False

def find_blocked_reactions(model):
    """Blocked reactions 찾기"""
    print("\n" + "="*70)
    print("2-3. Blocked Reactions 분석")
    print("="*70)
    
    print("Blocked reactions 분석 중... (시간이 걸릴 수 있습니다)")
    
    try:
        # 모든 exchange를 열어서 테스트
        for rxn in model.exchanges:
            rxn.lower_bound = -1000
            rxn.upper_bound = 1000
        
        blocked = cobra.flux_analysis.find_blocked_reactions(model, open_exchanges=True)
        
        print(f"\n[결과]")
        print(f"  전체 반응 수: {len(model.reactions)}")
        print(f"  Blocked reactions: {len(blocked)}")
        print(f"  활성 반응: {len(model.reactions) - len(blocked)}")
        print(f"  Blocked 비율: {len(blocked)/len(model.reactions)*100:.2f}%")
        
        if len(blocked) > 0:
            print(f"\nBlocked reactions 샘플 (최대 20개):")
            for i, rxn_id in enumerate(list(blocked)[:20]):
                try:
                    rxn = model.reactions.get_by_id(rxn_id)
                    print(f"  {i+1}. {rxn_id}: {rxn.name}")
                except:
                    print(f"  {i+1}. {rxn_id}")
        
        return blocked
        
    except Exception as e:
        print(f"[ERROR] Blocked reactions 분석 중 오류: {e}")
        return set()

def calculate_model_statistics(model):
    """모델 통계 계산"""
    stats = {
        'total_reactions': len(model.reactions),
        'total_metabolites': len(model.metabolites),
        'total_genes': len(model.genes),
        'compartments': len(model.compartments),
        'exchange_reactions': len(model.exchanges),
    }
    
    # GPR 통계
    reactions_with_gpr = 0
    reactions_without_gpr = 0
    
    for rxn in model.reactions:
        if len(rxn.genes) > 0:
            reactions_with_gpr += 1
        else:
            reactions_without_gpr += 1
    
    stats['reactions_with_gpr'] = reactions_with_gpr
    stats['reactions_without_gpr'] = reactions_without_gpr
    stats['gpr_coverage'] = reactions_with_gpr / len(model.reactions) * 100 if len(model.reactions) > 0 else 0
    
    # 가역성 통계
    reversible = sum(1 for rxn in model.reactions if rxn.reversibility)
    stats['reversible_reactions'] = reversible
    stats['irreversible_reactions'] = len(model.reactions) - reversible
    
    return stats

def print_qc_summary(model, biomass_flux, energy_loop_ok, blocked_reactions, stats):
    """QC 요약 출력"""
    print("\n" + "="*70)
    print("QC 요약 리포트")
    print("="*70)
    
    print("\n[모델 통계]")
    print(f"  반응 수: {stats['total_reactions']}")
    print(f"  대사물질 수: {stats['total_metabolites']}")
    print(f"  유전자 수: {stats['total_genes']}")
    print(f"  구획 수: {stats['compartments']}")
    print(f"  교환 반응 수: {stats['exchange_reactions']}")
    print(f"  GPR 연결 반응: {stats['reactions_with_gpr']} ({stats['gpr_coverage']:.1f}%)")
    print(f"  가역 반응: {stats['reversible_reactions']}")
    print(f"  비가역 반응: {stats['irreversible_reactions']}")
    print(f"  Blocked reactions: {len(blocked_reactions)} ({len(blocked_reactions)/stats['total_reactions']*100:.1f}%)")
    
    print("\n[성장 테스트 결과]")
    if biomass_flux > 0:
        print(f"  Acetate minimal medium: [PASS] 성장 가능 (flux = {biomass_flux:.6f} 1/h)")
    else:
        print(f"  Acetate minimal medium: [FAIL] 성장 불가")
    
    print(f"  에너지 루프 체크: {'[PASS]' if energy_loop_ok else '[FAIL]'} {'정상' if energy_loop_ok else '에너지 루프 발견'}")
    
    print("\n[Paper B 통과 기준]")
    if biomass_flux > 0 and energy_loop_ok:
        print("  [OK] Paper B 통과 가능: 기본 QC 통과")
    else:
        print("  [FAIL] Paper B 통과 불가: 모델 수정 필요")
        if biomass_flux <= 0:
            print("    - Acetate에서 성장 불가 → 갭필/경로 확인 필요")
        if not energy_loop_ok:
            print("    - 에너지 루프 발견 → 반응 방향성/스토이키오메트리 확인 필요")
    
    print("\n" + "="*70)

def main():
    """메인 함수"""
    # 모델 로드
    model = load_model("BaseModel.xml")
    
    # Biomass 반응 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응을 찾을 수 없습니다!")
        return
    
    print(f"[OK] Biomass 반응 발견: {biomass_rxn.id}\n")
    
    # 2-1. Acetate minimal medium 성장 테스트
    model_acetate = setup_acetate_minimal_medium(model.copy(), biomass_rxn)
    growth_ok, biomass_flux = test_growth_on_acetate(model_acetate, biomass_rxn)
    
    # 2-2. 에너지 루프 체크
    energy_loop_ok = test_energy_loop(model.copy(), biomass_rxn)
    
    # 2-3. Blocked reactions
    blocked_reactions = find_blocked_reactions(model)
    
    # 통계 계산
    stats = calculate_model_statistics(model)
    
    # QC 요약
    print_qc_summary(model, biomass_flux, energy_loop_ok, blocked_reactions, stats)
    
    # 결과를 파일로 저장
    results = {
        'biomass_flux': biomass_flux,
        'growth_ok': growth_ok,
        'energy_loop_ok': energy_loop_ok,
        'blocked_count': len(blocked_reactions),
        **stats
    }
    
    df_results = pd.DataFrame([results])
    df_results.to_csv('qc_results.csv', index=False)
    print("\n[OK] 결과가 qc_results.csv에 저장되었습니다.")

if __name__ == "__main__":
    main()

