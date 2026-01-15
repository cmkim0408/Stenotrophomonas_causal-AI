#!/usr/bin/env python
"""
성장 실패 원인 종합 리포트
"""

import cobra
import pandas as pd

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    return model

def generate_comprehensive_report():
    """종합 리포트 생성"""
    print("="*70)
    print("성장 실패 원인 종합 리포트")
    print("="*70)
    
    model = load_model("BaseModel.xml")
    
    print("\n[1. 확인된 사실]")
    print("  ✓ Exchange reaction (EX_ac_e) 추가 완료")
    print("  ✓ Transport reaction (ACt) 추가 완료")
    print("  ✓ Acetate 경로 (ACS, TCA, Glyoxylate) 반응 존재")
    print("  ✓ Biomass reaction (Growth) 존재")
    print("  ✓ 에너지 루프 없음 (정상)")
    
    print("\n[2. 발견된 문제들]")
    
    print("\n  2-1. ATP 생산 불가")
    print("    - 원인: Acetate만으로는 ATP를 생성할 수 없음")
    print("    - 영향: 모든 에너지 의존 반응 (ACS 포함) 차단")
    print("    - 증상: ATP 생산 테스트 실패 (objective = 0)")
    
    print("\n  2-2. CoA 생산 불가")
    print("    - 원인: Acetyl-CoA가 없어서 CoA를 재생할 수 없음")
    print("    - 영향: ACS 반응 작동 불가")
    print("    - 증상: CoA 생산 테스트 실패 (objective = 0)")
    
    print("\n  2-3. OAA Bootstrap 문제")
    print("    - 원인: 초기 OAA 생성 경로 없음")
    print("    - 순환 의존성:")
    print("      * CS는 OAA 필요")
    print("      * OAA는 MDH (Malate → OAA)로 생성")
    print("      * Malate는 MALS (Glyoxylate shunt)로 생성")
    print("      * Glyoxylate는 ICL (Isocitrate → Glyoxylate)로 생성")
    print("      * Isocitrate는 ACONT (Citrate → Isocitrate)로 생성")
    print("      * Citrate는 CS (Acetyl-CoA + OAA → Citrate)로 생성")
    print("      → 순환! 초기 OAA 없이는 시작 불가")
    
    print("\n  2-4. PEP 생산 불가")
    print("    - 원인: ATP가 없어서 PEP 생산 경로 작동 불가")
    print("    - 영향: PPC (PEP → OAA) 경로 사용 불가")
    
    print("\n[3. 근본 원인: Bootstrap Problem]")
    print("  모델이 '닭과 달걀' 문제에 직면:")
    print("    - ATP가 필요하지만 ATP 생산에 Acetyl-CoA 필요")
    print("    - Acetyl-CoA 생산에 ATP 필요 (ACS 반응)")
    print("    - 순환 의존성으로 인해 경로 시작 불가")
    
    print("\n[4. 해결 방안]")
    print("\n  방안 1: 초기 Bootstrap 대사물질 제공")
    print("    - 소량의 ATP 제공 (예: 0.1 mmol/gDW/h)")
    print("    - 또는 소량의 CoA 제공")
    print("    - 또는 소량의 OAA 제공")
    print("    → 이렇게 하면 경로가 시작될 수 있음")
    
    print("\n  방안 2: 대체 경로 추가")
    print("    - ATP 생산을 위한 다른 경로 (예: 광합성, 발효)")
    print("    - OAA 초기화 경로 (예: Pyruvate carboxylase 활성화)")
    print("    - CoA 초기화 경로")
    
    print("\n  방안 3: 반응 방향성 조정")
    print("    - 일부 반응의 가역성 확인")
    print("    - 경계 조건 조정")
    
    print("\n[5. 권장 조치사항]")
    print("  1. Bootstrap 테스트 수행")
    print("     - 소량의 ATP/CoA/OAA를 제공하고 성장 가능 여부 확인")
    print("  2. Gap analysis 상세 수행")
    print("     - 누락된 반응 식별")
    print("  3. 모델 수정")
    print("     - 필요한 반응 추가")
    print("     - 반응 방향성 수정")
    
    print("\n[6. 다음 단계]")
    print("  → Bootstrap 대사물질을 제공한 상태에서 FBA 재실행")
    print("  → 성장 가능하면: 모델이 정상, 단지 초기화가 필요함")
    print("  → 성장 불가하면: 추가 gap filling 필요")

def test_with_bootstrap(model):
    """Bootstrap 테스트"""
    print("\n" + "="*70)
    print("Bootstrap 테스트: 소량의 ATP 제공")
    print("="*70)
    
    # Medium 설정
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
    
    # 소량의 ATP 제공 (bootstrap)
    try:
        from cobra import Reaction, Metabolite
        atp_c = model.metabolites.get_by_id('atp_c')
        
        # ATP 제공 반응 생성
        atp_bootstrap = Reaction('EX_atp_bootstrap')
        atp_bootstrap.name = 'ATP bootstrap (minimal)'
        atp_bootstrap.lower_bound = -0.1  # 소량만 제공
        atp_bootstrap.upper_bound = 0
        atp_bootstrap.add_metabolites({atp_c: -1})
        model.add_reactions([atp_bootstrap])
        
        print("[OK] ATP bootstrap 반응 추가: -0.1 mmol/gDW/h")
        
        # Biomass objective
        model.objective = 'Growth'
        
        # FBA 실행
        solution = model.optimize()
        
        print(f"\nBootstrap FBA 결과:")
        print(f"  상태: {solution.status}")
        print(f"  Biomass flux: {solution.objective_value:.6f}")
        
        if solution.status == 'optimal' and solution.objective_value > 1e-6:
            print(f"\n  [SUCCESS] Bootstrap으로 성장 가능!")
            print(f"    → 모델은 정상이며, 초기화만 필요합니다")
            
            print(f"\n  주요 플럭스:")
            key_rxns = ['EX_ac_e', 'ACt', 'ACS', 'CS', 'ICL', 'MALS', 'Growth']
            for rxn_id in key_rxns:
                flux = solution.fluxes.get(rxn_id, 0)
                if abs(flux) > 1e-6:
                    print(f"    {rxn_id}: {flux:.4f}")
        else:
            print(f"\n  [FAIL] Bootstrap으로도 성장 불가")
            print(f"    → 추가 gap filling 필요")
        
        # Bootstrap 반응 제거
        model.remove_reactions([atp_bootstrap])
        
    except Exception as e:
        print(f"[ERROR] Bootstrap 테스트 실패: {e}")
        import traceback
        traceback.print_exc()

def main():
    model = load_model("BaseModel.xml")
    
    generate_comprehensive_report()
    test_with_bootstrap(model)
    
    print("\n" + "="*70)
    print("종합 리포트 완료")
    print("="*70)

if __name__ == "__main__":
    main()

