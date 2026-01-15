#!/usr/bin/env python
"""
ATPM=0일 때 성장 불가 원인 진단
"""

import cobra
from pathlib import Path

def load_model(model_path):
    model = cobra.io.read_sbml_model(model_path)
    return model

def setup_acetate_medium(model):
    """Acetate 미디어 설정"""
    # 모든 exchange 차단
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    # Acetate 허용
    ex_ac = model.reactions.get_by_id('EX_ac_e')
    ex_ac.upper_bound = 1000
    ex_ac.lower_bound = -1000
    
    # 필수 무기염
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

def diagnose_atpm0_issue(model):
    """ATPM=0일 때 문제 진단"""
    print("="*80)
    print("ATPM=0일 때 성장 불가 원인 진단")
    print("="*80)
    
    model.objective = 'Growth'
    atpm_rxn = model.reactions.get_by_id('ATPM')
    biomass_rxn = model.reactions.get_by_id('Growth')
    
    # ATPM=0으로 설정
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 0
    
    model = setup_acetate_medium(model)
    
    # FBA 수행
    print("\n[FBA 수행]")
    solution = model.optimize()
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    
    if solution.status == 'optimal' and solution.objective_value < 1e-6:
        print("\n[문제] ATPM=0일 때 성장률이 0입니다.")
        
        # 1. ATPM 제약을 제거하고 테스트
        print("\n[테스트 1] ATPM 제약 없이 (하한 0, 상한 1000)")
        atpm_rxn.lower_bound = 0
        atpm_rxn.upper_bound = 1000
        solution2 = model.optimize()
        print(f"  성장률: {solution2.objective_value:.6f}")
        print(f"  ATPM 플럭스: {solution2.fluxes.get('ATPM', 0):.6f}")
        
        if solution2.objective_value > 1e-6:
            print("  -> ATPM 제약을 완화하면 성장 가능!")
        
        # 2. Biomass 반응 분석
        print("\n[테스트 2] Biomass 반응 분석")
        print(f"  Biomass 반응 ID: {biomass_rxn.id}")
        print(f"  Biomass 반응식: {biomass_rxn.reaction[:200]}...")
        print(f"  Biomass 하한: {biomass_rxn.lower_bound}")
        print(f"  Biomass 상한: {biomass_rxn.upper_bound}")
        
        # 3. Exchange 반응 확인
        print("\n[테스트 3] 주요 Exchange 반응 플럭스 (ATPM=0)")
        key_exchanges = ['EX_ac_e', 'EX_o2_e', 'EX_co2_e', 'EX_nh4_e', 'EX_h2o_e']
        for ex_id in key_exchanges:
            try:
                ex_rxn = model.reactions.get_by_id(ex_id)
                flux = solution.fluxes.get(ex_id, 0.0)
                print(f"  {ex_id}: {flux:.6f} (하한: {ex_rxn.lower_bound}, 상한: {ex_rxn.upper_bound})")
            except KeyError:
                print(f"  {ex_id}: 반응 없음")
        
        # 4. 주요 대사 경로 플럭스 확인
        print("\n[테스트 4] 주요 대사 경로 플럭스 (ATPM=0)")
        key_reactions = ['CS', 'ICL', 'MALS', 'ICDHx', 'ACS', 'SUCOAACTr']
        for rxn_id in key_reactions:
            try:
                rxn = model.reactions.get_by_id(rxn_id)
                flux = solution.fluxes.get(rxn_id, 0.0)
                print(f"  {rxn_id}: {flux:.6f}")
            except KeyError:
                print(f"  {rxn_id}: 반응 없음")
        
        # 5. ATP 생성/소모 확인
        print("\n[테스트 5] ATP 관련 반응 확인")
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
                    if coeff > 0:  # 생성
                        atp_producing.append((rxn.id, flux * coeff))
                    elif coeff < 0:  # 소모
                        atp_consuming.append((rxn.id, abs(flux * coeff)))
            
            print(f"  ATP 생성 반응: {len(atp_producing)}개")
            for rxn_id, net_flux in sorted(atp_producing, key=lambda x: x[1], reverse=True)[:5]:
                print(f"    {rxn_id}: {net_flux:.6f}")
            print(f"  ATP 소모 반응: {len(atp_consuming)}개")
            for rxn_id, net_flux in sorted(atp_consuming, key=lambda x: x[1], reverse=True)[:5]:
                print(f"    {rxn_id}: {net_flux:.6f}")
        except KeyError:
            print("  atp_c 메타볼라이트 없음")
        
        # 6. 제약 조건 확인
        print("\n[테스트 6] 제약 조건 확인")
        print(f"  ATPM 하한: {atpm_rxn.lower_bound}")
        print(f"  ATPM 상한: {atpm_rxn.upper_bound}")
        
        # 7. 레퍼런스 모델과 비교 제안
        print("\n[해결 방안]")
        print("  1. 레퍼런스 모델과 직접 비교 필요")
        print("  2. 누락된 필수 반응 확인 필요")
        print("  3. 미디어 설정 확인 필요")
        print("  4. 모델 제약 조건 확인 필요")

def compare_with_reference():
    """레퍼런스 모델과 비교"""
    print("\n" + "="*80)
    print("레퍼런스 모델 테스트")
    print("="*80)
    
    base_path = Path(__file__).parent.parent
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    
    if ref_model_path.exists():
        try:
            ref_model = load_model(str(ref_model_path))
            ref_model.objective = 'Growth'
            
            # ATPM 찾기
            atpm_rxn = None
            for rxn in ref_model.reactions:
                if 'ATPM' in rxn.id.upper():
                    atpm_rxn = rxn
                    break
            
            if atpm_rxn:
                print(f"\n[레퍼런스 모델] ATPM 반응: {atpm_rxn.id}")
                atpm_rxn.lower_bound = 0
                atpm_rxn.upper_bound = 0
                
                # 미디어는 레퍼런스 모델의 원래 설정 사용
                # (이미 설정되어 있을 수 있음)
                
                solution = ref_model.optimize()
                print(f"  상태: {solution.status}")
                print(f"  성장률: {solution.objective_value:.6f}")
                print(f"  ICL: {solution.fluxes.get('ICL', 0):.6f}")
                print(f"  MALS: {solution.fluxes.get('MALS', 0):.6f}")
                
                if solution.objective_value > 1e-6:
                    print("  -> 레퍼런스 모델은 ATPM=0에서도 성장 가능!")
        except Exception as e:
            print(f"  레퍼런스 모델 테스트 실패: {e}")
    else:
        print(f"\n[레퍼런스 모델 파일 없음] {ref_model_path}")

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    
    model = load_model(str(model_path))
    
    # 진단
    diagnose_atpm0_issue(model)
    
    # 레퍼런스 모델과 비교
    compare_with_reference()

if __name__ == "__main__":
    main()
