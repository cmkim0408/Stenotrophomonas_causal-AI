#!/usr/bin/env python
"""
Acetate exchange 추가 후 성장 테스트
"""

import cobra
from cobra import Reaction

def load_model(model_path="BaseModel.xml"):
    """모델 로드"""
    print(f"모델 로드 중: {model_path}")
    model = cobra.io.read_sbml_model(model_path)
    print(f"[OK] 모델 로드 완료: {model.id}\n")
    return model

def ensure_acetate_exchange(model):
    """Acetate exchange 반응 확인 및 추가"""
    try:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        print(f"[OK] EX_ac_e 이미 존재: {ex_ac.reaction}")
    except KeyError:
        print("[INFO] EX_ac_e 없음, 추가 중...")
        try:
            ac_e = model.metabolites.get_by_id('ac_e')
        except KeyError:
            from cobra import Metabolite
            ac_c = model.metabolites.get_by_id('ac_c')
            ac_e = Metabolite(
                id='ac_e',
                formula=ac_c.formula,
                name='Acetate',
                compartment='C_e',
                charge=ac_c.charge
            )
            model.add_metabolites([ac_e])
            print(f"[OK] ac_e 생성")
        
        ex_ac = Reaction('EX_ac_e')
        ex_ac.name = 'Acetate exchange'
        ex_ac.lower_bound = -1000
        ex_ac.upper_bound = 1000
        ex_ac.add_metabolites({ac_e: -1})
        model.add_reactions([ex_ac])
        print(f"[OK] EX_ac_e 추가: {ex_ac.reaction}")
    
    return model

def setup_acetate_medium(model):
    """Acetate minimal medium 설정"""
    # 모든 exchange 차단
    for rxn in model.exchanges:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    
    # Acetate 허용
    try:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.lower_bound = -1000
        print(f"[OK] Acetate 허용: {ex_ac.id}")
    except KeyError:
        print("[ERROR] EX_ac_e를 찾을 수 없습니다!")
        return None
    
    # 필수 무기염/영양소
    essentials = {
        'EX_nh4_e': 'Ammonium',
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
    
    print("\n필수 영양소 설정:")
    for ex_id, name in essentials.items():
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            if ex_id == 'EX_co2_e':
                ex_rxn.lower_bound = -1000
                ex_rxn.upper_bound = 1000
            elif ex_id == 'EX_o2_e':
                ex_rxn.lower_bound = -1000
                ex_rxn.upper_bound = 1000
            else:
                ex_rxn.lower_bound = -1000
            print(f"  [OK] {name}: {ex_id}")
        except KeyError:
            pass
    
    return model

def test_growth(model):
    """성장 테스트"""
    print("\n" + "="*70)
    print("FBA 최적화 수행")
    print("="*70)
    
    # Biomass 반응 찾기
    biomass_rxn = None
    for name in ['Growth', 'BIOMASS', 'BIOMASS_Ecoli_core']:
        try:
            biomass_rxn = model.reactions.get_by_id(name)
            break
        except KeyError:
            continue
    
    if not biomass_rxn:
        print("[ERROR] Biomass 반응을 찾을 수 없습니다!")
        return False, 0
    
    model.objective = biomass_rxn.id
    print(f"[OK] Objective: {biomass_rxn.id}")
    
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
                return False, 0
        else:
            print(f"\n[FAIL] 최적화 실패: {solution.status}")
            return False, 0
            
    except Exception as e:
        print(f"\n[ERROR] FBA 최적화 중 오류: {e}")
        import traceback
        traceback.print_exc()
        return False, 0

def main():
    """메인 함수"""
    # 모델 로드
    model = load_model("BaseModel.xml")
    
    # Acetate exchange 확인 및 추가
    model = ensure_acetate_exchange(model)
    
    # 모델 저장
    try:
        cobra.io.write_sbml_model(model, "BaseModel.xml")
        print("\n[OK] BaseModel.xml 업데이트 완료")
    except Exception as e:
        print(f"\n[WARNING] 모델 저장 실패: {e}")
    
    # Acetate medium 설정
    model = setup_acetate_medium(model)
    
    # 성장 테스트
    growth_ok, biomass_flux = test_growth(model)
    
    print("\n" + "="*70)
    if growth_ok:
        print("[SUCCESS] Acetate minimal medium에서 성장 가능!")
    else:
        print("[FAIL] Acetate minimal medium에서 성장 불가")
    print("="*70)

if __name__ == "__main__":
    main()


