#!/usr/bin/env python
"""
레퍼런스 모델에만 있는 Transport 반응 추가 및 테스트
"""

import cobra
from pathlib import Path

def load_model(model_path):
    model = cobra.io.read_sbml_model(model_path)
    return model

def setup_acetate_medium(model):
    """Acetate 미디어 설정"""
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

def add_missing_transports_from_reference(model, ref_model):
    """레퍼런스 모델에만 있는 Transport 반응 추가"""
    print("="*80)
    print("레퍼런스 모델에만 있는 Transport 반응 추가")
    print("="*80)
    
    added = []
    
    # 1. T_coa_e_to_coa_c: coa_e --> coa_c
    if 'T_coa_e_to_coa_c' not in [r.id for r in model.reactions]:
        try:
            coa_e = model.metabolites.get_by_id('coa_e')
            coa_c = model.metabolites.get_by_id('coa_c')
            
            # 레퍼런스 모델에서 반응식 확인
            ref_rxn = ref_model.reactions.get_by_id('T_coa_e_to_coa_c')
            
            trans = cobra.Reaction('T_coa_e_to_coa_c')
            trans.name = 'CoA transport (e to c)'
            trans.lower_bound = -1000
            trans.upper_bound = 1000
            trans.add_metabolites({coa_e: -1.0, coa_c: 1.0})
            model.add_reactions([trans])
            added.append('T_coa_e_to_coa_c')
            print(f"\n[추가] T_coa_e_to_coa_c: {trans.reaction}")
        except (KeyError, AttributeError) as e:
            print(f"\n[건너뜀] T_coa_e_to_coa_c: {e}")
    
    # 2. T_pnto__R_e_to_c: pnto__R_e <=> pnto__Rc
    if 'T_pnto__R_e_to_c' not in [r.id for r in model.reactions]:
        try:
            # 레퍼런스 모델에서 반응식 확인
            ref_rxn = ref_model.reactions.get_by_id('T_pnto__R_e_to_c')
            
            # 메타볼라이트 확인
            pnto_e_id = None
            pnto_c_id = None
            
            for met, coeff in ref_rxn.metabolites.items():
                if coeff < 0:  # 소모
                    pnto_e_id = met.id
                elif coeff > 0:  # 생성
                    pnto_c_id = met.id
            
            if pnto_e_id and pnto_c_id:
                try:
                    pnto_e = model.metabolites.get_by_id(pnto_e_id)
                    pnto_c = model.metabolites.get_by_id(pnto_c_id)
                    
                    trans = cobra.Reaction('T_pnto__R_e_to_c')
                    trans.name = 'Pantothenate transport'
                    trans.lower_bound = -1000
                    trans.upper_bound = 1000
                    trans.add_metabolites({pnto_e: -1.0, pnto_c: 1.0})
                    model.add_reactions([trans])
                    added.append('T_pnto__R_e_to_c')
                    print(f"\n[추가] T_pnto__R_e_to_c: {trans.reaction}")
                except KeyError as e:
                    print(f"\n[건너뜀] T_pnto__R_e_to_c: 메타볼라이트 없음 - {e}")
        except (KeyError, AttributeError) as e:
            print(f"\n[건너뜀] T_pnto__R_e_to_c: {e}")
    
    # 3. EX_pnto__R_e 추가 (Exchange)
    if 'EX_pnto__R_e' not in [r.id for r in model.exchanges]:
        try:
            pnto_e = model.metabolites.get_by_id('pnto__R_e')
            
            ex_pnto = cobra.Reaction('EX_pnto__R_e')
            ex_pnto.name = 'Pantothenate exchange'
            ex_pnto.lower_bound = -1000
            ex_pnto.upper_bound = 1000
            ex_pnto.add_metabolites({pnto_e: -1.0})
            model.add_reactions([ex_pnto])
            added.append('EX_pnto__R_e')
            print(f"\n[추가] EX_pnto__R_e: {ex_pnto.reaction}")
        except KeyError as e:
            print(f"\n[건너뜀] EX_pnto__R_e: 메타볼라이트 없음 - {e}")
    
    return added

def test_with_all_fixes(model):
    """모든 수정사항 적용 후 테스트"""
    print("\n" + "="*80)
    print("모든 수정사항 적용 후 테스트")
    print("="*80)
    
    model = setup_acetate_medium(model)
    
    # Pantothenate 부트스트랩
    try:
        ex_pnto = model.reactions.get_by_id('EX_pnto__R_e')
        ex_pnto.lower_bound = -0.001
        ex_pnto.upper_bound = 1000
        print(f"\n[Pantothenate 부트스트랩] EX_pnto__R_e: 하한=-0.001")
    except KeyError:
        print(f"\n[경고] EX_pnto__R_e 없음")
    
    # ATPM=0 설정
    atpm_rxn = model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    model.objective = 'Growth'
    solution = model.optimize()
    
    print(f"\n[FBA 결과]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    
    acs_flux = solution.fluxes.get('ACS', 0.0)
    cs_flux = solution.fluxes.get('CS', 0.0)
    adk1_flux = solution.fluxes.get('ADK1', 0.0)
    icl_flux = solution.fluxes.get('ICL', 0.0)
    mals_flux = solution.fluxes.get('MALS', 0.0)
    
    print(f"\n[주요 반응 플럭스]")
    print(f"  ACS: {acs_flux:.6f}")
    print(f"  CS: {cs_flux:.6f}")
    print(f"  ADK1: {adk1_flux:.6f}")
    print(f"  ICL: {icl_flux:.6f}")
    print(f"  MALS: {mals_flux:.6f}")
    
    # Transport 플럭스 확인
    print(f"\n[Transport 플럭스]")
    transports = ['T_coa_e_to_coa_c', 'T_pnto__R_e_to_c', 'T_pnto__R_e_to_pnto__R_c']
    for trans_id in transports:
        try:
            flux = solution.fluxes.get(trans_id, 0.0)
            if abs(flux) > 1e-6:
                print(f"  {trans_id}: {flux:.6f}")
        except:
            pass
    
    # Exchange 플럭스 확인
    print(f"\n[Exchange 플럭스]")
    exchanges = ['EX_pnto__R_e', 'EX_coa_c']
    for ex_id in exchanges:
        try:
            flux = solution.fluxes.get(ex_id, 0.0)
            if abs(flux) > 1e-6:
                print(f"  {ex_id}: {flux:.6f}")
        except:
            pass
    
    if abs(acs_flux) > 1e-6:
        print(f"\n[성공] ACS 작동!")
        if solution.objective_value > 1e-6:
            print(f"  성장도 가능! (성장률: {solution.objective_value:.6f})")
        return True
    else:
        print(f"\n[실패] ACS 여전히 작동 안 함")
        return False

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    
    model = load_model(str(model_path))
    
    if not ref_model_path.exists():
        print(f"[오류] 레퍼런스 모델 파일 없음: {ref_model_path}")
        return
    
    ref_model = cobra.io.read_sbml_model(str(ref_model_path))
    
    print("="*80)
    print("레퍼런스 모델 Transport 반응 추가 및 테스트")
    print("="*80)
    
    # 레퍼런스 모델에만 있는 Transport 추가
    added = add_missing_transports_from_reference(model, ref_model)
    
    print(f"\n[추가된 반응] {len(added)}개: {', '.join(added)}")
    
    # 모든 수정사항 적용 후 테스트
    success = test_with_all_fixes(model)
    
    print("\n" + "="*80)
    print("최종 결과")
    print("="*80)
    
    if success:
        print(f"\n[성공] Transport 반응 추가로 문제 해결!")
        print(f"  -> ACS가 작동함")
        print(f"  -> 경로가 활성화됨")
    else:
        print(f"\n[실패] Transport 반응 추가로도 해결되지 않음")
        print(f"  -> 추가 조사 필요")

if __name__ == "__main__":
    main()
