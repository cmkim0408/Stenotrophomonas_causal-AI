#!/usr/bin/env python
"""
레퍼런스 모델에서 메타볼라이트와 Transport 반응 추가
"""

import cobra
from pathlib import Path

def load_model(model_path):
    model = cobra.io.read_sbml_model(model_path)
    return model

def add_missing_metabolites_from_reference(model, ref_model):
    """레퍼런스 모델에서 누락된 메타볼라이트 추가"""
    print("="*80)
    print("레퍼런스 모델에서 누락된 메타볼라이트 추가")
    print("="*80)
    
    added = []
    
    # 필요한 메타볼라이트 목록
    required_mets = {
        'coa_e': ('CoA', 'e'),
        'pnto__R_e': ('Pantothenate', 'e'),
        'pnto__Rc': ('Pantothenate', 'c'),  # 레퍼런스 모델의 실제 ID 확인 필요
    }
    
    for met_id, (name, comp) in required_mets.items():
        if met_id not in [m.id for m in model.metabolites]:
            try:
                # 레퍼런스 모델에서 메타볼라이트 정보 가져오기
                ref_met = ref_model.metabolites.get_by_id(met_id)
                
                new_met = cobra.Metabolite(
                    met_id,
                    name=ref_met.name if ref_met.name else name,
                    compartment=ref_met.compartment if hasattr(ref_met, 'compartment') else comp,
                    formula=ref_met.formula if hasattr(ref_met, 'formula') else None,
                    charge=ref_met.charge if hasattr(ref_met, 'charge') else None
                )
                model.add_metabolites([new_met])
                added.append(met_id)
                print(f"  추가: {met_id} ({name}, compartment={new_met.compartment})")
            except KeyError:
                # 레퍼런스 모델에도 없으면 새로 생성
                new_met = cobra.Metabolite(met_id, name=name, compartment=comp)
                model.add_metabolites([new_met])
                added.append(met_id)
                print(f"  생성: {met_id} ({name}, compartment={comp})")
    
    # pnto__Rc 대신 pnto__R_c 확인
    if 'pnto__Rc' not in [m.id for m in model.metabolites]:
        # 레퍼런스 모델에서 실제 ID 확인
        for met in ref_model.metabolites:
            if 'pnto' in met.id.lower() and met.compartment == 'c':
                if met.id not in [m.id for m in model.metabolites]:
                    try:
                        new_met = cobra.Metabolite(
                            met.id,
                            name=met.name if met.name else 'Pantothenate',
                            compartment=met.compartment,
                            formula=met.formula if hasattr(met, 'formula') else None,
                            charge=met.charge if hasattr(met, 'charge') else None
                        )
                        model.add_metabolites([new_met])
                        added.append(met.id)
                        print(f"  추가: {met.id} (레퍼런스에서 발견)")
                        break
                    except:
                        pass
    
    return added

def add_transports_with_metabolites(model, ref_model):
    """메타볼라이트 추가 후 Transport 반응 추가"""
    print("\n" + "="*80)
    print("Transport 반응 추가")
    print("="*80)
    
    added = []
    
    # 1. T_coa_e_to_coa_c
    if 'T_coa_e_to_coa_c' not in [r.id for r in model.reactions]:
        try:
            coa_e = model.metabolites.get_by_id('coa_e')
            coa_c = model.metabolites.get_by_id('coa_c')
            
            trans = cobra.Reaction('T_coa_e_to_coa_c')
            trans.name = 'CoA transport (e to c)'
            trans.lower_bound = -1000
            trans.upper_bound = 1000
            trans.add_metabolites({coa_e: -1.0, coa_c: 1.0})
            model.add_reactions([trans])
            added.append('T_coa_e_to_coa_c')
            print(f"  추가: T_coa_e_to_coa_c: {trans.reaction}")
        except KeyError as e:
            print(f"  건너뜀: T_coa_e_to_coa_c - {e}")
    
    # 2. T_pnto__R_e_to_c (레퍼런스 모델에서 실제 반응식 확인)
    if 'T_pnto__R_e_to_c' not in [r.id for r in model.reactions]:
        try:
            ref_rxn = ref_model.reactions.get_by_id('T_pnto__R_e_to_c')
            
            # 메타볼라이트 확인
            pnto_e_id = None
            pnto_c_id = None
            
            for met, coeff in ref_rxn.metabolites.items():
                if coeff < 0:
                    pnto_e_id = met.id
                elif coeff > 0:
                    pnto_c_id = met.id
            
            if pnto_e_id and pnto_c_id:
                pnto_e = model.metabolites.get_by_id(pnto_e_id)
                pnto_c = model.metabolites.get_by_id(pnto_c_id)
                
                trans = cobra.Reaction('T_pnto__R_e_to_c')
                trans.name = 'Pantothenate transport'
                trans.lower_bound = -1000
                trans.upper_bound = 1000
                trans.add_metabolites({pnto_e: -1.0, pnto_c: 1.0})
                model.add_reactions([trans])
                added.append('T_pnto__R_e_to_c')
                print(f"  추가: T_pnto__R_e_to_c: {trans.reaction}")
        except (KeyError, AttributeError) as e:
            print(f"  건너뜀: T_pnto__R_e_to_c - {e}")
    
    # 3. EX_pnto__R_e
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
            print(f"  추가: EX_pnto__R_e: {ex_pnto.reaction}")
        except KeyError as e:
            print(f"  건너뜀: EX_pnto__R_e - {e}")
    
    return added

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

def test_final(model):
    """최종 테스트"""
    print("\n" + "="*80)
    print("최종 테스트")
    print("="*80)
    
    model = setup_acetate_medium(model)
    
    # Pantothenate 부트스트랩
    try:
        ex_pnto = model.reactions.get_by_id('EX_pnto__R_e')
        ex_pnto.lower_bound = -0.001
        ex_pnto.upper_bound = 1000
        print(f"  Pantothenate 부트스트랩: EX_pnto__R_e 하한=-0.001")
    except KeyError:
        print(f"  경고: EX_pnto__R_e 없음")
    
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
    
    print(f"\n[주요 반응 플럭스]")
    print(f"  ACS: {acs_flux:.6f}")
    print(f"  CS: {cs_flux:.6f}")
    print(f"  ADK1: {adk1_flux:.6f}")
    
    if abs(acs_flux) > 1e-6:
        print(f"\n[성공] ACS 작동!")
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
    print("메타볼라이트 및 Transport 반응 추가")
    print("="*80)
    
    # 메타볼라이트 추가
    added_mets = add_missing_metabolites_from_reference(model, ref_model)
    
    # Transport 추가
    added_trans = add_transports_with_metabolites(model, ref_model)
    
    print(f"\n[추가 요약]")
    print(f"  메타볼라이트: {len(added_mets)}개")
    print(f"  Transport 반응: {len(added_trans)}개")
    
    # 최종 테스트
    success = test_final(model)
    
    print("\n" + "="*80)
    print("최종 결과")
    print("="*80)
    
    if success:
        print(f"\n[성공] 문제 해결!")
    else:
        print(f"\n[실패] 추가 조사 필요")

if __name__ == "__main__":
    main()
