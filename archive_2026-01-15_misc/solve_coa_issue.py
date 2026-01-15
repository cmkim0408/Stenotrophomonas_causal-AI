#!/usr/bin/env python
"""
CoA 문제 해결: 2번 → 1번 → 3번 순서로 수행
1. CoA 부트스트랩 제공 (2번)
2. CoA 합성 경로 확인 및 수정 (1번)
3. 레퍼런스 모델과 비교 (3번)
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

def solution_2_coa_bootstrap(model):
    """해결 방안 2번: CoA 부트스트랩 제공"""
    print("="*80)
    print("해결 방안 2번: CoA 부트스트랩 제공 (임시)")
    print("="*80)
    
    model = setup_acetate_medium(model)
    
    # CoA 부트스트랩 추가
    try:
        coa_c = model.metabolites.get_by_id('coa_c')
        ex_coa_id = 'EX_coa_c'
        
        if ex_coa_id not in [r.id for r in model.exchanges]:
            ex_coa = cobra.Reaction(ex_coa_id)
            ex_coa.name = 'CoA exchange'
            ex_coa.lower_bound = -0.001  # 아주 약간의 부트스트랩
            ex_coa.upper_bound = 1000
            
            coa_e_id = 'coa_e'
            try:
                coa_e = model.metabolites.get_by_id(coa_e_id)
            except KeyError:
                coa_e = cobra.Metabolite(coa_e_id, name='CoA', compartment='e')
                model.add_metabolites([coa_e])
            
            ex_coa.add_metabolites({coa_e: -1.0})
            model.add_reactions([ex_coa])
            print(f"\n[CoA 부트스트랩 추가]")
            print(f"  EX_coa_c: 하한={ex_coa.lower_bound}")
        else:
            ex_coa = model.reactions.get_by_id(ex_coa_id)
            ex_coa.lower_bound = -0.001
            ex_coa.upper_bound = 1000
            print(f"\n[CoA 부트스트랩 설정]")
            print(f"  EX_coa_c: 하한={ex_coa.lower_bound}")
    except KeyError:
        print(f"\n[오류] coa_c 메타볼라이트 없음")
        return None
    
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
    
    if abs(acs_flux) > 1e-6:
        print(f"\n[성공] CoA 부트스트랩으로 ACS 작동!")
        return True, solution
    else:
        print(f"\n[실패] CoA 부트스트랩으로도 ACS 작동 안 함")
        return False, solution

def solution_1_coa_synthesis_pathway(model):
    """해결 방안 1번: CoA 합성 경로 확인 및 수정"""
    print("\n" + "="*80)
    print("해결 방안 1번: CoA 합성 경로 확인 및 수정")
    print("="*80)
    
    model = setup_acetate_medium(model)
    
    # CoA 합성 경로 찾기
    try:
        coa_c = model.metabolites.get_by_id('coa_c')
        
        print(f"\n[CoA 합성 관련 반응 찾기]")
        coa_synthesis_keywords = ['COA', 'PNTO', 'panto', 'CoA.*synth', 'PAN', 'PPAT', 'DPCK']
        
        coa_reactions = []
        for rxn in model.reactions:
            rxn_id_upper = rxn.id.upper()
            rxn_name_upper = (rxn.name or "").upper()
            
            for keyword in coa_synthesis_keywords:
                if keyword in rxn_id_upper or keyword in rxn_name_upper:
                    try:
                        if 'coa_c' in [met.id for met in rxn.metabolites]:
                            coa_reactions.append(rxn)
                            break
                    except:
                        pass
        
        print(f"  총 {len(coa_reactions)}개 반응 발견")
        
        # CoA를 생성하는 반응 찾기
        coa_producing = []
        for rxn in coa_c.reactions:
            coeff = rxn.metabolites.get(coa_c, 0)
            if coeff > 0:  # 생성
                coa_producing.append((rxn.id, coeff, rxn.reaction))
        
        print(f"\n[CoA 생성 반응] {len(coa_producing)}개")
        for rxn_id, coeff, reaction in coa_producing[:10]:
            print(f"  {rxn_id}: 계수={coeff}")
            print(f"    반응식: {reaction[:100]}")
        
        # CoA 전구체 확인 (pantothenate)
        print(f"\n[CoA 전구체 확인]")
        panto_keywords = ['pnto', 'panto', 'PAN']
        panto_reactions = []
        for rxn in model.reactions:
            rxn_id_upper = rxn.id.upper()
            for keyword in panto_keywords:
                if keyword in rxn_id_upper:
                    panto_reactions.append(rxn)
                    break
        
        print(f"  Pantothenate 관련 반응: {len(panto_reactions)}개")
        for rxn in panto_reactions[:5]:
            print(f"    {rxn.id}: {rxn.reaction[:80]}")
        
        # Pantothenate exchange 확인
        panto_exchanges = ['EX_pnto__R_e', 'EX_pantothenate_e', 'EX_panto__R_e']
        print(f"\n[Pantothenate Exchange 확인]")
        for ex_id in panto_exchanges:
            if ex_id in [r.id for r in model.exchanges]:
                ex_rxn = model.reactions.get_by_id(ex_id)
                print(f"  {ex_id}: 존재 (bounds: [{ex_rxn.lower_bound}, {ex_rxn.upper_bound}])")
            else:
                print(f"  {ex_id}: 없음")
        
        # Pantothenate exchange 추가 (부트스트랩)
        print(f"\n[Pantothenate 부트스트랩 추가]")
        try:
            pnto_c = model.metabolites.get_by_id('pnto__R_c')
            ex_pnto_id = 'EX_pnto__R_e'
            
            if ex_pnto_id not in [r.id for r in model.exchanges]:
                ex_pnto = cobra.Reaction(ex_pnto_id)
                ex_pnto.name = 'Pantothenate exchange'
                ex_pnto.lower_bound = -0.001
                ex_pnto.upper_bound = 1000
                
                pnto_e_id = 'pnto__R_e'
                try:
                    pnto_e = model.metabolites.get_by_id(pnto_e_id)
                except KeyError:
                    pnto_e = cobra.Metabolite(pnto_e_id, name='Pantothenate', compartment='e')
                    model.add_metabolites([pnto_e])
                
                ex_pnto.add_metabolites({pnto_e: -1.0})
                model.add_reactions([ex_pnto])
                print(f"  추가: {ex_pnto_id} (하한: -0.001)")
            else:
                ex_pnto = model.reactions.get_by_id(ex_pnto_id)
                ex_pnto.lower_bound = -0.001
                ex_pnto.upper_bound = 1000
                print(f"  설정: {ex_pnto_id} (하한: -0.001)")
        except KeyError:
            print(f"  pnto__R_c 메타볼라이트 없음")
        
        # Transport 확인
        print(f"\n[Pantothenate Transport 확인]")
        pnto_transports = ['T_pnto__R_e_to_pnto__R_c', 'PANTOTe', 'PNTOt']
        for trans_id in pnto_transports:
            if trans_id in [r.id for r in model.reactions]:
                trans_rxn = model.reactions.get_by_id(trans_id)
                print(f"  {trans_id}: 존재 - {trans_rxn.reaction[:80]}")
            else:
                print(f"  {trans_id}: 없음")
        
        # Transport 추가 (필요시)
        if 'T_pnto__R_e_to_pnto__R_c' not in [r.id for r in model.reactions]:
            try:
                pnto_e = model.metabolites.get_by_id('pnto__R_e')
                pnto_c = model.metabolites.get_by_id('pnto__R_c')
                
                trans = cobra.Reaction('T_pnto__R_e_to_pnto__R_c')
                trans.name = 'Pantothenate transport'
                trans.lower_bound = -1000
                trans.upper_bound = 1000
                trans.add_metabolites({pnto_e: -1.0, pnto_c: 1.0})
                model.add_reactions([trans])
                print(f"\n  Transport 추가: T_pnto__R_e_to_pnto__R_c")
            except KeyError:
                print(f"\n  Transport 추가 실패: 메타볼라이트 없음")
        
    except KeyError:
        print(f"  coa_c 메타볼라이트 없음")
        return None
    
    # 테스트
    atpm_rxn = model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    model.objective = 'Growth'
    solution = model.optimize()
    
    print(f"\n[FBA 결과 (Pantothenate 부트스트랩)]")
    print(f"  상태: {solution.status}")
    print(f"  성장률: {solution.objective_value:.6f}")
    
    acs_flux = solution.fluxes.get('ACS', 0.0)
    print(f"  ACS: {acs_flux:.6f}")
    
    if abs(acs_flux) > 1e-6:
        print(f"\n[성공] Pantothenate 부트스트랩으로 ACS 작동!")
        return True, solution
    else:
        print(f"\n[실패] Pantothenate 부트스트랩으로도 ACS 작동 안 함")
        return False, solution

def solution_3_compare_with_reference(model):
    """해결 방안 3번: 레퍼런스 모델과 비교"""
    print("\n" + "="*80)
    print("해결 방안 3번: 레퍼런스 모델과 비교")
    print("="*80)
    
    base_path = Path(__file__).parent.parent
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    
    if not ref_model_path.exists():
        print(f"\n[레퍼런스 모델 파일 없음] {ref_model_path}")
        return None
    
    try:
        ref_model = cobra.io.read_sbml_model(str(ref_model_path))
        print(f"\n[레퍼런스 모델 로드 완료]")
        
        # CoA 합성 관련 반응 비교
        print(f"\n[CoA 합성 관련 반응 비교]")
        
        coa_keywords = ['COA', 'PNTO', 'panto', 'PAN', 'PPAT', 'DPCK']
        
        new_coa_rxns = set()
        for rxn in model.reactions:
            rxn_id_upper = rxn.id.upper()
            for keyword in coa_keywords:
                if keyword in rxn_id_upper:
                    try:
                        if 'coa_c' in [met.id for met in rxn.metabolites]:
                            new_coa_rxns.add(rxn.id)
                            break
                    except:
                        pass
        
        ref_coa_rxns = set()
        for rxn in ref_model.reactions:
            rxn_id_upper = rxn.id.upper()
            for keyword in coa_keywords:
                if keyword in rxn_id_upper:
                    try:
                        if 'coa_c' in [met.id for met in rxn.metabolites]:
                            ref_coa_rxns.add(rxn.id)
                            break
                    except:
                        pass
        
        print(f"  신규 모델 CoA 관련 반응: {len(new_coa_rxns)}개")
        print(f"  레퍼런스 모델 CoA 관련 반응: {len(ref_coa_rxns)}개")
        
        # 레퍼런스에만 있는 반응
        ref_only = ref_coa_rxns - new_coa_rxns
        if ref_only:
            print(f"\n[레퍼런스에만 있는 CoA 관련 반응] {len(ref_only)}개")
            for rxn_id in sorted(ref_only)[:10]:
                try:
                    rxn = ref_model.reactions.get_by_id(rxn_id)
                    print(f"  {rxn_id}: {rxn.reaction[:80]}")
                except:
                    pass
        
        # Pantothenate 관련 비교
        print(f"\n[Pantothenate 관련 반응 비교]")
        
        new_panto = set([r.id for r in model.reactions if 'PNTO' in r.id.upper() or 'PANTO' in r.id.upper()])
        ref_panto = set([r.id for r in ref_model.reactions if 'PNTO' in r.id.upper() or 'PANTO' in r.id.upper()])
        
        print(f"  신규 모델: {len(new_panto)}개")
        print(f"  레퍼런스 모델: {len(ref_panto)}개")
        
        ref_only_panto = ref_panto - new_panto
        if ref_only_panto:
            print(f"\n[레퍼런스에만 있는 Pantothenate 관련 반응] {len(ref_only_panto)}개")
            for rxn_id in sorted(ref_only_panto)[:10]:
                try:
                    rxn = ref_model.reactions.get_by_id(rxn_id)
                    print(f"  {rxn_id}: {rxn.reaction[:80]}")
                except:
                    pass
        
        # Exchange 비교
        print(f"\n[Exchange 비교]")
        new_ex = set([r.id for r in model.exchanges if 'PNTO' in r.id.upper() or 'PANTO' in r.id.upper()])
        ref_ex = set([r.id for r in ref_model.exchanges if 'PNTO' in r.id.upper() or 'PANTO' in r.id.upper()])
        
        print(f"  신규 모델 Exchange: {new_ex}")
        print(f"  레퍼런스 모델 Exchange: {ref_ex}")
        
        return ref_only, ref_only_panto
        
    except Exception as e:
        print(f"\n[오류] 레퍼런스 모델 비교 실패: {e}")
        return None

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    
    model = load_model(str(model_path))
    
    print("="*80)
    print("CoA 문제 해결: 2번 → 1번 → 3번 순서로 수행")
    print("="*80)
    
    # 2번: CoA 부트스트랩 제공
    success_2, solution_2 = solution_2_coa_bootstrap(model)
    
    # 1번: CoA 합성 경로 확인 및 수정
    success_1, solution_1 = solution_1_coa_synthesis_pathway(model)
    
    # 3번: 레퍼런스 모델과 비교
    ref_comparison = solution_3_compare_with_reference(model)
    
    print("\n" + "="*80)
    print("최종 결과")
    print("="*80)
    
    print(f"\n[2번: CoA 부트스트랩]")
    if success_2:
        print(f"  성공: ACS 작동!")
    else:
        print(f"  실패: ACS 작동 안 함")
    
    print(f"\n[1번: CoA 합성 경로]")
    if success_1:
        print(f"  성공: Pantothenate 부트스트랩으로 ACS 작동!")
    else:
        print(f"  실패: Pantothenate 부트스트랩으로도 ACS 작동 안 함")
    
    print(f"\n[3번: 레퍼런스 비교]")
    if ref_comparison:
        ref_only, ref_only_panto = ref_comparison
        if ref_only or ref_only_panto:
            print(f"  레퍼런스에만 있는 반응 발견!")
            print(f"  -> 추가 필요할 수 있음")
        else:
            print(f"  차이 없음")
    
    print(f"\n[권장 사항]")
    if success_1:
        print(f"  -> Pantothenate 부트스트랩 사용 권장")
    elif success_2:
        print(f"  -> CoA 부트스트랩 사용 (임시)")
    else:
        print(f"  -> 레퍼런스 모델의 누락된 반응 추가 필요")

if __name__ == "__main__":
    main()
