#!/usr/bin/env python
"""
CoA 합성 경로 완전성 확인
- 레퍼런스 모델의 CoA 합성 경로 확인
- 신규 모델과 비교
- 누락된 반응 찾기
"""

import cobra
from pathlib import Path

def load_model(model_path):
    model = cobra.io.read_sbml_model(model_path)
    return model

def find_coa_synthesis_pathway(model):
    """CoA 합성 경로 찾기"""
    try:
        coa_c = model.metabolites.get_by_id('coa_c')
        
        # CoA를 생성하는 반응 찾기
        coa_producing = []
        for rxn in coa_c.reactions:
            coeff = rxn.metabolites.get(coa_c, 0)
            if coeff > 0:  # 생성
                coa_producing.append(rxn)
        
        return coa_producing
    except KeyError:
        return []

def analyze_coa_synthesis(model, model_name):
    """CoA 합성 경로 분석"""
    print(f"\n{'='*80}")
    print(f"{model_name} 모델 - CoA 합성 경로 분석")
    print(f"{'='*80}")
    
    coa_producing = find_coa_synthesis_pathway(model)
    print(f"\n[CoA 생성 반응] {len(coa_producing)}개")
    
    # 주요 CoA 합성 반응 찾기
    key_keywords = ['PPAT', 'DPCK', 'PAN', 'panto', 'COA.*synth']
    
    key_reactions = []
    for rxn in coa_producing:
        rxn_id_upper = rxn.id.upper()
        rxn_name_upper = (rxn.name or "").upper()
        
        for keyword in key_keywords:
            if keyword in rxn_id_upper or keyword in rxn_name_upper:
                key_reactions.append(rxn)
                break
    
    print(f"\n[주요 CoA 합성 반응] {len(key_reactions)}개")
    for rxn in key_reactions[:10]:
        print(f"  {rxn.id}: {rxn.reaction[:100]}")
    
    # Pantothenate 관련 반응
    panto_reactions = []
    for rxn in model.reactions:
        rxn_id_upper = rxn.id.upper()
        if 'PNTO' in rxn_id_upper or 'PANTO' in rxn_id_upper or 'PAN' in rxn_id_upper:
            panto_reactions.append(rxn)
    
    print(f"\n[Pantothenate 관련 반응] {len(panto_reactions)}개")
    for rxn in panto_reactions[:10]:
        print(f"  {rxn.id}: {rxn.reaction[:80]}")
    
    return key_reactions, panto_reactions

def compare_coa_pathways(new_model, ref_model):
    """CoA 합성 경로 비교"""
    print(f"\n{'='*80}")
    print("CoA 합성 경로 비교")
    print(f"{'='*80}")
    
    new_key, new_panto = analyze_coa_synthesis(new_model, "신규")
    ref_key, ref_panto = analyze_coa_synthesis(ref_model, "레퍼런스")
    
    # 반응 ID 비교
    new_key_ids = set([r.id for r in new_key])
    ref_key_ids = set([r.id for r in ref_key])
    
    ref_only = ref_key_ids - new_key_ids
    
    if ref_only:
        print(f"\n[레퍼런스에만 있는 주요 CoA 합성 반응] {len(ref_only)}개")
        for rxn_id in sorted(ref_only):
            try:
                rxn = ref_model.reactions.get_by_id(rxn_id)
                print(f"  {rxn_id}: {rxn.reaction[:100]}")
            except:
                pass
    
    # Pantothenate 반응 비교
    new_panto_ids = set([r.id for r in new_panto])
    ref_panto_ids = set([r.id for r in ref_panto])
    
    ref_only_panto = ref_panto_ids - new_panto_ids
    
    if ref_only_panto:
        print(f"\n[레퍼런스에만 있는 Pantothenate 관련 반응] {len(ref_only_panto)}개")
        for rxn_id in sorted(ref_only_panto)[:10]:
            try:
                rxn = ref_model.reactions.get_by_id(rxn_id)
                print(f"  {rxn_id}: {rxn.reaction[:80]}")
            except:
                pass

def test_coa_demand_with_pantothenate(model):
    """Pantothenate 제공 시 CoA 생성 테스트"""
    print(f"\n{'='*80}")
    print("Pantothenate 제공 시 CoA 생성 테스트")
    print(f"{'='*80}")
    
    # 미디어 설정
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    # Pantothenate 제공
    try:
        ex_pnto = model.reactions.get_by_id('EX_pnto__R_e')
        ex_pnto.lower_bound = -0.001
        ex_pnto.upper_bound = 1000
        print(f"  Pantothenate 제공: EX_pnto__R_e 하한=-0.001")
    except KeyError:
        print(f"  경고: EX_pnto__R_e 없음")
        return None
    
    # CoA demand
    try:
        coa_c = model.metabolites.get_by_id('coa_c')
        
        coa_demand = cobra.Reaction('DM_coa_c')
        coa_demand.name = 'CoA demand'
        coa_demand.lower_bound = 0
        coa_demand.upper_bound = 1000
        coa_demand.add_metabolites({coa_c: -1.0})
        model.add_reactions([coa_demand])
        
        model.objective = 'DM_coa_c'
        solution = model.optimize()
        
        print(f"\n[CoA 생성 가능 여부]")
        print(f"  상태: {solution.status}")
        print(f"  CoA 최대 생산량: {solution.objective_value:.6f}")
        
        if solution.objective_value > 1e-6:
            print(f"  -> Pantothenate 제공 시 CoA 생성 가능!")
            
            # CoA 생성 반응 확인
            coa_producing = []
            for rxn in coa_c.reactions:
                if 'DM_' in rxn.id:
                    continue
                flux = solution.fluxes.get(rxn.id, 0.0)
                if abs(flux) > 1e-6:
                    coeff = rxn.metabolites.get(coa_c, 0)
                    if coeff > 0:
                        coa_producing.append((rxn.id, flux * coeff))
            
            print(f"\n[CoA 생성 반응 (플럭스 > 0)]")
            for rxn_id, net_flux in sorted(coa_producing, key=lambda x: x[1], reverse=True)[:5]:
                print(f"  {rxn_id}: {net_flux:.6f}")
        else:
            print(f"  -> Pantothenate 제공해도 CoA 생성 불가!")
            print(f"  -> CoA 합성 경로에 문제가 있음")
        
        model.remove_reactions([coa_demand])
        
        return solution.objective_value > 1e-6
        
    except KeyError:
        print(f"  coa_c 메타볼라이트 없음")
        return None

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    
    new_model = load_model(str(model_path))
    
    if not ref_model_path.exists():
        print(f"[오류] 레퍼런스 모델 파일 없음: {ref_model_path}")
        return
    
    ref_model = cobra.io.read_sbml_model(str(ref_model_path))
    
    print("="*80)
    print("CoA 합성 경로 완전성 확인")
    print("="*80)
    
    # CoA 합성 경로 비교
    compare_coa_pathways(new_model, ref_model)
    
    # Pantothenate 제공 시 CoA 생성 테스트
    can_produce_coa = test_coa_demand_with_pantothenate(new_model)
    
    print("\n" + "="*80)
    print("결론")
    print("="*80)
    
    if can_produce_coa:
        print(f"\n[OK] Pantothenate 제공 시 CoA 생성 가능!")
        print(f"  -> CoA 합성 경로는 정상")
        print(f"  -> 다른 문제가 있을 수 있음")
    else:
        print(f"\n[문제] Pantothenate 제공해도 CoA 생성 불가")
        print(f"  -> CoA 합성 경로에 문제가 있음")
        print(f"  -> 누락된 반응 추가 필요")

if __name__ == "__main__":
    main()
