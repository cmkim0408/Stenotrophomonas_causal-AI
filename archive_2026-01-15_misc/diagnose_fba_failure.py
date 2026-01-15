#!/usr/bin/env python
"""
FBA 실패 원인 진단
Acetate에서 성장하지 못하는 이유 분석
"""

import cobra
from cobra.flux_analysis import find_blocked_reactions

def load_model(model_path="BaseModel.xml"):
    """모델 로드"""
    print(f"모델 로드 중: {model_path}")
    model = cobra.io.read_sbml_model(model_path)
    print(f"[OK] 모델 로드 완료: {model.id}\n")
    return model

def find_biomass_reaction(model):
    """Biomass 반응 찾기"""
    for name in ['Growth', 'BIOMASS', 'BIOMASS_Ecoli_core']:
        try:
            return model.reactions.get_by_id(name)
        except KeyError:
            continue
    return None

def setup_acetate_medium(model):
    """Acetate medium 설정"""
    for rxn in model.exchanges:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    
    try:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.lower_bound = -1000
        ex_ac.upper_bound = 1000
    except KeyError:
        pass
    
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
    
    return model

def check_pathway_connectivity(model):
    """경로 연결성 확인"""
    print("="*70)
    print("경로 연결성 확인")
    print("="*70)
    
    # Acetate → Acetyl-CoA
    try:
        ac_e = model.metabolites.get_by_id('ac_e')
        ac_c = model.metabolites.get_by_id('ac_c')
        
        # Exchange 확인
        ac_ex_rxns = [r for r in ac_e.reactions if 'EX' in r.id]
        print(f"\n[Acetate Exchange]")
        if ac_ex_rxns:
            for rxn in ac_ex_rxns:
                print(f"  {rxn.id}: LB={rxn.lower_bound}, UB={rxn.upper_bound}")
        else:
            print(f"  [ERROR] Acetate exchange 반응 없음")
        
        # Transport 확인
        transport_rxns = set(ac_e.reactions) & set(ac_c.reactions)
        print(f"\n[Acetate Transport]")
        if transport_rxns:
            for rxn in transport_rxns:
                print(f"  {rxn.id}: {rxn.reaction}")
                print(f"    LB={rxn.lower_bound}, UB={rxn.upper_bound}")
        else:
            print(f"  [ERROR] Acetate transport 반응 없음")
        
    except KeyError as e:
        print(f"[ERROR] {e} metabolite 없음")
    
    # Acetate → Acetyl-CoA 반응 확인
    try:
        accoa_c = model.metabolites.get_by_id('accoa_c')
        acs_rxns = [r for r in ac_c.reactions if accoa_c in r.metabolites]
        print(f"\n[Acetate → Acetyl-CoA]")
        if acs_rxns:
            for rxn in acs_rxns:
                print(f"  {rxn.id}: {rxn.reaction}")
                print(f"    LB={rxn.lower_bound}, UB={rxn.upper_bound}")
                print(f"    유전자: {[g.id for g in rxn.genes]}")
        else:
            print(f"  [ERROR] ACS 반응 없음")
    except KeyError:
        pass
    
    # TCA/Glyoxylate 경로 확인
    tca_key_rxns = ['CS', 'ICL', 'MALS', 'MDH']
    print(f"\n[TCA/Glyoxylate 핵심 반응]")
    for rxn_id in tca_key_rxns:
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            print(f"  {rxn_id}: LB={rxn.lower_bound}, UB={rxn.upper_bound}")
        except KeyError:
            print(f"  {rxn_id}: [MISSING]")

def check_blocked_reactions(model):
    """블록된 반응 확인"""
    print("\n" + "="*70)
    print("블록된 반응 확인 (핵심 경로)")
    print("="*70)
    
    key_reactions = ['EX_ac_e', 'ACt', 'R_ACS', 'ACS', 'CS', 'ICL', 'MALS', 
                     'ICDHx', 'AKGDH', 'SUCD', 'FUM', 'MDH', 'Growth']
    
    blocked = find_blocked_reactions(model)
    blocked_ids = [rxn.id if hasattr(rxn, 'id') else rxn for rxn in blocked]
    
    print(f"\n핵심 반응 블록 상태:")
    for rxn_id in key_reactions:
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            if rxn_id in blocked_ids:
                print(f"  {rxn_id}: [BLOCKED] X")
            else:
                print(f"  {rxn_id}: [OK] O")
        except KeyError:
            print(f"  {rxn_id}: [MISSING] X")

def test_metabolite_production(model, metabolite_id):
    """특정 대사물질 생산 가능 여부 확인"""
    try:
        met = model.metabolites.get_by_id(metabolite_id)
        
        # 임시 exchange 반응 생성
        test_rxn = cobra.Reaction(f'TEST_{metabolite_id}')
        test_rxn.add_metabolites({met: 1})
        test_rxn.lower_bound = 0
        test_rxn.upper_bound = 1000
        
        model.add_reactions([test_rxn])
        model.objective = test_rxn.id
        
        solution = model.optimize()
        
        # 테스트 반응 제거
        model.remove_reactions([test_rxn])
        
        return solution.status == 'optimal' and solution.objective_value > 1e-6
    except:
        return False

def diagnose_metabolite_production(model):
    """주요 대사물질 생산 가능 여부 확인"""
    print("\n" + "="*70)
    print("주요 대사물질 생산 가능 여부 확인")
    print("="*70)
    
    metabolites = {
        'ac_c': 'Acetate (세포질)',
        'accoa_c': 'Acetyl-CoA',
        'cit_c': 'Citrate',
        'icit_c': 'Isocitrate',
        'glx_c': 'Glyoxylate',
        'mal__L_c': 'Malate',
        'oaa_c': 'Oxaloacetate'
    }
    
    for met_id, met_name in metabolites.items():
        can_produce = test_metabolite_production(model, met_id)
        status = "[OK] 생산 가능" if can_produce else "[FAIL] 생산 불가"
        print(f"  {met_name} ({met_id}): {status}")

def check_biomass_components(model, biomass_rxn):
    """Biomass 구성 요소 확인"""
    print("\n" + "="*70)
    print("Biomass 구성 요소 확인")
    print("="*70)
    
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"\nBiomass 반응: {biomass_rxn.id}")
    print(f"반응식: {biomass_rxn.reaction}")
    
    # 주요 대사물질 생산 가능 여부 확인
    metabolites = biomass_rxn.metabolites
    print(f"\n주요 구성 요소 (절대값 상위 10개):")
    
    sorted_mets = sorted(metabolites.items(), key=lambda x: abs(x[1]), reverse=True)
    for met, coeff in sorted_mets[:10]:
        can_produce = test_metabolite_production(model, met.id)
        status = "[OK]" if can_produce else "[FAIL]"
        print(f"  {status} {met.id}: {coeff:.4f}")

def main():
    """메인 함수"""
    print("="*70)
    print("FBA 실패 원인 진단")
    print("="*70)
    
    # 1. 모델 로드
    model = load_model("BaseModel.xml")
    
    # 2. Biomass 반응 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응을 찾을 수 없습니다!")
        return
    
    print(f"[OK] Biomass 반응: {biomass_rxn.id}\n")
    
    # 3. Acetate medium 설정
    model = setup_acetate_medium(model)
    
    # 4. 경로 연결성 확인
    check_pathway_connectivity(model)
    
    # 5. 블록된 반응 확인
    check_blocked_reactions(model)
    
    # 6. 대사물질 생산 가능 여부 확인
    diagnose_metabolite_production(model)
    
    # 7. Biomass 구성 요소 확인
    check_biomass_components(model, biomass_rxn)
    
    print("\n" + "="*70)
    print("진단 완료")
    print("="*70)

if __name__ == "__main__":
    main()
