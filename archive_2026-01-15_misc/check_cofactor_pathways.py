#!/usr/bin/env python
"""
보조인자 생산 경로 확인
- CoA 생산 경로 (Pantothenate -> CoA)
- NAD+/NADP+ 생산 경로
- FAD 생산 경로
- Folate 관련 경로
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

def check_coa_pathway(model):
    """CoA 생산 경로 확인"""
    print("\n" + "="*70)
    print("3단계: 보조인자 생산 경로 확인")
    print("="*70)
    print("\n[1. CoA 생산 경로]")
    
    # CoA 생합성 경로
    # Pantothenate -> 4'-Phosphopantothenate -> 4'-Phosphopantothenoylcysteine -> 
    # 4'-Phosphopantetheine -> Dephospho-CoA -> CoA
    
    coa_pathway_reactions = {
        'PNTO': 'Pantothenate kinase',
        'PPCS': '4-phosphopantothenoylcysteine synthetase',
        'PPCDC': '4-phosphopantothenoylcysteine decarboxylase',
        'PPAT': 'Phosphopantetheine adenylyltransferase',
        'DPCOAK': 'Dephospho-CoA kinase'
    }
    
    print("\nCoA 생합성 경로 반응:")
    found_coa = []
    for rxn_id, rxn_name in coa_pathway_reactions.items():
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            found_coa.append(rxn_id)
            print(f"  [OK] {rxn_id}: {rxn.name}")
        except KeyError:
            print(f"  [MISSING] {rxn_id}: {rxn_name}")
    
    # CoA 관련 대사물질 확인
    coa_metabolites = ['pnto_c', '4ppan_c', 'pppcs_c', 'pppan_c', 'dcoa_c', 'coa_c']
    
    print("\nCoA 생합성 중간 대사물질:")
    coa_met_status = []
    
    for met_id in coa_metabolites:
        try:
            met = model.metabolites.get_by_id(met_id)
            reactions = [r.id for r in met.reactions]
            
            # 생산 반응
            producing = [r for r in met.reactions if met in r.products]
            
            coa_met_status.append({
                'Metabolite': met_id,
                'Exists': True,
                'Reaction_Count': len(reactions),
                'Producing_Reactions': len(producing)
            })
            
            status = "[OK]" if producing else "[WARNING]"
            print(f"  {status} {met_id}: {len(reactions)}개 반응, {len(producing)}개 생산 반응")
        except KeyError:
            coa_met_status.append({
                'Metabolite': met_id,
                'Exists': False,
                'Reaction_Count': 0,
                'Producing_Reactions': 0
            })
            print(f"  [MISSING] {met_id}")
    
    # CoA 생산 가능 여부 테스트
    print("\nCoA 생산 가능 여부 테스트:")
    try:
        coa_c = model.metabolites.get_by_id('coa_c')
        
        # 포도당 설정
        for rxn in model.exchanges:
            rxn.upper_bound = 0
            rxn.lower_bound = 0
        
        try:
            ex_glc = model.reactions.get_by_id('EX_glc__D_e')
            ex_glc.lower_bound = -100
            ex_glc.upper_bound = 1000
        except KeyError:
            pass
        
        # 필수 영양소
        essentials = ['EX_nh4_e', 'EX_h2o_e', 'EX_h_e', 'EX_pi_e', 'EX_so4_e',
                      'EX_k_e', 'EX_na1_e', 'EX_mg2_e', 'EX_fe2_e', 'EX_co2_e', 'EX_o2_e']
        
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
        
        # Pantothenate 허용 (비타민 B5)
        try:
            ex_pnto = model.reactions.get_by_id('EX_pnto_e')
            ex_pnto.lower_bound = -1
            ex_pnto.upper_bound = 1000
            print("  [OK] Pantothenate exchange 허용")
        except KeyError:
            print("  [WARNING] EX_pnto_e 없음")
        
        # CoA 생산 테스트
        test_rxn = cobra.Reaction('TEST_coa_c')
        test_rxn.add_metabolites({coa_c: 1})
        test_rxn.lower_bound = 0
        test_rxn.upper_bound = 1000
        
        model.add_reactions([test_rxn])
        model.objective = test_rxn.id
        
        solution = model.optimize()
        
        model.remove_reactions([test_rxn])
        
        can_produce = solution.status == 'optimal' and solution.objective_value > 1e-6
        
        if can_produce:
            print(f"  [SUCCESS] CoA 생산 가능 (최대 플럭스: {solution.objective_value:.6f})")
        else:
            print(f"  [FAIL] CoA 생산 불가 ({solution.status})")
        
        return can_produce, coa_met_status
        
    except KeyError:
        print("  [ERROR] coa_c metabolite 없음")
        return False, coa_met_status

def check_nad_pathway(model):
    """NAD+/NADP+ 생산 경로 확인"""
    print("\n[2. NAD+/NADP+ 생산 경로]")
    
    # NAD+ 생합성 경로
    nad_pathway_reactions = {
        'NADS1': 'NAD synthase (nh3)',
        'NADS2': 'NAD synthase (glutamine)',
        'NADK': 'NAD kinase',
        'NADPPPS': 'NADP phosphatase',
        'NADDPPPS': 'NAD diphosphatase'
    }
    
    print("\nNAD+/NADP+ 생합성 반응:")
    found_nad = []
    for rxn_id, rxn_name in nad_pathway_reactions.items():
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            found_nad.append(rxn_id)
            print(f"  [OK] {rxn_id}: {rxn.name}")
        except KeyError:
            print(f"  [MISSING] {rxn_id}: {rxn_name}")
    
    # NAD+/NADP+ 대사물질 확인
    nad_metabolites = ['nad_c', 'nadh_c', 'nadp_c', 'nadph_c', 'nicotinamide_c', 'nicotinate_c']
    
    print("\nNAD+/NADP+ 관련 대사물질:")
    for met_id in nad_metabolites:
        try:
            met = model.metabolites.get_by_id(met_id)
            reactions = [r.id for r in met.reactions]
            
            producing = [r for r in met.reactions if met in r.products]
            
            status = "[OK]" if producing else "[WARNING]"
            print(f"  {status} {met_id}: {len(reactions)}개 반응, {len(producing)}개 생산 반응")
        except KeyError:
            print(f"  [MISSING] {met_id}")
    
    return found_nad

def check_fad_pathway(model):
    """FAD 생산 경로 확인"""
    print("\n[3. FAD 생산 경로]")
    
    # FAD 생합성 경로
    fad_pathway_reactions = {
        'RBFK': 'Riboflavin kinase',
        'FMNAT': 'FMN adenylyltransferase'
    }
    
    print("\nFAD 생합성 반응:")
    found_fad = []
    for rxn_id, rxn_name in fad_pathway_reactions.items():
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            found_fad.append(rxn_id)
            print(f"  [OK] {rxn_id}: {rxn.name}")
        except KeyError:
            print(f"  [MISSING] {rxn_id}: {rxn_name}")
    
    # FAD 관련 대사물질 확인
    fad_metabolites = ['ribflv_c', 'fmn_c', 'fad_c', 'fadh2_c']
    
    print("\nFAD 관련 대사물질:")
    for met_id in fad_metabolites:
        try:
            met = model.metabolites.get_by_id(met_id)
            reactions = [r.id for r in met.reactions]
            
            producing = [r for r in met.reactions if met in r.products]
            
            status = "[OK]" if producing else "[WARNING]"
            print(f"  {status} {met_id}: {len(reactions)}개 반응, {len(producing)}개 생산 반응")
        except KeyError:
            print(f"  [MISSING] {met_id}")
    
    return found_fad

def check_folate_pathway(model):
    """Folate 관련 경로 확인"""
    print("\n[4. Folate 관련 경로]")
    
    # Folate 관련 대사물질
    folate_metabolites = ['thf_c', 'mlthf_c', '10fthf_c', 'dhf_c']
    
    print("\nFolate 관련 대사물질:")
    folate_status = []
    
    for met_id in folate_metabolites:
        try:
            met = model.metabolites.get_by_id(met_id)
            reactions = [r.id for r in met.reactions]
            
            producing = [r for r in met.reactions if met in r.products]
            
            folate_status.append({
                'Metabolite': met_id,
                'Exists': True,
                'Reaction_Count': len(reactions),
                'Producing_Reactions': len(producing)
            })
            
            status = "[OK]" if producing else "[WARNING]"
            print(f"  {status} {met_id}: {len(reactions)}개 반응, {len(producing)}개 생산 반응")
        except KeyError:
            folate_status.append({
                'Metabolite': met_id,
                'Exists': False,
                'Reaction_Count': 0,
                'Producing_Reactions': 0
            })
            print(f"  [MISSING] {met_id}")
    
    return folate_status

def test_cofactor_production(model):
    """보조인자 생산 가능 여부 테스트"""
    print("\n" + "="*70)
    print("보조인자 생산 가능 여부 테스트")
    print("="*70)
    
    # 포도당 설정
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    try:
        ex_glc = model.reactions.get_by_id('EX_glc__D_e')
        ex_glc.lower_bound = -100
        ex_glc.upper_bound = 1000
    except KeyError:
        pass
    
    # 필수 영양소
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
    
    # 보조인자 생산 테스트
    cofactors = {
        'coa_c': 'CoA',
        'nad_c': 'NAD+',
        'nadp_c': 'NADP+',
        'fad_c': 'FAD',
        'thf_c': 'Tetrahydrofolate'
    }
    
    print("\n보조인자 생산 가능 여부:")
    production_results = []
    
    for met_id, name in cofactors.items():
        try:
            met = model.metabolites.get_by_id(met_id)
            
            # 임시 생산 반응
            test_rxn = cobra.Reaction(f'TEST_{met_id}')
            test_rxn.add_metabolites({met: 1})
            test_rxn.lower_bound = 0
            test_rxn.upper_bound = 1000
            
            model.add_reactions([test_rxn])
            model.objective = test_rxn.id
            
            solution = model.optimize()
            
            model.remove_reactions([test_rxn])
            
            can_produce = solution.status == 'optimal' and solution.objective_value > 1e-6
            
            production_results.append({
                'Cofactor': name,
                'Metabolite_ID': met_id,
                'Can_Produce': can_produce,
                'Status': solution.status,
                'Max_Flux': solution.objective_value if solution.status == 'optimal' else 0
            })
            
            status_icon = "[OK] 생산 가능" if can_produce else "[FAIL] 생산 불가"
            print(f"  {status_icon} {name} ({met_id})")
            
        except KeyError:
            production_results.append({
                'Cofactor': name,
                'Metabolite_ID': met_id,
                'Can_Produce': False,
                'Status': 'Metabolite missing',
                'Max_Flux': 0
            })
            print(f"  [MISSING] {name} ({met_id})")
    
    return production_results

def main():
    print("="*70)
    print("3단계: 보조인자 생산 경로 확인")
    print("="*70)
    
    # 모델 로드
    model = load_model("BaseModel.xml")
    
    # Biomass 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass: {biomass_rxn.id}")
    
    # 1. CoA 경로
    coa_can_produce, coa_met_status = check_coa_pathway(model)
    
    # 2. NAD+/NADP+ 경로
    nad_pathway = check_nad_pathway(model)
    
    # 3. FAD 경로
    fad_pathway = check_fad_pathway(model)
    
    # 4. Folate 경로
    folate_status = check_folate_pathway(model)
    
    # 5. 생산 가능 여부 테스트
    production_results = test_cofactor_production(model)
    
    # 결과 저장
    if coa_met_status:
        df_coa = pd.DataFrame(coa_met_status)
        df_coa.to_csv('coa_pathway_status.csv', index=False)
        print(f"\n[OK] CoA 경로 상태 저장: coa_pathway_status.csv")
    
    if folate_status:
        df_folate = pd.DataFrame(folate_status)
        df_folate.to_csv('folate_pathway_status.csv', index=False)
        print(f"[OK] Folate 경로 상태 저장: folate_pathway_status.csv")
    
    if production_results:
        df_prod = pd.DataFrame(production_results)
        df_prod.to_csv('cofactor_production_status.csv', index=False)
        print(f"[OK] 보조인자 생산 상태 저장: cofactor_production_status.csv")
    
    # 요약
    can_produce_count = sum(1 for s in production_results if s['Can_Produce'])
    total_count = len(production_results)
    
    print("\n" + "="*70)
    print("보조인자 생산 경로 분석 완료")
    print("="*70)
    print(f"\n생산 가능한 보조인자: {can_produce_count}/{total_count}")
    
    if can_produce_count < total_count:
        print("\n[문제 발견] 생산 불가능한 보조인자:")
        for status in production_results:
            if not status['Can_Produce']:
                print(f"  - {status['Cofactor']} ({status['Metabolite_ID']}): {status['Status']}")
    
    print("="*70)

if __name__ == "__main__":
    main()
