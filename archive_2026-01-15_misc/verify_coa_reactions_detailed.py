#!/usr/bin/env python
"""
CoA 경로 반응 상세 확인
실제로 존재하는 반응들이 올바른 대사물질을 사용하는지 확인
"""

import cobra
import pandas as pd

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    return model

def check_coa_reactions_detailed(model):
    """CoA 경로 반응 상세 확인"""
    print("="*70)
    print("CoA 생산 경로 반응 상세 확인")
    print("="*70)
    
    # 발견된 반응들 확인
    reactions_to_check = {
        'PNTK': {
            'expected': 'atp_c + pnto_c --> adp_c + 4ppan_c',
            'found_equation': 'atp_c + pnto__R_c --> 4ppan_c + adp_c + h_c',
            'description': 'Pantothenate kinase'
        },
        'PPCDC': {
            'expected': 'pppcs_c --> co2_c + pppan_c',
            'found_equation': '4ppcys_c + h_c --> co2_c + pan4p_c',
            'description': 'Phosphopantothenoylcysteine decarboxylase'
        },
        'DPCOAK': {
            'expected': 'atp_c + dcoa_c --> adp_c + coa_c',
            'found_equation': 'atp_c + dpcoa_c --> adp_c + coa_c + h_c',
            'description': 'Dephospho-CoA kinase'
        }
    }
    
    print("\n1. PNTO vs PNTK 반응 확인:")
    print("-" * 70)
    
    try:
        pntk = model.reactions.get_by_id('PNTK')
        print(f"[FOUND] PNTK 반응:")
        print(f"  반응식: {pntk.reaction}")
        print(f"  이름: {pntk.name}")
        print(f"  GPR: {pntk.gene_reaction_rule}")
        
        # 대사물질 확인
        pntk_reactants = [str(m) for m, c in pntk.metabolites.items() if c < 0]
        pntk_products = [str(m) for m, c in pntk.metabolites.items() if c > 0]
        
        print(f"  반응물: {pntk_reactants}")
        print(f"  생성물: {pntk_products}")
        
        # pnto_c vs pnto__R_c 확인
        pnto_c_exists = any('pnto_c' in str(m) for m in pntk.metabolites)
        pnto_R_c_exists = any('pnto__R_c' in str(m) for m in pntk.metabolites)
        
        print(f"\n  대사물질 확인:")
        print(f"    pnto_c 사용: {pnto_c_exists}")
        print(f"    pnto__R_c 사용: {pnto_R_c_exists}")
        
        if pnto_R_c_exists:
            try:
                pnto_R_c = model.metabolites.get_by_id('pnto__R_c')
                print(f"    pnto__R_c metabolite 존재: {pnto_R_c.name}")
            except KeyError:
                print(f"    [WARNING] pnto__R_c metabolite 없음")
            
            try:
                pnto_c = model.metabolites.get_by_id('pnto_c')
                print(f"    pnto_c metabolite 존재: {pnto_c.name}")
            except KeyError:
                print(f"    pnto_c metabolite 없음 (pnto__R_c만 존재)")
        
    except KeyError:
        print("[NOT FOUND] PNTK 반응 없음")
    
    print("\n2. PPCDC 반응 확인:")
    print("-" * 70)
    
    try:
        ppcdc = model.reactions.get_by_id('PPCDC')
        print(f"[FOUND] PPCDC 반응:")
        print(f"  반응식: {ppcdc.reaction}")
        print(f"  이름: {ppcdc.name}")
        print(f"  GPR: {ppcdc.gene_reaction_rule}")
        
        # 대사물질 확인
        ppcdc_mets = [str(m) for m in ppcdc.metabolites]
        print(f"  사용 대사물질: {ppcdc_mets}")
        
        # 예상 대사물질 확인
        expected_mets = ['pppcs_c', 'pppan_c', '4ppcys_c', 'pan4p_c']
        for met_id in expected_mets:
            try:
                met = model.metabolites.get_by_id(met_id)
                reactions = [r.id for r in met.reactions]
                print(f"    {met_id}: 존재 ({len(reactions)}개 반응)")
            except KeyError:
                print(f"    {met_id}: 없음")
        
        # 4ppcys_c vs pppcs_c 확인
        if '4ppcys_c' in ppcdc.reaction:
            print(f"\n  [NOTE] 모델은 4ppcys_c 사용 (pppcs_c 아님)")
            print(f"    이 둘이 같은 대사물질일 수 있음 (확인 필요)")
        
    except KeyError:
        print("[NOT FOUND] PPCDC 반응 없음")
    
    print("\n3. DPCOAK 반응 확인:")
    print("-" * 70)
    
    try:
        dpcoak = model.reactions.get_by_id('DPCOAK')
        print(f"[FOUND] DPCOAK 반응:")
        print(f"  반응식: {dpcoak.reaction}")
        print(f"  이름: {dpcoak.name}")
        print(f"  GPR: {dpcoak.gene_reaction_rule}")
        
        # dcoa_c vs dpcoa_c 확인
        dcoa_c_exists = any('dcoa_c' in str(m) for m in dpcoak.metabolites)
        dpcoa_c_exists = any('dpcoa_c' in str(m) for m in dpcoak.metabolites)
        
        print(f"\n  대사물질 확인:")
        print(f"    dcoa_c 사용: {dcoa_c_exists}")
        print(f"    dpcoa_c 사용: {dpcoa_c_exists}")
        
        if dpcoa_c_exists:
            try:
                dpcoa_c = model.metabolites.get_by_id('dpcoa_c')
                print(f"    dpcoa_c metabolite 존재: {dpcoa_c.name}")
            except KeyError:
                print(f"    [WARNING] dpcoa_c metabolite 없음")
            
            try:
                dcoa_c = model.metabolites.get_by_id('dcoa_c')
                print(f"    dcoa_c metabolite 존재: {dcoa_c.name}")
            except KeyError:
                print(f"    dcoa_c metabolite 없음 (dpcoa_c만 존재)")
        
    except KeyError:
        print("[NOT FOUND] DPCOAK 반응 없음")
    
    print("\n4. 누락된 반응 확인:")
    print("-" * 70)
    
    # PPCS 확인
    try:
        ppcs = model.reactions.get_by_id('PPCS')
        print(f"[FOUND] PPCS 반응:")
        print(f"  반응식: {ppcs.reaction}")
        print(f"  이름: {ppcs.name}")
        print(f"  GPR: {ppcs.gene_reaction_rule}")
    except KeyError:
        print("[NOT FOUND] PPCS 반응 없음")
        
        # 유사한 반응 찾기
        print("  유사한 반응 검색 중...")
        for rxn in model.reactions:
            rxn_str = rxn.reaction.lower()
            if '4ppan' in rxn_str or 'pppan' in rxn_str or 'pppcs' in rxn_str or '4ppcys' in rxn_str:
                if 'cys' in rxn_str or 'cysteine' in rxn_str.lower():
                    print(f"    - {rxn.id}: {rxn.name}")
                    print(f"      {rxn.reaction}")
    
    # PPAT 확인
    try:
        ppat = model.reactions.get_by_id('PPAT')
        print(f"\n[FOUND] PPAT 반응:")
        print(f"  반응식: {ppat.reaction}")
        print(f"  이름: {ppat.name}")
        print(f"  GPR: {ppat.gene_reaction_rule}")
    except KeyError:
        print("\n[NOT FOUND] PPAT 반응 없음")
        
        # 유사한 반응 찾기
        print("  유사한 반응 검색 중...")
        for rxn in model.reactions:
            rxn_str = rxn.reaction.lower()
            if ('pppan' in rxn_str or 'pan4p' in rxn_str) and ('dcoa' in rxn_str or 'dpcoa' in rxn_str):
                print(f"    - {rxn.id}: {rxn.name}")
                print(f"      {rxn.reaction}")
    
    print("\n5. CoA 생산 경로 연결성 확인:")
    print("-" * 70)
    
    # 경로상 대사물질들 확인
    coa_pathway_mets = {
        'pnto_c': 'Pantothenate',
        'pnto__R_c': 'Pantothenate (R)',
        '4ppan_c': '4\'-Phosphopantothenate',
        'pppcs_c': '4\'-Phosphopantothenoylcysteine',
        '4ppcys_c': '4\'-Phosphopantothenoylcysteine (alternative)',
        'pppan_c': '4\'-Phosphopantetheine',
        'pan4p_c': '4\'-Phosphopantetheine (alternative)',
        'dcoa_c': 'Dephospho-CoA',
        'dpcoa_c': 'Dephospho-CoA (alternative)',
        'coa_c': 'Coenzyme A'
    }
    
    print("\n경로상 대사물질 존재 여부:")
    met_status = []
    
    for met_id, met_name in coa_pathway_mets.items():
        try:
            met = model.metabolites.get_by_id(met_id)
            reactions = [r.id for r in met.reactions]
            
            # 생산 반응
            producing = [r for r in met.reactions if met in r.products]
            # 소비 반응
            consuming = [r for r in met.reactions if met in r.reactants]
            
            met_status.append({
                'Metabolite_ID': met_id,
                'Name': met_name,
                'Exists': True,
                'Total_Reactions': len(reactions),
                'Producing_Reactions': len(producing),
                'Consuming_Reactions': len(consuming)
            })
            
            print(f"  [OK] {met_id} ({met_name}):")
            print(f"    총 반응: {len(reactions)}개, 생산: {len(producing)}개, 소비: {len(consuming)}개")
            
        except KeyError:
            met_status.append({
                'Metabolite_ID': met_id,
                'Name': met_name,
                'Exists': False,
                'Total_Reactions': 0,
                'Producing_Reactions': 0,
                'Consuming_Reactions': 0
            })
            print(f"  [MISSING] {met_id} ({met_name})")
    
    # 결과 저장
    if met_status:
        df_mets = pd.DataFrame(met_status)
        df_mets.to_csv('coa_pathway_metabolites_status.csv', index=False)
        print(f"\n[OK] CoA 경로 대사물질 상태 저장: coa_pathway_metabolites_status.csv")
    
    print("\n" + "="*70)
    print("결론")
    print("="*70)
    
    print("\n1. PNTO 문제:")
    print("   - PNTK 반응은 존재하지만 pnto__R_c 사용")
    print("   - pnto_c와 pnto__R_c가 같은지 확인 필요")
    print("   - 또는 pnto_c -> pnto__R_c 변환 반응 필요")
    
    print("\n2. PPCS 문제:")
    print("   - PPCS 반응이 없음 (확인 필요)")
    print("   - 4ppan_c + cys__L_c -> pppcs_c 또는 4ppcys_c 반응 필요")
    
    print("\n3. PPCDC:")
    print("   - 반응은 존재하지만 4ppcys_c -> pan4p_c 사용")
    print("   - pppcs_c와 4ppcys_c가 같은지, pppan_c와 pan4p_c가 같은지 확인 필요")
    
    print("\n4. PPAT 문제:")
    print("   - PPAT 반응이 없음")
    print("   - pppan_c (또는 pan4p_c) -> dcoa_c (또는 dpcoa_c) 반응 필요")
    
    print("\n5. DPCOAK:")
    print("   - 반응은 존재하지만 dpcoa_c 사용")
    print("   - dcoa_c와 dpcoa_c가 같은지 확인 필요")

def main():
    model = load_model("BaseModel.xml")
    check_coa_reactions_detailed(model)

if __name__ == "__main__":
    main()
