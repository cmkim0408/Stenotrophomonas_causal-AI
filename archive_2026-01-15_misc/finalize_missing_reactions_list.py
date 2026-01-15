#!/usr/bin/env python
"""
누락된 반응 최종 정리
실제로 모델에 없는 반응만 정확히 리스트업
"""

import cobra
import pandas as pd

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    return model

def check_reaction_exists(model, rxn_id):
    """반응이 정확히 존재하는지 확인"""
    try:
        model.reactions.get_by_id(rxn_id)
        return True
    except KeyError:
        return False

def finalize_missing_reactions_list():
    """최종 누락 반응 리스트 정리"""
    
    # 모델 로드
    model = load_model("BaseModel.xml")
    
    # 진단 결과 기반 누락 반응 리스트
    missing_reactions = [
        # CoA 생산 경로
        {
            'Category': 'CoA Pathway',
            'Reaction_ID': 'PNTO',
            'Alternative_IDs': 'PNTOt2rpp (transport만 있음)',
            'Name': 'Pantothenate kinase',
            'Equation': 'atp_c + pnto_c --> adp_c + 4ppan_c',
            'Pathway': 'Pantothenate -> 4\'-Phosphopantothenate',
            'Notes': 'PNTOt2rpp는 transport만 있고 kinase 없음',
            'Priority': 'High'
        },
        {
            'Category': 'CoA Pathway',
            'Reaction_ID': 'PPCS',
            'Alternative_IDs': '없음 (유사: PPCSCT는 다른 반응)',
            'Name': '4-phosphopantothenoylcysteine synthetase',
            'Equation': 'atp_c + 4ppan_c + cys__L_c --> adp_c + pi_c + pppcs_c',
            'Pathway': '4\'-Phosphopantothenate -> 4\'-Phosphopantothenoylcysteine',
            'Notes': 'PPCSCT는 Propanoyl-CoA transferase로 다른 반응',
            'Priority': 'High'
        },
        {
            'Category': 'CoA Pathway',
            'Reaction_ID': 'PPAT',
            'Alternative_IDs': 'APPAT (유사하지만 pan4p_c 사용)',
            'Name': 'Phosphopantetheine adenylyltransferase',
            'Equation': 'atp_c + pppan_c --> ppi_c + dcoa_c',
            'Pathway': '4\'-Phosphopantetheine -> Dephospho-CoA',
            'Notes': 'APPAT는 pan4p_c를 사용하는 다른 반응',
            'Priority': 'High'
        },
        # 아미노산 생합성 경로
        {
            'Category': 'Amino Acid (Serine)',
            'Reaction_ID': 'SER',
            'Alternative_IDs': '없음 (PSP_L는 pser__L_c -> ser__L_c 존재)',
            'Name': 'Serine synthase / Phosphoserine phosphatase',
            'Equation': 'pser__L_c + h2o_c --> ser__L_c + pi_c',
            'Pathway': 'Serine pathway (Phosphoserine -> Serine)',
            'Notes': 'PSP_L 반응 확인 필요 - 실제로 존재할 수 있음',
            'Priority': 'Medium'
        },
        {
            'Category': 'Amino Acid (Serine)',
            'Reaction_ID': 'SERD',
            'Alternative_IDs': 'SERD_L (존재할 수 있음)',
            'Name': 'Serine dehydratase',
            'Equation': 'ser__L_c --> nh3_c + pyr_c + h2o_c',
            'Pathway': 'Serine pathway (Serine -> Pyruvate)',
            'Notes': 'SERD_L이 같은 반응일 수 있음 - 확인 필요',
            'Priority': 'Low'
        },
        {
            'Category': 'Amino Acid (Glycine)',
            'Reaction_ID': 'GCY',
            'Alternative_IDs': '없음',
            'Name': 'Glycine synthase / Glycine cleavage system',
            'Equation': 'co2_c + nh4_c + 5fthf_c + nad_c --> gly_c + nadh_c + mlthf_c',
            'Pathway': 'Glycine pathway (CO2 + NH4 -> Glycine)',
            'Notes': 'Glycine 생합성 경로 확인 필요',
            'Priority': 'Medium'
        },
        {
            'Category': 'Amino Acid (Glycine)',
            'Reaction_ID': 'SHMT',
            'Alternative_IDs': '없음',
            'Name': 'Serine hydroxymethyltransferase',
            'Equation': 'ser__L_c + thf_c --> gly_c + mlthf_c + h2o_c',
            'Pathway': 'Glycine pathway (Serine -> Glycine)',
            'Notes': 'Glycine 생합성 주요 경로',
            'Priority': 'High'
        },
        {
            'Category': 'Amino Acid (Aromatic)',
            'Reaction_ID': 'DAHPS',
            'Alternative_IDs': '없음 (AROAT, AROH는 다른 단계)',
            'Name': '3-deoxy-7-phosphoheptulonate synthase',
            'Equation': 'e4p_c + pep_c --> dah7p_c + pi_c + h2o_c',
            'Pathway': 'Aromatic AA (Erythrose-4-P + PEP -> DAHP)',
            'Notes': 'Aromatic AA 생합성 첫 단계',
            'Priority': 'High'
        },
        {
            'Category': 'Amino Acid (Aromatic)',
            'Reaction_ID': 'TRPS',
            'Alternative_IDs': 'TRPS1, TRPS2 (존재함)',
            'Name': 'Tryptophan synthase',
            'Equation': 'indole_c + ser__L_c --> trp__L_c + h2o_c',
            'Pathway': 'Aromatic AA (Indole + Serine -> Tryptophan)',
            'Notes': 'TRPS1, TRPS2로 이미 존재 - 문제 없음',
            'Priority': 'None'
        },
        {
            'Category': 'Amino Acid (Aromatic)',
            'Reaction_ID': 'TYRS',
            'Alternative_IDs': '없음',
            'Name': 'Tyrosine synthase',
            'Equation': 'prephenate_c --> tyr__L_c + co2_c + h2o_c',
            'Pathway': 'Aromatic AA (Prephenate -> Tyrosine)',
            'Notes': 'AROH 경로 확인 필요',
            'Priority': 'Medium'
        },
        {
            'Category': 'Amino Acid (Aromatic)',
            'Reaction_ID': 'PHES',
            'Alternative_IDs': '없음',
            'Name': 'Phenylalanine synthase',
            'Equation': 'prephenate_c --> phe__L_c + co2_c + h2o_c',
            'Pathway': 'Aromatic AA (Prephenate -> Phenylalanine)',
            'Notes': 'AROH 경로 확인 필요',
            'Priority': 'Medium'
        },
        {
            'Category': 'Amino Acid (Aspartate family)',
            'Reaction_ID': 'HOM',
            'Alternative_IDs': '없음',
            'Name': 'Homoserine dehydrogenase',
            'Equation': 'aspartate_semialdehyde_c + nadph_c + h_c --> hom__L_c + nadp_c + h2o_c',
            'Pathway': 'Aspartate family (Aspartate semialdehyde -> Homoserine)',
            'Notes': 'Homoserine 생합성 경로',
            'Priority': 'High'
        },
        {
            'Category': 'Amino Acid (Aspartate family)',
            'Reaction_ID': 'LYSS',
            'Alternative_IDs': '없음',
            'Name': 'Lysine synthase / Diaminopimelate pathway',
            'Equation': 'meso_dap_c --> lys__L_c + succ_c',
            'Pathway': 'Aspartate family (meso-Diaminopimelate -> Lysine)',
            'Notes': 'Lysine 생합성 경로',
            'Priority': 'High'
        },
        # 뉴클레오티드 생합성
        {
            'Category': 'Nucleotide',
            'Reaction_ID': 'PRPPS',
            'Alternative_IDs': 'PRPPS (존재함)',
            'Name': 'Phosphoribosyl pyrophosphate synthetase',
            'Equation': 'atp_c + r5p_c --> prpp_c + amp_c',
            'Pathway': 'Purine/Pyrimidine synthesis (R5P -> PRPP)',
            'Notes': '이미 존재 - 문제 없음',
            'Priority': 'None'
        }
    ]
    
    # 실제로 존재하는지 재확인
    print("="*70)
    print("누락된 반응 최종 확인")
    print("="*70)
    
    final_missing = []
    exists_list = []
    
    for rxn_info in missing_reactions:
        rxn_id = rxn_info['Reaction_ID']
        exists = check_reaction_exists(model, rxn_id)
        
        if exists:
            exists_list.append(rxn_info)
            print(f"\n[EXISTS] {rxn_id}: {rxn_info['Name']}")
        else:
            final_missing.append(rxn_info)
            print(f"\n[MISSING] {rxn_id}: {rxn_info['Name']}")
            print(f"  카테고리: {rxn_info['Category']}")
            print(f"  경로: {rxn_info['Pathway']}")
            print(f"  반응식: {rxn_info['Equation']}")
            if rxn_info['Alternative_IDs']:
                print(f"  대체 ID: {rxn_info['Alternative_IDs']}")
            if rxn_info['Notes']:
                print(f"  참고: {rxn_info['Notes']}")
    
    # 우선순위별 정리
    high_priority = [r for r in final_missing if r['Priority'] == 'High']
    medium_priority = [r for r in final_missing if r['Priority'] == 'Medium']
    low_priority = [r for r in final_missing if r['Priority'] == 'Low']
    
    print("\n" + "="*70)
    print("최종 누락 반응 요약")
    print("="*70)
    
    print(f"\n총 누락 반응: {len(final_missing)}개")
    print(f"  - High Priority: {len(high_priority)}개")
    print(f"  - Medium Priority: {len(medium_priority)}개")
    print(f"  - Low Priority: {len(low_priority)}개")
    
    print(f"\n이미 존재하는 반응: {len(exists_list)}개")
    for rxn in exists_list:
        print(f"  - {rxn['Reaction_ID']}: {rxn['Name']}")
    
    # 카테고리별 정리
    print("\n카테고리별 누락 반응:")
    by_category = {}
    for rxn in final_missing:
        category = rxn['Category']
        if category not in by_category:
            by_category[category] = []
        by_category[category].append(rxn)
    
    for category, rxns in by_category.items():
        print(f"\n  [{category}]: {len(rxns)}개")
        for rxn in rxns:
            priority_icon = '***' if rxn['Priority'] == 'High' else '**' if rxn['Priority'] == 'Medium' else '*'
            print(f"    {priority_icon} {rxn['Reaction_ID']}: {rxn['Name']}")
    
    # 결과 저장
    if final_missing:
        df_final = pd.DataFrame(final_missing)
        df_final.to_csv('final_missing_reactions_list.csv', index=False)
        print(f"\n[OK] 최종 누락 반응 리스트 저장: final_missing_reactions_list.csv")
        
        # 우선순위별로 별도 저장
        if high_priority:
            df_high = pd.DataFrame(high_priority)
            df_high.to_csv('missing_reactions_high_priority.csv', index=False)
            print(f"[OK] High Priority 누락 반응 저장: missing_reactions_high_priority.csv")
    
    return final_missing, exists_list

def main():
    final_missing, exists_list = finalize_missing_reactions_list()
    
    print("\n" + "="*70)
    print("최종 정리 완료")
    print("="*70)
    
    print("\n[다음 단계]")
    print("  1. High Priority 반응부터 추가 검토")
    print("  2. 실제로 대체 반응이 있는지 확인")
    print("  3. 필요한 반응만 Gap-filling 수행")

if __name__ == "__main__":
    main()
