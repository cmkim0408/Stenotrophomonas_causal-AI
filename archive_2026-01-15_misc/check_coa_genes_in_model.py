#!/usr/bin/env python
"""
CoA 생산 경로 유전자 확인
KEGG에서 확인된 유전자들이 모델에 연결되어 있는지 확인
"""

import cobra
import pandas as pd

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    print(f"[OK] 모델 로드: {model.id}")
    return model

def find_gene_by_pattern(model, patterns):
    """패턴으로 유전자 찾기"""
    found_genes = []
    for pattern in patterns:
        for gene in model.genes:
            if pattern.lower() in gene.id.lower() or pattern.lower() in gene.name.lower():
                if gene not in found_genes:
                    found_genes.append(gene)
    return found_genes

def check_coa_genes_in_model(model):
    """CoA 생산 경로 유전자 확인"""
    print("="*70)
    print("CoA 생산 경로 유전자 확인")
    print("="*70)
    
    # KEGG에서 확인된 유전자들
    coa_genes = [
        {
            'gene_id': 'Smlt0283',
            'name': 'pantothenate kinase',
            'ec': 'EC:2.7.1.33',
            'reaction': 'PNTO',
            'enzyme_name': 'pantothenate kinase',
            'expected_reaction': 'atp_c + pnto_c --> adp_c + 4ppan_c'
        },
        {
            'gene_id': 'Smlt0401',
            'name': 'dfp (pantothenate biosynthesis protein)',
            'ec': 'EC:4.1.1.36 6.3.2.5',
            'reaction': 'PPCS, PPCDC',
            'enzyme_name': 'phosphopantothenate---cysteine ligase / phosphopantothenoylcysteine decarboxylase',
            'expected_reaction_1': 'atp_c + 4ppan_c + cys__L_c --> adp_c + pi_c + pppcs_c',
            'expected_reaction_2': 'pppcs_c --> co2_c + pppan_c'
        },
        {
            'gene_id': 'Smlt1811',
            'name': 'coaD (phosphopantetheine adenylyltransferase)',
            'ec': 'EC:2.7.7.3',
            'reaction': 'PPAT',
            'enzyme_name': 'pantetheine-phosphate adenylyltransferase',
            'expected_reaction': 'atp_c + pppan_c --> ppi_c + dcoa_c'
        },
        {
            'gene_id': 'Smlt3761',
            'name': 'coaE (dephospho-CoA kinase)',
            'ec': 'EC:2.7.1.24',
            'reaction': 'DPCOAK',
            'enzyme_name': 'dephospho-CoA kinase',
            'expected_reaction': 'atp_c + dcoa_c --> adp_c + coa_c'
        }
    ]
    
    print("\nKEGG에서 확인된 유전자들이 모델에 있는지 확인:")
    print("="*70)
    
    gene_status = []
    
    for gene_info in coa_genes:
        gene_id = gene_info['gene_id']
        print(f"\n[유전자] {gene_id}: {gene_info['name']}")
        print(f"  EC: {gene_info['ec']}")
        print(f"  예상 반응: {gene_info.get('reaction', 'N/A')}")
        
        # 유전자 찾기
        try:
            gene = model.genes.get_by_id(gene_id)
            print(f"  상태: [EXISTS] 모델에 존재")
            print(f"    이름: {gene.name}")
            
            # 이 유전자가 연결된 반응 확인
            reactions = gene.reactions
            print(f"    연결된 반응 수: {len(reactions)}개")
            
            if reactions:
                print(f"    연결된 반응:")
                for rxn in reactions:
                    print(f"      - {rxn.id}: {rxn.name}")
                    print(f"        {rxn.reaction}")
                    
                    # GPR 확인
                    gpr = rxn.gene_reaction_rule
                    print(f"        GPR: {gpr}")
            else:
                print(f"    [WARNING] 반응에 연결되지 않음!")
            
            gene_status.append({
                'Gene_ID': gene_id,
                'Name': gene_info['name'],
                'EC': gene_info['ec'],
                'Expected_Reaction': gene_info['reaction'],
                'Exists': True,
                'Connected_Reactions': len(reactions),
                'Reaction_IDs': ', '.join([r.id for r in reactions])
            })
            
        except KeyError:
            # 패턴으로 찾기
            patterns = [gene_id.lower(), gene_info['gene_id'].replace('Smlt', '').lower()]
            found = find_gene_by_pattern(model, patterns)
            
            if found:
                print(f"  상태: [SIMILAR] 유사한 유전자 발견")
                for g in found:
                    print(f"    - {g.id}: {g.name}")
            else:
                print(f"  상태: [MISSING] 모델에 없음")
            
            gene_status.append({
                'Gene_ID': gene_id,
                'Name': gene_info['name'],
                'EC': gene_info['ec'],
                'Expected_Reaction': gene_info['reaction'],
                'Exists': False,
                'Connected_Reactions': 0,
                'Reaction_IDs': ''
            })
    
    # 반응 존재 여부 확인
    print("\n" + "="*70)
    print("예상 반응들이 모델에 있는지 확인:")
    print("="*70)
    
    expected_reactions = {
        'PNTO': {
            'equation': 'atp_c + pnto_c --> adp_c + 4ppan_c',
            'gene': 'Smlt0283'
        },
        'PPCS': {
            'equation': 'atp_c + 4ppan_c + cys__L_c --> adp_c + pi_c + pppcs_c',
            'gene': 'Smlt0401'
        },
        'PPCDC': {
            'equation': 'pppcs_c --> co2_c + pppan_c',
            'gene': 'Smlt0401'
        },
        'PPAT': {
            'equation': 'atp_c + pppan_c --> ppi_c + dcoa_c',
            'gene': 'Smlt1811'
        },
        'DPCOAK': {
            'equation': 'atp_c + dcoa_c --> adp_c + coa_c',
            'gene': 'Smlt3761'
        }
    }
    
    reaction_status = []
    
    for rxn_id, rxn_info in expected_reactions.items():
        print(f"\n[반응] {rxn_id}")
        print(f"  반응식: {rxn_info['equation']}")
        print(f"  예상 유전자: {rxn_info['gene']}")
        
        # 반응 찾기
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            print(f"  상태: [EXISTS] {rxn_id}")
            print(f"    이름: {rxn.name}")
            print(f"    실제 반응식: {rxn.reaction}")
            print(f"    GPR: {rxn.gene_reaction_rule}")
            
            # 유전자 확인
            if rxn_info['gene'] in rxn.gene_reaction_rule:
                print(f"    [OK] 예상 유전자 {rxn_info['gene']}가 GPR에 포함됨")
            else:
                print(f"    [WARNING] 예상 유전자 {rxn_info['gene']}가 GPR에 없음")
                print(f"      실제 GPR: {rxn.gene_reaction_rule}")
            
            reaction_status.append({
                'Reaction_ID': rxn_id,
                'Expected_Equation': rxn_info['equation'],
                'Expected_Gene': rxn_info['gene'],
                'Exists': True,
                'Actual_Equation': rxn.reaction,
                'GPR': rxn.gene_reaction_rule,
                'Gene_Match': rxn_info['gene'] in rxn.gene_reaction_rule
            })
            
        except KeyError:
            # 유사한 반응 찾기
            print(f"  상태: [MISSING] 반응 없음")
            
            # 반응식으로 유사한 반응 찾기
            similar_found = []
            reactants = rxn_info['equation'].split(' --> ')[0].split(' + ')
            products = rxn_info['equation'].split(' --> ')[1].split(' + ')
            
            # 간단한 패턴 매칭
            for rxn in model.reactions:
                rxn_str = rxn.reaction.lower()
                # 주요 대사물질 포함 여부 확인
                key_mets = ['pnto', '4ppan', 'pppcs', 'pppan', 'dcoa', 'coa']
                if any(met in rxn_str for met in key_mets if met in rxn_id.lower()):
                    if rxn not in similar_found:
                        similar_found.append(rxn)
            
            if similar_found:
                print(f"    유사한 반응:")
                for rxn in similar_found[:3]:
                    print(f"      - {rxn.id}: {rxn.name}")
                    print(f"        {rxn.reaction}")
            
            reaction_status.append({
                'Reaction_ID': rxn_id,
                'Expected_Equation': rxn_info['equation'],
                'Expected_Gene': rxn_info['gene'],
                'Exists': False,
                'Actual_Equation': '',
                'GPR': '',
                'Gene_Match': False
            })
    
    # 결과 저장
    if gene_status:
        df_genes = pd.DataFrame(gene_status)
        df_genes.to_csv('coa_genes_status.csv', index=False)
        print(f"\n[OK] 유전자 상태 저장: coa_genes_status.csv")
    
    if reaction_status:
        df_rxns = pd.DataFrame(reaction_status)
        df_rxns.to_csv('coa_reactions_status.csv', index=False)
        print(f"[OK] 반응 상태 저장: coa_reactions_status.csv")
    
    # 요약
    print("\n" + "="*70)
    print("요약")
    print("="*70)
    
    existing_genes = sum(1 for g in gene_status if g['Exists'])
    existing_reactions = sum(1 for r in reaction_status if r['Exists'])
    
    print(f"\n유전자:")
    print(f"  존재: {existing_genes}/{len(gene_status)}개")
    print(f"  누락: {len(gene_status) - existing_genes}개")
    
    print(f"\n반응:")
    print(f"  존재: {existing_reactions}/{len(reaction_status)}개")
    print(f"  누락: {len(reaction_status) - existing_reactions}개")
    
    # 유전자는 있지만 반응에 연결되지 않은 경우
    print(f"\n유전자-반응 연결:")
    for g in gene_status:
        if g['Exists'] and g['Connected_Reactions'] == 0:
            print(f"  [WARNING] {g['Gene_ID']}: 반응에 연결되지 않음")
    
    # 반응은 있지만 예상 유전자가 없는 경우
    for r in reaction_status:
        if r['Exists'] and not r['Gene_Match']:
            print(f"  [WARNING] {r['Reaction_ID']}: 예상 유전자 {r['Expected_Gene']}가 GPR에 없음")
    
    print("\n" + "="*70)

def main():
    model = load_model("BaseModel.xml")
    
    print(f"\n총 유전자 수: {len(model.genes)}")
    print(f"총 반응 수: {len(model.reactions)}")
    
    check_coa_genes_in_model(model)

if __name__ == "__main__":
    main()
