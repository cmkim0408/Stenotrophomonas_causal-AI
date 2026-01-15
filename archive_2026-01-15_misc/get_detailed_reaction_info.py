#!/usr/bin/env python
"""
레퍼런스 모델에서 누락된 반응들의 상세 정보 추출 (이름, 유전자 등)
"""

import cobra
import pandas as pd
from pathlib import Path

def load_model(model_path):
    print(f"모델 로드 중: {model_path}")
    model = cobra.io.read_sbml_model(model_path)
    print(f"[OK] 모델 로드 완료: {model.id}")
    return model

def get_reaction_full_info(model, rxn_id):
    """반응 전체 정보 가져오기"""
    try:
        rxn = model.reactions.get_by_id(rxn_id)
        genes = [g.id for g in rxn.genes]
        gpr = str(rxn.gene_reaction_rule) if hasattr(rxn, 'gene_reaction_rule') else ''
        
        # Pseudo reaction 확인 (SUPPLY_, DM_ 등으로 시작하거나 특정 패턴)
        is_pseudo = False
        if rxn_id.startswith('SUPPLY_') or rxn_id.startswith('DM_') or rxn_id.startswith('BIOMASS'):
            is_pseudo = True
        # 메타볼라이트가 하나만 있고 계수가 1인 경우 (공급 반응)
        if len(rxn.metabolites) == 1:
            met = list(rxn.metabolites.keys())[0]
            coeff = list(rxn.metabolites.values())[0]
            if coeff > 0 and len(rxn.reactants) == 0:
                is_pseudo = True
        
        return {
            'exists': True,
            'id': rxn.id,
            'name': rxn.name if rxn.name else '',
            'equation': rxn.reaction,
            'genes': ', '.join(genes) if genes else '',
            'gpr': gpr,
            'gene_count': len(genes),
            'subsystem': rxn.subsystem if hasattr(rxn, 'subsystem') else '',
            'is_pseudo': is_pseudo
        }
    except KeyError:
        return {
            'exists': False,
            'id': rxn_id,
            'name': '',
            'equation': '',
            'genes': '',
            'gpr': '',
            'gene_count': 0,
            'subsystem': '',
            'is_pseudo': False
        }

def main():
    base_path = Path(__file__).parent.parent
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    input_csv = base_path / "Stenotrophomonas-causal AI" / "missing_metabolic_reactions_only.csv"
    output_csv = base_path / "Stenotrophomonas-causal AI" / "missing_metabolic_reactions_detailed.csv"
    output_md = base_path / "Stenotrophomonas-causal AI" / "missing_metabolic_reactions_detailed.md"
    
    # 모델 로드
    ref_model = load_model(str(ref_model_path))
    
    # CSV 읽기
    df = pd.read_csv(input_csv)
    
    # 상세 정보 추가
    detailed_info = []
    
    for _, row in df.iterrows():
        rxn_id = row['reaction_id']
        info = get_reaction_full_info(ref_model, rxn_id)
        
        detailed_info.append({
            'reaction_id': rxn_id,
            'name': info['name'],
            'equation': info['equation'],
            'pathway': row['pathway'],
            'priority': row['priority'],
            'is_active': row['is_active'],
            'max_flux': row['max_flux'] if row['is_active'] else 0.0,
            'genes': info['genes'],
            'gene_reaction_rule': info['gpr'],
            'gene_count': info['gene_count'],
            'subsystem': info['subsystem'],
            'is_pseudo': info['is_pseudo'],
            'original_name': row.get('name', '')  # CSV에 있던 이름
        })
    
    # DataFrame 생성
    df_detailed = pd.DataFrame(detailed_info)
    
    # 정렬
    priority_order = {'HIGH': 1, 'MEDIUM': 2, 'LOW': 3}
    df_detailed['priority_num'] = df_detailed['priority'].map(priority_order)
    df_detailed = df_detailed.sort_values(['priority_num', 'pathway', 'reaction_id'])
    df_detailed = df_detailed.drop('priority_num', axis=1)
    
    # CSV 저장
    df_detailed.to_csv(output_csv, index=False, encoding='utf-8-sig')
    print(f"\n[OK] CSV 저장: {output_csv}")
    
    # Markdown 파일 생성
    with open(output_md, 'w', encoding='utf-8') as f:
        f.write("# 레퍼런스 모델에 있지만 신규 모델에 없는 대사 경로 반응 리스트\n")
        f.write("## (Transport 및 Exchange 제외)\n\n")
        f.write(f"총 **{len(df_detailed)}개** 반응\n\n")
        f.write("---\n\n")
        
        # Pseudo reaction 제외한 실제 반응
        df_real = df_detailed[~df_detailed['is_pseudo']]
        f.write(f"## 실제 대사 반응: {len(df_real)}개\n\n")
        f.write(f"### Pseudo 반응 (제외): {len(df_detailed[df_detailed['is_pseudo']])}개\n\n")
        
        for _, row in df_detailed[df_detailed['is_pseudo']].iterrows():
            f.write(f"- **{row['reaction_id']}**: {row['equation']} (Pseudo reaction)\n")
        
        f.write("\n---\n\n")
        
        # HIGH 우선순위
        df_high = df_real[df_real['priority'] == 'HIGH']
        if len(df_high) > 0:
            f.write("## HIGH 우선순위 반응 (실제 사용됨)\n\n")
            for i, (_, row) in enumerate(df_high.iterrows(), 1):
                f.write(f"### {i}. {row['reaction_id']}\n")
                if row['name']:
                    f.write(f"- **이름**: {row['name']}\n")
                f.write(f"- **반응식**: `{row['equation']}`\n")
                f.write(f"- **경로**: {row['pathway']}\n")
                if row['is_active']:
                    f.write(f"- **실제 사용됨**: Yes (플럭스: {row['max_flux']:.6f})\n")
                if row['genes']:
                    f.write(f"- **유전자**: {row['genes']}\n")
                if row['gene_reaction_rule']:
                    f.write(f"- **GPR**: {row['gene_reaction_rule']}\n")
                f.write("\n")
        
        # MEDIUM 우선순위
        df_medium = df_real[df_real['priority'] == 'MEDIUM']
        if len(df_medium) > 0:
            f.write("## MEDIUM 우선순위 반응\n\n")
            for i, (_, row) in enumerate(df_medium.iterrows(), 1):
                f.write(f"### {i}. {row['reaction_id']}\n")
                if row['name']:
                    f.write(f"- **이름**: {row['name']}\n")
                f.write(f"- **반응식**: `{row['equation']}`\n")
                f.write(f"- **경로**: {row['pathway']}\n")
                if row['genes']:
                    f.write(f"- **유전자**: {row['genes']}\n")
                if row['gene_reaction_rule']:
                    f.write(f"- **GPR**: {row['gene_reaction_rule']}\n")
                f.write("\n")
        
        # LOW 우선순위
        df_low = df_real[df_real['priority'] == 'LOW']
        if len(df_low) > 0:
            f.write("## LOW 우선순위 반응\n\n")
            for i, (_, row) in enumerate(df_low.iterrows(), 1):
                f.write(f"### {i}. {row['reaction_id']}\n")
                if row['name']:
                    f.write(f"- **이름**: {row['name']}\n")
                f.write(f"- **반응식**: `{row['equation']}`\n")
                f.write(f"- **경로**: {row['pathway']}\n")
                if row['genes']:
                    f.write(f"- **유전자**: {row['genes']}\n")
                if row['gene_reaction_rule']:
                    f.write(f"- **GPR**: {row['gene_reaction_rule']}\n")
                f.write("\n")
    
    print(f"[OK] Markdown 저장: {output_md}")
    
    # 요약 출력
    print("\n" + "="*70)
    print("요약")
    print("="*70)
    
    total = len(df_detailed)
    pseudo = len(df_detailed[df_detailed['is_pseudo']])
    real = total - pseudo
    with_genes = len(df_detailed[df_detailed['genes'] != ''])
    
    print(f"\n총 {total}개 반응")
    print(f"  실제 대사 반응: {real}개")
    print(f"  Pseudo 반응: {pseudo}개 (제외 대상)")
    print(f"  유전자 정보 있음: {with_genes}개")
    
    print("\nPseudo 반응 목록:")
    for _, row in df_detailed[df_detailed['is_pseudo']].iterrows():
        print(f"  - {row['reaction_id']}: {row['equation']}")

if __name__ == "__main__":
    main()
