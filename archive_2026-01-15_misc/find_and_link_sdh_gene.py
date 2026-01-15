#!/usr/bin/env python
"""
SDH 유전자 찾기 및 SUCD 반응에 연결
Smlt1796 (sdhC) 유전자를 모델에서 찾아서 SUCD 반응에 연결
"""

import cobra
import re

def load_model(model_path="BaseModel.xml"):
    """모델 로드"""
    print(f"모델 로드 중: {model_path}")
    model = cobra.io.read_sbml_model(model_path)
    print(f"[OK] 모델 로드 완료: {model.id}\n")
    return model

def search_sdh_genes(model):
    """SDH 관련 유전자 찾기"""
    print("="*70)
    print("SDH 관련 유전자 검색 (Smlt1796, sdhC)")
    print("="*70)
    
    # 다양한 검색 패턴
    search_patterns = [
        r'1796',
        r'sdh',
        r'SDH',
        r'Smlt1796',
        r'sdhC',
        r'sdh_c'
    ]
    
    found_genes = []
    
    print("\n전체 유전자 검색 중...")
    for gene in model.genes:
        gene_id = gene.id.upper()
        gene_upper = gene.id
        
        # 각 패턴으로 검색
        for pattern in search_patterns:
            if re.search(pattern, gene_id, re.IGNORECASE):
                if gene not in found_genes:
                    found_genes.append(gene)
                    print(f"\n  발견: {gene.id}")
                    
                    # 유전자와 연결된 반응 확인
                    if gene.reactions:
                        print(f"    연결된 반응 수: {len(gene.reactions)}개")
                        # SDH 관련 반응 찾기
                        sdh_related = []
                        for rxn in gene.reactions:
                            rxn_name_lower = rxn.name.lower()
                            if any(term in rxn_name_lower for term in ['succinate', 'fumarate', 'dehydrogenase', 'sdh']):
                                sdh_related.append(rxn)
                        
                        if sdh_related:
                            print(f"    SDH 관련 반응:")
                            for rxn in sdh_related[:3]:
                                print(f"      - {rxn.id}: {rxn.name}")
                    break
    
    # SUCD 반응 찾기
    print("\n" + "-"*70)
    print("SUCD 반응 확인")
    print("-"*70)
    try:
        sucd = model.reactions.get_by_id('SUCD')
        print(f"[OK] SUCD 반응 존재: {sucd.id}")
        print(f"  반응식: {sucd.reaction}")
        print(f"  현재 GPR: {sucd.gene_reaction_rule}")
        print(f"  현재 연결된 유전자: {[g.id for g in sucd.genes]}")
    except KeyError:
        print("[ERROR] SUCD 반응을 찾을 수 없습니다")
        return None, []
    
    return sucd, found_genes

def link_gene_to_sucd(model, sucd, sdh_genes):
    """SDH 유전자를 SUCD 반응에 연결"""
    print("\n" + "="*70)
    print("SUCD 반응에 유전자 연결")
    print("="*70)
    
    if not sdh_genes:
        print("[WARNING] SDH 관련 유전자를 찾을 수 없습니다.")
        print("  수동으로 Smlt1796 유전자를 추가 시도...")
        
        # 모델에 새로운 유전자 추가 시도
        try:
            # Smlt1796 형식으로 추가 시도
            new_gene = cobra.Gene('Smlt1796')
            new_gene.name = 'succinate dehydrogenase cytochrome b-556 subunit'
            model.genes.append(new_gene)
            sdh_genes = [new_gene]
            print(f"[OK] 새로운 유전자 추가: {new_gene.id}")
        except Exception as e:
            print(f"[ERROR] 유전자 추가 실패: {e}")
            # 모델에 이미 같은 ID가 있을 수 있음
            try:
                existing_gene = model.genes.get_by_id('Smlt1796')
                sdh_genes = [existing_gene]
                print(f"[OK] 기존 유전자 발견: {existing_gene.id}")
            except KeyError:
                print("[ERROR] Smlt1796 유전자를 찾거나 추가할 수 없습니다")
                return False
    
    # 가장 적합한 유전자 선택
    # sdhC는 cytochrome b subunit이므로 이름에 'b' 또는 'cytochrome'이 있는 것 우선
    selected_gene = None
    
    for gene in sdh_genes:
        gene_name_lower = (gene.name or '').lower()
        gene_id_lower = gene.id.lower()
        
        if '1796' in gene_id_lower or 'sdhc' in gene_id_lower or 'cytochrome' in gene_name_lower or 'b' in gene_name_lower:
            selected_gene = gene
            print(f"\n[선택] {gene.id} 유전자 선택")
            if gene.name:
                print(f"  이름: {gene.name}")
            break
    
    if not selected_gene and sdh_genes:
        selected_gene = sdh_genes[0]
        print(f"\n[선택] 첫 번째 발견된 유전자 선택: {selected_gene.id}")
    
    if not selected_gene:
        print("\n[ERROR] 연결할 유전자를 찾을 수 없습니다")
        return False
    
    # GPR 규칙 설정
    if len(sdh_genes) == 1:
        sucd.gene_reaction_rule = selected_gene.id
    else:
        # 여러 유전자가 있으면 AND 규칙 (복합체)
        gene_ids = [g.id for g in sdh_genes]
        sucd.gene_reaction_rule = ' and '.join(gene_ids)
    
    print(f"\n[OK] GPR 규칙 설정 완료")
    print(f"  규칙: {sucd.gene_reaction_rule}")
    print(f"  연결된 유전자: {[g.id for g in sucd.genes]}")
    
    return True

def verify_gene_linkage(model):
    """유전자 연결 검증"""
    print("\n" + "="*70)
    print("유전자 연결 검증")
    print("="*70)
    
    try:
        sucd = model.reactions.get_by_id('SUCD')
        genes = list(sucd.genes)
        
        if genes:
            print(f"[OK] SUCD 반응에 {len(genes)}개 유전자 연결됨")
            for gene in genes:
                print(f"  - {gene.id}")
                if gene.name:
                    print(f"    이름: {gene.name}")
                print(f"    연결된 반응 수: {len(gene.reactions)}개")
            
            print(f"\n  GPR 규칙: {sucd.gene_reaction_rule}")
            return True
        else:
            print("[WARNING] SUCD 반응에 유전자가 연결되어 있지 않습니다")
            return False
            
    except KeyError:
        print("[ERROR] SUCD 반응을 찾을 수 없습니다")
        return False

def main():
    model = load_model("BaseModel.xml")
    
    # 1. SDH 유전자 검색
    sucd, sdh_genes = search_sdh_genes(model)
    
    if not sucd:
        print("\n[ERROR] SUCD 반응을 찾을 수 없습니다. 먼저 add_sucd_reaction.py를 실행하세요.")
        return
    
    # 2. 유전자 연결
    success = link_gene_to_sucd(model, sucd, sdh_genes)
    
    if success:
        # 3. 검증
        verify_gene_linkage(model)
        
        # 4. 모델 저장
        print("\n" + "="*70)
        print("모델 저장")
        print("="*70)
        
        output_path = "BaseModel.xml"
        cobra.io.write_sbml_model(model, output_path)
        print(f"[OK] 모델 저장 완료: {output_path}")
        
        print("\n" + "="*70)
        print("작업 완료")
        print("="*70)
    else:
        print("\n[WARNING] 유전자 연결에 실패했습니다. 수동 확인이 필요할 수 있습니다.")

if __name__ == "__main__":
    main()
