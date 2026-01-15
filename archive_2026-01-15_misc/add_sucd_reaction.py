#!/usr/bin/env python
"""
SUCD (Succinate Dehydrogenase) 반응 추가
- 유전자: sdhC (Smlt1796)
- 반응식: succ_c + fad_c --> fum_c + fadh2_c
"""

import cobra
from cobra import Reaction, Metabolite

def load_model(model_path="BaseModel.xml"):
    """모델 로드"""
    print(f"모델 로드 중: {model_path}")
    model = cobra.io.read_sbml_model(model_path)
    print(f"[OK] 모델 로드 완료: {model.id}\n")
    return model

def find_sdh_genes(model):
    """sdh 관련 유전자 찾기"""
    print("="*70)
    print("SDH 관련 유전자 검색")
    print("="*70)
    
    sdh_genes = []
    search_terms = ['Smlt1796', '1796', 'sdh', 'SDH']
    
    for gene in model.genes:
        gene_str = gene.id.upper()
        if any(term.upper() in gene_str for term in search_terms):
            sdh_genes.append(gene)
            print(f"  발견: {gene.id}")
            # 이 유전자와 연결된 반응 확인
            if gene.reactions:
                print(f"    연결된 반응: {len(gene.reactions)}개")
                for rxn in list(gene.reactions)[:3]:
                    print(f"      - {rxn.id}: {rxn.name}")
    
    return sdh_genes

def check_metabolites_exist(model):
    """필수 대사물질 확인"""
    print("\n" + "="*70)
    print("필수 대사물질 확인")
    print("="*70)
    
    metabolites = {
        'succ_c': 'Succinate',
        'fum_c': 'Fumarate',
        'fad_c': 'FAD',
        'fadh2_c': 'FADH2'
    }
    
    found = {}
    missing = []
    
    for met_id, met_name in metabolites.items():
        try:
            met = model.metabolites.get_by_id(met_id)
            found[met_id] = met
            print(f"  [OK] {met_name} ({met_id})")
        except KeyError:
            missing.append(met_id)
            print(f"  [MISSING] {met_name} ({met_id})")
    
    return found, missing

def add_sucd_reaction(model, sdh_genes):
    """SUCD 반응 추가"""
    print("\n" + "="*70)
    print("SUCD 반응 추가")
    print("="*70)
    
    # 이미 SUCD 반응이 있는지 확인
    try:
        existing = model.reactions.get_by_id('SUCD')
        print(f"[INFO] SUCD 반응이 이미 존재합니다: {existing.id}")
        print(f"  반응식: {existing.reaction}")
        print(f"  유전자: {[g.id for g in existing.genes]}")
        return model, existing
    except KeyError:
        pass
    
    # 대사물질 확인
    found, missing = check_metabolites_exist(model)
    
    if missing:
        print(f"\n[ERROR] 필수 대사물질 누락: {', '.join(missing)}")
        print("  SUCD 반응을 추가할 수 없습니다.")
        return model, None
    
    # 대사물질 가져오기
    succ_c = found['succ_c']
    fum_c = found['fum_c']
    fad_c = found['fad_c']
    fadh2_c = found['fadh2_c']
    
    # SUCD 반응 생성
    sucd = Reaction('SUCD')
    sucd.name = 'Succinate dehydrogenase'
    sucd.lower_bound = 0  # 비가역
    sucd.upper_bound = 1000
    
    # 반응식: succ_c + fad_c --> fum_c + fadh2_c
    sucd.add_metabolites({
        succ_c: -1,
        fad_c: -1,
        fum_c: 1,
        fadh2_c: 1
    })
    
    # 유전자 연결
    if sdh_genes:
        # GPR 규칙 설정
        gene_ids = [g.id for g in sdh_genes]
        if len(gene_ids) == 1:
            sucd.gene_reaction_rule = gene_ids[0]
        else:
            # 여러 유전자가 있으면 AND 또는 OR 규칙 설정
            # 일반적으로 SDH는 복합체이므로 AND 규칙
            sucd.gene_reaction_rule = ' and '.join(gene_ids)
        
        print(f"[OK] GPR 규칙 설정: {sucd.gene_reaction_rule}")
    else:
        print("[WARNING] SDH 관련 유전자를 찾을 수 없습니다.")
        print("  유전자 없이 반응을 추가합니다.")
        # 기본 유전자 ID 시도
        try:
            gene = model.genes.get_by_id('Smlt1796')
            sucd.gene_reaction_rule = 'Smlt1796'
            print(f"[OK] Smlt1796 유전자 연결")
        except KeyError:
            try:
                # FAGFNPBA 형식 검색
                for g in model.genes:
                    if '1796' in g.id or ('sdh' in g.id.lower() and 'c' in g.id.lower()):
                        sucd.gene_reaction_rule = g.id
                        print(f"[OK] {g.id} 유전자 연결")
                        break
            except:
                pass
    
    # 모델에 반응 추가
    model.add_reactions([sucd])
    print(f"\n[OK] SUCD 반응 추가 완료")
    print(f"  반응식: {sucd.reaction}")
    print(f"  이름: {sucd.name}")
    print(f"  유전자 규칙: {sucd.gene_reaction_rule}")
    
    return model, sucd

def verify_sucd_reaction(model):
    """SUCD 반응 검증"""
    print("\n" + "="*70)
    print("SUCD 반응 검증")
    print("="*70)
    
    try:
        sucd = model.reactions.get_by_id('SUCD')
        print(f"[OK] SUCD 반응 존재: {sucd.id}")
        print(f"  반응식: {sucd.reaction}")
        print(f"  유전자: {[g.id for g in sucd.genes]}")
        print(f"  GPR: {sucd.gene_reaction_rule}")
        
        # 반응식 검증
        expected_reactants = {'succ_c', 'fad_c'}
        expected_products = {'fum_c', 'fadh2_c'}
        
        reactants = {m.id for m in sucd.reactants}
        products = {m.id for m in sucd.products}
        
        if expected_reactants.issubset(reactants) and expected_products.issubset(products):
            print(f"  [OK] 반응식이 올바릅니다")
        else:
            print(f"  [WARNING] 반응식 검증 실패")
            print(f"    예상 반응물: {expected_reactants}")
            print(f"    실제 반응물: {reactants}")
            print(f"    예상 생성물: {expected_products}")
            print(f"    실제 생성물: {products}")
        
        return True
    except KeyError:
        print(f"[ERROR] SUCD 반응을 찾을 수 없습니다")
        return False

def main():
    model = load_model("BaseModel.xml")
    
    # 1. SDH 관련 유전자 찾기
    sdh_genes = find_sdh_genes(model)
    
    # 2. SUCD 반응 추가
    model, sucd = add_sucd_reaction(model, sdh_genes)
    
    # 3. 검증
    if sucd:
        verify_sucd_reaction(model)
        
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

if __name__ == "__main__":
    main()
