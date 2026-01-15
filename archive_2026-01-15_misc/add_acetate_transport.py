#!/usr/bin/env python
"""
Acetate transport reaction 추가
ac_e -> ac_c transport 반응 추가
"""

import cobra
from cobra import Reaction, Metabolite

def load_model(model_path="BaseModel.xml"):
    """모델 로드"""
    print(f"모델 로드 중: {model_path}")
    model = cobra.io.read_sbml_model(model_path)
    print(f"[OK] 모델 로드 완료: {model.id}\n")
    return model

def add_acetate_transport(model):
    """
    Acetate transport reaction 추가
    """
    print("="*70)
    print("Acetate Transport Reaction 추가")
    print("="*70)
    
    # ac_e와 ac_c 확인
    try:
        ac_e = model.metabolites.get_by_id('ac_e')
        print(f"[OK] ac_e 존재: {ac_e.id}")
    except KeyError:
        print("[ERROR] ac_e metabolite가 없습니다!")
        return None
    
    try:
        ac_c = model.metabolites.get_by_id('ac_c')
        print(f"[OK] ac_c 존재: {ac_c.id}")
    except KeyError:
        print("[ERROR] ac_c metabolite가 없습니다!")
        return None
    
    # Transport 반응 확인
    transport_ids = ['ACt', 'ACtex', 'ACt2', 'ACt2r', 'ac_transport']
    transport_rxn = None
    
    for trans_id in transport_ids:
        try:
            transport_rxn = model.reactions.get_by_id(trans_id)
            print(f"[OK] {trans_id} transport 반응이 이미 존재합니다: {transport_rxn.reaction}")
            break
        except KeyError:
            continue
    
    if not transport_rxn:
        # 기존 transport 반응 패턴 확인
        for rxn in model.reactions:
            if ac_e in rxn.metabolites and ac_c in rxn.metabolites:
                transport_rxn = rxn
                print(f"[OK] Transport 반응 발견: {rxn.id}: {rxn.reaction}")
                break
        
        if not transport_rxn:
            print("[INFO] Acetate transport 반응이 없습니다. 생성 중...")
            
            # 간단한 diffusion transport 추가: ac_e <=> ac_c
            transport_rxn = Reaction('ACt')
            transport_rxn.name = 'Acetate transport via diffusion'
            transport_rxn.lower_bound = -1000  # 가역적
            transport_rxn.upper_bound = 1000
            
            # ac_e <=> ac_c
            transport_rxn.add_metabolites({
                ac_e: -1,
                ac_c: 1
            })
            
            model.add_reactions([transport_rxn])
            print(f"[OK] ACt transport 반응 생성: {transport_rxn.reaction}")
            print(f"  하한: {transport_rxn.lower_bound}, 상한: {transport_rxn.upper_bound}")
    
    return model

def save_model(model, output_path="BaseModel.xml"):
    """모델 저장"""
    print(f"\n모델 저장 중: {output_path}")
    try:
        cobra.io.write_sbml_model(model, output_path)
        print(f"[OK] 모델 저장 완료: {output_path}")
    except Exception as e:
        print(f"[ERROR] 모델 저장 실패: {e}")
        raise

def main():
    """메인 함수"""
    # 모델 로드
    model = load_model("BaseModel.xml")
    
    # Acetate transport 추가
    model = add_acetate_transport(model)
    
    if model:
        # 모델 저장
        save_model(model, "BaseModel.xml")
        
        print("\n" + "="*70)
        print("Acetate transport reaction 추가 완료!")
        print("="*70)
    else:
        print("\n[ERROR] Transport 반응 추가 실패")

if __name__ == "__main__":
    main()


