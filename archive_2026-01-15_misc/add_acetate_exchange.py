#!/usr/bin/env python
"""
Acetate exchange reaction 추가 스크립트
"""

import cobra
from cobra import Reaction, Metabolite

def load_model(model_path="BaseModel.xml"):
    """모델 로드"""
    print(f"모델 로드 중: {model_path}")
    model = cobra.io.read_sbml_model(model_path)
    print(f"[OK] 모델 로드 완료: {model.id}\n")
    return model

def add_acetate_exchange(model):
    """
    Acetate exchange reaction 추가
    """
    print("="*70)
    print("Acetate Exchange Reaction 추가")
    print("="*70)
    
    # ac_e metabolite 확인
    try:
        ac_e = model.metabolites.get_by_id('ac_e')
        print(f"[OK] ac_e metabolite 존재: {ac_e.id}")
    except KeyError:
        print("[INFO] ac_e metabolite가 없습니다. 생성 중...")
        # ac_c에서 정보 가져오기
        ac_c = model.metabolites.get_by_id('ac_c')
        
        # ac_e 생성 (extracellular compartment)
        ac_e = Metabolite(
            id='ac_e',
            formula=ac_c.formula,
            name='Acetate',
            compartment='C_e',
            charge=ac_c.charge
        )
        model.add_metabolites([ac_e])
        print(f"[OK] ac_e metabolite 생성: {ac_e.id}")
    
    # EX_ac_e 반응 확인
    try:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        print(f"[OK] EX_ac_e 반응이 이미 존재합니다: {ex_ac.id}")
        print(f"  반응식: {ex_ac.reaction}")
        print(f"  하한: {ex_ac.lower_bound}, 상한: {ex_ac.upper_bound}")
    except KeyError:
        print("[INFO] EX_ac_e 반응이 없습니다. 생성 중...")
        
        # Exchange reaction 생성
        ex_ac = Reaction('EX_ac_e')
        ex_ac.name = 'Acetate exchange'
        ex_ac.lower_bound = -1000  # Uptake 허용
        ex_ac.upper_bound = 1000   # Secretion 허용
        
        # ac_e <=> (exchange reaction)
        ex_ac.add_metabolites({ac_e: -1})
        
        model.add_reactions([ex_ac])
        print(f"[OK] EX_ac_e 반응 생성: {ex_ac.id}")
        print(f"  반응식: {ex_ac.reaction}")
        print(f"  하한: {ex_ac.lower_bound}, 상한: {ex_ac.upper_bound}")
    
    return model

def save_model(model, output_path="BaseModel_with_acetate.xml"):
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
    
    # Acetate exchange 추가
    model = add_acetate_exchange(model)
    
    # 모델 저장 (BaseModel.xml 업데이트)
    save_model(model, "BaseModel.xml")
    
    print("\n" + "="*70)
    print("Acetate exchange reaction 추가 완료!")
    print("="*70)

if __name__ == "__main__":
    main()

