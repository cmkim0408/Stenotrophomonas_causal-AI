#!/usr/bin/env python
"""
Genome-Scale Metabolic Model (GSM) 구축 스크립트
BaseModel.xml을 기반으로 GSM을 로드하고 분석합니다.
"""

import cobra
from cobra import Model
import pandas as pd
import numpy as np

def load_model(model_path="BaseModel.xml"):
    """
    SBML 파일에서 대사 모델을 로드합니다.
    
    Parameters:
    -----------
    model_path : str
        SBML 모델 파일 경로
    
    Returns:
    --------
    model : cobra.Model
        로드된 대사 모델
    """
    print(f"모델 로드 중: {model_path}")
    try:
        model = cobra.io.read_sbml_model(model_path)
        print(f"[OK] 모델 로드 완료: {model.id}")
        return model
    except Exception as e:
        print(f"[ERROR] 모델 로드 실패: {e}")
        raise

def print_model_summary(model):
    """
    모델의 기본 정보를 출력합니다.
    
    Parameters:
    -----------
    model : cobra.Model
        대사 모델
    """
    print("\n" + "="*60)
    print("GSM 모델 요약 정보")
    print("="*60)
    print(f"모델 ID: {model.id}")
    print(f"모델 이름: {model.name}")
    print(f"\n반응(Reactions) 수: {len(model.reactions)}")
    print(f"대사물질(Metabolites) 수: {len(model.metabolites)}")
    print(f"유전자(Genes) 수: {len(model.genes)}")
    print(f"구획(Compartments) 수: {len(model.compartments)}")
    
    # 구획 정보
    print(f"\n구획 목록: {list(model.compartments.keys())}")
    
    # 경계 조건이 있는 대사물질 수 (exchange reaction과 연결된 대사물질)
    exchange_reactions = [r for r in model.reactions if r.id.startswith('EX_')]
    print(f"교환 반응(Exchange reactions) 수: {len(exchange_reactions)}")
    
    # 반응 타입별 통계
    reversible_reactions = [r for r in model.reactions if r.reversibility]
    print(f"가역 반응(Reversible) 수: {len(reversible_reactions)}")
    print(f"비가역 반응(Irreversible) 수: {len(model.reactions) - len(reversible_reactions)}")
    
    # 유전자-반응 연결 통계
    reactions_with_genes = [r for r in model.reactions if len(r.genes) > 0]
    print(f"유전자와 연결된 반응 수: {len(reactions_with_genes)}")
    
    print("="*60 + "\n")

def check_model_consistency(model):
    """
    모델의 일관성을 검사합니다.
    
    Parameters:
    -----------
    model : cobra.Model
        대사 모델
    """
    print("\n" + "="*60)
    print("모델 일관성 검사")
    print("="*60)
    
    # 기본 검증
    try:
        model.optimize()
        print("[OK] FBA 최적화 성공")
    except Exception as e:
        print(f"[ERROR] FBA 최적화 실패: {e}")
    
    # 반응과 대사물질 연결 확인
    reactions_with_no_metabolites = [r for r in model.reactions if len(r.metabolites) == 0]
    if reactions_with_no_metabolites:
        print(f"[WARNING] 대사물질이 없는 반응 수: {len(reactions_with_no_metabolites)}")
    else:
        print("[OK] 모든 반응에 대사물질이 연결되어 있습니다")
    
    print("="*60 + "\n")

def find_acs_reactions(model):
    """
    AcsA 관련 반응을 찾습니다.
    
    Parameters:
    -----------
    model : cobra.Model
        대사 모델
    
    Returns:
    --------
    acs_reactions : list
        AcsA 관련 반응 리스트
    """
    acs_reactions = []
    for reaction in model.reactions:
        # 반응 ID나 이름에 ACS가 포함된 경우
        if 'ACS' in reaction.id.upper() or 'ACS' in reaction.name.upper():
            acs_reactions.append(reaction)
        
        # Smlt4623 유전자와 연결된 반응
        for gene in reaction.genes:
            if 'Smlt4623' in gene.id or '4623' in gene.id:
                if reaction not in acs_reactions:
                    acs_reactions.append(reaction)
                break
    
    return acs_reactions

def print_acs_reactions(acs_reactions):
    """
    AcsA 관련 반응 정보를 출력합니다.
    
    Parameters:
    -----------
    acs_reactions : list
        AcsA 관련 반응 리스트
    """
    if acs_reactions:
        print("\n" + "="*60)
        print("AcsA 관련 반응 정보")
        print("="*60)
        for rxn in acs_reactions:
            print(f"\n반응 ID: {rxn.id}")
            print(f"반응 이름: {rxn.name}")
            print(f"반응식: {rxn.reaction}")
            print(f"유전자: {[gene.id for gene in rxn.genes]}")
            print(f"가역성: {rxn.reversibility}")
            print(f"하한: {rxn.lower_bound}, 상한: {rxn.upper_bound}")
        print("="*60 + "\n")
    else:
        print("\n[WARNING] AcsA 관련 반응을 찾을 수 없습니다.\n")

def main():
    """
    메인 함수
    """
    # 모델 로드
    model = load_model("BaseModel.xml")
    
    # 모델 요약 정보 출력
    print_model_summary(model)
    
    # AcsA 관련 반응 확인
    acs_reactions = find_acs_reactions(model)
    print_acs_reactions(acs_reactions)
    
    # 모델 일관성 검사
    check_model_consistency(model)
    
    # 모델 저장 (선택사항)
    # cobra.io.write_sbml_model(model, "GSM_model.xml")
    
    print("GSM 구축 완료!")
    return model

if __name__ == "__main__":
    model = main()

