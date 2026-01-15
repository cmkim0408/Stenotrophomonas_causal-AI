#!/usr/bin/env python
"""
레퍼런스 모델의 실제 FBA 설정 확인
- fba_flux_gradient_acid.csv가 어떻게 생성되었는지 확인
"""

import pandas as pd
from pathlib import Path

def check_reference_fba_results():
    """레퍼런스 FBA 결과 확인"""
    print("="*80)
    print("레퍼런스 모델 FBA 결과 확인")
    print("="*80)
    
    base_path = Path(__file__).parent.parent
    fba_file = base_path / "Stenotrophomonas" / "fba_flux_gradient_acid.csv"
    
    if not fba_file.exists():
        print(f"[오류] 파일 없음: {fba_file}")
        return None
    
    df = pd.read_csv(fba_file)
    
    print(f"\n[파일 정보]")
    print(f"  총 반응 수: {len(df)}")
    print(f"  컬럼: {list(df.columns)}")
    
    # ACS_ADP 플럭스 확인
    acs_adp_row = df[df.iloc[:,0] == 'ACS_ADP']
    if not acs_adp_row.empty:
        print(f"\n[ACS_ADP 플럭스]")
        for col in df.columns[1:]:
            flux = acs_adp_row.iloc[0][col]
            print(f"  {col}: {flux:.6f}")
    
    # Growth 플럭스 확인
    growth_row = df[df.iloc[:,0] == 'Growth']
    if not growth_row.empty:
        print(f"\n[Growth 플럭스]")
        for col in df.columns[1:]:
            flux = growth_row.iloc[0][col]
            print(f"  {col}: {flux:.6f}")
    
    # EX_pnto__R_e 플럭스 확인
    ex_pnto_row = df[df.iloc[:,0] == 'EX_pnto__R_e']
    if not ex_pnto_row.empty:
        print(f"\n[EX_pnto__R_e 플럭스]")
        for col in df.columns[1:]:
            flux = ex_pnto_row.iloc[0][col]
            print(f"  {col}: {flux:.6f}")
    
    # ICL, MALS 플럭스 확인
    icl_row = df[df.iloc[:,0] == 'ICL']
    mals_row = df[df.iloc[:,0] == 'MALS']
    
    if not icl_row.empty:
        print(f"\n[ICL 플럭스]")
        for col in df.columns[1:]:
            flux = icl_row.iloc[0][col]
            print(f"  {col}: {flux:.6f}")
    
    if not mals_row.empty:
        print(f"\n[MALS 플럭스]")
        for col in df.columns[1:]:
            flux = mals_row.iloc[0][col]
            print(f"  {col}: {flux:.6f}")
    
    return df

def find_fba_script():
    """FBA 스크립트 찾기"""
    print("\n" + "="*80)
    print("FBA 생성 스크립트 찾기")
    print("="*80)
    
    base_path = Path(__file__).parent.parent
    sten_path = base_path / "Stenotrophomonas"
    
    # fba_flux_gradient 관련 스크립트 찾기
    fba_scripts = list(sten_path.glob("*flux_gradient*.py"))
    fba_scripts += list(sten_path.glob("*fba*gradient*.py"))
    
    print(f"\n[발견된 스크립트] {len(fba_scripts)}개")
    for script in fba_scripts:
        print(f"  {script.name}")
    
    return fba_scripts

def main():
    print("="*80)
    print("레퍼런스 모델 FBA 설정 확인")
    print("="*80)
    
    # FBA 결과 확인
    df = check_reference_fba_results()
    
    # 스크립트 찾기
    scripts = find_fba_script()
    
    print("\n" + "="*80)
    print("결론")
    print("="*80)
    
    if df is not None:
        print(f"\n[확인] 레퍼런스 FBA 결과 파일 분석 완료")
        print(f"  -> 실제 사용된 설정 확인 필요")

if __name__ == "__main__":
    main()
