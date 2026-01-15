#!/usr/bin/env python
"""
레퍼런스 FBA 결과 파일의 실제 값 확인
- fba_flux_gradient_acid.csv에서 ATPM_0일 때 실제 플럭스 확인
- 레퍼런스 모델이 실제로 어떻게 작동했는지 확인
"""

import pandas as pd
from pathlib import Path

def check_reference_fba_results():
    """레퍼런스 FBA 결과 확인"""
    base_path = Path(__file__).parent.parent
    fba_file = base_path / "Stenotrophomonas" / "fba_flux_gradient_acid.csv"
    
    if not fba_file.exists():
        print(f"[오류] 파일 없음: {fba_file}")
        return None
    
    df = pd.read_csv(fba_file)
    
    print("="*80)
    print("레퍼런스 FBA 결과 파일 분석")
    print("="*80)
    
    print(f"\n[파일 정보]")
    print(f"  총 반응 수: {len(df)}")
    print(f"  컬럼: {list(df.columns)}")
    
    # ATPM_0 컬럼 확인
    if 'ATPM_0' in df.columns:
        print(f"\n[ATPM_0 컬럼 확인]")
        
        # 주요 반응 플럭스 확인
        key_reactions = ['ACS_ADP', 'ACS', 'CS', 'ICL', 'MALS', 'SUCDi', 'Growth', 'EX_ac_e', 'ATPS4rpp']
        
        print(f"\n[ATPM=0일 때 주요 반응 플럭스]")
        print(f"{'반응':<20} {'플럭스':<15}")
        print("-" * 35)
        
        for rxn_id in key_reactions:
            row = df[df.iloc[:,0] == rxn_id]
            if not row.empty:
                flux = row['ATPM_0'].iloc[0]
                if abs(flux) > 1e-6:
                    print(f"{rxn_id:<20} {flux:>15.6f}")
        
        # ATPM_5, ATPM_10도 확인
        print(f"\n[ATPM=5일 때 주요 반응 플럭스]")
        if 'ATPM_5' in df.columns:
            for rxn_id in key_reactions:
                row = df[df.iloc[:,0] == rxn_id]
                if not row.empty:
                    flux = row['ATPM_5'].iloc[0]
                    if abs(flux) > 1e-6:
                        print(f"{rxn_id:<20} {flux:>15.6f}")
        
        print(f"\n[ATPM=10일 때 주요 반응 플럭스]")
        if 'ATPM_10' in df.columns:
            for rxn_id in key_reactions:
                row = df[df.iloc[:,0] == rxn_id]
                if not row.empty:
                    flux = row['ATPM_10'].iloc[0]
                    if abs(flux) > 1e-6:
                        print(f"{rxn_id:<20} {flux:>15.6f}")
    
    # Exchange 플럭스 확인
    print(f"\n[Exchange 플럭스 확인 (ATPM=0)]")
    if 'ATPM_0' in df.columns:
        ex_reactions = ['EX_ac_e', 'EX_o2_e', 'EX_hco3_e', 'EX_nh4_e', 'EX_pi_e']
        for ex_id in ex_reactions:
            row = df[df.iloc[:,0] == ex_id]
            if not row.empty:
                flux = row['ATPM_0'].iloc[0]
                print(f"  {ex_id:20s}: {flux:>12.6f}")
    
    return df

def main():
    df = check_reference_fba_results()
    
    print("\n" + "="*80)
    print("결론")
    print("="*80)
    print(f"\n[레퍼런스 FBA 결과 분석]")
    print(f"  -> 실제 플럭스 값 확인 완료")
    print(f"  -> 이 값들을 기준으로 신규 모델과 비교 필요")

if __name__ == "__main__":
    main()
