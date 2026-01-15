#!/usr/bin/env python
"""
종합 진단 보고서 - 1, 2, 3단계 결과 통합
"""

import pandas as pd

def create_final_report():
    """최종 종합 보고서 생성"""
    
    print("="*70)
    print("모델 종합 진단 보고서")
    print("1, 2, 3단계 분석 결과 통합")
    print("="*70)
    
    print("\n[진단 완료 항목]")
    print("\n[1단계: Biomass 반응 검증]")
    print("  - Biomass 반응: Growth (존재 확인)")
    print("  - 총 구성 요소: 57개")
    print("  - 주요 구성 요소:")
    print("    * ATP: -54.12 (대량 필요)")
    print("    * 아미노산: 20종")
    print("    * 기타 뉴클레오티드: GTP, UTP, CTP, dNTP")
    print("    * 보조인자: CoA, NAD+, NADP+, FAD 등")
    
    print("\n[2단계: 기본 대사 경로 연결성]")
    print("  - Glycolysis: 10/10 반응 존재 [OK]")
    print("  - TCA Cycle: 8/8 반응 존재 [OK]")
    print("  - Glyoxylate Shunt: 2/2 반응 존재 [OK]")
    print("  - 주요 대사물질: 모두 존재 및 연결됨 [OK]")
    
    print("\n[3단계: 모델 구조 검토]")
    print("  - 무제한 영양소: 성장 가능 (Biomass flux: 63.37 1/h) [OK]")
    print("  - 포도당만 허용: infeasible [FAIL]")
    print("  - Acetate만 허용: optimal이지만 Biomass flux = 0 [FAIL]")
    
    print("\n" + "="*70)
    print("[1단계 상세 분석: 뉴클레오티드 생합성 경로]")
    print("="*70)
    
    print("\n1. Purine 생합성 경로:")
    print("  - 발견된 Purine 관련 반응: 68개")
    print("  - 주요 대사물질 상태:")
    print("    * IMP_c: 6개 반응, 2개 생산 반응 [OK]")
    print("    * AMP_c: 112개 반응, 106개 생산 반응 [OK]")
    print("    * GMP_c: 10개 반응, 8개 생산 반응 [OK]")
    print("    * ADP_c: 270개 반응, 263개 생산 반응 [OK]")
    print("    * ATP_c: 367개 반응, 3개 생산 반응 [WARNING]")
    print("    * GDP_c: 18개 반응, 10개 생산 반응 [OK]")
    print("    * GTP_c: 23개 반응, 4개 생산 반응 [OK]")
    
    print("\n2. Pyrimidine 생합성 경로:")
    print("  - 발견된 Pyrimidine 관련 반응: 109개")
    print("  - 주요 대사물질 상태:")
    print("    * UMP_c: 7개 반응, 5개 생산 반응 [OK]")
    print("    * CMP_c: 24개 반응, 21개 생산 반응 [OK]")
    print("    * UDP_c: 11개 반응, 7개 생산 반응 [OK]")
    print("    * CDP_c: 9개 반응, 5개 생산 반응 [OK]")
    
    print("\n3. dNTP 합성 경로:")
    print("  - dNTP 관련 반응은 존재하나 생산 경로 확인 필요")
    print("  - Ribonucleotide reductase (RNR) 경로 확인 필요")
    
    print("\n4. 뉴클레오티드 생산 가능 여부 (포도당 기반):")
    print("  - ATP, GTP, UTP, CTP 생산 경로 문제")
    print("  - 포도당만으로는 생산 불가능")
    
    print("\n" + "="*70)
    print("[2단계 상세 분석: 아미노산 생합성 경로]")
    print("="*70)
    
    print("\n1. 20종 아미노산 생합성 경로:")
    print("  - 대부분의 아미노산에 생산 반응 존재 [OK]")
    print("  - 그러나 포도당 기반 생산 테스트 결과 infeasible")
    
    print("\n2. 주요 아미노산 생합성 경로 상태:")
    print("  - Serine pathway: 1/3 반응 존재 [PARTIAL]")
    print("    * [MISSING] SER, SERD")
    print("    * [OK] PGCD")
    print("  - Glycine pathway: 0/2 반응 존재 [MISSING]")
    print("    * [MISSING] GCY, SHMT")
    print("  - Branched-chain AA (Val/Leu/Ile): 5/5 반응 존재 [OK]")
    print("    * [OK] ALS, 2OXOVISO, IPMS, IPMD, LEUTA")
    print("  - Proline pathway: 2/2 반응 존재 [OK]")
    print("    * [OK] P5CS, P5CD")
    print("  - Aromatic AA (Phe/Tyr/Trp): 1/5 반응 존재 [PARTIAL]")
    print("    * [OK] CHORS")
    print("    * [MISSING] DAHPS, TRPS, TYRS, PHES")
    print("  - Aspartate family: 4/6 반응 존재 [PARTIAL]")
    print("    * [OK] ASPTA, ASPK, THRS, METS")
    print("    * [MISSING] HOM, LYSS")
    
    print("\n3. 아미노산 생산 가능 여부 (포도당 기반):")
    print("  - 대부분의 아미노산 생산 반응은 존재하나")
    print("  - 포도당만으로는 생산 불가능 (infeasible)")
    print("  - 추가 영양소 또는 경로 필요")
    
    print("\n" + "="*70)
    print("[3단계 상세 분석: 보조인자 생산 경로]")
    print("="*70)
    
    print("\n1. CoA 생산 경로:")
    print("  - 경로 상태: 2/5 반응 존재 [PARTIAL]")
    print("    * [MISSING] PNTO: Pantothenate kinase")
    print("    * [MISSING] PPCS: 4-phosphopantothenoylcysteine synthetase")
    print("    * [OK] PPCDC: Phosphopantothenoylcysteine decarboxylase")
    print("    * [MISSING] PPAT: Phosphopantetheine adenylyltransferase")
    print("    * [OK] DPCOAK: Dephospho-CoA kinase")
    print("  - 중간 대사물질 상태:")
    print("    * [MISSING] pnto_c (Pantothenate)")
    print("    * [OK] 4ppan_c (2개 반응)")
    print("    * [MISSING] pppcs_c, pppan_c, dcoa_c")
    print("    * [OK] coa_c (191개 반응, 74개 생산 반응)")
    print("  - 생산 가능 여부: [FAIL] infeasible")
    print("  - 문제: Pantothenate -> CoA 경로 불완전")
    
    print("\n2. NAD+/NADP+ 생산 경로:")
    print("  - 경로 상태: 3/5 반응 존재 [PARTIAL]")
    print("    * [OK] NADS1: NAD synthase (nh3)")
    print("    * [OK] NADS2: Nicotinate-mononucleotide adenylyltransferase")
    print("    * [OK] NADK: NAD kinase")
    print("    * [MISSING] NADPPPS, NADDPPPS")
    print("  - 대사물질 상태:")
    print("    * [OK] nad_c, nadh_c, nadp_c, nadph_c 모두 존재")
    print("  - 생산 가능 여부: [FAIL] infeasible")
    
    print("\n3. FAD 생산 경로:")
    print("  - 경로 상태: 2/2 반응 존재 [OK]")
    print("    * [OK] RBFK: Riboflavin kinase")
    print("    * [OK] FMNAT: FMN adenylyltransferase")
    print("  - 대사물질 상태:")
    print("    * [OK] ribflv_c, fmn_c, fad_c, fadh2_c 모두 존재")
    print("  - 생산 가능 여부: [FAIL] infeasible")
    
    print("\n4. Folate 관련 경로:")
    print("  - 대사물질 상태: [OK]")
    print("    * [OK] thf_c, mlthf_c, 10fthf_c, dhf_c 모두 존재")
    print("  - 생산 가능 여부: [FAIL] infeasible")
    
    print("\n" + "="*70)
    print("[핵심 문제 요약]")
    print("="*70)
    
    print("\n[문제 1: 뉴클레오티드 생산 경로]")
    print("  - ATP 생산 반응: 3개만 존재 (367개 반응 중)")
    print("  - 포도당만으로는 ATP 생산 불가능")
    print("  - GTP, UTP, CTP도 생산 불가능")
    print("  - 뉴클레오티드 De novo synthesis 경로 불완전 가능성")
    
    print("\n[문제 2: 아미노산 생합성 경로]")
    print("  - 일부 아미노산 생합성 경로 불완전:")
    print("    * Serine pathway: SER, SERD 누락")
    print("    * Glycine pathway: GCY, SHMT 누락")
    print("    * Aromatic AA: DAHPS, TRPS, TYRS, PHES 누락")
    print("    * Aspartate family: HOM, LYSS 누락")
    print("  - 포도당만으로는 아미노산 생산 불가능")
    
    print("\n[문제 3: 보조인자 생산 경로]")
    print("  - CoA 생산 경로 불완전:")
    print("    * PNTO (Pantothenate kinase) 누락")
    print("    * PPCS (4-phosphopantothenoylcysteine synthetase) 누락")
    print("    * PPAT (Phosphopantetheine adenylyltransferase) 누락")
    print("    * Pantothenate -> CoA 경로 차단")
    print("  - 모든 보조인자 (CoA, NAD+, NADP+, FAD, Folate) 생산 불가능")
    print("  - 포도당만으로는 생산 불가능 (infeasible)")
    
    print("\n" + "="*70)
    print("[종합 결론]")
    print("="*70)
    
    print("\n[모델 상태]")
    print("  1. 모델 구조: 정상 (무제한 영양소로 성장 가능)")
    print("  2. 탄소원 처리 경로: 정상 (Glycolysis, TCA, Glyoxylate 모두 존재)")
    print("  3. 뉴클레오티드 생합성: 부분적 (경로는 있으나 생산 불가)")
    print("  4. 아미노산 생합성: 부분적 (일부 경로 누락)")
    print("  5. 보조인자 생산: 불완전 (특히 CoA 경로)")
    
    print("\n[핵심 문제]")
    print("  - 포도당만으로는 Biomass 구성 요소 생산 불가능")
    print("  - 뉴클레오티드, 아미노산, 보조인자 모두 생산 불가능")
    print("  - 무제한 영양소로는 성장 가능 -> 구조는 정상")
    print("  - -> 생합성 경로가 불완전하거나 부트스트랩 문제")
    
    print("\n[권장 해결 방안]")
    print("  1. 누락된 반응 추가:")
    print("     - CoA 생산 경로: PNTO, PPCS, PPAT")
    print("     - 아미노산 생합성: SER, SERD, GCY, SHMT 등")
    print("     - 뉴클레오티드 생합성 경로 보완")
    print("  2. 부트스트랩 영양소 추가:")
    print("     - 소량의 CoA 또는 Pantothenate")
    print("     - 소량의 뉴클레오티드 (ATP, GTP, UTP, CTP)")
    print("     - 소량의 NAD+/NADP+")
    print("  3. Gap-filling 수행:")
    print("     - 생산 불가능한 대사물질에 대한 경로 추가")
    print("     - 부족한 반응 식별 및 추가")
    
    print("\n" + "="*70)
    print("진단 완료")
    print("="*70)

def main():
    create_final_report()
    
    # 결과 요약 CSV 저장
    summary = {
        'Category': [
            '모델 구조',
            'Glycolysis',
            'TCA Cycle',
            'Glyoxylate Shunt',
            '무제한 영양소 성장',
            'Purine 생합성',
            'Pyrimidine 생합성',
            'ATP 생산',
            '아미노산 생합성 (전체)',
            'Serine pathway',
            'Glycine pathway',
            'Branched-chain AA',
            'Proline pathway',
            'Aromatic AA',
            'Aspartate family',
            'CoA 생산 경로',
            'NAD+/NADP+ 경로',
            'FAD 경로',
            '보조인자 생산 (전체)'
        ],
        'Status': [
            '정상',
            '10/10 반응 존재',
            '8/8 반응 존재',
            '2/2 반응 존재',
            '가능 (63.37 1/h)',
            '부분적 (경로 존재)',
            '부분적 (경로 존재)',
            '생산 불가',
            '부분적 (일부 경로 누락)',
            '1/3 반응 존재',
            '0/2 반응 존재',
            '5/5 반응 존재',
            '2/2 반응 존재',
            '1/5 반응 존재',
            '4/6 반응 존재',
            '2/5 반응 존재',
            '3/5 반응 존재',
            '2/2 반응 존재',
            '생산 불가'
        ],
        'Can_Produce_Glucose': [
            '가능',
            'N/A',
            'N/A',
            'N/A',
            '가능',
            '불가능',
            '불가능',
            '불가능',
            '불가능',
            '불가능',
            '불가능',
            '불가능',
            '불가능',
            '불가능',
            '불가능',
            '불가능',
            '불가능',
            '불가능',
            '불가능'
        ]
    }
    
    df_summary = pd.DataFrame(summary)
    df_summary.to_csv('comprehensive_diagnosis_summary.csv', index=False)
    print(f"\n[OK] 종합 진단 요약 저장: comprehensive_diagnosis_summary.csv")

if __name__ == "__main__":
    main()
