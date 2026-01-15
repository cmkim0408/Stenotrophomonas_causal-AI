#!/usr/bin/env python
"""
뉴클레오티드 생합성 경로 확인
- De novo purine synthesis
- De novo pyrimidine synthesis
- Nucleotide salvage pathway
- dNTP 합성 경로
"""

import cobra
import pandas as pd

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    print(f"[OK] 모델 로드: {model.id}")
    return model

def find_biomass_reaction(model):
    for name in ['Growth', 'BIOMASS']:
        try:
            return model.reactions.get_by_id(name)
        except KeyError:
            continue
    return None

def check_purine_synthesis(model):
    """Purine 생합성 경로 확인"""
    print("\n" + "="*70)
    print("1. Purine 생합성 경로 확인 (De novo)")
    print("="*70)
    
    # Purine 생합성 주요 반응
    purine_reactions = {
        'PRPPS': 'Phosphoribosyl pyrophosphate synthetase (R5P + ATP -> PRPP)',
        'GAR': 'Glycinamide ribotide synthetase',
        'FGAMS': 'Formylglycinamidine ribotide synthetase',
        'AIR': '5-Aminoimidazole ribotide synthetase',
        'ADSL': 'Adenylosuccinate lyase',
        'ADSL2r': 'Adenylosuccinate lyase (alternate)',
        'ADSS': 'Adenylosuccinate synthetase',
        'IMP': 'Inosine monophosphate 관련',
        'ADK': 'Adenosine kinase',
        'ADPRT': 'Adenine phosphoribosyltransferase',
        'GMP': 'Guanosine monophosphate 관련',
        'GMPR': 'GMP reductase',
        'GUAPRT': 'Guanine phosphoribosyltransferase'
    }
    
    # 실제 반응 ID 패턴 검색
    print("\nPurine 생합성 반응 검색:")
    found_reactions = []
    missing_reactions = []
    
    # 패턴 검색
    search_terms = ['PRPP', 'purine', 'adenine', 'guanine', 'IMP', 'AMP', 'GMP',
                    'adenylosuccinate', 'glycinamide', 'aminoimidazole', 'formylglycinamide']
    
    for term in search_terms:
        for rxn in model.reactions:
            if term.lower() in rxn.id.lower() or term.lower() in rxn.name.lower():
                if rxn not in found_reactions:
                    found_reactions.append(rxn)
    
    print(f"\n발견된 Purine 관련 반응: {len(found_reactions)}개")
    
    # IMP, AMP, GMP 관련 반응 찾기
    important_purines = ['imp', 'amp', 'gmp', 'adp', 'atp', 'gdp', 'gtp']
    
    print("\n주요 Purine 대사물질 및 반응:")
    for purine in important_purines:
        met_id = f'{purine}_c'
        try:
            met = model.metabolites.get_by_id(met_id)
            reactions = [r.id for r in met.reactions]
            print(f"\n  {purine.upper()}_c: {len(reactions)}개 반응")
            # Purine 생산 관련 반응 찾기
            producing_rxns = []
            for rxn in met.reactions:
                if met in rxn.products:
                    producing_rxns.append(rxn)
            
            if producing_rxns:
                print(f"    생산 반응: {len(producing_rxns)}개")
                for rxn in producing_rxns[:5]:
                    print(f"      - {rxn.id}: {rxn.name}")
                    print(f"        {rxn.reaction}")
            else:
                print(f"    [WARNING] 생산 반응 없음!")
        except KeyError:
            print(f"  {purine.upper()}_c: [MISSING]")
    
    return found_reactions

def check_pyrimidine_synthesis(model):
    """Pyrimidine 생합성 경로 확인"""
    print("\n" + "="*70)
    print("2. Pyrimidine 생합성 경로 확인 (De novo)")
    print("="*70)
    
    # Pyrimidine 생합성 주요 반응
    pyrimidine_reactions = {
        'CARPS': 'Carbamoyl phosphate synthetase',
        'ASPCT': 'Aspartate carbamoyltransferase',
        'DHORTS': 'Dihydroorotate synthase',
        'DHORD': 'Dihydroorotate dehydrogenase',
        'ORPT': 'Orotate phosphoribosyltransferase',
        'OMPDC': 'Orotidine-5-phosphate decarboxylase',
        'UMP': 'Uridine monophosphate 관련',
        'CMP': 'Cytidine monophosphate 관련',
        'CTPS': 'CTP synthetase',
        'NDPK': 'Nucleoside diphosphate kinase'
    }
    
    # 패턴 검색
    search_terms = ['pyrimidine', 'uracil', 'cytosine', 'thymine', 'orotate', 
                    'dihydroorotate', 'carbamoyl', 'UMP', 'CMP', 'UDP', 'CDP', 'UTP', 'CTP']
    
    found_reactions = []
    
    for term in search_terms:
        for rxn in model.reactions:
            if term.lower() in rxn.id.lower() or term.lower() in rxn.name.lower():
                if rxn not in found_reactions:
                    found_reactions.append(rxn)
    
    print(f"\n발견된 Pyrimidine 관련 반응: {len(found_reactions)}개")
    
    # UMP, CMP, UTP, CTP 관련 반응 찾기
    important_pyrimidines = ['ump', 'cmp', 'udp', 'cdp', 'utp', 'ctp', 'dtmp', 'dtdp', 'dttp', 'dcdp', 'dctp']
    
    print("\n주요 Pyrimidine 대사물질 및 반응:")
    for pyr in important_pyrimidines:
        met_id = f'{pyr}_c'
        try:
            met = model.metabolites.get_by_id(met_id)
            reactions = [r.id for r in met.reactions]
            print(f"\n  {pyr.upper()}_c: {len(reactions)}개 반응")
            
            # 생산 반응 찾기
            producing_rxns = []
            for rxn in met.reactions:
                if met in rxn.products:
                    producing_rxns.append(rxn)
            
            if producing_rxns:
                print(f"    생산 반응: {len(producing_rxns)}개")
                for rxn in producing_rxns[:5]:
                    print(f"      - {rxn.id}: {rxn.name}")
                    print(f"        {rxn.reaction}")
            else:
                print(f"    [WARNING] 생산 반응 없음!")
        except KeyError:
            print(f"  {pyr.upper()}_c: [MISSING]")
    
    return found_reactions

def check_dntp_synthesis(model):
    """dNTP 합성 경로 확인"""
    print("\n" + "="*70)
    print("3. dNTP 합성 경로 확인")
    print("="*70)
    
    # dNTP 관련 대사물질
    dntp_list = ['datp', 'dgtp', 'dctp', 'dttp', 'dadp', 'dgdp', 'dcdp', 'dtdp']
    
    print("\ndNTP 대사물질 확인:")
    dntp_status = []
    
    for dntp in dntp_list:
        met_id = f'{dntp}_c'
        try:
            met = model.metabolites.get_by_id(met_id)
            reactions = [r.id for r in met.reactions]
            
            # 생산 반응 찾기
            producing_rxns = []
            for rxn in met.reactions:
                if met in rxn.products:
                    producing_rxns.append(rxn)
            
            dntp_status.append({
                'dNTP': dntp.upper(),
                'Exists': True,
                'Reaction_Count': len(reactions),
                'Producing_Reactions': len(producing_rxns)
            })
            
            status = "[OK]" if producing_rxns else "[WARNING]"
            print(f"  {status} {dntp.upper()}_c: {len(reactions)}개 반응, {len(producing_rxns)}개 생산 반응")
            
            if producing_rxns:
                for rxn in producing_rxns[:3]:
                    print(f"      - {rxn.id}: {rxn.name}")
        except KeyError:
            dntp_status.append({
                'dNTP': dntp.upper(),
                'Exists': False,
                'Reaction_Count': 0,
                'Producing_Reactions': 0
            })
            print(f"  [MISSING] {dntp.upper()}_c")
    
    # Ribonucleotide reductase (RNR) 확인
    print("\nRibonucleotide Reductase (RNR) 확인:")
    rnr_patterns = ['RNR', 'ribonucleotide', 'reductase', 'NDPK']
    
    rnr_found = []
    for pattern in rnr_patterns:
        for rxn in model.reactions:
            if pattern.lower() in rxn.id.lower() or pattern.lower() in rxn.name.lower():
                if 'ribonucleotide' in rxn.name.lower() or 'NDPK' in rxn.id:
                    if rxn not in rnr_found:
                        rnr_found.append(rxn)
    
    if rnr_found:
        print(f"  발견된 RNR 관련 반응: {len(rnr_found)}개")
        for rxn in rnr_found[:5]:
            print(f"    - {rxn.id}: {rxn.name}")
            print(f"      {rxn.reaction}")
    else:
        print("  [WARNING] RNR 관련 반응 없음")
    
    return dntp_status, rnr_found

def test_nucleotide_production(model, biomass_rxn):
    """뉴클레오티드 생산 가능 여부 테스트"""
    print("\n" + "="*70)
    print("4. 뉴클레오티드 생산 가능 여부 테스트")
    print("="*70)
    
    # 포도당 설정
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    try:
        ex_glc = model.reactions.get_by_id('EX_glc__D_e')
        ex_glc.lower_bound = -100
        ex_glc.upper_bound = 1000
    except KeyError:
        pass
    
    # 필수 영양소
    essentials = ['EX_nh4_e', 'EX_h2o_e', 'EX_h_e', 'EX_pi_e', 'EX_so4_e',
                  'EX_k_e', 'EX_na1_e', 'EX_mg2_e', 'EX_ca2_e', 'EX_fe2_e',
                  'EX_mn2_e', 'EX_zn2_e', 'EX_co2_e', 'EX_o2_e']
    
    for ex_id in essentials:
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            if ex_id in ['EX_co2_e', 'EX_o2_e']:
                ex_rxn.lower_bound = -1000
                ex_rxn.upper_bound = 1000
            else:
                ex_rxn.lower_bound = -1000
        except KeyError:
            pass
    
    # 주요 뉴클레오티드 생산 테스트
    nucleotides = ['atp_c', 'gtp_c', 'utp_c', 'ctp_c', 'datp_c', 'dgtp_c', 'dctp_c', 'dttp_c']
    
    print("\n뉴클레오티드 생산 가능 여부 (포도당 기반):")
    production_status = []
    
    for nuc_id in nucleotides:
        try:
            met = model.metabolites.get_by_id(nuc_id)
            
            # 임시 생산 반응 생성
            test_rxn = cobra.Reaction(f'TEST_{nuc_id}')
            test_rxn.add_metabolites({met: 1})
            test_rxn.lower_bound = 0
            test_rxn.upper_bound = 1000
            
            model.add_reactions([test_rxn])
            model.objective = test_rxn.id
            
            solution = model.optimize()
            
            model.remove_reactions([test_rxn])
            
            can_produce = solution.status == 'optimal' and solution.objective_value > 1e-6
            
            production_status.append({
                'Nucleotide': nuc_id,
                'Can_Produce': can_produce,
                'Status': solution.status,
                'Max_Flux': solution.objective_value if solution.status == 'optimal' else 0
            })
            
            status_icon = "[OK] 생산 가능" if can_produce else "[FAIL] 생산 불가"
            print(f"  {status_icon} {nuc_id}")
            if solution.status == 'optimal':
                print(f"    최대 플럭스: {solution.objective_value:.6f}")
            
        except KeyError:
            production_status.append({
                'Nucleotide': nuc_id,
                'Can_Produce': False,
                'Status': 'Metabolite missing',
                'Max_Flux': 0
            })
            print(f"  [MISSING] {nuc_id}")
    
    return production_status

def main():
    print("="*70)
    print("1단계: 뉴클레오티드 생합성 경로 확인")
    print("="*70)
    
    # 모델 로드
    model = load_model("BaseModel.xml")
    
    # Biomass 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass: {biomass_rxn.id}")
    
    # 1. Purine 생합성 경로
    purine_rxns = check_purine_synthesis(model)
    
    # 2. Pyrimidine 생합성 경로
    pyrimidine_rxns = check_pyrimidine_synthesis(model)
    
    # 3. dNTP 합성 경로
    dntp_status, rnr_found = check_dntp_synthesis(model)
    
    # 4. 생산 가능 여부 테스트
    production_status = test_nucleotide_production(model, biomass_rxn)
    
    # 결과 저장
    if dntp_status:
        df_dntp = pd.DataFrame(dntp_status)
        df_dntp.to_csv('dntp_synthesis_status.csv', index=False)
        print(f"\n[OK] dNTP 상태 저장: dntp_synthesis_status.csv")
    
    if production_status:
        df_prod = pd.DataFrame(production_status)
        df_prod.to_csv('nucleotide_production_status.csv', index=False)
        print(f"[OK] 뉴클레오티드 생산 상태 저장: nucleotide_production_status.csv")
    
    # 요약
    can_produce_count = sum(1 for s in production_status if s['Can_Produce'])
    total_count = len(production_status)
    
    print("\n" + "="*70)
    print("뉴클레오티드 생합성 경로 분석 완료")
    print("="*70)
    print(f"\n생산 가능한 뉴클레오티드: {can_produce_count}/{total_count}")
    
    if can_produce_count < total_count:
        print("\n[문제 발견] 일부 뉴클레오티드 생산 불가능")
        for status in production_status:
            if not status['Can_Produce']:
                print(f"  - {status['Nucleotide']}: {status['Status']}")
    
    print("="*70)

if __name__ == "__main__":
    main()
