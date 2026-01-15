#!/usr/bin/env python
"""
아미노산 생합성 경로 확인
20종 아미노산 생합성 경로 확인
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

def check_amino_acid_synthesis(model):
    """아미노산 생합성 경로 확인"""
    print("\n" + "="*70)
    print("2단계: 아미노산 생합성 경로 확인")
    print("="*70)
    
    # 20종 아미노산
    amino_acids = {
        'ala__L_c': {'name': 'Alanine', 'pathway': 'Pyruvate'},
        'arg__L_c': {'name': 'Arginine', 'pathway': 'Glutamate'},
        'asn__L_c': {'name': 'Asparagine', 'pathway': 'Aspartate'},
        'asp__L_c': {'name': 'Aspartic acid', 'pathway': 'Oxaloacetate'},
        'cys__L_c': {'name': 'Cysteine', 'pathway': 'Serine'},
        'gln__L_c': {'name': 'Glutamine', 'pathway': 'Glutamate'},
        'glu__L_c': {'name': 'Glutamic acid', 'pathway': 'α-Ketoglutarate'},
        'gly_c': {'name': 'Glycine', 'pathway': 'Serine'},
        'his__L_c': {'name': 'Histidine', 'pathway': 'PRPP'},
        'ile__L_c': {'name': 'Isoleucine', 'pathway': 'Threonine'},
        'leu__L_c': {'name': 'Leucine', 'pathway': 'Pyruvate'},
        'lys__L_c': {'name': 'Lysine', 'pathway': 'Aspartate'},
        'met__L_c': {'name': 'Methionine', 'pathway': 'Aspartate'},
        'phe__L_c': {'name': 'Phenylalanine', 'pathway': 'Chorismate'},
        'pro__L_c': {'name': 'Proline', 'pathway': 'Glutamate'},
        'ser__L_c': {'name': 'Serine', 'pathway': '3-Phosphoglycerate'},
        'thr__L_c': {'name': 'Threonine', 'pathway': 'Aspartate'},
        'trp__L_c': {'name': 'Tryptophan', 'pathway': 'Chorismate'},
        'tyr__L_c': {'name': 'Tyrosine', 'pathway': 'Chorismate'},
        'val__L_c': {'name': 'Valine', 'pathway': 'Pyruvate'}
    }
    
    print("\n20종 아미노산 생합성 경로 확인:")
    print(f"{'아미노산':<20} {'이름':<20} {'경로':<20} {'생산 반응':<15} {'상태':<15}")
    print("-" * 95)
    
    amino_acid_status = []
    
    for aa_id, info in amino_acids.items():
        try:
            met = model.metabolites.get_by_id(aa_id)
            reactions = [r.id for r in met.reactions]
            
            # 생산 반응 찾기
            producing_rxns = []
            for rxn in met.reactions:
                if met in rxn.products:
                    producing_rxns.append(rxn)
            
            amino_acid_status.append({
                'Amino_Acid_ID': aa_id,
                'Name': info['name'],
                'Pathway': info['pathway'],
                'Exists': True,
                'Total_Reactions': len(reactions),
                'Producing_Reactions': len(producing_rxns),
                'Can_Produce': len(producing_rxns) > 0
            })
            
            status = "[OK]" if producing_rxns else "[WARNING]"
            print(f"{aa_id:<20} {info['name']:<20} {info['pathway']:<20} {len(producing_rxns):<15} {status:<15}")
            
            if producing_rxns:
                # 주요 생산 반응 출력
                for rxn in producing_rxns[:3]:
                    print(f"      - {rxn.id}: {rxn.name}")
                    print(f"        {rxn.reaction}")
            
        except KeyError:
            amino_acid_status.append({
                'Amino_Acid_ID': aa_id,
                'Name': info['name'],
                'Pathway': info['pathway'],
                'Exists': False,
                'Total_Reactions': 0,
                'Producing_Reactions': 0,
                'Can_Produce': False
            })
            print(f"{aa_id:<20} {info['name']:<20} {info['pathway']:<20} {'0':<15} {'[MISSING]':<15}")
    
    return amino_acid_status

def test_amino_acid_production(model, biomass_rxn):
    """아미노산 생산 가능 여부 테스트"""
    print("\n" + "="*70)
    print("아미노산 생산 가능 여부 테스트 (포도당 기반)")
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
    
    # 아미노산 생산 테스트
    amino_acids = ['ala__L_c', 'arg__L_c', 'asn__L_c', 'asp__L_c', 'cys__L_c',
                   'gln__L_c', 'glu__L_c', 'gly_c', 'his__L_c', 'ile__L_c',
                   'leu__L_c', 'lys__L_c', 'met__L_c', 'phe__L_c', 'pro__L_c',
                   'ser__L_c', 'thr__L_c', 'trp__L_c', 'tyr__L_c', 'val__L_c']
    
    print("\n아미노산 생산 가능 여부:")
    production_results = []
    
    for aa_id in amino_acids:
        try:
            met = model.metabolites.get_by_id(aa_id)
            
            # 임시 생산 반응
            test_rxn = cobra.Reaction(f'TEST_{aa_id}')
            test_rxn.add_metabolites({met: 1})
            test_rxn.lower_bound = 0
            test_rxn.upper_bound = 1000
            
            model.add_reactions([test_rxn])
            model.objective = test_rxn.id
            
            solution = model.optimize()
            
            model.remove_reactions([test_rxn])
            
            can_produce = solution.status == 'optimal' and solution.objective_value > 1e-6
            
            production_results.append({
                'Amino_Acid': aa_id,
                'Can_Produce': can_produce,
                'Status': solution.status,
                'Max_Flux': solution.objective_value if solution.status == 'optimal' else 0
            })
            
            status_icon = "[OK] 생산 가능" if can_produce else "[FAIL] 생산 불가"
            print(f"  {status_icon} {aa_id}")
            
        except KeyError:
            production_results.append({
                'Amino_Acid': aa_id,
                'Can_Produce': False,
                'Status': 'Metabolite missing',
                'Max_Flux': 0
            })
            print(f"  [MISSING] {aa_id}")
    
    return production_results

def check_key_amino_acid_pathways(model):
    """주요 아미노산 생합성 경로 확인"""
    print("\n" + "="*70)
    print("주요 아미노산 생합성 경로 상세 확인")
    print("="*70)
    
    # 중요한 아미노산 생합성 반응
    key_pathways = {
        'Serine pathway': ['SER', 'SERD', 'PGCD'],
        'Glycine pathway': ['GCY', 'SHMT'],
        'Branched-chain AA (Val/Leu/Ile)': ['ALS', '2OXOVISO', 'IPMS', 'IPMD', 'LEUTA'],
        'Proline pathway': ['P5CS', 'P5CD'],
        'Aromatic AA (Phe/Tyr/Trp)': ['CHORS', 'DAHPS', 'TRPS', 'TYRS', 'PHES'],
        'Aspartate family': ['ASPTA', 'ASPK', 'HOM', 'THRS', 'METS', 'LYSS']
    }
    
    for pathway_name, rxn_ids in key_pathways.items():
        print(f"\n[{pathway_name}]")
        found_count = 0
        for rxn_id in rxn_ids:
            try:
                rxn = model.reactions.get_by_id(rxn_id)
                found_count += 1
                print(f"  [OK] {rxn_id}: {rxn.name}")
            except KeyError:
                print(f"  [MISSING] {rxn_id}")
        print(f"  존재율: {found_count}/{len(rxn_ids)}")

def main():
    print("="*70)
    print("2단계: 아미노산 생합성 경로 확인")
    print("="*70)
    
    # 모델 로드
    model = load_model("BaseModel.xml")
    
    # Biomass 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응 없음")
        return
    
    print(f"[OK] Biomass: {biomass_rxn.id}")
    
    # 아미노산 생합성 경로 확인
    aa_status = check_amino_acid_synthesis(model)
    
    # 주요 경로 상세 확인
    check_key_amino_acid_pathways(model)
    
    # 생산 가능 여부 테스트
    production_results = test_amino_acid_production(model, biomass_rxn)
    
    # 결과 저장
    if aa_status:
        df_aa = pd.DataFrame(aa_status)
        df_aa.to_csv('amino_acid_synthesis_status.csv', index=False)
        print(f"\n[OK] 아미노산 상태 저장: amino_acid_synthesis_status.csv")
    
    if production_results:
        df_prod = pd.DataFrame(production_results)
        df_prod.to_csv('amino_acid_production_status.csv', index=False)
        print(f"[OK] 아미노산 생산 상태 저장: amino_acid_production_status.csv")
    
    # 요약
    can_produce_count = sum(1 for s in production_results if s['Can_Produce'])
    total_count = len(production_results)
    
    print("\n" + "="*70)
    print("아미노산 생합성 경로 분석 완료")
    print("="*70)
    print(f"\n생산 가능한 아미노산: {can_produce_count}/{total_count}")
    
    if can_produce_count < total_count:
        print("\n[문제 발견] 생산 불가능한 아미노산:")
        for status in production_results:
            if not status['Can_Produce']:
                print(f"  - {status['Amino_Acid']}: {status['Status']}")
    
    print("="*70)

if __name__ == "__main__":
    main()
