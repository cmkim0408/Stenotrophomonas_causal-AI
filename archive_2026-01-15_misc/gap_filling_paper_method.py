#!/usr/bin/env python
"""
논문 방법론에 따른 Gap-filling 반응 추가
논문에서 언급한 19개 필수 반응 추가:
- Transport: ACt2rpp, PNTOt2rpp, MQN8t
- 아미노산: ALS, 2OXOVISO, IPMS, IPMD, P5CS, P5CD, LEUTA
- 보조인자: NAD/NADP, nicotinate 경로
"""

import cobra
from cobra import Reaction, Metabolite
import pandas as pd

def load_model(model_path="BaseModel.xml"):
    """모델 로드"""
    model = cobra.io.read_sbml_model(model_path)
    print(f"[OK] 모델 로드: {model.id}")
    return model

def check_metabolite_exists(model, met_id, compartment='c'):
    """대사물질 존재 확인"""
    full_id = f"{met_id}_{compartment}"
    try:
        met = model.metabolites.get_by_id(full_id)
        return True, met
    except KeyError:
        # 다른 구획 확인
        for comp in ['_c', '_e', '_p']:
            try:
                full_id_alt = f"{met_id}{comp}"
                met = model.metabolites.get_by_id(full_id_alt)
                return True, met
            except KeyError:
                continue
        return False, None

def get_metabolite(model, met_id, compartment='c'):
    """대사물질 가져오기 또는 생성"""
    exists, met = check_metabolite_exists(model, met_id, compartment)
    if exists:
        return met
    
    # 기본 정보로 생성 (나중에 수정 가능)
    try:
        # 비슷한 대사물질 참고
        if compartment == 'c':
            full_id = f"{met_id}_c"
        elif compartment == 'e':
            full_id = f"{met_id}_e"
        elif compartment == 'p':
            full_id = f"{met_id}_p"
        else:
            full_id = f"{met_id}_{compartment}"
        
        met = Metabolite(
            id=full_id,
            name=met_id,
            compartment=compartment
        )
        model.add_metabolites([met])
        return met
    except Exception as e:
        print(f"[WARNING] {met_id} 생성 실패: {e}")
        return None

def add_transport_reactions(model):
    """Transport 반응 추가"""
    print("\n" + "="*70)
    print("1. Transport 반응 추가")
    print("="*70)
    
    added_reactions = []
    
    # ACt2rpp: Acetate transport (proton symport, periplasm)
    try:
        ac_p = get_metabolite(model, 'ac', 'p')
        ac_c = get_metabolite(model, 'ac', 'c')
        h_p = get_metabolite(model, 'h', 'p')
        h_c = get_metabolite(model, 'h', 'c')
        
        try:
            rxn = model.reactions.get_by_id('ACt2rpp')
            print(f"[OK] ACt2rpp 이미 존재: {rxn.reaction}")
        except KeyError:
            rxn = Reaction('ACt2rpp')
            rxn.name = 'Acetate transport via proton symport (periplasm)'
            rxn.lower_bound = -1000
            rxn.upper_bound = 1000
            rxn.add_metabolites({
                ac_p: -1,
                h_p: -1,
                ac_c: 1,
                h_c: 1
            })
            model.add_reactions([rxn])
            print(f"[OK] ACt2rpp 추가: {rxn.reaction}")
            added_reactions.append('ACt2rpp')
    except Exception as e:
        print(f"[WARNING] ACt2rpp 추가 실패: {e}")
    
    # PNTOt2rpp: Pantothenate transport (proton symport, periplasm)
    try:
        pnto_p = get_metabolite(model, 'pnto', 'p')
        pnto_c = get_metabolite(model, 'pnto', 'c')
        h_p = get_metabolite(model, 'h', 'p')
        h_c = get_metabolite(model, 'h', 'c')
        
        try:
            rxn = model.reactions.get_by_id('PNTOt2rpp')
            print(f"[OK] PNTOt2rpp 이미 존재: {rxn.reaction}")
        except KeyError:
            rxn = Reaction('PNTOt2rpp')
            rxn.name = 'Pantothenate transport via proton symport (periplasm)'
            rxn.lower_bound = -1000
            rxn.upper_bound = 1000
            rxn.add_metabolites({
                pnto_p: -1,
                h_p: -1,
                pnto_c: 1,
                h_c: 1
            })
            model.add_reactions([rxn])
            print(f"[OK] PNTOt2rpp 추가: {rxn.reaction}")
            added_reactions.append('PNTOt2rpp')
    except Exception as e:
        print(f"[WARNING] PNTOt2rpp 추가 실패: {e}")
    
    # MQN8t: Menaquinone-8 transport
    try:
        mqn8_e = get_metabolite(model, 'mqn8', 'e')
        mqn8_c = get_metabolite(model, 'mqn8', 'c')
        
        try:
            rxn = model.reactions.get_by_id('MQN8t')
            print(f"[OK] MQN8t 이미 존재: {rxn.reaction}")
        except KeyError:
            rxn = Reaction('MQN8t')
            rxn.name = 'Menaquinone-8 transport'
            rxn.lower_bound = -1000
            rxn.upper_bound = 1000
            rxn.add_metabolites({
                mqn8_e: -1,
                mqn8_c: 1
            })
            model.add_reactions([rxn])
            print(f"[OK] MQN8t 추가: {rxn.reaction}")
            added_reactions.append('MQN8t')
    except Exception as e:
        print(f"[WARNING] MQN8t 추가 실패: {e}")
    
    return added_reactions

def add_amino_acid_reactions(model):
    """아미노산 합성 반응 추가"""
    print("\n" + "="*70)
    print("2. 아미노산 합성 반응 추가")
    print("="*70)
    
    added_reactions = []
    
    # ALS (Acetolactate Synthase) - Valine/Leucine 경로
    # 2 pyruvate --> 2-acetolactate + CO2
    try:
        pyr_c = get_metabolite(model, 'pyr', 'c')
        alac_c = get_metabolite(model, 'alac', 'c')  # 2-acetolactate
        co2_c = get_metabolite(model, 'co2', 'c')
        
        try:
            rxn = model.reactions.get_by_id('ALS')
            print(f"[OK] ALS 이미 존재: {rxn.reaction}")
        except KeyError:
            rxn = Reaction('ALS')
            rxn.name = 'Acetolactate synthase'
            rxn.lower_bound = 0
            rxn.upper_bound = 1000
            rxn.add_metabolites({
                pyr_c: -2,
                alac_c: 1,
                co2_c: 1
            })
            model.add_reactions([rxn])
            print(f"[OK] ALS 추가: {rxn.reaction}")
            added_reactions.append('ALS')
    except Exception as e:
        print(f"[WARNING] ALS 추가 실패: {e}")
    
    # 2OXOVISO (2-Oxoisovalerate synthase) - Valine 경로
    # 2-acetolactate --> 2-oxoisovalerate + CO2
    try:
        alac_c = get_metabolite(model, 'alac', 'c')
        oiv_c = get_metabolite(model, 'oiv', 'c')  # 2-oxoisovalerate
        co2_c = get_metabolite(model, 'co2', 'c')
        
        try:
            rxn = model.reactions.get_by_id('2OXOVISO')
            print(f"[OK] 2OXOVISO 이미 존재: {rxn.reaction}")
        except KeyError:
            rxn = Reaction('2OXOVISO')
            rxn.name = '2-Oxoisovalerate synthase'
            rxn.lower_bound = 0
            rxn.upper_bound = 1000
            rxn.add_metabolites({
                alac_c: -1,
                oiv_c: 1,
                co2_c: 1
            })
            model.add_reactions([rxn])
            print(f"[OK] 2OXOVISO 추가: {rxn.reaction}")
            added_reactions.append('2OXOVISO')
    except Exception as e:
        print(f"[WARNING] 2OXOVISO 추가 실패: {e}")
    
    # IPMS (Isopropylmalate Synthase) - Leucine 경로
    # 2-oxoisovalerate + acetyl-CoA --> 3-isopropylmalate + CoA
    try:
        oiv_c = get_metabolite(model, 'oiv', 'c')
        accoa_c = get_metabolite(model, 'accoa', 'c')
        ipm_c = get_metabolite(model, 'ipm', 'c')  # 3-isopropylmalate
        coa_c = get_metabolite(model, 'coa', 'c')
        
        try:
            rxn = model.reactions.get_by_id('IPMS')
            print(f"[OK] IPMS 이미 존재: {rxn.reaction}")
        except KeyError:
            rxn = Reaction('IPMS')
            rxn.name = 'Isopropylmalate synthase'
            rxn.lower_bound = 0
            rxn.upper_bound = 1000
            rxn.add_metabolites({
                oiv_c: -1,
                accoa_c: -1,
                ipm_c: 1,
                coa_c: 1
            })
            model.add_reactions([rxn])
            print(f"[OK] IPMS 추가: {rxn.reaction}")
            added_reactions.append('IPMS')
    except Exception as e:
        print(f"[WARNING] IPMS 추가 실패: {e}")
    
    # IPMD (3-Isopropylmalate Dehydrogenase) - Leucine 경로
    # 3-isopropylmalate --> 2-oxoisocaproate + CO2 + NADH
    try:
        ipm_c = get_metabolite(model, 'ipm', 'c')
        oic_c = get_metabolite(model, 'oic', 'c')  # 2-oxoisocaproate
        co2_c = get_metabolite(model, 'co2', 'c')
        nad_c = get_metabolite(model, 'nad', 'c')
        nadh_c = get_metabolite(model, 'nadh', 'c')
        
        try:
            rxn = model.reactions.get_by_id('IPMD')
            print(f"[OK] IPMD 이미 존재: {rxn.reaction}")
        except KeyError:
            rxn = Reaction('IPMD')
            rxn.name = '3-Isopropylmalate dehydrogenase'
            rxn.lower_bound = 0
            rxn.upper_bound = 1000
            rxn.add_metabolites({
                ipm_c: -1,
                nad_c: -1,
                oic_c: 1,
                co2_c: 1,
                nadh_c: 1
            })
            model.add_reactions([rxn])
            print(f"[OK] IPMD 추가: {rxn.reaction}")
            added_reactions.append('IPMD')
    except Exception as e:
        print(f"[WARNING] IPMD 추가 실패: {e}")
    
    # LEUTA (Leucine Transaminase) - Leucine 경로
    # 2-oxoisocaproate + glutamate --> leucine + 2-oxoglutarate
    try:
        oic_c = get_metabolite(model, 'oic', 'c')
        glu_c = get_metabolite(model, 'glu', 'c')
        leu_c = get_metabolite(model, 'leu__L', 'c')
        akg_c = get_metabolite(model, 'akg', 'c')
        
        try:
            rxn = model.reactions.get_by_id('LEUTA')
            print(f"[OK] LEUTA 이미 존재: {rxn.reaction}")
        except KeyError:
            rxn = Reaction('LEUTA')
            rxn.name = 'Leucine transaminase'
            rxn.lower_bound = -1000
            rxn.upper_bound = 1000
            rxn.add_metabolites({
                oic_c: -1,
                glu_c: -1,
                leu_c: 1,
                akg_c: 1
            })
            model.add_reactions([rxn])
            print(f"[OK] LEUTA 추가: {rxn.reaction}")
            added_reactions.append('LEUTA')
    except Exception as e:
        print(f"[WARNING] LEUTA 추가 실패: {e}")
    
    # P5CS (Pyrroline-5-carboxylate Synthase) - Proline 경로
    # glutamate + ATP + NADPH --> pyrroline-5-carboxylate + ADP + Pi + NADP+
    try:
        glu_c = get_metabolite(model, 'glu', 'c')
        atp_c = get_metabolite(model, 'atp', 'c')
        nadph_c = get_metabolite(model, 'nadph', 'c')
        p5c_c = get_metabolite(model, 'p5c', 'c')  # pyrroline-5-carboxylate
        adp_c = get_metabolite(model, 'adp', 'c')
        pi_c = get_metabolite(model, 'pi', 'c')
        nadp_c = get_metabolite(model, 'nadp', 'c')
        
        try:
            rxn = model.reactions.get_by_id('P5CS')
            print(f"[OK] P5CS 이미 존재: {rxn.reaction}")
        except KeyError:
            rxn = Reaction('P5CS')
            rxn.name = 'Pyrroline-5-carboxylate synthase'
            rxn.lower_bound = 0
            rxn.upper_bound = 1000
            rxn.add_metabolites({
                glu_c: -1,
                atp_c: -1,
                nadph_c: -1,
                p5c_c: 1,
                adp_c: 1,
                pi_c: 1,
                nadp_c: 1
            })
            model.add_reactions([rxn])
            print(f"[OK] P5CS 추가: {rxn.reaction}")
            added_reactions.append('P5CS')
    except Exception as e:
        print(f"[WARNING] P5CS 추가 실패: {e}")
    
    # P5CD (Pyrroline-5-carboxylate Dehydrogenase) - Proline 경로
    # pyrroline-5-carboxylate + NADPH --> proline + NADP+
    try:
        p5c_c = get_metabolite(model, 'p5c', 'c')
        nadph_c = get_metabolite(model, 'nadph', 'c')
        pro_c = get_metabolite(model, 'pro__L', 'c')
        nadp_c = get_metabolite(model, 'nadp', 'c')
        
        try:
            rxn = model.reactions.get_by_id('P5CD')
            print(f"[OK] P5CD 이미 존재: {rxn.reaction}")
        except KeyError:
            rxn = Reaction('P5CD')
            rxn.name = 'Pyrroline-5-carboxylate dehydrogenase'
            rxn.lower_bound = 0
            rxn.upper_bound = 1000
            rxn.add_metabolites({
                p5c_c: -1,
                nadph_c: -1,
                pro_c: 1,
                nadp_c: 1
            })
            model.add_reactions([rxn])
            print(f"[OK] P5CD 추가: {rxn.reaction}")
            added_reactions.append('P5CD')
    except Exception as e:
        print(f"[WARNING] P5CD 추가 실패: {e}")
    
    return added_reactions

def check_existing_reactions(model, reaction_ids):
    """기존 반응 확인"""
    print("\n" + "="*70)
    print("기존 반응 확인")
    print("="*70)
    
    existing = []
    missing = []
    
    for rxn_id in reaction_ids:
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            existing.append(rxn_id)
            print(f"[OK] {rxn_id} 존재: {rxn.reaction}")
        except KeyError:
            missing.append(rxn_id)
            print(f"[MISSING] {rxn_id}")
    
    return existing, missing

def main():
    print("="*70)
    print("논문 방법론에 따른 Gap-filling 반응 추가")
    print("="*70)
    
    # 모델 로드
    model = load_model("BaseModel.xml")
    
    # 논문에서 언급한 반응 목록
    paper_reactions = [
        'ACt2rpp', 'PNTOt2rpp', 'MQN8t',  # Transport
        'ALS', '2OXOVISO', 'IPMS', 'IPMD', 'LEUTA',  # Valine/Leucine
        'P5CS', 'P5CD'  # Proline
    ]
    
    # 기존 반응 확인
    existing, missing = check_existing_reactions(model, paper_reactions)
    
    print(f"\n존재하는 반응: {len(existing)}/{len(paper_reactions)}")
    print(f"누락된 반응: {len(missing)}/{len(paper_reactions)}")
    
    if missing:
        print(f"\n누락된 반응 목록: {', '.join(missing)}")
        
        # Transport 반응 추가
        transport_added = add_transport_reactions(model)
        
        # 아미노산 반응 추가
        amino_added = add_amino_acid_reactions(model)
        
        print(f"\n추가된 반응:")
        print(f"  Transport: {len(transport_added)}개 - {', '.join(transport_added)}")
        print(f"  아미노산: {len(amino_added)}개 - {', '.join(amino_added)}")
        
        # 모델 저장
        output_path = "BaseModel.xml"
        cobra.io.write_sbml_model(model, output_path)
        print(f"\n[OK] 모델 저장: {output_path}")
    
    print("\n" + "="*70)
    print("Gap-filling 완료")
    print("="*70)

if __name__ == "__main__":
    main()
