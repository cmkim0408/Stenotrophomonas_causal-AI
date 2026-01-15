#!/usr/bin/env python
"""
Acetate 기반 FBA (Flux Balance Analysis) 계산 및 분석
- Acetate minimal medium 설정
- Biomass 최대화 FBA 수행
- 주요 경로 플럭스 분석 (TCA, Glyoxylate shunt 등)
"""

import cobra
import pandas as pd
import numpy as np

def load_model(model_path="BaseModel.xml"):
    """모델 로드"""
    print(f"모델 로드 중: {model_path}")
    model = cobra.io.read_sbml_model(model_path)
    print(f"[OK] 모델 로드 완료: {model.id}\n")
    return model

def find_biomass_reaction(model):
    """Biomass 반응 찾기"""
    # Objective로 설정된 반응 확인
    if model.objective:
        try:
            obj_vars = list(model.objective.variables)
            for var in obj_vars:
                try:
                    rxn_id = var.name
                    rxn = model.reactions.get_by_id(rxn_id)
                    if 'growth' in rxn.id.lower() or 'biomass' in rxn.id.lower():
                        return rxn
                except (KeyError, AttributeError):
                    continue
        except AttributeError:
            pass
    
    # 이름으로 찾기
    for name in ['Growth', 'BIOMASS', 'BIOMASS_Ecoli_core', 'R_BIOMASS']:
        try:
            return model.reactions.get_by_id(name)
        except KeyError:
            continue
    
    # 패턴으로 찾기
    for rxn in model.reactions:
        if 'biomass' in rxn.id.lower() or 'growth' in rxn.id.lower():
            return rxn
    
    return None

def setup_acetate_medium(model):
    """Acetate minimal medium 설정"""
    print("="*70)
    print("Acetate Minimal Medium 설정")
    print("="*70)
    
    # 모든 exchange 차단
    for rxn in model.exchanges:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    
    # Acetate 허용 (탄소원)
    try:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        ex_ac.lower_bound = -1000  # 무제한 uptake
        ex_ac.upper_bound = 1000
        print(f"[OK] Acetate (탄소원): {ex_ac.id} (-1000 ~ 1000)")
    except KeyError:
        print("[WARNING] EX_ac_e를 찾을 수 없습니다. 다른 이름으로 검색 중...")
        for rxn in model.exchanges:
            if 'ac' in rxn.id.lower() and '_e' in rxn.id:
                rxn.lower_bound = -1000
                rxn.upper_bound = 1000
                print(f"[OK] Acetate (탄소원): {rxn.id} (-1000 ~ 1000)")
                break
    
    # 필수 무기염 및 영양소
    essential_exchanges = {
        'EX_nh4_e': 'Ammonium (질소원)',
        'EX_h2o_e': 'Water',
        'EX_h_e': 'Proton',
        'EX_pi_e': 'Phosphate',
        'EX_so4_e': 'Sulfate',
        'EX_k_e': 'Potassium',
        'EX_na1_e': 'Sodium',
        'EX_mg2_e': 'Magnesium',
        'EX_ca2_e': 'Calcium',
        'EX_fe2_e': 'Iron (II)',
        'EX_mn2_e': 'Manganese',
        'EX_zn2_e': 'Zinc',
        'EX_co2_e': 'CO2 (양방향)',
        'EX_o2_e': 'Oxygen (양방향)'
    }
    
    print("\n필수 영양소 설정:")
    for ex_id, name in essential_exchanges.items():
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            if ex_id in ['EX_co2_e', 'EX_o2_e']:
                ex_rxn.lower_bound = -1000  # 양방향 허용
                ex_rxn.upper_bound = 1000
            else:
                ex_rxn.lower_bound = -1000  # Uptake만 허용
                ex_rxn.upper_bound = 0
            print(f"  [OK] {name}: {ex_id}")
        except KeyError:
            pass
    
    return model

def run_fba(model, biomass_rxn):
    """FBA 최적화 수행"""
    print("\n" + "="*70)
    print("FBA 최적화 수행")
    print("="*70)
    
    model.objective = biomass_rxn.id
    print(f"Objective: {biomass_rxn.id}")
    
    try:
        solution = model.optimize()
        
        print(f"\n[최적화 결과]")
        print(f"  상태: {solution.status}")
        print(f"  Objective value: {solution.objective_value:.6f}")
        
        if solution.status == 'optimal':
            biomass_flux = solution.fluxes.get(biomass_rxn.id, 0)
            print(f"  Biomass flux: {biomass_flux:.6f} 1/h")
            
            if biomass_flux > 0:
                print(f"\n  [OK] 성장 가능!")
                if 0.001 <= biomass_flux <= 2.0:
                    print(f"    합리적인 성장률 범위 (0.001 ~ 2.0 1/h)")
                elif biomass_flux < 0.001:
                    print(f"    [WARNING] 매우 낮은 성장률 (< 0.001 1/h)")
            else:
                print(f"\n  [FAIL] 성장 불가 (Biomass flux = 0)")
        
        return solution
        
    except Exception as e:
        print(f"\n[ERROR] FBA 최적화 실패: {e}")
        import traceback
        traceback.print_exc()
        return None

def analyze_exchange_fluxes(model, solution):
    """Exchange 플럭스 분석"""
    print("\n" + "="*70)
    print("Exchange 플럭스 분석")
    print("="*70)
    
    exchange_data = []
    for rxn in model.exchanges:
        flux = solution.fluxes.get(rxn.id, 0)
        if abs(flux) > 1e-6:
            # 방향 결정
            if flux < 0:
                direction = "Uptake (섭취)"
                flux_abs = abs(flux)
            else:
                direction = "Secretion (분비)"
                flux_abs = flux
            
            exchange_data.append({
                '반응 ID': rxn.id,
                '이름': rxn.name,
                '플럭스': flux,
                '절대값': flux_abs,
                '방향': direction
            })
    
    if exchange_data:
        df = pd.DataFrame(exchange_data)
        df = df.sort_values('절대값', ascending=False)
        
        print(f"\n활성 Exchange 반응: {len(df)}개\n")
        print(df.to_string(index=False))
        
        # 주요 통계
        uptake_total = df[df['방향'] == 'Uptake (섭취)']['절대값'].sum()
        secretion_total = df[df['방향'] == 'Secretion (분비)']['절대값'].sum()
        
        print(f"\n총 섭취 플럭스: {uptake_total:.6f}")
        print(f"총 분비 플럭스: {secretion_total:.6f}")
    else:
        print("\n[WARNING] 활성 Exchange 반응이 없습니다")

def analyze_tca_glyoxylate_fluxes(model, solution):
    """TCA cycle 및 Glyoxylate shunt 플럭스 분석"""
    print("\n" + "="*70)
    print("TCA Cycle 및 Glyoxylate Shunt 플럭스 분석")
    print("="*70)
    
    # TCA cycle 반응
    tca_reactions = {
        'CS': 'Citrate Synthase',
        'ACONT': 'Aconitase',
        'ICDHx': 'Isocitrate Dehydrogenase (NAD+)',
        'ICDHyr': 'Isocitrate Dehydrogenase (NADP+)',
        'AKGDH': 'α-Ketoglutarate Dehydrogenase',
        'SUCOAS': 'Succinyl-CoA Synthetase',
        'SUCD': 'Succinate Dehydrogenase',
        'FUM': 'Fumarase',
        'MDH': 'Malate Dehydrogenase'
    }
    
    # Glyoxylate shunt 반응
    glyoxylate_reactions = {
        'ICL': 'Isocitrate Lyase',
        'MALS': 'Malate Synthase'
    }
    
    all_reactions = {**tca_reactions, **glyoxylate_reactions}
    
    flux_data = []
    for rxn_id, rxn_name in all_reactions.items():
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            flux = solution.fluxes.get(rxn_id, 0)
            
            pathway = 'TCA' if rxn_id in tca_reactions else 'Glyoxylate'
            
            flux_data.append({
                '경로': pathway,
                '반응 ID': rxn_id,
                '이름': rxn_name,
                '플럭스': flux,
                '절대값': abs(flux),
                '활성': 'Yes' if abs(flux) > 1e-6 else 'No'
            })
        except KeyError:
            flux_data.append({
                '경로': 'TCA' if rxn_id in tca_reactions else 'Glyoxylate',
                '반응 ID': rxn_id,
                '이름': rxn_name,
                '플럭스': 0,
                '절대값': 0,
                '활성': 'Missing'
            })
    
    df = pd.DataFrame(flux_data)
    
    print("\n[TCA Cycle 반응]")
    tca_df = df[df['경로'] == 'TCA'].copy()
    print(tca_df.to_string(index=False))
    
    print("\n[Glyoxylate Shunt 반응]")
    gly_df = df[df['경로'] == 'Glyoxylate'].copy()
    print(gly_df.to_string(index=False))
    
    # 활성 반응 요약
    active_tca = len(tca_df[tca_df['활성'] == 'Yes'])
    active_gly = len(gly_df[gly_df['활성'] == 'Yes'])
    
    print(f"\n활성 TCA cycle 반응: {active_tca}/{len(tca_reactions)}")
    print(f"활성 Glyoxylate shunt 반응: {active_gly}/{len(glyoxylate_reactions)}")
    
    # Glyoxylate shunt 활성 여부 확인
    if active_gly > 0:
        icl_flux = df[df['반응 ID'] == 'ICL']['플럭스'].values[0]
        mals_flux = df[df['반응 ID'] == 'MALS']['플럭스'].values[0]
        print(f"\n[Glyoxylate Shunt 활성]")
        print(f"  ICL 플럭스: {icl_flux:.6f}")
        print(f"  MALS 플럭스: {mals_flux:.6f}")
        if abs(icl_flux) > 1e-6 and abs(mals_flux) > 1e-6:
            print(f"  → Glyoxylate shunt가 활성화되어 탄소 보존 경로 사용 중")

def analyze_acetate_pathway_fluxes(model, solution):
    """Acetate 경로 플럭스 분석"""
    print("\n" + "="*70)
    print("Acetate → Biomass 경로 플럭스 분석")
    print("="*70)
    
    pathway_reactions = {
        'EX_ac_e': 'Acetate Exchange',
        'ACt': 'Acetate Transport',
        'R_ACS': 'Acetyl-CoA Synthetase',
        'ACS': 'Acetyl-CoA Synthetase (alternate)',
        'CS': 'Citrate Synthase',
        'ICL': 'Isocitrate Lyase',
        'MALS': 'Malate Synthase',
        'MDH': 'Malate Dehydrogenase'
    }
    
    print("\n주요 경로 반응:")
    pathway_data = []
    for rxn_id, rxn_name in pathway_reactions.items():
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            flux = solution.fluxes.get(rxn_id, 0)
            pathway_data.append({
                '단계': rxn_name,
                '반응 ID': rxn_id,
                '플럭스': flux,
                '절대값': abs(flux),
                '활성': 'Yes' if abs(flux) > 1e-6 else 'No'
            })
            if abs(flux) > 1e-6:
                print(f"  {rxn_id} ({rxn_name}): {flux:.6f}")
        except KeyError:
            pass
    
    if pathway_data:
        df = pd.DataFrame(pathway_data)
        print(f"\n활성 경로 반응: {len(df[df['활성'] == 'Yes'])}/{len(pathway_data)}")
    
    # Acetate uptake rate
    try:
        ex_ac = model.reactions.get_by_id('EX_ac_e')
        acetate_uptake = abs(solution.fluxes.get('EX_ac_e', 0))
        if acetate_uptake > 0:
            print(f"\nAcetate 섭취 속도: {acetate_uptake:.6f} mmol/gDW/h")
    except KeyError:
        pass

def analyze_energy_production(model, solution):
    """에너지 생성 분석"""
    print("\n" + "="*70)
    print("에너지 생성 분석")
    print("="*70)
    
    # ATP 생성/소비 반응
    atp_related = []
    for rxn in model.reactions:
        if abs(solution.fluxes.get(rxn.id, 0)) > 1e-6:
            try:
                metabolites = rxn.metabolites
                atp_c = model.metabolites.get_by_id('atp_c')
                adp_c = model.metabolites.get_by_id('adp_c')
                
                atp_production = 0
                atp_consumption = 0
                
                if atp_c in metabolites:
                    coeff = metabolites[atp_c]
                    flux = solution.fluxes.get(rxn.id, 0)
                    if coeff > 0:
                        atp_production = coeff * flux
                    else:
                        atp_consumption = abs(coeff * flux)
                
                if abs(atp_production) > 1e-6 or abs(atp_consumption) > 1e-6:
                    atp_related.append({
                        '반응 ID': rxn.id,
                        '이름': rxn.name,
                        '플럭스': flux,
                        'ATP 생성': atp_production if atp_production > 0 else 0,
                        'ATP 소비': atp_consumption if atp_consumption > 0 else 0
                    })
            except (KeyError, AttributeError):
                continue
    
    if atp_related:
        df = pd.DataFrame(atp_related)
        df = df.sort_values('ATP 생성', ascending=False)
        
        print("\n주요 ATP 생성 반응 (상위 10개):")
        top_producers = df[df['ATP 생성'] > 0].head(10)
        if len(top_producers) > 0:
            print(top_producers.to_string(index=False))
        
        total_atp_production = df['ATP 생성'].sum()
        total_atp_consumption = df['ATP 소비'].sum()
        
        print(f"\n총 ATP 생성: {total_atp_production:.6f} mmol/gDW/h")
        print(f"총 ATP 소비: {total_atp_consumption:.6f} mmol/gDW/h")
        print(f"순 ATP 생성: {total_atp_production - total_atp_consumption:.6f} mmol/gDW/h")

def main():
    """메인 함수"""
    print("="*70)
    print("Acetate 기반 FBA 계산 및 분석")
    print("="*70)
    
    # 1. 모델 로드
    model = load_model("BaseModel.xml")
    
    # 2. Biomass 반응 찾기
    biomass_rxn = find_biomass_reaction(model)
    if not biomass_rxn:
        print("[ERROR] Biomass 반응을 찾을 수 없습니다!")
        return
    
    print(f"[OK] Biomass 반응: {biomass_rxn.id}")
    
    # 3. Acetate medium 설정
    model = setup_acetate_medium(model)
    
    # 4. FBA 수행
    solution = run_fba(model, biomass_rxn)
    
    if not solution or solution.status != 'optimal':
        print("\n[ERROR] FBA 최적화 실패")
        return
    
    # 5. Exchange 플럭스 분석
    analyze_exchange_fluxes(model, solution)
    
    # 6. TCA/Glyoxylate 플럭스 분석
    analyze_tca_glyoxylate_fluxes(model, solution)
    
    # 7. Acetate 경로 플럭스 분석
    analyze_acetate_pathway_fluxes(model, solution)
    
    # 8. 에너지 생성 분석
    analyze_energy_production(model, solution)
    
    # 최종 요약
    print("\n" + "="*70)
    print("최종 요약")
    print("="*70)
    
    biomass_flux = solution.fluxes.get(biomass_rxn.id, 0)
    acetate_uptake = abs(solution.fluxes.get('EX_ac_e', 0))
    
    print(f"\nBiomass flux: {biomass_flux:.6f} 1/h")
    print(f"Acetate uptake: {acetate_uptake:.6f} mmol/gDW/h")
    
    if acetate_uptake > 0:
        yield_biomass = biomass_flux / acetate_uptake if acetate_uptake > 0 else 0
        print(f"Biomass yield: {yield_biomass:.6f} gDW/mmol acetate")
    
    print("\n[OK] FBA 분석 완료!")
    print("="*70)

if __name__ == "__main__":
    main()
