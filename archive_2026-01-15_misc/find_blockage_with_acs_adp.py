#!/usr/bin/env python
"""
실제 막힌 지점 찾기 (ACS_ADP 포함)
"""

import cobra
from pathlib import Path
import csv

def load_model(model_path):
    model = cobra.io.read_sbml_model(model_path)
    return model

def add_acs_adp_if_missing(model, ref_model):
    """ACS_ADP 추가 (없으면)"""
    if 'ACS_ADP' in [r.id for r in model.reactions]:
        return False
    
    try:
        ref_acs_adp = ref_model.reactions.get_by_id('ACS_ADP')
        
        new_acs_adp = cobra.Reaction('ACS_ADP')
        new_acs_adp.name = ref_acs_adp.name
        new_acs_adp.lower_bound = ref_acs_adp.lower_bound
        new_acs_adp.upper_bound = ref_acs_adp.upper_bound
        
        metabolites_dict = {}
        for met, coeff in ref_acs_adp.metabolites.items():
            if met.id in [m.id for m in model.metabolites]:
                metabolites_dict[model.metabolites.get_by_id(met.id)] = coeff
            else:
                new_met = cobra.Metabolite(
                    met.id,
                    name=met.name,
                    compartment=met.compartment,
                    formula=met.formula if hasattr(met, 'formula') else None,
                    charge=met.charge if hasattr(met, 'charge') else None
                )
                model.add_metabolites([new_met])
                metabolites_dict[new_met] = coeff
        
        new_acs_adp.add_metabolites(metabolites_dict)
        model.add_reactions([new_acs_adp])
        return True
    except:
        return False

def apply_media_from_tsv(model, tsv_path):
    """TSV 파일에서 미디어 적용"""
    with open(tsv_path, "r", encoding="utf-8-sig", newline="") as f:
        raw_lines = f.readlines()
    
    lines = []
    for ln in raw_lines:
        ln_stripped = ln.strip()
        if ln_stripped and not ln_stripped.startswith("#"):
            lines.append(ln)
    
    rdr = csv.DictReader(lines, delimiter="\t")
    fieldnames = list(rdr.fieldnames)
    
    if len(fieldnames) >= 3 and fieldnames[1].replace(".", "").replace("-", "").isdigit():
        col_ex = fieldnames[0]
        col_lb = fieldnames[1]
        col_ub = fieldnames[2]
    else:
        col_ex = fieldnames[0]
        col_lb = fieldnames[1] if len(fieldnames) > 1 else fieldnames[0]
        col_ub = fieldnames[2] if len(fieldnames) > 2 else fieldnames[1]
    
    applied = 0
    for row in rdr:
        ex_id = row[col_ex].strip()
        if not ex_id or ex_id.lower() == "exchange_id":
            continue
        
        try:
            lb = float(row[col_lb])
            ub = float(row[col_ub])
            
            if ex_id in model.reactions:
                rxn = model.reactions.get_by_id(ex_id)
                rxn.lower_bound = lb
                rxn.upper_bound = ub
                applied += 1
        except (KeyError, ValueError, TypeError):
            continue
    
    return applied

def check_max_production(model, met_id):
    """메타볼라이트 최대 생산량 확인"""
    try:
        met = model.metabolites.get_by_id(met_id)
        demand = cobra.Reaction(f'DM_{met_id}')
        demand.lower_bound = 0
        demand.upper_bound = 1000
        demand.add_metabolites({met: -1.0})
        
        model_test = model.copy()
        model_test.add_reactions([demand])
        model_test.objective = demand.id
        
        solution = model_test.optimize()
        return solution.objective_value if solution.status == 'optimal' else 0.0
    except:
        return None

def check_acetate_pathway_step_by_step(model, media_path):
    """Acetate 경로를 단계별로 확인"""
    print("="*80)
    print("Acetate 경로 단계별 확인")
    print("="*80)
    
    # 미디어 설정
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    if media_path.exists():
        apply_media_from_tsv(model, media_path)
    
    # ATPM=0 설정
    atpm_rxn = model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    # 1단계: Acetate uptake 확인
    print(f"\n[1단계] Acetate uptake 확인")
    ex_ac = model.reactions.get_by_id('EX_ac_e')
    ex_ac.lower_bound = -1.0
    ex_ac.upper_bound = 1000
    
    model.objective = 'EX_ac_e'
    ex_ac.objective_coefficient = -1.0
    solution = model.optimize()
    
    ac_flux = solution.fluxes.get('EX_ac_e', 0.0)
    print(f"  EX_ac_e 플럭스: {ac_flux:.6f}")
    
    if abs(ac_flux) < 1e-6:
        print(f"  [문제] Acetate uptake 불가")
        return False
    else:
        print(f"  [OK] Acetate uptake 가능")
    
    # 2단계: ATP 생성 확인
    print(f"\n[2단계] ATP 생성 확인")
    atp_max = check_max_production(model, 'atp_c')
    print(f"  ATP 최대 생산량: {atp_max:.6f}")
    
    if atp_max < 1e-6:
        print(f"  [문제] ATP 생성 불가")
        return False
    else:
        print(f"  [OK] ATP 생성 가능")
    
    # 3단계: CoA 생성 확인
    print(f"\n[3단계] CoA 생성 확인")
    coa_max = check_max_production(model, 'coa_c')
    print(f"  CoA 최대 생산량: {coa_max:.6f}")
    
    if coa_max < 1e-6:
        print(f"  [문제] CoA 생성 불가")
    else:
        print(f"  [OK] CoA 생성 가능")
    
    # 4단계: CoA 부트스트랩 제공 시 ACS_ADP 작동 확인
    print(f"\n[4단계] CoA 부트스트랩 제공 시 ACS_ADP 작동 확인")
    
    # CoA exchange 추가
    try:
        coa_c = model.metabolites.get_by_id('coa_c')
        ex_coa_id = 'EX_coa_c'
        
        if ex_coa_id not in [r.id for r in model.exchanges]:
            ex_coa = cobra.Reaction(ex_coa_id)
            ex_coa.name = 'CoA exchange'
            ex_coa.lower_bound = -0.001
            ex_coa.upper_bound = 1000
            
            coa_e_id = 'coa_e'
            try:
                coa_e = model.metabolites.get_by_id(coa_e_id)
            except KeyError:
                coa_e = cobra.Metabolite(coa_e_id, name='CoA', compartment='e')
                model.add_metabolites([coa_e])
            
            ex_coa.add_metabolites({coa_e: -1.0})
            model.add_reactions([ex_coa])
        
        # Transport 추가
        trans_id = 'T_coa_e_to_coa_c'
        if trans_id not in [r.id for r in model.reactions]:
            trans = cobra.Reaction(trans_id)
            trans.name = 'CoA transport'
            trans.lower_bound = -1000
            trans.upper_bound = 1000
            trans.add_metabolites({coa_e: -1.0, coa_c: 1.0})
            model.add_reactions([trans])
    except KeyError:
        print(f"  [오류] coa_c 메타볼라이트 없음")
        return False
    
    model.objective = 'Growth'
    solution = model.optimize()
    
    acs_adp_flux = solution.fluxes.get('ACS_ADP', 0.0)
    acs_flux = solution.fluxes.get('ACS', 0.0)
    
    print(f"  ACS_ADP 플럭스: {acs_adp_flux:.6f}")
    print(f"  ACS 플럭스: {acs_flux:.6f}")
    
    if abs(acs_adp_flux) < 1e-6 and abs(acs_flux) < 1e-6:
        print(f"  [문제] CoA 제공해도 ACS/ACS_ADP 작동 안 함")
        
        # ACS_ADP 반응식 확인
        try:
            acs_adp = model.reactions.get_by_id('ACS_ADP')
            print(f"\n  ACS_ADP 반응식: {acs_adp.reaction}")
            
            # 필요한 메타볼라이트 확인
            print(f"  필요한 메타볼라이트:")
            for met, coeff in acs_adp.metabolites.items():
                if coeff < 0:  # 소모
                    print(f"    {met.id} (계수: {coeff})")
                    met_max = check_max_production(model, met.id)
                    if met_max is not None:
                        print(f"      최대 생산량: {met_max:.6f}")
                        if met_max < 1e-6:
                            print(f"      [문제] {met.id} 생성 불가!")
        except:
            pass
        
        return False
    else:
        print(f"  [OK] CoA 제공 시 ACS/ACS_ADP 작동")
        return True
    
    # 5단계: Acetyl-CoA -> Citrate 확인
    print(f"\n[5단계] Acetyl-CoA -> Citrate 확인")
    
    # Acetyl-CoA와 OAA 부트스트랩
    try:
        accoa_c = model.metabolites.get_by_id('accoa_c')
        oaa_c = model.metabolites.get_by_id('oaa_c')
        
        for met, met_id in [(accoa_c, 'accoa_c'), (oaa_c, 'oaa_c')]:
            ex_id = f'EX_{met_id}'
            if ex_id not in [r.id for r in model.exchanges]:
                ex_rxn = cobra.Reaction(ex_id)
                ex_rxn.name = f'{met_id} exchange'
                ex_rxn.lower_bound = -0.001
                ex_rxn.upper_bound = 1000
                
                met_e_id = met_id.replace('_c', '_e')
                try:
                    met_e = model.metabolites.get_by_id(met_e_id)
                except KeyError:
                    met_e = cobra.Metabolite(met_e_id, name=met.name, compartment='e')
                    model.add_metabolites([met_e])
                
                ex_rxn.add_metabolites({met_e: -1.0})
                model.add_reactions([ex_rxn])
                
                # Transport 추가
                trans_id = f'T_{met_e_id}_to_{met_id}'
                if trans_id not in [r.id for r in model.reactions]:
                    trans = cobra.Reaction(trans_id)
                    trans.name = f'{met_id} transport'
                    trans.lower_bound = -1000
                    trans.upper_bound = 1000
                    trans.add_metabolites({met_e: -1.0, met: 1.0})
                    model.add_reactions([trans])
    except KeyError as e:
        print(f"  [오류] 메타볼라이트 없음: {e}")
        return False
    
    model.objective = 'Growth'
    solution = model.optimize()
    
    cs_flux = solution.fluxes.get('CS', 0.0)
    print(f"  CS 플럭스: {cs_flux:.6f}")
    
    if abs(cs_flux) < 1e-6:
        print(f"  [문제] Acetyl-CoA와 OAA 제공해도 CS 작동 안 함")
        return False
    else:
        print(f"  [OK] Acetyl-CoA와 OAA 제공 시 CS 작동")
        return True

def main():
    base_path = Path(__file__).parent.parent
    model_path = base_path / "Stenotrophomonas-causal AI" / "BaseModel.xml"
    ref_model_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5_clean" / "model_YE0p5.xml"
    media_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5" / "Acetate_YE0p5__nocmnt__normalized.tsv"
    
    model = load_model(str(model_path))
    
    if not ref_model_path.exists():
        print(f"[오류] 레퍼런스 모델 파일 없음: {ref_model_path}")
        return
    
    ref_model = cobra.io.read_sbml_model(str(ref_model_path))
    
    print("="*80)
    print("실제 막힌 지점 찾기")
    print("="*80)
    
    # ACS_ADP 추가
    added = add_acs_adp_if_missing(model, ref_model)
    if added:
        print(f"\n[ACS_ADP 추가 완료]")
    
    # 단계별 확인
    result = check_acetate_pathway_step_by_step(model, media_path)
    
    print("\n" + "="*80)
    print("결론")
    print("="*80)
    
    if result:
        print(f"\n[성공] 경로가 작동함")
    else:
        print(f"\n[실패] 경로가 막혀있음")
        print(f"  -> 위의 단계별 확인 결과를 참고하여 막힌 지점 확인")

if __name__ == "__main__":
    main()
