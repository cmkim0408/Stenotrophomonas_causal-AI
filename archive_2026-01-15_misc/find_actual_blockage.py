#!/usr/bin/env python
"""
실제 막힌 지점 찾기
- CoA를 제공해도 작동하지 않으므로 다른 곳이 막혀있음
- 단계별로 확인
"""

import cobra
from pathlib import Path
import csv

def load_model(model_path):
    model = cobra.io.read_sbml_model(model_path)
    return model

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

def test_with_bootstrap(model, media_path, bootstrap_mets=None):
    """부트스트랩 제공 후 테스트"""
    # 미디어 설정
    for rxn in model.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    if media_path.exists():
        apply_media_from_tsv(model, media_path)
    
    # 부트스트랩 제공
    if bootstrap_mets:
        print(f"\n[부트스트랩 제공]")
        for met_id, amount in bootstrap_mets.items():
            try:
                met = model.metabolites.get_by_id(met_id)
                ex_id = f'EX_{met_id}'
                
                if ex_id not in [r.id for r in model.exchanges]:
                    ex_rxn = cobra.Reaction(ex_id)
                    ex_rxn.name = f'{met_id} exchange'
                    ex_rxn.lower_bound = -amount
                    ex_rxn.upper_bound = 1000
                    
                    met_e_id = met_id.replace('_c', '_e')
                    try:
                        met_e = model.metabolites.get_by_id(met_e_id)
                    except KeyError:
                        met_e = cobra.Metabolite(met_e_id, name=met.name, compartment='e')
                        model.add_metabolites([met_e])
                    
                    ex_rxn.add_metabolites({met_e: -1.0})
                    model.add_reactions([ex_rxn])
                else:
                    ex_rxn = model.reactions.get_by_id(ex_id)
                    ex_rxn.lower_bound = -amount
                    ex_rxn.upper_bound = 1000
                
                print(f"  {ex_id}: 하한={-amount}")
            except KeyError:
                print(f"  {met_id}: 메타볼라이트 없음")
    
    # ATPM=0 설정
    atpm_rxn = model.reactions.get_by_id('ATPM')
    atpm_rxn.lower_bound = 0
    atpm_rxn.upper_bound = 1000
    
    model.objective = 'Growth'
    solution = model.optimize()
    
    return solution

def check_acetate_pathway_step_by_step(model, media_path):
    """Acetate 경로를 단계별로 확인"""
    print("="*80)
    print("Acetate 경로 단계별 확인")
    print("="*80)
    
    # 1단계: Acetate uptake 확인
    print(f"\n[1단계] Acetate uptake 확인")
    model_test = model.copy()
    for rxn in model_test.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    if media_path.exists():
        apply_media_from_tsv(model_test, media_path)
    
    ex_ac = model_test.reactions.get_by_id('EX_ac_e')
    ex_ac.lower_bound = -1.0
    ex_ac.upper_bound = 1000
    
    # Acetate가 세포 내로 들어오는지 확인
    try:
        ac_c = model_test.metabolites.get_by_id('ac_c')
        ac_e = model_test.metabolites.get_by_id('ac_e')
        
        # Transport 반응 찾기
        ac_transports = []
        for rxn in ac_c.reactions:
            if ac_e in rxn.metabolites:
                ac_transports.append(rxn.id)
        
        print(f"  Acetate transport 반응: {ac_transports}")
        
        # Transport가 작동하는지 확인
        model_test.objective = 'EX_ac_e'
        model_test.reactions.get_by_id('EX_ac_e').objective_coefficient = -1.0
        solution = model_test.optimize()
        
        ac_flux = solution.fluxes.get('EX_ac_e', 0.0)
        print(f"  EX_ac_e 플럭스: {ac_flux:.6f}")
        
        if abs(ac_flux) < 1e-6:
            print(f"  [문제] Acetate uptake 불가")
            return False
        else:
            print(f"  [OK] Acetate uptake 가능")
    except KeyError:
        print(f"  [문제] ac_c 또는 ac_e 메타볼라이트 없음")
        return False
    
    # 2단계: ATP 생성 확인
    print(f"\n[2단계] ATP 생성 확인")
    model_test = model.copy()
    for rxn in model_test.exchanges:
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    if media_path.exists():
        apply_media_from_tsv(model_test, media_path)
    
    # ATP demand
    try:
        atp_c = model_test.metabolites.get_by_id('atp_c')
        atp_demand = cobra.Reaction('DM_atp_c')
        atp_demand.name = 'ATP demand'
        atp_demand.lower_bound = 0
        atp_demand.upper_bound = 1000
        atp_demand.add_metabolites({atp_c: -1.0})
        model_test.add_reactions([atp_demand])
        
        model_test.objective = 'DM_atp_c'
        solution = model_test.optimize()
        
        atp_production = solution.objective_value
        print(f"  ATP 최대 생산량: {atp_production:.6f}")
        
        if atp_production < 1e-6:
            print(f"  [문제] ATP 생성 불가")
            return False
        else:
            print(f"  [OK] ATP 생성 가능")
    except KeyError:
        print(f"  [문제] atp_c 메타볼라이트 없음")
        return False
    
    # 3단계: CoA 제공 시 ACS_ADP 작동 확인
    print(f"\n[3단계] CoA 제공 시 ACS_ADP 작동 확인")
    bootstrap = {'coa_c': 0.001}
    solution = test_with_bootstrap(model.copy(), media_path, bootstrap_mets=bootstrap)
    
    acs_adp_flux = solution.fluxes.get('ACS_ADP', 0.0)
    print(f"  ACS_ADP 플럭스: {acs_adp_flux:.6f}")
    
    if abs(acs_adp_flux) < 1e-6:
        print(f"  [문제] CoA 제공해도 ACS_ADP 작동 안 함")
        
        # ACS_ADP가 blocked인지 확인
        try:
            acs_adp = model.reactions.get_by_id('ACS_ADP')
            blocked = cobra.flux_analysis.find_blocked_reactions(model, open_exchanges=True)
            if acs_adp.id in blocked:
                print(f"  [원인] ACS_ADP가 blocked 상태")
            else:
                print(f"  [원인] ACS_ADP는 blocked 아님, 다른 문제")
                
                # ACS_ADP 반응식 확인
                print(f"  ACS_ADP 반응식: {acs_adp.reaction}")
                
                # 필요한 메타볼라이트 확인
                for met, coeff in acs_adp.metabolites.items():
                    if coeff < 0:  # 소모
                        print(f"    필요: {met.id} (계수: {coeff})")
                        # 최대 생산량 확인
                        try:
                            met_demand = cobra.Reaction(f'DM_{met.id}')
                            met_demand.lower_bound = 0
                            met_demand.upper_bound = 1000
                            met_demand.add_metabolites({met: -1.0})
                            model_test2 = model.copy()
                            model_test2.add_reactions([met_demand])
                            model_test2.objective = met_demand.id
                            sol_test = model_test2.optimize()
                            max_prod = sol_test.objective_value
                            print(f"      최대 생산량: {max_prod:.6f}")
                        except:
                            pass
        except:
            pass
        
        return False
    else:
        print(f"  [OK] CoA 제공 시 ACS_ADP 작동")
        return True
    
    # 4단계: Acetyl-CoA -> Citrate 확인
    print(f"\n[4단계] Acetyl-CoA -> Citrate 확인")
    bootstrap = {'accoa_c': 0.001, 'oaa_c': 0.001}
    solution = test_with_bootstrap(model.copy(), media_path, bootstrap_mets=bootstrap)
    
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
    media_path = base_path / "Stenotrophomonas" / "scenarios" / "YE0p5" / "Acetate_YE0p5__nocmnt__normalized.tsv"
    
    model = load_model(str(model_path))
    
    # ACS_ADP 추가 확인
    if 'ACS_ADP' not in [r.id for r in model.reactions]:
        print("[경고] ACS_ADP가 모델에 없습니다. 먼저 추가해주세요.")
        return
    
    print("="*80)
    print("실제 막힌 지점 찾기")
    print("="*80)
    
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
