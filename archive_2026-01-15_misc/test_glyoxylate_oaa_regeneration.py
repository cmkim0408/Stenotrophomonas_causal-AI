#!/usr/bin/env python
"""
Glyoxylate shunt를 통한 OAA 재생성 확인
아세트산만으로 성장하려면 OAA를 재생성할 수 있어야 함
"""

import cobra
from cobra import Reaction

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    return model

def setup_medium(model, bootstrap=None):
    """Medium 설정"""
    for rxn in model.exchanges:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    
    for ex_id, bounds in [
        ('EX_ac_e', (-1000, 1000)),
        ('EX_nh4_e', (-1000, 0)),
        ('EX_h2o_e', (-1000, 0)),
        ('EX_h_e', (-1000, 0)),
        ('EX_pi_e', (-1000, 0)),
        ('EX_so4_e', (-1000, 0)),
        ('EX_k_e', (-1000, 0)),
        ('EX_mg2_e', (-1000, 0)),
        ('EX_fe2_e', (-1000, 0)),
        ('EX_o2_e', (-1000, 1000))
    ]:
        try:
            ex_rxn = model.reactions.get_by_id(ex_id)
            ex_rxn.lower_bound = bounds[0]
            ex_rxn.upper_bound = bounds[1]
        except KeyError:
            pass
    
    # Bootstrap 대사물질 제공
    boot_rxns = []
    if bootstrap:
        for met_id, amount in bootstrap.items():
            met = model.metabolites.get_by_id(met_id)
            boot = Reaction(f'EX_{met_id}_boot')
            boot.lower_bound = -amount
            boot.add_metabolites({met: -1})
            model.add_reactions([boot])
            boot_rxns.append(boot)
    
    return model, boot_rxns

def analyze_glyoxylate_shunt(model):
    """Glyoxylate shunt 분석"""
    print("="*70)
    print("Glyoxylate Shunt 분석")
    print("="*70)
    
    try:
        icl = model.reactions.get_by_id('ICL')
        mals = model.reactions.get_by_id('MALS')
        
        print(f"ICL (Isocitrate Lyase): {icl.reaction}")
        print(f"MALS (Malate Synthase): {mals.reaction}")
        
        print(f"\n[분석]")
        print(f"  1. ICL: Isocitrate → Glyoxylate + Succinate")
        print(f"  2. MALS: Glyoxylate + Acetyl-CoA → Malate")
        print(f"  3. MDH: Malate → OAA")
        print(f"  → 결과: Isocitrate + Acetyl-CoA → Succinate + OAA")
        print(f"  → OAA가 재생성됨!")
        
        # 필요한 대사물질 확인
        print(f"\n필요한 초기 대사물질:")
        print(f"  - Acetyl-CoA (ACS로 생성 가능)")
        print(f"  - Isocitrate (Acetyl-CoA + OAA → Citrate → Isocitrate)")
        print(f"  → 초기 OAA가 필요함!")
        
    except KeyError as e:
        print(f"[ERROR] {e}")

def test_oaa_regeneration(model):
    """OAA 재생성 테스트"""
    print("\n" + "="*70)
    print("OAA 재생성 테스트 (초기 OAA 제공)")
    print("="*70)
    
    # 초기 OAA 제공
    model, boot_rxns = setup_medium(model, bootstrap={'oaa_c': 0.01, 'coa_c': 0.01})
    
    # OAA 생산 테스트 (순환 가능 여부)
    oaa_c = model.metabolites.get_by_id('oaa_c')
    test_oaa = Reaction('TEST_oaa')
    test_oaa.lower_bound = 0
    test_oaa.upper_bound = 1000
    test_oaa.add_metabolites({oaa_c: -1})
    model.add_reactions([test_oaa])
    model.objective = 'TEST_oaa'
    
    solution = model.optimize()
    
    print(f"초기 OAA + CoA 제공 후 OAA 재생성:")
    print(f"  상태: {solution.status}")
    print(f"  OAA 생산: {solution.objective_value:.6f}")
    
    if solution.objective_value > 0.01:
        print(f"  [SUCCESS] OAA 재생성 가능!")
        print(f"\n  활성화된 경로:")
        for rxn_id in ['EX_ac_e', 'ACt', 'ACS', 'CS', 'ACONT', 'ICL', 'MALS', 'MDH']:
            flux = solution.fluxes.get(rxn_id, 0)
            if abs(flux) > 1e-6:
                print(f"    {rxn_id}: {flux:.4f}")
    else:
        print(f"  [FAIL] OAA 재생성 불가")
        print(f"\n  원인 분석:")
        for rxn_id in ['EX_ac_e', 'ACt', 'ACS', 'CS', 'ICL', 'MALS', 'MDH']:
            flux = solution.fluxes.get(rxn_id, 0)
            print(f"    {rxn_id}: {flux:.6f}")
    
    model.remove_reactions([test_oaa] + boot_rxns)

def test_complete_acetate_growth(model):
    """완전한 아세트산 성장 테스트"""
    print("\n" + "="*70)
    print("완전한 아세트산 성장 테스트 (초기 OAA + CoA 제공)")
    print("="*70)
    
    # 최소한의 bootstrap: OAA와 CoA만 제공
    model, boot_rxns = setup_medium(model, bootstrap={'oaa_c': 0.01, 'coa_c': 0.01})
    
    model.objective = 'Growth'
    solution = model.optimize()
    
    print(f"결과:")
    print(f"  상태: {solution.status}")
    print(f"  Biomass: {solution.objective_value:.6f}")
    
    if solution.objective_value > 1e-6:
        print(f"  [SUCCESS] 아세트산만으로 성장 가능!")
        print(f"\n  주요 경로 플럭스:")
        pathway = [
            ('EX_ac_e', 'Acetate uptake'),
            ('ACt', 'Acetate transport'),
            ('ACS', 'Acetyl-CoA synthetase'),
            ('CS', 'Citrate synthase'),
            ('ICL', 'Isocitrate lyase'),
            ('MALS', 'Malate synthase'),
            ('MDH', 'Malate dehydrogenase'),
            ('ME1', 'Malic enzyme 1'),
            ('ME2', 'Malic enzyme 2'),
            ('Growth', 'Biomass')
        ]
        for rxn_id, name in pathway:
            flux = solution.fluxes.get(rxn_id, 0)
            if abs(flux) > 1e-6:
                print(f"    {rxn_id:10s} ({name:25s}): {flux:.6f}")
        
        print(f"\n  [결론] 초기 OAA + CoA만 제공하면 순환이 시작되어")
        print(f"         아세트산만으로 지속적인 성장 가능!")
    else:
        print(f"  [FAIL] 성장 불가")
        print(f"\n  상세 분석:")
        for rxn_id in ['EX_ac_e', 'ACt', 'ACS', 'CS', 'ICL', 'MALS', 'MDH', 'Growth']:
            flux = solution.fluxes.get(rxn_id, 0)
            print(f"    {rxn_id}: {flux:.6f}")
    
    model.remove_reactions(boot_rxns)

def main():
    model = load_model("BaseModel.xml")
    
    analyze_glyoxylate_shunt(model)
    test_oaa_regeneration(model)
    test_complete_acetate_growth(model)
    
    print("\n" + "="*70)
    print("최종 결론")
    print("="*70)
    print("아세트산만으로 ATP를 만들 수 있는가?")
    print("→ 이론적으로는 가능하지만, 초기 조건이 필요:")
    print("  1. CoA (초기화)")
    print("  2. OAA (초기화 또는 재생성 경로)")
    print("")
    print("Glyoxylate shunt를 통해 OAA를 재생성할 수 있지만,")
    print("처음에는 초기 OAA가 필요합니다.")
    print("")
    print("실제 생물학적으로는 세포에 초기 대사물질이 존재하므로,")
    print("모델에서도 초기 OAA와 CoA를 제공하면 아세트산만으로")
    print("지속적인 ATP 생산 및 성장이 가능해야 합니다.")
    print("="*70)

if __name__ == "__main__":
    main()

