#!/usr/bin/env python
"""Bootstrap 테스트: 소량 ATP 제공"""

import cobra
from cobra import Reaction

model = cobra.io.read_sbml_model('BaseModel.xml')

# Medium 설정
for rxn in model.exchanges:
    rxn.lower_bound = 0
    rxn.upper_bound = 0

# Acetate 및 필수 영양소
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

# ATP bootstrap
atp_c = model.metabolites.get_by_id('atp_c')
atp_boot = Reaction('EX_atp_boot')
atp_boot.name = 'ATP bootstrap'
atp_boot.lower_bound = -0.1
atp_boot.upper_bound = 0
atp_boot.add_metabolites({atp_c: -1})
model.add_reactions([atp_boot])

model.objective = 'Growth'
sol = model.optimize()

print(f"Bootstrap 테스트 결과:")
print(f"  상태: {sol.status}")
print(f"  Biomass flux: {sol.objective_value:.6f}")

if sol.objective_value > 1e-6:
    print(f"\n[SUCCESS] Bootstrap으로 성장 가능!")
    print(f"주요 플럭스:")
    for rxn_id in ['EX_ac_e', 'ACt', 'ACS', 'CS', 'ICL', 'MALS', 'Growth']:
        flux = sol.fluxes.get(rxn_id, 0)
        if abs(flux) > 1e-6:
            print(f"  {rxn_id}: {flux:.4f}")
else:
    print(f"\n[FAIL] Bootstrap으로도 성장 불가")

