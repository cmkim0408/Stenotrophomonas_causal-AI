#!/usr/bin/env python
import cobra

m = cobra.io.read_sbml_model('BaseModel.xml')
rxns = ['ICL', 'MALS', 'ICDHx', 'CS']
for rxn_id in rxns:
    if rxn_id in [r.id for r in m.reactions]:
        r = m.reactions.get_by_id(rxn_id)
        print(f'{rxn_id}: bounds=[{r.lower_bound}, {r.upper_bound}]')
        print(f'  reaction: {r.reaction[:100]}')
    else:
        print(f'{rxn_id}: 반응 없음')
