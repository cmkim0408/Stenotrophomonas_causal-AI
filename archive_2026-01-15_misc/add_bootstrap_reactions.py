#!/usr/bin/env python
"""
부트스트랩 반응 추가
Acetate 기반 성장을 위한 최소 부트스트랩
"""

import cobra

def load_model(model_path="BaseModel.xml"):
    model = cobra.io.read_sbml_model(model_path)
    return model

def add_bootstrap_reactions(model):
    """부트스트랩 반응 추가"""
    print("부트스트랩 반응 추가 중...")
    
    bootstrap_components = {'atp_c': -0.05, 'coa_c': -0.01}
    
    added_reactions = []
    
    for met_id, supply_rate in bootstrap_components.items():
        try:
            met = model.metabolites.get_by_id(met_id)
            
            # 기존 demand 반응 확인
            dm_id = f'DM_{met_id}'
            try:
                dm_rxn = model.reactions.get_by_id(dm_id)
                print(f"  [UPDATE] {dm_id}: lower_bound를 {supply_rate}로 설정")
                dm_rxn.lower_bound = supply_rate
            except KeyError:
                # 새 반응 추가
                dm_rxn = cobra.Reaction(dm_id)
                dm_rxn.name = f'{met_id} bootstrap demand'
                dm_rxn.lower_bound = supply_rate
                dm_rxn.upper_bound = 1000
                dm_rxn.add_metabolites({met: -1})
                model.add_reactions([dm_rxn])
                print(f"  [ADD] {dm_id}: {supply_rate} mmol/gDCW/h")
                added_reactions.append(dm_id)
        except KeyError:
            print(f"  [SKIP] {met_id} metabolite 없음")
    
    print(f"\n총 {len(added_reactions)}개 반응 추가됨")
    return added_reactions

def main():
    model = load_model("BaseModel.xml")
    
    print("="*70)
    print("부트스트랩 반응 추가")
    print("="*70)
    
    added = add_bootstrap_reactions(model)
    
    # 모델 저장
    output_path = "BaseModel_with_bootstrap.xml"
    cobra.io.write_sbml_model(model, output_path)
    print(f"\n[OK] 모델 저장: {output_path}")
    
    print("="*70)

if __name__ == "__main__":
    main()
