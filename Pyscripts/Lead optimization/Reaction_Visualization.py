import os
import requests
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import json
import pandas as pd

Git_repo_url =  "https://raw.githubusercontent.com/durrantlab/autogrow4/master/autogrow/operators/mutation/smiles_click_chem/reaction_libraries/all_rxns/All_Rxns_rxn_library.json"

# Make an image directory
Mol_img_dir = "reaction_images"
os.makedirs(Mol_img_dir, exist_ok=True)
json_file = "reaction_dir.json"
csv_file = "reaction_dir.csv"

# Fetch the data
response = requests.get(Git_repo_url)
if response.status_code == 200:
    rxn_data = response.json()
else:
    print("Failed to retrieve the reaction library")
    exit()

Rxn_list = []

def generate_rxn_image(reaction_name, smirks, reactants, idx, Mol_img_dir):
    try:
        reaction = AllChem.ReactionFromSmarts(smirks)
        if reaction is None:
            print(f"Failed to generate reaction object for {reaction_name}")
            return None
        
        reactant_mols = [Chem.MolFromSmiles(smiles) for smiles in reactants]
        if None in reactant_mols:
            print(f"Failed to parse reactants for {reaction_name}")
            return None
        
        products = reaction.RunReactants(tuple(reactant_mols))
        if not products or not products[0]:
            print(f"Failed to run reactants to get product for {reaction_name}")
            return None
        
        product_mols = products[0]
        full_reaction = AllChem.ChemicalReaction()
        for mol in reactant_mols:
            full_reaction.AddReactantTemplate(mol)
        for mol in product_mols:
            full_reaction.AddProductTemplate(mol)
        full_reaction.Initialize()
        
        img = Draw.ReactionToImage(full_reaction, subImgSize=(400, 400))
        img_path = os.path.join(Mol_img_dir, f"reaction_{idx+1}.png")
        img.save(img_path)
        print(f"Processed {reaction_name} and saved image to {img_path}")
        return img_path
    except Exception as e:
        print(f"Failed to process {reaction_name}: {e}")
        return None

for idx, (reaction_name, reaction_info) in enumerate(rxn_data.items()):
    try:
        smirks = reaction_info.get('reaction_string')
        default_reactants = reaction_info.get('example_rxn_reactants', [])
        
        print(f"\nProcessing {reaction_name}")
        user_input = input("Enter reactant SMILES (comma-separated) or press enter to use default: ").strip()
        reactants = user_input.split(',') if user_input else default_reactants
        
        if not smirks or not reactants:
            print(f"Skipping {reaction_name}: Missing SMIRKS or reactants")
            continue
        
        img_path = generate_rxn_image(reaction_name, smirks, reactants, idx, Mol_img_dir)
        if img_path:
            Rxn_list.append({
                "reaction_name": reaction_name,
                "smirks": smirks,
                "reactants": reactants,
                "image_path": img_path
            })
    except Exception as e:
        print(f"Failed to process {reaction_name}: {e}")

if Rxn_list:
    with open(json_file, "w") as f:
        json.dump(Rxn_list, f, indent=4)
    df = pd.DataFrame(Rxn_list)
    df.to_csv(csv_file, index=False)
    print(f"Saved reaction data to {json_file} and {csv_file}")
else:
    print("No valid reactions were processed.")
