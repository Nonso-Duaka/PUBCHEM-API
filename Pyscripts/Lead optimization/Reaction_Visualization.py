import os
import requests
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import json
import pandas as pd
from PIL import Image  # For saving images

# GitHub repository URL containing the reaction library
Git_repo_url = "https://raw.githubusercontent.com/durrantlab/autogrow4/master/autogrow/operators/mutation/smiles_click_chem/reaction_libraries/all_rxns/All_Rxns_rxn_library.json"

# Create a directory to store reaction images
Mol_img_dir = "reaction_images"
os.makedirs(Mol_img_dir, exist_ok=True)

# File names for storing the results in JSON and CSV formats
json_file = "reaction_dir.json"
csv_file = "reaction_dir.csv"

# Fetch the reaction data from the GitHub repository
response = requests.get(Git_repo_url)
if response.status_code == 200:
    rxn_data = response.json()
else:
    print("Failed to retrieve the reaction library")
    exit()

# Initialize an empty list to store reaction data
Rxn_list = []

def generate_rxn_image(reaction_name, smirks, reactants, idx, Mol_img_dir):
    """
    Generates and saves an image for a chemical reaction based on SMIRKS and reactants.
    Parameters:
    - reaction_name (str): The name of the reaction.
    - smirks (str): The SMIRKS pattern representing the reaction.
    - reactants (list): A list of SMILES strings representing the reactants.
    - idx (int): The index of the reaction in the list.
    - Mol_img_dir (str): The directory path to save the generated reaction images.
    Returns:
    - str: The file path where the reaction image was saved, or None if the image could not be generated.
    """
    try:
        print(f"Generating reaction image for {reaction_name} with SMIRKS: {smirks} and reactants: {reactants}")

        # Generate reaction object from SMIRKS
        reaction = AllChem.ReactionFromSmarts(smirks)
        if reaction is None:
            print(f"ERROR: Failed to generate reaction object for {reaction_name}")
            return None

        # Parse reactant SMILES into molecules
        reactant_mols = [Chem.MolFromSmiles(smiles) for smiles in reactants]
        if None in reactant_mols:
            print(f"ERROR: Failed to parse reactants for {reaction_name}. Check input: {reactants}")
            return None

        # Run the reaction on the reactants
        products = reaction.RunReactants(tuple(reactant_mols))
        if not products or not products[0]:
            print(f"ERROR: Failed to run reactants for {reaction_name}. Reactants: {reactants}")
            return None

        # Create a full reaction object
        full_reaction = AllChem.ChemicalReaction()
        for mol in reactant_mols:
            full_reaction.AddReactantTemplate(mol)
        for mol in products[0]:
            full_reaction.AddProductTemplate(mol)
        full_reaction.Initialize()

        # Generate and save the reaction image
        img = Draw.ReactionToImage(full_reaction, subImgSize=(400, 400))
        img_path = os.path.join(Mol_img_dir, f"reaction_{idx+1}.png")
        img.save(img_path)
        print(f"Processed {reaction_name} and saved image to {img_path}")
        return img_path

    except Exception as e:
        print(f"ERROR: Failed to process {reaction_name}: {e}")
        return None

# Process each reaction in the fetched reaction data
for idx, (reaction_name, reaction_info) in enumerate(rxn_data.items()):
    try:
        smirks = reaction_info.get('reaction_string')
        default_reactants = reaction_info.get('example_rxn_reactants', [])

        print(f"\nProcessing {reaction_name}")

        # Debug: Print available default reactants
        print(f"Default reactants for {reaction_name}: {default_reactants}")

        # Allow the user to input reactant SMILES or use default ones
        user_input = input("Enter reactant SMILES (comma-separated) or press enter to use default: ").strip()
        reactants = user_input.split(',') if user_input else default_reactants

        # Debug: Show selected reactants
        print(f"Using reactants: {reactants}")

        # Validate SMIRKS and reactants
        if not smirks:
            print(f"Skipping {reaction_name}: Missing SMIRKS")
            continue
        if not reactants:
            print(f"Skipping {reaction_name}: Missing reactants")
            continue

        # Generate the reaction image and save it
        img_path = generate_rxn_image(reaction_name, smirks, reactants, idx, Mol_img_dir)
        if img_path:
            # Append the reaction data to the list
            Rxn_list.append({
                "reaction_name": reaction_name,
                "smirks": smirks,
                "reactants": reactants,
                "image_path": img_path
            })

    except Exception as e:
        
        print(f"ERROR: Failed to process {reaction_name}: {e}")

# Save the results to JSON and CSV files
if Rxn_list:
    with open(json_file, "w") as f:
        json.dump(Rxn_list, f, indent=4)
    df = pd.DataFrame(Rxn_list)
    df.to_csv(csv_file, index=False)
    print(f"Saved reaction data to {json_file} and {csv_file}")
else:
    print("No valid reactions were processed.")