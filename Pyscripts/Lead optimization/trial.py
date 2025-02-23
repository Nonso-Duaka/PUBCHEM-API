import json
import requests
from urllib.parse import quote_plus  # To handle encoding of SMARTS string

PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
GITHUB_JSON_URL = "https://raw.githubusercontent.com/durrantlab/autogrow4/master/autogrow/operators/mutation/smiles_click_chem/reaction_libraries/all_rxns/All_Rxns_rxn_library.json"

def fetch_substructures(smarts, max_results=2):
    """
    Fetch compounds from PubChem that match the given SMARTS, returning CIDs.
    Handles JSON decode errors gracefully.
    """
    encoded_smarts = quote_plus(smarts)
    search_url = f"{PUBCHEM_BASE_URL}/compound/fastsubstructure/smarts/{encoded_smarts}/cids/JSON"
    
    try:
        response = requests.get(search_url, timeout=10)  # Set a timeout to prevent hanging requests
        response.raise_for_status()  # Raise an error for HTTP errors (404, 500, etc.)

        # ✅ Check if response is JSON before parsing
        if "application/json" not in response.headers.get("Content-Type", ""):
            print(f"Warning: Non-JSON response for SMARTS '{smarts}'. Response content:\n{response.text[:500]}")
            return []  # Skip this SMARTS

        data = response.json()  # Safe JSON decoding
        cids = data.get("IdentifierList", {}).get("CID", [])
        return cids[:max_results]  # Return up to `max_results` CIDs

    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON for SMARTS '{smarts}'. Response content:\n{response.text[:500]}")
        return []  # Skip this SMARTS and continue execution

    except requests.RequestException as e:
        print(f"Error fetching substructures for SMARTS '{smarts}': {e}")
        return []  # Skip and continue

def fetch_reaction_data():
    """
    Fetch the reaction JSON file from GitHub.
    """
    try:
        response = requests.get(GITHUB_JSON_URL)
        response.raise_for_status()
        return response.json()
    except requests.RequestException as e:
        print(f"Error fetching reaction data: {e}")
        return {}

def process_reactions():
    """
    Fetch reactions from GitHub, extract reactants (SMARTS), and trigger substructure searches.
    """
    reactions = fetch_reaction_data()
    reaction_data = {}
    
    print("Processing reactions...")

    for reaction_name, details in reactions.items():
        reaction_string = details.get("reaction_string", "")
        if not reaction_string:
            print(f"Skipping reaction '{reaction_name}' due to missing reaction string.")
            continue

        try:
            reactants, _ = reaction_string.split(">>")
            reactants = reactants.split(".")
        except ValueError as e:
            print(f"Error parsing reaction string for '{reaction_name}': {e}")
            continue
        
        print(f"Processing reaction: {reaction_name}")
        
        substructure_results = {}
        for reactant in reactants:
            print(f"  - Fetching substructures for reactant: {reactant}")
            substructures = fetch_substructures(reactant)
            substructure_results[reactant] = substructures
        
        reaction_data[reaction_name] = {
            "reaction_string": reaction_string,
            "substructures": substructure_results
        }
    
    # Save results to file
    with open("filtered_reactions.json", "w") as output_file:
        json.dump(reaction_data, output_file, indent=4)
    
    print("Processing complete. Filtered reactions saved in 'filtered_reactions.json'.")
    return reaction_data

if __name__ == "__main__":
    results = process_reactions()
