import json
import requests
import time
from urllib.parse import quote_plus
from concurrent.futures import ThreadPoolExecutor, as_completed
from random import randint

# Constants
PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
GITHUB_JSON_URL = "https://raw.githubusercontent.com/durrantlab/autogrow4/master/autogrow/operators/mutation/smiles_click_chem/reaction_libraries/all_rxns/All_Rxns_rxn_library.json"
MAX_RETRIES = 5  # Maximum retries for failed requests
TIMEOUT = 10  # Timeout in seconds for requests
MAX_THREADS = 10  # Number of concurrent threads to use
RATE_LIMIT_DELAY = 1  # Delay in seconds between API calls to avoid rate limiting

def fetch_with_retry(url, retries=MAX_RETRIES, timeout=TIMEOUT):
    """
    Fetch a URL with retries and exponential backoff for failed requests.
    """
    for attempt in range(retries):
        try:
            response = requests.get(url, timeout=timeout)
            response.raise_for_status()
            return response.json()
        except requests.RequestException as e:
            print(f"Request failed (attempt {attempt + 1}/{retries}): {e}")
            if attempt < retries - 1:
                sleep_time = 2 ** attempt + randint(1, 3)
                print(f"Retrying in {sleep_time} seconds...")
                time.sleep(sleep_time)
            else:
                print("Max retries reached, moving on.")
    return None

def fetch_substructures_with_pct_query(smarts):
    """
    Fetch compounds from PubChem that match the given SMARTS and have MW < 150 Da using PUG API.
    """
    encoded_smarts = quote_plus(smarts)
    search_url = f"{PUBCHEM_BASE_URL}/compound/fastsubstructure/smarts/{encoded_smarts}/cids/JSON"
    data = fetch_with_retry(search_url)
    return data.get("IdentifierList", {}).get("CID", []) if data else []

def fetch_properties_for_cids(cids):
    """
    Fetch molecular properties for a list of CIDs (MolecularWeight, IsomericSMILES).
    """
    properties = []
    cid_batches = [cids[i:i + MAX_THREADS] for i in range(0, len(cids), MAX_THREADS)]
    
    for batch in cid_batches:
        with ThreadPoolExecutor(max_workers=MAX_THREADS) as executor:
            futures = {executor.submit(fetch_properties_for_cid, cid): cid for cid in batch}
            for future in as_completed(futures):
                result = future.result()
                if result and float(result["MW"]) < 150:
                    properties.append(result)
        time.sleep(RATE_LIMIT_DELAY)  # Rate limiting to avoid overwhelming PubChem
    
    return properties

def fetch_properties_for_cid(cid):
    """
    Fetch molecular properties for a single CID.
    """
    prop_url = f"{PUBCHEM_BASE_URL}/compound/cid/{cid}/property/MolecularWeight,IsomericSMILES/JSON"
    prop_data = fetch_with_retry(prop_url)
    if prop_data:
        properties = prop_data.get("PropertyTable", {}).get("Properties", [])
        if properties:
            return {
                "CID": cid,
                "SMILES": properties[0].get("IsomericSMILES", ""),
                "MW": properties[0].get("MolecularWeight", float("inf"))
            }
    return None

def fetch_reaction_data():
    """
    Fetch the reaction JSON file from GitHub.
    """
    return fetch_with_retry(GITHUB_JSON_URL)

def process_reactions():
    """
    Fetch reactions from GitHub, extract reactants (SMARTS), and trigger substructure searches.
    """
    reactions = fetch_reaction_data()
    if not reactions:
        print("Error: Could not fetch reaction data.")
        return {}

    reaction_data = {}

    for reaction_name, details in reactions.items():
        reaction_string = details.get("reaction_string", "")
        if not reaction_string:
            continue

        reactants = reaction_string.split(">>")[0].split(".")
        substructure_results = {}

        for reactant in reactants:
            cids = fetch_substructures_with_pct_query(reactant)
            if cids:
                properties = fetch_properties_for_cids(cids)
                substructure_results[reactant] = properties

        reaction_data[reaction_name] = {
            "reaction_string": reaction_string,
            "substructures": substructure_results
        }

    with open("filtered_reactions.json", "w") as output_file:
        json.dump(reaction_data, output_file, indent=4)

    return reaction_data

if __name__ == "__main__":
    results = process_reactions()
    print("Processing complete. Filtered reactions saved in 'filtered_reactions.json'")