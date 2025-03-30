import requests
import json
from time import sleep
from urllib.parse import quote
import os

class PubChemSMARTSQuery:
    def __init__(self):
        self.base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        self.delay = 1.5  # Delay between requests to avoid rate limiting
        self.max_results = 50  # Limit number of compounds per query
        self.timeout = 30  # Timeout for requests in seconds
        
    def get_reaction_data(self):
        """Fetch reaction data from AutoGrow4 repository"""
        url = "https://raw.githubusercontent.com/durrantlab/autogrow4/master/autogrow/operators/mutation/smiles_click_chem/reaction_libraries/all_rxns/All_Rxns_rxn_library.json"
        print("Fetching reaction data from GitHub...")
        try:
            response = requests.get(url, timeout=self.timeout)
            if response.status_code != 200:
                raise Exception(f"Failed to fetch reaction data: Status {response.status_code}")
            return response.json()
        except requests.exceptions.RequestException as e:
            raise Exception(f"Error fetching reaction data: {str(e)}")
    
    def search_by_smarts(self, smarts_pattern):
        """Search PubChem for compounds matching a SMARTS pattern using POST"""
        url = f"{self.base_url}/compound/fastsubstructure/smarts/cids/JSON"
        
        print(f"\nSearching for pattern: {smarts_pattern}")
        
        # Use POST request with data parameter
        data = {'smarts': smarts_pattern}
        headers = {'Content-Type': 'application/x-www-form-urlencoded'}
        
        try:
            response = requests.post(url, data=data, headers=headers, timeout=self.timeout)
            sleep(self.delay)  # Sleep to avoid rate limiting
            
            if response.status_code != 200:
                print(f"Error: Status code {response.status_code}")
                print(f"Response text: {response.text[:500]}")  # Print first 500 chars of response
                return []
                
            try:
                result = response.json()
                cids = result.get('IdentifierList', {}).get('CID', [])
                print(f"Found {len(cids)} matching compounds")
                
                # Limit results
                return cids[:self.max_results] if len(cids) > self.max_results else cids
            except json.JSONDecodeError as e:
                print(f"Error decoding JSON response: {str(e)}")
                print(f"Response content: {response.text[:500]}...")  # Print first 500 chars
                return []
                
        except requests.exceptions.RequestException as e:
            print(f"Request failed: {str(e)}")
            return []
    
    def get_smiles_for_cids(self, cids):
        """Get SMILES for a list of compound IDs"""
        if not cids:
            return []
            
        cid_str = ",".join(map(str, cids))
        url = f"{self.base_url}/compound/cid/{cid_str}/property/CanonicalSMILES/JSON"
        
        try:
            response = requests.get(url, timeout=self.timeout)
            sleep(self.delay)  # Sleep to avoid rate limiting
            
            if response.status_code != 200:
                print(f"Error: Status code {response.status_code}")
                print(f"Response text: {response.text[:500]}")
                return []
                
            results = []
            try:
                for prop in response.json().get('PropertyTable', {}).get('Properties', []):
                    results.append({
                        'CID': prop.get('CID'),
                        'SMILES': prop.get('CanonicalSMILES', '')
                    })
                return results
            except Exception as e:
                print(f"Error processing SMILES data: {str(e)}")
                print(f"Response content: {response.text[:500]}...")
                return []
                
        except requests.exceptions.RequestException as e:
            print(f"Request failed: {str(e)}")
            return []
    
    def query_all_reactions(self, output_file="reaction_fragments.json"):
        """Query PubChem for all reaction SMARTS patterns"""
        # Get reaction data
        try:
            reaction_data = self.get_reaction_data()
        except Exception as e:
            print(f"Fatal error: {str(e)}")
            return {}
        
        results = {}
        total_reactions = len(reaction_data)
        processed = 0
        
        # Process each reaction
        for rxn_name, rxn_info in reaction_data.items():
            processed += 1
            if 'group_smarts' not in rxn_info or not rxn_info['group_smarts']:
                print(f"\n[{processed}/{total_reactions}] Skipping {rxn_name}: No group_smarts found")
                continue
                
            print(f"\n[{processed}/{total_reactions}] Processing reaction: {rxn_name}")
            
            results[rxn_name] = {
                'reaction_name': rxn_info.get('reaction_name', rxn_name),
                'group_smarts': rxn_info['group_smarts'],
                'reactants': []
            }
            
            # Process each SMARTS pattern for this reaction
            for i, smarts in enumerate(rxn_info['group_smarts']):
                print(f"\nProcessing reactant {i+1}/{len(rxn_info['group_smarts'])}")
                print(f"SMARTS pattern: {smarts}")
                
                # Search for compounds
                cids = self.search_by_smarts(smarts)
                
                if cids:
                    # Get SMILES for these CIDs
                    compounds = self.get_smiles_for_cids(cids)
                    
                    results[rxn_name]['reactants'].append({
                        'smarts': smarts,
                        'compounds': compounds,
                        'num_compounds': len(compounds)
                    })
                    
                    print(f"Found {len(compounds)} compounds for this reactant")
                else:
                    results[rxn_name]['reactants'].append({
                        'smarts': smarts,
                        'compounds': [],
                        'num_compounds': 0
                    })
                    print(f"No compounds found for this reactant")
                
                # Save intermediate results after each reactant
                self._save_intermediate_results(results, output_file)
        
        # Final save to JSON file
        self._save_final_results(results, output_file)
        return results
    
    def _save_intermediate_results(self, results, output_file):
        """Save intermediate results to a temporary file"""
        temp_file = f"temp_{output_file}"
        try:
            with open(temp_file, 'w') as f:
                json.dump(results, f, indent=2)
        except Exception as e:
            print(f"Warning: Could not save intermediate results: {str(e)}")
    
    def _save_final_results(self, results, output_file):
        """Save final results to JSON file"""
        try:
            with open(output_file, 'w') as f:
                json.dump(results, f, indent=2)
            print(f"\nResults saved to {os.path.abspath(output_file)}")
        except Exception as e:
            print(f"Error saving final results: {str(e)}")

# Run the script
if __name__ == "__main__":
    print("Starting PubChem query using reaction SMARTS patterns...")
    querier = PubChemSMARTSQuery()
    querier.query_all_reactions()