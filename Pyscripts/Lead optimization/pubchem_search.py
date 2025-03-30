import requests
import json
from time import sleep
import re
import os
from urllib.parse import quote

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
        """Search PubChem for compounds matching a SMARTS pattern using POST.
        Returns ALL CIDs found, even if JSON parsing fails."""
        url = f"{self.base_url}/compound/fastsubstructure/smarts/cids/JSON"
        
        print(f"\nSearching for pattern: {smarts_pattern}")
        
        # Use POST request with data parameter
        data = {'smarts': smarts_pattern}
        headers = {'Content-Type': 'application/x-www-form-urlencoded'}
        
        try:
            response = requests.post(url, data=data, headers=headers, timeout=self.timeout)
            sleep(self.delay)
            
            if response.status_code != 200:
                print(f"Error: Status code {response.status_code}")
                print(f"Response text: {response.text[:500]}")
                return []
                
            # Save raw response for debugging
            with open("last_response.txt", "w") as f:
                f.write(response.text)
            
            # List to store ALL found CIDs (from both JSON and text extraction)
            all_cids = []
            
            # 1. First try to parse as complete JSON
            try:
                result = response.json()
                json_cids = result.get('IdentifierList', {}).get('CID', [])
                print(f"Found {len(json_cids)} CIDs via JSON parse")
                all_cids.extend(json_cids)
            except json.JSONDecodeError:
                print("Warning: JSON parse failed (will rely on text extraction)")
            
            # 2. ALWAYS try to extract CIDs from text (even if JSON parse worked)
            text_cids = self._extract_cids_from_text(response.text)
            print(f"Extracted {len(text_cids)} CIDs from text response")
            all_cids.extend(text_cids)
            
            # Remove duplicates while preserving order
            seen = set()
            unique_cids = [x for x in all_cids if not (x in seen or seen.add(x))]
            
            print(f"Total unique CIDs found: {len(unique_cids)}")
            return unique_cids[:self.max_results] if len(unique_cids) > self.max_results else unique_cids
                
        except requests.exceptions.RequestException as e:
            print(f"Request failed: {str(e)}")
            return []
    
    def _extract_cids_from_text(self, response_text):
        """Extract CIDs from raw response text using multiple patterns"""
        cids = []
        # Multiple patterns to catch different formats
        patterns = [
            r'"CID":\s*\[([\d,\s]+)',      # Standard JSON format
            r'"CID"\s*:\s*\[([\d\s,]+)',   # Variant with different spacing
            r'\[(\d+(?:\s*,\s*\d+)*)\]',    # Bare arrays
            r'\b\d{3,}\b'                   # Standalone large numbers (likely CIDs)
        ]
        
        for pattern in patterns:
            matches = re.findall(pattern, response_text)
            for match in matches:
                numbers = re.findall(r'\d+', str(match))
                cids.extend([int(cid) for cid in numbers])
        
        return cids
    
    def get_smiles_for_cids(self, cids):
        """Get SMILES for a list of compound IDs"""
        if not cids:
            return []
            
        cid_str = ",".join(map(str, cids))
        url = f"{self.base_url}/compound/cid/{cid_str}/property/CanonicalSMILES/JSON"
        
        try:
            response = requests.get(url, timeout=self.timeout)
            sleep(self.delay)
            
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
        try:
            reaction_data = self.get_reaction_data()
        except Exception as e:
            print(f"Fatal error: {str(e)}")
            return {}
        
        results = {}
        total_reactions = len(reaction_data)
        processed = 0
        
        for rxn_name, rxn_info in reaction_data.items():
            processed += 1
            if 'group_smarts' not in rxn_info or not rxn_info['group_smarts']:
                print(f"\n[{processed}/{total_reactions}] Skipping {rxn_name}: No group_smarts")
                continue
                
            print(f"\n[{processed}/{total_reactions}] Processing reaction: {rxn_name}")
            
            results[rxn_name] = {
                'reaction_name': rxn_info.get('reaction_name', rxn_name),
                'group_smarts': rxn_info['group_smarts'],
                'reactants': []
            }
            
            for i, smarts in enumerate(rxn_info['group_smarts']):
                print(f"\nProcessing reactant {i+1}/{len(rxn_info['group_smarts'])}")
                print(f"SMARTS pattern: {smarts}")
                
                cids = self.search_by_smarts(smarts)
                
                if cids:
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
                    print("No compounds found for this reactant")
                
                self._save_intermediate_results(results, output_file)
        
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

if __name__ == "__main__":
    print("Starting PubChem query using reaction SMARTS patterns...")
    querier = PubChemSMARTSQuery()
    querier.query_all_reactions()