import requests
import json
from time import sleep
import re
import os
from urllib.parse import quote

class PubChemSMARTSQuery:
    def __init__(self):
        self.base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        self.delay = 1.5
        self.timeout = 30
        self.max_mw = 150
        self.page_size = 100
        
    def get_reaction_data(self):
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
        url = f"{self.base_url}/compound/fastsubstructure/smarts/cids/JSON"
        print(f"\nSearching for pattern: {smarts_pattern}")
        data = {'smarts': smarts_pattern}
        headers = {'Content-Type': 'application/x-www-form-urlencoded'}
        try:
            count_url = f"{self.base_url}/compound/fastsubstructure/smarts/cids/TXT?smarts={quote(smarts_pattern)}"
            count_response = requests.get(count_url, timeout=self.timeout)
            sleep(self.delay)
            if count_response.status_code != 200:
                print(f"Error getting count: Status code {count_response.status_code}")
                print(f"Response text: {count_response.text[:500]}")
                return []
            all_cids = self._extract_cids_from_text(count_response.text)
            total_cids = len(all_cids)
            print(f"Total compounds found: {total_cids}")
            if not all_cids:
                return []
            filtered_cids = []
            for i in range(0, len(all_cids), self.page_size):
                page_cids = all_cids[i:i + self.page_size]
                print(f"Processing page {i//self.page_size + 1} with {len(page_cids)} compounds...")
                properties = self._get_compound_properties(page_cids)
                for cid, mw in properties:
                    try:
                        mw_float = float(mw) if isinstance(mw, str) else mw
                        if mw_float is not None and mw_float <= self.max_mw:
                            filtered_cids.append(cid)
                    except (ValueError, TypeError) as e:
                        print(f"Warning: Could not parse MW for CID {cid}: {str(e)}")
                        continue
                print(f"After MW filtering: {len(filtered_cids)} compounds so far")
            print(f"Final count after MW filtering: {len(filtered_cids)}/{total_cids} compounds")
            return filtered_cids
        except requests.exceptions.RequestException as e:
            print(f"Request failed: {str(e)}")
            return []
    
    def _get_compound_properties(self, cids):
        if not cids:
            return []
        cid_str = ",".join(map(str, cids))
        url = f"{self.base_url}/compound/cid/{cid_str}/property/MolecularWeight/JSON"
        try:
            response = requests.get(url, timeout=self.timeout)
            sleep(self.delay)
            if response.status_code != 200:
                print(f"Error getting properties: Status code {response.status_code}")
                print(f"Response text: {response.text[:500]}")
                return [(cid, None) for cid in cids]
            results = []
            try:
                data = response.json()
                if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
                    for prop in data['PropertyTable']['Properties']:
                        cid = prop.get('CID')
                        mw = prop.get('MolecularWeight')
                        if mw is not None:
                            try:
                                mw = float(mw)
                            except (ValueError, TypeError):
                                mw = None
                        results.append((cid, mw))
                return results
            except Exception as e:
                print(f"Error processing property data: {str(e)}")
                print(f"Response content: {response.text[:500]}...")
                return [(cid, None) for cid in cids]
        except requests.exceptions.RequestException as e:
            print(f"Request failed: {str(e)}")
            return [(cid, None) for cid in cids]
    
    def _extract_cids_from_text(self, response_text):
        cids = []
        patterns = [
            r'"CID":\s*\[([\d,\s]+)',
            r'"CID"\s*:\s*\[([\d\s,]+)',
            r'\[(\d+(?:\s*,\s*\d+)*)\]',
            r'\b\d{3,}\b'
        ]
        for pattern in patterns:
            matches = re.findall(pattern, response_text)
            for match in matches:
                numbers = re.findall(r'\d+', str(match))
                cids.extend([int(cid) for cid in numbers])
        seen = set()
        unique_cids = [x for x in cids if not (x in seen or seen.add(x))]
        return unique_cids
    
    def get_smiles_for_cids(self, cids):
        if not cids:
            return []
        results = []
        for i in range(0, len(cids), self.page_size):
            page_cids = cids[i:i + self.page_size]
            cid_str = ",".join(map(str, page_cids))
            url = f"{self.base_url}/compound/cid/{cid_str}/property/CanonicalSMILES,IsomericSMILES/JSON"
            try:
                response = requests.get(url, timeout=self.timeout)
                sleep(self.delay)
                if response.status_code != 200:
                    print(f"Error: Status code {response.status_code}")
                    print(f"Response text: {response.text[:500]}")
                    continue
                try:
                    data = response.json()
                    if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
                        for prop in data['PropertyTable']['Properties']:
                            results.append({
                                'CID': prop.get('CID'),
                                'CanonicalSMILES': prop.get('CanonicalSMILES', ''),
                                'IsomericSMILES': prop.get('IsomericSMILES', '')
                            })
                except Exception as e:
                    print(f"Error processing SMILES data: {str(e)}")
                    print(f"Response content: {response.text[:500]}...")
            except requests.exceptions.RequestException as e:
                print(f"Request failed: {str(e)}")
        return results
    
    def query_all_reactions(self, output_file="reaction_fragments.json"):
        try:
            reaction_data = self.get_reaction_data()
        except Exception as e:
            print(f"Fatal error: {str(e)}")
            return {}
        
        temp_file = f"temp_{output_file}"
        if os.path.exists(temp_file):
            with open(temp_file, 'r') as f:
                results = json.load(f)
            print(f"Resuming from saved progress in {temp_file}...")
        else:
            results = {}
        
        total_reactions = len(reaction_data)
        processed = 0
        
        try:
            for rxn_name, rxn_info in reaction_data.items():
                processed += 1
                if rxn_name in results:
                    print(f"\n[{processed}/{total_reactions}] Skipping already processed reaction: {rxn_name}")
                    continue
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
                        print(f"Found {len(compounds)} compounds for this reactant (MW ≤ {self.max_mw})")
                    else:
                        results[rxn_name]['reactants'].append({
                            'smarts': smarts,
                            'compounds': [],
                            'num_compounds': 0
                        })
                        print("No compounds found for this reactant")
                    self._save_intermediate_results(results, output_file)
        except KeyboardInterrupt:
            print("\nInterrupted by user. Progress saved.")
            self._save_final_results(results, output_file)
            return results
        
        self._save_final_results(results, output_file)
        return results

    def _save_intermediate_results(self, results, output_file):
        temp_file = f"temp_{output_file}"
        try:
            with open(temp_file, 'w') as f:
                json.dump(results, f, indent=2)
        except Exception as e:
            print(f"Warning: Could not save intermediate results: {str(e)}")
    
    def _save_final_results(self, results, output_file):
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
