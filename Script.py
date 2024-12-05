import requests
import json
import os
from requests.exceptions import Timeout, RequestException

def fetch_CID(SMILES):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{SMILES}/cids/JSON"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()
        
        cids = data.get("IdentifierList", {}).get("CID", [])
        if cids:
            return str(cids[0])  
        else:
            return "Error: Invalid SMILES string or no CID found for the given SMILES."
    except Timeout:
        return "Network Error: Request timed out. Please try again."
    except RequestException as e:
        return f"Network Error: {e}"
    except (json.JSONDecodeError, KeyError) as e:
        return f"Data Error: Unexpected API response format - {e}"
    

def fetch_molecule_details(cid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/description/JSON"
    
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()
        information_list = data.get("InformationList", {}).get("Information", [])
        results = []

        for info in information_list:
            details = {
                "CID": info.get("CID", "N/A")
            }
            description = info.get("Description")
            source = info.get("DescriptionSourceName")
            url = info.get("DescriptionURL")
            
            if description:
                details["Description"] = description
            if source:
                details["Source"] = source
            if url:
                details["URL"] = url

            results.append(details)

        return results if results else [{"error": "No information available for the given CID."}]

    except Timeout:
        return {"error": "Request timed out while fetching data."}
    except RequestException as e:
        return {"error": f"Network issue occurred: {e}"}
    except (KeyError, IndexError, ValueError) as e:
        return {"error": f"Unexpected data format: {e}"}


def fetch_patents(cid):
   
    try:
        # Request for patent-related properties (PatentCount, PatentFamilyCount)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/PatentCount,PatentFamilyCount/JSON"
        response = requests.get(url)
        
        if response.status_code == 200:
            data = response.json()
            properties = data.get("PropertyTable", {}).get("Properties", [])
            
            # If properties are found, extract patent information
            if properties:
                result = properties[0]  # As we request multiple properties, take the first result
                patent_count = result.get("PatentCount", "Not available")
                patent_family_count = result.get("PatentFamilyCount", "Not available")
                
                return {
                    "PatentCount": patent_count,
                    "PatentFamilyCount": patent_family_count
                }
            else:
                return {"CID": cid, "Error": "No patent data available"}
        else:
            return {"CID": cid, "Error": f"Failed to fetch data. Status code: {response.status_code}"}
    
    except requests.exceptions.RequestException as e:
        return {"CID": cid, "Error": f"An error occurred: {str(e)}"}


"""def fetch_patents_by_cid(cid):
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/classification/hnid/4501233/patents/JSON"
    response = requests.get(url)
    
    # Check if the request was successful
    if response.status_code != 200:
        return {"Error": f"Failed to fetch data. Status code: {response.status_code}"}
    
    # Attempt to parse the JSON response
    try:
        data = response.json()
    except ValueError:
        return {"Error": "Invalid JSON response"}
    
    # Extract patent data from the parsed JSON
    patents = data.get("Patents", [])
    if not patents:
        return {"Error": "No patent data available"}
    
    # Filter patents by the given CID
    filtered_patents = [patent for patent in patents if patent.get("CID") == str(cid)]
    if not filtered_patents:
        return {"CID": cid, "Error": "No patent data found for this CID"}
    
    return {"CID": cid, "Patents": filtered_patents}


def fetch_suppliers(cid):
    # Update the URL to fetch information from the 'substance' source table
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/sourcetable/substance/JSON?cid={cid}"

    try:
        # Send the GET request with a timeout to avoid hanging
        response = requests.get(url, timeout=10)
        response.raise_for_status()  # Raise exception for 4xx/5xx responses
        data = response.json()
        
        # Extract supplier information from the response
        sources = data.get('InformationList', {}).get('Information', [])
        suppliers_info = []
        
        for source in sources:
            if "SourceName" in source:
                suppliers_info.append({
                    "Supplier": source.get("SourceName", "N/A"),
                    "URL": source.get("SourceURL", "N/A")
                })

        # Return supplier info or an error message if no suppliers found
        if suppliers_info:
            return suppliers_info
        else:
            return {"error": "No supplier information available"}

    except requests.exceptions.RequestException as e:
        # Handle network errors, including timeouts and invalid responses
        return {"error": f"Failed to retrieve supplier information due to network issue: {e}"}"""


def fetch_synonyms(cid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()
        
        descriptions = data.get("InformationList", {}).get("Information", [{}])
        if descriptions:
            synonyms = descriptions[0]
            result = {"Synonyms": synonyms.get("Synonym", "No Synonyms available")}
            return result
        else:
            return {"error": "No synonyms available"}
    except Timeout:
        return {"error": "Failed to retrieve synonyms due to request timeout."}
    except RequestException as e:
        return {"error": f"Failed to retrieve synonyms due to network issue: {e}"}
    except (json.JSONDecodeError, KeyError) as e:
        return {"error": f"Data Error: Unexpected API response format - {e}"}



def fetch_compounds_properties(cid):
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/"
        "MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES,InChI,InChIKey,"
        "IUPACName,Title,XLogP,ExactMass,MonoisotopicMass,TPSA,Complexity,Charge,"
        "HBondDonorCount,HBondAcceptorCount,RotatableBondCount,HeavyAtomCount,"
        "IsotopeAtomCount,AtomStereoCount,DefinedAtomStereoCount,UndefinedAtomStereoCount,"
        "BondStereoCount,DefinedBondStereoCount,UndefinedBondStereoCount,CovalentUnitCount,"
        "PatentCount,PatentFamilyCount,LiteratureCount,Volume3D,XStericQuadrupole3D,"
        "YStericQuadrupole3D,ZStericQuadrupole3D,FeatureCount3D,FeatureAcceptorCount3D,"
        "FeatureDonorCount3D,FeatureAnionCount3D,FeatureCationCount3D,FeatureRingCount3D,"
        "FeatureHydrophobeCount3D,ConformerModelRMSD3D,EffectiveRotorCount3D,"
        "ConformerCount3D,Fingerprint2D/JSON"
    )
    
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()
        
        properties = data.get('PropertyTable', {}).get('Properties', [])
        if properties:
            Property = properties[0]
            result = {
                "Molecular Formula": Property.get("MolecularFormula", "N/A"),
                "Molecular Weight": Property.get("MolecularWeight", "N/A"),
                "Canonical SMILES": Property.get("CanonicalSMILES", "N/A"),
                "Isomeric SMILES": Property.get("IsomericSMILES", "N/A"),
                "InChI": Property.get("InChI", "N/A"),
                "InChIKey": Property.get("InChIKey", "N/A"),
                "IUPAC Name": Property.get("IUPACName", "N/A"),
                "Title": Property.get("Title", "N/A"),
                "XLogP": Property.get("XLogP", "N/A"),
                "Exact Mass": Property.get("ExactMass", "N/A"),
                "Monoisotopic Mass": Property.get("MonoisotopicMass", "N/A"),
                "Topological Polar Surface Area (TPSA)": Property.get("TPSA", "N/A"),
                "Complexity": Property.get("Complexity", "N/A"),
                "Charge": Property.get("Charge", "N/A"),
                "Hydrogen Bond Donor Count": Property.get("HBondDonorCount", "N/A"),
                "Hydrogen Bond Acceptor Count": Property.get("HBondAcceptorCount", "N/A"),
                "Rotatable Bond Count": Property.get("RotatableBondCount", "N/A"),
                "Heavy Atom Count": Property.get("HeavyAtomCount", "N/A"),
                "Isotope Atom Count": Property.get("IsotopeAtomCount", "N/A"),
                "Atom Stereo Count": Property.get("AtomStereoCount", "N/A"),
                "Defined Atom Stereo Count": Property.get("DefinedAtomStereoCount", "N/A"),
                "Undefined Atom Stereo Count": Property.get("UndefinedAtomStereoCount", "N/A"),
                "Bond Stereo Count": Property.get("BondStereoCount", "N/A"),
                "Defined Bond Stereo Count": Property.get("DefinedBondStereoCount", "N/A"),
                "Undefined Bond Stereo Count": Property.get("UndefinedBondStereoCount", "N/A"),
                "Covalent Unit Count": Property.get("CovalentUnitCount", "N/A"),
                "Patent Count": Property.get("PatentCount", "N/A"),
                "Patent Family Count": Property.get("PatentFamilyCount", "N/A"),
                "Literature Count": Property.get("LiteratureCount", "N/A"),
                "Volume (3D)": Property.get("Volume3D", "N/A"),
                "X Steric Quadrupole (3D)": Property.get("XStericQuadrupole3D", "N/A"),
                "Y Steric Quadrupole (3D)": Property.get("YStericQuadrupole3D", "N/A"),
                "Z Steric Quadrupole (3D)": Property.get("ZStericQuadrupole3D", "N/A"),
                "Feature Count (3D)": Property.get("FeatureCount3D", "N/A"),
                "Feature Acceptor Count (3D)": Property.get("FeatureAcceptorCount3D", "N/A"),
                "Feature Donor Count (3D)": Property.get("FeatureDonorCount3D", "N/A"),
                "Feature Anion Count (3D)": Property.get("FeatureAnionCount3D", "N/A"),
                "Feature Cation Count (3D)": Property.get("FeatureCationCount3D", "N/A"),
                "Feature Ring Count (3D)": Property.get("FeatureRingCount3D", "N/A"),
                "Feature Hydrophobe Count (3D)": Property.get("FeatureHydrophobeCount3D", "N/A"),
                "Conformer Model RMSD (3D)": Property.get("ConformerModelRMSD3D", "N/A"),
                "Effective Rotor Count (3D)": Property.get("EffectiveRotorCount3D", "N/A"),
                "Conformer Count (3D)": Property.get("ConformerCount3D", "N/A"),
                "Fingerprint (2D)": Property.get("Fingerprint2D", "N/A")
            }
            return result
        else:
            return {"error": "No properties found for the provided CID"}
    except Timeout:
        return {"error": "Failed to retrieve properties due to request timeout."}
    except RequestException as e:
        return {"error": f"Failed to retrieve properties due to network issue: {e}"}
    except (json.JSONDecodeError, KeyError) as e:
        return {"error": f"Data Error: Unexpected API response format - {e}"}


    

def fetch_hazard_information(cid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON/?response_type=display&heading=GHS%20Classification"
    
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()
     
        safety_section = next(
            section for section in data["Record"]["Section"]
            if section["TOCHeading"] == "Safety and Hazards"
        )
        hazards_section = next(
            section for section in safety_section["Section"]
            if section["TOCHeading"] == "Hazards Identification"
        )
        ghs_section = next(
            section for section in hazards_section["Section"]
            if section["TOCHeading"] == "GHS Classification"
        )
      
        pictograms = [
            markup["URL"]
            for info in ghs_section["Information"]
            if info["Name"] == "Pictogram(s)"
            for markup in info["Value"]["StringWithMarkup"][0]["Markup"]
        ]
      
        signal = next(
            info["Value"]["StringWithMarkup"][0]["String"]
            for info in ghs_section["Information"]
            if info["Name"] == "Signal"
        )
       
        hazard_statements = [
            statement["String"]
            for info in ghs_section["Information"]
            if info["Name"] == "GHS Hazard Statements"
            for statement in info["Value"]["StringWithMarkup"]
        ]
        
        result = {
            "Pictograms": pictograms,
            "Signal": signal,
            "Hazard Statements": hazard_statements
        }
        return result

    except Timeout:
        return {"error": "Failed to retrieve hazard information due to request timeout."}
    except RequestException as e:
        return {"error": f"Failed to retrieve hazard information due to network issue: {e}"}
    except (KeyError, StopIteration, IndexError, json.JSONDecodeError) as e:
        return {"error": f"Data Error: Unexpected API response format - {e}"}



def fetch_active_assays(cid):
  
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/assaysummary/JSON"
    
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()

       
        table = data.get("Table", {})
        columns = table.get("Columns", {}).get("Column", [])
        rows = table.get("Row", [])

     
        try:
            outcome_index = columns.index("Activity Outcome")
        except ValueError:
            return {"error": "Activity Outcome column not found in data."}

    
        active_assays = []
        for row in rows:
            cells = row.get("Cell", [])
            if len(cells) > outcome_index and cells[outcome_index] == "Active":
                assay = {columns[i]: cells[i] for i in range(len(cells))}
                active_assays.append(assay)

        if not active_assays:
            return {"error": "No active assays found for the provided CID."}

        return {"ActiveAssays": active_assays}

    except requests.exceptions.Timeout:
        return {"error": "Request timed out while fetching assay summary data."}
    except requests.exceptions.RequestException as e:
        return {"error": f"Network issue occurred: {e}"}
    except (ValueError, KeyError) as e:
        return {"error": f"Unexpected API response format: {e}"}



def save_file(result, filename="compound_data.json"):
    full_path = os.path.abspath(filename)
    with open(filename, "w") as jf:
        json.dump(result, jf, indent=4)
    print(f"Data saved successfully to {filename} at {full_path}.")


def fetch_similar_compounds(smiles, threshold=95, max_records=100):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/{smiles}/cids/JSON?Threshold={threshold}&MaxRecords={max_records}" 
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()
        cids = data.get("IdentifierList", {}).get("CID", [])
        if not cids:
            return {"error": "No similar compounds found. Please adjust the threshold or check the SMILES input."}

        compound_data = []
        for cid in cids:
            prop_url = (
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/"
                "CanonicalSMILES/JSON"
            )
            try:
                prop_response = requests.get(prop_url, timeout=5)
                prop_response.raise_for_status()
                prop_data = prop_response.json()
                properties = prop_data.get('PropertyTable', {}).get('Properties', [])
                if properties:
                    smiles_string = properties[0].get("CanonicalSMILES", "N/A")
                    compound_data.append({"CID": cid, "SMILES": smiles_string})
            except Exception as e:
                continue
        if compound_data:
            return {"Similar Compounds": compound_data}
        else:
            return {"error": "Failed to retrieve SMILES strings for the similar compounds."}
    except Timeout:
        return {"error": "Request timed out while fetching similar compounds."}
    except RequestException as e:
        return {"error": f"Network issue occurred: {e}"}
    except (KeyError, json.JSONDecodeError) as e:
        return {"error": f"Unexpected data format: {e}"}


def fetch_substructure_compounds(smiles, max_records=100):
    
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsubstructure/smiles/{smiles}/cids/JSON?MatchIsotopes=true&MaxRecords={max_records}"
    
    
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()

        cids = data.get("IdentifierList", {}).get("CID", [])
        if cids:
            return {"Substructure Compounds": cids}
        else:
            return {"error": "No substructure compounds found. Please check the SMILES input."}
    except Timeout:
        return {"error": "Request timed out while fetching substructure compounds."}
    except RequestException as e:
        return {"error": f"Network issue occurred: {e}"}
    except (KeyError, json.JSONDecodeError) as e:
        return {"error": f"Unexpected data format: {e}"}


def fetch_superstructure_compounds(smiles, max_records=100):
    
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsuperstructure/smiles/{smiles}/cids/JSON?MatchIsotopes=true&MaxRecords={max_records}"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()

        cids = data.get("IdentifierList", {}).get("CID", [])
        if cids:
            return {"Superstructure Compounds": cids}
        else:
            return {"error": "No superstructure compounds found. Please check the SMILES input."}
    except Timeout:
        return {"error": "Request timed out while fetching superstructure compounds."}
    except RequestException as e:
        return {"error": f"Network issue occurred: {e}"}
    except (KeyError, json.JSONDecodeError) as e:
        return {"error": f"Unexpected data format: {e}"}

def main():
    SMILES = input("Enter the SMILES of your choice: ").strip()
    threshold = input("Enter the Tanimoto coefficient threshold (default is 95): ").strip()

    if not threshold.isdigit():
        threshold = 95
    else:
        threshold = int(threshold)

    cid = fetch_CID(SMILES)
    
    if "Error" in cid:
        print(cid)
        return
    
    description = fetch_molecule_details(cid)
    property = fetch_compounds_properties(cid)
    synonyms = fetch_synonyms(cid)
    safety_info = fetch_hazard_information(cid)
    patents_count = fetch_patents(cid)
    bioassay_data = fetch_active_assays(cid)
    #patent_numbers =  fetch_patents_by_cid(cid)
    #suppliers = fetch_suppliers(cid)
    similar_compounds = fetch_similar_compounds(SMILES, threshold)  
    substructure_compounds = fetch_substructure_compounds(SMILES)
    superstructure_compounds = fetch_superstructure_compounds(SMILES)
   
    combined_data = {
        "Description": description,
        "Properties": property,
        "Synonyms": synonyms,
        "Safety Information": safety_info,
        "Patents": patents_count,
        "Bioassay Data": bioassay_data,
        #"Patent Numbers": patent_numbers,
        #"Suppliers": suppliers,
        "Similar Compounds": similar_compounds,
        "Substructure Compounds":substructure_compounds,
        "Superstructure Compounds":superstructure_compounds
    }
    
    save_file(combined_data)

if __name__ == "__main__":
    main()
