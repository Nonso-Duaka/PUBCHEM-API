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

def fetch_patents_by_cid(cid):
    # PubChem API URL to get patent data for the given CID
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/patents/JSON"
    
    # Send GET request to fetch patent data
    response = requests.get(url)
    
    if response.status_code == 200:
        data = response.json()
        
        # Check if patents are present in the response
        if "Patent" in data:
            patents = data["Patent"]
            patent_info = []
            for patent in patents:
                # Extract patent number and title (if available)
                patent_number = patent.get("PatentNumber", "No patent number")
                patent_title = patent.get("Title", "No title available")
                
                patent_info.append({
                    "Patent Number": patent_number,
                    "Title": patent_title
                })
            return patent_info
        else:
            return {"Error": "No patents found for this CID."}
    else:
        return {"Error": f"Failed to fetch patent data. Status code: {response.status_code}"}


def fetch_molecule_details(cid):
   
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/description/JSON"
    
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()

        # Extract the InformationList from the response
        information_list = data.get("InformationList", {}).get("Information", [])

        # Initialize the results list
        results = []

        if information_list:
            for info in information_list:
                # Create a dictionary for each information entry
                details = {
                    "CID": info.get("CID", "N/A"),
                    "Description": info.get("Description", "No Description available"),
                    "Source": info.get("DescriptionSourceName", "N/A"),
                    "URL": info.get("DescriptionURL", "N/A")
                }


                results.append(details)
        else:
            # Add a default entry only if no information is available
            results.append({
                "Description": "No Description available",
                "Source": "N/A",
                "URL": "N/A"
            })

        return results

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
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularWeight,MolecularFormula,IsomericSMILES,CanonicalSMILES,IUPACName,Charge/JSON"
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
                "Isomeric SMILES": Property.get("IsomericSMILES", "N/A"),
                "IUPAC Name": Property.get("IUPACName", "N/A")
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
        
        # Navigate to the desired section
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
        
        # Extract pictograms
        pictograms = [
            markup["URL"]
            for info in ghs_section["Information"]
            if info["Name"] == "Pictogram(s)"
            for markup in info["Value"]["StringWithMarkup"][0]["Markup"]
        ]
        
        # Extract signal word
        signal = next(
            info["Value"]["StringWithMarkup"][0]["String"]
            for info in ghs_section["Information"]
            if info["Name"] == "Signal"
        )
        
        # Extract hazard statements
        hazard_statements = [
            statement["String"]
            for info in ghs_section["Information"]
            if info["Name"] == "GHS Hazard Statements"
            for statement in info["Value"]["StringWithMarkup"]
        ]
        
        # Combine results into a structured dictionary
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
        if cids:
            return {"Similar Compounds": cids}
        else:
            return {"error": "No similar compounds found. Please adjust the threshold or check the SMILES input."}
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
    
    # Default threshold if user input is empty or invalid
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
        #"Patent Numbers": patent_numbers,
        #"Suppliers": suppliers,
        "Similar Compounds": similar_compounds,
        "Substructure Compounds":substructure_compounds,
        "Superstructure Compounds":superstructure_compounds
    }
    
    save_file(combined_data)

if __name__ == "__main__":
    main()
