import requests
import json
import os
from requests.exceptions import Timeout, RequestException

def fetch_CID(SMILES):
    """
    Fetches the PubChem Compound ID (CID) for a given SMILES (Simplified Molecular Input Line Entry System) string.
    This function queries the PubChem API to retrieve the corresponding CID for a given SMILES string.
    If a valid CID is found, it returns the first CID from the response. Otherwise, it returns an error message.
    Args:
        SMILES (str): The SMILES representation of a chemical compound.
    Returns:
        str: The CID of the compound if found, otherwise an error message.
    """
    # Construct the PubChem API request URL using the provided SMILES string
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{SMILES}/cids/JSON"
    try:
        # Send a GET request to the PubChem API with a timeout of 10 seconds
        response = requests.get(url, timeout=10)
        response.raise_for_status()  # Raise an HTTPError if the request returned an error status
        # Parse the response JSON data
        data = response.json()
        # Extract the list of CIDs from the response
        cids = data.get("IdentifierList", {}).get("CID", [])
        if cids:
            return str(cids[0])  # Return the first CID as a string
        else:
            return "Error: Invalid SMILES string or no CID found for the given SMILES."
    # Handle request timeout errors
    except Timeout:
        return "Network Error: Request timed out. Please try again."
    # Handle general network-related errors
    except RequestException as e:
        return f"Network Error: {e}"
    # Handle cases where the API response format is unexpected
    except (json.JSONDecodeError, KeyError) as e:
        return f"Data Error: Unexpected API response format - {e}"
 
def fetch_molecule_details(cid):
    """
    Fetches a compound's description and related metadata from PubChem using its CID.
    This function queries the PubChem API to retrieve a detailed description of a compound
    based on its Compound ID (CID). The retrieved details may include:
    - Description of the compound
    - Source of the description
    - URL for additional information
    Args:
        cid (str): The PubChem Compound ID.
    Returns:
        list: A list of dictionaries containing molecule details (CID, Description, Source, URL).
              If no information is found, it returns an error message.
    """
    # Construct the PubChem API request URL using the provided CID
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/description/JSON"
    try:
        # Send a GET request to the PubChem API with a timeout of 10 seconds
        response = requests.get(url, timeout=10)
        response.raise_for_status()  # Raise an exception for HTTP errors
        # Parse the JSON response
        data = response.json()
        # Extract the information list containing compound details
        information_list = data.get("InformationList", {}).get("Information", [])
        results = []
        # Iterate through the extracted information and structure the details
        for info in information_list:
            details = {
                "CID": info.get("CID", "N/A")  # Ensure CID is present in the output
            }
            description = info.get("Description")
            source = info.get("DescriptionSourceName")
            url = info.get("DescriptionURL")
            # Add available details to the dictionary
            if description:
                details["Description"] = description
            if source:
                details["Source"] = source
            if url:
                details["URL"] = url
            # Append structured details to the results list
            results.append(details)
        # Return results if any information is found; otherwise, return an error message
        return results if results else [{"error": "No information available for the given CID."}]
    # Handle request timeout errors
    except Timeout:
        return {"error": "Request timed out while fetching data."}
    # Handle general network-related exceptions
    except RequestException as e:
        return {"error": f"Network issue occurred: {e}"}
    # Handle unexpected response formats and missing keys
    except (KeyError, IndexError, ValueError) as e:
        return {"error": f"Unexpected data format: {e}"}

def fetch_patents(cid):
    """
    Fetches patent-related properties for a given compound from PubChem.
    This function queries the PubChem API to retrieve patent-related data for a given 
    Compound ID (CID), specifically:
    - The number of patents associated with the compound (PatentCount)
    - The number of patent families associated with the compound (PatentFamilyCount)
    Args:
        cid (str): The PubChem Compound ID.
    Returns:
        dict: A dictionary containing:
            - "PatentCount": The number of patents linked to the compound.
            - "PatentFamilyCount": The number of patent families.
            - If no data is found, an error message is returned.
    """
    try:
        # Construct the PubChem API request URL for patent-related properties
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/PatentCount,PatentFamilyCount/JSON"
        # Send a GET request to the PubChem API
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for HTTP errors
        # Parse the response JSON data
        data = response.json()
        # Extract the list of properties from the response
        properties = data.get("PropertyTable", {}).get("Properties", [])
        # If properties exist, extract patent details
        if properties:
            result = properties[0]  # Get the first result (since multiple properties were requested)
            patent_count = result.get("PatentCount", "Not available")
            patent_family_count = result.get("PatentFamilyCount", "Not available")
            return {
                "PatentCount": patent_count,
                "PatentFamilyCount": patent_family_count
            }
        else:
            return {"CID": cid, "Error": "No patent data available"}
    # Handle network-related errors
    except requests.exceptions.RequestException as e:
        return {"CID": cid, "Error": f"An error occurred: {str(e)}"}

def fetch_synonyms(cid):
    """
    Fetches synonyms for a given compound from PubChem.
    This function queries the PubChem API to retrieve a list of synonyms for a given 
    Compound ID (CID). Synonyms can include common names, trade names, and alternative
    identifiers used in chemical databases.
    Args:
        cid (str): The PubChem Compound ID.
    Returns:
        dict: A dictionary containing:
            - "Synonyms": A list of synonyms for the compound.
            - If no synonyms are found, an error message is returned.
    """
    # Construct the PubChem API request URL to fetch synonyms
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
    try:
        # Send a GET request to the PubChem API with a 10-second timeout
        response = requests.get(url, timeout=10)
        response.raise_for_status()  # Raise an exception for HTTP errors

        # Parse the response JSON data
        data = response.json()

        # Extract the list of synonyms from the response
        descriptions = data.get("InformationList", {}).get("Information", [{}])

        if descriptions:
            synonyms = descriptions[0]
            return {"Synonyms": synonyms.get("Synonym", "No Synonyms available")}
        else:
            return {"error": "No synonyms available"}
    # Handle request timeout errors
    except Timeout:
        return {"error": "Failed to retrieve synonyms due to request timeout."}
    # Handle general network-related errors
    except RequestException as e:
        return {"error": f"Failed to retrieve synonyms due to network issue: {e}"}
    # Handle unexpected response formats and missing keys
    except (json.JSONDecodeError, KeyError) as e:
        return {"error": f"Data Error: Unexpected API response format - {e}"}

def fetch_compounds_properties(cid):
    """
    Fetches a comprehensive set of chemical properties for a given compound from PubChem.
    This function queries the PubChem API to retrieve various physicochemical and 
    structural properties of a compound using its Compound ID (CID). The retrieved 
    properties include molecular weight, formula, SMILES notation, InChIKey, 
    hydrogen bond counts, 3D structural data, and more.
    Args:
        cid (str): The PubChem Compound ID.
    Returns:
        dict: A dictionary containing various chemical properties.
              If no properties are found, an error message is returned.
    """
    # Construct the PubChem API request URL with multiple chemical properties
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
        # Send a GET request to the PubChem API with a 10-second timeout
        response = requests.get(url, timeout=10)
        response.raise_for_status()  # Raise an exception for HTTP errors
        # Parse the response JSON data
        data = response.json()
        # Extract the list of chemical properties
        properties = data.get('PropertyTable', {}).get('Properties', [])
        # If properties exist, extract and format them
        if properties:
            prop = properties[0]
            result = {
                "Molecular Formula": prop.get("MolecularFormula", "N/A"),
                "Molecular Weight": prop.get("MolecularWeight", "N/A"),
                "Canonical SMILES": prop.get("CanonicalSMILES", "N/A"),
                "Isomeric SMILES": prop.get("IsomericSMILES", "N/A"),
                "InChI": prop.get("InChI", "N/A"),
                "InChIKey": prop.get("InChIKey", "N/A"),
                "IUPAC Name": prop.get("IUPACName", "N/A"),
                "Title": prop.get("Title", "N/A"),
                "XLogP": prop.get("XLogP", "N/A"),
                "Exact Mass": prop.get("ExactMass", "N/A"),
                "Monoisotopic Mass": prop.get("MonoisotopicMass", "N/A"),
                "Topological Polar Surface Area (TPSA)": prop.get("TPSA", "N/A"),
                "Complexity": prop.get("Complexity", "N/A"),
                "Charge": prop.get("Charge", "N/A"),
                "Hydrogen Bond Donor Count": prop.get("HBondDonorCount", "N/A"),
                "Hydrogen Bond Acceptor Count": prop.get("HBondAcceptorCount", "N/A"),
                "Rotatable Bond Count": prop.get("RotatableBondCount", "N/A"),
                "Heavy Atom Count": prop.get("HeavyAtomCount", "N/A"),
                "Isotope Atom Count": prop.get("IsotopeAtomCount", "N/A"),
                "Atom Stereo Count": prop.get("AtomStereoCount", "N/A"),
                "Defined Atom Stereo Count": prop.get("DefinedAtomStereoCount", "N/A"),
                "Undefined Atom Stereo Count": prop.get("UndefinedAtomStereoCount", "N/A"),
                "Bond Stereo Count": prop.get("BondStereoCount", "N/A"),
                "Defined Bond Stereo Count": prop.get("DefinedBondStereoCount", "N/A"),
                "Undefined Bond Stereo Count": prop.get("UndefinedBondStereoCount", "N/A"),
                "Covalent Unit Count": prop.get("CovalentUnitCount", "N/A"),
                "Patent Count": prop.get("PatentCount", "N/A"),
                "Patent Family Count": prop.get("PatentFamilyCount", "N/A"),
                "Literature Count": prop.get("LiteratureCount", "N/A"),
                "Volume (3D)": prop.get("Volume3D", "N/A"),
                "X Steric Quadrupole (3D)": prop.get("XStericQuadrupole3D", "N/A"),
                "Y Steric Quadrupole (3D)": prop.get("YStericQuadrupole3D", "N/A"),
                "Z Steric Quadrupole (3D)": prop.get("ZStericQuadrupole3D", "N/A"),
                "Feature Count (3D)": prop.get("FeatureCount3D", "N/A"),
                "Feature Acceptor Count (3D)": prop.get("FeatureAcceptorCount3D", "N/A"),
                "Feature Donor Count (3D)": prop.get("FeatureDonorCount3D", "N/A"),
                "Feature Anion Count (3D)": prop.get("FeatureAnionCount3D", "N/A"),
                "Feature Cation Count (3D)": prop.get("FeatureCationCount3D", "N/A"),
                "Feature Ring Count (3D)": prop.get("FeatureRingCount3D", "N/A"),
                "Feature Hydrophobe Count (3D)": prop.get("FeatureHydrophobeCount3D", "N/A"),
                "Conformer Model RMSD (3D)": prop.get("ConformerModelRMSD3D", "N/A"),
                "Effective Rotor Count (3D)": prop.get("EffectiveRotorCount3D", "N/A"),
                "Conformer Count (3D)": prop.get("ConformerCount3D", "N/A"),
                "Fingerprint (2D)": prop.get("Fingerprint2D", "N/A")
            }
            return result
        else:
            return {"error": "No properties found for the provided CID"}
    # Handle request timeout errors
    except Timeout:
        return {"error": "Failed to retrieve properties due to request timeout."}
    # Handle general network-related errors
    except RequestException as e:
        return {"error": f"Failed to retrieve properties due to network issue: {e}"}
    # Handle unexpected response formats and missing keys
    except (json.JSONDecodeError, KeyError) as e:
        return {"error": f"Data Error: Unexpected API response format - {e}"}

def fetch_hazard_information(cid):
    """
    Fetches the hazard classification and safety data for a given compound from PubChem.
    This function queries the PubChem API to retrieve Global Harmonized System (GHS) 
    hazard classification details for a given Compound ID (CID). The retrieved 
    information includes:
    - Pictograms (hazard symbols)
    - Signal words (e.g., "Danger", "Warning")
    - Hazard statements (e.g., "Causes skin irritation")
    Args:
        cid (str): The PubChem Compound ID.
    Returns:
        dict: A dictionary containing:
            - "Pictograms": List of URLs to hazard pictograms.
            - "Signal": Signal word indicating hazard severity.
            - "Hazard Statements": List of hazard descriptions.
        If no hazard data is found, an error message is returned.
    """
    # Construct the PubChem API request URL to fetch GHS classification data
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON/?response_type=display&heading=GHS%20Classification"
    try:
        # Send a GET request to the PubChem API with a 10-second timeout
        response = requests.get(url, timeout=10)
        response.raise_for_status()  # Raise an exception for HTTP errors
        # Parse the JSON response
        data = response.json()
        # Navigate the JSON structure to extract hazard-related sections
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
        # Extract hazard pictograms (URLs to images)
        pictograms = [
            markup["URL"]
            for info in ghs_section["Information"]
            if info["Name"] == "Pictogram(s)"
            for markup in info["Value"]["StringWithMarkup"][0]["Markup"]
        ]
        # Extract the signal word (e.g., "Danger" or "Warning")
        signal = next(
            info["Value"]["StringWithMarkup"][0]["String"]
            for info in ghs_section["Information"]
            if info["Name"] == "Signal"
        )
        # Extract hazard statements (e.g., "Causes skin irritation")
        hazard_statements = [
            statement["String"]
            for info in ghs_section["Information"]
            if info["Name"] == "GHS Hazard Statements"
            for statement in info["Value"]["StringWithMarkup"]
        ]
        # Structure the extracted data into a dictionary
        result = {
            "Pictograms": pictograms,
            "Signal": signal,
            "Hazard Statements": hazard_statements
        }
        return result

    # Handle request timeout errors
    except Timeout:
        return {"error": "Failed to retrieve hazard information due to request timeout."}

    # Handle general network-related errors
    except RequestException as e:
        return {"error": f"Failed to retrieve hazard information due to network issue: {e}"}

    # Handle unexpected response formats and missing keys
    except (KeyError, StopIteration, IndexError, json.JSONDecodeError) as e:
        return {"error": f"Data Error: Unexpected API response format - {e}"}

def fetch_active_assays(cid):
    """
    Fetches bioassay activity data for a given compound from PubChem.
    This function queries the PubChem API to retrieve bioassay summary data for a given 
    Compound ID (CID). It filters and returns only those assays where the compound 
    has shown active results.
    Args:
        cid (str): The PubChem Compound ID.
    Returns:
        dict: A dictionary containing:
            - "ActiveAssays": A list of dictionaries representing bioassays where 
              the compound was found to be active.
            - If no active assays are found, an error message is returned.
    """
    # Construct the PubChem API request URL for assay summary data
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/assaysummary/JSON"
    try:
        # Send a GET request to the PubChem API with a 10-second timeout
        response = requests.get(url, timeout=10)
        response.raise_for_status()  # Raise an exception for HTTP errors
        # Parse the response JSON data
        data = response.json()
        # Extract assay summary table data
        table = data.get("Table", {})
        columns = table.get("Columns", {}).get("Column", [])  # List of column names
        rows = table.get("Row", [])  # List of assay results
        # Identify the index of the "Activity Outcome" column
        try:
            outcome_index = columns.index("Activity Outcome")
        except ValueError:
            return {"error": "Activity Outcome column not found in data."}
        # Filter the rows where the "Activity Outcome" is "Active"
        active_assays = []
        for row in rows:
            cells = row.get("Cell", [])
            if len(cells) > outcome_index and cells[outcome_index] == "Active":
                # Convert row data into a dictionary with column names as keys
                assay = {columns[i]: cells[i] for i in range(len(cells))}
                active_assays.append(assay)
        # Return results if any active assays are found, otherwise return an error
        if not active_assays:
            return {"error": "No active assays found for the provided CID."}
        return {"ActiveAssays": active_assays}
    # Handle request timeout errors
    except requests.exceptions.Timeout:
        return {"error": "Request timed out while fetching assay summary data."}
    # Handle general network-related errors
    except requests.exceptions.RequestException as e:
        return {"error": f"Network issue occurred: {e}"}
    # Handle unexpected response formats and missing keys
    except (ValueError, KeyError) as e:
        return {"error": f"Unexpected API response format: {e}"}

def save_file(result, filename="compound_data.json"):
    """
    Saves the retrieved data to a JSON file.
    This function takes a dictionary of retrieved compound data and saves it 
    as a formatted JSON file in the specified location.
    Args:
        result (dict): The data to save.
        filename (str, optional): The name of the JSON file (default is "compound_data.json").
    Returns:
        None: The function does not return a value but prints a confirmation message.
    """
    # Get the absolute file path
    full_path = os.path.abspath(filename)
    # Open the file in write mode and save the JSON data
    with open(filename, "w") as jf:
        json.dump(result, jf, indent=4)  # Pretty-print JSON with indentation
    # Print a confirmation message with the file path
    print(f"Data saved successfully to {filename} at {full_path}.")

def fetch_similar_compounds(smiles, threshold=95, max_records=100):
    """
    Fetches similar compounds from PubChem based on a given SMILES string.
    This function makes a request to the PubChem API to retrieve compounds with a 2D similarity above
    the given threshold relative to the provided SMILES string. It returns the list of compound IDs 
    (CIDs) and their corresponding canonical SMILES strings. If no similar compounds are found or if 
    an error occurs during the request, an appropriate error message is returned.
    Parameters:
    smiles (str): A SMILES string representing the compound whose similar compounds are to be found.
    threshold (int, optional): The minimum similarity percentage required for compounds to be considered similar (default is 95%).
    max_records (int, optional): The maximum number of similar compounds to return (default is 100).
    Returns:
    dict: A dictionary containing a list of similar compound CIDs with their SMILES or an error message.
          Example:
          - Success: {"Similar Compounds": [{"CID": CID1, "SMILES": SMILES1}, {"CID": CID2, "SMILES": SMILES2}, ...]}
          - Failure: {"error": "Error message"}
    """
    # Construct the PubChem API URL with the provided SMILES, threshold, and max_records parameters
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/{smiles}/cids/JSON?Threshold={threshold}&MaxRecords={max_records}"
    try:
        # Send a GET request to the PubChem API with a 10-second timeout
        response = requests.get(url, timeout=10)
        # Raise an exception for non-2xx HTTP responses (e.g., 404 or 500)
        response.raise_for_status()
        # Parse the JSON response from the API
        data = response.json()
        # Extract the list of compound IDs (CIDs) from the response
        cids = data.get("IdentifierList", {}).get("CID", [])
        # If no CIDs are found, return an error message
        if not cids:
            return {"error": "No similar compounds found. Please adjust the threshold or check the SMILES input."}
        compound_data = []
        # For each CID, fetch its canonical SMILES string
        for cid in cids:
            prop_url = (
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/"
                "CanonicalSMILES/JSON"
            )
            try:
                # Send a GET request to fetch the SMILES string for the compound
                prop_response = requests.get(prop_url, timeout=5)
                prop_response.raise_for_status()
                # Parse the JSON response containing the SMILES
                prop_data = prop_response.json()
                # Extract the SMILES string if available
                properties = prop_data.get('PropertyTable', {}).get('Properties', [])
                if properties:
                    smiles_string = properties[0].get("CanonicalSMILES", "N/A")
                    compound_data.append({"CID": cid, "SMILES": smiles_string})
            except Exception as e:
                # If an error occurs while fetching properties for a CID, continue with the next one
                continue
        # Return the list of similar compounds if found
        if compound_data:
            return {"Similar Compounds": compound_data}
        else:
            return {"error": "Failed to retrieve SMILES strings for the similar compounds."}
    except Timeout:
        # Handle request timeout errors
        return {"error": "Request timed out while fetching similar compounds."}
    except RequestException as e:
        # Handle network-related issues (e.g., connection errors)
        return {"error": f"Network issue occurred: {e}"}
    except (KeyError, json.JSONDecodeError) as e:
        # Handle unexpected data formats (e.g., missing keys or invalid JSON)
        return {"error": f"Unexpected data format: {e}"}

def fetch_substructure_compounds(smiles, max_records=100):
    """
    Fetches substructure compounds from PubChem using a given SMILES string.
    This function makes a request to the PubChem API to retrieve compounds that contain a substructure 
    matching the provided SMILES string. It returns a list of compound IDs (CIDs) in the response. If no 
    compounds are found or if an error occurs during the request, an appropriate error message is returned.
    Parameters:
    smiles (str): A SMILES string representing the substructure whose matching compounds are to be found.
    max_records (int, optional): The maximum number of records to return (default is 100).
    Returns:
    dict: A dictionary containing a list of substructure compound CIDs or an error message.
          Example:
          - Success: {"Substructure Compounds": [CID1, CID2, ...]}
          - Failure: {"error": "Error message"}
    """
    # Construct the PubChem API URL with the provided SMILES and max_records parameters
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsubstructure/smiles/{smiles}/cids/JSON?MatchIsotopes=true&MaxRecords={max_records}"
    try:
        # Send a GET request to the PubChem API with a 10-second timeout
        response = requests.get(url, timeout=10)
        # Raise an exception for non-2xx HTTP responses (e.g., 404 or 500)
        response.raise_for_status()
        # Parse the JSON response from the API
        data = response.json()
        # Extract the list of compound IDs (CIDs) from the response
        cids = data.get("IdentifierList", {}).get("CID", [])
        # If CIDs are found, return them in a dictionary
        if cids:
            return {"Substructure Compounds": cids}
        else:
            # If no CIDs are found, return an error message
            return {"error": "No substructure compounds found. Please check the SMILES input."}
    except Timeout:
        # Handle request timeout errors
        return {"error": "Request timed out while fetching substructure compounds."}   
    except RequestException as e:
        # Handle network-related issues (e.g., connection errors)
        return {"error": f"Network issue occurred: {e}"}
    except (KeyError, json.JSONDecodeError) as e:
        # Handle unexpected data formats (e.g., missing keys or invalid JSON)
        return {"error": f"Unexpected data format: {e}"}

def fetch_superstructure_compounds(smiles, max_records=100):
    """
    Fetches superstructure compounds from PubChem using a given SMILES string.
    This function makes a request to the PubChem API to retrieve compounds that have a superstructure 
    matching the provided SMILES string. It returns a list of compound IDs (CIDs) in the response. If no 
    compounds are found or if an error occurs during the request, an appropriate error message is returned.
    Parameters:
    smiles (str): A SMILES string representing the compound whose superstructure compounds are to be found.
    max_records (int, optional): The maximum number of records to return (default is 100).
    Returns:
    dict: A dictionary containing a list of superstructure compound CIDs or an error message.
          Example:
          - Success: {"Superstructure Compounds": [CID1, CID2, ...]}
          - Failure: {"error": "Error message"}
    """
    # Construct the PubChem API URL with the provided SMILES and max_records parameters
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsuperstructure/smiles/{smiles}/cids/JSON?MatchIsotopes=true&MaxRecords={max_records}"
    try:
        # Send a GET request to the PubChem API with a 10-second timeout
        response = requests.get(url, timeout=10)      
        # Raise an exception for non-2xx HTTP responses (e.g., 404 or 500)
        response.raise_for_status()    
        # Parse the JSON response from the API
        data = response.json()
        # Extract the list of compound IDs (CIDs) from the response
        cids = data.get("IdentifierList", {}).get("CID", [])
        # If CIDs are found, return them in a dictionary
        if cids:
            return {"Superstructure Compounds": cids}
        else:
            # If no CIDs are found, return an error message
            return {"error": "No superstructure compounds found. Please check the SMILES input."}
    except Timeout:
        # Handle request timeout errors
        return {"error": "Request timed out while fetching superstructure compounds."}
    except RequestException as e:
        # Handle network-related issues (e.g., connection errors)
        return {"error": f"Network issue occurred: {e}"}
    except (KeyError, json.JSONDecodeError) as e:
        # Handle unexpected data formats (e.g., missing keys or invalid JSON)
        return {"error": f"Unexpected data format: {e}"}

def main():
    """
    Main function to fetch and display various details about a compound based on its SMILES string.
    This function allows the user to input a SMILES string, then fetches various details about the 
    compound associated with that SMILES, including:
    - Description
    - Properties
    - Synonyms
    - Safety Information
    - Patents
    - Bioassay Data
    - Similar Compounds
    - Substructure Compounds
    - Superstructure Compounds
    The function combines the data from these various fetch operations and saves it into a file.
    """
    # Prompt the user for a SMILES string and strip any leading or trailing spaces
    SMILES = input("Enter the SMILES of your choice: ").strip()
    # Prompt the user for the Tanimoto coefficient threshold (default is 95)
    threshold = input("Enter the Tanimoto coefficient threshold (default is 95): ").strip()
    # If the threshold is not a valid number, default to 95
    if not threshold.isdigit():
        threshold = 95
    else:
        threshold = int(threshold)
    # Fetch the CID for the provided SMILES string
    cid = fetch_CID(SMILES)
    # If there is an error in fetching the CID, print the error and return
    if "Error" in cid:
        print(cid)
        return
    # Fetch the other details using the CID
    description = fetch_molecule_details(cid)
    property = fetch_compounds_properties(cid)
    synonyms = fetch_synonyms(cid)
    safety_info = fetch_hazard_information(cid)
    patents_count = fetch_patents(cid)
    bioassay_data = fetch_active_assays(cid)
    # Fetch similar, substructure, and superstructure compounds using the SMILES string and threshold
    similar_compounds = fetch_similar_compounds(SMILES, threshold)  
    substructure_compounds = fetch_substructure_compounds(SMILES)
    superstructure_compounds = fetch_superstructure_compounds(SMILES)
   
    # Combine all the fetched data into a dictionary
    combined_data = {
        "Description": description,
        "Properties": property,
        "Synonyms": synonyms,
        "Safety Information": safety_info,
        "Patents": patents_count,
        "Bioassay Data": bioassay_data,
        "Similar Compounds": similar_compounds,
        "Substructure Compounds": substructure_compounds,
        "Superstructure Compounds": superstructure_compounds
    }
    # Save the combined data into a file (the `save_file` function must be defined elsewhere)
    save_file(combined_data)
if __name__ == "__main__":
    main()
