import json
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
import sys
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
import functools

def process_single_compound(compound):
    """
    Process a single compound with RDKit filtering
    Returns tuple: (processed_compound_dict or None, skip_reason or None)
    """
    smiles = compound.get('SMILES', '')
    cid = compound.get('CID', '')
    
    # Parse SMILES with RDKit
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        return None, f"Could not parse SMILES for CID {cid}: {smiles}"
    
    # Check for single atom fragments
    if mol.GetNumAtoms() <= 1:
        return None, f"Single atom fragment CID {cid}"
    
    # Define unwanted elements (metals, metalloids, and specific non-metals)
    unwanted_elements = {
        # 'B',   # Boron. I would like to remove it, but two of the reactions
        # I'm interested in involve boronic acid. aryl_boronic_acid_robust and
        # boronic_acid_robust

        'Si',  # Silicon
        'Se',  # Selenium
        'Al',  # Aluminum
        'As',  # Arsenic
        # Common metals
        'Li', 'Na', 'K', 'Rb', 'Cs', 'Fr',  # Alkali metals
        'Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra',  # Alkaline earth metals
        'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',  # Transition metals (3d)
        'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',  # Transition metals (4d)
        'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',   # Transition metals (5d)
        'Ga', 'In', 'Tl', 'Ge', 'Sn', 'Pb', 'Sb', 'Bi', 'Po',       # Post-transition metals
        'Te',  # Tellurium
        # Lanthanides and Actinides
        'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
        'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr'
    }
    
    # Check for unwanted elements
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in unwanted_elements:
            return None, f"Contains unwanted element {atom.GetSymbol()} in CID {cid}"
    
    # Neutralize each ionizable group independently
    try:
        # Custom neutralization function that neutralizes each atom independently
        def neutralize_atoms(mol):
            """
            Neutralize ionizable atoms independently rather than maintaining overall charge balance
            """
            pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
            at_matches = mol.GetSubstructMatches(pattern)
            at_matches_list = [y[0] for y in at_matches]
            
            if len(at_matches_list) > 0:
                for at_idx in at_matches_list:
                    atom = mol.GetAtomWithIdx(at_idx)
                    chg = atom.GetFormalCharge()
                    hcount = atom.GetTotalNumHs()
                    
                    if chg < 0:
                        # Add hydrogen to neutralize negative charge
                        atom.SetFormalCharge(0)
                        atom.SetNumExplicitHs(hcount + abs(chg))
                    elif chg > 0:
                        # Remove hydrogen to neutralize positive charge (if possible)
                        if hcount >= chg:
                            atom.SetFormalCharge(0)
                            atom.SetNumExplicitHs(hcount - chg)
                        else:
                            # Can't neutralize this positive charge by removing H
                            pass
            
            return mol
        
        # Apply custom neutralization
        mol = neutralize_atoms(mol)
        
        # Sanitize the molecule
        Chem.SanitizeMol(mol)
        
    except Exception as e:
        return None, f"Could not neutralize CID {cid}: {str(e)}"
    
    # Calculate accurate molecular weight
    accurate_mw = Descriptors.MolWt(mol)
    
    # Exclude if MW > 250
    if accurate_mw > 250:
        return None, f"MW {accurate_mw:.2f} > 250 for CID {cid}"
    
    # Canonicalize SMILES
    canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
    
    # Return processed compound
    processed_compound = {
        "CID": cid,
        "SMILES": canonical_smiles,
        "MW": round(accurate_mw, 2)
    }
    
    return processed_compound, None

def process_molecules(input_file, output_file, num_processes=None):
    """
    Process molecules from JSON file using RDKit with multiprocessing:
    - Calculate accurate molecular weight
    - Exclude compounds with MW > 250
    - Exclude fragments containing unwanted elements (metals, metalloids, etc.)
    - Exclude single atom fragments
    - Neutralize molecules
    - Canonicalize SMILES
    """
    
    # Load input JSON
    try:
        with open(input_file, 'r') as f:
            data = json.load(f)
    except FileNotFoundError:
        print(f"Error: File {input_file} not found")
        return
    except json.JSONDecodeError:
        print(f"Error: Invalid JSON in {input_file}")
        return
    
    compounds = data.get('compounds', [])
    
    # Determine number of processes
    if num_processes is None:
        num_processes = min(cpu_count(), len(compounds))
    
    print(f"Processing {len(compounds)} compounds using {num_processes} processes...")
    
    # Process compounds with multiprocessing
    processed_compounds = []
    skip_reasons = []
    
    with Pool(processes=num_processes) as pool:
        # Use tqdm with multiprocessing
        results = list(tqdm(
            pool.imap(process_single_compound, compounds),
            total=len(compounds),
            desc="Processing compounds",
            unit="compound"
        ))
    
    # Collect results and remove duplicates
    seen_smiles = set()
    duplicate_count = 0
    
    for result, skip_reason in results:
        if result is not None:
            canonical_smiles = result['SMILES']
            if canonical_smiles not in seen_smiles:
                seen_smiles.add(canonical_smiles)
                processed_compounds.append(result)
            else:
                duplicate_count += 1
        else:
            skip_reasons.append(skip_reason)
    
    # Create output data
    output_data = {
        "compounds": processed_compounds
    }
    
    # Write output JSON
    try:
        with open(output_file, 'w') as f:
            json.dump(output_data, f, indent=2)
        
        print(f"\nProcessing complete!")
        print(f"Original compounds: {len(compounds)}")
        print(f"Processed compounds: {len(processed_compounds)}")
        print(f"Skipped compounds: {len(skip_reasons)}")
        print(f"Duplicate SMILES removed: {duplicate_count}")
        print(f"Unique compounds retained: {len(processed_compounds)}")
        
        # Show breakdown of skip reasons
        if skip_reasons:
            print("\nSkip reasons breakdown:")
            skip_types = {}
            for reason in skip_reasons:
                if "Could not parse" in reason:
                    skip_types["Invalid SMILES"] = skip_types.get("Invalid SMILES", 0) + 1
                elif "Single atom" in reason:
                    skip_types["Single atom"] = skip_types.get("Single atom", 0) + 1
                elif "unwanted element" in reason:
                    skip_types["Contains unwanted elements"] = skip_types.get("Contains unwanted elements", 0) + 1
                elif "neutralize" in reason:
                    skip_types["Neutralization failed"] = skip_types.get("Neutralization failed", 0) + 1
                elif "MW" in reason:
                    skip_types["MW > 250"] = skip_types.get("MW > 250", 0) + 1
            
            for skip_type, count in skip_types.items():
                print(f"  {skip_type}: {count}")
        
        print(f"Output saved to: {output_file}")
        
    except Exception as e:
        print(f"Error writing output file: {e}")

if __name__ == "__main__":
    input_file = "1.small_molecules_1st_pass.json"
    output_file = "2.small_molecules_2nd_pass.json"
    
    # Use all available CPU cores, or specify a number
    process_molecules(input_file, output_file, num_processes=None)