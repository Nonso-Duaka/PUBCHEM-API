import json
from rdkit import Chem
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import sys
from collections import defaultdict

def match_compound_to_patterns(args):
    """
    Match a single compound against all functional group patterns
    Returns tuple: (compound, list_of_matching_pattern_names)
    """
    compound, patterns = args
    
    smiles = compound.get('SMILES', '')
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        return compound, []
    
    matching_patterns = []
    
    for pattern_name, smarts_pattern in patterns.items():
        try:
            # Create pattern molecule from SMARTS
            pattern_mol = Chem.MolFromSmarts(smarts_pattern)
            
            if pattern_mol is None:
                print(f"Warning: Could not parse SMARTS pattern for {pattern_name}: {smarts_pattern}")
                continue
            
            # Check if pattern matches the molecule
            if mol.HasSubstructMatch(pattern_mol):
                matching_patterns.append(pattern_name)
                
        except Exception as e:
            print(f"Error matching pattern {pattern_name}: {e}")
            continue
    
    return compound, matching_patterns

def process_functional_groups(molecules_file, patterns_file, output_file, num_processes=None):
    """
    Process molecules and group them by functional group patterns
    """
    
    # Load molecules
    try:
        with open(molecules_file, 'r') as f:
            molecules_data = json.load(f)
        compounds = molecules_data.get('compounds', [])
        print(f"Loaded {len(compounds)} compounds from {molecules_file}")
    except FileNotFoundError:
        print(f"Error: File {molecules_file} not found")
        return
    except json.JSONDecodeError:
        print(f"Error: Invalid JSON in {molecules_file}")
        return
    
    # Load functional group patterns
    try:
        with open(patterns_file, 'r') as f:
            patterns = json.load(f)
        print(f"Loaded {len(patterns)} functional group patterns from {patterns_file}")
    except FileNotFoundError:
        print(f"Error: File {patterns_file} not found")
        return
    except json.JSONDecodeError:
        print(f"Error: Invalid JSON in {patterns_file}")
        return
    
    # Determine number of processes
    if num_processes is None:
        num_processes = min(cpu_count(), len(compounds))
    
    print(f"Processing with {num_processes} processes...")
    
    # Prepare arguments for multiprocessing
    # Each compound needs access to all patterns
    args_list = [(compound, patterns) for compound in compounds]
    
    # Process compounds with multiprocessing
    with Pool(processes=num_processes) as pool:
        results = list(tqdm(
            pool.imap(match_compound_to_patterns, args_list),
            total=len(compounds),
            desc="Matching functional groups",
            unit="compound"
        ))
    
    # Organize results by functional group
    functional_group_matches = defaultdict(list)
    
    print("Organizing results by functional groups...")
    for compound, matching_patterns in tqdm(results, desc="Organizing results"):
        for pattern_name in matching_patterns:
            functional_group_matches[pattern_name].append(compound)
    
    # Initialize all pattern names (even if no matches)
    output_data = {}
    for pattern_name in patterns.keys():
        output_data[pattern_name] = functional_group_matches.get(pattern_name, [])
    
    # Write output JSON
    try:
        with open(output_file, 'w') as f:
            json.dump(output_data, f, indent=2)
        
        print(f"\nProcessing complete!")
        print(f"Total compounds processed: {len(compounds)}")
        
        # Show statistics
        print("\nFunctional group match statistics:")
        total_matches = 0
        non_empty_groups = 0
        
        for pattern_name, matched_compounds in output_data.items():
            count = len(matched_compounds)
            if count > 0:
                non_empty_groups += 1
                total_matches += count
                print(f"  {pattern_name}: {count} compounds")
        
        print(f"\nSummary:")
        print(f"  Functional groups with matches: {non_empty_groups}/{len(patterns)}")
        print(f"  Total pattern matches: {total_matches}")
        print(f"  Average matches per compound: {total_matches/len(compounds):.2f}")
        print(f"Output saved to: {output_file}")
        
    except Exception as e:
        print(f"Error writing output file: {e}")

def main():
    molecules_file = "2.small_molecules_2nd_pass.json"
    patterns_file = "functional_groups.json"
    output_file = "functional_group_matches.json"
    
    # Use all available CPU cores, or specify a number
    process_functional_groups(molecules_file, patterns_file, output_file, num_processes=None)

if __name__ == "__main__":
    main()