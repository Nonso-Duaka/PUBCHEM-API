#!/usr/bin/env python3
"""
PubChem Small Molecule Extractor - Fast Version (No RDKit)

Downloads PubChem SDF files one by one from FTP and extracts CID + SMILES
for compounds with molecular weight < 150 Da. Uses simple atom counting
for desalting instead of RDKit for much faster performance.

Requirements:
    pip install requests tqdm
"""

import gzip
import requests
import re
import json
import csv
from io import BytesIO
from tqdm import tqdm
import os
from typing import List, Dict, Tuple, Set, Any, Optional


class PubChemSmallMoleculeExtractor:
    def __init__(self, max_mw: float = 150.0) -> None:
        """
        Initialize the extractor.
        
        Args:
            max_mw (float): Maximum molecular weight threshold
        """
        self.max_mw = max_mw
        self.base_ftp_url = "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/"
        self.results = []
        self.seen_smiles = set()  # Track unique SMILES to avoid duplicates
        
    def count_atoms_in_smiles_part(self, smiles_part: str) -> int:
        """
        Count atoms in a SMILES string part using regex.
        This method finds atoms in brackets (e.g., [Na+]), two-letter
        elements (e.g., Cl), and single-letter elements (e.g., C, c).
        Args:
            smiles_part (str): Part of SMILES string.
        Returns:
            int: The count of atoms.
        """
        # Regex to find atoms: matches bracketed expressions, two-letter
        # elements (e.g., Cl, Br), or single-letter elements (e.g., C, c).
        atom_pattern = r'\[[^\]]+\]|[A-Z][a-z]?|[a-z]'
        return len(re.findall(atom_pattern, smiles_part))
    
    def desalt_smiles(self, smiles):
        """
        Simple desalting: if SMILES contains dots, keep the part with most atoms.
        
        Args:
            smiles (str): Input SMILES string (potentially with salts)
            
        Returns:
            str: Desalted SMILES string (largest component)
        """
        # Only desalt if there's a dot (indicating multiple components)
        if '.' not in smiles:
            return smiles.strip()
        
        # Split on dots to get components
        components = smiles.split('.')
        
        # Find component with most atoms
        max_atoms = 0
        largest_component = components[0]
        
        for component in components:
            atom_count = self.count_atoms_in_smiles_part(component)
            if atom_count > max_atoms:
                max_atoms = atom_count
                largest_component = component
        
        return largest_component.strip()
    
    def estimate_molecular_weight(self, smiles: str) -> float:
        """
        Calculate molecular weight from SMILES using simple atom counting.
        This version correctly handles aromatic (lowercase) atoms and does
        not include an estimate for hydrogen atoms.
        Args:
            smiles (str): SMILES string
            
        Returns:
            float: Molecular weight in Da
            
        Raises:
            ValueError: If unknown atoms are found or calculation fails
        """
        # Simple atomic weights for common elements
        atomic_weights = {
            'C': 12.01, 'N': 14.01, 'O': 16.00, 'S': 32.06,
            'P': 30.97, 'F': 19.00, 'Cl': 35.45, 'Br': 79.90,
            'I': 126.90, 'H': 1.008, 'B': 10.81, 'Si': 28.09,
            'Na': 22.99, 'K': 39.10, 'Ca': 40.08, 'Mg': 24.31,
            'Al': 26.98, 'Fe': 55.85, 'Zn': 65.38, 'Cu': 63.55,
            'Mn': 54.94, 'Co': 58.93, 'Ni': 58.69, 'Cr': 51.996,
            'Li': 6.94, 'Be': 9.01, 'Se': 78.96, 'As': 74.92,
            'Sn': 118.71, 'Pb': 207.2, 'Ag': 107.87, 'Au': 196.97,
            'Pt': 195.08, 'Pd': 106.42, 'Rh': 102.91, 'Ru': 101.07,
            'Os': 190.23, 'Ir': 192.22, 'Re': 186.21, 'W': 183.84,
            'Mo': 95.96, 'Tc': 98.91, 'V': 50.94, 'Ti': 47.87,
            'Sc': 44.96, 'Y': 88.91, 'La': 138.91, 'Ce': 140.12,
            'Pr': 140.91, 'Nd': 144.24, 'Pm': 145.0, 'Sm': 150.36,
            'Eu': 151.96, 'Gd': 157.25, 'Tb': 158.93, 'Dy': 162.50,
            'Ho': 164.93, 'Er': 167.26, 'Tm': 168.93, 'Yb': 173.05,
            'Lu': 174.97, 'Hf': 178.49, 'Ta': 180.95, 'Ra': 226.03,
            'Ac': 227.03, 'Th': 232.04, 'Pa': 231.04, 'U': 238.03,
            'Np': 237.05, 'Pu': 244.06, 'Am': 243.06, 'Cm': 247.07,
            'Bk': 247.07, 'Cf': 251.08, 'Es': 252.08, 'Fm': 257.10,
            'Md': 258.10, 'No': 259.10, 'Lr': 262.11, 'Rf': 267.12,
            'Db': 268.13, 'Sg': 271.13, 'Bh': 270.13, 'Hs': 277.15,
            'Mt': 276.15, 'Ds': 281.16, 'Rg': 280.16, 'Cn': 285.17
        }
        total_weight = 0.0
        unknown_atoms = set()
        # Find atoms inside square brackets, e.g., [Na+], [O-], extracting the element symbol.
        bracket_atoms = re.findall(r'\[([A-Z][a-z]?)', smiles)
        # Remove square bracket sections to process the rest of the string.
        simple_smiles = re.sub(r'\[[^\]]*\]', '', smiles)
        # Find regular atoms (e.g., C, Cl) and aromatic atoms (e.g., c, n).
        other_atoms = re.findall(r'[A-Z][a-z]?|[a-z]', simple_smiles)
        all_atom_symbols = bracket_atoms + other_atoms
        for symbol in all_atom_symbols:
            # Normalize symbol to element (e.g., 'c' becomes 'C') for lookup.
            element = symbol.capitalize()
            if element in atomic_weights:
                total_weight += atomic_weights[element]
            else:
                unknown_atoms.add(symbol)
        # Check for unknown atoms and raise error if found
        if unknown_atoms:
            raise ValueError(f"Unknown atoms found in SMILES '{smiles}': {', '.join(sorted(unknown_atoms))}")
        # Per user request, the inaccurate hydrogen estimation has been removed.
        return total_weight
        
    def parse_sdf_record(self, record: str, debug_count: int = 0) -> Optional[Dict[str, Any]]:
        """
        Parse a single SDF record to extract CID, SMILES, and molecular weight.
        Automatically desalts SMILES by keeping the largest component.
        
        Args:
            record (str): Single SDF record
            debug_count (int): Record number for debugging
            
        Returns:
            dict or None: Compound data if MW < threshold, None otherwise
        """
        try:
            # Extract CID
            cid_match = re.search(r'> <PUBCHEM_COMPOUND_CID>\n(\d+)', record)
            if not cid_match:
                return None
            cid = int(cid_match.group(1))
            
            # Extract canonical SMILES
            smiles_match = re.search(r'> <PUBCHEM_OPENEYE_CAN_SMILES>\n([^\n]+)', record)
            if not smiles_match:
                # Try alternative SMILES field
                smiles_match = re.search(r'> <PUBCHEM_OPENEYE_ISO_SMILES>\n([^\n]+)', record)
            
            if not smiles_match:
                # This was an `assert False` which would crash the script.
                # Replaced with a non-crashing warning to improve robustness.
                assert False, f"Warning: SMILES not found for CID {cid}. Skipping record."
                return None
            
            raw_smiles = smiles_match.group(1).strip()
            
            # Desalt: keep the largest component (only if dots present)
            smiles = self.desalt_smiles(raw_smiles)
            if not smiles:  # If desalting resulted in empty string
                return None
            
            # Estimate molecular weight from the desalted SMILES
            try:
                mw = self.estimate_molecular_weight(smiles)
            except Exception:
                # If MW estimation fails, skip this compound
                return None
            # Round the MW *before* checking against the threshold. This ensures
            # that a compound with MW 149.9998, which rounds to 150.000,
            # is correctly excluded.
            rounded_mw = round(mw, 3)
            # Check molecular weight threshold on the rounded value
            if rounded_mw >= self.max_mw:
                return None
            
            return {
                'CID': cid,
                'SMILES': smiles,
                'MW': rounded_mw
            }
            
        except (ValueError, AttributeError) as e:
            return None
    
    def process_sdf_stream(self, sdf_stream: Any) -> Tuple[List[Dict[str, Any]], int, int]:
        """
        Process an SDF stream and extract small molecules.
        
        Args:
            sdf_stream: File-like object containing SDF data
            
        Returns:
            tuple: (list of new compound dictionaries, count of new compounds, count of duplicates)
        """
        new_compounds = []
        duplicate_count = 0
        current_record = ""
        
        for line in sdf_stream:
            if isinstance(line, bytes):
                # Use 'replace' to handle decoding errors without crashing or ignoring them silently.
                line = line.decode('utf-8', errors='replace')
            current_record += line
            
            # End of record marker
            if line.strip() == "$$$$":
                compound = self.parse_sdf_record(current_record)
                if compound:
                    smiles = compound['SMILES']
                    if smiles not in self.seen_smiles:
                        # New unique compound
                        self.seen_smiles.add(smiles)
                        new_compounds.append(compound)
                    else:
                        # Duplicate SMILES
                        duplicate_count += 1
                current_record = ""
        return new_compounds, len(new_compounds), duplicate_count
    
    def download_and_process_file(self, filename: str) -> Tuple[List[Dict[str, Any]], int, int]:
        """
        Download and process a single SDF file.
        
        Args:
            filename (str): Name of the SDF file to download
            
        Returns:
            tuple: (list of new compounds, count of new compounds, count of duplicates)
        """
        url = f"{self.base_ftp_url}{filename}"
        
        try:
            response = requests.get(url, stream=True, timeout=300)
            response.raise_for_status()
            
            # Get total file size for progress bar
            total_size = int(response.headers.get('content-length', 0))
            
            # Download with progress bar
            downloaded_data = BytesIO()
            with tqdm(total=total_size, unit='B', unit_scale=True, desc=f"Downloading {filename}") as pbar:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        downloaded_data.write(chunk)
                        pbar.update(len(chunk))
            
            # Reset stream position
            downloaded_data.seek(0)
            
            # Decompress gzipped content
            with gzip.GzipFile(fileobj=downloaded_data) as gz_file:
                new_compounds, new_count, duplicate_count = self.process_sdf_stream(gz_file)
            
            print(f"Found {new_count} NEW compounds with MW < {self.max_mw} in {filename}")
            if duplicate_count > 0:
                print(f"Skipped {duplicate_count} duplicate SMILES in {filename}")
            
            return new_compounds, new_count, duplicate_count
            
        except Exception as e:
            print(f"Error processing {filename}: {e}")
            return [], 0, 0
    
    def get_file_list(self):
        """
        Get list of SDF files from PubChem FTP directory.
        
        Returns:
            list: List of SDF filenames
        """
        try:
            print("Fetching file list from PubChem FTP...")
            response = requests.get(self.base_ftp_url, timeout=30)
            response.raise_for_status()
            
            # Print first 1000 characters of response for debugging
            print("First 1000 characters of FTP directory listing:")
            print(response.text[:1000])
            print("=" * 50)
            
            # Updated regex pattern for current PubChem file naming convention
            # Looking for patterns like: Compound_000000001_000500000.sdf.gz
            filenames = re.findall(r'Compound_\d{9}_\d{9}\.sdf\.gz', response.text)
            
            if not filenames:
                # Try alternative patterns
                alt_patterns = [
                    r'Compound_\d+_\d+\.sdf\.gz',  # More flexible digit matching
                    r'href="(Compound_[^"]+\.sdf\.gz)"',  # Look for href links
                ]
                
                for pattern in alt_patterns:
                    filenames = re.findall(pattern, response.text)
                    if filenames:
                        print(f"Found files using pattern: {pattern}")
                        break
            
            # Remove duplicates while preserving order
            seen = set()
            unique_filenames = []
            for filename in filenames:
                if filename not in seen:
                    seen.add(filename)
                    unique_filenames.append(filename)
            filenames = unique_filenames
            
            if filenames:
                print(f"Found {len(filenames)} SDF files")
                print("First 5 files:", filenames[:5])
            else:
                print("No SDF files found! Checking for any .gz files...")
                all_gz = re.findall(r'href="([^"]+\.gz)"', response.text)
                print(f"Found {len(all_gz)} .gz files:", all_gz[:10])
                
            return sorted(filenames)
            
        except Exception as e:
            print(f"Error getting file list: {e}")
            print(f"Response status code: {getattr(response, 'status_code', 'N/A')}")
            return []
    
    def extract_all_small_molecules(self, output_file: str = "small_molecules.json",
                                        resume: bool = True, max_files: Optional[int] = None) -> List[Dict[str, Any]]:
        """
        Download and process all PubChem SDF files to extract small molecules.
        
        Args:
            output_file (str): Output filename for results
            resume (bool): Whether to resume from existing progress
            max_files (int, optional): Maximum number of files to process (for testing)
        Returns:
            List[Dict[str, Any]]: A list of new compounds found during this run.
        """
        # Load existing progress if resuming
        processed_files = set()
        if resume and os.path.exists(output_file):
            try:
                with open(output_file, 'r') as f:
                    existing_data = json.load(f)
                    self.results = existing_data.get('compounds', [])
                    processed_files = set(existing_data.get('processed_files', []))
                    
                    # Rebuild the seen_smiles set from existing data
                    self.seen_smiles = {compound['SMILES'] for compound in self.results}
                    
                print(f"Resuming: loaded {len(self.results)} existing compounds")
                print(f"Already processed {len(processed_files)} files")
                print(f"Tracking {len(self.seen_smiles)} unique SMILES")
            except Exception as e:
                print(f"Error loading existing data: {e}")
                # Reset if loading failed
                self.results = []
                self.seen_smiles = set()
                processed_files = set()
        
        # Get list of files to process
        all_files = self.get_file_list()
        if not all_files:
            print("No files found! Check the file list method.")
            return []
        # Filter out already processed files
        files_to_process = [f for f in all_files if f not in processed_files]
        
        if max_files:
            files_to_process = files_to_process[:max_files]
        
        print(f"Found {len(all_files)} total files")
        print(f"Processing {len(files_to_process)} files")
        # Track totals and compounds found specifically in this run
        compounds_this_run = []
        total_new_compounds = 0
        total_duplicates = 0
        
        # Process files with progress bar
        for i, filename in enumerate(files_to_process):
            print(f"\nProcessing file {i+1}/{len(files_to_process)}: {filename}")
            new_compounds, new_count, duplicate_count = self.download_and_process_file(filename)
            
            self.results.extend(new_compounds)
            compounds_this_run.extend(new_compounds)
            processed_files.add(filename)
            
            total_new_compounds += new_count
            total_duplicates += duplicate_count
            
            # Save progress periodically
            if len(processed_files) % 10 == 0:
                self._save_progress(output_file, processed_files)
                print(f"Progress saved. Total unique compounds so far: {len(self.results)}")
        
        # Final save
        self._save_progress(output_file, processed_files)
        
        print(f"\n" + "="*60)
        print(f"EXTRACTION COMPLETE!")
        print(f"Total NEW compounds found in this run: {total_new_compounds}")
        print(f"Total duplicates skipped this run: {total_duplicates}")
        print(f"Total unique compounds in cumulative database: {len(self.results)}")
        print(f"Cumulative results saved to {output_file}")
        print(f"="*60)
        return compounds_this_run
    
    def _save_progress(self, output_file, processed_files):
        """Save current progress to file."""
        try:
            data = {
                'compounds': self.results,
                'processed_files': list(processed_files),
                'total_compounds': len(self.results),
                'total_unique_smiles': len(self.seen_smiles),
                'max_mw': self.max_mw
            }
            
            with open(output_file, 'w') as f:
                json.dump(data, f, indent=2)
                
            # Also save as CSV for easy viewing
            csv_file = output_file.replace('.json', '.csv')
            with open(csv_file, 'w', newline='') as f:
                if self.results:
                    writer = csv.DictWriter(f, fieldnames=['CID', 'SMILES', 'MW'])
                    writer.writeheader()
                    writer.writerows(self.results)
                    
        except Exception as e:
            print(f"Error saving progress: {e}")
    
    def export_smiles_file(self, compounds: List[Dict[str, Any]], output_file: str = "small_molecules.smi") -> None:
        """
        Export results as a simple SMILES file.
        
        Args:
            compounds (List[Dict[str, Any]]): The list of compounds to export.
            output_file (str): Output SMILES filename
        """
        try:
            with open(output_file, 'w') as f:
                for compound in compounds:
                    f.write(f"{compound['SMILES']} {compound['CID']}\n")
            print(f"SMILES file with {len(compounds)} entries saved to {output_file}")
        except Exception as e:
            print(f"Error saving SMILES file: {e}")
    def analyze_results(self, compounds: List[Dict[str, Any]]) -> None:
        """
        Print analysis of the provided list of compounds.
        Args:
            compounds (List[Dict[str, Any]]): The list of compounds to analyze.
        """
        if not compounds:
            print("No compounds to analyze for this run.")
            return
        
        print(f"\n" + "="*50)
        print(f"ANALYSIS OF NEWLY FOUND COMPOUNDS")
        print(f"="*50)
        print(f"Total unique compounds found in this run: {len(compounds)}")
        # Molecular weight distribution
        mws = [compound['MW'] for compound in compounds]
        print(f"Molecular weight range: {min(mws):.2f} - {max(mws):.2f} Da")
        print(f"Average molecular weight: {sum(mws)/len(mws):.2f} Da")
        # Show some example compounds
        print(f"\nExample compounds from this run:")
        for i, compound in enumerate(compounds[:5]):
            print(f"  {i+1}. CID {compound['CID']}: {compound['SMILES']} (MW: {compound['MW']})")
        if len(compounds) > 5:
            print(f"  ... and {len(compounds) - 5} more compounds")


def main() -> None:
    """Main function to run the extraction."""
    print("PubChem Small Molecule Extractor - Fast Version (No RDKit)")
    print("=" * 60)
    
    # Initialize extractor
    extractor = PubChemSmallMoleculeExtractor(max_mw=150.0)
    
    # Test desalting functionality
    print("\nTesting desalting functionality...")
    test_cases = [
        "CCO",  # ethanol (no salt)
        "C(=O)(C(=O)[O-])N.[Na+]",  # sodium salt example
        "CC(C)CC(C(=O)[O-])N.[Na+]",  # another salt
        "CC(=O)O.CCN",  # mixture
        "C1=CC=CC=C1.[Cl-]",  # benzene with chloride
        "c1ccccc1Cl" # aromatic example with two-letter element
    ]
    
    for smiles in test_cases:
        desalted = extractor.desalt_smiles(smiles)
        if desalted:
            try:
                mw = extractor.estimate_molecular_weight(desalted)
                atom_count = extractor.count_atoms_in_smiles_part(desalted)
                print(f"  {smiles}")
                print(f" -> {desalted} (MW: {mw:.2f}, atoms: {atom_count})")
            except Exception as e:
                print(f" -> {desalted} (MW calculation failed: {e})")
        else:
            print(f"  {smiles} -> FAILED")
    
    # Test with just getting the file list first
    print("\nTesting file list retrieval...")
    files = extractor.get_file_list()
    
    if not files:
        print("ERROR: No files found. Please check the FTP directory structure.")
        return
    
    print(f"SUCCESS: Found {len(files)} files")
    
    # For testing, limit the number of files
    print("\nStarting extraction with first 2 files for testing...")
    newly_found_compounds = extractor.extract_all_small_molecules(max_files=2)
    # Analyze the results from THIS RUN
    extractor.analyze_results(newly_found_compounds)
    # Export a SMILES file containing only the compounds from THIS RUN
    extractor.export_smiles_file(
        newly_found_compounds,
        output_file="newly_found_molecules.smi"
    )
    # The full, cumulative database is always stored in small_molecules.json/.csv
    # To run a full extraction, you can remove the `max_files` argument.
    # Example for a full run (will take a long time):
    # print("\nStarting full extraction...")
    # all_new_compounds = extractor.extract_all_small_molecules()
    # extractor.analyze_results(all_new_compounds)
    # extractor.export_smiles_file(all_new_compounds, "all_new_compounds.smi")

if __name__ == "__main__":
    main()