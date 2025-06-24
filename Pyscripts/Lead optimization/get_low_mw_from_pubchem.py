#!/usr/bin/env python3
"""
PubChem Small Molecule Extractor

Downloads PubChem SDF files one by one from FTP and extracts CID + SMILES
for compounds with molecular weight < 150 Da. Does not save SDF files permanently.
"""

import gzip
import requests
import re
import json
import csv
from io import BytesIO
from tqdm import tqdm
import os

class PubChemSmallMoleculeExtractor:
    def __init__(self, max_mw=150.0):
        """
        Initialize the extractor.
        
        Args:
            max_mw (float): Maximum molecular weight threshold
        """
        self.max_mw = max_mw
        self.base_ftp_url = "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/"
        self.results = []
        self.seen_smiles = set()  # Track unique SMILES to avoid duplicates
        
    def parse_sdf_record(self, record, debug_count=0):
        """
        Parse a single SDF record to extract CID, SMILES, and molecular weight.
        
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
                # if debug_count <= 3:  # Debug first few records
                #     print(f"  Debug record {debug_count}: No CID found")
                return None
            cid = int(cid_match.group(1))
            
            # Extract molecular weight
            mw_match = re.search(r'> <PUBCHEM_MOLECULAR_WEIGHT>\n([\d.]+)', record)
            if not mw_match:
                # if debug_count <= 3:
                #     print(f"  Debug record {debug_count} (CID {cid}): No MW found")
                return None
            mw = float(mw_match.group(1))
            
            # Debug first few records regardless of MW
            # if debug_count <= 3:
            #     print(f"  Debug record {debug_count} (CID {cid}): MW = {mw}")
            
            # Check molecular weight threshold
            if mw >= self.max_mw:
                return None
            
            # Extract canonical SMILES
            smiles_match = re.search(r'> <PUBCHEM_OPENEYE_CAN_SMILES>\n([^\n]+)', record)
            if not smiles_match:
                # Try alternative SMILES field
                smiles_match = re.search(r'> <PUBCHEM_OPENEYE_ISO_SMILES>\n([^\n]+)', record)
            
            if not smiles_match:
                # if debug_count <= 3:
                #     print(f"  Debug record {debug_count} (CID {cid}, MW {mw}): No SMILES found")
                return None
            
            smiles = smiles_match.group(1).strip()
            
            # if debug_count <= 3:
            #     print(f"  Debug record {debug_count} (CID {cid}, MW {mw}): Found SMILES = {smiles[:50]}...")
            
            return {
                'CID': cid,
                'SMILES': smiles,
                'MW': mw
            }
            
        except (ValueError, AttributeError) as e:
            # if debug_count <= 3:
            #     print(f"  Debug record {debug_count}: Parse error = {e}")
            return None
    
    def process_sdf_stream(self, sdf_stream):
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
                line = line.decode('utf-8', errors='ignore')
            
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
    
    def download_and_process_file(self, filename):
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
    
    def extract_all_small_molecules(self, output_file="small_molecules.json", 
                                  resume=True, max_files=None):
        """
        Download and process all PubChem SDF files to extract small molecules.
        
        Args:
            output_file (str): Output filename for results
            resume (bool): Whether to resume from existing progress
            max_files (int, optional): Maximum number of files to process (for testing)
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
            return
        
        # Filter out already processed files
        files_to_process = [f for f in all_files if f not in processed_files]
        
        if max_files:
            files_to_process = files_to_process[:max_files]
        
        print(f"Found {len(all_files)} total files")
        print(f"Processing {len(files_to_process)} files")
        
        # Track totals across all files
        total_new_compounds = 0
        total_duplicates = 0
        
        # Process files with progress bar
        for i, filename in enumerate(files_to_process):
            print(f"\nProcessing file {i+1}/{len(files_to_process)}: {filename}")
            new_compounds, new_count, duplicate_count = self.download_and_process_file(filename)
            
            self.results.extend(new_compounds)
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
        print(f"Total NEW compounds with MW < {self.max_mw}: {total_new_compounds}")
        print(f"Total duplicates skipped: {total_duplicates}")
        print(f"Total unique compounds in database: {len(self.results)}")
        print(f"Results saved to {output_file}")
        print(f"="*60)
    
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
    
    def export_smiles_file(self, output_file="small_molecules.smi"):
        """
        Export results as a simple SMILES file.
        
        Args:
            output_file (str): Output SMILES filename
        """
        try:
            with open(output_file, 'w') as f:
                for compound in self.results:
                    f.write(f"{compound['SMILES']} {compound['CID']}\n")
            print(f"SMILES file saved to {output_file}")
        except Exception as e:
            print(f"Error saving SMILES file: {e}")


def main():
    """Main function to run the extraction."""
    print("PubChem Small Molecule Extractor")
    print("=" * 40)
    
    # Initialize extractor
    extractor = PubChemSmallMoleculeExtractor(max_mw=150.0)
    
    # Test with just getting the file list first
    print("Testing file list retrieval...")
    files = extractor.get_file_list()
    
    if not files:
        print("ERROR: No files found. Please check the FTP directory structure.")
        return
    
    print(f"SUCCESS: Found {len(files)} files")
    
    # For testing, limit the number of files
    print("\nStarting extraction with first 2 files for testing...")
    extractor.extract_all_small_molecules(max_files=2)
    
    # Uncomment for full extraction (will take several hours)
    # extractor.extract_all_small_molecules()
    
    # Export as SMILES file
    extractor.export_smiles_file()


if __name__ == "__main__":
    main()