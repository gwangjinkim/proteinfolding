'''
brew install pymol

conda create --name proteinfolding
conda activate proteinfolding

conda install -c conda-forge python ipython
pip install reqests biopython pytest

git init



'''

import os
import requests
from pathlib import Path
from Bio.PDB import PDBParser
import numpy as np
import matplotlib.pyplot as plt
from plum import dispatch
import Bio.PDB


class ProteinAnalyzer:
    def __init__(self, uniprot_accession, out_dir="/Users/josephus/Downloads/pdb_results/"):
        """
        Initialize the ProteinAnalyzer with a UniProt accession number.

        Parameters:
            uniprot_accession (str): The UniProt accession number of the protein.
        """
        self.uniprot_accession = uniprot_accession
        self.pdb_ids = []
        self.out_dir = Path(out_dir)
        self.pdb_files = []
        self.structure = None
        self.residue_b_factors = []
        self.residue_ids = []
        self.high_flex_residues = []
        self.threshold = None
        self.pdb_info = {}

        # ensure existence of output directory
        os.makedirs(self.out_dir, exist_ok=True)

    def fetch_pdb_ids(self):
        """
        Fetch PDB IDs associated with the UniProt accession number.
        """
        # Uniprot's API
        response = requests.get(f'https://rest.uniprot.org/uniprotkb/{self.uniprot_accession}.json')

        if response.status_code == 200:
            data = response.json()
            self.pdb_ids = [ref['id'] for ref in data['uniProtKBCrossReferences'] if ref['database'] == 'PDB']
            print(f"PDB IDs associated with UniProt accession {self.uniprot_accession}: {self.pdb_ids}")
        else:
            print(f"Failed to retrieve data for {self.uniprot_accession}")
            self.pdb_ids = []
        
    def download_pdb_files(self):
        """
        Download PDB files for the fetched PDB IDs.
        """
        if not self.pdb_ids:
            print("No PDB IDs to download.")
            return
         
        for pdb_id in self.pdb_ids:
            pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            pdb_response = requests.get(pdb_url)
            
            if pdb_response.status_code == 200:
                # Save the PDB file
                filename = Path(f"{self.out_dir}/{pdb_id}.pdb")
                with open(filename, 'wb') as file:
                    file.write(pdb_response.content)
                self.pdb_files.append(filename)
                print(f"Downloaded PDB file for {pdb_id}")
            else:
                print(f"Failed to download PDB file for {pdb_id}")

    def load_pdb_files_from_directory(self):
        """
        Load all PDB files from self.out_dir and update self.pdb_files and self.pdb_ids.
        """
        # Find all files in self.out_dir ending with '.pdb'
        pdb_files = list(self.out_dir.glob('*.pdb'))
        if pdb_files:
            self.pdb_files = pdb_files
            # Extract PDB IDs from filenames
            self.pdb_ids = [pdb_file.stem.upper() for pdb_file in pdb_files]
            print(f"Found {len(pdb_files)} PDB files in {self.out_dir}")
            print(f"PDB IDs: {self.pdb_ids}")
        else:
            print(f"No PDB files found in {self.out_dir}")
            
    def parse_pdb_file(self, pdb_file):
        """
        Parse a PDB file and store its structure.

        Parameters:
            pdb_file (Path): The path to the PDB file.
        """
        parser = PDBParser(QUIET=True)
        self.structure = parser.get_structure(pdb_file.stem, pdb_file)

    def analyze_flexibility(self):
        """
        Analyze the parsed PDB structure to identify regions of high flexibility based on B-factors.
        """
        if not self.structure:
            print("No structure to analyze.")
            return
        self.residue_b_factors = []
        self.residue_ids = []
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    # Skip hetero atoms and water molecules
                    if residue.id[0] != ' ':
                        continue
                    b_factors = [atom.get_bfactor() for atom in residue]
                    avg_b_factor = np.mean(b_factors)
                    self.residue_b_factors.append(avg_b_factor)
                    # Create a unique residue identifier
                    res_id = f"{chain.id}_{residue.id[1]}_{residue.get_resname()}"
                    self.residue_ids.append(res_id)
        # Identify residues with high B-factors
        mean_b = np.mean(self.residue_b_factors)
        std_b = np.std(self.residue_b_factors)
        self.threshold = mean_b + std_b
        self.high_flex_residues = [
            (res_id, b_factor) for res_id, b_factor in zip(self.residue_ids, self.residue_b_factors) if b_factor > self.threshold
        ]
        print("Residues with high flexibility:")
        for res_id, b_factor in self.high_flex_residues:
            print(f"{res_id}: B-factor = {b_factor:.2f}")

    def create_save_path(self, pdb_file):
        save_dir = f"{self.out_dir}/plots/"
        os.makedirs(save_dir, exist_ok=True)
        save_path = f"{save_dir}{pdb_file.stem}.png"
        return save_path
        
    def plot_b_factors(self, save_path=None):
        """
        Plot the average B-factors along the sequence and indicate the high flexibility threshold.
        """
        if not self.residue_b_factors:
            print("No B-factors to plot. Please run analyze_flexibility() first.")
            return
        plt.figure(figsize=(10, 4))
        plt.plot(self.residue_b_factors, label='Average B-factor')
        if self.threshold is not None:
            plt.axhline(y=self.threshold, color='r', linestyle='--', label='High flexibility threshold')
        plt.xlabel('Residue Index')
        plt.ylabel('Average B-factor')
        plt.title(f'B-factors along the sequence of {self.structure.id}')
        plt.legend()
        if save_path:
            plt.savefig(save_path)
            print(f"Plot saved to {save_path}")
            plt.close()
        else:
            plt.show()

    def parse_pdb_header(self, pdb_file):
        """
        Parse the PDB header to extract resolution and R-factors.

        Parameters:
            pdb_file (Path): The path to the PDB file.

        Returns:
            dict: A dictionary containing resolution, R-factor, and free R-factor.
        """
        header_info = {'resolution': None, 'r_value': None, 'free_r_value': None}
        with open(pdb_file, 'r') as file:
            for line in file:
                if line.startswith('REMARK   2 RESOLUTION.'):
                    match = re.search(r"([\d\.]+)\s+ANGSTROMS\.", line)
                    if match:
                        header_info['resolution'] = float(match.group(1))
                elif line.startswith('REMARK   3   R VALUE'):
                    if 'WORKING SET' in line:
                        match = re.search(r"R VALUE\s+\(WORKING SET\)\s+:\s+([\d\.]+)", line)
                        if match:
                            header_info['r_value'] = float(match.group(1))
                    elif 'FREE R VALUE' in line:
                        match = re.search(r"FREE R VALUE\s+:\s+([\d\.]+)", line)
                        if match:
                            header_info['free_r_value'] = float(match.group(1))
                elif line.startswith('END'):
                    break  # Stop parsing after header
        return header_info

    def get_longest_contiguous_chain_length(self, pdb_file):
        """
        Calculate the length of the longest contiguous chain in the PDB file.

        Parameters:
            pdb_file (Path): The path to the PDB file.

        Returns:
            tuple: Chain ID and length of the longest contiguous chain.
        """
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_file.stem, pdb_file)
        max_length = 0
        max_chain_id = None
        for model in structure:
            for chain in model:
                # Get sorted list of residue IDs
                residue_numbers = [residue.id[1] for residue in chain if residue.id[0] == ' ']
                if not residue_numbers:
                    continue
                residue_numbers.sort()
                # Identify contiguous segments
                contiguous_segments = []
                current_segment = [residue_numbers[0]]
                for i in range(1, len(residue_numbers)):
                    if residue_numbers[i] == residue_numbers[i - 1] + 1:
                        current_segment.append(residue_numbers[i])
                    else:
                        contiguous_segments.append(current_segment)
                        current_segment = [residue_numbers[i]]
                contiguous_segments.append(current_segment)
                # Find the longest contiguous segment
                for segment in contiguous_segments:
                    segment_length = len(segment)
                    if segment_length > max_length:
                        max_length = segment_length
                        max_chain_id = chain.id
        return max_chain_id, max_length

    def assess_pdb_files(self):
        """
        Assess all PDB files to find the one with the longest contiguous chain and best quality.
        """
        self.pdb_files = list(self.out_dir.glob('*.pdb'))
        if not self.pdb_files:
            print("No PDB files to assess.")
            return
        for pdb_file in self.pdb_files:
            pdb_id = pdb_file.stem.upper()
            header_info = self.parse_pdb_header(pdb_file)
            chain_id, chain_length = self.get_longest_contiguous_chain_length(pdb_file)
            self.pdb_info[pdb_id] = {
                'pdb_file': pdb_file,
                'resolution': header_info['resolution'],
                'r_value': header_info['r_value'],
                'free_r_value': header_info['free_r_value'],
                'chain_id': chain_id,
                'chain_length': chain_length
            }
        # Display the collected information
        print("\nPDB File Assessment:")
        for pdb_id, info in self.pdb_info.items():
            print(f"PDB ID: {pdb_id}")
            print(f"  Chain ID: {info['chain_id']}")
            print(f"  Chain Length: {info['chain_length']}")
            print(f"  Resolution: {info['resolution']} Å")
            print(f"  R-value: {info['r_value']}")
            print(f"  Free R-value: {info['free_r_value']}")

    def select_best_pdb(self):
        """
        Select the PDB file with the longest contiguous chain and best quality.

        Returns:
            str: The PDB ID of the best PDB file.
        """
        if not self.pdb_info:
            print("No PDB information available to select the best file.")
            return None
        # Sort PDBs based on chain length and resolution
        sorted_pdbs = sorted(
            self.pdb_info.items(),
            key=lambda item: (
                -item[1]['chain_length'],      # Longer chains first
                item[1]['resolution'] if item[1]['resolution'] is not None else float('inf'),  # Lower resolution values (better quality) first
                item[1]['free_r_value'] if item[1]['free_r_value'] is not None else float('inf')  # Lower free R-value
            )
        )
        best_pdb_id = sorted_pdbs[0][0]
        print(f"\nSelected best PDB: {best_pdb_id}")
        return best_pdb_id






import os
import requests

class AlphaFoldStructureFetcher:
    def __init__(self, uniprot_id, output_dir='structures'):
        self.uniprot_id = uniprot_id
        self.base_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
        self.output_dir = output_dir
        self.pdb_file = None

    def fetch_structure(self):
        """Fetches the AlphaFold structure (PDB format) for the given UniProt ID."""
        # Ensure the output directory exists
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        # Define the output PDB file path
        self.pdb_file = os.path.join(self.output_dir, f"{self.uniprot_id}.pdb")

        # Check if the PDB file already exists locally
        if os.path.exists(self.pdb_file):
            print(f"PDB file for {self.uniprot_id} already exists: {self.pdb_file}")
            return self.pdb_file

        # Fetch the PDB structure from AlphaFold
        response = requests.get(self.base_url)

        if response.status_code == 200:
            with open(self.pdb_file, 'wb') as f:
                f.write(response.content)
            print(f"Downloaded PDB file for {self.uniprot_id}: {self.pdb_file}")
        else:
            print(f"Failed to retrieve the PDB file for UniProt ID {self.uniprot_id}. Status code: {response.status_code}")
            self.pdb_file = None

        return self.pdb_file

class AlphaFoldPDBAnalyzer:
    def __init__(self, pdb_file):
        self.pdb_file = pdb_file
        self.structure = None
        self.pLDDT_scores = []
        self.rigid_regions = []
        self.flexible_regions = []
        self.load_structure()

    def load_structure(self):
        """Loads the PDB structure using BioPython."""
        parser = Bio.PDB.PDBParser(QUIET=True)
        self.structure = parser.get_structure("protein", self.pdb_file)
        self.extract_pLDDT_scores()

    def extract_pLDDT_scores(self):
        """Extracts the pLDDT scores from B-factors in the PDB file."""
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        self.pLDDT_scores.append(atom.bfactor)  # pLDDT is stored in B-factor field
                        break  # Only extract the B-factor from one atom per residue

    def classify_regions(self, rigid_threshold=90, flexible_threshold=70):
        """Classifies regions into rigid and flexible based on pLDDT scores."""
        for i, score in enumerate(self.pLDDT_scores):
            if score >= rigid_threshold:
                self.rigid_regions.append((i, score))
            elif score <= flexible_threshold:
                self.flexible_regions.append((i, score))

    def plot_pLDDT_distribution(self):
        """Plots the distribution of pLDDT scores."""
        plt.hist(self.pLDDT_scores, bins=50, color='skyblue', edgecolor='black')
        plt.title('pLDDT Score Distribution')
        plt.xlabel('pLDDT Score')
        plt.ylabel('Residue Count')
        plt.show()

    def get_rigid_regions(self):
        """Returns rigid regions based on pLDDT classification."""
        return self.rigid_regions

    def get_flexible_regions(self):
        """Returns flexible regions based on pLDDT classification."""
        return self.flexible_regions

class ProteinVisualizer:
    def __init__(self, pdb_file):
        self.pdb_file = pdb_file

    def visualize_by_pLDDT(self, output_file="colored_by_pLDDT.pse"):
        """Visualizes the PDB structure by coloring regions based on pLDDT scores."""
        try:
            import pymol
            from pymol import cmd

            pymol.finish_launching()
            cmd.load(self.pdb_file)

            # Set coloring by B-factor (which contains pLDDT scores)
            cmd.spectrum("b", "blue_white_red", minimum=0, maximum=100)

            # Save the visualization session
            cmd.save(output_file)
            print(f"Visualization saved to {output_file}")
        except ImportError:
            print("PyMOL is not installed or not available in the environment.")
        except Exception as e:
            print(f"Error during PyMOL visualization: {e}")

class ProteinVisualizer:
    def __init__(self, pdb_file):
        self.pdb_file = pdb_file

    def visualize_by_pLDDT(self, output_file="colored_by_pLDDT.pse"):
        """Visualizes the PDB structure by coloring regions based on pLDDT scores."""
        try:
            import pymol
            from pymol import cmd

            pymol.finish_launching(['pymol', '-cq'])  # Launch PyMOL in headless mode (no GUI)
            cmd.load(self.pdb_file)

            # Set coloring by B-factor (which contains pLDDT scores)
            cmd.spectrum("b", "blue_white_red", minimum=0, maximum=100)

            # Save the visualization session
            cmd.save(output_file)
            print(f"Visualization saved to {output_file}")
        except ImportError:
            print("PyMOL is not installed or not available in the environment.")
        except Exception as e:
            print(f"Error during PyMOL visualization: {e}")

'''  
# Create an instance of the ProteinAnalyzer class
analyzer = ProteinAnalyzer("P0DOX5")

# Step 1: Fetch PDB IDs associated with the UniProt accession
analyzer.fetch_pdb_ids()

# Step 2: Download the PDB files for the fetched PDB IDs
analyzer.download_pdb_files()

# Step 3: Assess all PDB files
analyzer.assess_pdb_files()

# Step 4: Select the best PDB file
best_pdb_id = analyzer.select_best_pdb()

# Step 5: Perform further analysis on the best PDB file
if best_pdb_id:
    best_pdb_file = analyzer.pdb_info[best_pdb_id]['pdb_file']
    analyzer.parse_pdb_file(best_pdb_file)
    analyzer.analyze_flexibility()
    # Save the B-factors plot
    plot_filename = analyzer.out_dir / f"{analyzer.structure.id}_b_factors.png"
    analyzer.plot_b_factors(save_path=plot_filename)

## /Users/josephus/Downloads/pdb_results/6B70_b_factors.png
'''

'''
uniprot_id = "P0DOX5"
fetcher = AlphaFoldStructureFetcher(uniprot_id)
pdb_file = fetcher.fetch_structure()

if pdb_file:
    print(f"Structure successfully retrieved: {pdb_file}")
else:
    print("Failed to retrieve the structure.")


analyzer = AlphaFoldPDBAnalyzer(pdb_file)

# Classify regions
analyzer.classify_regions()

# Plot plDDT distribution
analyzer.plot_pLDDT_distribution()

# Retrieve rigit and flexible regions
rigid_regions = analyzer.get_rigid_regions()
flexible_regions = analyzer.get_flexible_regions()

print(f"Rigid regions: {rigid_regions[:10]}") # first 10 rigit residues
print(f"Flexible regions: {flexible_regions[:10]}") # first 10 flexible residues

# Visualize structure colored by pLDDT
visualizer = ProteinVisualizer(pdb_file)
visualizer.visualize_by_pLDDT()

# open from ipython with:
# !pymol colored_by_pLDDT.pse

# to see the differences better, you can set:
# spectrum b, blue_white_red, minimum=85, maximum=100
'''


