# ProteinFolding Package

`ProteinFolding` is a small Python package that allows you to fetch and analyze protein structures, particularly those predicted by AlphaFold. This package provides functionality to:
- Fetch protein structures from the AlphaFold database by UniProt ID.
- Analyze predicted local distance difference test (pLDDT) scores to identify flexible and rigid regions.
- Visualize proteins using PyMOL, with options for custom coloring based on pLDDT scores.

## Table of Contents

- [Installation](#installation)
- [Classes Overview](#classes-overview)
- [Usage](#usage)
  - [Fetching Protein Structures](#fetching-protein-structures)
  - [Analyzing Protein Flexibility and Rigidity](#analyzing-protein-flexibility-and-rigidity)
  - [Visualizing Proteins with PyMOL](#visualizing-proteins-with-pymol)
  - [Changing Coloring in PyMOL](#changing-coloring-in-pymol)
- [Examples](#examples)
- [Contributing](#contributing)
- [License](#license)

## Installation

You can install this package by directly installing it with `pip` from this github repository:
```bash
pip install -e git+https://github.com/gwangjinkim/proteinfolding.git#egg=proteinfolding
```
Ensure you have PyMOL installed to use the visualization features:
- Install PyMOL: [PyMOL installation guide](https://pymol.org/2/#download)
```
# in MacOS
brew install pymol
```

## Classes Overview

### `ProteinAnalyzer`
An interface for analyzing protein structures in general. It can be extended or customized to handle protein-specific analyses, including those beyond AlphaFold predictions.

### `AlphaFoldStructureFetcher`
This class fetches protein structures (in PDB format) from the AlphaFold database using a given UniProt ID.

#### Methods:
- `fetch_structure()`: Downloads the structure for a given UniProt ID from AlphaFold and saves it locally in PDB format.

### `AlphaFoldPDBAnalyzer`
Analyzes protein structures from AlphaFold, focusing on identifying flexible and rigid regions using pLDDT scores.

#### Methods:
- `load_structure()`: Loads the PDB structure using BioPython.
- `extract_pLDDT_scores()`: Extracts pLDDT scores from the PDB's B-factor field.
- `classify_regions()`: Classifies residues as either rigid or flexible based on pLDDT thresholds.
- `plot_pLDDT_distribution()`: Plots the distribution of pLDDT scores.
- `get_rigid_regions()`: Returns a list of rigid regions.
- `get_flexible_regions()`: Returns a list of flexible regions.

### `ProteinVisualizer`
Visualizes protein structures using PyMOL. This class can color proteins based on pLDDT scores.

#### Methods:
- `visualize_by_pLDDT(output_file)`: Colors the protein structure based on pLDDT scores and saves the PyMOL session file.





## ProteinAnalyzer Class

The `ProteinAnalyzer` class is designed to analyze protein structures by fetching and downloading PDB files from the UniProt and RCSB PDB databases, parsing the PDB files, and analyzing flexibility based on B-factors.


### Key Methods:

- **`fetch_pdb_ids()`**: Fetches PDB IDs associated with a given UniProt accession number using the UniProt API.
- **`download_pdb_files()`**: Downloads PDB files for the fetched PDB IDs.
- **`load_pdb_files_from_directory()`**: Loads all PDB files from the specified output directory.
- **`parse_pdb_file(pdb_file)`**: Parses a PDB file and stores its structure.
- **`analyze_flexibility()`**: Analyzes the structure to identify regions of high flexibility based on B-factors.
- **`plot_b_factors(save_path=None)`**: Plots the average B-factors along the sequence and indicates the high-flexibility threshold.
- **`get_longest_contiguous_chain_length(pdb_file)`**: Calculates the length of the longest contiguous chain in the PDB file.
- **`select_best_pdb()`**: Selects the best PDB file based on chain length and resolution.

### Example Usage:

1. **Initialize the `ProteinAnalyzer`** with a UniProt accession number:

```python
from proteinfolding import ProteinAnalyzer

analyzer = ProteinAnalyzer("P0DOX5")  # Example UniProt accession
```

2. **Fetch PDB IDs** associated with the protein:

```python
analyzer.fetch_pdb_ids()
```

3. **Download the PDB files** for the fetched PDB IDs:

```python
analyzer.download_pdb_files()
```

4. **Load PDB files** from a directory (optional, if you've already downloaded the files):

```python
analyzer.load_pdb_files_from_directory()
```

5. **Parse a PDB file**:

```python
pdb_file = analyzer.pdb_files[0]  # Use the first PDB file for example
analyzer.parse_pdb_file(pdb_file)
```

6. **Analyze the flexibility** of the protein structure:

```python
analyzer.analyze_flexibility()
```

7. **Plot B-factors** to visualize flexible and rigid regions:

```python
save_path = analyzer.create_save_path(pdb_file)
analyzer.plot_b_factors(save_path=save_path)
```

8. **Assess all PDB files** to identify the best one based on structure quality and chain length:

```python
analyzer.assess_pdb_files()
best_pdb_id = analyzer.select_best_pdb()
```

### Example Output:

After running the above steps, the console will display the PDB IDs, downloaded files, and any highly flexible regions found in the protein. The plots will be saved as images, and the best PDB file will be selected based on chain length and resolution.

### Dependencies:

The `ProteinAnalyzer` class depends on the following Python libraries:
- `requests`
- `BioPython`
- `numpy`
- `matplotlib`

Make sure these dependencies are installed before running the analysis. You can install them using:

```bash
pip install requests biopython numpy matplotlib
```

### Customization:

- **Output Directory**: The output directory for downloaded PDB files and plots can be customized when initializing the `ProteinAnalyzer`.
- **Threshold Adjustment**: The threshold for high flexibility is automatically calculated based on the mean and standard deviation of the B-factors, but you can modify the method to adjust this threshold if necessary.






## Usage

### Fetching Protein Structures

To fetch a protein structure from AlphaFold using its UniProt ID:

```python
from proteinfolding import AlphaFoldStructureFetcher

# Fetch protein structure using UniProt ID
uniprot_id = "P0DOX5"
fetcher = AlphaFoldStructureFetcher(uniprot_id)
pdb_file = fetcher.fetch_structure()
```

This will download the structure to the `structures/` folder and return the path to the PDB file.

### Analyzing Protein Flexibility and Rigidity

Once you've downloaded the PDB file, you can analyze the flexible and rigid regions of the protein based on pLDDT scores:

```python
from proteinfolding import AlphaFoldPDBAnalyzer

# Analyze the downloaded PDB file
analyzer = AlphaFoldPDBAnalyzer(pdb_file)
analyzer.classify_regions()

# Get regions classified as rigid or flexible
rigid_regions = analyzer.get_rigid_regions()
flexible_regions = analyzer.get_flexible_regions()

# Plot pLDDT distribution
analyzer.plot_pLDDT_distribution()

print("Rigid Regions:", rigid_regions)
print("Flexible Regions:", flexible_regions)
```

### Visualizing Proteins with PyMOL

To visualize the protein structure in PyMOL with custom coloring based on pLDDT scores:

```python
from proteinfolding import ProteinVisualizer

visualizer = ProteinVisualizer(pdb_file)
visualizer.visualize_by_pLDDT(output_file="colored_by_pLDDT.pse")
```

This saves the PyMOL session file (`.pse`) with the structure colored according to the pLDDT scores.

### Changing Coloring in PyMOL

If you want to change the coloring spectrum (for example, to set a minimum of 50 and a maximum of 100), you can do so within PyMOL using the following commands:

1. Open PyMOL and load your session file:

    ```bash
    pymol colored_by_pLDDT.pse
    ```

2. In the PyMOL command line, type:

    ```bash
    spectrum b, blue_white_red, minimum=50, maximum=100
    ```

   This command will recolor the protein based on pLDDT values:

- **Blue**: Represents regions with pLDDT scores close to 50 (more flexible).
   - **Red**: Represents regions with pLDDT scores close to 100 (more rigid).

3. To save this modified session, type:

    ```bash
    save recolored_by_pLDDT.pse
    ```

This will save the newly colored session to `recolored_by_pLDDT.pse`.

## Examples

Here’s a complete example of fetching, analyzing, and visualizing a protein structure using the classes provided in the `ProteinFolding` package:

```python
from proteinfolding.main import AlphaFoldStructureFetcher, AlphaFoldPDBAnalyzer, ProteinVisualizer

# Step 1: Fetch the PDB structure from AlphaFold using the UniProt ID
uniprot_id = "P0DOX5"
fetcher = AlphaFoldStructureFetcher(uniprot_id)
pdb_file = fetcher.fetch_structure()

# Step 2: Analyze the flexible and rigid regions based on pLDDT scores
analyzer = AlphaFoldPDBAnalyzer(pdb_file)
analyzer.classify_regions()
rigid_regions = analyzer.get_rigid_regions()
flexible_regions = analyzer.get_flexible_regions()

# Step 3: Visualize the structure using PyMOL and save the session
visualizer = ProteinVisualizer(pdb_file)
visualizer.visualize_by_pLDDT(output_file="colored_by_pLDDT.pse")
```

## Contributing

Contributions are welcome! Please fork the repository and submit a pull request for any improvements or new features.

## License

This package is licensed under the MIT License. See `LICENSE` for details.


