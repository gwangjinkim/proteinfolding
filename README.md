# ProteinFolding Package

`ProteinFolding` is a Python package that allows you to fetch and analyze protein structures, particularly those predicted by AlphaFold. This package provides functionality to:
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

You can install this package by cloning the repository and installing its dependencies.

```bash
git clone https://github.com/yourusername/proteinfolding.git
cd proteinfolding
pip install -r requirements.txt
```

Ensure you have PyMOL installed to use the visualization features:
- Install PyMOL: [PyMOL installation guide](https://pymol.org/2/#download)

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

## Usage

### Fetching Protein Structures

To fetch a protein structure from AlphaFold using its UniProt ID:

```python
from proteinfolding.main import AlphaFoldStructureFetcher

# Fetch protein structure using UniProt ID
uniprot_id = "P0DOX5"
fetcher = AlphaFoldStructureFetcher(uniprot_id)
pdb_file = fetcher.fetch_structure()
```

This will download the structure to the `structures/` folder and return the path to the PDB file.

### Analyzing Protein Flexibility and Rigidity

Once you've downloaded the PDB file, you can analyze the flexible and rigid regions of the protein based on pLDDT scores:

```python
from proteinfolding.main import AlphaFoldPDBAnalyzer

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
from proteinfolding.main import ProteinVisualizer

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
