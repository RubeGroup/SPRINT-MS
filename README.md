# SPRINT-MS
Repository containing code base for the Pool Solver, Pool Designer and notebooks for recreating results from [SPRINT-MS manuscript](https://www.biorxiv.org/content/10.1101/2025.11.21.689547v1).
This repository contains code for simulating dataset and for processing data from actual experiments. 

# To make this work:
 1. Clone the github repository `SPRINT-MS`
 2. Add the following lines (with updated absolute) to `~/.profile`:
 `export IPMS_DIR="/path/to/git/SPRINT-MS"`
 `export PYTHONPATH="$IPMS_DIR:$PYTHONPATH"`
 3. Install Python>=3.0, scipy, numpy, matplotlib and scikit-learn. 

 
# Conda environment:

For example to make a conda environment to run SPRINT-MS:

- `conda create --name sprintms` - this creates conda environment called sprintms
- `source activate sprintms` - this activates environment
- `conda install scipy numpy matplotlib scikit-learn` 


# Data included:

- `rawData/` -  `20250404_PooledIP_15x30/` contains mixing scheme and output data from pooled antibody experiments
and `20250512_Pooled_Lysate_IPs/` contains mixing scheme and output data from pooled lysate experiments. 
- `standardizedData/` - contains data in a standard form, used as input for code in `poolData/` and `poolSolver/`.
`20250404_PooledIP_15x30/` for the Ab experiments and `20250512_Pooled_Lysate_IPs/` for lysate experiments, and 
and `README` for generating this standardized data.
- `referenceData/`- contains PPI data from various literature sources 
 `20250813_human_literature_interactions_UniProtIDs.tsv` for pooled AB experiments and 
 `20250813_yeast_literature_interactions_UniProtIDs.tsv` for pooled lysate experiments.

# Code organization:

- `poolData/` - utility funtions in `protein.py`, `pooledDataset.py` and `proteinProteinMatrix.py` for processing and plotting data from experiments,
`interactionList.py` contains code to process data from literature (`referenceData/knownPPIs`).
- `poolSolver/` - different functions for versions of Pool Solver algorithms (preliminary versions are `leastSquaresSolver.py`, `nnlsSolver.py`, `correlationSolver.py` and `subsetSelectionProteinSolvers.py`, final Pool Solver is in `bestSubsetSelectionPoolSolver.py`)
- `poolDesigner/` - code for generating mixing matrix (Pool Designer algorithm)
- `poolSimulator/` - code for simulating experiments, with different number of PPIs, noise levels etc. 



# Examples provided:

The following jupyter notebooks include some examples to run the platform:

- `Simulated_Experiment` - code for generating figure 1 and supplementary figure 2 in manuscript, 
which illustrate outputs from simulated experiments.
- `Pooled_Ab_Experiments` - code for generating figure 2 in manuscript, showing data from pooled Ab experiments.
- `Pooled_Lysate_Experiments` - code for generating figure 3 in manuscript, showing data from pooled lysate experiments.
- `Pool_Designer` - code for generating mixing matrices with different pools and antibodies


 


