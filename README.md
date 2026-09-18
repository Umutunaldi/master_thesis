# Master thesis analysis

This repository contains the Python and R scripts used for the structural
analysis of isoform-specific protein interactions in the master thesis. It also
contains the selected MODELLER structures, their model scores, the FoldX
interaction-energy outputs, and the final tab-separated result tables.

## Main files

- `thesis_paper.py`: complete Python analysis
- `environment.yml`: conda environment used for Python, MODELLER, MMseqs2 and R
- `structman_output/`: StructMAn protein-protein interaction result used as input
- `modeller_complex_runs/`: selected (best-ranked) MODELLER structure and score table for each job
- `foldx_results/`: saved FoldX `AnalyseComplex` output for each selected model
- `summary_results/`: combined result tables used by the R scripts
- `interface_analysis/`: mapped residue contacts used for the networks and contact maps
- `supplementary_tables/`: final supplementary tables
- `r_scripts/`: R scripts used to generate the quantitative figures

## Conda environment

environment.yml

conda env create -f environment.yml
conda activate thesis

## Python script

Python script should be run in the same folder.

python thesis_paper.py



The selected (best-ranked) MODELLER structures and FoldX result files are already included. To run MODELLER and FoldX again, need to remove the modeller_complex_runs and foldx_results directories.


## AlphaFold templates

The AlphaFold and AlphaFold-Multimer complex template files are listed in
`required_af_templates.tsv`. They need to be downloaded. This repo doesn't provide that since that might require a permission.


For a complete remodelling run, save the five PDB files in a directory named
`af_models` and use the exact template filenames listed in the manifest. The
already selected MODELLER models allow the downstream interface and FoldX
analyses to be reproduced without redistributing these template files.

## FoldX

FoldX is not included in this repository or in `environment.yml`. It must be
obtained separately from the official FoldX website:

<https://foldxsuite.crg.eu/>

The existing 52 FoldX interaction-energy files are included, so the script
does not call FoldX when those outputs are present. To calculate a new or
missing job, FoldX needs to be separately installed and its executable should be available in the
system path under the command name `foldx`. The script uses `RepairPDB` followed
by `AnalyseComplex` with the two model-chain identifiers.

## R figures

Run the R scripts from the repository root after the Python result tables are
present:

```bash
Rscript r_scripts/interaction_networks.R
Rscript r_scripts/contact_maps.R
Rscript r_scripts/summary_figures.R
```

The figures are written to the `figures` directory. The final analysis uses an
interactor-chain C-alpha RMSD threshold of 2 Å.

## Notes

The supplied StructMAn table was produced with the development version used in
the thesis. StructMAn itself is available from <https://tools.helmholtz-hips.de/structman> and <https://github.com/kalininalab/d-StructMAn>.
summary_results contains generated outputs included for reproducibility. It could be deleted and the script would still generate it.
