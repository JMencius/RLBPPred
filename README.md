# RLBPPred

`RLBPPred` is a command-line tool for ranking proteins as candidate R-loop-binding proteins. It scores each protein with a pretrained SaProt-650M-AF2 encoder and a 300-member MLP ensemble. Start with a UniProt FASTA file, a local protein structure, or precomputed SaProt embeddings.

## Installation
`RLBPPred` only support local installation, but will suport PyPI or Bioconda online installation in the future.

Use Python 3.10–3.13. Install [v0.3.0](https://github.com/callmest/RLBPPred/releases/tag/v0.3.0) from source to run the bundled examples:

```
conda create -n rlbppred python=3.10
git clone --branch v0.3.0 --depth 1 https://github.com/callmest/RLBPPred.git
cd RLBPPred
pip install '.[fasta]'
rlbppred setup
```
The `fasta` extra supports both sequence and structure inputs. Setup downloads and verifies about 3.6 GB once. The bundled Foldseek executable supports modern macOS and Linux x86_64. Prediction runs on **CPU** by default.

## Usage
Full usage of `RLBPPred predict`
```

```


## Examples
1. If your target protein got structure in [AlphaFoldDB](https://alphafold.ebi.ac.uk/)
```
rlbppred predict --fasta examples/P49917.fasta --output-dir results_lig4
```
2. User-defined structure
```bash
rlbppred predict --pdb examples/1UBQ.pdb --chain A --output-dir results_1ubq
rlbppred predict --pdb examples/1UBQ.mmcif --chain A --output-dir results_1ubq
```

## Expected results
Open `results_lig4/predictions.csv`. With AlphaFold DB structure version 6, the example produced this row:

| protein_id | probability | classification_at_0_5 | above_screening_threshold |
| :---: | :---: | :---: | :---: |
| P49917 | 0.991215 | 1 | 1 |

The score is rounded here. `probability` is the mean ensemble score, not a calibrated probability of binding. The last two columns flag scores **≥0.5** and **≥0.98** (the default screening cutoff), respectively. A high score prioritizes a protein for follow-up; it does not establish binding.

## Documentation

- [Model card](MODEL_CARD.md) and [validation results](VALIDATION.md) — training data, screening cutoff, and checks against reference predictions.
- [Installation and offline setup](docs/assets.md) — model downloads, versions, checksums, and cache locations.
- [Example inputs](examples/README.md) — file sources and requirements for local structures.
- [RLBPPred-paper-code](https://github.com/callmest/RLBPPred-paper-code) — training and analysis code for reproducing the paper.


## Acknoledgement
`RLBPPred` is dependent on the following libraries, we are grateful to all the developers/maintainers:
- [SaProt](https://github.com/westlake-repl/Saprot): Structure-aware protein language model
- [Pytorch](https://github.com/pytorch/pytorch): Model training
- [BioPython](https://github.com/biopython/biopython): IO


