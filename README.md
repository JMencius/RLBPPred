# RLBPPred

RLBPPred is a command-line tool for ranking proteins as candidate R-loop-binding proteins. It scores each protein with a pretrained SaProt-650M-AF2 encoder and a 300-member MLP ensemble. Start with a UniProt FASTA file, a local protein structure, or precomputed SaProt embeddings.

## Installation

Use Python 3.10–3.13. Install [v0.3.0](https://github.com/callmest/RLBPPred/releases/tag/v0.3.0) from source to run the bundled examples:

```bash
git clone --branch v0.3.0 --depth 1 https://github.com/callmest/RLBPPred.git
cd RLBPPred
python -m venv .venv
source .venv/bin/activate
python -m pip install '.[fasta]'
rlbppred setup
```

The `fasta` extra supports both sequence and structure inputs. Setup downloads and verifies about 3.6 GB once. The bundled Foldseek executable supports macOS and Linux x86-64 with AVX2. Prediction runs on CPU by default.

## Usage

### UniProt FASTA

Run the included human DNA ligase 4 (LIG4) example from the repository root:

```bash
rlbppred predict --fasta examples/P49917.fasta --output-dir results_lig4
```

Open `results_lig4/predictions.csv`. With AlphaFold DB structure version 6, the example produced this row:

| protein_id | probability | classification_at_0_5 | above_screening_threshold |
| --- | ---: | ---: | ---: |
| P49917 | 0.991215 | 1 | 1 |

The score is rounded here. `probability` is the mean ensemble score, not a calibrated probability of binding. The last two columns flag scores **≥0.5** and **≥0.98** (the default screening cutoff), respectively. A high score prioritizes a protein for follow-up; it does not establish binding.

For your own FASTA file, include a UniProt accession in each header (for example, `>P49917`). RLBPPred retrieves the full-length structure from AlphaFold DB and checks for an exact sequence match. This requires internet access; the tool does not fold a sequence for which no matching structure is available. Use a new output directory for each run.

### Local structure (PDB or mmCIF)

To score chain A of the included human ubiquitin crystal structure, run:

```bash
rlbppred predict --pdb examples/1UBQ.pdb --chain A --output-dir results_1ubq
```

This writes `results_1ubq/predictions.csv`. No FASTA file or internet connection is needed after setup. For a structure with several protein chains, select one with `--chain`. The score applies only to the residues present in that chain. If the source is not recognized, specify `--structure-source af2` or `--structure-source experimental` according to how the structure was obtained.

The model was trained using AlphaFold DB structures. For AF2 input, structural features are masked where pLDDT is below 70; experimental B-factors are not treated as pLDDT. Neither predictive accuracy nor the 0.98 screening cutoff has been validated for experimental structures.

### Precomputed embeddings (NPZ)

If you already have SaProt-650M-AF2 embeddings, skip Foldseek and the encoder download. From the cloned repository, use the minimal installation below instead of the full setup. Your NPZ file must contain `protein_id` (one ID per row) and `X` (raw, unstandardized `float32` vectors of shape `N × 1280`).

```bash
python -m pip install .
rlbppred setup --components model
rlbppred predict --embeddings features.npz --output predictions.csv
```

## Documentation

- [Model card](MODEL_CARD.md) and [validation results](VALIDATION.md) — training data, screening cutoff, and checks against reference predictions.
- [Installation and offline setup](docs/assets.md) — model downloads, versions, checksums, and cache locations.
- [Example inputs](examples/README.md) — file sources and requirements for local structures.
- [RLBPPred-paper-code](https://github.com/callmest/RLBPPred-paper-code) — training and analysis code for reproducing the paper.

RLBPPred is distributed under the [MIT license](LICENSE). See [third-party notices](THIRD_PARTY_NOTICES.md) for the tools, models, and data it uses.
