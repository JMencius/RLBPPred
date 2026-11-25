# RLBPPred

Score proteins for potential R-loop binding using the pretrained **SaProt-650M-AF2 + 300-member MLP ensemble**. Start with a UniProt FASTA, a local structure, or previously computed embeddings. No training is required.

## Install and download the model

From the root of **this source checkout**, in a fresh Python 3.10–3.13 environment:

```bash
python -m pip install '.[fasta]'
rlbppred setup
```

The `fasta` extra installs dependencies for **both** FASTA and PDB inputs. Setup downloads and verifies about 3.6 GB once (SaProt, Foldseek and the MLP). The bundled Foldseek download supports macOS and Linux x86-64 with AVX2. The commands below describe the development version, not the published v0.2.1 release.

## Predict: choose one input

**UniProt FASTA** — the included example is human LIG4. The program retrieves a full-length AlphaFold DB structure that matches the FASTA sequence exactly. This route requires internet access and a UniProt accession in each header. It does not fold arbitrary sequences or accept local structures alongside a FASTA.

```bash
rlbppred predict --fasta examples/P49917.fasta --output-dir results_fasta
```

The older `predict-fasta` command remains available for this AFDB-only route.

**Local structure** — try the included human ubiquitin crystal structure ([PDB 1UBQ](https://www.rcsb.org/structure/1UBQ)). The program reads the protein sequence from the coordinates; no FASTA or structure download is needed. Use `--chain` to select one protein chain when several are present.

```bash
rlbppred predict --pdb examples/1UBQ.pdb --chain A --output-dir results_pdb
```

Local structures may be AF2 predictions or experimental coordinates in PDB or mmCIF format. If the source cannot be identified automatically, specify `--structure-source af2` or `--structure-source experimental` **only if you know its origin**. Scores describe the selected, modeled chain or fragment, not necessarily the full-length protein.

**Interpretation:** The model was trained using AFDB structures. For AF2 input, RLBPPred masks structural features at residues with pLDDT below 70. For experimental structures it does **not** use B-factors as pLDDT. Accuracy on experimental input has not been established, and the screening threshold has not been calibrated for this input type.

**Existing embeddings** — provide raw, unstandardized SaProt-AF2 features in an NPZ file with `protein_id` and float32 `X` (N × 1,280):

```bash
rlbppred predict --embeddings features.npz --output predictions.csv
```

## Read the result

For the PDB example, open `results_pdb/predictions.csv`; for the FASTA example, open `results_fasta/predictions.csv`. Each row gives the protein ID, the ensemble score (`probability`), the classification at 0.5 and a screening flag at 0.98 by default. A high score is not proof of binding. A structure-based run also records the input source and selected chain in `structures/provenance.json` and saves a `run.json` receipt. Use a new output directory for each run.

[Model card](MODEL_CARD.md) · [Asset versions and offline setup](docs/assets.md) · [Example provenance](examples/README.md) · [Validation](VALIDATION.md) · [Training and paper reproduction](https://github.com/callmest/RLBPPred-paper-code) · [Third-party notices](THIRD_PARTY_NOTICES.md)
