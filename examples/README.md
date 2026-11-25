# Example inputs

Run these from the root of the source checkout after installation and `rlbppred setup`.

| Input | Command | Result |
|---|---|---|
| Human DNA ligase 4, UniProt [P49917](https://www.uniprot.org/uniprotkb/P49917/entry) | `rlbppred predict --fasta examples/P49917.fasta --output-dir results_fasta` | `results_fasta/predictions.csv` |
| Human ubiquitin crystal structure, chain A, [PDB 1UBQ](https://www.rcsb.org/structure/1UBQ) | `rlbppred predict --pdb examples/1UBQ.pdb --chain A --output-dir results_pdb` | `results_pdb/predictions.csv` |

The FASTA example retrieves its structure from AlphaFold DB. The PDB example uses only the included coordinates and cached model assets; it does not need a FASTA or a network connection. For this example the output row identifies `1UBQ:A`, the observed ubiquitin chain, rather than a UniProt protein inferred from the filename. This demonstrates the software interface, not experimental evidence for R-loop binding.

For another structure, use one model and an unambiguous protein chain. The reader rejects unsupported modifications, alternate conformations and breaks in residue numbering rather than silently joining or dropping them. A continuous fragment can be scored, but it is not automatically a full-length protein. Inspect the selected chain and observed residue range in `structures/provenance.json` inside the output directory.

## Where the example structure came from

`1UBQ.pdb` is the unmodified [RCSB PDB coordinate file](https://files.rcsb.org/download/1UBQ.pdb) for entry 1UBQ (archive revision 1.3, 14 February 2024). It contains one 76-residue human ubiquitin chain solved by X-ray diffraction at 1.80 Å; all 76 deposited protein residues have coordinates.

- File SHA-256: `d4a6812d8951cf6594e6a0763f089e35f5a80b62acb3c117b2c5565228a7b161`.
- Primary citation: Vijay-Kumar, S., Bugg, C. E. & Cook, W. J. Structure of ubiquitin refined at 1.8 Å resolution. *J. Mol. Biol.* **194**, 531–544 (1987). [DOI](https://doi.org/10.1016/0022-2836%2887%2990679-6).
- PDB archive data are distributed under the [wwPDB CC0 policy](https://www.wwpdb.org/about/usage-policies).
