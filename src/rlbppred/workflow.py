"""Two independent inputs, one fixed encoder and prediction ensemble."""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sys


def predict_fasta(fasta, output_dir, *, saprot_model=None, foldseek=None, model_dir=None,
                  assets_dir=None, device='cpu', batch_size=256, threshold=.98, max_residues=None):
    """Score UniProt FASTA records using sequence-matched structures from AFDB."""
    return _run('fasta', fasta, output_dir, saprot_model=saprot_model, foldseek=foldseek,
                model_dir=model_dir, assets_dir=assets_dir, device=device, batch_size=batch_size,
                threshold=threshold, max_residues=max_residues)


def predict_pdb(pdb, output_dir, *, chain=None, structure_source='auto', saprot_model=None,
                foldseek=None, model_dir=None, assets_dir=None, device='cpu', batch_size=256,
                threshold=.98, max_residues=None):
    """Score the observed protein chain in a local coordinate file, without FASTA or network."""
    return _run('pdb', pdb, output_dir, chain=chain, structure_source=structure_source,
                saprot_model=saprot_model, foldseek=foldseek, model_dir=model_dir,
                assets_dir=assets_dir, device=device, batch_size=batch_size, threshold=threshold,
                max_residues=max_residues)


def _run(mode, input_path, output_dir, *, chain=None, structure_source=None, saprot_model,
         foldseek, model_dir, assets_dir, device, batch_size, threshold, max_residues):
    from . import __version__, encoder, io, runtime, structure
    from .cli import default_assets_dir, default_model_dir

    if not isinstance(batch_size, int) or isinstance(batch_size, bool) or batch_size < 1:
        raise ValueError('batch-size must be positive')
    if not math.isfinite(threshold) or not 0 <= threshold <= 1:
        raise ValueError('threshold must be finite and in [0, 1]')
    if max_residues is not None and (isinstance(max_residues, bool) or not isinstance(max_residues, int) or max_residues < 1):
        raise ValueError('max-residues must be positive')
    destination = Path(output_dir).expanduser()
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(f'Run directory already exists: {destination}')
    input_path = Path(input_path).expanduser()
    if not input_path.is_file():
        raise ValueError(f'Input is not a regular file: {input_path}')
    fasta_mode = mode == 'fasta'
    records = structure.read_fasta(input_path) if fasta_mode else []
    if fasta_mode and max_residues is not None and any(len(seq) > max_residues for _, seq in records):
        raise ValueError('A FASTA sequence exceeds max-residues; no sequences are truncated')
    cache = Path(assets_dir) if assets_dir is not None else default_assets_dir()
    checkpoint = Path(saprot_model) if saprot_model is not None else cache / 'SaProt_650M_AF2.pt'
    binary = Path(foldseek) if foldseek is not None else cache / 'foldseek/bin/foldseek'
    bundle = Path(model_dir) if model_dir is not None else default_model_dir()
    for name, path in [('SaProt checkpoint', checkpoint), ('Foldseek', binary),
                       ('MLP manifest', bundle / 'model.json'), ('MLP weights', bundle / 'model.safetensors')]:
        if not path.is_file():
            raise ValueError(f'{name} not found: {path}. Run rlbppred setup first or supply its path.')
    with input_path.open('rb') as handle:
        digest = hashlib.sha256()
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    destination.mkdir(parents=True)
    receipt = {'rlbppred_version': __version__, 'started_utc': datetime.now(timezone.utc).isoformat(),
               'input_mode': mode, f'{mode}_sha256': digest.hexdigest(),
               'proteins': len(records) if fasta_mode else 1, 'device': device, 'screening_threshold': threshold,
               'operator': '>=', 'status': 'running'}
    receipt_path = destination / 'run.json'
    def save_receipt():
        receipt_path.write_text(json.dumps(receipt, indent=2, allow_nan=False) + '\n')
    save_receipt()
    try:
        if fasta_mode:
            print('Retrieving sequence-matched full-length AF2 structures from AlphaFold DB...', file=sys.stderr)
            resolved = structure.resolve_structures(records, destination)
        else:
            print('Reading the observed protein chain from the local structure...', file=sys.stderr)
            item = structure.resolve_pdb(input_path, chain=chain, source=structure_source, workdir=destination)
            resolved = [item]
            meta = item['provenance']
            receipt.update(structure_source=meta['structure_source'], chain_id=item['chain_id'],
                           sequence_scope=meta['sequence_scope'], scored_residues=len(item['sequence']),
                           structure_provenance='structures/provenance.json')
            if meta['structure_sha256'] != receipt['pdb_sha256']:
                raise ValueError('Structure file changed while preparing this run')
            if max_residues is not None and len(item['sequence']) > max_residues:
                raise ValueError('The selected chain exceeds max-residues; it was not truncated')
            if meta['structure_source'] == 'experimental':
                receipt['warning'] = ('Experimental coordinates: B-factors are not used as pLDDT. '
                                      'Benchmark performance and the screening threshold were not validated '
                                      'for this input mode. The score describes the observed chain only.')
                print(receipt['warning'], file=sys.stderr)
        print('Generating structure-aware AA+3Di tokens...', file=sys.stderr)
        tokens, token_metadata = structure.make_tokens(resolved, binary, destination, plddt_threshold=70.)
        print('Extracting fixed SaProt-650M-AF2 global-mean features...', file=sys.stderr)
        ids, matrix, feature_metadata = encoder.encode_tokens(tokens, checkpoint, device=device,
                                                               output=destination / 'features.npz',
                                                               max_residues=max_residues)
        print('Scoring with the frozen 300-member strict MLP...', file=sys.stderr)
        model = runtime.load_model(bundle, device=device)
        scores = runtime.predict_embeddings(model, matrix, batch_size=batch_size)
        io.write_predictions(destination / 'predictions.csv', ids, scores, threshold)
        receipt.update(status='complete', finished_utc=datetime.now(timezone.utc).isoformat(),
                       model_id=runtime.MODEL_ID, foldseek_version=token_metadata['foldseek_version'],
                       encoder_sha256=feature_metadata['checkpoint_sha256'],
                       above_threshold=int((scores >= threshold).sum()),
                       note='High scores are not automatically novel candidates; deployment FPR/FDR is not guaranteed.')
        save_receipt()
    except Exception as exc:
        receipt.update(status='failed', error=f'{type(exc).__name__}: {exc}')
        save_receipt()
        raise
    return receipt
