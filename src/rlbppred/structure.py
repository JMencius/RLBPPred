"""Prepare Foldseek 3Di from accession-matched AFDB or standalone coordinates."""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlparse
from urllib.request import Request, urlopen

from . import __version__


AMINO_ACIDS = frozenset('ACDEFGHIKLMNPQRSTVWY')
STRUCTURE_ALPHABET = frozenset('pynwrqhgdlvtmfsaeikc#')
AFDB_API = 'https://alphafold.ebi.ac.uk/api/prediction/'
_ID = re.compile(r'[A-Za-z0-9][A-Za-z0-9_.:|+\-]*\Z')
_UNIPROT = re.compile(
    r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-[1-9][0-9]*)?\Z')
_AF2 = re.compile(r'ALPHAFOLD(?:[\s_-]*2\b|[\s_-]+(?:MONOMER|MULTIMER)[\s_-]+V?2(?:\.|\b))')
_THREE_TO_ONE = dict(zip(
    'ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR'.split(),
    'ACDEFGHIKLMNPQRSTVWY'))


def _digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def _sequence_digest(sequence):
    return hashlib.sha256(sequence.encode('ascii')).hexdigest()


def _validate_records(records):
    if not isinstance(records, (list, tuple)) or not records:
        raise ValueError('FASTA must contain at least one protein')
    seen = set()
    for protein, sequence in records:
        if not isinstance(protein, str) or not _ID.fullmatch(protein):
            raise ValueError('Protein identifiers must be safe nonempty ASCII identifiers')
        if protein in seen:
            raise ValueError(f'Duplicate FASTA identifier: {protein}')
        seen.add(protein)
        if not isinstance(sequence, str) or not sequence or set(sequence) - AMINO_ACIDS:
            raise ValueError(f'{protein}: expected a nonempty sequence containing the 20 standard amino acids')


def read_fasta(path):
    """Read standard protein FASTA; sp|ACCESSION|NAME headers become ACCESSION.

    Simple first-token identifiers are preserved. Sequence letters are
    uppercased; whitespace is ignored, but X, stop characters and gaps are not.
    """
    records = []
    protein = None
    pieces = []
    with Path(path).open(encoding='utf-8-sig') as stream:
        for number, raw in enumerate(stream, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if protein is not None:
                    records.append((protein, ''.join(pieces).upper()))
                header = line[1:].split()
                if not header:
                    raise ValueError(f'FASTA line {number}: empty header')
                protein = header[0]
                fields = protein.split('|')
                if fields[0] in ('sp', 'tr'):
                    if len(fields) != 3 or not fields[1]:
                        raise ValueError(f'FASTA line {number}: malformed UniProt header')
                    protein = fields[1]
                pieces = []
            else:
                if protein is None:
                    raise ValueError(f'FASTA line {number}: sequence appears before a header')
                pieces.append(''.join(line.split()))
    if protein is not None:
        records.append((protein, ''.join(pieces).upper()))
    _validate_records(records)
    return records


def _new_subdir(workdir, name):
    workdir = Path(workdir).expanduser()
    if workdir.is_symlink():
        raise ValueError('Run directory must not be a symbolic link')
    workdir.mkdir(parents=True, exist_ok=True)
    target = workdir / name
    if target.exists() or target.is_symlink():
        raise FileExistsError(f'Output directory already exists: {target}')
    target.mkdir()
    return target


def _checked_url(url, api=False):
    parsed = urlparse(url)
    prefix = '/api/prediction/' if api else '/files/'
    if (parsed.scheme != 'https' or parsed.hostname != 'alphafold.ebi.ac.uk'
            or parsed.username or parsed.password or parsed.port not in (None, 443)
            or not parsed.path.startswith(prefix)):
        raise ValueError('Expected an official HTTPS AlphaFold DB URL')
    return url


def _fetch_bytes(url, *, api=False, max_bytes=100_000_000):
    """Bounded public download; no user sequences or credentials are submitted."""
    _checked_url(url, api=api)
    request = Request(url, headers={'User-Agent': f'RLBPPred/{__version__} structure retrieval',
                                    'Accept': 'application/json' if api else '*/*'})
    try:
        with urlopen(request, timeout=45) as response:
            _checked_url(response.geturl(), api=api)
            content = response.read(max_bytes + 1)
    except (HTTPError, URLError, OSError, TimeoutError) as exc:
        raise ValueError(f'AlphaFold DB retrieval failed: {exc}') from exc
    if len(content) > max_bytes:
        raise ValueError('AlphaFold DB response exceeds the safety size limit')
    return content


_EXPERIMENT = re.compile(r'\b(?:DIFFRACTION|NMR|MICROSCOPY|SCATTERING|CRYSTALLOGRAPHY|TRANSFER|TOMOGRAPHY)\b')
_OTHER_PREDICTOR = re.compile(r'ALPHAFOLD[\s_-]*3\b|\bESMFOLD\b|\bROSETTAFOLD\w*\b|\bBOLTZ\b')


def _local_identity(path, source):
    """Recognize declared structure origin, never infer it from numeric B-values."""
    if source not in ('auto', 'af2', 'experimental'):
        raise ValueError('structure-source must be auto, af2 or experimental')
    if path.suffix.lower() == '.pdb':
        titles, methods = [], []
        with path.open(encoding='ascii', errors='replace') as inp:
            for line in inp:
                tag = line[:6].strip()
                if tag in ('ATOM', 'HETATM', 'MODEL'):
                    break
                if tag in ('TITLE', 'HEADER', 'EXPDTA'):
                    (methods if tag == 'EXPDTA' else titles).append(line[6:].strip().upper())
    else:
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict
        meta = MMCIF2Dict(str(path))
        methods = [v.upper() for v in meta.get('_exptl.method', ())]
        titles = [v.upper() for v in meta.get('_struct.title', ())]
        names = meta.get('_software.name', ())
        versions = meta.get('_software.version', ())
        titles.extend((name + ' ' + (versions[i] if i < len(versions) else '')).upper()
                      for i, name in enumerate(names))
    header = ' '.join(titles + methods)
    recognized = None
    if any(_EXPERIMENT.search(m) for m in methods):
        recognized = 'experimental'
    elif _OTHER_PREDICTOR.search(header):
        raise ValueError('Unsupported prediction method: expected AlphaFold2 or an experimental structure')
    elif _AF2.search(header):
        recognized = 'af2'
    theoretical = any('THEORETICAL' in m or 'COMPUT' in m for m in methods)
    metadata = {}
    declared = None
    sidecar = Path(str(path) + '.json')
    if sidecar.is_file() and not sidecar.is_symlink():
        declaration = json.loads(sidecar.read_text(encoding='utf-8'))
        if not isinstance(declaration, dict):
            raise ValueError('Expected a JSON object in the structure sidecar')
        origin = declaration.get('structure_source')
        if origin in ('alphafold2', 'af2') and declaration.get('bfactor_is_plddt') is True:
            declared = 'af2'
        elif origin == 'experimental' and declaration.get('bfactor_is_plddt', False) is False:
            declared = 'experimental'
        else:
            raise ValueError('Invalid structure sidecar source or pLDDT declaration')
        metadata['sidecar_sha256'] = _digest(sidecar)
    origins = {s for s in (recognized, declared, source if source != 'auto' else None) if s}
    if len(origins) > 1 or (theoretical and 'experimental' in origins):
        raise ValueError('Conflicting structure source in coordinate metadata, sidecar or --structure-source')
    if not origins:
        raise ValueError('Unknown structure source. Set --structure-source af2 only for AF2 with pLDDT '
                         'in B-factors, or --structure-source experimental for experimental coordinates.')
    kind = origins.pop()
    return {**metadata, 'structure_source': kind, 'bfactor_is_plddt': kind == 'af2',
            'source_evidence': 'structure metadata' if recognized else ('sidecar' if declared else 'command line')}


def _model(path):
    try:
        from Bio.PDB import MMCIFParser, PDBParser
        from Bio.PDB.PDBExceptions import PDBConstructionException
    except ImportError as exc:
        raise RuntimeError('Biopython is required; install the fasta package extra') from exc
    parser = MMCIFParser(QUIET=True) if path.suffix.lower() in ('.cif', '.mmcif') else PDBParser(PERMISSIVE=False, QUIET=True)
    try:
        structure = parser.get_structure('input', str(path))
        models = list(structure.get_models())
        if len(models) != 1:
            raise ValueError('Exactly one structure model is required')
        return models[0]
    except (KeyError, TypeError, OSError, PDBConstructionException) as exc:
        raise ValueError(f'{path.name}: cannot parse a valid structure: {exc}') from exc


def _protein_like(residue):
    from Bio.PDB.Polypeptide import is_aa
    return is_aa(residue) or all(key in residue for key in ('CA', 'N', 'C'))


def _read_chain(chain, af2):
    import numpy as np
    if not re.fullmatch(r'[A-Za-z0-9]+', chain.id):
        raise ValueError('Blank or special-character chain identifiers are unsupported')
    sequence, confidence, residue_ids = [], [], []
    for residue in chain:
        flag, number, ins = residue.id
        name = residue.get_resname().upper()
        if flag.strip() and not _protein_like(residue):
            continue  # solvent, nucleotides in hetero records, bound ligands
        if flag.strip() or name not in _THREE_TO_ONE:
            raise ValueError(f'Nonstandard protein residue {name} in chain {chain.id}; cannot silently skip it')
        if residue.is_disordered():
            raise ValueError('Alternate locations/ambiguous residues: select a single conformer before scoring')
        atoms = list(residue.get_atoms())
        if 'CA' not in residue or any(atom.is_disordered() for atom in atoms):
            raise ValueError('Missing C-alpha or alternate-location atoms are unsupported')
        if residue_ids and number not in (residue_ids[-1][0], residue_ids[-1][0]+1):
            raise ValueError(f'Noncontiguous residue numbering in chain {chain.id}: scoring across a '
                             'possible missing internal segment is unsupported')
        residue_ids.append((number, ins.strip()))
        sequence.append(_THREE_TO_ONE[name])
        if af2:
            values = [float(atom.get_bfactor()) for atom in atoms]
            if any(not math.isfinite(v) or not 0 <= v <= 100 for v in values):
                raise ValueError('pLDDT values must be finite and in [0, 100]')
            confidence.append(float(np.mean(values)))
    if not sequence:
        raise ValueError(f'No observed protein residues in chain {chain.id}')
    return ''.join(sequence), confidence if af2 else None, [f'{n}{ins}' for n, ins in residue_ids]


def _matching_chain(path, expected_sequence):
    matches, lengths = [], []
    for chain in _model(path):
        if not any(_protein_like(res) for res in chain):
            continue
        observed, confidence, _ = _read_chain(chain, af2=True)
        lengths.append(len(observed))
        if observed == expected_sequence:
            matches.append((chain.id, confidence))
    if len(matches) != 1:
        if matches:
            raise ValueError('Multiple chains exactly match FASTA; chain choice is ambiguous')
        raise ValueError(f'No complete structure chain exactly matches FASTA '
                         f'(FASTA length {len(expected_sequence)}, structure chain lengths {lengths})')
    return matches[0]


def _select_afdb_record(payload, protein, sequence):
    if not isinstance(payload, list) or not payload:
        raise ValueError(f'{protein}: no AlphaFold DB prediction records')
    matches = []
    for row in payload:
        if not isinstance(row, dict):
            continue
        if (row.get('uniprotAccession') != protein or row.get('uniprotSequence') != sequence
                or row.get('sequence', sequence) != sequence):
            continue
        start = row.get('uniprotStart', row.get('sequenceStart'))
        end = row.get('uniprotEnd', row.get('sequenceEnd'))
        try:
            complete = int(start) == 1 and int(end) == len(sequence)
            version = int(row.get('latestVersion'))
        except (ValueError, TypeError):
            continue
        tool = row.get('toolUsed', '')
        entry = row.get('entryId', row.get('modelEntityId'))
        if (complete and version >= 1 and isinstance(tool, str) and _AF2.search(tool.upper())
                and entry == f'AF-{protein}-F1' and row.get('isComplex', False) is False):
            matches.append(row)
    if len(matches) != 1:
        raise ValueError(f'{protein}: expected exactly one full-length AFDB record with the exact FASTA sequence; '
                         'fragments, isoform mismatches and ambiguous records are not silently selected')
    return matches[0]


def resolve_structures(records, workdir):
    """Remote-only FASTA route: require one exact full-length AFDB AF2 entry per accession."""
    _validate_records(records)
    destination = _new_subdir(workdir, 'structures')
    resolved = []
    for protein, sequence in records:
        if not _UNIPROT.fullmatch(protein):
            raise ValueError(f'{protein}: AFDB retrieval requires a UniProt accession; '
                             'use --pdb instead to score a local structure')
        api_url = AFDB_API + quote(protein, safe='')
        raw = _fetch_bytes(api_url, api=True, max_bytes=20_000_000)
        try:
            row = _select_afdb_record(json.loads(raw), protein, sequence)
        except (json.JSONDecodeError, UnicodeError) as exc:
            raise ValueError(f'{protein}: invalid AlphaFold DB metadata JSON') from exc
        with (destination / (protein + '.afdb.json')).open('xb') as target:
            target.write(raw)
        url = row.get('pdbUrl') or row.get('cifUrl')
        if not isinstance(url, str):
            raise ValueError(f'{protein}: AFDB did not provide a structure URL')
        _checked_url(url)
        suffix = Path(urlparse(url).path).suffix.lower()
        if suffix not in ('.pdb', '.cif', '.mmcif'):
            raise ValueError(f'{protein}: unsupported structure file format')
        entry_id = row.get('entryId', row.get('modelEntityId'))
        expected_name = f'{entry_id}-model_v{int(row["latestVersion"])}{suffix}'
        if Path(urlparse(url).path).name != expected_name:
            raise ValueError(f'{protein}: structure URL does not match AFDB entry/version metadata')
        content = _fetch_bytes(url)
        path = destination / (protein + suffix)
        with path.open('xb') as target:
            target.write(content)
        chain, confidence = _matching_chain(path, sequence)
        provenance = {'source': 'alphafold_database', 'structure_source': 'af2', 'sequence_scope': 'full_length_fasta',
                      'accession': protein,
                      'entry_id': entry_id, 'version': int(row['latestVersion']),
                      'structure_method': row['toolUsed'],
                      'metadata_url': api_url, 'structure_url': url,
                      'retrieved_utc': datetime.now(timezone.utc).isoformat(),
                      'bfactor_is_plddt': True}
        provenance.update(structure_sha256=_digest(path), sequence_sha256=_sequence_digest(sequence),
                          length=len(sequence), chain_id=chain, complete_sequence_match=True)
        resolved.append({'protein_id': protein, 'sequence': sequence, 'path': str(path),
                         'chain_id': chain, 'plddt': confidence, 'provenance': provenance})
    with (destination / 'provenance.json').open('x', encoding='utf-8') as target:
        json.dump([r['provenance'] for r in resolved], target, indent=2, allow_nan=False)
        target.write('\n')
    return resolved


def resolve_pdb(path, *, chain=None, source='auto', workdir=None):
    """Select one observed protein chain from a standalone PDB/mmCIF; no FASTA or network.

    The coordinates alone do not establish full-length canonical protein coverage.
    With a workdir, retain an immutable run snapshot with its source declaration.
    """
    original = Path(path).expanduser().resolve(strict=True)
    suffix = original.suffix.lower()
    if not original.is_file() or suffix not in ('.pdb', '.cif', '.mmcif'):
        raise ValueError('Expected an uncompressed .pdb, .cif or .mmcif coordinate file')
    stem = re.sub(r'[^A-Za-z0-9_.-]+', '_', original.stem).strip('._-') or 'structure'
    destination = None
    path = original
    if workdir is not None:
        destination = _new_subdir(workdir, 'structures')
        path = destination / (stem + suffix)
        shutil.copyfile(original, path)
        sidecar = Path(str(original) + '.json')
        if sidecar.is_file() and not sidecar.is_symlink():
            shutil.copyfile(sidecar, str(path) + '.json')
    model = _model(path)
    identity = _local_identity(path, source)
    choices = {c.id: c for c in model if any(_protein_like(r) for r in c)}
    if not choices:
        raise ValueError('No protein chains found in the coordinate file')
    available = ', '.join(repr(c) for c in choices)
    if chain is None:
        if len(choices) != 1:
            raise ValueError(f'Multiple protein chains: {available}; specify --chain')
        chain = next(iter(choices))
    if chain not in choices:
        raise ValueError(f'Unknown protein chain {chain!r}; available protein chains: {available}')
    sequence, confidence, residue_ids = _read_chain(choices[chain], af2=identity['structure_source']=='af2')
    protein = f'{stem}:{chain}'
    provenance = {'source': 'local_structure', **identity, 'source_file': original.name,
                  'structure_sha256': _digest(path), 'sequence_sha256': _sequence_digest(sequence),
                  'protein_id': protein, 'length': len(sequence), 'chain_id': chain,
                  'sequence_scope': 'observed_chain', 'full_length_verified': False,
                  'first_observed_residue': residue_ids[0], 'last_observed_residue': residue_ids[-1],
                  'model_serial': model.serial_num}
    if destination is not None:
        with (destination / 'provenance.json').open('x', encoding='utf-8') as f:
            json.dump([provenance], f, indent=2, allow_nan=False)
            f.write('\n')
    return {'protein_id': protein, 'sequence': sequence, 'path': str(path),
            'chain_id': chain, 'plddt': confidence, 'provenance': provenance}


def make_tokens(resolved, foldseek, workdir, plddt_threshold=70.0):
    """Return (AA+3Di records, provenance) without implicit trimming or skipping."""
    if not math.isfinite(plddt_threshold) or not 0 <= plddt_threshold <= 100:
        raise ValueError('pLDDT threshold must be finite and in [0, 100]')
    if not resolved:
        raise ValueError('No resolved structures were supplied')
    executable = shutil.which(str(foldseek)) or str(Path(foldseek).expanduser())
    binary = Path(executable).resolve(strict=True)
    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise ValueError('Foldseek must be an executable regular file')
    try:
        version = subprocess.run([str(binary), 'version'], capture_output=True, text=True, timeout=30, check=False)
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise ValueError(f'Cannot execute Foldseek on this platform: {exc}') from exc
    if version.returncode or not version.stdout.strip():
        raise ValueError('Foldseek version command failed; supply a compatible, versioned binary')
    destination = _new_subdir(workdir, 'tokens')
    records, rows = [], []
    for item in resolved:
        protein, sequence, chain = item['protein_id'], item['sequence'], item['chain_id']
        _validate_records([(protein, sequence)])
        confidence = item['plddt']
        kind = item['provenance']['structure_source']
        if kind not in ('af2', 'experimental'):
            raise ValueError(f'{protein}: unknown structure source')
        if kind == 'experimental' and confidence is not None:
            raise ValueError(f'{protein}: experimental B-factors must not be supplied as pLDDT')
        if kind == 'af2' and (confidence is None or len(confidence) != len(sequence) or
                             any(not math.isfinite(v) or not 0 <= v <= 100 for v in confidence)):
            raise ValueError(f'{protein}: pLDDT coverage does not match the complete sequence')
        path = Path(item['path']).resolve(strict=True)
        if _digest(path) != item['provenance']['structure_sha256']:
            raise ValueError(f'{protein}: structure changed after sequence/confidence validation')
        with tempfile.TemporaryDirectory(prefix='foldseek-', dir=destination) as temporary:
            descriptor = Path(temporary) / 'descriptor.tsv'
            command = [str(binary), 'structureto3didescriptor', '-v', '0', '--threads', '1',
                       '--chain-name-mode', '1', str(path), str(descriptor)]
            try:
                result = subprocess.run(command, capture_output=True, text=True, timeout=300, check=False)
            except (OSError, subprocess.TimeoutExpired) as exc:
                raise ValueError(f'{protein}: Foldseek execution failed: {exc}') from exc
            if result.returncode or not descriptor.is_file():
                raise ValueError(f'{protein}: Foldseek descriptor generation failed '
                                 f'(exit {result.returncode}): {result.stderr[-1000:]}')
            matches = []
            for line in descriptor.read_text().splitlines():
                if not line.strip():
                    continue
                columns = line.split('\t')
                if len(columns) < 3:
                    raise ValueError(f'{protein}: malformed Foldseek descriptor row')
                name, amino_acids, structural = columns[:3]
                label = name.split()[0]
                if '_' not in label:
                    raise ValueError(f'{protein}: Foldseek output lacks explicit chain identifiers')
                prefix, result_chain = label.rsplit('_', 1)
                if result_chain == chain and Path(prefix).name in (path.name, path.stem):
                    matches.append((amino_acids, structural.lower()))
            if len(matches) != 1:
                raise ValueError(f'{protein}: Foldseek output did not identify exactly the selected structure chain')
            amino_acids, structural = matches[0]
            if amino_acids != sequence or len(structural) != len(sequence):
                raise ValueError(f'{protein}: Foldseek AA/3Di does not exactly cover the selected sequence')
            if set(structural) - STRUCTURE_ALPHABET:
                raise ValueError(f'{protein}: Foldseek returned unsupported 3Di symbols')
            masked = (structural if confidence is None else
                      ''.join('#' if score < plddt_threshold else symbol
                              for symbol, score in zip(structural, confidence)))
            combined = ''.join(aa + symbol for aa, symbol in zip(sequence, masked))
        records.append((protein, combined))
        rows.append({'protein_id': protein, 'length': len(sequence), 'chain_id': chain,
                     'structure_sha256': item['provenance']['structure_sha256'],
                     'token_sha256': hashlib.sha256(combined.encode('ascii')).hexdigest(),
                     'structure_source': kind,
                     'masked_residues': 0 if confidence is None else sum(score < plddt_threshold for score in confidence),
                     'plddt_threshold': None if confidence is None else plddt_threshold})
    provenance = {'foldseek_version': version.stdout.strip(), 'foldseek_sha256': _digest(binary),
                  'mode': 'structureto3didescriptor', 'chain_name_mode': 1, 'threads': 1,
                  'plddt_policy': 'AF2: mask 3Di below the mean per-residue pLDDT cutoff; experimental: no B-factor masking',
                  'records': rows}
    with (destination / 'tokens.fasta').open('x', encoding='ascii') as target:
        for protein, combined in records:
            target.write(f'>{protein}\n{combined}\n')
    with (destination / 'provenance.json').open('x', encoding='utf-8') as target:
        json.dump(provenance, target, indent=2, allow_nan=False)
        target.write('\n')
    return records, provenance
