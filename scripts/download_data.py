"""Restore pinned owner mirrors. These hashes do not authenticate organizer provenance."""
import concurrent.futures as cf
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1048576), b''): h.update(b)
    return h.hexdigest()

def fetch(repo, ref, path, dest):
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + '.partial')
    try:
        with tmp.open('wb') as f:
            subprocess.run(['gh','api',f'repos/{repo}/contents/{path}?ref={ref}',
                            '-H','Accept: application/vnd.github.raw'], stdout=f, check=True, timeout=300)
        tmp.replace(dest)
    finally: tmp.unlink(missing_ok=True)

def restore(e):
    dest = ROOT / 'data' / e['dest']
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists() or sha(dest) != e['sha256']:
        if 'parts' in e:
            paths = [ROOT/'data'/Path(p).name for p in e['parts']]
            with cf.ThreadPoolExecutor(max_workers=3) as pool:
                list(pool.map(lambda pair: fetch(e['repo'],e['ref'],pair[0],pair[1]),zip(e['parts'],paths)))
            tmp = dest.with_suffix('.partial')
            with tmp.open('wb') as f:
                for p in paths:
                    with p.open('rb') as inp:
                        for b in iter(lambda: inp.read(1048576),b''): f.write(b)
                    p.unlink()
            tmp.replace(dest)
        else: fetch(e['repo'], e['ref'], e['path'], dest)
    if sha(dest) != e['sha256']:
        dest.unlink()
        raise ValueError(f"Checksum mismatch: {e['id']}")
    print('Verified', e['id'], dest.stat().st_size, flush=True)
    return dict(id=e['id'], sha256=e['sha256'], bytes=dest.stat().st_size, provenance=e['provenance'])

if __name__ == '__main__':
    manifest = json.loads((ROOT/'research/prior_data_manifest.json').read_text())
    entries = [e for e in manifest['files'] if e['group']=='core' and e['id']!='existing_faults']
    receipts = [restore(e) for e in entries]
    fetch('buffedlizard55-lab/GEMSDOE32','de3fa9270bdf3a9c7a8829a8df9e42f610ab681e',
          'docs/downloads/gemsdoe32-h33-h33-2-b2-20261004T220000Z-e5eb6e7e-zeros.tif',ROOT/'data/incumbent32.tif')
    (ROOT/'evidence/data_receipt.json').write_text(json.dumps(receipts,indent=2)+'\n')
