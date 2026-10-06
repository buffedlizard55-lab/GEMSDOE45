"""Independent release checks; uses shipped artifacts, no raw training-data dependency."""
import hashlib
import json
from pathlib import Path
import zipfile
from urllib.parse import unquote,urlsplit
import numpy as np
import rasterio
from bs4 import BeautifulSoup
from core import ROOT,sha256

def verify():
    docs=ROOT/'docs'; audit=json.loads((docs/'downloads/submission-audit.json').read_text())
    tif=docs/'downloads'/audit['file']
    if sha256(tif)!=audit['sha256']: raise ValueError('Artifact SHA mismatch')
    with rasterio.open(tif) as s:
        p=s.read(1); mask=s.read_masks(1)>0
        assert s.count==1 and s.dtypes==('float32',)
        assert s.shape==tuple(audit['shape'])==(3730,3292)
        assert s.crs.to_epsg()==32611 and list(s.transform)==audit['transform']
        assert np.isfinite(p).all() and p.min()>=0 and p.max()<=1 and s.nodata is None
        assert not p[~mask].any() and mask.sum()==audit['footprint_pixels']
        assert hashlib.sha256(p.tobytes()).hexdigest()==audit['pixel_sha256']
        assert (p>0).sum()==audit['positive_pixels']
    with zipfile.ZipFile(tif.with_suffix('.zip')) as z:
        assert z.namelist()==[tif.name]
        assert hashlib.sha256(z.read(tif.name)).hexdigest()==audit['sha256']
    report=json.loads((ROOT/'evidence/holdout.json').read_text())
    assert report['preregistration_sha256']==sha256(ROOT/'research/hypotheses.md')==audit['preregistration_sha256']
    delta=[]
    for fold in report['folds']:
        for method in ['baseline','phase']:
            r=fold[method]['test']
            np.testing.assert_allclose(r['dti'],r['tp']/(r['tp']+.2*r['fp']+.8*r['fn']+1e-12),atol=1e-12)
            np.testing.assert_allclose(r['tp']+r['fn'],r['truth'],atol=1e-10)
            assert fold[method]['chosen_fraction']==max(fold[method]['calibration'],key=lambda x:x['dti'])['fraction']
        d=fold['phase']['test']['dti']-fold['baseline']['test']['dti']
        np.testing.assert_allclose(d,fold['paired_delta'],atol=1e-12); delta.append(d)
    np.testing.assert_allclose(np.mean(delta),report['mean_delta'],atol=1e-12)
    assert report['local_proxy_pass']==all(d>0 for d in delta)
    assert report['slot_approved'] is False and audit['slot_approved'] is False
    prior=json.loads((ROOT/'evidence/prior_fingerprints.json').read_text())['files']
    comparable=[r for r in prior if r['status']=='fingerprinted']
    assert len(comparable)>=1
    assert all(r['canonical_pixel_sha256']!=audit['pixel_sha256'] for r in comparable)
    uniqueness=json.loads((ROOT/'evidence/uniqueness.json').read_text())
    assert uniqueness['compared']==len(comparable) and uniqueness['identical_predictions']==0
    checked=0
    for path in docs.glob('*.html'):
        soup=BeautifulSoup(path.read_text(),'html.parser')
        assert soup.find('h1') and soup.find('main',id='main')
        for tag in soup.find_all(['a','img','script','link']):
            url=tag.get('href') or tag.get('src')
            if not url: continue
            parsed=urlsplit(url)
            if parsed.scheme or parsed.netloc: continue
            dest=path.parent/unquote(parsed.path) if parsed.path else path
            if not dest.exists(): raise ValueError(f'Broken local link {path.name}: {url}')
            if parsed.fragment and dest.suffix=='.html':
                target=BeautifulSoup(dest.read_text(),'html.parser')
                if not target.find(id=unquote(parsed.fragment)): raise ValueError(f'Broken anchor {url}')
            checked+=1
    for pth in (ROOT/'evidence').glob('*.json'):
        assert (docs/'data'/pth.name).read_bytes()==pth.read_bytes(),f'Stale docs data: {pth}'
    result=dict(format=True,zip=True,preregistration=True,holdout_arithmetic=True,unique_compared=len(comparable),local_links_checked=checked,slot_approved=False)
    print(json.dumps(result,indent=2)); return result
if __name__=='__main__': verify()
