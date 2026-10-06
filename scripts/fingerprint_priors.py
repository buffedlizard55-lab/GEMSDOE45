"""Pixel-level comparison inventory, not just renamed-file checks. Downloaded priors are NOT model inputs."""
import concurrent.futures as cf
import hashlib
import json
from pathlib import Path
import numpy as np
import rasterio
from rasterio.windows import Window
from download_data import fetch
from core import ROOT

def fingerprint(item):
    repo,ref,t=item; path=ROOT/'data/prior_tmp'/f'{t["blob"]}.tif'
    row=dict(repo=repo,commit=ref,path=t['path'],blob=t['blob'])
    try:
        fetch('buffedlizard55-lab/'+repo,ref,t['path'],path)
        with rasterio.open(ROOT/'data/sample_submission.tif') as template,rasterio.open(path) as s:
            if s.count!=1 or s.shape!=template.shape or s.crs!=template.crs or s.transform!=template.transform:
                row['status']='different_grid_or_bands'; return row
            h=hashlib.sha256(); positive=0; invalid=0
            for y in range(0,s.height,128):
                win=Window(0,y,s.width,min(128,s.height-y))
                a=s.read(1,window=win).astype('float32')
                valid=np.isfinite(template.read(1,window=win)) & (template.read_masks(1,window=win)>0)
                invalid+=int(np.count_nonzero(valid & (~np.isfinite(a)|(a<0)|(a>1))))
                a[~valid]=0
                # Do NOT collapse within-footprint invalid values to legit predictions.
                h.update(a.tobytes()); positive+=int(np.count_nonzero(a>0))
            row.update(status='fingerprinted',canonical_pixel_sha256=h.hexdigest(),positive_pixels=positive,invalid_in_footprint=invalid)
    except Exception as e: row.update(status='unavailable',error=str(e)[:250])
    finally: path.unlink(missing_ok=True)
    return row

if __name__=='__main__':
    inventory=json.loads((ROOT/'evidence/prior_inventory.json').read_text())
    items={t['blob']:(r['repo'],r['commit'],t) for r in inventory['repositories'] for t in r.get('tifs',[])}
    with cf.ThreadPoolExecutor(max_workers=3) as pool:
        rows=[]
        for i,row in enumerate(pool.map(fingerprint,items.values())):
            rows.append(row)
            print(i+1,len(items),row['repo'],row['status'],flush=True)
    (ROOT/'evidence/prior_fingerprints.json').write_text(json.dumps(dict(scope='All unique Git blobs ending .tif in paths containing download, submission or artifact, from accessible prior main trees; not all releases or unlinked files',files=rows),indent=2)+'\n')
