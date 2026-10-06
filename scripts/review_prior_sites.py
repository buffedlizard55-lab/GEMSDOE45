"""Read the published landing-page source of each supplied prior repository via GitHub.
A source scan is not a complete audit of every historical implementation.
"""
import concurrent.futures as cf
import datetime as dt
import json
import subprocess
from pathlib import Path
from bs4 import BeautifulSoup
ROOT = Path(__file__).resolve().parents[1]
REPOS = ['GEMSDOE','GEMSDOE2','GEMSDOE3','GEMSDOE4','5GEMSDOE','6GEMSDOE','7GEMSDOE','8GEMSDOE','GEMSDOE9','GEMSDOE10','11GEMSDOE','12GEMSDOE','13GEMSDOE','14GEMSDOE','15GEMSDOE','16GEMSDOE','17GEMSDOE','18GEMSDOE','19GEMSDOE','20GEMSDOE','GEMSDOE21'] + ['GEMSDOE'+str(i) for i in range(22,42)] + ['42GEMSDOE','43GEMSDOE','44GEMSDOE','GEMSDOE42','GEMSDOE43','GEMSDOE44']
def api(path):
    r = subprocess.run(['gh','api',path],capture_output=True,text=True,timeout=60)
    if r.returncode: raise RuntimeError(r.stderr.strip()[:180])
    return json.loads(r.stdout)
def scan(repo):
    base='repos/buffedlizard55-lab/'+repo
    try:
        tree=api(base+'/git/trees/main?recursive=1')
        paths={x['path']:x for x in tree['tree']}
        path=next((p for p in ['docs/index.html','index.html','docs/index.md','index.md','README.md'] if p in paths),None)
        if not path: raise ValueError('No landing source')
        r=subprocess.run(['gh','api',base+'/contents/'+path+'?ref='+tree['sha'],'-H','Accept: application/vnd.github.raw'],capture_output=True,text=True,check=True,timeout=60)
        soup=BeautifulSoup(r.stdout,'html.parser')
        for e in soup(['style','script']): e.decompose()
        text=soup.get_text(' ',strip=True)
        dest=ROOT/'research/prior_summaries'/f'{repo}.txt'; dest.parent.mkdir(parents=True,exist_ok=True)
        dest.write_text(text[:24000])
        tifs=[{'path':p,'bytes':x.get('size'), 'blob':x['sha']} for p,x in paths.items() if p.lower().endswith('.tif') and any(k in p for k in ['download','submission','artifact'])]
        return dict(repo=repo,status='reviewed landing source',commit=tree['sha'],path=path,url=f'https://github.com/buffedlizard55-lab/{repo}/blob/{tree["sha"]}/{path}',tifs=tifs,characters=len(text),truncated=len(text)>24000)
    except Exception as e: return dict(repo=repo,status='unavailable',error=str(e))
if __name__=='__main__':
    with cf.ThreadPoolExecutor(max_workers=6) as pool: results=list(pool.map(scan,REPOS))
    (ROOT/'evidence/prior_inventory.json').write_text(json.dumps(dict(checked_at=dt.datetime.now(dt.timezone.utc).isoformat(),scope='Landing sources; not exhaustive source-code novelty proof',repositories=results),indent=2)+'\n')
    for x in results: print(x['repo'],x['status'],len(x.get('tifs',[])),flush=True)
