"""Non-authenticated source-health feed; failures are surfaced, never guessed into scores."""
import concurrent.futures as cf
import datetime as dt
import hashlib
import json
from pathlib import Path
import re
import requests
from bs4 import BeautifulSoup
ROOT=Path(__file__).resolve().parents[1]

def leaderboard_top(html):
    soup=BeautifulSoup(html,'html.parser')
    if not soup.find('h1') or 'leaderboard' not in soup.find('h1').get_text().lower():
        raise ValueError('Not a leaderboard page')
    for table in soup.find_all('table'):
        if 'Tversky' not in table.get_text(): continue
        for row in table.find_all('tr'):
            cells=row.find_all('td')
            if not cells: continue
            for cell in cells:
                value=cell.get_text(' ',strip=True)
                if re.fullmatch(r'0\.\d{4}',value): return float(value)
    raise ValueError('Leaderboard schema changed or no score rows')

def check(source):
    row=dict(id=source['id'],url=source['url'],checked_at=dt.datetime.now(dt.timezone.utc).isoformat())
    try:
        r=requests.get(source['url'],timeout=(10,30),headers={'User-Agent':'GEMSDOE45-source-health/1.0'})
        row.update(http_status=r.status_code,final_url=r.url)
        r.raise_for_status()
        if '/login/' in r.url:
            row['status']='login_required'; return row
        row.update(status='reachable',sha256=hashlib.sha256(r.content).hexdigest())
        if source['id']=='leaderboard':
            row['top_score']=leaderboard_top(r.text); row['status']='parsed'
    except Exception as e: row.update(status='error',error=str(e)[:350])
    return row

if __name__=='__main__':
    sources=json.loads((ROOT/'evidence/official_sources.json').read_text())['sources']
    with cf.ThreadPoolExecutor(max_workers=4) as pool: rows=list(pool.map(check,sources))
    target=ROOT/'docs/data/source_health.json'; target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(dict(checked_at=dt.datetime.now(dt.timezone.utc).isoformat(),sources=rows,note='Reachability is not independent verification of every source claim. Failed checks do not replace the dated research snapshot.'),indent=2)+'\n')
    print(json.dumps({r['id']:r['status'] for r in rows},indent=2))
