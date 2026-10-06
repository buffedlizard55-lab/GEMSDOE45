import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import pytest
from refresh_sources import leaderboard_top
from verify_repo import verify

def test_published_release():
    assert verify()['format']

def test_feed_parses_only_named_leaderboard_table():
    fixture='<h1>Leaderboard</h1><table><tr><th>DW-Tversky</th></tr><tr><td>#1</td><td>23 submissions</td><td>0.3345</td></tr><tr><td>0.3262</td></tr></table>'
    assert leaderboard_top(fixture)==.3345

@pytest.mark.parametrize('html',['<h1>Login</h1><p>0.9999</p>','<h1>Leaderboard</h1><p>0.9999</p>','<h1>Leaderboard</h1><table><tr><td>0.9999</td></tr></table>'])
def test_feed_fails_closed_on_schema_change(html):
    with pytest.raises(ValueError): leaderboard_top(html)
