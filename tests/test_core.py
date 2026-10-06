import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin
from core import dti,pack,phase_features,quadrant_masks,write_tif,validate_tif
from run_experiment import local_gate

def brute(p,g):
    truth=np.argwhere(g); pixels=np.argwhere(p>0)
    tp=sum(max([p[tuple(x)]*max(1-np.linalg.norm(x-t)/3,0) for x in pixels]+[0]) for t in truth)
    fp=sum(p[tuple(x)]*(1-max([max(1-np.linalg.norm(x-t)/3,0) for t in truth]+[0])) for x in pixels)
    fn=len(truth)-tp
    return tp/(tp+.2*fp+.8*fn+1e-12)

@pytest.mark.parametrize('binary',[False,True])
@pytest.mark.parametrize('seed',range(5))
def test_metric_against_independent_double_loop(binary,seed):
    rng=np.random.default_rng(seed); p=rng.random((9,11)).astype('float32'); g=rng.random(p.shape)<.1
    if binary: p=(p>.7).astype('float32')
    assert dti(p,g)['dti']==pytest.approx(brute(p,g),abs=1e-7)

def test_metric_empty_perfect_kernel():
    z=np.zeros((11,11),np.float32); g=z.astype(bool); g[5,5]=True
    assert dti(z,g)['dti']==0
    z[5,5]=1; assert dti(z,g)['dti']==pytest.approx(1)
    assert dti(z,np.zeros_like(g))['fp']==1
    z[5,5]=0; z[5,8]=1; assert dti(z,g)['tp']==0
    z[5,8]=0; z[5,6]=1; assert dti(z,g)['tp']==pytest.approx(2/3)

def test_metric_mask():
    p=np.zeros((12,12),np.float32); g=p.astype(bool); g[6,6]=True; p[6,6]=1
    mask=np.ones_like(g); mask[6,6]=False
    assert dti(p,g,mask)['truth']==0

@pytest.mark.parametrize('bad',[np.nan,np.inf,-1,1.1])
def test_invalid_metric_rejected(bad):
    p=np.full((4,4),bad,np.float32)
    with pytest.raises(ValueError): dti(p,np.ones_like(p,bool))

def test_pack_spacing_budget_determinism():
    rng=np.random.default_rng(17); a=rng.random((50,50)).astype('float32'); v=np.ones_like(a,bool); v[:3]=False
    p=pack(a,v,.016); points=np.argwhere(p)
    assert p.sum()==round(v.sum()*.016)
    for i in range(len(points)):
        for j in range(i): assert np.linalg.norm(points[i]-points[j])>=2.8
    np.testing.assert_array_equal(p,pack(a,v,.016)); assert not p[~v].any()

def test_phase_bounds_polarity():
    rng=np.random.default_rng(45); a=rng.normal(size=(100,100)).astype('float32')
    odd,even,f,signed=phase_features(a,2)
    other=phase_features(-a,2)
    assert np.all((f>=0)&(f<=1)) and np.all(np.abs(signed)<=1.00001)
    np.testing.assert_allclose(f,other[2]); np.testing.assert_allclose(signed,-other[3])

def test_step_is_more_odd_than_symmetric_ridge_at_center():
    x=np.arange(201)-100
    step=np.tile((x>=0).astype('float32'),(101,1))
    ridge=np.tile(np.exp(-x*x/50).astype('float32'),(101,1))
    assert phase_features(step,2)[2][50,100]>phase_features(ridge,2)[2][50,100]

def test_spatial_embargo_and_gate():
    q,safe=quadrant_masks((300,300),np.ones((300,300),bool))
    assert not safe[121:180].any(); assert not safe[:,121:180].any()
    assert local_gate([.1,.1,.1,.1]); assert not local_gate([.1,.1,-.1,.1])
    assert not local_gate([np.nan,.1,.1,.1]); assert not local_gate([])

def template(tmp_path):
    path=tmp_path/'template.tif'; a=np.zeros((20,30),np.float32); a[:3]=np.nan
    with rasterio.open(path,'w',driver='GTiff',width=30,height=20,count=1,dtype='float32',crs='EPSG:32611',transform=from_origin(243350,4508550,100,100),nodata=np.nan) as s: s.write(a,1)
    return path

def test_tif_roundtrip_and_internal_mask(tmp_path):
    t=template(tmp_path); out=tmp_path/'out.tif'; p=np.ones((20,30),np.float32)
    receipt=write_tif(out,p,t)
    assert receipt['local_format_pass'] and receipt['nonfinite']==0 and receipt['positive_pixels']==510
    with rasterio.open(out) as s:
        assert s.read(1,masked=True).mask[:3].all()
        assert not s.read(1)[:3].any()
    assert not out.with_suffix('.tif.msk').exists()

@pytest.mark.parametrize('bad',[np.nan,np.inf,-1,1.1])
def test_writer_rejects_not_clips(tmp_path,bad):
    with pytest.raises(ValueError): write_tif(tmp_path/'bad.tif',np.full((20,30),bad,np.float32),template(tmp_path))

def test_grid_mismatch(tmp_path):
    t=template(tmp_path); out=tmp_path/'out.tif'
    write_tif(out,np.zeros((20,30),np.float32),t)
    with rasterio.open(out,'r+') as s: s.transform=from_origin(1,1,100,100)
    with pytest.raises(ValueError,match='Grid'): validate_tif(out,t)
