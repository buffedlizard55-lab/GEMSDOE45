"""Original H45 numerical primitives. DTI implements the published formula, not its server."""
import hashlib
from pathlib import Path
import numpy as np
from scipy import ndimage as ndi
import rasterio

ROOT = Path(__file__).resolve().parents[1]

def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for chunk in iter(lambda:f.read(1048576),b''): h.update(chunk)
    return h.hexdigest()

def dti(pred, truth, mask=None):
    p=np.asarray(pred,dtype=np.float32)
    g=np.asarray(truth,dtype=bool)
    if p.shape!=g.shape or p.ndim!=2: raise ValueError('Matching 2D grids required')
    if not np.isfinite(p).all() or np.any((p<0)|(p>1)): raise ValueError('Invalid probability')
    if mask is not None:
        if mask.shape!=p.shape: raise ValueError('Mask shape')
        p=np.where(mask,p,0); g=g & mask
    if not g.any():
        return dict(dti=0.0,tp=0.0,fp=float(p.sum(dtype=np.float64)),fn=0.0,truth=0)
    distance=ndi.distance_transform_edt(~g)
    fp=float(np.sum(p*np.minimum(distance/3,1),dtype=np.float64))
    if np.all((p==0)|(p==1)):
        cover=np.maximum(1-ndi.distance_transform_edt(p==0)/3,0) if p.any() else np.zeros_like(p)
    else:
        cover=np.zeros_like(p)
        h,w=p.shape
        for dy in range(-2,3):
            for dx in range(-2,3):
                weight=max(1-np.hypot(dx,dy)/3,0)
                ys=slice(max(0,dy),min(h,h+dy)); xs=slice(max(0,dx),min(w,w+dx))
                yt=slice(max(0,-dy),min(h,h-dy)); xt=slice(max(0,-dx),min(w,w-dx))
                np.maximum(cover[ys,xs],p[yt,xt]*weight,out=cover[ys,xs])
    tp=float(cover[g].sum(dtype=np.float64)); fn=float(g.sum()-tp)
    return dict(dti=tp/(tp+.2*fp+.8*fn+1e-12),tp=tp,fp=fp,fn=fn,truth=int(g.sum()))

def phase_features(a, sigma):
    """Odd/even differential energy, bounded signed curvature, scale normalized."""
    gx=ndi.gaussian_filter(a,sigma,order=(0,1))*sigma
    gy=ndi.gaussian_filter(a,sigma,order=(1,0))*sigma
    xx=ndi.gaussian_filter(a,sigma,order=(0,2))*sigma**2
    yy=ndi.gaussian_filter(a,sigma,order=(2,0))*sigma**2
    xy=ndi.gaussian_filter(a,sigma,order=(1,1))*sigma**2
    odd=np.hypot(gx,gy)
    even=np.sqrt(xx*xx+yy*yy+2*xy*xy)
    denom=odd+even+np.float32(1e-6)
    return odd,even,odd/denom,(xx+yy)/(np.sqrt(2)*denom)

def quadrant_masks(shape, valid, gap=30):
    h,w=shape; yy,xx=np.ogrid[:h,:w]
    q=((yy>=h//2)*2+(xx>=w//2)).astype(np.int8)
    safe=valid & (np.abs(yy-h//2)>=gap) & (np.abs(xx-w//2)>=gap)
    return q,safe

def pack(score, valid, fraction):
    """Deterministic greedy d=2.8 separation. No coordinate-lattice emission."""
    if not 0<fraction<=1: raise ValueError('Fraction outside (0,1]')
    if score.shape!=valid.shape or not np.isfinite(score[valid]).all(): raise ValueError('Score grid')
    h,w=score.shape; ids=np.flatnonzero(valid)
    order=ids[np.argsort(-score.ravel()[ids],kind='stable')]
    budget=int(round(fraction*len(ids)))
    occupied=np.zeros((h,w),bool); result=np.zeros((h,w),np.float32)
    offsets=[(dy,dx) for dy in range(-2,3) for dx in range(-2,3) if dy*dy+dx*dx<2.8**2]
    count=0
    for idx in order:
        if count>=budget: break
        y,x=divmod(int(idx),w)
        if occupied[y,x]: continue
        result[y,x]=1; count+=1
        for dy,dx in offsets:
            if 0<=y+dy<h and 0<=x+dx<w: occupied[y+dy,x+dx]=True
    return result

def write_tif(path, pred, template):
    with rasterio.open(template) as t:
        profile=t.profile.copy(); valid=np.isfinite(t.read(1)) & (t.read_masks(1)>0)
        if pred.shape!=t.shape: raise ValueError('Prediction shape')
    if not np.isfinite(pred).all() or np.any((pred<0)|(pred>1)): raise ValueError('Invalid predictions; refusing to silently repair')
    profile.update(driver='GTiff',dtype='float32',count=1,nodata=None,compress='deflate',predictor=3,tiled=True,blockxsize=256,blockysize=256)
    with rasterio.Env(GDAL_TIFF_INTERNAL_MASK=True):
        with rasterio.open(path,'w',**profile) as out:
            out.write(np.where(valid,pred,0).astype('float32'),1)
            out.write_mask(valid.astype('uint8')*255)
            out.update_tags(method='H45-PHASE',status='RESEARCH_ONLY_NOT_SLOT_APPROVED')
    return validate_tif(path,template)

def validate_tif(path,template):
    with rasterio.open(template) as t,rasterio.open(path) as s:
        if s.count!=1 or s.dtypes!=('float32',): raise ValueError('Requires one float32 band')
        if s.crs!=t.crs or s.shape!=t.shape or s.transform!=t.transform: raise ValueError('Grid mismatch')
        p=s.read(1); valid=np.isfinite(t.read(1)) & (t.read_masks(1)>0)
        if not np.isfinite(p).all() or np.any((p<0)|(p>1)): raise ValueError('Outside [0,1] or nonfinite')
        if s.nodata is not None: raise ValueError('Must not carry nodata sentinel')
        if np.any(p[~valid]!=0): raise ValueError('Nonzero outside footprint')
        if not np.array_equal(s.read_masks(1)>0,valid): raise ValueError('Footprint mask mismatch')
        return dict(file=Path(path).name,sha256=sha256(path),pixel_sha256=hashlib.sha256(p.tobytes()).hexdigest(),shape=list(s.shape),crs=s.crs.to_string(),transform=list(s.transform),count=s.count,dtype=s.dtypes[0],min=float(p.min()),max=float(p.max()),nonfinite=int((~np.isfinite(p)).sum()),positive_pixels=int((p>0).sum()),footprint_pixels=int(valid.sum()),nodata=s.nodata,internal_mask=True,bytes=Path(path).stat().st_size,local_format_pass=True,portal_acceptance='unverified',slot_approved=False)
