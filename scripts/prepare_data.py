"""Build frozen feature ablation; template used only for geometry and finite footprint."""
import json
import numpy as np
import rasterio
from scipy import ndimage as ndi
from core import ROOT,phase_features,sha256

def main():
    data=ROOT/'data'
    with rasterio.open(data/'sample_submission.tif') as t:
        valid=np.isfinite(t.read(1)) & (t.read_masks(1)>0); shape=t.shape; transform=t.transform; crs=t.crs
    np.save(data/'valid.npy',valid)
    with rasterio.open(data/'labels.tif') as s:
        assert s.shape==shape and s.transform==transform and s.crs==crs
        labels=s.read(1)==1
    np.save(data/'truth.npy',labels & valid)
    stack=np.lib.format.open_memmap(data/'features.npy',mode='w+',dtype='float32',shape=(46,shape[0]*shape[1]))
    names=[]; metadata=[]
    with rasterio.open(data/'training_features.tif') as s:
        assert s.shape==shape and s.transform==transform and s.crs==crs and s.count==19
        for i in range(1,20):
            a=s.read(i,masked=True).filled(np.nan)
            a[~np.isfinite(a)]=np.nan
            stack[i-1]=a.ravel(); names.append(s.tags(i).get('band_name',f'band{i}'))
            metadata.append(dict(band=i,tags=s.tags(i),invalid=int((~np.isfinite(a)).sum())))
        derived=[]
        for k,band in enumerate([12,2,13]):
            a=s.read(band,masked=True).filled(np.nan)
            good=np.isfinite(a)
            if not good.any(): raise ValueError('Entire layer invalid')
            idx=ndi.distance_transform_edt(~good,return_distances=False,return_indices=True)
            a=a[tuple(idx)]; del idx
            # Constant offsets are irrelevant to derivatives, but subtracting the median
            # reduces discrete Gaussian derivative DC leakage. This uses no labels.
            a=a-np.median(a[good]); prior=None
            for j,sigma in enumerate([2,5]):
                odd,even,fraction,signed=phase_features(a,sigma)
                stack[19+k*4+j*2]=odd.ravel(); stack[20+k*4+j*2]=even.ravel()
                names.extend([f'b{band}_odd_s{sigma}',f'b{band}_even_s{sigma}'])
                stack[31+k*5+j*2]=fraction.ravel(); stack[32+k*5+j*2]=signed.ravel()
                derived.extend([f'b{band}_odd_fraction_s{sigma}',f'b{band}_signed_phase_s{sigma}'])
                if prior is not None:
                    stack[35+k*5]=(prior*signed).ravel(); derived.append(f'b{band}_signed_phase_agreement')
                prior=signed
                print('features',band,sigma,flush=True)
            stack.flush()
    receipt=dict(shape=list(shape),channels=names+derived,baseline_channels=31,candidate_channels=46,band_metadata=metadata,preregistration_sha256=sha256(ROOT/'research/hypotheses.md'),feature_bytes=(data/'features.npy').stat().st_size,template_warning='Mirror template has positive labels; values never used as predictions or model inputs.')
    (ROOT/'evidence/features.json').write_text(json.dumps(receipt,indent=2)+'\n')
if __name__=='__main__': main()
