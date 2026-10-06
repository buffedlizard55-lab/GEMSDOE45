"""Frozen four-quadrant ablation, then full-data inference. Never uploads to competition."""
import os
os.environ.setdefault('OMP_NUM_THREADS','2')
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
import datetime as dt
import json
import zipfile
import numpy as np
import rasterio
from sklearn.ensemble import HistGradientBoostingClassifier
from threadpoolctl import threadpool_limits
from core import ROOT,dti,pack,quadrant_masks,write_tif,sha256

FRACTIONS=[.004,.008,.016]

def fit(stack,truth,domain,nfeatures,seed):
    rng=np.random.default_rng(seed)
    y=truth.ravel(); mask=domain.ravel()
    pos=np.flatnonzero(mask & y); neg=np.flatnonzero(mask & ~y)
    if not len(pos) or not len(neg): raise ValueError('Training requires both classes')
    ids=np.concatenate([rng.choice(pos,min(30000,len(pos)),replace=False),rng.choice(neg,min(120000,len(neg)),replace=False)])
    rng.shuffle(ids)
    x=np.array(stack[:nfeatures,ids].T,order='C')
    model=HistGradientBoostingClassifier(max_iter=120,max_leaf_nodes=15,learning_rate=.08,l2_regularization=10,min_samples_leaf=50,early_stopping=False,random_state=seed)
    with threadpool_limits(limits=2): model.fit(x,y[ids])
    return model,dict(positive=int(y[ids].sum()),unlabelled=int((~y[ids]).sum()),seed=seed)

def predict(model,stack,mask,nfeatures):
    ids=np.flatnonzero(mask); out=np.zeros(mask.size,np.float32)
    with threadpool_limits(limits=2):
        for i in range(0,len(ids),100000):
            chunk=ids[i:i+100000]
            out[chunk]=model.predict_proba(np.array(stack[:nfeatures,chunk].T,order='C'))[:,1]
    return out.reshape(mask.shape)

def local_gate(deltas):
    return len(deltas)==4 and all(np.isfinite(d) and d>0 for d in deltas)

def main():
    data=ROOT/'data'; evidence=ROOT/'evidence'; evidence.mkdir(exist_ok=True)
    stack=np.load(data/'features.npy',mmap_mode='r'); valid=np.load(data/'valid.npy'); truth=np.load(data/'truth.npy')
    q,safe=quadrant_masks(valid.shape,valid)
    folds=[]
    for fold in range(4):
        cal=(fold+1)%4
        train=safe & (q!=fold) & (q!=cal); calibration=safe & (q==cal); test=safe & (q==fold)
        row=dict(fold=fold,calibration_quadrant=cal,train_pixels=int(train.sum()),calibration_pixels=int(calibration.sum()),test_pixels=int(test.sum()),test_truth=int(truth[test].sum()))
        for method,n in [('baseline',31),('phase',46)]:
            print('FIT',fold,method,flush=True)
            model,sample=fit(stack,truth,train,n,4506+fold)
            scores=predict(model,stack,calibration|test,n)
            trials=[]
            for fraction in FRACTIONS:
                pred=pack(scores,calibration,fraction)
                trials.append(dict(fraction=fraction,**dti(pred,truth,calibration)))
            chosen=max(trials,key=lambda x:x['dti'])['fraction']
            pred=pack(scores,test,chosen)
            row[method]=dict(chosen_fraction=chosen,calibration=trials,test=dti(pred,truth,test),samples=sample,emitted=int(pred.sum()))
            print('RESULT',fold,method,row[method]['test']['dti'],chosen,flush=True)
        row['paired_delta']=row['phase']['test']['dti']-row['baseline']['test']['dti']
        folds.append(row)
        (evidence/'holdout_partial.json').write_text(json.dumps(folds,indent=2)+'\n')
    deltas=[r['paired_delta'] for r in folds]
    report=dict(created_at=dt.datetime.now(dt.timezone.utc).isoformat(),protocol='research/hypotheses.md',preregistration_sha256=sha256(ROOT/'research/hypotheses.md'),folds=folds,mean_delta=float(np.mean(deltas)),delta_range=[min(deltas),max(deltas)],baseline_mean=float(np.mean([r['baseline']['test']['dti'] for r in folds])),phase_mean=float(np.mean([r['phase']['test']['dti'] for r in folds])),local_proxy_pass=local_gate(deltas),slot_approved=False,slot_reasons=['No independent off-catalogue labels','No comparable inherited holdout-best detector; H33 was pruned against visible catalogue'],score_type='Catalogue-transfer spatial proxy, NOT leaderboard score')
    (evidence/'holdout.json').write_text(json.dumps(report,indent=2)+'\n')
    (evidence/'holdout_partial.json').unlink(missing_ok=True)
    fraction=float(np.median([r['phase']['chosen_fraction'] for r in folds]))
    print('FULL FIT',fraction,flush=True)
    model,sample=fit(stack,truth,valid,46,4506)
    scores=predict(model,stack,valid,46)
    # This follows the official ALL-fault target: no catalogue exclusion or copied dots.
    pred=pack(scores,valid,fraction)
    pixel_id=__import__('hashlib').sha256(pred.tobytes()).hexdigest()[:12]
    name=f'gems45-phase-parity-20261006-{pixel_id}-research'
    dest=ROOT/'docs/downloads'/f'{name}.tif'; dest.parent.mkdir(parents=True,exist_ok=True)
    receipt=write_tif(dest,pred,data/'sample_submission.tif')
    receipt.update(method='H45-PHASE',fraction=fraction,training_samples=sample,unique_name='GEMS45-PHASE-'+pixel_id,note='H45 differential phase parity; 46 channels; nested geographic calibration; all-fault target. Research only; off-catalogue gate not passed.',preregistration_sha256=report['preregistration_sha256'])
    with rasterio.open(data/'incumbent32.tif') as s:
        old=s.read(1)
        receipt['incumbent32_pixel_differences']=int(np.count_nonzero(pred!=old))
        receipt['incumbent32_sha256']=sha256(data/'incumbent32.tif')
        receipt['incumbent32_jaccard']=float(np.count_nonzero((pred>0)&(old>0))/max(1,np.count_nonzero((pred>0)|(old>0))))
    if receipt['incumbent32_pixel_differences']==0: raise ValueError('Candidate duplicates incumbent')
    (ROOT/'docs/downloads/submission-audit.json').write_text(json.dumps(receipt,indent=2)+'\n')
    with zipfile.ZipFile(dest.with_suffix('.zip'),'w',zipfile.ZIP_DEFLATED) as z: z.write(dest,arcname=dest.name)
    np.save(data/'prediction.npy',pred)
    print('DONE',dest,receipt,flush=True)
if __name__=='__main__': main()
