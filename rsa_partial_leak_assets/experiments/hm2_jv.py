# Pure-Python/fpylll implementation of jvdsn Herrmann-May multivariate for this challenge.
import sys, time, math, itertools, argparse
import numpy as np
from fpylll import IntegerMatrix, LLL, BKZ
sys.path.append('/mnt/data')
from rsa_partial_lib import N, ct, e, mask, leak, decrypt_from_p

XS,XL,YS,YL = 265,155,600,230
X,Y = 1<<XL, 1<<YL
A0,B0 = 1<<XS, 1<<YS
base = leak & ~(((1<<XL)-1)<<XS) & ~(((1<<YL)-1)<<YS)

def poly_mul(p,q):
    r={}
    for (i,j),a in p.items():
        for (k,l),b in q.items():
            r[(i+k,j+l)] = r.get((i+k,j+l),0) + a*b
    return {m:c for m,c in r.items() if c}

def poly_pow(p,k):
    r={(0,0):1}
    b=p
    while k:
        if k&1: r=poly_mul(r,b)
        k//=2
        if k: b=poly_mul(b,b)
    return r

def poly_scale_shift(p, coeff, sx, sy):
    if coeff==0: return {}
    return {(i+sx,j+sy): c*coeff for (i,j),c in p.items()}

def build(C,m,t,lead='x'):
    # Build jvdsn shared.small_roots.herrmann_may_multivariate shifts for two variables.
    if lead=='x':
        inv=pow(A0,-1,N)
        f={(1,0):1, (0,1):(B0*inv)%N, (0,0):(C*inv)%N}  # monic in x
        bounds={(1,0):X,(0,1):Y}
        nonlead='y'
    else:
        inv=pow(B0,-1,N)
        f={(0,1):1, (1,0):(A0*inv)%N, (0,0):(C*inv)%N}  # monic in y
        bounds={(1,0):X,(0,1):Y}
        nonlead='x'
    fpows=[{(0,0):1}]
    for k in range(1,m+1):
        fpows.append(poly_mul(fpows[-1], f))       # NO coefficient reduction after f % N
    shifts=[]
    meta=[]
    for k in range(m+1):
        coeff = pow(N, max(t-k,0))
        # _get_shifts starts at j=1, so only powers of the non-leading variable, exponent < m+1-k
        for s in range(m+1-k):
            sx,sy = ((0,s) if nonlead=='y' else (s,0))
            shifts.append(poly_scale_shift(fpows[k], coeff, sx, sy))
            meta.append((k,s))
    mons=sorted(set().union(*(p.keys() for p in shifts)), key=lambda z:(z[0]+z[1],z[0],z[1]))
    col={mon:i for i,mon in enumerate(mons)}
    scales=[(X**i)*(Y**j) for i,j in mons]
    M=IntegerMatrix(len(shifts), len(mons))
    for r,p in enumerate(shifts):
        for mon,c in p.items():
            M[r,col[mon]] = int(c * scales[col[mon]])
    return M, mons, scales

def unscale_row(M,r,mons,scales):
    d={}; norm2=0; w=0
    for j,(mon,s) in enumerate(zip(mons,scales)):
        v=int(M[r,j])
        if v:
            norm2 += v*v; w += 1
            if v % s: return None, norm2, w
            d[mon] = v//s
    return d,norm2,w

def poly_eval_mod(p,a,b,mod):
    # p small prime, a,b residues
    max_i=max((i for i,j in p), default=0); max_j=max((j for i,j in p), default=0)
    xp=[1]*(max_i+1); yp=[1]*(max_j+1)
    for i in range(1,max_i+1): xp[i]=xp[i-1]*a%mod
    for j in range(1,max_j+1): yp[j]=yp[j-1]*b%mod
    return sum((c%mod)*xp[i]*yp[j] for (i,j),c in p.items())%mod

_pow_cache={}
def _grid(pr,mi,mj):
    key=(pr,mi,mj)
    if key in _pow_cache: return _pow_cache[key]
    a=np.repeat(np.arange(pr,dtype=np.int64), pr)
    b=np.tile(np.arange(pr,dtype=np.int64), pr)
    xp=[np.ones(pr*pr,dtype=np.int64)]
    yp=[np.ones(pr*pr,dtype=np.int64)]
    for i in range(1,mi+1): xp.append((xp[-1]*a)%pr)
    for j in range(1,mj+1): yp.append((yp[-1]*b)%pr)
    _pow_cache[key]=(a,b,xp,yp)
    return a,b,xp,yp

def roots_mod(polys,pr):
    if not polys: return []
    mi=max(i for p in polys for i,j in p); mj=max(j for p in polys for i,j in p)
    a,b,xp,yp=_grid(pr,mi,mj)
    ok=np.ones(pr*pr,dtype=bool)
    for p in polys:
        idx=np.nonzero(ok)[0]
        if len(idx)==0: return []
        s=np.zeros(len(idx),dtype=np.int64)
        for (i,j),c in p.items():
            cc=c%pr
            if cc: s=(s + cc*((xp[i][idx]*yp[j][idx])%pr))%pr
        ok[idx] &= (s==0)
        if ok.sum()>50 and p is polys[-1]:
            # caller will abort large root sets; avoid list building cost for huge sets
            pass
    idx=np.nonzero(ok)[0]
    if len(idx)>10000: return [(0,0)]*10001
    return list(zip(a[idx].tolist(), b[idx].tolist()))

def crt_states(states,mod,R,pr,cap=1000):
    inv=pow(mod,-1,pr); out=[]
    for x,y in states:
        xm=x%pr; ym=y%pr
        for a,b in R:
            nx=x + mod*(((a-xm)*inv)%pr)
            ny=y + mod*(((b-ym)*inv)%pr)
            # range pruning once mod gets large-ish
            if nx < X and ny < Y:
                out.append((nx,ny))
                if len(out)>cap: return [],0
    return out, mod*pr

PRIMES=[257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541]

def recover_candidates(polys, max_polys=6, verbose=False):
    # Try increasing subsets; too many bad polynomials can kill the true root.
    for r in range(2, min(max_polys,len(polys))+1):
        P=polys[:r]
        states=[(0,0)]; mod=1
        good=True
        for pr in PRIMES:
            R=roots_mod(P,pr)
            if verbose: print('  r',r,'prime',pr,'roots',len(R),'states',len(states),'modbits',mod.bit_length(), flush=True)
            if not R or len(R)>50:
                good=False; break
            states,mod=crt_states(states,mod,R,pr)
            if not states:
                good=False; break
            if mod > X and mod > Y:
                return states
        if good and mod > X and mod > Y:
            return states
    return []

def try_one(cid,m,t,lead,delta=0.99,bkz=0, verbose=False):
    lo,hi=cid&15,cid>>4
    C=base | (lo<<150) | (hi<<920)
    M,mons,scales=build(C,m,t,lead)
    t0=time.time()
    LLL.reduction(M, delta=delta, method='proved', float_type='mpfr', precision=256)
    if bkz:
        BKZ.reduction(M, BKZ.Param(block_size=bkz))
    dt=time.time()-t0
    rows=[]
    for r in range(M.nrows):
        p,n2,w=unscale_row(M,r,mons,scales)
        if not p or len(p)<=1: continue
        # jvdsn filters by norm^2 * w < N^2
        if w and n2*w < N*N:
            rows.append((n2.bit_length(), w, p))
    rows.sort(key=lambda z:z[0])
    print(f'cid={cid:03d} lo={lo:x} hi={hi:x} lead={lead} m={m} t={t} dim={M.nrows}x{M.ncols} LLL={dt:.2f}s good={len(rows)} bits={[b for b,_,_ in rows[:5]]}', flush=True)
    if not rows: return False
    # remove duplicates
    uniq=[]; seen=set()
    for _,_,p in rows:
        key=tuple(sorted(p.items()))
        if key not in seen:
            uniq.append(p); seen.add(key)
    # Try CRT on selected polynomials. Also try sliding windows because first rows can be dependent/bad.
    candidate_sets=[]
    candidate_sets.append(uniq[:12])
    for off in range(0,min(4,len(uniq))): candidate_sets.append(uniq[off:off+12])
    for polys in candidate_sets:
        cands=recover_candidates(polys, max_polys=min(6,len(polys)), verbose=verbose)
        if verbose: print('cands',len(cands), flush=True)
        for x,y in cands:
            p=C+(x<<XS)+(y<<YS)
            if N%p==0 and (p&mask)==leak:
                print('FOUND')
                print('p =',hex(p))
                print('q =',hex(N//p))
                print('m =',decrypt_from_p(p))
                print('m.hex =',decrypt_from_p(p).hex())
                return True
    return False

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('-m',type=int,default=8); ap.add_argument('-t',type=int,default=None)
    ap.add_argument('--lead',choices=['x','y','both'],default='both')
    ap.add_argument('--a',type=int,default=0); ap.add_argument('--b',type=int,default=256)
    ap.add_argument('--verbose',action='store_true')
    args=ap.parse_args(); t=args.t if args.t is not None else int((1-math.sqrt(0.5))*args.m)
    leads=['x','y'] if args.lead=='both' else [args.lead]
    for cid in range(args.a,args.b):
        for lead in leads:
            if try_one(cid,args.m,t,lead,verbose=args.verbose): sys.exit(0)
    print('not found')
