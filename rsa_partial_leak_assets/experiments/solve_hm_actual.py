import sys, time, itertools, sympy as sp
sys.path.append('/mnt/data')
from constants import N,e,ct,MASK,LEAK
from hm_multivar import hm_shifts, create_lattice, reduce_lattice, reconstruct_row, eval_poly, row_norm_bits
from root_recover import recover_crt

A=1<<265; B=1<<600; X=1<<155; Y=1<<230
WX=((1<<155)-1)<<265; WY=((1<<230)-1)<<600
PRIMES=list(sp.primerange(1009,10000))[:80]

def dec(p):
    q=N//p
    phi=(p-1)*(q-1)
    d=pow(e,-1,phi)
    m=pow(ct,d,N)
    return m.to_bytes((m.bit_length()+7)//8,'big')

def try_candidate(lo,hi,m=8,t=3,root=True,method='lll'):
    base=LEAK | (lo<<150) | (hi<<920)
    p0=base & ~WX & ~WY
    invA=pow(A,-1,N)
    C=(B*invA)%N; D=(p0*invA)%N
    f={(1,0):1,(0,1):C,(0,0):D}
    shifts=hm_shifts(f,N,m,t)
    M,mons=create_lattice(shifts,X,Y)
    t0=time.time()
    reduce_lattice(M,method)
    print(f'LLL lo={lo:x} hi={hi:x} m={m} t={t} dim={M.nrows} time={time.time()-t0:.1f}s', flush=True)
    polys=[]
    # take first rows with low norm. For actual no true eval; include nonzero rows with not too huge norm.
    for r in range(M.nrows):
        pol=reconstruct_row(M,r,mons,X,Y)
        if pol:
            polys.append(pol)
    # root recovery using first 12? We don't know. Try all but root_recover only first few pairs.
    if not root: return None
    outs=recover_crt(polys, X, Y, PRIMES, verbose=True)
    print('outs',len(outs),outs[:2], flush=True)
    for xx,yy in outs:
        p=p0 + A*xx + B*yy
        if N%p==0 and (p&MASK)==LEAK:
            print('FOUND',lo,hi)
            print('p=',hex(p))
            print('q=',hex(N//p))
            print('m=',dec(p))
            return p
    return None

if __name__=='__main__':
    lo=int(sys.argv[1],0) if len(sys.argv)>1 else 0
    hi=int(sys.argv[2],0) if len(sys.argv)>2 else 0
    m=int(sys.argv[3]) if len(sys.argv)>3 else 8
    t=int(sys.argv[4]) if len(sys.argv)>4 else 3
    try_candidate(lo,hi,m,t)
