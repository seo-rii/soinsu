import sys, time, math, itertools
sys.path.append('/mnt/data')
from rsa_partial_lib import N,ct,mask,leak,e,x_start,x_len,y_start,y_len,X,Y,C_for,decrypt_from_p
from fpylll import IntegerMatrix, LLL, BKZ

Acoef = 1<<x_start
Bcoef = 1<<y_start

# Bivariate polynomial dict (ix,iy)->int over ZZ. For normalized basepolys coefficients are residues 0..N-1.
def p_add(p,q):
    r=p.copy()
    for m,c in q.items():
        if c:
            r[m]=r.get(m,0)+c
            if r[m]==0: del r[m]
    return r

def p_mul(p,q):
    r={}
    for (i,j),a in p.items():
        if not a: continue
        for (k,l),b in q.items():
            if b:
                m=(i+k,j+l)
                r[m]=r.get(m,0)+a*b
    return {m:c for m,c in r.items() if c}

def p_pow(p,k):
    r={(0,0):1}
    b=p
    while k:
        if k&1: r=p_mul(r,b)
        k//=2
        if k: b=p_mul(b,b)
    return r

def p_monom(sx,sy,coef=1):
    return {(sx,sy):coef} if coef else {}

def p_scale(p,c):
    if c==1: return p.copy()
    if c==0: return {}
    return {m:v*c for m,v in p.items() if v*c}

def p_shift(p,sx,sy,coef=1):
    if coef==0: return {}
    return {(i+sx,j+sy):v*coef for (i,j),v in p.items() if v*coef}

def norm_base(C, varidx):
    # varidx 0 normalize coefficient of x to 1, varidx 1 normalize y to 1, modulo N representatives 0..N-1.
    if varidx==0:
        inv=pow(Acoef,-1,N)
        return {(0,0):(C*inv)%N, (1,0):1, (0,1):(Bcoef*inv)%N}
    else:
        inv=pow(Bcoef,-1,N)
        return {(0,0):(C*inv)%N, (1,0):(Acoef*inv)%N, (0,1):1}

def build_hm_linear_rows(C,m,t,beta=0.499, include_dups=True):
    shifts=[]
    for varidx in (0,1):
        base=norm_base(C,varidx)
        # precompute powers
        powers=[{(0,0):1}]
        for k in range(1,m+1):
            powers.append(p_mul(powers[-1],base))
        for k in range(m+1):
            # j total degree among variables except varidx, from 0..m-k
            for j in range(m-k+1):
                sx,sy=(0,j) if varidx==0 else (j,0)
                Npow=pow(N, max(t-k,0))
                poly=p_shift(powers[k], sx, sy, Npow)
                shifts.append((varidx,k,j,poly))
    # Optionally remove exact duplicate polynomials? Could happen for k=0, j=0 for both vars = N^t
    if not include_dups:
        seen=set(); new=[]
        for varidx,k,j,p in shifts:
            key=tuple(sorted(p.items()))
            if key not in seen:
                seen.add(key); new.append((varidx,k,j,p))
        shifts=new
    mons=sorted(set().union(*(set(p.keys()) for _,_,_,p in shifts)), key=lambda z:(z[0]+z[1],z[0],z[1]))
    col={m:i for i,m in enumerate(mons)}
    scales=[(1<< (x_len*i + y_len*j)) for i,j in mons]
    rows=[]; meta=[]
    for varidx,k,j,p in shifts:
        row=[0]*len(mons)
        for mon,c in p.items():
            row[col[mon]]=int(c)*scales[col[mon]]
        rows.append(row); meta.append((varidx,k,j))
    return rows, mons, scales, meta

def reduce_rows(rows, delta=0.99, precision=256, bkz_block=None):
    M=len(rows); n=len(rows[0])
    A=IntegerMatrix(M,n)
    for i,row in enumerate(rows):
        for j,v in enumerate(row): A[i,j]=int(v)
    t0=time.time()
    LLL.reduction(A, delta=delta, method='proved', float_type='mpfr', precision=precision)
    if bkz_block:
        par=BKZ.Param(block_size=bkz_block, strategies=BKZ.DEFAULT_STRATEGY, max_loops=1)
        BKZ.reduction(A, par)
    dt=time.time()-t0
    out=[[int(A[i,j]) for j in range(n)] for i in range(M)]
    return out, dt

def unscale_row(row, mons, scales):
    d={}
    for val,mon,sc in zip(row,mons,scales):
        if val:
            if val % sc != 0:
                return None
            d[mon]=val//sc
    return d

def row_l1_bits(row):
    s=0
    for v in row: s += abs(v)
    return s.bit_length() if s else 0

def row_max_bits(row):
    vals=[abs(v).bit_length() for v in row if v]
    return max(vals) if vals else 0

def get_polys(C,m,t,beta=0.499,limit=40,verbose=True,bkz_block=None,filter_bound=False, include_dups=False):
    rows,mons,scales,meta=build_hm_linear_rows(C,m,t,beta,include_dups=include_dups)
    R,dt=reduce_rows(rows,precision=256,bkz_block=bkz_block)
    if verbose:
        print(f'HM rows {len(rows)} cols {len(mons)} m={m} t={t} LLL {dt:.2f}s', flush=True)
        print('first maxbits', [row_max_bits(r) for r in R[:12]], flush=True)
        print('first l1bits ', [row_l1_bits(r) for r in R[:12]], flush=True)
        print('howgrave bit approx', int(math.log2(N)*beta*t), flush=True)
    hbits=math.log2(N)*beta*t
    polys=[]
    for r in R:
        if row_max_bits(r)==0: continue
        if filter_bound and row_l1_bits(r) >= hbits: continue
        c=unscale_row(r,mons,scales)
        if c is None: continue
        # skip constants/monomials? keep
        mb=max(abs(v).bit_length() for v in c.values()) if c else 0
        deg=max((i+j for i,j in c), default=0)
        terms=len(c)
        polys.append((row_l1_bits(r), row_max_bits(r), mb, deg, terms, c))
    polys.sort(key=lambda z:z[0])
    return polys[:limit]

# Fast modular resultant/root via brute x,y for small primes. For low degree and p~1000 ok.
def poly_y_at_x(coeffs,xv,p):
    maxj=max((j for i,j in coeffs.keys()), default=0)
    arr=[0]*(maxj+1)
    max_i=max((i for i,j in coeffs.keys()), default=0)
    xp=[1]*(max_i+1)
    for i in range(1,max_i+1): xp[i]=xp[i-1]*xv%p
    for (i,j),c in coeffs.items():
        arr[j]=(arr[j]+(c%p)*xp[i])%p
    while arr and arr[-1]==0: arr.pop()
    return arr

def poly_trim(a):
    while a and a[-1]==0: a.pop()
    return a

def poly_mod_rem(a,b,p):
    a=a[:]; b=poly_trim([x%p for x in b])
    if not b: raise ZeroDivisionError
    inv=pow(b[-1],-1,p)
    while len(a)>=len(b) and a:
        coef=a[-1]*inv%p
        if coef:
            off=len(a)-len(b)
            for j,bj in enumerate(b):
                a[off+j]=(a[off+j]-coef*bj)%p
        poly_trim(a)
    return a

def poly_mod_gcd(a,b,p):
    a=poly_trim([x%p for x in a]); b=poly_trim([x%p for x in b])
    while b:
        a,b=b,poly_mod_rem(a,b,p)
    if a:
        inv=pow(a[-1],-1,p); a=[x*inv%p for x in a]
    return a

def poly_eval_y(poly,y,p):
    s=0
    for c in reversed(poly): s=(s*y+c)%p
    return s

def roots_pair_mod_fast(c1,c2,p,max_roots=1000):
    roots=[]
    for xv in range(p):
        py1=poly_y_at_x(c1,xv,p); py2=poly_y_at_x(c2,xv,p)
        if not py1 and not py2: continue
        if not py1: g=py2
        elif not py2: g=py1
        else:
            try: g=poly_mod_gcd(py1,py2,p)
            except Exception: continue
        if len(g)>1:
            # find y roots by brute, p should be small
            for yv in range(p):
                if poly_eval_y(g,yv,p)==0:
                    roots.append((xv,yv))
                    if len(roots)>=max_roots: return roots
    return roots

def crt_pair(a,m,b,n):
    t=((b-a)%n)*pow(m%n,-1,n)%n
    return (a+m*t)%(m*n), m*n

def recover_pair(c1,c2,primes,verbose=False,max_roots_per_prime=30,max_cand=1000):
    cand=[(0,0,1)]
    for pr in primes:
        roots=roots_pair_mod_fast(c1,c2,pr,max_roots=max_roots_per_prime+1)
        if verbose: print(' prime',pr,'roots',len(roots),roots[:5], flush=True)
        if not roots or len(roots)>max_roots_per_prime:
            return []
        new=[]
        for xr,yr in roots:
            for ax,ay,M in cand:
                nx,NM=crt_pair(ax,M,xr,pr)
                ny,_=crt_pair(ay,M,yr,pr)
                if NM>X and nx>=X: continue
                if NM>Y and ny>=Y: continue
                new.append((nx,ny,NM))
                if len(new)>max_cand: return []
        cand=new
        if verbose: print(' cand',len(cand),'Mbits',cand[0][2].bit_length() if cand else 0, flush=True)
        if not cand: return []
        if cand[0][2] > max(X,Y): return cand
    return cand

def verify(C,xv,yv):
    if not(0<=xv<X and 0<=yv<Y): return False,None
    p=C + (xv<<x_start) + (yv<<y_start)
    if (p & mask)!=leak: return False,None
    if N%p==0: return True,p
    return False,None

def scan_candidate(low,high,m,t,limit=30,beta=0.499,verbose=True,primes=None,filter_bound=False):
    C=C_for(low,high)
    polys=get_polys(C,m,t,beta=beta,limit=limit,verbose=verbose,filter_bound=filter_bound)
    if verbose:
        print('polys meta', [(a,b,c,d,e) for a,b,c,d,e,_ in polys[:12]], flush=True)
    if primes is None:
        primes=[101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409]
    # skip constant or monomials maybe; try pairs with enough terms and not too sparse
    for ia in range(min(len(polys),limit)):
        for ib in range(ia+1,min(len(polys),limit)):
            c1=polys[ia][-1]; c2=polys[ib][-1]
            # need at least 2 variables involved? but try
            if len(c1)<2 or len(c2)<2: continue
            if verbose: print('TRY',ia,ib,polys[ia][:5],polys[ib][:5], flush=True)
            cand=recover_pair(c1,c2,primes,verbose=verbose,max_roots_per_prime=20,max_cand=200)
            for xv,yv,M in cand:
                ok,p=verify(C,xv,yv)
                if verbose: print(' cand Mbits',M.bit_length(),'xbits',xv.bit_length(),'ybits',yv.bit_length(),'ok',ok, flush=True)
                if ok:
                    print('FOUND low high',low,high)
                    print('p=',hex(p))
                    print('plaintext=',decrypt_from_p(p))
                    return p
    return None

if __name__=='__main__':
    low=int(sys.argv[1]) if len(sys.argv)>1 else 0
    high=int(sys.argv[2]) if len(sys.argv)>2 else 0
    m=int(sys.argv[3]) if len(sys.argv)>3 else 7
    t=int(sys.argv[4]) if len(sys.argv)>4 else 2
    scan_candidate(low,high,m,t,limit=30,verbose=True)
