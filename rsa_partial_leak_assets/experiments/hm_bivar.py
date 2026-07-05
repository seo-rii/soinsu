from fpylll import IntegerMatrix, LLL
from math import prod
from itertools import combinations
import sympy as sp

# Polynomial dict: monomial (i,j) -> int coeff

def poly_clean(p):
    return {m:c for m,c in p.items() if c}

def poly_add(a,b):
    r=a.copy()
    for m,c in b.items():
        r[m]=r.get(m,0)+c
        if r[m]==0: del r[m]
    return r

def poly_mul(a,b):
    r={}
    for (i,j),c in a.items():
        for (k,l),d in b.items():
            m=(i+k,j+l)
            r[m]=r.get(m,0)+c*d
    return poly_clean(r)

def poly_pow(a,n):
    r={(0,0):1}
    b=a
    while n:
        if n&1: r=poly_mul(r,b)
        n//=2
        if n: b=poly_mul(b,b)
    return r

def poly_shift(a, sx, sy, coeff=1):
    return {(i+sx,j+sy):c*coeff for (i,j),c in a.items()}

def iiter(m,n):
    arr=[0]*n; s=0; stop=False
    while True:
        yield tuple(arr)
        # next like sage code
        done=True
        for idx in range(n-1,-1,-1):
            if s==m or arr[idx]==m:
                s-=arr[idx]; arr[idx]=0; continue
            arr[idx]+=1; s+=1; done=False; break
        if done: break

def hm_polys_bivar(N, c0, cy, X, Y, m=8, t=3, lll_rows=None, verbose=True):
    # f = x + cy*y + c0 over ZZ residues
    f={(1,0):1,(0,1):cy%N,(0,0):c0%N}
    g=[]; monoms=[]; Xmuls=[]
    for ii in iiter(m,2):
        k=ii[0]; jy=ii[1]
        pk=poly_pow(f,k)
        mult=pow(N, max(t-k,0))
        gi=poly_shift(pk, 0, jy, mult)
        g.append(gi); monoms.append((k,jy)); Xmuls.append((X**k)*(Y**jy))
    n=len(g)
    M=[[0]*n for _ in range(n)]
    for i,gi in enumerate(g):
        for j in range(i+1):
            mon=monoms[j]
            coeff=gi.get(mon,0)
            # for constant j=0 Xmuls[0]=1
            M[i][j]=coeff*Xmuls[j]
    if verbose:
        print('LLL dim',n,'m,t',m,t,'maxbits',max(abs(v).bit_length() for row in M for v in row if v), flush=True)
    B=IntegerMatrix.from_matrix(M)
    LLL.reduction(B, delta=0.99, eta=0.51, method='proved', float_type='mpfr', precision=128)
    if verbose: print('LLL done', flush=True)
    rows=lll_rows or n
    hs=[]
    for i in range(min(rows,n)):
        h={}
        for j,mon in enumerate(monoms):
            val=int(B[i,j])
            if val:
                xm=Xmuls[j]
                if val % xm != 0:
                    # Should not happen with exact LLL? may happen due to combinations? Actually lattice columns scaled, entries remain multiples.
                    pass
                coeff=val//xm
                h[mon]=h.get(mon,0)+coeff
        hs.append(poly_clean(h))
    return hs, monoms

def eval_poly(p,x,y):
    s=0
    xpows={0:1}; ypows={0:1}
    for i,j in p:
        if i not in xpows: xpows[i]=x**i
        if j not in ypows: ypows[j]=y**j
    for (i,j),c in p.items(): s += c*xpows[i]*ypows[j]
    return s

def poly_to_sympy_mod(p, mod=None):
    x,y=sp.symbols('x y')
    expr=0
    if mod is None:
        for (i,j),c in p.items(): expr += c*x**i*y**j
        return sp.Poly(expr,x,y,domain=sp.ZZ)
    else:
        for (i,j),c in p.items(): expr += (c%mod)*x**i*y**j
        return sp.Poly(expr,x,y, modulus=mod)

def roots_mod_pair(p1,p2,prime):
    x,y=sp.symbols('x y')
    P1=poly_to_sympy_mod(p1, prime)
    P2=poly_to_sympy_mod(p2, prime)
    try:
        R=sp.resultant(P1.as_expr(), P2.as_expr(), y)
        Rpoly=sp.Poly(R, x, modulus=prime)
    except Exception as e:
        return []
    if Rpoly.is_zero:
        return []
    # factor/roots mod prime. brute if prime small? choose primes ~2^15, degree <=m, can use ground_roots
    try:
        xrts=[int(r)%prime for r,mul in Rpoly.ground_roots().items()]
    except Exception:
        # brute fallback for small prime only
        xrts=[]
        for xv in range(prime):
            if Rpoly.eval(xv)%prime==0: xrts.append(xv)
    out=[]
    for xv in xrts:
        # gcd of univariate polynomials in y
        Q1=sp.Poly(P1.as_expr().subs(x,xv), y, modulus=prime)
        Q2=sp.Poly(P2.as_expr().subs(x,xv), y, modulus=prime)
        G=sp.gcd(Q1,Q2)
        if G.is_zero: continue
        try:
            yrts=[int(r)%prime for r,mul in G.ground_roots().items()]
        except Exception:
            yrts=[]
            for yv in range(prime):
                if Q1.eval(yv)%prime==0 and Q2.eval(yv)%prime==0: yrts.append(yv)
        for yv in yrts:
            out.append((xv,yv))
    return list(set(out))

def crt_pair(a,m,b,n):
    # m,n coprime
    inv=pow(m,-1,n)
    t=((b-a)%n)*inv % n
    return a + m*t, m*n

def recover_roots_crt(hs, X, Y, verify_func, primes=None, max_pairs=20, verbose=True):
    # Try pairs of first hs; use CRT with unique-ish roots over primes.
    if primes is None:
        primes=[10007,10009,10037,10039,10061,10067,10069,10079,10091,10093,10099,10103,10111,10133,10139,10141,10151,10159,10163,10169,10177,10181,10193,10211,10223,10243,10247,10253,10259,10267]
    pairs=list(combinations(range(min(len(hs),10)),2))[:max_pairs]
    for ia,ib in pairs:
        p1,p2=hs[ia],hs[ib]
        candidates=[(0,0,1)]  # xr, yr, mod
        if verbose: print('pair',ia,ib,'terms',len(p1),len(p2), flush=True)
        for pr in primes:
            roots=roots_mod_pair(p1,p2,pr)
            if verbose: print(' prime',pr,'roots',len(roots),'cand',len(candidates), flush=True)
            if not roots:
                candidates=[]; break
            new=[]
            for xr,yr,mod in candidates:
                for xmod,ymod in roots:
                    nx,nmod=crt_pair(xr,mod,xmod,pr)
                    ny,_=crt_pair(yr,mod,ymod,pr)
                    # canonical positive root; if nmod exceeds 2*bound, can check exact direct
                    if nmod > 2*max(X,Y):
                        xval=nx if nx<X else nx-nmod
                        yval=ny if ny<Y else ny-nmod
                        if 0 <= xval < X and 0 <= yval < Y:
                            res=verify_func(xval,yval)
                            if res: return res
                    # prune if residues impossible? keep canonical positive, root assumed nonnegative
                    # if mod > bound, residue must equal value, so if nx>=X impossible unless nx-nmod in range but root nonneg no
                    if nmod > X and nx >= X: continue
                    if nmod > Y and ny >= Y: continue
                    new.append((nx,ny,nmod))
            # Dedup and cap
            seen=set(); candidates2=[]
            for c in new:
                key=(c[0],c[1],c[2])
                if key not in seen:
                    seen.add(key); candidates2.append(c)
            candidates=candidates2[:10000]
            if not candidates: break
            if candidates and candidates[0][2] > max(X,Y):
                for xr,yr,mod in candidates:
                    if 0 <= xr < X and 0 <= yr < Y:
                        res=verify_func(xr,yr)
                        if res: return res
                break
    return None

def mon_div(a,b):
    if a[0] >= b[0] and a[1] >= b[1]:
        return (a[0]-b[0], a[1]-b[1])
    return None

def poly_monomials_power_linear(m):
    # monomials of (c0 + x + cy y)^m = all (i,j) with i+j<=m
    return {(i,j) for i in range(m+1) for j in range(m+1-i)}

def generic_shift_polys_bivar(N,c0,cy,X,Y,m=2,d=4,lll_rows=None,verbose=True):
    # f leading monomial x under lex x>y -> l=(1,0); f monic in x
    f={(1,0):1,(0,1):cy%N,(0,0):c0%N}
    l=(1,0)
    M=[]
    fm_mons=poly_monomials_power_linear(m)
    for k in range(m+1):
        Mk=set(); T=poly_monomials_power_linear(m-k)
        lk=(k,0)
        for mon in fm_mons:
            div=mon_div(mon,lk)
            if div is not None and div in T:
                for ex in range(d):
                    for ey in range(d):
                        Mk.add((mon[0]+ex, mon[1]+ey))
        M.append(Mk)
    M.append(set())
    shifts=[]
    for k in range(m+1):
        fk=poly_pow(f,k)
        for mon in sorted(M[k]-M[k+1]):
            div=mon_div(mon,(k,0))
            if div is None: continue
            shifts.append(poly_shift(fk, div[0], div[1], pow(N,m-k)))
    # collect monomials sorted by total then x/y maybe
    monoms=sorted(set().union(*[set(s.keys()) for s in shifts]), key=lambda ij:(ij[0]+ij[1],ij[0],ij[1]))
    factors=[(X**i)*(Y**j) for i,j in monoms]
    nr=len(shifts); nc=len(monoms)
    Mmat=[[0]*nc for _ in range(nr)]
    for i,sh in enumerate(shifts):
        for j,mon in enumerate(monoms):
            Mmat[i][j]=sh.get(mon,0)*factors[j]
    if verbose:
        print('Generic shift LLL rows/cols',nr,nc,'m,d',m,d,'maxbits',max(abs(v).bit_length() for row in Mmat for v in row if v),flush=True)
    B=IntegerMatrix.from_matrix(Mmat)
    LLL.reduction(B, delta=0.99, eta=0.51, method='proved', float_type='mpfr', precision=128)
    if verbose: print('LLL done',flush=True)
    rows=lll_rows or min(nr,nc)
    hs=[]
    for i in range(min(rows,B.nrows)):
        h={}
        for j,mon in enumerate(monoms):
            val=int(B[i,j])
            if val:
                factor=factors[j]
                coeff=val//factor if val%factor==0 else int(round(val/factor))
                h[mon]=h.get(mon,0)+coeff
        hs.append(poly_clean(h))
    return hs,monoms

def hm_polys_bivar_balanced(N, c0, cy, X, Y, m=8, t=4, lll_rows=None, verbose=True):
    def bal(a):
        a%=N
        return a-N if a>N//2 else a
    f={(1,0):1,(0,1):bal(cy),(0,0):bal(c0)}
    g=[]; monoms=[]; Xmuls=[]
    for ii in iiter(m,2):
        k=ii[0]; jy=ii[1]
        pk=poly_pow(f,k)
        mult=pow(N, max(t-k,0))
        gi=poly_shift(pk, 0, jy, mult)
        g.append(gi); monoms.append((k,jy)); Xmuls.append((X**k)*(Y**jy))
    n=len(g)
    M=[[0]*n for _ in range(n)]
    for i,gi in enumerate(g):
        for j in range(i+1):
            M[i][j]=gi.get(monoms[j],0)*Xmuls[j]
    if verbose:
        print('BAL LLL dim',n,'m,t',m,t,'maxbits',max(abs(v).bit_length() for row in M for v in row if v), flush=True)
    B=IntegerMatrix.from_matrix(M)
    LLL.reduction(B, delta=0.99, eta=0.51, method='proved', float_type='mpfr', precision=128)
    if verbose: print('LLL done', flush=True)
    hs=[]
    for i in range(min(lll_rows or n,n)):
        h={}
        for j,mon in enumerate(monoms):
            val=int(B[i,j])
            if val: h[mon]=h.get(mon,0)+val//Xmuls[j]
        hs.append(poly_clean(h))
    return hs,monoms

def hm_polys_bivar_bkz(N, c0, cy, X, Y, m=8, t=4, block=20, lll_rows=None, verbose=True):
    f={(1,0):1,(0,1):cy%N,(0,0):c0%N}
    g=[]; monoms=[]; Xmuls=[]
    for ii in iiter(m,2):
        k,jy=ii; g.append(poly_shift(poly_pow(f,k),0,jy,pow(N,max(t-k,0)))); monoms.append(ii); Xmuls.append(X**k*Y**jy)
    n=len(g); M=[[0]*n for _ in range(n)]
    for i,gi in enumerate(g):
        for j in range(i+1): M[i][j]=gi.get(monoms[j],0)*Xmuls[j]
    if verbose: print('BKZ LLL dim',n,'block',block,flush=True)
    B=IntegerMatrix.from_matrix(M)
    LLL.reduction(B,delta=0.99,eta=0.51,method='proved',float_type='mpfr',precision=128)
    if block:
        from fpylll import BKZ
        par=BKZ.Param(block_size=block, max_loops=2)
        BKZ.reduction(B,par)
    if verbose: print('red done',flush=True)
    hs=[]
    for i in range(min(lll_rows or n,n)):
        h={}
        for j,mon in enumerate(monoms):
            val=int(B[i,j])
            if val: h[mon]=h.get(mon,0)+val//Xmuls[j]
        hs.append(poly_clean(h))
    return hs,monoms

def linear_bivar_multilead(N, p0, ax, ay, X, Y, m=6, t=3, lll_rows=None, verbose=True):
    # Build kionactf-style Herrmann-May using both variables as leading.
    # f = p0 + ax*x + ay*y mod N; normalize by ax then by ay.
    invx=pow(ax,-1,N); invy=pow(ay,-1,N)
    fx={(1,0):1,(0,1):(ay*invx)%N,(0,0):(p0*invx)%N}
    fy={(0,1):1,(1,0):(ax*invy)%N,(0,0):(p0*invy)%N}
    shifts=[]
    for f,other in [(fx,'y'),(fy,'x')]:
        for k in range(m+1):
            fk=poly_pow(f,k); mult=pow(N,max(t-k,0))
            for j in range(m-k+1):
                if other=='y': shifts.append(poly_shift(fk,0,j,mult))
                else: shifts.append(poly_shift(fk,j,0,mult))
    monoms=sorted(set().union(*[set(s.keys()) for s in shifts]), key=lambda ij:(ij[0]+ij[1],ij[0],ij[1]))
    factors=[X**i*Y**j for i,j in monoms]
    nr=len(shifts); nc=len(monoms)
    M=[[0]*nc for _ in range(nr)]
    for i,sh in enumerate(shifts):
        for j,mon in enumerate(monoms): M[i][j]=sh.get(mon,0)*factors[j]
    if verbose: print('linear multilead LLL rows/cols',nr,nc,'m,t',m,t,'maxbits',max(abs(v).bit_length() for row in M for v in row if v),flush=True)
    B=IntegerMatrix.from_matrix(M)
    LLL.reduction(B,delta=0.99,eta=0.51,method='proved',float_type='mpfr',precision=128)
    if verbose: print('LLL done',flush=True)
    rows=lll_rows or min(nr,nc)
    hs=[]
    for i in range(min(rows,B.nrows)):
        h={}
        for j,mon in enumerate(monoms):
            val=int(B[i,j])
            if val: h[mon]=h.get(mon,0)+val//factors[j]
        hs.append(poly_clean(h))
    return hs,monoms
