from fpylll import IntegerMatrix, LLL
from itertools import combinations
import sympy as sp

def clean(p): return {m:c for m,c in p.items() if c}
def mul(a,b):
    n=len(next(iter(a))) if a else len(next(iter(b)))
    r={}
    for ma,ca in a.items():
      for mb,cb in b.items():
        m=tuple(ma[i]+mb[i] for i in range(n)); r[m]=r.get(m,0)+ca*cb
    return clean(r)
def powpoly(a,e):
    n=len(next(iter(a)))
    r={(0,)*n:1}; b=a
    while e:
      if e&1: r=mul(r,b)
      e//=2
      if e: b=mul(b,b)
    return r
def shift(a, exps, coeff=1):
    return {tuple(m[i]+exps[i] for i in range(len(exps))): c*coeff for m,c in a.items()}
def iiter(m,n):
    arr=[0]*n; s=0
    while True:
      yield tuple(arr)
      done=True
      for idx in range(n-1,-1,-1):
        if s==m or arr[idx]==m:
          s-=arr[idx]; arr[idx]=0; continue
        arr[idx]+=1; s+=1; done=False; break
      if done: break

def hm_polys(N,c0,cs,bounds,m,t,rows=None,prec=128,verbose=True):
    n=1+len(cs)
    f={(1,)+((0,)*(n-1)):1,(0,)*n:c0%N}
    for idx,c in enumerate(cs, start=1):
        mon=[0]*n; mon[idx]=1; f[tuple(mon)]=c%N
    gs=[]; mons=[]; factors=[]
    for ii in iiter(m,n):
        k=ii[0]
        gi=shift(powpoly(f,k), (0,)+ii[1:], pow(N,max(t-k,0)))
        gs.append(gi); mons.append(ii)
        fac=1
        for b,e in zip(bounds,ii): fac*=b**e
        factors.append(fac)
    dim=len(gs); M=[[0]*dim for _ in range(dim)]
    for i,gi in enumerate(gs):
        for j in range(i+1): M[i][j]=gi.get(mons[j],0)*factors[j]
    if verbose: print('HM n',n,'dim',dim,'m,t',m,t,'maxbits',max(abs(v).bit_length() for row in M for v in row if v),flush=True)
    B=IntegerMatrix.from_matrix(M)
    LLL.reduction(B,delta=0.99,eta=0.51,method='proved',float_type='mpfr',precision=prec)
    if verbose: print('LLL done',flush=True)
    hs=[]
    for i in range(min(rows or dim,dim)):
        h={}
        for j,mon in enumerate(mons):
            val=int(B[i,j])
            if val:
                fac=factors[j]
                if val%fac: pass
                h[mon]=h.get(mon,0)+val//fac
        hs.append(clean(h))
    return hs, mons

def evalp(p, vals):
    s=0; n=len(vals); powcache=[{0:1} for _ in vals]
    for m,c in p.items():
      term=c
      for i,e in enumerate(m):
        if e not in powcache[i]: powcache[i][e]=vals[i]**e
        term*=powcache[i][e]
      s+=term
    return s

def to_sympy_mod(p,mod,vars):
    expr=0
    for m,c in p.items():
        term=c%mod
        for v,e in zip(vars,m):
            if e: term*=v**e
        expr+=term
    return sp.Poly(expr,*vars,modulus=mod)

def solve_mod_groebner(polys,prime,nvars):
    vs=sp.symbols('x0:'+str(nvars))
    exprs=[]
    for p in polys:
        P=to_sympy_mod(p,prime,vs)
        if not P.is_zero: exprs.append(P.as_expr())
    if len(exprs)<nvars: return []
    try:
        G=sp.groebner(exprs,*vs,modulus=prime,order='lex')
    except Exception as e:
        return []
    # zero-dimensional? Use triangular form if possible; else brute last variable? We'll use solve_poly_system? No modulus.
    # Extract univariate polys for vars and recursively brute roots using triangular basis.
    basis=[sp.Poly(g,*vs,modulus=prime) for g in G.polys]
    # recursive assignment from first variable? Lex usually has univar in last? inspect all univariate.
    sols=[{}]
    # Find variables one by one by any polynomial becoming univariate under assignments.
    for _ in range(nvars):
        progress=False; new_sols=[]
        for sol in sols:
            # if complete
            if len(sol)==nvars:
                new_sols.append(sol); continue
            found=False
            for var in reversed(vs):
                if var in sol: continue
                # build gcd/product of polys involving only this var after substituting knowns
                candidates=[]
                for P in basis:
                    expr=P.as_expr().subs(sol)
                    rem_vars=[v for v in vs if v not in sol and v!=var]
                    if any(expr.has(v) for v in rem_vars): continue
                    Q=sp.Poly(expr,var,modulus=prime)
                    if Q.is_zero or Q.degree()<1: continue
                    roots=list(Q.ground_roots().keys())
                    if roots:
                        candidates=[int(r)%prime for r in roots]; break
                if candidates:
                    for r in candidates:
                        ns=sol.copy(); ns[var]=r; new_sols.append(ns)
                    progress=True; found=True; break
            if not found and len(sol)<nvars:
                # fallback do not keep
                pass
        sols=new_sols
        if not progress: break
    out=[]
    for sol in sols:
        if len(sol)==nvars:
            vals=tuple(int(sol[v])%prime for v in vs)
            ok=True
            for e in exprs:
                if int(e.subs({vs[i]:vals[i] for i in range(nvars)}))%prime:
                    ok=False; break
            if ok: out.append(vals)
    return list(set(out))

def crt(a,m,b,n):
    t=((b-a)%n)*pow(m,-1,n)%n
    return a+m*t, m*n

def recover_crt_groebner(hs,bounds,verify,primes=None,npolys=None,verbose=True):
    n=len(bounds)
    if primes is None: primes=[10007,10009,10037,10039,10061,10067,10069,10079,10091,10093,10099,10103,10111]
    # Try subsets of first rows; initially first npolys rows.
    subsets=[]
    r=npolys or min(len(hs), n+2)
    subsets.append(tuple(range(r)))
    # also combinations if needed
    for sub in subsets:
        cands=[((0,)*n,1)]
        if verbose: print('subset',sub,flush=True)
        for pr in primes:
            roots=solve_mod_groebner([hs[i] for i in sub], pr, n)
            if verbose: print(' prime',pr,'roots',len(roots),'cand',len(cands),flush=True)
            if not roots: cands=[]; break
            new=[]
            for residues,mod in cands:
                for root in roots:
                    vals=[]; nmod=None; ok=True
                    for idx in range(n):
                        nv,nm=crt(residues[idx],mod,root[idx],pr)
                        if nm>bounds[idx] and nv>=bounds[idx]: ok=False; break
                        vals.append(nv); nmod=nm
                    if ok:
                        if nmod>max(bounds):
                            res=verify(tuple(vals))
                            if res: return res
                        new.append((tuple(vals),nmod))
            # dedup cap
            seen=set(); c2=[]
            for vals,mod in new:
                key=(vals,mod)
                if key not in seen: seen.add(key); c2.append((vals,mod))
            cands=c2[:10000]
            if not cands: break
    return None
