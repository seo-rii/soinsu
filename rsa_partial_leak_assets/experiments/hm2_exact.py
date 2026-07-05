import sys, math, time, itertools
sys.path.append('/mnt/data')
from fpylll import IntegerMatrix, LLL, BKZ
from rsa_partial_lib import N,ct,mask,leak,e,C_for,decrypt_from_p,x_start,x_len,y_start,y_len,X,Y

# f = x + A*y + C mod p, p|N. X=2^155,Y=2^230. Shifts: y^i f^k N^(max(t-k,0)), i<=m-k.
A = ((1<<y_start) * pow(1<<x_start, -1, N)) % N

def f_poly(C):
    return {(1,0):1,(0,1):A,(0,0):(C*pow(1<<x_start,-1,N))%N}

def mul(P,Q):
    R={}
    for (i,j),a in P.items():
      for (k,l),b in Q.items():
        R[(i+k,j+l)] = R.get((i+k,j+l),0)+a*b
    return {m:c for m,c in R.items() if c}

def shift(P, sx, sy, coef=1):
    return {(i+sx,j+sy):c*coef for (i,j),c in P.items() if c*coef}

def build(C,m,t):
    F=f_poly(C); P=[{(0,0):1}]
    for k in range(1,m+1): P.append(mul(P[-1],F))
    polys=[]
    for k in range(m+1):
      for i in range(m-k+1):
        polys.append(shift(P[k],0,i,pow(N,max(t-k,0))))
    mons=sorted(set().union(*[set(p) for p in polys]), key=lambda z:(z[0]+z[1],z[0],z[1]))
    col={m:i for i,m in enumerate(mons)}
    sc=[(X**i)*(Y**j) for i,j in mons]
    B=[]
    for P0 in polys:
      row=[0]*len(mons)
      for mon,c in P0.items(): row[col[mon]]=int(c)*sc[col[mon]]
      B.append(row)
    return B,mons,sc

def lll(B,delta=0.99):
    M=IntegerMatrix(len(B),len(B[0]))
    for i,r in enumerate(B):
      for j,v in enumerate(r): M[i,j]=int(v)
    st=time.time(); LLL.reduction(M,delta=delta,method='proved',float_type='mpfr',precision=256)
    return [[int(M[i,j]) for j in range(M.ncols)] for i in range(M.nrows)], time.time()-st

def unscale(row,mons,sc):
    P={}
    for v,mon,s in zip(row,mons,sc):
      if v:
        if v%s: return None
        P[mon]=v//s
    return P

def evalP(P,x,y): return sum(c*pow(x,i)*pow(y,j) for (i,j),c in P.items())
def bits_l1(r): return (sum(abs(x) for x in r)).bit_length()
def bits_max(r): return max([abs(x).bit_length() for x in r if x] or [0])

def trim(a):
    while a and a[-1]%MOD==0: a.pop()
    return [x%MOD for x in a]
def prem(a,b):
    a=trim(a[:]); b=trim(b[:]); inv=pow(b[-1],-1,MOD)
    while len(a)>=len(b) and a:
      c=a[-1]*inv%MOD
      if c:
        off=len(a)-len(b)
        for j,bj in enumerate(b): a[off+j]=(a[off+j]-c*bj)%MOD
      trim(a)
    return a
def pgcd(a,b):
    a=trim(a); b=trim(b)
    while b: a,b=b,prem(a,b)
    if a:
      inv=pow(a[-1],-1,MOD); a=[x*inv%MOD for x in a]
    return a
def peval(a,y):
    s=0
    for c in reversed(a): s=(s*y+c)%MOD
    return s

def ypoly(P,x,p):
    mi=max([i for i,j in P] or [0]); mj=max([j for i,j in P] or [0])
    xp=[1]*(mi+1)
    for i in range(1,mi+1): xp[i]=xp[i-1]*x%p
    a=[0]*(mj+1)
    for (i,j),c in P.items(): a[j]=(a[j]+(c%p)*xp[i])%p
    return a

def roots_mod(polys,p):
    global MOD; MOD=p
    roots=[]
    P0=polys[0]
    for x in range(p):
      g=None
      for P in polys:
        a=ypoly(P,x,p)
        if not trim(a): continue
        g=a if g is None else pgcd(g,a)
        if not g: break
      if g and len(g)>1 and len(g)<=5:
        for y in range(p):
          if peval(g,y)==0 and all(sum((c%p)*pow(x,i,p)*pow(y,j,p) for (i,j),c in P.items())%p==0 for P in polys):
            roots.append((x,y))
    return roots

def crt(a,m,b,n):
    t=((b-a)%n)*pow(m%n,-1,n)%n
    return (a+m*t)%(m*n)

def lift_roots(polys, primes=(31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107)):
    cand=[(0,0,1)]
    for p in primes:
      rs=roots_mod(polys,p)
      print('  mod',p,'roots',len(rs),rs[:6],flush=True)
      if not rs or len(rs)>20: return []
      nc=[]
      for x,y in rs:
        for a,b,M in cand:
          nc.append((crt(a,M,x,p),crt(b,M,y,p),M*p))
      cand=nc[:2000]
      if not cand: return []
      if cand[0][2] > max(X,Y): return cand
    return cand

def solve(low,high,m,t,limit=12):
    C=C_for(low,high); B,mons,sc=build(C,m,t); R,dt=lll(B)
    print(f'cand low={low} high={high} m={m} t={t} dim={len(B)} time={dt:.1f}s',flush=True)
    print(' bits', [bits_l1(r) for r in R[:10]], 'target', int(0.5*t*N.bit_length()), flush=True)
    # take short nontrivial polynomials; avoid constants/monomials
    polys=[]
    for r in R:
      P=unscale(r,mons,sc)
      if not P: continue
      if len(P)>=2 and max(i+j for i,j in P)>0:
        polys.append(P)
      if len(polys)>=limit: break
    # try subsets increasing
    for k in range(2,min(6,len(polys))+1):
      for idxs in itertools.combinations(range(min(len(polys),8)), k):
        print(' try',idxs,flush=True)
        for x,y,M in lift_roots([polys[i] for i in idxs]):
          if x<X and y<Y:
            Pval=C + (x<<x_start) + (y<<y_start)
            if (Pval&mask)==leak and N%Pval==0:
              print('FOUND',low,high,hex(Pval)); print(decrypt_from_p(Pval)); return Pval
    return None

if __name__=='__main__':
    low=int(sys.argv[1]); high=int(sys.argv[2]); m=int(sys.argv[3]); t=int(sys.argv[4])
    solve(low,high,m,t)
