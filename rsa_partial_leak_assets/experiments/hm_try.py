from math import sqrt, gcd
from fpylll import IntegerMatrix, LLL
import sympy as sp, time, sys

def H(s): return int(''.join(s.split()),16)
MASK=H('''ffffffffffffffff fffffffff0ffffff ffffffffffffffff c00000000000fffe
0000000000000000 000003ffe0000000 0000000000ffffff ffffffffffffffff
ffffffffffffffff fffffff000000000 000003ffe0000000 00000000000001ff
ffffffffffffffff fffffffffc3fffff ffffffffffffffff ffffffffffffffff''')
PAND=H('''ffa360d46885c534 d186538170633faf c2c0548a2e24a2c1 c0000000000039e2
0000000000000000 000000a520000000 00000000003e2de4 c436d2ca740a6246
99e1a1af94045c63 261323c000000000 000003bba0000000 00000000000000e5
0b0bc2461fcbac07 26360c2c0809450a 9a892cbf1d98ceee 48827591ccc593c9''')
N=H('''e505004fb5d34eb7 12d48ff4bbe8d27f c388133c6c0e7340 01061c0ee0a4edc6
37c04fe8dd376185 de8ba04d0ccdbabb 93ab7c371b88d92e 865eec42b028c61d
d7004ebf2ebb5d69 d0a09142be5c9de4 da16e514eea31817 2ecda6cd192073eb
afb1e02d522ec053 34590ea6d75960c4 937bf64f9700db17 7a4aa3da6aae6807
e5e32c0d0e428a0d b68d299f20c235d8 4ef459b0cf118286 59c31663c9ea8204
4b28152c89a9c36c 3ec4303bd36664fd 77fb02c58340bdae 21120326d83fc017
34bc90048dec9fe3 5f08c8fdc523abf8 4a91ec430f495672 37c3153a2035ff62
5613b6dc3e6cb14d 50e18b8a79b25d67 8465b3ad02f5b7d8 18a1e2d635a0baf1''')
ct=H('''8919342826ef3821 5af31e00c9290c4c 50ef9ff9e1afc591 47fab5b096361035
e85f5fc95b73b069 7813b57b831a807d 41bcbecde5b9e663 9e2845b14e395ed0
e5d995e63709ac0c 5ee2337228ee76bc bad857b14904aa2e 8e9997671908a634
d0d1dda1d062ce7f 2e3293ddec8f5cce 26029292d594a062 dcf317d2a8380f43
d72551889efceb87 6c8945a50382272e 76ed6b6fcdff1603 44e9e948e2b6e740
e78bedf25f30e2c7 eeb5f74686c8eadc 29cea04ff08cfd86 dfd3d2a1632bf04a
d5cfa369892a2da4 0f0dc0098ce6b731 d841aab3d0c8b78e b69c4625c47c4ad7
158d49bb5d879581 e02bc525abe47f39 f699864bc5ce1de7 19430dae7aa5480b''')
# dict polynomial in x,y
def add(a,b):
 r=a.copy()
 for m,c in b.items():
  r[m]=r.get(m,0)+c
  if r[m]==0: del r[m]
 return r
def mul(a,b,mod=None):
 r={}
 for (i,j),c in a.items():
  for (k,l),d in b.items():
   m=(i+k,j+l); v=r.get(m,0)+c*d
   if mod: v%=mod
   r[m]=v
 return {m:c for m,c in r.items() if c}
def powmod_poly(f,k,mod):
 r={(0,0):1}
 base={m:c%mod for m,c in f.items()}
 while k:
  if k&1: r=mul(r,base,mod)
  k//=2
  if k: base=mul(base,base,mod)
 return r
def rowpoly(row, mons, X, Y):
 d={}
 for idx,(i,j) in enumerate(mons):
  v=int(row[idx])
  if v:
   den=(X**i)*(Y**j)
   assert v%den==0
   d[(i,j)]=v//den
 return d
def tos(d):
 x,y=sp.symbols('x y')
 return sum(sp.Integer(c)*x**i*y**j for (i,j),c in d.items())
def hm_lattice(f,X,Y,N,m,t):
 polys=[]
 for k in range(m+1):
  fk=powmod_poly(f,k,N**k if k else 1) if k else {(0,0):1}
  mult=N**max(t-k,0)
  for i in range(m-k+1): # y-shifts only
   g={(a,b+i):c*mult for (a,b),c in fk.items()}
   polys.append(g)
 # order as HM: k then i; monomials discovered: x^k y^{i+?}. just include all sorted with row-leading order maybe.
 mons=sorted(set(mm for p in polys for mm in p), key=lambda u:(u[0]+u[1],u[0],u[1]))
 # But for lower triangular maybe sort by (xdeg, ydeg) corresponding? Let's just use all; LLL invariant to col order.
 M=IntegerMatrix(len(polys), len(mons))
 for r,p in enumerate(polys):
  for c,(a,b) in enumerate(mons):
   M[r,c]=p.get((a,b),0)*(X**a)*(Y**b)
 return M,mons

def roots_from_polys(hs,X,Y,limit_pairs=50,verbose=False):
 x,y=sp.symbols('x y')
 sy=[tos(h) for h in hs]
 # remove duplicates / proportional? maybe
 for i in range(min(len(sy),limit_pairs)):
  for j in range(i+1,min(len(sy),limit_pairs)):
   try:
    R=sp.resultant(sy[i], sy[j], x)
    if R==0: continue
    P=sp.Poly(R,y)
    # numeric? factor linear roots only maybe degree manageable.
    # use sp.ground_roots for rational roots
    gr=sp.polys.polytools.ground_roots(P)
    cand=[]
    for rr,mult in gr.items():
     if rr.is_Integer:
      yy=int(rr)
      if 0<=yy<Y: cand.append(yy)
    if verbose and cand: print('cand y',cand,'pair',i,j,'deg',P.degree())
    for yy in cand:
     Pi=sp.Poly(sy[i].subs(y,yy),x)
     grx=sp.polys.polytools.ground_roots(Pi)
     for rr,mult in grx.items():
      if rr.is_Integer:
       xx=int(rr)
       if 0<=xx<X: yield xx,yy
   except Exception as e:
    if verbose: print('res err',i,j,type(e),e)

def try_cand(low,high,m=8,t=None,rows=8,verbose=False):
 if t is None: t=max(1,round((1-sqrt(0.5))*m))
 xmask=((1<<155)-1)<<265; ymask=((1<<230)-1)<<600
 const=(PAND | (low<<150) | (high<<920)) & ~xmask & ~ymask
 a1=1<<265; a2=1<<600
 inv=pow(a1,-1,N)
 A=(a2*inv)%N; B=(const*inv)%N
 # choose signed reps to reduce coefficients?
 if A>N//2: A-=N
 if B>N//2: B-=N
 f={(1,0):1,(0,1):A,(0,0):B}
 X=1<<155; Y=1<<230
 st=time.time(); M,mons=hm_lattice(f,X,Y,N,m,t)
 if verbose: print('lat',M.nrows,M.ncols,'m,t',m,t,'build',time.time()-st)
 st=time.time(); LLL.reduction(M, method='proved', float_type='mpfr', precision=256)
 if verbose: print('lll',time.time()-st)
 hs=[rowpoly([M[r,c] for c in range(M.ncols)], mons, X,Y) for r in range(min(rows,M.nrows))]
 if verbose:
  for idx,h in enumerate(hs[:4]): print('h',idx,'terms',len(h),'deg',max(sum(mm) for mm in h),'normbits',max(abs(c).bit_length() for c in h.values()))
 for xx,yy in roots_from_polys(hs,X,Y,verbose=verbose):
  p=const+(xx<<265)+(yy<<600)
  if (p&MASK)==PAND and p>1 and N%p==0:
   return p,xx,yy
 return None

if __name__=='__main__':
 print('bounds exp bits',155+230)
 for m,t in [(8,2),(9,3),(10,3),(12,4)]:
  print('PARAM',m,t,flush=True)
  for high in range(16):
   for low in range(16):
    print('try',m,t,high,low,flush=True)
    ans=try_cand(low,high,m,t,rows=12,verbose=False)
    if ans:
     p,xx,yy=ans; print('FOUND',m,t,high,low,hex(p));
     q=N//p; phi=(p-1)*(q-1); d=pow(65537,-1,phi); mm=pow(ct,d,N); print(mm.to_bytes((mm.bit_length()+7)//8,'big')); sys.exit()
 print('none')
