from fpylll import IntegerMatrix, LLL
import sympy as sp, random

def poly_mul(a,b):
 r={}
 for i,c in a.items():
  for j,d in b.items(): r[i+j]=r.get(i+j,0)+c*d
 return {i:c for i,c in r.items() if c}
def poly_pow(f,k):
 r={0:1}
 for _ in range(k): r=poly_mul(r,f)
 return r
def build(f,N,X,m,t):
 rows=[]; fp=[poly_pow(f,k) for k in range(m+1)]
 for k in range(m+1):
  base={i:c*pow(N,m-k) for i,c in fp[k].items()}
  maxshift = t if k==m else t # maybe m-k?
  for j in range(maxshift+1): rows.append({i+j:c for i,c in base.items()})
 mons=sorted(set().union(*[set(r) for r in rows]))
 M=IntegerMatrix(len(rows),len(mons))
 for r,row in enumerate(rows):
  for c,mon in enumerate(mons):
   val=row.get(mon,0)
   if val: M[r,c]=val*pow(X,mon)
 return M,mons

def rec(v,mons,X):
 p={}
 for val,mon in zip(v,mons):
  if val:
   den=pow(X,mon); assert val%den==0
   p[mon]=int(val//den)
 return p

def ev(p,x): return sum(c*pow(x,i) for i,c in p.items())
for bits in [128,192,256]:
 p=int(sp.randprime(1<<(bits//2-1),1<<(bits//2))); q=int(sp.randprime(1<<(bits//2-1),1<<(bits//2))); N=p*q
 known=bits//4+5 # known lsb bits, unknown p high bits length bits/4-5 perhaps
 b=known; x0=p>>b; p0=p & ((1<<b)-1); X=1<<(bits//2-b)
 f={0:p0,1:1<<b}
 print('bits',bits,'x bits',x0.bit_length(),'known',b)
 for m,t in [(2,1),(3,1),(4,1),(3,2),(4,2),(5,2),(6,2)]:
  M,mons=build(f,N,X,m,t)
  LLL.reduction(M,method='proved',float_type='mpfr',precision=128)
  vals=[]
  for r in range(min(M.nrows,10)):
   pol=rec([M[r,c] for c in range(M.ncols)],mons,X)
   val=ev(pol,x0)
   vals.append((0 if val==0 else abs(val).bit_length(), val==0, len(pol)))
  print('m,t,dim',m,t,M.nrows,M.ncols,vals[:5])
