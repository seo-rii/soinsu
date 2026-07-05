import sys,math,time,itertools
sys.path.append('/mnt/data')
from fpylll import IntegerMatrix, LLL
from rsa_partial_lib import N,leak,mask,C_for
chunks=[(265,84),(362,58),(600,69),(682,87),(784,46)]
# use leading first var; clear variable bit ranges already done? C_for clears grouped not individual, but okay create base clearing chunks

def clear(v,s,l): return v & ~(((1<<l)-1)<<s)
base=leak
for s,l in chunks: base=clear(base,s,l)
def C(low,high): return base | (low<<150)|(high<<920)

def mul(P,Q):
 R={}
 for a,ca in P.items():
  for b,cb in Q.items():
   m=tuple(x+y for x,y in zip(a,b)); R[m]=R.get(m,0)+ca*cb
 return {m:c for m,c in R.items() if c}

def sh(P,exp,coef):
 return {tuple(a+b for a,b in zip(m,exp)):c*coef for m,c in P.items() if c*coef}

def build(C0,m,t):
 n=len(chunks); a=[pow(1<<chunks[0][0],-1,N)*(1<<s)%N for s,l in chunks]
 c=C0*pow(1<<chunks[0][0],-1,N)%N
 F={(0,)*n:c}; F[(1,)+(0,)*(n-1)]=1
 for j in range(1,n):
  e=[0]*n; e[j]=1; F[tuple(e)]=a[j]
 P=[{(0,)*n:1}]
 for k in range(1,m+1): P.append(mul(P[-1],F))
 polys=[]
 for k in range(m+1):
  # shifts over variables 2..n total <= m-k
  for exps_tail in itertools.product(range(m+1), repeat=n-1):
   if sum(exps_tail)<=m-k:
    polys.append(sh(P[k], (0,)+exps_tail, pow(N,max(t-k,0))))
 mons=sorted(set().union(*[set(p) for p in polys]), key=lambda z:(sum(z),z))
 col={mon:i for i,mon in enumerate(mons)}; bounds=[1<<l for s,l in chunks]
 sc=[]
 for mon in mons:
  v=1
  for B,e in zip(bounds,mon): v*=B**e
  sc.append(v)
 rows=[]
 for P0 in polys:
  r=[0]*len(mons)
  for mon,c0 in P0.items(): r[col[mon]]=int(c0)*sc[col[mon]]
  rows.append(r)
 return rows
rows=build(C(0,0),4,1); print('dim',len(rows),len(rows[0]))
A=IntegerMatrix(len(rows),len(rows[0]))
for i,r in enumerate(rows):
 for j,v in enumerate(r): A[i,j]=int(v)
st=time.time(); LLL.reduction(A,delta=0.99,method='proved',float_type='mpfr',precision=192); print('time',time.time()-st)
def l1(i): return sum(abs(int(A[i,j])) for j in range(A.ncols)).bit_length()
print([l1(i) for i in range(12)], 'target',1024)
