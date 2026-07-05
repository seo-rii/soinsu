import random, sys, time, math
from fpylll import IntegerMatrix,LLL
# constants mask/leak from rsa_partial_lib
sys.path.append('/mnt/data')
from rsa_partial_lib import mask,leak
XS,XL,YS,YL=265,155,600,230; X=1<<XL; Y=1<<YL; A=1<<XS; B=1<<YS
# generate p matching leak
base=leak & ~(((1<<XL)-1)<<XS) & ~(((1<<YL)-1)<<YS)
lo=random.randrange(16); hi=random.randrange(16)
x0=random.randrange(X); y0=random.randrange(Y)
p=base | (lo<<150)|(hi<<920) | (x0<<XS)|(y0<<YS)
p |= (1<<1023)|1
q=random.getrandbits(1024)|(1<<1023)|1
N=p*q
C=base | (lo<<150)|(hi<<920)
print('lo hi',lo,hi,'Nbits',N.bit_length())

def mul(P,Q,mod=0):
 r={}
 for (i,j),a in P.items():
  for (k,l),b in Q.items():
   r[(i+k,j+l)]=r.get((i+k,j+l),0)+a*b
 return {m:(c%mod if mod else c) for m,c in r.items() if (c%mod if mod else c)}

def build(m,t,lead='y'):
 if lead=='x': inv=pow(A,-1,N); f={(1,0):1,(0,1):(B*inv)%N,(0,0):(C*inv)%N}; sh=lambda k,i:(0,i)
 else: inv=pow(B,-1,N); f={(0,1):1,(1,0):(A*inv)%N,(0,0):(C*inv)%N}; sh=lambda k,i:(i,0)
 pw=[{(0,0):1}]
 for k in range(1,m+1): pw.append(mul(pw[-1],f,pow(N,k)))
 polys=[]
 for k in range(m+1):
  for i in range(m-k+1):
   sx,sy=sh(k,i); mult=pow(N,max(t-k,0))
   polys.append({(a+sx,b+sy):coef*mult for (a,b),coef in pw[k].items()})
 mons=sorted(set().union(*(P.keys() for P in polys)),key=lambda z:(z[0]+z[1],z[0],z[1])); sc=[X**i*Y**j for i,j in mons]; col={u:i for i,u in enumerate(mons)}
 M=IntegerMatrix(len(polys),len(mons))
 for r,P in enumerate(polys):
  for mon,coef in P.items(): M[r,col[mon]]=int(coef*sc[col[mon]])
 return M,mons,sc

def unscale(row,mons,sc):
 d={}
 for v,mon,s in zip(row,mons,sc):
  v=int(v)
  if v:
   if v%s: return None
   d[mon]=v//s
 return d

def evalp(P,x,y):
 return sum(c*(x**i)*(y**j) for (i,j),c in P.items())
for m,t in [(6,2),(8,2),(8,3),(10,3)]:
 M,mons,sc=build(m,t,'y')
 t0=time.time(); LLL.reduction(M,delta=.99,method='proved',float_type='mpfr',precision=256); print('m,t',m,t,'dim',M.nrows,M.ncols,'LLL',time.time()-t0)
 good=[]; bits=[]
 for r in range(min(M.nrows,20)):
  row=[int(M[r,j]) for j in range(M.ncols)]
  bits.append(sum(abs(v) for v in row).bit_length() if any(row) else 0)
  P=unscale(row,mons,sc)
  if P and evalp(P,x0,y0)==0:
   good.append(r)
 print('bits',bits[:10],'good rows',good[:10])
