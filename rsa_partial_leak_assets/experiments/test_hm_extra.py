import random, sys, time, math
from fpylll import IntegerMatrix,LLL
sys.path.append('/mnt/data')
from rsa_partial_lib import mask,leak
XS,XL,YS,YL=265,155,600,230; X=1<<XL; Y=1<<YL; A=1<<XS; B=1<<YS
base=leak & ~(((1<<XL)-1)<<XS) & ~(((1<<YL)-1)<<YS)
lo=random.randrange(16); hi=random.randrange(16); x0=random.randrange(X); y0=random.randrange(Y)
p=base | (lo<<150)|(hi<<920) | (x0<<XS)|(y0<<YS); p|=(1<<1023)|1
q=random.getrandbits(1024)|(1<<1023)|1; N=p*q; C=base|(lo<<150)|(hi<<920)
print('lo hi',lo,hi,'Nbits',N.bit_length(),flush=True)

def mul(P,Q,mod=0):
 r={}
 for (i,j),a in P.items():
  for (k,l),b in Q.items():
   r[(i+k,j+l)]=r.get((i+k,j+l),0)+a*b
 return {m:(c%mod if mod else c) for m,c in r.items() if (c%mod if mod else c)}

def build(m,t,s):
 inv=pow(B,-1,N); f={(0,1):1,(1,0):(A*inv)%N,(0,0):(C*inv)%N} # monic y, shift x
 pw=[{(0,0):1}]
 for k in range(1,m+1): pw.append(mul(pw[-1],f,pow(N,k)))
 polys=[]
 for k in range(m+1):
  max_i=m-k+s
  for i in range(max_i+1):
   mult=pow(N,max(t-k,0)); polys.append({(a+i,b):coef*mult for (a,b),coef in pw[k].items()})
 mons=sorted(set().union(*(P.keys() for P in polys)),key=lambda z:(z[1]+z[0],z[1],z[0]))
 sc=[X**i*Y**j for i,j in mons]; col={u:i for i,u in enumerate(mons)}
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

def ev(P): return sum(c*(x0**i)*(y0**j) for (i,j),c in P.items())
for m,t,s in [(6,2,0),(6,2,2),(6,2,4),(7,2,4),(8,2,4),(8,3,4),(8,3,8),(10,3,6)]:
 M,mons,sc=build(m,t,s); print('param',m,t,s,'dim',M.nrows,M.ncols,flush=True)
 t0=time.time();
 try: LLL.reduction(M,delta=.99) # fast
 except Exception as e: print('LLLerr',e,flush=True); continue
 print('LLL',time.time()-t0,flush=True)
 bits=[]; good=[]
 for r in range(min(M.nrows,30)):
  row=[int(M[r,j]) for j in range(M.ncols)]; bits.append(sum(abs(v) for v in row).bit_length() if any(row) else 0)
  P=unscale(row,mons,sc)
  if P and ev(P)==0: good.append(r)
 print('bits',bits[:8],'good',good[:10],flush=True)
