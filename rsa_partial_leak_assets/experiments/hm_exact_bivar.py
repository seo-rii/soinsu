import sys,math,time,itertools
sys.path.append('/mnt/data')
from fpylll import IntegerMatrix, LLL, BKZ
from rsa_partial_lib import N,ct,mask,leak,e,decrypt_from_p
# constants for p=C+2^265*x+2^600*y, after guessing p[150..153],p[920..923]
XS,YS=265,600; XL,YL=155,230; X=1<<XL; Y=1<<YL; A0=1<<XS; B0=1<<YS
base=leak & ~(((1<<XL)-1)<<XS) & ~(((1<<YL)-1)<<YS)

def C(low,high): return base | (low<<150) | (high<<920)
def padd(p,q):
 r=p.copy()
 for m,c in q.items():
  c+=r.get(m,0)
  if c: r[m]=c
  elif m in r: del r[m]
 return r
def pmul(p,q,mod=None):
 r={}
 for (i,j),a in p.items():
  for (k,l),b in q.items():
   m=(i+k,j+l); r[m]=r.get(m,0)+a*b
 if mod: r={m:c%mod for m,c in r.items() if c%mod}
 return {m:c for m,c in r.items() if c}
def build(c,m,t,lead='x'):
 # monic f modulo N in lead variable; use representatives 0..N-1
 if lead=='x':
  inv=pow(A0,-1,N); f={(1,0):1,(0,1):(B0*inv)%N,(0,0):(c*inv)%N}
  sh=lambda k,i:(0,i)
 else:
  inv=pow(B0,-1,N); f={(0,1):1,(1,0):(A0*inv)%N,(0,0):(c*inv)%N}
  sh=lambda k,i:(i,0)
 pw=[{(0,0):1}]
 for k in range(1,m+1): pw.append(pmul(pw[-1],f, pow(N,k)))
 polys=[]
 for k in range(m+1):
  for i in range(m-k+1):
   sx,sy=sh(k,i); mult=pow(N,max(t-k,0))
   polys.append({(a+sx,b+sy):coef*mult for (a,b),coef in pw[k].items()})
 mons=sorted(set().union(*[p.keys() for p in polys]), key=lambda z:(z[0]+z[1],z[0],z[1]))
 sc=[(X**i)*(Y**j) for i,j in mons]; col={u:i for i,u in enumerate(mons)}
 M=IntegerMatrix(len(polys),len(mons))
 for r,p in enumerate(polys):
  for mon,coef in p.items(): M[r,col[mon]]=int(coef*sc[col[mon]])
 return M,mons,sc

def reduce(M,delta=.99,bkz=0):
 t0=time.time(); LLL.reduction(M,delta=delta,method='proved',float_type='mpfr',precision=256)
 if bkz:
  BKZ.reduction(M,BKZ.Param(block_size=bkz,max_loops=1))
 return time.time()-t0

def unscale(row,mons,sc):
 d={}
 for v,mon,s in zip(row,mons,sc):
  v=int(v)
  if v:
   if v%s: return None
   d[mon]=v//s
 return d

def evalp(p,x,y,mod=None):
 s=0
 for (i,j),c in p.items(): s += c*pow(x,i,mod or 1<<100000)*pow(y,j,mod or 1<<100000) if mod else c*(x**i)*(y**j)
 return s%mod if mod else s

def roots_mod(polys,pr):
 good=[]
 # precompute powers
 xp=[[1] for _ in range(pr)]; yp=[[1] for _ in range(pr)]
 max_i=max([i for p in polys for i,j in p] or [0]); max_j=max([j for p in polys for i,j in p] or [0])
 for a in range(pr):
  xp[a]=[1]*(max_i+1)
  for i in range(1,max_i+1): xp[a][i]=xp[a][i-1]*a%pr
 for b in range(pr):
  yp[b]=[1]*(max_j+1)
  for j in range(1,max_j+1): yp[b][j]=yp[b][j-1]*b%pr
 for a in range(pr):
  for b in range(pr):
   ok=True
   for p in polys:
    s=0
    for (i,j),c in p.items(): s=(s+(c%pr)*xp[a][i]*yp[b][j])%pr
    if s: ok=False; break
   if ok: good.append((a,b))
 return good

def crt_pair(state,mod,roots,pr,cap=20000):
 inv=pow(mod,-1,pr); out=[]
 for x,y in state:
  for a,b in roots:
   tx=((a-x)%pr)*inv%pr; ty=((b-y)%pr)*inv%pr
   out.append((x+mod*tx, y+mod*ty))
   if len(out)>cap: return [],0,False
 return out,mod*pr,True

def recover(polys,verbose=False):
 primes=[257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467]
 st=[(0,0)]; mod=1
 for pr in primes:
  roots=roots_mod(polys[:min(4,len(polys))],pr)
  if verbose: print('p',pr,'roots',len(roots),'state',len(st),'modbits',mod.bit_length(),flush=True)
  if not roots: return []
  st,mod,ok=crt_pair(st,mod,roots,pr)
  if not ok: return []
  if mod>X and mod>Y:
   return [(x,y) for x,y in st if x<X and y<Y]
 return []

def try_one(low,high,m=14,t=None,lead='x',bkz=0,lim=12):
 if t is None: t=max(1,round((1-math.sqrt(.5))*m))
 c=C(low,high); M,mons,sc=build(c,m,t,lead); dt=reduce(M,bkz=bkz)
 rows=[]
 for r in range(M.nrows):
  row=[int(M[r,j]) for j in range(M.ncols)]
  if not any(row): continue
  pol=unscale(row,mons,sc)
  if pol and len(pol)>1: rows.append((sum(abs(v) for v in row).bit_length(),len(pol),max(i+j for i,j in pol),pol))
 rows.sort(key=lambda z:z[0]); polys=[z[3] for z in rows[:lim]]
 print('cand',low,high,'m,t',m,t,'dim',M.nrows,M.ncols,'LLL',round(dt,2),'bits',[z[0] for z in rows[:6]],flush=True)
 cands=recover(polys)
 for x,y in cands:
  p=c + (x<<XS)+(y<<YS)
  if p>1 and N%p==0 and (p&mask)==leak:
   print('FOUND',p); print(decrypt_from_p(p)); return p
 return None
if __name__=='__main__':
 m=int(sys.argv[1]) if len(sys.argv)>1 else 14; t=int(sys.argv[2]) if len(sys.argv)>2 else None
 for high in range(16):
  for low in range(16):
   if try_one(low,high,m,t): raise SystemExit
