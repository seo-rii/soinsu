# Exact Herrmann-May bivariate attempt for problem 7.  Requires: pip install fpylll
# Usage: python3 hm2_exact_compact.py 8 3 0 256   # scans candidate ids [0,256)
import sys,math,time
from fpylll import IntegerMatrix,LLL
N=int('e505004fb5d34eb712d48ff4bbe8d27fc388133c6c0e734001061c0ee0a4edc637c04fe8dd376185de8ba04d0ccdbabb93ab7c371b88d92e865eec42b028c61dd7004ebf2ebb5d69d0a09142be5c9de4da16e514eea318172ecda6cd192073ebafb1e02d522ec05334590ea6d75960c4937bf64f9700db177a4aa3da6aae6807e5e32c0d0e428a0db68d299f20c235d84ef459b0cf11828659c31663c9ea82044b28152c89a9c36c3ec4303bd36664fd77fb02c58340bdae21120326d83fc01734bc90048dec9fe35f08c8fdc523abf84a91ec430f49567237c3153a2035ff625613b6dc3e6cb14d50e18b8a79b25d678465b3ad02f5b7d818a1e2d635a0baf1',16)
ct=int('8919342826ef38215af31e00c9290c4c50ef9ff9e1afc59147fab5b096361035e85f5fc95b73b0697813b57b831a807d41bcbecde5b9e6639e2845b14e395ed0e5d995e63709ac0c5ee2337228ee76bcbad857b14904aa2e8e9997671908a634d0d1dda1d062ce7f2e3293ddec8f5cce26029292d594a062dcf317d2a8380f43d72551889efceb876c8945a50382272e76ed6b6fcdff160344e9e948e2b6e740e78bedf25f30e2c7eeb5f74686c8eadc29cea04ff08cfd86dfd3d2a1632bf04ad5cfa369892a2da40f0dc0098ce6b731d841aab3d0c8b78eb69c4625c47c4ad7158d49bb5d879581e02bc525abe47f39f699864bc5ce1de719430dae7aa5480b',16)
mask=int('fffffffffffffffffffffffff0ffffffffffffffffffffffc00000000000fffe0000000000000000000003ffe00000000000000000fffffffffffffffffffffffffffffffffffffff000000000000003ffe00000000000000000000001fffffffffffffffffffffffffc3fffffffffffffffffffffffffffffffff',16)
leak=int('ffa360d46885c534d186538170633fafc2c0548a2e24a2c1c0000000000039e20000000000000000000000a52000000000000000003e2de4c436d2ca740a624699e1a1af94045c63261323c000000000000003bba0000000000000000000000e50b0bc2461fcbac0726360c2c0809450a9a892cbf1d98ceee48827591ccc593c9',16)
e=65537; XS,XL,YS,YL=265,155,600,230; X,Y=1<<XL,1<<YL; A,B=1<<XS,1<<YS
base=leak & ~(((1<<XL)-1)<<XS) & ~(((1<<YL)-1)<<YS)
def mul(p,q,mod=0):
 r={}
 for (i,j),a in p.items():
  for (k,l),b in q.items(): r[(i+k,j+l)]=r.get((i+k,j+l),0)+a*b
 return {m:(c%mod if mod else c) for m,c in r.items() if (c%mod if mod else c)}
def build(C,m,t):                 # f monic in y; shifts in smaller x
 inv=pow(B,-1,N); f={(0,1):1,(1,0):(A*inv)%N,(0,0):(C*inv)%N}; pw=[{(0,0):1}]
 for k in range(1,m+1): pw.append(mul(pw[-1],f,N**k))
 P=[]
 for k in range(m+1):
  for i in range(m-k+1): P.append({(a+i,b):c*(N**max(t-k,0)) for (a,b),c in pw[k].items()})
 Mns=sorted(set().union(*(p.keys() for p in P)),key=lambda z:(z[0]+z[1],z[0])); sc=[X**i*Y**j for i,j in Mns]; col={u:i for i,u in enumerate(Mns)}
 M=IntegerMatrix(len(P),len(Mns))
 for r,p in enumerate(P):
  for mon,c in p.items(): M[r,col[mon]]=int(c*sc[col[mon]])
 return M,Mns,sc
def unsc(row,mn,sc):
 d={}
 for v,mon,s in zip(row,mn,sc):
  v=int(v)
  if v:
   if v%s: return None
   d[mon]=v//s
 return d
def roots_mod(P,pr):
 mi=max(i for p in P for i,j in p); mj=max(j for p in P for i,j in p); xp=[[1]*(mi+1) for _ in range(pr)]; yp=[[1]*(mj+1) for _ in range(pr)]
 for a in range(pr):
  for i in range(1,mi+1): xp[a][i]=xp[a][i-1]*a%pr
 for b in range(pr):
  for j in range(1,mj+1): yp[b][j]=yp[b][j-1]*b%pr
 R=[]
 for a in range(pr):
  for b in range(pr):
   if all(sum((c%pr)*xp[a][i]*yp[b][j] for (i,j),c in p.items())%pr==0 for p in P): R.append((a,b))
 return R
def crt(st,mod,R,pr,cap=30000):
 out=[]; inv=pow(mod,-1,pr)
 for x,y in st:
  for a,b in R:
   out.append((x+mod*(((a-x)%pr)*inv%pr), y+mod*(((b-y)%pr)*inv%pr)))
   if len(out)>cap: return [],0
 return out,mod*pr
def rec(P):
 st=[(0,0)]; mod=1
 for pr in [257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467]:
  R=roots_mod(P[:4],pr)
  if not R: return []
  st,mod=crt(st,mod,R,pr)
  if not st: return []
  if mod>X and mod>Y: return [(x,y) for x,y in st if x<X and y<Y]
 return []
def dec(p):
 q=N//p; d=pow(e,-1,(p-1)*(q-1)); z=pow(ct,d,N); return z.to_bytes((z.bit_length()+7)//8,'big')
def one(cid,m,t):
 lo,hi=cid&15,cid>>4; C=base|(lo<<150)|(hi<<920); M,mn,sc=build(C,m,t); t0=time.time(); LLL.reduction(M,delta=.99,method='proved',float_type='mpfr',precision=256); dt=time.time()-t0
 rows=[]
 for r in range(M.nrows):
  row=[int(M[r,j]) for j in range(M.ncols)]; p=unsc(row,mn,sc)
  if p and len(p)>1: rows.append((sum(abs(v) for v in row).bit_length(),p))
 rows=sorted(rows,key=lambda z:z[0])[:10]; print(cid,lo,hi,'dim',M.nrows,'LLL',round(dt,2),'bits',[b for b,_ in rows[:5]],flush=True)
 for x,y in rec([p for _,p in rows]):
  p=C+(x<<XS)+(y<<YS)
  if N%p==0 and (p&mask)==leak: print('FOUND p=',p,'\nmsg=',dec(p)); return True
 return False
if __name__=='__main__':
 m=int(sys.argv[1]) if len(sys.argv)>1 else 8; t=int(sys.argv[2]) if len(sys.argv)>2 else max(1,round((1-math.sqrt(.5))*m)); a=int(sys.argv[3]) if len(sys.argv)>3 else 0; b=int(sys.argv[4]) if len(sys.argv)>4 else 256
 for cid in range(a,b):
  if one(cid,m,t): break
