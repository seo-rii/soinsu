import sys,itertools,time
sys.path.append('/mnt/data')
from hm2_exact import build,lll,unscale,bits_l1,roots_mod,crt,X,Y,C_for,N,mask,leak,decrypt_from_p,x_start,y_start

def lift(polys, primes=(31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257)):
    cand=[(0,0,1)]
    for p in primes:
        rs=roots_mod(polys,p)
        if not rs or len(rs)>12: return []
        nc=[]
        for x,y in rs:
          for a,b,M in cand:
            NM=M*p; nx=crt(a,M,x,p); ny=crt(b,M,y,p)
            if NM>X and nx>=X: continue
            if NM>Y and ny>=Y: continue
            nc.append((nx,ny,NM))
            if len(nc)>200: return []
        cand=nc
        if not cand: return []
        if cand[0][2] > max(X,Y): return cand
    return cand

def check(low,high,m=8,t=2):
    C=C_for(low,high); B,mons,sc=build(C,m,t); R,dt=lll(B,0.99)
    polys=[]
    for r in R:
      P=unscale(r,mons,sc)
      if P and len(P)>=2 and max(i+j for i,j in P)>0:
        polys.append(P)
      if len(polys)>=10: break
    pairs=list(itertools.combinations(range(min(len(polys),8)),2))
    # try also triples from first 6
    sets=[(i,j) for i,j in pairs] + list(itertools.combinations(range(min(len(polys),6)),3))
    for S in sets:
      cand=lift([polys[i] for i in S])
      for x,y,M in cand:
        if x<X and y<Y:
          p=C+(x<<x_start)+(y<<y_start)
          if (p&mask)==leak and N%p==0:
            return p,dt,S
    return None,dt,None

st=time.time(); total=0
for high in range(16):
  for low in range(16):
    p,dt,S=check(low,high,8,2); total+=dt
    print('done',high,low,'lll',round(dt,2),'elapsed',round(time.time()-st,1),'found',bool(p),flush=True)
    if p:
      print('FOUND',low,high,S,hex(p)); print(decrypt_from_p(p)); raise SystemExit
print('not found',time.time()-st,'lll_total',total)
