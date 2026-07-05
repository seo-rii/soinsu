import sys,time,math
sys.path.append('/mnt/data')
import hm_try
from fpylll import LLL
from sympy import primerange

def get_hs(low,high,m=8,t=2,rows=6):
    xmask=((1<<155)-1)<<265; ymask=((1<<230)-1)<<600
    const=(hm_try.PAND | (low<<150) | (high<<920)) & ~xmask & ~ymask
    inv=pow(1<<265,-1,hm_try.N); A=((1<<600)*inv)%hm_try.N; B=(const*inv)%hm_try.N
    # use residues 0..N-1 per HM; signed maybe not? Use same signed as before? try both maybe
    if A>hm_try.N//2: A-=hm_try.N
    if B>hm_try.N//2: B-=hm_try.N
    f={(1,0):1,(0,1):A,(0,0):B}; X=1<<155; Y=1<<230
    M,mons=hm_try.hm_lattice(f,X,Y,hm_try.N,m,t)
    LLL.reduction(M)
    hs=[hm_try.rowpoly([M[r,c] for c in range(M.ncols)], mons, X,Y) for r in range(min(rows,M.nrows))]
    return const,hs

def roots_mod_brute(polys,p,max_roots=10000):
    # reduce coefficients
    ps=[]; maxx=maxy=0
    for h in polys:
        d={}
        for (i,j),c in h.items():
            cc=c%p
            if cc: d[(i,j)]=cc; maxx=max(maxx,i); maxy=max(maxy,j)
        ps.append(d)
    xpows=[[1]*p for _ in range(maxx+1)]
    ypows=[[1]*p for _ in range(maxy+1)]
    for i in range(1,maxx+1):
        prev=xpows[i-1]; xpows[i]=[(prev[x]*x)%p for x in range(p)]
    for j in range(1,maxy+1):
        prev=ypows[j-1]; ypows[j]=[(prev[y]*y)%p for y in range(p)]
    roots=[]
    for y in range(p):
        # precompute each polynomial values over x
        ok=[True]*p
        for h in ps:
            vals=[0]*p
            for (i,j),c in h.items():
                cy=(c*ypows[j][y])%p
                if cy:
                    xi=xpows[i]
                    for x in range(p): vals[x]=(vals[x]+cy*xi[x])%p
            for x,v in enumerate(vals):
                if v: ok[x]=False
        for x in range(p):
            if ok[x]:
                roots.append((x,y))
                if len(roots)>max_roots: return roots
    return roots

def crt_pair(a,m,b,n):
    # moduli coprime
    return (a + ((b-a)*pow(m,-1,n)%n)*m)%(m*n), m*n

def lift_roots(polys, primes, max_states=200):
    states=[(0,0,1)]
    for p in primes:
        roots=roots_mod_brute(polys,p,max_roots=1000)
        if not roots: return []
        new=[]
        for xr,yr,mod in states:
            invm=pow(mod,-1,p)
            for a,b in roots:
                xnew=(xr+((a-xr)%p)*invm%p*mod)%(mod*p)
                ynew=(yr+((b-yr)%p)*invm%p*mod)%(mod*p)
                new.append((xnew,ynew,mod*p))
                if len(new)>max_states: break
            if len(new)>max_states: break
        states=new
        # print('p',p,'roots',len(roots),'states',len(states),'modbits',states[0][2].bit_length() if states else 0)
    return states

if __name__=='__main__':
    const,hs=get_hs(0,0,8,2,rows=5)
    for p in [101,103,107,109,113,127,131,137,139,149]:
        st=time.time(); r=roots_mod_brute(hs[:3],p,max_roots=1000); print(p,len(r),time.time()-st,r[:5])
