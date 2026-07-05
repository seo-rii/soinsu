import sys, math, time, itertools, sympy as sp
from fpylll import IntegerMatrix, LLL

def H(s): return int(''.join(s.split()),16)
N=H('''e505004fb5d34eb7 12d48ff4bbe8d27f c388133c6c0e7340 01061c0ee0a4edc6
37c04fe8dd376185 de8ba04d0ccdbabb 93ab7c371b88d92e 865eec42b028c61d
d7004ebf2ebb5d69 d0a09142be5c9de4 da16e514eea31817 2ecda6cd192073eb
afb1e02d522ec053 34590ea6d75960c4 937bf64f9700db17 7a4aa3da6aae6807
e5e32c0d0e428a0d b68d299f20c235d8 4ef459b0cf118286 59c31663c9ea8204
4b28152c89a9c36c 3ec4303bd36664fd 77fb02c58340bdae 21120326d83fc017
34bc90048dec9fe3 5f08c8fdc523abf8 4a91ec430f495672 37c3153a2035ff62
5613b6dc3e6cb14d 50e18b8a79b25d67 8465b3ad02f5b7d8 18a1e2d635a0baf1''')
MASK=H('''ffffffffffffffff fffffffff0ffffff ffffffffffffffff c00000000000fffe
0000000000000000 000003ffe0000000 0000000000ffffff ffffffffffffffff
ffffffffffffffff fffffff000000000 000003ffe0000000 00000000000001ff
ffffffffffffffff fffffffffc3fffff ffffffffffffffff ffffffffffffffff''')
LEAK=H('''ffa360d46885c534 d186538170633faf c2c0548a2e24a2c1 c0000000000039e2
0000000000000000 000000a520000000 00000000003e2de4 c436d2ca740a6246
99e1a1af94045c63 261323c000000000 000003bba0000000 00000000000000e5
0b0bc2461fcbac07 26360c2c0809450a 9a892cbf1d98ceee 48827591ccc593c9''')
ct=H('''8919342826ef3821 5af31e00c9290c4c 50ef9ff9e1afc591 47fab5b096361035
e85f5fc95b73b069 7813b57b831a807d 41bcbecde5b9e663 9e2845b14e395ed0
e5d995e63709ac0c 5ee2337228ee76bc bad857b14904aa2e 8e9997671908a634
d0d1dda1d062ce7f 2e3293ddec8f5cce 26029292d594a062 dcf317d2a8380f43
d72551889efceb87 6c8945a50382272e 76ed6b6fcdff1603 44e9e948e2b6e740
e78bedf25f30e2c7 eeb5f74686c8eadc 29cea04ff08cfd86 dfd3d2a1632bf04a
d5cfa369892a2da4 0f0dc0098ce6b731 d841aab3d0c8b78e b69c4625c47c4ad7
158d49bb5d879581 e02bc525abe47f39 f699864bc5ce1de7 19430dae7aa5480b''')
e=65537
XS,XL,YS,YL=265,155,600,230
X,Y=1<<XL,1<<YL; A,B=1<<XS,1<<YS
base=LEAK & ~(((1<<XL)-1)<<XS) & ~(((1<<YL)-1)<<YS)

def mul(p,q,mod=None):
    r={}
    for (i,j),a in p.items():
        for (k,l),b in q.items():
            c=a*b
            if mod: c%=mod
            if c:
                e=(i+k,j+l); r[e]=r.get(e,0)+c
    if mod: r={e:c%mod for e,c in r.items() if c%mod}
    return r

def build(C,m,t,both=True):
    polys=[]
    # base normalized in x and y, coefficients reduced modulo N
    bases=[]
    invA=pow(A,-1,N); bases.append({(1,0):1,(0,1):(B*invA)%N,(0,0):(C*invA)%N}) # monic x; shift y^j
    if both:
        invB=pow(B,-1,N); bases.append({(0,1):1,(1,0):(A*invB)%N,(0,0):(C*invB)%N}) # monic y; shift x^j
    for bi,f in enumerate(bases):
        pows=[{(0,0):1}]
        for k in range(1,m+1):
            # Important: kionactf/sage does not reduce f^k mod N; it changes mod coeffs to ZZ then powers.
            # Coefficients are representatives modulo N; multiplication over ZZ.
            pows.append(mul(pows[-1], f, None))
        for k in range(m+1):
            mult=pow(N, max(t-k,0))
            for j in range(m-k+1):
                d={}
                for (u,v),c in pows[k].items():
                    if bi==0: mon=(u,v+j)      # monic x: exclude x, shift y
                    else:     mon=(u+j,v)      # monic y: exclude y, shift x
                    d[mon]=d.get(mon,0)+c*mult
                polys.append(d)
    # dedup rows by dictionary contents (k=0 rows duplicate between bases)
    uniq=[]; seen=set()
    for p in polys:
        key=tuple(sorted(p.items()))
        if key not in seen:
            seen.add(key); uniq.append(p)
    polys=uniq
    mons=sorted(set(e for p in polys for e in p), key=lambda z:(z[0]+z[1],z[0],z[1]))
    scales=[(X**i)*(Y**j) for i,j in mons]
    col={e:i for i,e in enumerate(mons)}
    M=IntegerMatrix(len(polys), len(mons))
    for r,p in enumerate(polys):
        for mon,c in p.items():
            M[r,col[mon]]=int(c*scales[col[mon]])
    return M,mons,scales

def unscale_row(M,r,mons,scales):
    d={}
    for j,mon in enumerate(mons):
        v=int(M[r,j])
        if v:
            s=scales[j]
            # after LLL, divisibility may fail due to linearly combined rows? It should hold over ZZ if rows all scaled same cols.
            if v % s: return None
            d[mon]=v//s
    return d

def eval_poly(poly,xv,yv):
    return sum(c*pow(xv,i)*pow(yv,j) for (i,j),c in poly.items())

def eval_mod(poly,xv,yv,P):
    return sum((c%P)*pow(xv,i,P)*pow(yv,j,P) for (i,j),c in poly.items())%P

# fast brute roots mod small prime, useful when primes ~ 257..1009
_mod_cache={}

def eval_mod(poly,xv,yv,P):
    return sum((c%P)*pow(xv,i,P)*pow(yv,j,P) for (i,j),c in poly.items())%P

def to_expr(poly,P,x,y):
    return sum((c%P)*x**i*y**j for (i,j),c in poly.items())

def lin_roots(expr,P,var):
    try:
        facs=sp.factor_list(expr, modulus=P)[1]
    except Exception:
        return []
    out=[]
    for fac,mul in facs:
        fp=sp.Poly(fac,var,modulus=P)
        if fp.degree()==1:
            a=int(fp.nth(1))%P; b=int(fp.nth(0))%P
            if a: out.append((-b*pow(a,-1,P))%P)
    return list(dict.fromkeys(out))

def roots_mod_pair(p1,p2,P,validators=()):
    x,y=sp.symbols('x y')
    f=sp.Poly(to_expr(p1,P,x,y),x,y,modulus=P)
    g=sp.Poly(to_expr(p2,P,x,y),x,y,modulus=P)
    try:
        R=sp.resultant(f.as_expr(),g.as_expr(),y)
        Rx=sp.Poly(R,x,modulus=P)
        if Rx.is_zero: return []
        xs=lin_roots(Rx.as_expr(),P,x)
    except Exception as e:
        return []
    out=[]
    for xv in xs[:200]:
        try:
            fy=sp.Poly(f.as_expr().subs(x,xv), y, modulus=P)
            gy=sp.Poly(g.as_expr().subs(x,xv), y, modulus=P)
            h=sp.gcd(fy,gy)
            ys=lin_roots(h.as_expr(),P,y) if h.degree()>0 else []
        except Exception:
            ys=[]
        for yv in ys:
            if eval_mod(p1,xv,yv,P)==0 and eval_mod(p2,xv,yv,P)==0 and all(eval_mod(v,xv,yv,P)==0 for v in validators):
                out.append((xv,yv))
    return list(dict.fromkeys(out))

def crt_pair(a,m,b,n):
    return (a + m*(((b-a)%n)*pow(m%n,-1,n)%n))%(m*n), m*n

def recover_from_polys(polys,C,maxpair=28):
    xmask=(MASK>>XS)&(X-1); xleak=(LEAK>>XS)&(X-1)
    ymask=(MASK>>YS)&(Y-1); yleak=(LEAK>>YS)&(Y-1)
    primes=[1000003,1000033,1000037,1000039,1000081,1000099,1000117,1000121,1000133,1000151,1000159,1000171,1000183,1000187,1000193,1000199,1000211,1000213,1000231,1000249,1000253,1000273,1000289]
    L=min(len(polys),10); pairs=list(itertools.combinations(range(L),2))[:maxpair]
    for ia,ib in pairs:
        vals=[polys[k] for k in range(L) if k not in (ia,ib)][:4]
        states=[(0,1,0,1)]
        for P in primes:
            R=roots_mod_pair(polys[ia],polys[ib],P,validators=vals[:2])
            if not R or len(R)>64:
                states=[]; break
            ns=[]
            for xr,mx,yr,my in states:
                for xp,yp in R:
                    nx,nmx=crt_pair(xr,mx,xp,P); ny,nmy=crt_pair(yr,my,yp,P)
                    if (nx<X or nmx<X) and (ny<Y or nmy<Y): ns.append((nx,nmx,ny,nmy))
                    if len(ns)>5000: break
                if len(ns)>5000: break
            states=list(dict.fromkeys(ns))
            if not states: break
            if states[0][1]>=X and states[0][3]>=Y:
                for xr,mx,yr,my in states:
                    if xr<X and yr<Y and (xr&xmask)==xleak and (yr&ymask)==yleak:
                        p=C+(xr<<XS)+(yr<<YS)
                        if (p&MASK)==LEAK and N%p==0: return p
                break
    return None

def one(cid,m,t,both=True,keep=28):
    lo,hi=cid&15,cid>>4; C=base | (lo<<150) | (hi<<920)
    M,mons,scales=build(C,m,t,both); st=time.time()
    LLL.reduction(M,delta=0.99,method='proved',float_type='mpfr',precision=256)
    dt=time.time()-st
    rows=[]
    for r in range(M.nrows):
        l1=sum(abs(int(M[r,j])) for j in range(M.ncols))
        if not l1: continue
        pol=unscale_row(M,r,mons,scales)
        if pol and len(pol)>1: rows.append((l1.bit_length(),pol))
    rows=sorted(rows,key=lambda z:z[0])[:keep]
    print(f'cid={cid:03d} lo={lo:x} hi={hi:x} dim={M.nrows}x{M.ncols} LLL={dt:.2f}s bits={[b for b,_ in rows[:6]]}',flush=True)
    return recover_from_polys([p for _,p in rows],C)

def dec(p):
    q=N//p; d=pow(e,-1,(p-1)*(q-1)); m=pow(ct,d,N)
    return m.to_bytes((m.bit_length()+7)//8,'big')
if __name__=='__main__':
    m=int(sys.argv[1]) if len(sys.argv)>1 else 8; t=int(sys.argv[2]) if len(sys.argv)>2 else 3
    a=int(sys.argv[3]) if len(sys.argv)>3 else 0; b=int(sys.argv[4]) if len(sys.argv)>4 else 256
    both=bool(int(sys.argv[5])) if len(sys.argv)>5 else True
    for cid in range(a,b):
        p=one(cid,m,t,both)
        if p:
            print('FOUND'); print('p=',hex(p)); print('q=',hex(N//p)); print('m=',dec(p)); print('mhex=',dec(p).hex()); sys.exit()
    print('not found')
