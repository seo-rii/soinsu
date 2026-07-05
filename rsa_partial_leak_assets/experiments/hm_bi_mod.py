from fpylll import IntegerMatrix, LLL
import sympy as sp, itertools, time, sys, math

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
PM=H('''ffa360d46885c534 d186538170633faf c2c0548a2e24a2c1 c0000000000039e2
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

def poly_mul(a,b):
    r={}
    for (i,j),ai in a.items():
        for (k,l),bj in b.items():
            e=(i+k,j+l); r[e]=r.get(e,0)+ai*bj
    return {e:c for e,c in r.items() if c}
def poly_pow(f,k):
    r={(0,0):1}
    for _ in range(k): r=poly_mul(r,f)
    return r

def build_polys(low,high,m=6,t=3,keep=8):
    p0=PM | (low<<150) | (high<<920)
    inv=pow(1<<265,-1,N); A=((1<<600)*inv)%N; C=(p0*inv)%N
    f={(1,0):1,(0,1):A,(0,0):C}
    shifts=[]
    for k in range(m+1):
        base=poly_pow(f,k); mult=pow(N,max(t-k,0))
        for iy in range(m+1-k): shifts.append({(i,j+iy):c*mult for (i,j),c in base.items()})
    mons=sorted(set(e for g in shifts for e in g), key=lambda e:(e[0]+e[1],e[1],e[0]))
    X,Y=1<<155,1<<230; scales=[pow(X,i)*pow(Y,j) for i,j in mons]
    B=IntegerMatrix(len(shifts),len(mons))
    for r,g in enumerate(shifts):
        for c,e in enumerate(mons): B[r,c]=g.get(e,0)*scales[c]
    LLL.reduction(B, method='proved', delta=0.99)
    polys=[]
    for r in range(min(keep,B.nrows)):
        d={}
        for c,e in enumerate(mons):
            v=int(B[r,c])
            if v: d[e]=v//scales[c]
        if len(d)>1: polys.append(d)
    return p0,polys

def eval_mod(poly,xv,yv,P):
    s=0
    for (i,j),c in poly.items(): s=(s+(c%P)*pow(xv,i,P)*pow(yv,j,P))%P
    return s

def to_expr(poly, P, x, y):
    return sum((c%P)*x**i*y**j for (i,j),c in poly.items())

def lin_roots_univar(poly,P,var):
    # factor and return roots of all linear factors over GF(P)
    roots=[]
    try:
        fl=sp.factor_list(poly, modulus=P)[1]
    except Exception:
        return roots
    for fac,mul in fl:
        fp=sp.Poly(fac,var, modulus=P)
        if fp.degree()==1:
            a=int(fp.nth(1))%P; b=int(fp.nth(0))%P
            if a: roots.append((-b*pow(a,-1,P))%P)
    return list(set(roots))

def roots_mod_pair(p1,p2,P):
    x,y=sp.symbols('x y')
    f=sp.Poly(to_expr(p1,P,x,y), x,y, modulus=P)
    g=sp.Poly(to_expr(p2,P,x,y), x,y, modulus=P)
    try:
        R=sp.resultant(f.as_expr(), g.as_expr(), y)
        Rx=sp.Poly(R, x, modulus=P)
        if Rx.is_zero: return []
        xs=lin_roots_univar(Rx.as_expr(),P,x)
    except Exception as e:
        return []
    out=[]
    for xv in xs[:50]:
        fy=sp.Poly(f.as_expr().subs(x,xv), y, modulus=P)
        gy=sp.Poly(g.as_expr().subs(x,xv), y, modulus=P)
        try:
            h=sp.gcd(fy,gy)
            ys=lin_roots_univar(h.as_expr(),P,y) if h.degree()>0 else []
        except Exception:
            ys=[]
        for yv in ys:
            if eval_mod(p1,xv,yv,P)==0 and eval_mod(p2,xv,yv,P)==0: out.append((xv,yv))
    return list(set(out))

def crt_pair(a,m,b,n):
    return (a + ((b-a)%n)*pow(m,-1,n)%n*m)%(m*n), m*n

def recover(low,high,m=6,t=3,keep=8):
    p0,polys=build_polys(low,high,m,t,keep)
    if len(polys)<2: return None
    X,Y=1<<155,1<<230
    # try several pairs; use 30-bit primes
    primes=[1000003,1000033,1000037,1000039,1000081,1000099,1000117,1000121,1000133,1000151,1000159,1000171,1000183,1000187,1000193,1000199]
    for ia,ib in itertools.combinations(range(min(len(polys),6)),2):
        states=[(0,1,0,1)] # x,mx,y,my
        ok=True
        for P in primes:
            rts=roots_mod_pair(polys[ia],polys[ib],P)
            if not rts:
                ok=False; break
            ns=[]
            for xr,mx,yr,my in states:
                for xP,yP in rts[:20]:
                    if math.gcd(mx,P)!=1 or math.gcd(my,P)!=1: continue
                    nx,nmx=crt_pair(xr,mx,xP,P); ny,nmy=crt_pair(yr,my,yP,P)
                    # keep if feasible modulo bounds (there exists value <bound with residue)
                    if nx<X or nmx<X:
                        if ny<Y or nmy<Y:
                            ns.append((nx,nmx,ny,nmy))
            # dedup and cap
            ded={ (a,b,c,d) for a,b,c,d in ns }
            states=list(ded)[:100]
            if not states: ok=False; break
            if states[0][1]>X and states[0][3]>Y:
                for xr,mx,yr,my in states:
                    if xr<X and yr<Y:
                        p=p0+(xr<<265)+(yr<<600)
                        if (p&MASK)==PM and N%p==0: return p
                break
        # maybe test final states
        for xr,mx,yr,my in states:
            if xr<X and yr<Y:
                p=p0+(xr<<265)+(yr<<600)
                if (p&MASK)==PM and N%p==0: return p
    return None

if __name__=='__main__':
    print('start')
    M=int(sys.argv[1]) if len(sys.argv)>1 else 6; T=int(sys.argv[2]) if len(sys.argv)>2 else 3
    for h in range(16):
      for l in range(16):
        st=time.time(); p=recover(l,h,M,T,8); print(h,l,round(time.time()-st,2),bool(p), flush=True)
        if p:
            print(hex(p)); q=N//p; print(hex(q));
            d=pow(65537,-1,(p-1)*(q-1)); msg=pow(ct,d,N).to_bytes(256,'big').lstrip(b'\0'); print(msg, msg.hex()); sys.exit()
